#!/usr/bin/env python3
"""
Complete QR Code Generator from Scratch

A comprehensive implementation of QR code generation without external libraries,
following the ISO/IEC 18004 specification.

This file accompanies the LaTeX tutorial document and includes:
- Galois Field GF(256) arithmetic
- Polynomial operations over GF(256)
- Reed-Solomon error correction encoding
- Data encoding (numeric, alphanumeric, byte modes)
- BCH code for format information
- QR code matrix construction with function patterns
- Data masking with penalty calculation
- Complete QR code generation

Author: Tutorial Implementation
Based on: ISO/IEC 18004, Thonky's QR Code Tutorial, and Wikiversity Reed-Solomon Guide
"""

from typing import List, Tuple, Optional, Callable

#==============================================================================
# GALOIS FIELD GF(256) ARITHMETIC
#==============================================================================

class GF256:
    """
    Galois Field GF(2^8) arithmetic for QR codes.
    
    Uses the primitive polynomial: x^8 + x^4 + x^3 + x^2 + 1 (0x11d = 285)
    with generator alpha = 2.
    
    References:
    - https://en.wikipedia.org/wiki/Finite_field_arithmetic
    - https://research.swtch.com/field
    """
    
    PRIMITIVE_POLY = 0x11d  # x^8 + x^4 + x^3 + x^2 + 1 = 285
    
    def __init__(self):
        """Initialize the field with precomputed exp and log tables."""
        self.exp_table = [0] * 512  # Extended for convenience
        self.log_table = [0] * 256
        self._build_tables()
    
    def _build_tables(self):
        """Build exponential and logarithm lookup tables using alpha = 2."""
        x = 1
        for i in range(255):
            self.exp_table[i] = x
            self.exp_table[i + 255] = x  # Duplicate for easy modulo
            self.log_table[x] = i
            
            # Multiply by alpha (2) with reduction
            x = self._multiply_no_table(x, 2)
        
        self.log_table[0] = -1  # log(0) is undefined
    
    def _multiply_no_table(self, a: int, b: int) -> int:
        """
        Multiply two GF(256) elements without using tables.
        Uses Russian peasant multiplication with polynomial reduction.
        """
        result = 0
        while b > 0:
            if b & 1:  # If lowest bit is set
                result ^= a  # Add (XOR) a to result
            b >>= 1
            a <<= 1
            if a & 0x100:  # If degree >= 8
                a ^= self.PRIMITIVE_POLY  # Reduce modulo primitive
        return result
    
    def add(self, a: int, b: int) -> int:
        """Addition in GF(256) is XOR."""
        return a ^ b
    
    def subtract(self, a: int, b: int) -> int:
        """Subtraction in GF(256) is the same as addition."""
        return a ^ b
    
    def multiply(self, a: int, b: int) -> int:
        """Multiply two GF(256) elements using log tables."""
        if a == 0 or b == 0:
            return 0
        return self.exp_table[self.log_table[a] + self.log_table[b]]
    
    def divide(self, a: int, b: int) -> int:
        """Divide a by b in GF(256)."""
        if b == 0:
            raise ZeroDivisionError("Division by zero in GF(256)")
        if a == 0:
            return 0
        return self.exp_table[(self.log_table[a] - self.log_table[b]) % 255]
    
    def power(self, a: int, n: int) -> int:
        """Raise a to the power n in GF(256)."""
        if a == 0:
            return 0 if n > 0 else 1
        return self.exp_table[(self.log_table[a] * n) % 255]
    
    def inverse(self, a: int) -> int:
        """Find multiplicative inverse of a in GF(256)."""
        if a == 0:
            raise ZeroDivisionError("No inverse for 0")
        # a^(-1) = a^254 since a^255 = 1
        return self.exp_table[255 - self.log_table[a]]


# Global GF256 instance
gf = GF256()


#==============================================================================
# POLYNOMIAL OPERATIONS OVER GF(256)
#==============================================================================

class Polynomial:
    """
    Polynomial with coefficients in GF(256).
    
    Coefficients are stored in ascending order of degree:
    coeffs[i] is the coefficient of x^i.
    """
    
    def __init__(self, coefficients: List[int], gf_instance: GF256 = None):
        """Initialize polynomial with coefficients."""
        self.gf = gf_instance or gf
        # Trim leading zeros (from the highest degree end)
        self.coeffs = list(coefficients)
        while len(self.coeffs) > 1 and self.coeffs[-1] == 0:
            self.coeffs.pop()
    
    @property
    def degree(self) -> int:
        """Return the degree of the polynomial."""
        return len(self.coeffs) - 1
    
    def __repr__(self) -> str:
        terms = []
        for i, c in enumerate(self.coeffs):
            if c != 0:
                if i == 0:
                    terms.append(f"{c}")
                elif i == 1:
                    terms.append(f"{c}x")
                else:
                    terms.append(f"{c}x^{i}")
        return " + ".join(terms) if terms else "0"
    
    def evaluate(self, x: int) -> int:
        """Evaluate polynomial at x using Horner's method."""
        result = 0
        for coeff in reversed(self.coeffs):
            result = self.gf.add(self.gf.multiply(result, x), coeff)
        return result
    
    def add(self, other: 'Polynomial') -> 'Polynomial':
        """Add two polynomials."""
        max_len = max(len(self.coeffs), len(other.coeffs))
        a = self.coeffs + [0] * (max_len - len(self.coeffs))
        b = other.coeffs + [0] * (max_len - len(other.coeffs))
        result = [self.gf.add(a[i], b[i]) for i in range(max_len)]
        return Polynomial(result, self.gf)
    
    def multiply(self, other: 'Polynomial') -> 'Polynomial':
        """Multiply two polynomials."""
        result = [0] * (len(self.coeffs) + len(other.coeffs) - 1)
        for i, a in enumerate(self.coeffs):
            for j, b in enumerate(other.coeffs):
                product = self.gf.multiply(a, b)
                result[i + j] = self.gf.add(result[i + j], product)
        return Polynomial(result, self.gf)
    
    def scale(self, scalar: int) -> 'Polynomial':
        """Multiply polynomial by a scalar."""
        result = [self.gf.multiply(c, scalar) for c in self.coeffs]
        return Polynomial(result, self.gf)


#==============================================================================
# REED-SOLOMON ERROR CORRECTION
#==============================================================================

class ReedSolomonEncoder:
    """
    Reed-Solomon encoder for QR code error correction.
    
    References:
    - https://en.wikipedia.org/wiki/Reed-Solomon_error_correction
    - https://en.wikiversity.org/wiki/Reed-Solomon_codes_for_coders
    """
    
    def __init__(self, gf_instance: GF256 = None):
        self.gf = gf_instance or gf
        self._generator_cache = {}
    
    def build_generator(self, num_ec_codewords: int) -> Polynomial:
        """
        Build generator polynomial for given number of EC codewords.
        
        g(x) = (x - alpha^0)(x - alpha^1)...(x - alpha^(n-1))
             = (x + alpha^0)(x + alpha^1)...(x + alpha^(n-1))
             
        In GF(256), subtraction equals addition.
        """
        if num_ec_codewords in self._generator_cache:
            return self._generator_cache[num_ec_codewords]
        
        # Start with g(x) = 1
        gen = Polynomial([1], self.gf)
        
        for i in range(num_ec_codewords):
            # Multiply by (x + alpha^i)
            # coeffs [alpha^i, 1] represents alpha^i + x
            factor = Polynomial([self.gf.exp_table[i], 1], self.gf)
            gen = gen.multiply(factor)
        
        self._generator_cache[num_ec_codewords] = gen
        return gen
    
    def encode(self, data: List[int], num_ec_codewords: int) -> List[int]:
        """
        Encode data bytes with Reed-Solomon error correction.
        
        Args:
            data: List of data bytes (integers 0-255)
            num_ec_codewords: Number of error correction codewords to generate
        
        Returns:
            List of error correction codewords
        """
        generator = self.build_generator(num_ec_codewords)
        return self._divide_for_remainder(data, generator, num_ec_codewords)
    
    def _divide_for_remainder(self, data: List[int], generator: Polynomial, 
                              num_ec: int) -> List[int]:
        """
        Compute remainder of message polynomial divided by generator.
        Uses the shift-register approach for efficiency.
        """
        result = list(data) + [0] * num_ec
        gen_coeffs = list(reversed(generator.coeffs))  # High degree first
        
        for i in range(len(data)):
            coeff = result[i]
            if coeff != 0:
                for j in range(len(gen_coeffs)):
                    result[i + j] ^= self.gf.multiply(gen_coeffs[j], coeff)
        
        return result[-num_ec:]


# Global encoder instance
rs_encoder = ReedSolomonEncoder()


#==============================================================================
# DATA ENCODING MODES
#==============================================================================

# Mode indicators (4-bit values)
MODE_NUMERIC = 0b0001
MODE_ALPHANUMERIC = 0b0010
MODE_BYTE = 0b0100
MODE_KANJI = 0b1000
MODE_TERMINATOR = 0b0000

# Alphanumeric character mapping
ALPHANUMERIC_TABLE = {
    '0': 0, '1': 1, '2': 2, '3': 3, '4': 4,
    '5': 5, '6': 6, '7': 7, '8': 8, '9': 9,
    'A': 10, 'B': 11, 'C': 12, 'D': 13, 'E': 14,
    'F': 15, 'G': 16, 'H': 17, 'I': 18, 'J': 19,
    'K': 20, 'L': 21, 'M': 22, 'N': 23, 'O': 24,
    'P': 25, 'Q': 26, 'R': 27, 'S': 28, 'T': 29,
    'U': 30, 'V': 31, 'W': 32, 'X': 33, 'Y': 34,
    'Z': 35, ' ': 36, '$': 37, '%': 38, '*': 39,
    '+': 40, '-': 41, '.': 42, '/': 43, ':': 44
}


def get_character_count_bits(version: int, mode: int) -> int:
    """Get the number of bits for the character count indicator."""
    if version <= 9:
        table = {MODE_NUMERIC: 10, MODE_ALPHANUMERIC: 9, 
                 MODE_BYTE: 8, MODE_KANJI: 8}
    elif version <= 26:
        table = {MODE_NUMERIC: 12, MODE_ALPHANUMERIC: 11, 
                 MODE_BYTE: 16, MODE_KANJI: 10}
    else:
        table = {MODE_NUMERIC: 14, MODE_ALPHANUMERIC: 13, 
                 MODE_BYTE: 16, MODE_KANJI: 12}
    return table.get(mode, 8)


def detect_mode(data: str) -> int:
    """Detect the most efficient encoding mode for the data."""
    if all(c.isdigit() for c in data):
        return MODE_NUMERIC
    if all(c in ALPHANUMERIC_TABLE for c in data):
        return MODE_ALPHANUMERIC
    return MODE_BYTE


def int_to_bits(value: int, length: int) -> List[int]:
    """Convert integer to list of bits with specified length."""
    return [(value >> (length - 1 - i)) & 1 for i in range(length)]


def bits_to_bytes(bits: List[int]) -> List[int]:
    """Convert list of bits to list of bytes."""
    # Pad to multiple of 8
    bits = list(bits)
    while len(bits) % 8 != 0:
        bits.append(0)
    
    bytes_list = []
    for i in range(0, len(bits), 8):
        byte = 0
        for j in range(8):
            byte = (byte << 1) | bits[i + j]
        bytes_list.append(byte)
    
    return bytes_list


def encode_numeric(data: str) -> List[int]:
    """Encode numeric data. Returns list of bits."""
    bits = []
    i = 0
    while i < len(data):
        if i + 3 <= len(data):
            value = int(data[i:i+3])
            bits.extend(int_to_bits(value, 10))
            i += 3
        elif i + 2 <= len(data):
            value = int(data[i:i+2])
            bits.extend(int_to_bits(value, 7))
            i += 2
        else:
            value = int(data[i])
            bits.extend(int_to_bits(value, 4))
            i += 1
    return bits


def encode_alphanumeric(data: str) -> List[int]:
    """Encode alphanumeric data. Returns list of bits."""
    bits = []
    i = 0
    while i < len(data):
        if i + 2 <= len(data):
            v1 = ALPHANUMERIC_TABLE[data[i]]
            v2 = ALPHANUMERIC_TABLE[data[i + 1]]
            value = 45 * v1 + v2
            bits.extend(int_to_bits(value, 11))
            i += 2
        else:
            value = ALPHANUMERIC_TABLE[data[i]]
            bits.extend(int_to_bits(value, 6))
            i += 1
    return bits


def encode_byte(data: str) -> List[int]:
    """Encode byte data. Returns list of bits."""
    bits = []
    data_bytes = data.encode('utf-8') if isinstance(data, str) else data
    for byte in data_bytes:
        bits.extend(int_to_bits(byte, 8))
    return bits


def encode_data(data: str, version: int, mode: int = None) -> List[int]:
    """
    Encode data for QR code.
    
    Returns: List of bits including mode indicator and character count.
    """
    if mode is None:
        mode = detect_mode(data)
    
    bits = []
    
    # Mode indicator (4 bits)
    bits.extend(int_to_bits(mode, 4))
    
    # Character count indicator
    count_bits = get_character_count_bits(version, mode)
    char_count = len(data.encode('utf-8') if mode == MODE_BYTE else data)
    bits.extend(int_to_bits(char_count, count_bits))
    
    # Data encoding
    if mode == MODE_NUMERIC:
        bits.extend(encode_numeric(data))
    elif mode == MODE_ALPHANUMERIC:
        bits.extend(encode_alphanumeric(data))
    else:  # MODE_BYTE
        bits.extend(encode_byte(data))
    
    return bits


#==============================================================================
# BCH CODE FOR FORMAT INFORMATION
#==============================================================================

# BCH generator polynomial: x^10 + x^8 + x^5 + x^4 + x^2 + x + 1
BCH_GENERATOR = 0b10100110111

# Format mask pattern
FORMAT_MASK = 0b101010000010010

# Error correction level bits
EC_LEVEL_BITS = {
    'L': 0b01,
    'M': 0b00,
    'Q': 0b11,
    'H': 0b10
}


def bch_encode(data_5bits: int) -> int:
    """
    Encode 5 data bits using (15,5) BCH code.
    
    Args:
        data_5bits: 5-bit integer (EC level 2 bits + mask pattern 3 bits)
    
    Returns:
        15-bit encoded format information (before final XOR)
    """
    dividend = data_5bits << 10
    remainder = dividend
    for i in range(14, 9, -1):
        if remainder & (1 << i):
            remainder ^= BCH_GENERATOR << (i - 10)
    return (data_5bits << 10) | remainder


def get_format_string(ec_level: str, mask_pattern: int) -> int:
    """Generate the complete 15-bit format string."""
    data_5bits = (EC_LEVEL_BITS[ec_level] << 3) | mask_pattern
    encoded = bch_encode(data_5bits)
    return encoded ^ FORMAT_MASK


def format_bits_to_list(format_int: int) -> List[int]:
    """Convert 15-bit integer to list of bits."""
    return [(format_int >> (14 - i)) & 1 for i in range(15)]


#==============================================================================
# QR CODE MATRIX CONSTRUCTION
#==============================================================================

class QRMatrix:
    """
    QR Code matrix construction and manipulation.
    
    References:
    - https://www.thonky.com/qr-code-tutorial/module-placement-matrix
    """
    
    # Alignment pattern positions for versions 2-40
    ALIGNMENT_POSITIONS = {
        1: [],
        2: [6, 18],
        3: [6, 22],
        4: [6, 26],
        5: [6, 30],
        6: [6, 34],
        7: [6, 22, 38],
        8: [6, 24, 42],
        9: [6, 26, 46],
        10: [6, 28, 50],
    }
    
    def __init__(self, version: int):
        self.version = version
        self.size = 4 * version + 17
        
        # Matrix values: None=unassigned, 0=white, 1=black
        self.matrix = [[None] * self.size for _ in range(self.size)]
        
        # Track which modules are function patterns
        self.is_function = [[False] * self.size for _ in range(self.size)]
        
        self._place_function_patterns()
    
    def _place_function_patterns(self):
        """Place all function patterns."""
        self._place_finder_patterns()
        self._place_separators()
        self._place_timing_patterns()
        self._place_alignment_patterns()
        self._place_dark_module()
        self._reserve_format_area()
        if self.version >= 7:
            self._reserve_version_area()
    
    def _place_finder_patterns(self):
        """Place the three finder patterns."""
        positions = [
            (0, 0),                          # Top-left
            (self.size - 7, 0),              # Top-right
            (0, self.size - 7)               # Bottom-left
        ]
        for (x, y) in positions:
            self._place_finder_pattern(x, y)
    
    def _place_finder_pattern(self, x: int, y: int):
        """Place a single finder pattern at position (x, y)."""
        for dy in range(7):
            for dx in range(7):
                if (dy == 0 or dy == 6 or dx == 0 or dx == 6 or
                    (2 <= dx <= 4 and 2 <= dy <= 4)):
                    value = 1
                else:
                    value = 0
                self.matrix[y + dy][x + dx] = value
                self.is_function[y + dy][x + dx] = True
    
    def _place_separators(self):
        """Place white separators around finder patterns."""
        # Horizontal
        for x in range(8):
            if x < self.size:
                self._set_function(x, 7, 0)
                self._set_function(self.size - 8 + x, 7, 0)
                self._set_function(x, self.size - 8, 0)
        
        # Vertical
        for y in range(8):
            if y < self.size:
                self._set_function(7, y, 0)
                self._set_function(self.size - 8, y, 0)
                self._set_function(7, self.size - 8 + y, 0)
    
    def _place_timing_patterns(self):
        """Place timing patterns (row 6 and column 6)."""
        for i in range(8, self.size - 8):
            value = (i + 1) % 2
            self._set_function(i, 6, value)
            self._set_function(6, i, value)
    
    def _place_alignment_patterns(self):
        """Place alignment patterns for version 2+."""
        if self.version < 2:
            return
        
        positions = self._get_alignment_positions()
        for row in positions:
            for col in positions:
                if self._overlaps_finder(row, col):
                    continue
                self._place_alignment_pattern(col, row)
    
    def _get_alignment_positions(self) -> List[int]:
        """Get alignment pattern center positions."""
        if self.version in self.ALIGNMENT_POSITIONS:
            return self.ALIGNMENT_POSITIONS[self.version]
        
        first = 6
        last = self.size - 7
        num_intervals = (self.version // 7) + 1
        step = (last - first) // num_intervals
        step = ((step + 1) // 2) * 2
        
        positions = [first]
        pos = last
        while pos > first + step:
            positions.insert(1, pos)
            pos -= step
        positions.append(last)
        return positions
    
    def _overlaps_finder(self, row: int, col: int) -> bool:
        """Check if alignment pattern would overlap finder patterns."""
        if row <= 8 and col <= 8:
            return True
        if row <= 8 and col >= self.size - 9:
            return True
        if row >= self.size - 9 and col <= 8:
            return True
        return False
    
    def _place_alignment_pattern(self, x: int, y: int):
        """Place a single alignment pattern centered at (x, y)."""
        for dy in range(-2, 3):
            for dx in range(-2, 3):
                if abs(dy) == 2 or abs(dx) == 2 or (dy == 0 and dx == 0):
                    value = 1
                else:
                    value = 0
                self._set_function(x + dx, y + dy, value)
    
    def _place_dark_module(self):
        """Place the dark module."""
        x, y = 8, 4 * self.version + 9
        self._set_function(x, y, 1)
    
    def _reserve_format_area(self):
        """Reserve space for format information."""
        for i in range(9):
            self.is_function[8][i] = True
            self.is_function[i][8] = True
        
        for i in range(8):
            self.is_function[8][self.size - 1 - i] = True
            self.is_function[self.size - 1 - i][8] = True
    
    def _reserve_version_area(self):
        """Reserve space for version information (version 7+)."""
        for i in range(6):
            for j in range(3):
                self.is_function[i][self.size - 11 + j] = True
                self.is_function[self.size - 11 + j][i] = True
    
    def _set_function(self, x: int, y: int, value: int):
        """Set a function pattern module."""
        if 0 <= x < self.size and 0 <= y < self.size:
            self.matrix[y][x] = value
            self.is_function[y][x] = True
    
    def place_data(self, data_bits: List[int]) -> int:
        """Place data bits in zigzag pattern."""
        bit_index = 0
        x = self.size - 1
        upward = True
        
        while x >= 0:
            if x == 6:
                x -= 1
            
            y_range = range(self.size - 1, -1, -1) if upward else range(self.size)
            for y in y_range:
                for dx in [0, -1]:
                    col = x + dx
                    if col < 0:
                        continue
                    if self.is_function[y][col]:
                        continue
                    
                    if bit_index < len(data_bits):
                        self.matrix[y][col] = data_bits[bit_index]
                        bit_index += 1
                    else:
                        self.matrix[y][col] = 0
            
            x -= 2
            upward = not upward
        
        return bit_index
    
    def to_string(self, border: int = 4) -> str:
        """Convert matrix to string with quiet zone border."""
        lines = []
        
        # Top border
        for _ in range(border):
            lines.append("  " * (self.size + 2 * border))
        
        for row in self.matrix:
            line = "  " * border  # Left border
            for cell in row:
                if cell == 1:
                    line += "██"
                else:
                    line += "  "
            line += "  " * border  # Right border
            lines.append(line)
        
        # Bottom border
        for _ in range(border):
            lines.append("  " * (self.size + 2 * border))
        
        return "\n".join(lines)


#==============================================================================
# DATA MASKING
#==============================================================================

MASK_PATTERNS: List[Callable[[int, int], bool]] = [
    lambda r, c: (r + c) % 2 == 0,
    lambda r, c: r % 2 == 0,
    lambda r, c: c % 3 == 0,
    lambda r, c: (r + c) % 3 == 0,
    lambda r, c: (r // 2 + c // 3) % 2 == 0,
    lambda r, c: (r * c) % 2 + (r * c) % 3 == 0,
    lambda r, c: ((r * c) % 2 + (r * c) % 3) % 2 == 0,
    lambda r, c: ((r + c) % 2 + (r * c) % 3) % 2 == 0,
]


def apply_mask(matrix: List[List[int]], is_function: List[List[bool]], 
               mask_num: int) -> List[List[int]]:
    """Apply mask pattern to data modules only."""
    size = len(matrix)
    result = [[0 if c is None else c for c in row] for row in matrix]  # Replace None with 0
    mask_func = MASK_PATTERNS[mask_num]
    
    for r in range(size):
        for c in range(size):
            if not is_function[r][c] and mask_func(r, c):
                result[r][c] ^= 1
    
    return result


def calculate_penalty(matrix: List[List[int]]) -> int:
    """Calculate total penalty score for a masked matrix."""
    size = len(matrix)
    penalty = 0
    penalty += _penalty_runs(matrix, size)
    penalty += _penalty_boxes(matrix, size)
    penalty += _penalty_finder_like(matrix, size)
    penalty += _penalty_balance(matrix, size)
    return penalty


def _penalty_runs(matrix: List[List[int]], size: int) -> int:
    """Penalty for runs of 5+ same-color modules."""
    penalty = 0
    
    for r in range(size):
        run_length = 1
        for c in range(1, size):
            curr = matrix[r][c] if matrix[r][c] is not None else 0
            prev = matrix[r][c-1] if matrix[r][c-1] is not None else 0
            if curr == prev:
                run_length += 1
            else:
                if run_length >= 5:
                    penalty += 3 + (run_length - 5)
                run_length = 1
        if run_length >= 5:
            penalty += 3 + (run_length - 5)
    
    for c in range(size):
        run_length = 1
        for r in range(1, size):
            curr = matrix[r][c] if matrix[r][c] is not None else 0
            prev = matrix[r-1][c] if matrix[r-1][c] is not None else 0
            if curr == prev:
                run_length += 1
            else:
                if run_length >= 5:
                    penalty += 3 + (run_length - 5)
                run_length = 1
        if run_length >= 5:
            penalty += 3 + (run_length - 5)
    
    return penalty


def _penalty_boxes(matrix: List[List[int]], size: int) -> int:
    """Penalty for 2x2 same-color boxes."""
    penalty = 0
    for r in range(size - 1):
        for c in range(size - 1):
            color = matrix[r][c] if matrix[r][c] is not None else 0
            c1 = matrix[r][c+1] if matrix[r][c+1] is not None else 0
            c2 = matrix[r+1][c] if matrix[r+1][c] is not None else 0
            c3 = matrix[r+1][c+1] if matrix[r+1][c+1] is not None else 0
            if c1 == color and c2 == color and c3 == color:
                penalty += 3
    return penalty


def _penalty_finder_like(matrix: List[List[int]], size: int) -> int:
    """Penalty for patterns similar to finder patterns."""
    penalty = 0
    pattern1 = [1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0]
    pattern2 = [0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1]
    
    for r in range(size):
        for c in range(size - 10):
            row_pattern = [matrix[r][c+i] if matrix[r][c+i] is not None else 0 for i in range(11)]
            if row_pattern in [pattern1, pattern2]:
                penalty += 40
    
    for c in range(size):
        for r in range(size - 10):
            col_pattern = [matrix[r+i][c] if matrix[r+i][c] is not None else 0 for i in range(11)]
            if col_pattern in [pattern1, pattern2]:
                penalty += 40
    
    return penalty


def _penalty_balance(matrix: List[List[int]], size: int) -> int:
    """Penalty based on dark/light module ratio."""
    dark_count = sum(sum(1 if c == 1 else 0 for c in row) for row in matrix)
    total = size * size
    percent = (dark_count * 100) // total
    
    prev_multiple = percent - (percent % 5)
    next_multiple = prev_multiple + 5
    
    penalty = min(
        abs(prev_multiple - 50) // 5,
        abs(next_multiple - 50) // 5
    ) * 10
    
    return penalty


def choose_best_mask(matrix: List[List[int]], 
                     is_function: List[List[bool]]) -> Tuple[int, int]:
    """Choose the mask pattern with lowest penalty."""
    best_mask = 0
    best_penalty = float('inf')
    
    for mask_num in range(8):
        masked = apply_mask(matrix, is_function, mask_num)
        penalty = calculate_penalty(masked)
        
        if penalty < best_penalty:
            best_penalty = penalty
            best_mask = mask_num
    
    return best_mask, best_penalty


#==============================================================================
# COMPLETE QR CODE GENERATOR
#==============================================================================

# Data capacity tables (simplified for versions 1-10)
EC_CODEWORDS = {
    1: {'L': 7, 'M': 10, 'Q': 13, 'H': 17},
    2: {'L': 10, 'M': 16, 'Q': 22, 'H': 28},
    3: {'L': 15, 'M': 26, 'Q': 36, 'H': 44},
    4: {'L': 20, 'M': 36, 'Q': 52, 'H': 64},
    5: {'L': 26, 'M': 48, 'Q': 72, 'H': 88},
    6: {'L': 36, 'M': 64, 'Q': 96, 'H': 112},
    7: {'L': 40, 'M': 72, 'Q': 108, 'H': 130},
    8: {'L': 48, 'M': 88, 'Q': 132, 'H': 156},
    9: {'L': 60, 'M': 110, 'Q': 160, 'H': 192},
    10: {'L': 72, 'M': 130, 'Q': 192, 'H': 224},
}

DATA_CAPACITY = {
    1: {'L': 19, 'M': 16, 'Q': 13, 'H': 9},
    2: {'L': 34, 'M': 28, 'Q': 22, 'H': 16},
    3: {'L': 55, 'M': 44, 'Q': 34, 'H': 26},
    4: {'L': 80, 'M': 64, 'Q': 48, 'H': 36},
    5: {'L': 108, 'M': 86, 'Q': 62, 'H': 46},
    6: {'L': 136, 'M': 108, 'Q': 76, 'H': 60},
    7: {'L': 156, 'M': 124, 'Q': 88, 'H': 66},
    8: {'L': 194, 'M': 154, 'Q': 110, 'H': 86},
    9: {'L': 232, 'M': 182, 'Q': 132, 'H': 100},
    10: {'L': 274, 'M': 216, 'Q': 154, 'H': 122},
}


class QRCodeGenerator:
    """Complete QR code generator."""
    
    def __init__(self, ec_level: str = 'M'):
        """
        Initialize generator with error correction level.
        
        Args:
            ec_level: 'L' (7%), 'M' (15%), 'Q' (25%), or 'H' (30%)
        """
        if ec_level not in ['L', 'M', 'Q', 'H']:
            raise ValueError(f"Invalid error correction level: {ec_level}")
        self.ec_level = ec_level
    
    def generate(self, data: str, version: int = None) -> List[List[int]]:
        """
        Generate a QR code for the given data.
        
        Args:
            data: String to encode
            version: QR version (1-40), or None to auto-detect
        
        Returns:
            2D list of 0s and 1s representing the QR code
        """
        # Step 1: Determine version
        if version is None:
            version = self._determine_version(data)
        
        if version > 10:
            raise ValueError("This simplified implementation supports versions 1-10")
        
        print(f"Generating Version {version} QR Code with EC Level {self.ec_level}")
        
        # Step 2: Encode data
        data_bits = encode_data(data, version)
        
        # Step 3: Add terminator and padding
        data_codewords = self._pad_data(data_bits, version)
        print(f"Data codewords ({len(data_codewords)}): {data_codewords[:10]}...")
        
        # Step 4: Generate error correction
        num_ec = EC_CODEWORDS[version][self.ec_level]
        ec_codewords = rs_encoder.encode(data_codewords, num_ec)
        print(f"EC codewords ({len(ec_codewords)}): {ec_codewords}")
        
        # Step 5: Combine message
        final_message = data_codewords + ec_codewords
        
        # Step 6: Convert to bits
        final_bits = []
        for byte in final_message:
            final_bits.extend(int_to_bits(byte, 8))
        
        # Step 7: Create matrix
        qr = QRMatrix(version)
        
        # Step 8: Place data
        bits_placed = qr.place_data(final_bits)
        print(f"Placed {bits_placed} bits in {qr.size}x{qr.size} matrix")
        
        # Step 9: Apply best mask
        best_mask, penalty = choose_best_mask(qr.matrix, qr.is_function)
        final_matrix = apply_mask(qr.matrix, qr.is_function, best_mask)
        print(f"Applied mask pattern {best_mask} (penalty: {penalty})")
        
        # Step 10: Add format information
        self._add_format_info(final_matrix, qr.is_function, best_mask, qr.size)
        
        return final_matrix
    
    def _determine_version(self, data: str) -> int:
        """Determine minimum version for the data."""
        mode = detect_mode(data)
        
        # Calculate data length based on mode
        if mode == MODE_BYTE:
            data_len = len(data.encode('utf-8'))
        else:
            data_len = len(data)
        
        for version in range(1, 11):
            capacity = DATA_CAPACITY.get(version, {}).get(self.ec_level, 0)
            if capacity >= data_len:
                return version
        
        raise ValueError("Data too long for versions 1-10")
    
    def _pad_data(self, data_bits: List[int], version: int) -> List[int]:
        """Pad data to required length."""
        bits = list(data_bits)
        capacity_bytes = DATA_CAPACITY[version][self.ec_level]
        capacity_bits = capacity_bytes * 8
        
        # Add terminator (up to 4 bits)
        term_bits = min(4, capacity_bits - len(bits))
        bits.extend([0] * term_bits)
        
        # Pad to byte boundary
        while len(bits) % 8 != 0:
            bits.append(0)
        
        # Convert to bytes
        codewords = bits_to_bytes(bits)
        
        # Add pad codewords (alternating 236, 17)
        pad_bytes = [236, 17]
        i = 0
        while len(codewords) < capacity_bytes:
            codewords.append(pad_bytes[i % 2])
            i += 1
        
        return codewords
    
    def _add_format_info(self, matrix: List[List[int]], 
                         is_function: List[List[bool]], 
                         mask: int, size: int):
        """Add format information to the matrix."""
        format_bits = format_bits_to_list(get_format_string(self.ec_level, mask))
        
        # Primary position: around top-left finder
        # Bits 0-5: row 8, columns 0-5
        # Bit 6: row 8, column 7
        # Bit 7: row 8, column 8
        # Bits 8-14: column 8, rows 7 down to 0 (skipping 6)
        
        for i in range(6):
            matrix[8][i] = format_bits[i]
        matrix[8][7] = format_bits[6]
        matrix[8][8] = format_bits[7]
        matrix[7][8] = format_bits[8]
        for i in range(6):
            matrix[5 - i][8] = format_bits[9 + i]
        
        # Secondary position
        # Bits 0-7: column 8, rows (size-1) down to (size-8)
        # Bits 8-14: row 8, columns (size-8) to (size-1)
        
        for i in range(7):
            matrix[size - 1 - i][8] = format_bits[i]
        for i in range(8):
            matrix[8][size - 8 + i] = format_bits[7 + i]


def matrix_to_image(matrix: List[List[int]], filename: str, 
                    scale: int = 10, border: int = 4):
    """
    Save QR code matrix as PNG image.
    
    Requires PIL/Pillow: pip install Pillow
    """
    try:
        from PIL import Image
    except ImportError:
        print("PIL not available. Install with: pip install Pillow")
        return
    
    size = len(matrix)
    img_size = (size + 2 * border) * scale
    
    img = Image.new('1', (img_size, img_size), 1)  # White background
    pixels = img.load()
    
    for y in range(size):
        for x in range(size):
            if matrix[y][x] == 1:
                # Fill scaled pixel area
                for dy in range(scale):
                    for dx in range(scale):
                        px = (border + x) * scale + dx
                        py = (border + y) * scale + dy
                        pixels[px, py] = 0  # Black
    
    img.save(filename)
    print(f"Saved QR code to {filename}")


#==============================================================================
# DEMONSTRATION
#==============================================================================

def demo():
    """Run demonstrations of all components."""
    print("=" * 70)
    print("QR CODE GENERATOR - COMPLETE DEMONSTRATION")
    print("=" * 70)
    
    # Demo 1: GF(256) Arithmetic
    print("\n" + "=" * 70)
    print("1. GALOIS FIELD GF(256) ARITHMETIC")
    print("=" * 70)
    
    a, b = 83, 202
    print(f"\na = {a}, b = {b}")
    print(f"a + b (XOR) = {gf.add(a, b)}")
    print(f"a * b = {gf.multiply(a, b)}")
    print(f"a / b = {gf.divide(a, b)}")
    print(f"a^10 = {gf.power(a, 10)}")
    print(f"a^(-1) = {gf.inverse(a)}")
    print(f"Verify: a * a^(-1) = {gf.multiply(a, gf.inverse(a))}")
    
    # Demo 2: Reed-Solomon
    print("\n" + "=" * 70)
    print("2. REED-SOLOMON ERROR CORRECTION")
    print("=" * 70)
    
    # Example data (Version 1-M has 16 data codewords, 10 EC)
    data = [16, 32, 12, 86, 97, 128, 236, 17, 236, 17, 236, 17, 236, 17, 236, 17]
    num_ec = 10
    
    print(f"\nData ({len(data)} bytes): {data}")
    gen = rs_encoder.build_generator(num_ec)
    print(f"Generator polynomial (degree {gen.degree})")
    ec = rs_encoder.encode(data, num_ec)
    print(f"EC codewords ({len(ec)}): {ec}")
    
    # Demo 3: Data Encoding
    print("\n" + "=" * 70)
    print("3. DATA ENCODING")
    print("=" * 70)
    
    test_data = "HELLO WORLD"
    print(f"\nEncoding: '{test_data}'")
    print(f"Detected mode: {'ALPHANUMERIC' if detect_mode(test_data) == MODE_ALPHANUMERIC else 'OTHER'}")
    bits = encode_data(test_data, 1, MODE_ALPHANUMERIC)
    print(f"Encoded bits ({len(bits)}): {''.join(map(str, bits[:40]))}...")
    
    # Demo 4: Complete QR Code Generation
    print("\n" + "=" * 70)
    print("4. COMPLETE QR CODE GENERATION")
    print("=" * 70)
    
    qr_gen = QRCodeGenerator(ec_level='M')
    
    # Generate a QR code
    test_message = "HELLO"
    print(f"\nGenerating QR code for: '{test_message}'")
    
    try:
        qr_matrix = qr_gen.generate(test_message, version=1)
        
        # Display the QR code
        print("\n" + "=" * 70)
        print("GENERATED QR CODE:")
        print("=" * 70)
        
        qr_obj = QRMatrix(1)  # Create fresh matrix for display helper
        qr_obj.matrix = qr_matrix
        print(qr_obj.to_string(border=2))
        
        # Try to save as image
        try:
            matrix_to_image(qr_matrix, "qr_hello.png")
        except Exception as e:
            print(f"(Could not save image: {e})")
        
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
    
    print("\n" + "=" * 70)
    print("DEMONSTRATION COMPLETE")
    print("=" * 70)


if __name__ == "__main__":
    demo()
