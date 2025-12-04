# QR Code Generator Tutorial

A complete implementation of QR code generation from scratch in Python, following the ISO/IEC 18004 specification. This project includes both the implementation and a comprehensive LaTeX tutorial document.

## Contents

- **`qr_code_generator.py`** - Complete QR code generator implementation
- **`qr_code_tutorial.tex`** - LaTeX source for the tutorial document
- **`qr_code_tutorial.pdf`** - Compiled tutorial document
- **`qr_hello.png`** - Example generated QR code

## Features

The implementation includes:

- Galois Field GF(256) arithmetic
- Polynomial operations over GF(256)
- Reed-Solomon error correction encoding
- Data encoding (numeric, alphanumeric, byte modes)
- BCH code for format information
- QR code matrix construction with function patterns
- Data masking with penalty calculation
- Complete QR code generation

## Requirements

- Python 3.x
- Pillow (optional, for image generation): `pip install Pillow`

## Usage

Run the demonstration:

```bash
python qr_code_generator.py
```

This will generate a QR code for "HELLO" and save it as `qr_hello.png`.

## Example

```python
from qr_code_generator import QRCodeGenerator, matrix_to_image

# Create a QR code generator
qr_gen = QRCodeGenerator(ec_level='M')

# Generate a QR code
matrix = qr_gen.generate("Hello, World!")

# Save as image (requires Pillow)
matrix_to_image(matrix, "output.png")
```

## Error Correction Levels

- **L** - ~7% error correction
- **M** - ~15% error correction (default)
- **Q** - ~25% error correction
- **H** - ~30% error correction

## License

This is an educational implementation based on ISO/IEC 18004, Thonky's QR Code Tutorial, and the Wikiversity Reed-Solomon Guide.

