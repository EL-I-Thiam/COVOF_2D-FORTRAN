#!/usr/bin/env python3
import sys
import pypdf

if len(sys.argv) < 2:
    print("Usage: check_pages.py <pdf_file> [min_pages]")
    sys.exit(1)

pdf_file = sys.argv[1]
min_pages = int(sys.argv[2]) if len(sys.argv) > 2 else 40

try:
    reader = pypdf.PdfReader(pdf_file)
    pages = len(reader.pages)
    status = "OK" if pages >= min_pages else "ATTENTION (Moins de {} pages)".format(min_pages)
    print(f"{pdf_file}: {pages} pages [{status}]")
    if pages < min_pages:
        sys.exit(2)
except Exception as e:
    print(f"Erreur de lecture de {pdf_file}: {e}")
    sys.exit(1)
