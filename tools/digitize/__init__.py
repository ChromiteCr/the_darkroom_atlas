"""Digitise curves from manufacturer datasheet PDFs into reference-curve CSV files.

See ``tools/digitize/README.md``. The core (geometry, calibration, path selection, CSV
output) needs only numpy, PyYAML and pydantic; reading PDFs needs pdfminer.six and
rendering needs pypdfium2, both in the ``digitize`` dependency group.
"""
