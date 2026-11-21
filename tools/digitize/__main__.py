"""Entry point: ``uv run --group digitize python -I -m tools.digitize``."""

from tools.digitize.cli import main

raise SystemExit(main())
