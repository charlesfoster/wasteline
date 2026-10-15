import sys
from pathlib import Path

BIN = Path(__file__).resolve().parents[2] / "bin"
if str(BIN) not in sys.path:
    sys.path.insert(0, str(BIN))
