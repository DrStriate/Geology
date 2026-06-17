from pathlib import Path
import sys

# Ensure the project root is on sys.path when running app.py from a subfolder.
def setRoot():
    ROOT = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(ROOT))