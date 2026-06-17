# The following code stub will make the 'run Pyrthon file from a subfolder' work,
#  but it is not recommended to use this in production code.
'''
# Ensure the project root is on sys.path when running app.py from a subfolder.
from pathlib import Path
import sys
ROOT_DIR = Path(__file__).resolve().parent.parent
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))
#end of the code stub
'''

from config_root import setRoot
setRoot()

from src.utils import custom_helper

def main():
    result = custom_helper()
    print(result)

if __name__ == "__main__":
    main()
