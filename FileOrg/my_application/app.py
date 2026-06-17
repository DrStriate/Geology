# The following code stub will make the 'run Pyrthon file from a subfolder' work,
# in VS Code but it is not recommended to use this in production code.
from config_root import setRoot
setRoot()

from src.utils import custom_helper

def main():
    result = custom_helper()
    print(result)

if __name__ == "__main__":
    main()
