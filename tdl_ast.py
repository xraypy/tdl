#!/usr/bin/python 

import sys
from astlib.shell import shell
shell(debug=True, scripts = sys.argv[1:]).cmdloop()

