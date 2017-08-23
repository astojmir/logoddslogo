#!/usr/bin/env python

import cgi
import cgitb; cgitb.enable()
from logoddslogolib._cgi import main

if __name__ == "__main__":

    main(__file__, proteins=False)
