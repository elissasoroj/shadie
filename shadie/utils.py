#!/usr/bin/env python

"""
Random utilities
"""

import sys
from typing import Union, List
from loguru import logger


# logger will show time, function name, and message.
LOGFORMAT = (
    "{time:hh:mm} | {level: <7} | "
    "<b><red>{function: <15}</red></b> | "
    "<level>{message}</level>"
)


# colorize the logger if stdout is IPython/Jupyter or a terminal (TTY)
try:
    import IPython
    TTY1 = bool(IPython.get_ipython())
except ImportError:
    TTY1 = False
TTY2 = sys.stdout.isatty()


def set_loglevel(loglevel="INFO"):
    """
    Set the loglevel for loguru logger. Using 'enable' here as 
    described in the loguru docs for logging inside of a library.
    This sets the level at which logger calls will be displayed 
    throughout the rest of the code.
    """
    config = {}
    config["handlers"] = [{
        "sink": sys.stdout,
        "format": LOGFORMAT,
        "level": loglevel,
        "colorize": TTY1 or TTY2,
    }]
    logger.configure(**config)
    logger.enable("shadie")
