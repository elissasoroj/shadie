#!/usr/bin/env python

"""
Random utilities
"""

import sys
from loguru import logger


LOGFORMAT = (
    "🌿 {level} | "
    "<level>{file}</level> | "
    "<black>{message}</black>"
)

def colorize():
    """colorize the logger if stderr is IPython/Jupyter or a terminal (TTY)"""
    try:
        import IPython        
        tty1 = bool(IPython.get_ipython())
    except ImportError:
        tty1 = False
    tty2 = sys.stderr.isatty()
    if tty1 or tty2:
        return True
    return False

LOGGERS = [0]
def set_log_level(log_level="INFO"):
    """Set the loglevel for loguru logger. 

    This removes default loguru handler, but leaves any others, 
    and adds a new one that will filter to only print logs from 
    shadie modules, which should use `logger.bind(name='shadie')`.
    """
    for idx in LOGGERS:
        try:
            logger.remove(idx)
        except ValueError:
            pass
    idx = logger.add(
        sink=sys.stderr,
        level=log_level,
        colorize=colorize(),
        format=LOGFORMAT,
        filter=lambda x: x['extra'].get("name") == "shadie",
    )
    LOGGERS.append(idx)
    logger.enable("shadie")
