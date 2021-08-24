'''
========================================================================

    File Name:      debug.py
    Author:         Jason Li
    Description:    Print tabbed output information.
    Usage:          prt(msg, tab=0)
                        Print tabbed output message.
                    prtList(lst, filename="test.dat")
                        Print the content of a list to GNUPlot-able \
                        data file.
                    timeStamp()
                        Print time stamp.
                    
========================================================================
'''

import datetime
from const import *

curTab = 0

def prt(msg, tab=0):
    """
    Print tabbed output message.
    Input:  msg --Message.
            tab --1 (Increase tab from next message), 0 (Keep current 
                  tab), -1 (Reduce tab from this message). (Default: 0)
    """
    global curTab
    if (tab>0):
        tabs = "\t"*curTab
        curTab += 1
    elif (tab<0):
        curTab -= 1
        tabs = "\t"*curTab
    else:
        tabs = "\t"*curTab
    print(tabs+msg)

def prtList(lst, filename="test.dat"):
    """
    Print the content of a list to GNUPlot-able data file.
    Input:  lst                 --List to print.
            filename            --Filename. (Default: test.dat)
    Output: Scratch/[filename]  --Output file.
    """
    with open(scratchPath+filename, "w") as f:
        [f.write(str(var)+"\n") for var in lst]
        

def timeStamp():
    """Print time stamp."""
    print(datetime.datetime.now())
