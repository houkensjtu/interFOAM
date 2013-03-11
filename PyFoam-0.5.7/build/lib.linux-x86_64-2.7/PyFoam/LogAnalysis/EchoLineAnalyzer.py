#  ICE Revision: $Id$ 
"""Echos a line"""

from LogLineAnalyzer import LogLineAnalyzer

class EchoLineAnalyzer(LogLineAnalyzer):
    """Test implementation. Simply echos every line it gets"""
    
    def __init__(self):
        LogLineAnalyzer.__init__(self)

    def doAnalysis(self,line):
        print "<"+self.parent.getTime()+">"+line+"<"
        
