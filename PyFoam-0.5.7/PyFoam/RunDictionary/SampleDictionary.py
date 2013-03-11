#  ICE Revision: $Id:$
"""Working with a directory of samples"""

from os import path,listdir

class SampleDirectory(object):
    """A directory of sampled times"""

    def __init__(self,case,dirName="samples"):
        """@param case: The case directory
        @param dirName: Name of the directory with the samples"""

        self.dir=path.join(case,dirName)
        self.times=[]
        
    
