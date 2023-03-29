########################################################################
#                    Part of the SNL NuMAD Toolbox                     #
#  Developed by Sandia National Laboratories Wind Energy Technologies  #
#              See license.txt for disclaimer information              #
########################################################################

class Stack:
    """A class definition for a stack of composite layers.

    Parameters
    ----------

    Attributes
    ----------
    name : string
        Name of the stack or composite material used by NuMAD, e.g. '000000_HP_LE_PANEL'
    indices : list
        Indices of stack, ``[in board station, out board station, 
        1st kepoint, 2nd keypoint]``, e.g. ``[ibSta,obSta,keypt1,keypt2]``
    plygroups : list
        List of ``ply`` dataclasses
        
    Example
    -------
        ``stack = StackDef();``
    
    See also ``xlsBlade``, ``BladeDef``, ``BladeDef.updateBOM``
    """

    name: str = None
    indices = None
    plygroups: list = []

    def addply(self,ply):
        """This method adds a Ply object to stack

        Parameters
        ----------
        ply : Ply object
            Ply object to be added
        Returns
        -------
        None

        Example
        -------
            ``stack.addply(ply)``
        """
        if self.plygroups:
            if ((ply.component == self.plygroups[-1].component) and
                    (ply.angle == self.plygroups[-1].angle)):
                self.plygroups[-1].nPlies += 1
            else:
                self.plygroups.append(ply)
        else:
            self.plygroups = []
            self.plygroups.append(ply)

            
