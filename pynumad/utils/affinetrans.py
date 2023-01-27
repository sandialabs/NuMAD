from scipy.spatial.transform import Rotation
import numpy as np

def rotation(axis, angle):
    """
    Designed to replace matlab's makehgtform
    """
    r = Rotation.from_euler(axis, angle)
    rmatrix = np.eye(4)
    rmatrix[0:3,0:3] = r.as_matrix()
    return rmatrix

def translation(xtrans,ytrans,ztrans):
    tmatrix = np.eye(4)
    tmatrix[0:3,3] = [xtrans,ytrans,ztrans]
    return tmatrix


if __name__ == "__main__":
    print(translation(2,3,4))
    print(rotation('z', 1))