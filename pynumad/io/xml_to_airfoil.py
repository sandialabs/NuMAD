import re
import numpy as np

def xml_to_airfoil(airfoil, filecontents = None): 
    """TODO docstring
    
    Parameters
    ----------

    Returns
    -------
    """
    # DEVNOTE: this works for the kinds of af files
    #          found in the BAR reference blades
    #          but not for general xmls
    #          this should be extended to be more general at some point
    # [coords, reference] = readAirfoilXML(filecontents)

    # The following regular expression pattern matches any number of characters
    # found between the opening and closing "reference" tags
    fulltext = ''.join(filecontents)
    
    pattern = '<reference>(.*)</reference>'
    t = re.search(pattern,fulltext)
    reference = t.group(1)

    for line in filecontents:
        #check if there is a tag
        if re.search('<',line):
            continue
        # otherwise, assume coordinate data
        else:
            if re.search('\t',line):
                line = line.replace('\t', ' ')
            x, y = line.split(' ')
            x = float(x)
            y = float(y)
            try:
                coords = np.append(coords, [[x,y]],axis=0)
            except UnboundLocalError:
                coords = np.array([[x,y]])
    airfoil.reference = reference
    airfoil.coordinates = coords