import numpy as np

from scipy.interpolate import interp1d, PchipInterpolator,\
    PPoly, CubicSpline


def interpolator_wrap(x,v,xq,method = 'linear',axis = 0, extrapolation = None):
    """This function is designed to emulate the arg structure and output
    of matlabs interp1d function.
    """
    if method == 'linear':
        interpolator = interp1d(x,v,'linear', axis,bounds_error=False, fill_value= 'extrapolate')
        vq = interpolator(xq)
        return vq
    if method == 'pchip':
        interpolator = PchipInterpolator(x,v,axis,extrapolate=True)
        vq = interpolator(xq)
        return vq
    if method == 'spline':
        interpolator = interp1d(x,v,'cubic', axis,bounds_error=False,fill_value='extrapolate')
        vq = interpolator(xq)
        return vq
    if method == 'pp':
        pass
    if method == 'v5cubic':
        raise Exception("Method error for interpolator_wrap. 'v5cubic' not implemented")
    if method == 'makima':
        raise Exception("Method error for interpolator_wrap. 'makima' not implemented")
    if method == 'nearest':
        raise Exception("Method error for interpolator_wrap. 'nearest' not implemented")
    if method == 'next':
        raise Exception("Method error for interpolator_wrap. 'next' not implemented")
    if method == 'previous':
        raise Exception("Method error for interpolator_wrap. 'previous' not implemented")

    
def calcGenLinePP(blade_struct: dict):
    """Calculate blade reference line piecewise polynomials
    blade_struct = calcGenLinePP(blade_struct) updates the piecewise
    polynomial representation of the blade's Presweep and Precurve
    reference lines. This function is called by NuMAD_genline.
    
    The fields PresweepRef and PrecurveRef are required in blade_struct.
    Each of these fields has the following data structure:
        method: 'normal' | 'shear'
                This field is not used by calcGenLinePP.
        table: N-by-3 matrix with columns span,offset,slope
                This table provides the offset and slope constraints of the
                reference line at specific spanwise locations along the
                blade. NaN may be used wherever a constraint is not
                desired.
        pptype: 'poly' | 'spline' | 'pchip' | 'linear' | 'disabled'
                This field selects the interpolation method to use to
                create the piecewise polynomial
                poly = minimum order polynomial which satisfies all constraints
                spline = cubic spline (offset constraints only)
                pchip = shape-preserving cubic spline (offset constraints only)
                linear = linear interpolation (offset constraints only)
                disabled = returns straight line
            pp: piecewise polynomial data created by this function
            dpp: piecewise polynomial data of reference line's derivative
    
       See also NuMAD_genline, PPoly, interp1.
    """
    
    # PresweepRef
    spline_type = blade_struct["PresweepRef"]["pptype"]
    PresweepRef = blade_struct["PresweepRef"]["table"]
    if spline_type in ['linear','spline','pchip']:
        if PresweepRef.shape[0] > 1:
            pp = CubicSpline(PresweepRef[:,0],PresweepRef[:,1])
        else:
            pp = CubicSpline([0,1], [0,0])
    
    blade_struct["PresweepRef"]["pp"] = pp
    # dc = np.diag(np.arange(pp.order - 1,1+- 1,- 1),1)
    
    # blade_struct["PresweepRef"]["dpp"] = PPoly(pp.breaks,pp.coefs * dc)
    
    # PrecurveRef
    spline_type = blade_struct["PrecurveRef"]["pptype"]
    PrecurveRef = blade_struct["PrecurveRef"]["table"]
    if spline_type in ['linear','spline','pchip']:
        if PrecurveRef.shape[0] > 1:
            pp = CubicSpline(PrecurveRef[0,:],PrecurveRef[1,:])
        else:
            pp = CubicSpline([0,1], [0,0])
    
    blade_struct["PrecurveRef"]["pp"] = pp
    # dc = np.diag(np.arange(pp.order - 1,1+- 1,- 1),1)
    
    # blade_struct["PrecurveRef"]["dpp"] = PPoly(pp.breaks,pp.coefs * dc)
    
    return blade_struct