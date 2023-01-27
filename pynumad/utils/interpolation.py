from scipy import interpolate


def interpolator_wrap(x,v,xq,method = 'linear',axis = 0, extrapolation = None):
    """This function is designed to emulate the arg structure and output
    of matlabs interp1d function.
    """
    if method == 'linear':
        interpolator = interpolate.interp1d(x,v,'linear', axis,bounds_error=False, fill_value= 'extrapolate')
        vq = interpolator(xq)
        return vq
    if method == 'pchip':
        interpolator = interpolate.PchipInterpolator(x,v,axis,extrapolate=True)
        vq = interpolator(xq)
        return vq
    if method == 'spline':
        interpolator = interpolate.interp1d(x,v,'cubic', axis,bounds_error=False,fill_value='extrapolate')
        vq = interpolator(xq)
        return vq
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