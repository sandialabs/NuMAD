function blade_struct = calcGenLinePP(blade_struct)
%CALCGENLINEPP  Calculate blade reference line piecewise polynomials
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   blade_struct = calcGenLinePP(blade_struct) updates the piecewise
%   polynomial representation of the blade's Presweep and Precurve 
%   reference lines. This function is called by NuMAD_genline.
%   
%   The fields PresweepRef and PrecurveRef are required in blade_struct.
%   Each of these fields has the following data structure:
%       method: 'normal' | 'shear'
%               This field is not used by calcGenLinePP.
%        table: N-by-3 matrix with columns span,offset,slope
%               This table provides the offset and slope constraints of the
%               reference line at specific spanwise locations along the
%               blade. NaN may be used wherever a constraint is not
%               desired.
%       pptype: 'poly' | 'spline' | 'pchip' | 'linear' | 'disabled'
%               This field selects the interpolation method to use to
%               create the piecewise polynomial
%               poly = minimum order polynomial which satisfies all constraints
%               spline = cubic spline (offset constraints only)
%               pchip = shape-preserving cubic spline (offset constraints only)
%               linear = linear interpolation (offset constraints only)
%               disabled = returns straight line
%           pp: piecewise polynomial data created by this function
%          dpp: piecewise polynomial data of reference line's derivative
%
%   See also NuMAD_genline, mkpp, interp1.

    spline_type = blade_struct.PresweepRef.pptype;
    PresweepRef = blade_struct.PresweepRef.table;
    switch spline_type
        case {'linear','spline','pchip'}
            if size(PresweepRef,1) > 1
                pp = interp1(PresweepRef(:,1),PresweepRef(:,2),spline_type,'pp');
            else
                pp = mkpp(PresweepRef(1,1)*[1 1], [0 PresweepRef(1,2)]);
            end
        case 'poly'
            keypts = PresweepRef;
            Nkp = size(keypts,1);
            indkp = find(isnan(keypts(:,2:3))==0)+Nkp;
            K = length(indkp);  % order of poly
            A = zeros(K);
            B = zeros(K,1);
            kp.x0 = min(keypts(:,1));
            for krow = 1:K
                kp.i = indkp(krow);
                [kp.r kp.c] = ind2sub(size(keypts),kp.i);
                kp.x = keypts(kp.r,1);
                kp.y = keypts(kp.r,2);
                kp.m = keypts(kp.r,3);
                if kp.c == 2  % position constraint
                    for kcol = 1:K
                        % note that mkpp() requires the x-coord to be
                        % relative to the start of each break
                        A(krow,kcol) = (kp.x-kp.x0).^(K-kcol);
                    end
                    B(krow,1) = kp.y;
                else % slope constraint
                    for kcol = 1:K
                        if isnan(kp.m)
                            error('logical error - poly slope');
                        end
                        if kcol < K
                            A(krow,kcol) = (K-kcol)*(kp.x-kp.x0).^(K-kcol-1);
                        else
                            A(krow,kcol) = 0;
                        end
                    end
                    B(krow,1) = kp.m;
                end
            end
            coefs = A \ B;
            pp = mkpp([kp.x0 max(keypts(:,1))], coefs');
        case 'disabled'
            pp = mkpp([0 0], [0 0]);
    end
    blade_struct.PresweepRef.pp = pp;
    dc = diag(pp.order-1:-1:1,1);  % derivative coefficient matrix
    blade_struct.PresweepRef.dpp = mkpp(pp.breaks,pp.coefs*dc); % create spline with derivative coefficents
    
    spline_type = blade_struct.PrecurveRef.pptype;
    PrecurveRef = blade_struct.PrecurveRef.table;
    switch spline_type
        case {'linear','spline','pchip'}
            if size(PrecurveRef,1) > 1
                pp = interp1(PrecurveRef(:,1),PrecurveRef(:,2),spline_type,'pp');
            else
                pp = mkpp(PrecurveRef(1,1)*[1 1], [0 PrecurveRef(1,2)]);
            end
        case 'poly'
            keypts = PrecurveRef;
            Nkp = size(keypts,1);
            indkp = find(isnan(keypts(:,2:3))==0)+Nkp;
            K = length(indkp);  % order of poly
            A = zeros(K);
            B = zeros(K,1);
            kp.x0 = min(keypts(:,1));
            for krow = 1:K
                kp.i = indkp(krow);
                [kp.r kp.c] = ind2sub(size(keypts),kp.i);
                kp.x = keypts(kp.r,1);
                kp.y = keypts(kp.r,2);
                kp.m = keypts(kp.r,3);
                if kp.c == 2  % position constraint
                    for kcol = 1:K
                        % note that mkpp() requires the x-coord to be
                        % relative to the start of each break
                        A(krow,kcol) = (kp.x-kp.x0).^(K-kcol);
                    end
                    B(krow,1) = kp.y;
                else % slope constraint
                    for kcol = 1:K
                        if isnan(kp.m)
                            error('logical error - poly slope');
                        end
                        if kcol < K
                            A(krow,kcol) = (K-kcol)*(kp.x-kp.x0).^(K-kcol-1);
                        else
                            A(krow,kcol) = 0;
                        end
                    end
                    B(krow,1) = kp.m;
                end
            end
            coefs = A \ B;
            pp = mkpp([kp.x0 max(keypts(:,1))], coefs');
        case 'disabled'
            pp = mkpp([0 0], [0 0]);
    end
    blade_struct.PrecurveRef.pp = pp;
    dc = diag(pp.order-1:-1:1,1);  % derivative coefficient matrix
    blade_struct.PrecurveRef.dpp = mkpp(pp.breaks,pp.coefs*dc); % create spline with derivative coefficents

end