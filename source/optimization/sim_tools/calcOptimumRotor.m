function  [mu,c_opt,phi,beta] = calcOptimumRotor(tsr,Cl,N,r,R,aoa_design,approach,varargin)
% CALCOPTIMUMROTOR Calculate chord and inflow angle (or blade twist) for
%      optimum rotor design, according to various theories.
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
% [mu,c_opt,phi,beta] =
%           calcOptimumRotor(tsr,Cl,N,r,R,aoa_design,approach,hubRad,ad)
%
%  mu = r/R;
%  c_opt = optimal chord distibution
%  phi = inflow angles associated with optimal chord distribution
%  beta = blade pitch schedule associated with phi
%
%  tsr = design tip speed ratio
%  Cl = design Cl
%  N = number of blades
%  r = [nx1] array of n local rotor radius points (m)
%  R = outer rotor radius (m)
%  aoa_design = design angle of attack (degrees); used for computation of blade twist
%     distribution once inflow angle, phi, is known.
%  hubRad = (required only for approach = 4).  hub radius (m)
%  ad = (required only for approach = 5).  AeroDyn input file data structure as read by readFastAD()
%  approach = integer; choose 1 thru 5:
%     1  Manwell, Betz optimum with no wake rotation.
%     2  Manwell, Betz optimum with wake rotation.
%     3  Burton, optimum rotor.
%     4  Burton, optimum with hub and tip losses.
%     5  Burton, linear chord taper; requires determination of design
%        angles of attack for specified Cl distribution using actual airfoil
%        data.
%
% References:
%   Manwell, J. F., McGowan, J. G., and Rogers, A. L. Wind Energy
%       Explained: Theory, Design and Application: John Wiley & Sons,
%       Ltd., 2008.
%   Burton, T., Sharpe, D., Jenkins, N., and Bossanyi, E. Wind Energy
%       Handbook. Chichester: John Wiley and Sons, Ltd., 2008.*
%
% * Note: Numerous errors in formulas from the Burton edition were found in
%         the course of coding this script.  Expect discrepancies between
%         the book formulas and formulas in this script.
%

if approach==4  % use hub radius only for tip/hub loss computations
    if isempty(varargin)
        error('hubRadius needs to be specified for approach=4');
        return;
    else
        hubRad=varargin{1};
        mur=hubRad/R;
    end
end

if approach==5  % use hub radius only for tip/hub loss computations
    if isempty(varargin)
        error('AeroDyn input data structure, ad, needs to be specified for approach=5');
        return;
    else
        ad=varargin{end};
    end
end

mu=r/R;
tsrr=tsr*mu;

disp(' ')
disp(' ')

switch approach
    
    case 1
        disp('Betz optimum with no wake rotation, using Manwell')
        
        phi=atand(2/3./tsrr);   % Eq.(3.6.5) Manwell
        c_opt=2*4*pi*r/3./tsrr/Cl/N.*sind(phi);  % Eq.(3.6.6) Manwell
        beta=phi-aoa_design;
        
    case 2 
        disp('Betz optimum with wake rotation, using Manwell')
        
        phi=2/3*atand(1./tsrr);   % Eq.(3.6.5) Manwell
        c_opt=2*4*pi*r/Cl/N.*(1-cosd(phi));  % Eq.(3.6.6) Manwell
        beta=phi-aoa_design;
        
    case 3 
        disp('Optimum rotor, using Burton')
        
        num=8/9;
        sqroot=(1-1/3)^2+tsr^2*mu.^2.*(1+2/9/tsr^2./mu.^2).^2;
        den=tsr*Cl*sqrt(sqroot);
        sigmar=num./den;  % Eq.(3.67a) Burton
        c_opt=sigmar.*(2*pi*R/N);
        num=1-1/3;
        den=tsr.*mu.*(1+2/3^2/tsr^2./mu.^2);
        phi=atand(num./den);  % Eq.(3.68a) Burton
        beta=phi-aoa_design;
        
    case 4 
        disp('Optimum with hub and tip losses, using Burton')
        
        % calculate tip and hub loss factors
        a=1/3;
        sqroot=1+(tsr*mu).^2/(1-a)^2;
        exponent1=(N/2)*(mu-1)./mu.*sqrt(sqroot);
        ft=2/pi*acos(exp(exponent1));  % Eq.(3.76) Burton, tip loss factor
        exponent2=(-N/2)*(mu-mur)./mu.*sqrt(sqroot);
        fr=2/pi*acos(exp(exponent2));  % Eq.(3.78) Burton, root loss factor
        f=ft.*fr;  % Eq.(3.79) Burton
        %     f=f./f;  % this line can be used to check that when losses are 1, you get the same curve as in approach=3
        %calculate chord and twist distributions
        num=4*a*(1-a);
        sqroot=(1-a./f).^2+(tsr*mu.*(1+a*(1-a./f)/tsr^2./mu.^2./f)).^2;
        den=tsr*Cl.*sqrt(sqroot);
        sigmar=num./den;  % Eq.(3.85) Burton
        c_opt=sigmar.*(2*pi*R/N);
        num=1-a./f;
        den=tsr.*mu.*(1+(a*(1-a./f))/tsr^2./mu.^2./f);
        phi=atand(num./den);  % Eq.(3.86) Burton
        beta=phi-aoa_design;
        
    case 5 
        disp('Linear chord taper, using Burton')
        
%         sp=0.70;
        sp=0.80;
        c_opt=(8/9/tsr/sp)*(2-tsr*mu/tsr/sp)*(2*pi/Cl/tsr/N)*R;  % Eq.(3.69), Burton
        numer=8/9;
        sqroot=(1-1/3)^2 + tsr^2.*mu.^2.*(1+2./9./tsr^2./mu.^2).^2;
        denom=(N.*c_opt.*tsr./2./pi).*sqrt(sqroot);
        Cl_opt=numer./denom*R;  % Following Eq.(3.69), Burton
        
        figure(500)
        plot(mu,Cl_opt)
        title({'C_l distribution',...
            sprintf('Design C_l=%4.2f, design TSR=%4.2f',Cl,tsr)})
        xlabel('r/R')
        ylabel('C_l')
        legend('Linear taper blade optimum')
        grid on
        
        % find the AoA's that enable the optimal Cl values for linear taper
        for i=1:ad.NumFoil
            fn=ad.FoilNm{i};
            fn=strrep(fn,'"','');
            af=readAirfoilData(fn);
            pointer1=find(af.AoA>-25);
            pointer2=find(af.AoA<30);
            pointer=intersect(pointer1,pointer2);
            LoverD=af.CL./af.CD;
            LoverDmax=max(LoverD);
            pointer_LoverDmax=find(LoverD==LoverDmax);
            pointer_LoverDmax=pointer_LoverDmax(1);
            Cldesign=af.CL(pointer_LoverDmax);
            Cddesign=af.CD(pointer_LoverDmax);
            AoAdesign=af.AoA(pointer_LoverDmax);
            
            figure(501)
            plot(af.AoA(pointer),LoverD(pointer),'k-o',af.AoA(pointer_LoverDmax),LoverD(pointer_LoverDmax),'r-v')
            xlabel('Angle of attack, deg')
            ylabel('C_l/C_D')
            
            Clmin=-1.5;
            Clmax=2;
            
            figure(502)
            subplot(1,2,1)
            plot(af.CD(pointer),af.CL(pointer),'k-o',...
                [0 Cddesign],[0 Cldesign],'r-v')
            axis([0 .05 Clmin Clmax])
            ylabel('C_l')
            xlabel('C_d')
            title(sprintf('Max C_l/C_d=%5.2f',LoverDmax))
            grid on
            legend('C_l/C_d curve','Max C_l/C_d','Location','East')
            subplot(1,2,2)
            plot(af.AoA(pointer),af.CL(pointer),'k-o',...
                AoAdesign,Cldesign,'r-v',...
                [0 20],[Cl_opt(i) Cl_opt(i)],'r')
            axis([-20 20 Clmin Clmax])
            xlabel('Angle of attack, deg')
            ylabel('C_l')
            title([sprintf('Cl_m_a_x_C_l_/_C_d=%5.2f',Cldesign),', ',...
                sprintf('AoA_m_a_x_C_l_/_C_d=%5.2f',AoAdesign)])
            grid on
            legend(sprintf('C_l data for %s',fn),'C-l for max C-l/C_d','Desired C_l for optimum linear taper blade','Location','Southeast')
            
            disp('Use the cursor to choose the design alpha for the plotted airfoil')
            [AoA(i),~] = ginput(1);
            AoA(i)
        end
        
        for i=1:ad.BldNodes
            AoA_BldNodes(i)=AoA(ad.NFoil(i));
        end
        
        % compute inflow angles
        tanphi=(1-1/3) ./ (tsr.*mu.*(1+2./3./tsr^2./mu.^2));
        phi=atand(tanphi);   % Eq.(3.68a), Burton
        beta=phi-AoA_BldNodes';
        
end

end
