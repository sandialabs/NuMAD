if 0  % set up parameters manually
    flutterInput.fstFile='Test13.fst';
    flutterInput.bladeFile='Baseline_Blade.dat';
%    flutterInput.bladeFile='FASTBlade_precomp.dat';
    flutterInput.aeroFile='Test13_AD.ipt';
    flutterInput.outFile='name.out';
    flutterInput.LCS=ones(4,1)*2*pi;  % one value for each airfoil in the AD file
    flutterInput.pitchAxisDomain=[0    0.7000    1.6000    2.5000    3.4000    4.3000    5.2000    6.1000    7.0000    8.7500   10.5000   12.2500   14.0000   15.7500   17.5000   19.2500   21.0000   22.7500   24.5000   26.2500   28.0000   29.7500   31.5000   33.2500];  % spanwise locations at which pitch axis is defined by pitchAxisVal's
    flutterInput.pitchAxisVal=[0.5000    0.5000    0.4768    0.4536    0.4304    0.4071    0.3839    0.3607    0.3375    0.3375    0.3375    0.3375    0.3375    0.3375    0.3375    0.3375    0.3375    0.3375    0.3375    0.3375    0.3375    0.3375    0.3375    0.3375];  % pitch axis values at locations defined by pitchAxisDomain
    flutterInput.numadBladeLen=33.25;
    flutterInput.OmegaArray=(0:2.5:50);  % array of rotor speeds (RPM) for flutter analysis
else % use flutter input parameters defined and saved by NuMAD
    load flutterInput
end

tol=0.001;

if 1  % peform classical flutter analysis
    [freq,damp]=feaAutoAllModes(flutterInput,tol,'F');
else  % visualize mode shapes at a given speed and excitation frequency
    omega=37.8;
    freq=6.3;
    fea(flutterInput,omega,freq,1,'F');
end

