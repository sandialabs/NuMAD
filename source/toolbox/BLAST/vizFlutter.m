function vizFlutter(inputfile,res,modeNum,sf,aviname)
%vizFlutter Creates animation of mode shape
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   vizFlutter(inputfile,res,modeNum,sf,aviname)
%                    
%   This function accepts an input file containing blade information and a
%   results file to create an animation of the mode shapes of a blade
%   projected onto a three-dimensional wireframe through the deformations
%   of beam theory.
%
%      input:
%      inputfile    = input filename for BLAST analysis job
%      res          = results filename from BLAST analysis
%      modeNum      = mode number to visualize
%      sf           = scale factor for mode shape deformations
%      aviname      = name of avi file to be generated
 
%      output:
%      None: There is no explicit output for this function. The specified
%      avi file is generated in the job directory.



% Create movie file with required parameters and set up window
  fps= 10;
  mov = VideoWriter(aviname);
  mov.FrameRate = fps;
  mov.Quality = 75;
  open(mov);
  
  close all;
    
  fig1 = figure;
  set(fig1, 'Units', 'normalized', 'Position', [0.25 0.25 0.75 0.75]);
  
  %read input files
   [bladefn,adfn,bladeLength] = readInputFile(inputfile);
  
   %create blade wireframe
   [B1,B2,B3l,B3u,ealoc2,ealoc3] = createBladeWireFrame(bladefn,adfn,bladeLength,0,0,0,[0 0 0],0);
   
   % read in phase (phase1) and out of phase (phase2) modeshapes
   [numStations,~]=size(B1);
   numModesInOutputFile = 10;
   [phase1list,phase2list] = readResults(res,numModesInOutputFile,numStations);
   
   %select modeshape of interest
   phase1 = phase1list(:,:,modeNum);
   phase2 = phase2list(:,:,modeNum);
   
   %process deformation of blade for in and out of phase shapes
   [defxu,defxl,defyu,defyl,defz] = processPointDeformation(B1,B2,B3u,B3l,ealoc2,ealoc3,phase1);
   [defxu2,defxl2,defyu2,defyl2,defz2] = processPointDeformation(B1,B2,B3u,B3l,ealoc2,ealoc3,phase2);
   
   %temporary deformation for getting axes limits
   sfac=1; cfac=1;
   B1defu = B1 + defxu.*sfac.*sf + defxu2.*cfac.*sf;
   B1defl = B1 + defxl.*sfac.*sf + defxl2.*cfac.*sf;
   B2defu = B2 + defyu.*sfac.*sf + defyu2.*cfac.*sf;
   B2defl = B2 + defyl.*sfac.*sf + defyl2.*cfac.*sf;
   B3defu = B3u + defz.*sfac.*sf + defz2.*cfac.*sf; 
   B3defl = B3l + defz.*sfac.*sf + defz2.*cfac.*sf; 
   
   [axdata] = calculateAxesBounds(B1defl,B1defu,B2defl,B2defu,B3defl,B3defu,1.5);
   
   %set deltaOmegat
   delOmegat = pi/10;
   omegat=0;
   for i=1:100
        %time stepping calculating cos and sin of omegat
        cfac=cos(omegat);
        sfac=sin(omegat);
        omegat=omegat+delOmegat;
   
        %calculate deformed blade shape at instance in time (sf=scale
        %factor)
        B1defu = B1 + defxu.*sfac.*sf + defxu2.*cfac.*sf;
        B1defl = B1 + defxl.*sfac.*sf + defxl2.*cfac.*sf;
        B2defu = B2 + defyu.*sfac.*sf + defyu2.*cfac.*sf;
        B2defl = B2 + defyl.*sfac.*sf + defyl2.*cfac.*sf;
        B3defu = B3u + defz.*sfac.*sf + defz2.*cfac.*sf; 
        B3defl = B3l + defz.*sfac.*sf + defz2.*cfac.*sf; 
      
        %plot undeformed mesh and surface of deformed blade.
        mesh(B1,B2,B3u,'EdgeColor','black');
        hold on;
        mesh(B1,B2,B3l,'EdgeColor','black');
        alpha(0.5);
        surf(B1defu,B2defu,B3defu,'FaceColor','cyan','EdgeColor','none');
        surf(B1defl,B2defl,B3defl,'FaceColor','cyan','EdgeColor','none');
        hold off;
        
        %set axis and lighting commands
        axis equal;
        axis(axdata);
        grid on;
        view(3);
        set(gcf,'Color',[1  1 1]);
        camlight right; lighting phong
    
        % put this plot in a movieframe
        F = getframe(gcf);
        writeVideo(mov,F);

   end

    %close movie file
    close(mov);
   
end

function [bladefn,adfn,bladeLength] = readInputFile(inputfile)

%   This function reads the main  blast input file
%
%   input:
%   inputfile   = BLAST input file
%
%   output:
%   bladefn     = blade data file
%   adfn        = aerodyn data file
%   bladeLength = length of blade

    fid = fopen(inputfile,'r');
    fstfn = fscanf(fid,'%s',1);
    bladefn = fscanf(fid,'%s',1);
    adfn = fscanf(fid,'%s',1);
    fclose(fid);
    
    fast=readFastMain(fstfn);
    bladeLength=fast.TurbConf.TipRad;
    
end

function [defxu,defxl,defyu,defyl,defz] = processPointDeformation(B1,B2,B3u,B3l,ealoc2,ealoc3,phase)

%This function  calculates the deformation of a point in the body
%
%    input:
%    B1       = 1 coordinates of body
%    B2       = 2 coordiantes of body
%    B3l      = 3 coordinates of lower body
%    B3u      = 3 coordinates of upper body
%    ealoc2   = 2 coordinates of elastic axis
%    ealoc3   = 3 doordinates of elastic axis
%    phase    = mode-shape of deformation

%    output:
%    defxu    = deformation in x of upper body
%    defxl    = deformation in x of lower body
%    defyu    = deformation in y of upper body
%    defyl    = deformation in y of lower body
%    defz     = deformation in z of body

    eax = phase(:,1);
    eay = phase(:,2);
    eaz = phase(:,3);
    eatx = phase(:,4);
    eaty = phase(:,5); %CHECK THESE FOR SIGN....
    eatz = phase(:,6);
    
    [l,w]=size(B1);
    
    for i=1:l
        for j=1:w
       
            defxu(i,j) = eax(i) - eaty(i)*(B2(i,j)-ealoc2(i)) + eatz(i)*(B3u(i,j) - ealoc3(i));
            defxl(i,j) = eax(i) - eaty(i)*((B2(i,j)-ealoc2(i))  + eatz(i)*(B3l(i,j) - ealoc3(i))); 
            defyu(i,j) = eay(i) - eatx(i)*(B3u(i,j) - ealoc3(i));
            defyl(i,j) = eay(i) - eatx(i)*(B3l(i,j) - ealoc3(i));
            defz(i,j) = eaz(i) + eatx(i)*(B2(i,j)-ealoc2(i)) ;
            
        end
    end
end

function [axdata] = calculateAxesBounds(B1l,B1u,B2l,B2u,B3l,B3u,axfac)

%This function calculates the axes bounds for creating an animation
%
%    input:
%    B1l      = 1 coordinates of lower body
%    B1u      = 1 coordinates of upper body
%    B2l      = 2 coordiantes of lower body
%    B2u      = 2 coordinates of upper body
%    B3l      = 3 coordinates of lower body
%    B3u      = 3 coordinates of upper body
%    axfac    = factor to extend axes by for padding

%    output:
%    axdata   = three dimensional axes limits

    %get bounds on axes
    maxx1 = max(max(abs(B1l)));
    maxx2 = max(max(abs(B1u)));

    maxy1 = max(max(abs(B2l)));
    maxy2 = max(max(abs(B2u)));

    maxz1 = max(max(abs(B3l)));
    maxz2 = max(max(abs(B3u)));

    if(maxx1 < maxx2)
        maxx=maxx2;
    else
        maxx=maxx1;
    end

    if(maxy1 < maxy2)
        maxy=maxy2;
    else
        maxy=maxy1;
    end

    if(maxz1 < maxz2)
        maxz=maxz2;
    else
        maxz=maxz1;
    end

    axdata = [0 maxx -maxy*axfac maxy*axfac -maxz*axfac maxz*axfac];

end
