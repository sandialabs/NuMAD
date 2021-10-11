function loadsFORad2ansys(z,Mi,Ti,theta,prebend,forcefilename)
    % Inputs
    % z - spanwise location where loads (M and T ar defined). Must include root
    % Mi - the input moment at corresponding to the spanwise location, z.
         % Pure edgewise moment when theta = 0. Must include root moment
    % Ti - input blade twisting moments at the spanwise locations
    % theta - Angle projection
    % forcefilename - name of outputfile to be created. ".forces" will be
    %                 appeded to forcefilename. 

    %Transform input resultant moment to the ANSYS coordinate system
    Mx=Mi*sind(theta); 
    My=Mi*cosd(theta);
    
    
    Mx(1)
    My(1)
    diffz=diff(z);

    fx=getForcesFromMoments(diffz,My);
    fy=getForcesFromMoments(diffz,Mx);

    
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % figure(3);
    % plot(z,Mi,'-ok')
    % hold on;
    % plot(z(1:end-1),mCheck,'*k--')
    % legend('Mi','mCheck')
    % xlabel('Spanwise coordinate, z [m]')
    % ylabel('Moment, M(z) [N m]')
    % hold on;
    % 
    
%     figure(4)
%     plot(z,fx,'ok-')
%     hold on;
%    % plot(z,fy,'xb-')
%     hold on;
            
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%

    % Output the results to a file that ad2ansys will read
    fid=fopen(strcat(forcefilename,'.forces'),'wt');
    fprintf(fid,'\tZ\t\t\t\t Fx\t\t\t\t Fy\t\t\t\t M\t\t\t\t Alpha\t\t\t x_off\t\t\t y_off\n');

    for kk=2:length(z) %skip first point (the root)
        fprintf(fid,'%.8E\t%.8E\t%.8E\t%.8E\t%.8E\t%.8E\t%.8E\n',z(kk), fx(kk),fy(kk),Ti(kk),0,0,prebend(kk-1));
    end

    fclose(fid);
end
% function f=relocateForces(f,z,zNew)
%     for i=2:numel(z) %skip first point (the root)
%         Mroot=z(i)*f(i); %Moment at root due to f(i)
%         f(i)=Mroot/zNew(i);
%     end
% 
% end
% 
% function z=reoderZ(z,L)
%     dz=zeros(size(z)); %Vector to hold the distances that each point needs to shift
%     IBedge=z(1);       %Assuming first poit is zero
%     halfDZ = z(2);  % initialize section width
%     for bk = 2:numel(z)-1% Skip the first point (root) and the last point
%     
% 
% %         if bk>max
% %            disp('')
% % 
% %         end         
%             
%         OBedge = z(bk) + halfDZ;  % outboard edge of section
%         IBedge=z(bk)-halfDZ;
%         alpha=z(bk+1) - OBedge;
%         if (bk~=numel(z)) &&  alpha<=0
%             if abs(alpha)<0.01
%                 OBedge = z(bk)+(z(bk+1)-z(bk))/2;
%             else
%                 OBedge = z(bk+1)+alpha;  % Reassign outboard edge of section
%             end
%             
%             %OBedge = 2*z(bk+1)- OBedge;  % outboard edge of section
%             dz(bk)=IBedge+(OBedge-IBedge)/2-z(bk);  %for updating force such that moment equivalency is obationed
%             z(bk)=IBedge+(OBedge-IBedge)/2;
%             
%             %Print (Can be deleteed)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             fprintf('\nStation: %i: IBedge = %f  Z=%f  OBedge=%f',bk, IBedge,z(bk),OBedge)
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             
% %             if abs(z(bk+1)-L)<0.0
%             halfDZ = z(bk+1) - OBedge;
%         else
%             %Print (Can be deleteed)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             fprintf('\nStation: %i: IBedge = %f  Z=%f  OBedge=%f',bk, IBedge,z(bk),OBedge)
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             
%             halfDZ = z(bk+1) - OBedge;
%         end
% 
%         
%     end
%     bk=numel(z);
%         
%     if z(bk)+alpha ~= L
%         IBedge=z(bk)-alpha;
%         halfDZ=(L-IBedge)/2;
%         dz(bk)=IBedge+halfDZ-z(bk);
%         z(bk)=IBedge+halfDZ;
%         OBedge=L;
%     end
%     fprintf('\nStation: %i: IBedge = %f  Z=%f  OBedge=%f',bk, IBedge,z(bk),OBedge)
% 
% end

    function ql=getForcesFromMoments(diffz,M)
        ct=0;
        ql=zeros(size(M)); %Forces to be output
        for i=length(diffz):-1:1
            ql(i+1)=M(i)/diffz(i);
            dz=diffz(i);
            for k=1:ct
                dz=dz+diffz(i+k);
                ql(i+1)=ql(i+1)-1/diffz(i)*ql(i+1+k)*dz;
            end
            ct=ct+1;
        end
        
        %Moment Check
        ct=0;
        mCheck=zeros(size(diffz));
        for i=length(diffz):-1:1
             dz=diffz(i);
            mCheck(i)=ql(i+1)*dz;
            for k=1:ct
                dz=dz+diffz(i+k);
                mCheck(i)=mCheck(i)+ql(i+1+k)*dz;
            end
            ct=ct+1;
        end
    end
