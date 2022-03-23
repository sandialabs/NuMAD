function designvar = ReadANSYS_DeflectionResults(blade, config, iLoad, deflectionFilename)
    nSpan=length(blade.ispan);
    data=zeros(nSpan,6);  %u1Avg,u2Avg,u3Avg,0,theta2,theta3
    for iSpan=1:nSpan
        fileName=[deflectionFilename '-' int2str(iSpan) '.out'];
        temp=txt2mat(fileName);  %node number, Displacements of the nodes on the crossection
        delete(fileName);

        %Displacement
        for k=1:3
            data(iSpan, k) =  mean(temp(:,k+4));
        end

        [nNode,~]=size(temp);  %Number of nodes
        [xmax,LE]=max(temp(:,2)); %Find the max x location (column 2) 
        [xmin,TE]=min(temp(:,2)); %Find the max x location (column 2) 
        [ymax,LP]=max(temp(:,3)); %Find the max y location (column 3) 
        [ymin,HP]=min(temp(:,3)); %Find the max y location (column 3) 
        %close all;
        %plot(temp(:,2),temp(:,3),'ok')
        %hold on;
        P = temp(LE,2:4); %x,y,z coordinates for point P at leading edge
        Q = temp(TE,2:4); %x,y,z coordinates for point Q at trailing edge
        PQ=P-Q;
        %quiver(Q(1),Q(2),PQ(1),PQ(2));
        %plot(temp(:,2)+temp(LE,2),temp(:,3)+temp(LE,3),'xb')
        %axis equal;
        R=P+temp(LE,2:4);
        S=Q+temp(TE,2:4);
        RS = R-S;
        %quiver(Q(1),Q(2),RS(1),RS(2));
        %data(iSpan, 5) =  180/pi* acos(dot(RS(1:2:3),PQ(1:2:3))/(vecnorm(RS(1:2:3))*vecnorm(PQ(1:2:3))));
        %data(iSpan, 6) =  180/pi* acos(dot(RS(1:2),PQ(1:2))/(vecnorm(RS(1:2))*vecnorm(PQ(1:2))));
        index = 1:2:3;
        a = RS(index(1))*PQ(index(1));
        b = RS(index(2))*PQ(index(2));
        c = sqrt(PQ(index(1))^2+PQ(index(2))^2);
        d = sqrt(RS(index(1))^2+RS(index(2))^2);
        data(iSpan, 5) =  180/pi* acos((a+b)/(c*d));

        index = 1:2;
        a = RS(index(1))*PQ(index(1));
        b = RS(index(2))*PQ(index(2));
        c = sqrt(PQ(index(1))^2+PQ(index(2))^2);
        d = sqrt(RS(index(1))^2+RS(index(2))^2);
                
        arg=(a+b)/(c*d);
        if  arg>1
            if round(arg,6) ==1
                data(iSpan, 6) =  180/pi* acos(round(arg,6));
            else
                data(iSpan, 6) =  180/pi* acos(arg); %Imaginary results
            end
        end

        T = temp(LP,2:4); %x,y,z coordinates for point T on suction side
        U = temp(HP,2:4); %x,y,z coordinates for point U on pressure side
        TU=T-U;
        
        V=T+temp(LP,2:4);
        W=U+temp(HP,2:4);
        VW = V-W;

        index = [2,3];
        a = VW(index(1))*TU(index(1));
        b = VW(index(2))*TU(index(2));
        c = sqrt(TU(index(1))^2+TU(index(2))^2);
        d = sqrt(VW(index(1))^2+VW(index(2))^2);
        data(iSpan, 4) =  180/pi* acos((a+b)/(c*d));
        %title(['ispan:' int2str(iSpan) ' theta:' num2str(data(iSpan, 6))])
    end
    designvar.deflection{iLoad} = [];      
    for jj=1:6
        designvar.deflection{iLoad}=[designvar.deflection{iLoad} data(:,jj)];
    end 
end
    