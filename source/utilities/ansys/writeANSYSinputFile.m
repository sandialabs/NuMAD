function writeANSYSinputFile(fid,mat,ansysSecNumber,coreMatName)
    %%%%%Find the face sheet%%%%
    cellMat={};
    for i=1:length(mat.layer)
        cellMat = [cellMat; {mat.layer(i).layerName}];   %Create a cell array to use "find"
    end
    kbalsa=find(strcmp(coreMatName,cellMat)==1);
    iLayer=1:(kbalsa-1); %Number of distinct materials in the top face

    % Find the number of layers in the face
    qty=0; %Quantity of layers counter
    for i=1:numel(iLayer)
        qty=qty+mat.layer(iLayer(i)).quantity;
    end

    %Loop through the top facesheet layers

    fprintf(fid, '!*************** ansysSecNumber = %i ***************\n',ansysSecNumber);
    fprintf(fid, '/POST1\n');
    fprintf(fid, '*DEL,iel\n');
    fprintf(fid, '*DEL,enum\n');
    fprintf(fid, '*DEL,nelTemp\n');
    fprintf(fid, 'RSYS, SOLU\n');
    fprintf(fid, 'ALLSEL\n');
    fprintf(fid, 'ESEL, S, SEC,,%i\n',ansysSecNumber);
    fprintf(fid, '*GET, enum, ELEM, 0, NUM, MIN, !  lowest element number in the selected set\n');
    fprintf(fid, '*get, nelTemp, ELEM,0,count\n');
    fprintf(fid, '*DIM, iel,ARRAY,nelTemp\n');
    fname=strcat('section-',int2str(ansysSecNumber),'-faceAvgStresses'); %Text file to store averages stresses to be computed by ansys

    if isfile(strcat(fname,'.txt'))
        delete(strcat(fname,'.txt')) %Check if file exisits and delete it because following ansys commands append to file.
    end
    fprintf(fid, '*CFOPEN, %s, txt,,APPEND\n',fname);

    %Create an array with the element numbers in the selected set
    fprintf(fid, '*DO, J, 1,nelTemp  !Loop through elements\n');
    fprintf(fid, 'iel(J)=enum\n');
    fprintf(fid, 'enum =ELNEXT(enum)  !Next higher element number above N in selected set\n');
    fprintf(fid, '*ENDDO\n');

    fprintf(fid, '\n');
    fprintf(fid, 'ALLSEL\n');
    fprintf(fid, '*DO, J, 1,nelTemp  !Loop through elements\n');
    fprintf(fid, '	S11a=0 !Initialize average stress variables for each element\n');
    fprintf(fid, '	S22a=0\n');
    fprintf(fid, '	S33a=0\n');
    fprintf(fid, '	S23a=0\n');
    fprintf(fid, '	S13a=0\n');
    fprintf(fid, '	S12a=0\n');
    fprintf(fid, '	*DO, I, 1,%i    !Loop through face layers\n',qty);
    fprintf(fid, '	    LAYER,I\n');
    fprintf(fid, '		SHELL,MID   !Stress result at midlayer\n');
    fprintf(fid, '		ESEL,S,ELEM,,iel(J)\n');
    fprintf(fid, '		ETABLE,ERAS !Each element gets a new element table\n');
    fprintf(fid, '	    ETABLE,S11,S,X,AVG !AVG - Store averaged element centroid value\n');
    fprintf(fid, '	    ETABLE,S22,S,Y,AVG\n');
    fprintf(fid, '		ETABLE,S33,S,Z,AVG\n');
    fprintf(fid, '		ETABLE,S23,S,YZ,AVG\n');
    fprintf(fid, '		ETABLE,S13,S,XZ,AVG\n');
    fprintf(fid, '        ETABLE,S12,S,XY,AVG\n');
    fprintf(fid, '	    *GET,tempS11, ELEM, iel(J), ETAB, S11\n');
    fprintf(fid, '		*GET,tempS22, ELEM, iel(J), ETAB, S22\n');
    fprintf(fid, '		*GET,tempS33, ELEM, iel(J), ETAB, S33\n');
    fprintf(fid, '		*GET,tempS23, ELEM, iel(J), ETAB, S23\n');
    fprintf(fid, '		*GET,tempS13, ELEM, iel(J), ETAB, S13\n');
    fprintf(fid, '		*GET,tempS12, ELEM, iel(J), ETAB, S12\n');
    fprintf(fid, '		S11a=S11a+tempS11\n');
    fprintf(fid, '		S22a=S22a+tempS22\n');
    fprintf(fid, '		S33a=S33a+tempS33\n');
    fprintf(fid, '		S23a=S23a+tempS23\n');
    fprintf(fid, '		S13a=S13a+tempS13\n');
    fprintf(fid, '		S12a=S12a+tempS12\n');
    fprintf(fid, '	*ENDDO\n');
    fprintf(fid, '	S11a=S11a/%i\n',qty);
    fprintf(fid, '	S22a=S22a/%i\n',qty);
    fprintf(fid, '	S33a=S33a/%i\n',qty);
    fprintf(fid, '	S23a=S23a/%i\n',qty);
    fprintf(fid, '	S13a=S13a/%i\n',qty);
    fprintf(fid, '	S12a=S12a/%i\n',qty);
    fprintf(fid, '	ELNO=iel(J)  !It is needed to refer to ELNO in the command below\n');
    fprintf(fid, '*VWRITE,ELNO,S11a,S22a,S33a,S23a,S13a,S12a\n');
    fprintf(fid, '(E20.12,E20.12,E20.12,E20.12,E20.12,E20.12,E20.12)\n');
    fprintf(fid, '*ENDDO\n');
    fprintf(fid, '*CFCLOS\n');
    fprintf(fid, '\n');
    fprintf(fid, '\n');
end