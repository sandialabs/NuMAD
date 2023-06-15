import numpy as np
import os
from pynumad.analysis.ansys.utility import txt2mat

def read_1_ANSYSoutput(filename):
    with open(filename, 'r') as fid:
        format = '%s %f %s'
        lines = fid.readlines()
        imax = len(lines)
        data = np.nan
        ct = 0
        
        for i in range(imax):
            tline = lines[i]
            try:
                formatted_line = np.loadtxt(tline,format)
                # will only work if "format" is readable
                # in that line.
                if not len(formatted_line[2])==0:
                    data = formatted_line[2]
                    # w/ a numerical argument is reported.
            finally:
                pass
    
    fid.close()
    return data


def read_nlist(filename): 
    if filename is None:
        filename = 'NLIST.lis'
    
    # Open the file and read the entire contents
    with open(filename, 'rb') as fid:
        # filecontents = np.transpose(fread(fid,inf,'uint8=>char'))
        filecontents = fid.read()
    nlist = []
    tbl_hdrs = regexp(filecontents,'NODE\s*X\s*Y\s*Z\s*')
    for kTbl in range(tbl_hdrs.size-1):
        tbl = filecontents(range(tbl_hdrs(kTbl),tbl_hdrs(kTbl + 1)))
        data = regexprep(tbl,'\s*NODE\s*X\s*Y\s*Z\s*','')
        nlist.appen(np.loadtxt(data,'%f %f %f %f'))
    
    # get the last table
    kTbl = 0
    tbl = filecontents(np.arange(tbl_hdrs(kTbl + 1),end()+1))
    data = regexprep(tbl,'\s*NODE\s*X\s*Y\s*Z\s*','')
    
    nlist.append(np.loadtxt(data,'%f %f %f %f'))
    nlist = np.array(nlist)
    return nlist

def readAnsysDeflections(blade, config, iLoad, deflectionFilename): 
    nSpan = len(blade.ispan)
    data = np.zeros((nSpan,6))
    
    for iSpan in range(nSpan):
        fileName = deflectionFilename+'-'+str(iSpan)+'.out'
        temp_results = txt2mat(fileName)
        os.remove(fileName)
        #Displacement
        for k in range(3):
            data[iSpan,k] = np.mean(temp_results[:,k + 4])
        nNode,__ = temp_results.shape
        xmax = np.amax(temp_results[:,1])
        LE = np.argmax(temp_results[:,1])
        xmin = np.amin(temp_results[:,1])
        TE = np.argmin(temp_results[:,1])
        ymax = np.amax(temp_results[:,2])
        LP = np.argmax(temp_results[:,2])
        ymin = np.amin(temp_results[:,2])
        HP = np.argmin(temp_results[:,2])
        #close all;
        #plot(temp(:,2),temp(:,3),'ok')
        #hold on;
        P = temp_results[LE,1:5]
        Q = temp_results[TE,1:5]
        PQ = P - Q
        #quiver(Q(1),Q(2),PQ(1),PQ(2));
        #plot(temp(:,2)+temp(LE,2),temp(:,3)+temp(LE,3),'xb')
        #axis equal;
        R = P + temp_results[LE,1:5]
        S = Q + temp_results[TE,1:5]
        RS = R - S
        #quiver(Q(1),Q(2),RS(1),RS(2));
        #data(iSpan, 5) =  180/pi* acos(dot(RS(1:2:3),PQ(1:2:3))/(vecnorm(RS(1:2:3))*vecnorm(PQ(1:2:3))));
        #data(iSpan, 6) =  180/pi* acos(dot(RS(1:2),PQ(1:2))/(vecnorm(RS(1:2))*vecnorm(PQ(1:2))));
        index = [0,2]
        a = RS[index[0]] * PQ[index[0]]
        b = RS[index[1]] * PQ[index[1]]
        c = np.sqrt(PQ[index[0]] ** 2 + PQ[index[1]] ** 2)
        d = np.sqrt(RS[index[0]] ** 2 + RS[index[1]] ** 2)
        data[iSpan,5] = 180 / np.pi * np.arccos((a + b) / (c * d))
        index = [0,1]
        a = RS[index[0]] * PQ[index[0]]
        b = RS[index[1]] * PQ[index[1]]
        c = np.sqrt(PQ[index[0]] ** 2 + PQ[index[1]] ** 2)
        d = np.sqrt(RS[index[0]] ** 2 + RS[index[1]] ** 2)
        arg = (a + b) / (c * d)
        if arg > 1:
            if np.round(arg,6) == 1:
                data[iSpan,5] = 180 / np.pi * np.arccos(np.round(arg,6))
            else:
                data[iSpan,5] = 180 / np.pi * np.arccos(arg)
        T = temp_results[LP,1:5]
        U = temp_results[HP,1:5]
        TU = T - U
        V = T + temp_results[LP,1:5]
        W = U + temp_results[HP,1:5]
        VW = V - W
        index = [1,2]
        a = VW[index[0]] * TU[index[0]]
        b = VW[index[1]] * TU[index[1]]
        c = np.sqrt(TU[index[0]] ** 2 + TU[index[1]] ** 2)
        d = np.sqrt(VW[index[0]] ** 2 + VW[index[1]] ** 2)
        data[iSpan,3] = 180 / np.pi * np.arccos((a + b) / (c * d))
        #title(['ispan:' int2str(iSpan) ' theta:' num2str(data(iSpan, 6))])
    
    deflections = []
    for jj in range(6):
        deflections.append(data[:,jj])
    
    return deflections

def readANSYSElementTable(filename, pat, NCOLS): 
    # readANSYSElementTable Read an ANSYS POST1 element table listing.
    # **********************************************************************
    # *                   Part of the SNL NuMAD Toolbox                    *
    # * Developed by Sandia National Laboratories Wind Energy Technologies *
    # *             See license.txt for disclaimer information             *
    # **********************************************************************
    # data = readANSYSElementTable(filename)
    # Read an ANSYS POST1 element table listing.
    # Usage: data = readANSYSElementTable(filename)
    #  where FILENAME is file name string, default 'Strains.txt'
    #        pat - pattern that repeated in the table e.g pat = 'ELEM\s*EPELX\s*EPELY\s*EPELZ\s*EPELXY\s*EPELYZ\s*EPELXZ\s*';
    #        NCOL - number of columns in the data table
    #        DATA is 7-column matrix [e.g.ELEM, EPELX, EPELY, EPELZ, EPELXY, EPELYZ, EPELXZ]
    
    defaultfn = 'Strains.txt'
    # hard-code filename if not specified
    if filename is None:
        filename = defaultfn
        pat = 'ELEM\s*EPELX\s*EPELY\s*EPELZ\s*EPELXY\s*EPELYZ\s*EPELXZ\s*'
        NCOLS = 7
    
    #     # user select filename if not specified
#     if ~exist('filename','var') || isempty(filename)
#         [fn,pn] = uigetfile( ...
#             {'*.txt','Text files(*.txt)'; ...
#             '*.*','All files (*.*)'},...
#             'Select ANSYS element list',defaultfn);
#         if isequal(fn,0) || isequal(pn,0)
#             disp('Operation canceled by user.')
#             return;
#         end
#         filename = fullfile(pn,fn);
#     end
    
    # Open the file and read the entire contents
    with open(filename, 'rb') as fid:
        filecontents = fid.read()
    #assignin('base','filecontents',filecontents);  #debugging
    
    # process the tables
    data = cell(1,NCOLS)
    tbl_hdrs = regexp(filecontents,pat)
    
    tbl_hdrs[end() + 1] = np.asarray(filecontents).size
    
    #assignin('base','tbl_hdrs',tbl_hdrs);  #debugging
    for kTbl in range(tbl_hdrs.size-1):
        tbl = filecontents(np.arange(tbl_hdrs(kTbl),tbl_hdrs(kTbl + 1) - 1+1))
        tbl = regexprep(tbl,pat,'')
        data = np.array([[data],[np.loadtxt(tbl,np.matlib.repmat(' %f',1,NCOLS))]])
    
    data = cell2mat(data)
    return data

def readAnsysFailure(fileName = None): 
    fid = open(fileName)
    terminate = 0
    elemFailure = []
    while (terminate == 0):

        fLine = fgetl(fid)
        if (fLine == - 1):
            terminate = 1
        else:
            fArray = str2num(fLine)
            if (len(fArray) > 2):
                elemFailure = np.array([[elemFailure],[fArray[1]]])
                #elemFailure = elemFailure + 60000*fArray(2)^2;

    
    fid.close()
    #     fid = fopen(fileName);
    #     terminate = 0;
    #     elemFailure = zeros(numEls,1);
    #     k = 1;
    #     while(terminate == 0)
    #         fLine = fgetl(fid);
    #         if(fLine == -1)
    #             terminate = 1;
    #         else
    #             fArray = str2num(fLine);
    #             if(length(fArray) > 2)
    #                 elemFailure(k) = fArray(2);
    #             end
    #         end
    #     end
    # fid.close();
    return elemFailure

def readAnsysFreq(fileName = None): 
    fid = open(fileName)
    fLine = fgetl(fid)
    while (not contains(fLine,'1') ):

        fLine = fgetl(fid)

    
    terminate = 0
    modeNum = 0
    Freq = []
    while (terminate == 0):

        if (fLine == - 1):
            terminate = 1
        else:
            lnLst = str2num(fLine)
            if (len(lnLst) >= 2):
                if (lnLst[1] > 0.0):
                    modeNum = modeNum + 1
                    Freq = np.array([Freq,lnLst[1]])
            fLine = fgetl(fid)

    
    fid.close()
    return Freq

def readAnsysLinearBuckling(blade = None,config = None,iLoad = None,fid = None,bucklingFilename = None): 
    fid = open(np.array([bucklingFilename,'.out']))
    for jj in np.arange(1,5+1).reshape(-1):
        tline = fgetl(fid)
    
    data = cell(1,5)
    while 1:

        tline = fgetl(fid)
        if not ischar(tline) :
            break
        data = np.array([[data],[np.loadtxt(tline,'%f %f %f %f %f')]])

    
    fid.close()
    print(' ')
    data = cell2mat(data)
    linearLoadFactors = data(np.arange(1,config.analysisFlags.globalBuckling+1),2)
    
    os.delete(np.array([bucklingFilename,'.out']))
    return linearLoadFactors

def readANSYSnforce(filename = None): 
    # readANSYSnforce  Read an ANSYS list of Elements.
    # **********************************************************************
    # *                   Part of the SNL NuMAD Toolbox                    *
    # * Developed by Sandia National Laboratories Wind Energy Technologies *
    # *             See license.txt for disclaimer information             *
    # **********************************************************************
    # Read an ANSYS list of Elements.
    # Usage: data = readANSYSnforce(FILENAME)
    #  where FILENAME is file name string, default 'Elements.txt'
    #        DATA is 3-column matrix [ELEM, MAT, SEC]
    
    defaultfn = 'nforce.txt'
    # hard-code filename if not specified
    if not ('filename' is not None) :
        filename = defaultfn
    
    #     # user select filename if not specified
    #     if ~exist('filename','var') || isempty(filename)
    #         [fn,pn] = uigetfile( ...
    #             {'*.txt','Text files(*.txt)'; ...
    #             '*.*','All files (*.*)'},...
    #             'Select ANSYS element list',defaultfn);
    #         if isequal(fn,0) || isequal(pn,0)
    #             disp('Operation canceled by user.')
    #             return;
    #         end
    #         filename = fullfile(pn,fn);
    #     end
    
    # Open the file and read the entire contents
    fid = open(filename)
    if (fid == - 1):
        raise Exception('Could not open file "%s"',filename)
    
    filecontents = np.transpose(fread(fid,inf,'uint8=>char'))
    fid.close()
    #assignin('base','filecontents',filecontents);  #debugging
    
    # process the tables
    NCOLS = 7
    data = cell(1,NCOLS)
    pat = 'NODE\s*FX\s*FY\s*FZ\s*MX\s*MY\s*MZ\s*'
    #pat='*** NOTE ***\s*CP =\s*\d*.\d*\s*TIME=\s*\d*:\d*:\d*\r\s*Use the SLIST command to list section data for element\s*\d*.\s*Section\s*\r\s*data overrides the real constant data.|*** NOTE ***\s*CP =\s*\d*.\d*\s*TIME=\s*\d*:\d*:\d*\r\s*Use the SLIST command to list section data for element\s*\d*.\s*Section data\s*\r\s*overrides the real constant data. ';
    tbl_hdrs = regexp(filecontents,pat)
    
    tbl_hdrs[end() + 1] = np.asarray(filecontents).size
    
    #assignin('base','tbl_hdrs',tbl_hdrs);  #debugging
    for kTbl in np.arange(1,np.asarray(tbl_hdrs).size - 1+1).reshape(-1):
        tbl = filecontents(np.arange(tbl_hdrs(kTbl),tbl_hdrs(kTbl + 1) - 1+1))
        tbl = regexprep(tbl,pat,'')
        data = np.array([[data],[np.loadtxt(tbl,np.matlib.repmat(' %f',1,NCOLS))]])
    
    data = cell2mat(data)
    # only columns 1,2,6 are needed now
    # Note: np.loadtxt can skip fields by using #*f in place of #f on the
    # columns to skip
    #data = data(:,[1 2 6]);
    
    return data


def readANSYSoutputs(filename = None,ncol = None): 
    fid = open(filename,'r')
    if (fid == - 1):
        raise Exception('Could not open file "%s"',filename)
    
    format = '%f' * ncol
    
    file = fileread(filename)
    
    lines = strsplit(file,'\n')
    
    imax = len(lines)
    
    tempdata = np.zeros((imax,ncol))
    
    #have been read, then the number of lines with numeric data will be known.
    #The remaining zeros will be truncated.
    
    ct = 0
    
    for i in range(imax):
        tline = fgetl(fid)
        try:
            a = np.loadtxt(tline,format)
            b = cell2mat(a)
            # in that line.
            if not len(b)==0 :
                ct = ct + 1
                tempdata[ct,:] = b
        finally:
            pass
        #     disp(tline)
        #     if ischar(tline)
        #         ct=ct+1
        #     end
        #data=[data; np.loadtxt(tline,'#f #f #f #f')];
    
    fid.close()
    data = tempdata[0:ct+1,:]
    
    return data

def readANSYSStrains(filename, flag): 
    # readANSYSStrains  Read an ANSYS POST1 element table listing.
    # **********************************************************************
    # *                   Part of the SNL NuMAD Toolbox                    *
    # * Developed by Sandia National Laboratories Wind Energy Technologies *
    # *             See license.txt for disclaimer information             *
    # **********************************************************************
    # data = readANSYSStrains(filename)
    # Read an ANSYS POST1 element table listing.
    # Usage: data = read_strains(filename,flag)
    #  where FILENAME is file name string, default 'Strains.txt'
    #  and flag is either NODE or ELEM
    #        DATA is 7-column matrix [ELEM, EPELX, EPELY, EPELZ, EPELXY, EPELYZ, EPELXZ]
    
    # hard-code filename if not specified
    if filename is None:
        filename = 'Strains.txt'
    
    #     # user select filename if not specified
    #     if ~exist('filename','var') || isempty(filename)
    #         [fn,pn] = uigetfile( ...
    #             {'*.txt','Text files(*.txt)'; ...
    #             '*.*','All files (*.*)'},...
    #             'Select ANSYS element list',defaultfn);
    #         if isequal(fn,0) || isequal(pn,0)
    #             disp('Operation canceled by user.')
    #             return;
    #         end
    #         filename = fullfile(pn,fn);
    #     end
    
    # Open the file and read the entire contents
    with open(filename, 'rb') as fid:
        if (fid == - 1):
            raise Exception('Could not open file "%s"',filename)
        # filecontents = np.transpose(fread(fid,inf,'uint8=>char'))
        filecontents = fid.read()

    #assignin('base','filecontents',filecontents);  #debugging
    
    # process the tables
    NCOLS = 7
    data = []
    pat = flag+'\s*EPELX\s*EPELY\s*EPELZ\s*EPELXY\s*EPELYZ\s*EPELXZ\s*'
    tbl_hdrs = regexp(filecontents,pat)
    
    tbl_hdrs[end() + 1] = np.asarray(filecontents).size
    
    #assignin('base','tbl_hdrs',tbl_hdrs);  #debugging
    for kTbl in range(tbl_hdrs.size-1):
        tbl = filecontents(np.arange(tbl_hdrs(kTbl),tbl_hdrs(kTbl + 1) - 1+1))
        tbl = regexprep(tbl,pat,'')
        data.append([np.loadtxt(tbl,np.matlib.repmat(' %f',1,NCOLS))])
    
    data = cell2mat(data)
    return data