nmdfn='WindPact15MWv1.3.nmd';
bmodes_path='C:\DesignCodes\BModes_v3.00.00\bmodes.exe';
NuMAD2BPE2FASTBlade(nmdfn,bmodes_path)
%     delete FASTBlade.dat
    delete utlsuite*
    delete BPE_SectionData.mat
    delete bmodes.out
    delete shell7bpe.src
%     delete displacement.txt
    delete output.txt
    disp('BPE-related files have been deleted.')


