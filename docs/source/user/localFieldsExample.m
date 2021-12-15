%First get ansysBladeMaterials
[isoorthoInModel,compsInModel,SkinAreas,app] = getMatrialLayerInfoWithOutGUI(blade);
bladeMatNames=cell(numel(blade.materials),1);

for iMat=1:numel(blade.materials)
    bladeMatNames{iMat}=blade.materials(iMat).name;
end

matPointer=zeros(numel(isoorthoInModel),1);

for iMat=1:numel(isoorthoInModel)
    ansysMPnumber = find(strcmp(isoorthoInModel(iMat),bladeMatNames)==1);
    matPointer(iMat)=ansysMPnumber;
end

ansysBladeMaterials=blade.materials(matPointer);
sections = readANSYSSections('NuMAD/Sections.txt');  %Path to where results are located
elements = readANSYSElem(['NuMAD/Elements.txt']);
myresult=extractFieldsThruThickness('plateStrains-all-1.txt',sections,elements,ansysBladeMaterials,elNo,coordSys)