% Generate QMLib dataset

%clc;
%clear;
inputdirname='MaxCut/Beasley/Beasley';
outputdirname='Beasley';
filenames=dir([inputdirname,'/be*'])';
nbprobs=length(filenames);
if exist(outputdirname,'dir')==0
    mkdir(outputdirname);
end

count=1;
for file=filenames
    fprintf('Loading problem %s ... %d / %d\n',file.name, count, nbprobs);
    [~,~,Q,l,s]=loadMaxCut([file.folder,'/',file.name]);
    N=size(Q,1);
    d=2;
    count=count+1;
    
    save([outputdirname,'/',file.name,'.mat'],'Q','l','s','N','d');
end
