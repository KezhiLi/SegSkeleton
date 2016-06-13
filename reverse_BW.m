%% add path of segworm algorithms
segworm_path = 'C:\Kezhi\WormTrackingSoftware\SegWorm-master\SegWorm-master';
addpath(genpath([segworm_path,'.']));

% current path and folder
folder = 'skeleton\';
path_from = ['X:\Kezhi\DataSet\AllFiles\OutSource_files\All_Label\Tif\'];
path_to = [path_from,folder];

root_folder = genpath([path_from,'.']);

file=dir([path_from,'*.tif']);
num_file = size(file,1);
isNormalized = 1;
skeleton = {};


for nf = 381;
    disp([num2str(nf),'/',num2str(num_file)])
    
    tif_file_ori = file(nf).name(1:end-4);
    
    %% uniform the name accordinly
    if ~strcmp(tif_file_ori(end-2:end),')-1')
        if strcmp(tif_file_ori(end),')')
            movefile([path_from,file(nf).name],[path_from,tif_file_ori,'-1.tif']);
            tif_file = [tif_file_ori,'-1'];
        elseif strcmp(tif_file_ori(end-1),'-')
            movefile([path_from,file(nf).name],[path_from,tif_file_ori(1:end-1),'1.tif']);
            tif_file = [tif_file_ori(1:end-1),'1'];
        elseif strcmp(tif_file_ori(end-2),'-')
            movefile([path_from,file(nf).name],[path_from,tif_file_ori(1:end-3),'-1.tif']);
            tif_file = [tif_file_ori(1:end-3),'-1'];
        end
    else
        tif_file = tif_file_ori;
    end

    fileRead=[path_from,tif_file,'.tif'];
    tif_path = [path_from,tif_file,'.tif'];
    tiff_info = imfinfo(tif_path);
  
    %% calculate the skeleton and save it
    fileWrite=[path_from,tif_file,'-1.tif'];
    gray = round(max(max(img1))/2);
    
    for ii=1:size(tiff_info, 1)
        img=imread(fileRead,ii);
        if nf == 381
            img_ori = img;
            img = zeros(size(img));
            img(img_ori==1)=1;
            img(img_ori==0)=255;
            img(img_ori==2)=0;
        end
        if ii==1
            imwrite(img,fileWrite,'tif');
        else
            imwrite(img,fileWrite,'tif','WriteMode','append');
        end
    end
end
