% Input: manually labeled binary tif
% Output: make .tif file with equal size for every frame: because it is
% found that sometimes one or several frames in a .tif file has different
% size to other frame. This script is to make sure that all frames are with
% size [480, 640].
%

%% add path of segworm algorithms
segworm_path = 'C:\Kezhi\MyCode!!!\ManualVideos\SegSkeleton';
addpath(genpath([segworm_path,'.']));

% current path and folder
Av_folder = 'skeleton_Ave\';
Ev_folder = 'skeleton_segworm\';
path_from = ['X:\Kezhi\DataSet\AllFiles\OutSource_files\All_Label\Tif\'];
path_Av = [path_from,Av_folder];
path_Ev = [path_from,Ev_folder];

root_folder = genpath([path_from,'.']);

file_tif=dir([path_from,'*.tif']);
num_file = size(file_tif,1);

for nf = 1:  num_file;
    disp([num2str(nf),'/',num2str(num_file)])
    % if the frame number is other than 60
    if size(result{nf,3},1)~=60      
        flag = 0;
        % obtain the filename
        tif_file_ori = file_tif(nf).name;
        tiff_info = imfinfo([path_from,tif_file_ori]);
        fileWrite = [path_from,tif_file_ori(1:end-4),'-resize.tif'];
        
        for ii = 1:size(tiff_info, 1)
            tiff = double(imread([path_from,tif_file_ori],ii));
            I = tiff/max(max(tiff));
            %         imshow(tiff,[]);
            %         frame=getframe;
            %         im=frame2im(frame);
            %         [I,map]=rgb2ind(im,256);
            if sum(abs(size(tiff)-[480, 640]))>0
                II =(zeros(480, 640));
                II(1:size(I,1),1:size(I,2)) = I;
                I = II;
                flag = 1;
            end
            if ii==1
                imwrite(I,fileWrite,'tif');
            else
                imwrite(I,fileWrite,'tif','WriteMode','append');
            end
            
        end
        if flag < 1
            delete(fileWrite);
        end
    end
end