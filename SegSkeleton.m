% Input: manually labeled binary tif
% Output: skeleton points in .mat format by using segworm 
% 

%% add path of segworm algorithms
segworm_path = 'C:\Kezhi\WormTrackingSoftware\SegWorm-master\SegWorm-master';
addpath(genpath([segworm_path,'.']));

% current path and folder
folder = 'skeleton_segworm\';
path_from = ['X:\Kezhi\DataSet\AllFiles\OutSource_files\All_Label\Tif\'];
path_to = [path_from,folder];

root_folder = genpath([path_from,'.']);

file=dir([path_from,'*.tif']);
num_file = size(file,1);
isNormalized = 1;
skeleton = {};


for nf = 1: num_file;
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
    
    %% fix image size problem (some images with size 641*481)
    img0=imread(fileRead,1);
    if size(img0)~= [640,480] & size(img0)~= [480,640]
        ImgWrite_new=[path_to,tif_file(1:end-1),'-1.tif'];
        
        % edges should be 0, otherwise there is a problem
        offset_mtx =(img0(1:4,1:4)<0.5).*[[1:4];[2:5];[3:6];[4:7]];
        offset_mtx(offset_mtx==0)=10;
        [r,c]=find(offset_mtx==min(min(offset_mtx)));
        if r~=1&c~=1
            for ii=1:size(tiff_info, 1)
                img_old=imread(fileRead,ii);
                if size(img_old,1)>size(img_old,2)
                    img1 = zeros(640,480);
                else
                    img1 = zeros(480,640);
                end
                img_part = img_old(r:end,c:end);
                img1(1:size(img_part,1),1:size(img_part,2)) = img_part;
                if ii==1
                    imwrite(img1,ImgWrite_new,'tif');
                else
                    imwrite(img1,ImgWrite_new,'tif','WriteMode','append');
                end
            end
        else
            error('img size is wrong and cannot be fixed: %s',tif_file);
        end
         movefile([path_from,tif_file,'.tif'],[path_from,'size_problem_file\',tif_file(1:end-1),'99.tif']);
         movefile(ImgWrite_new,[path_from,tif_file,'.tif']);
    end
  
    %% calculate the skeleton and save it
    fileWrite=[path_to,tif_file,'.mat'];
    gray = round(max(max(img1))/2);
    
    for ii=1:size(tiff_info, 1)
        img=imread(fileRead,ii);
        if max(max(img))<0.01
            skeleton{ii} = 0;
            % imshow(img);
        else
            try
                skeleton{ii} =  segWorm_ske(img, ii, isNormalized,0);
            catch ME
                skeleton{ii} = 1; % 1 means skeletonization failed and the worm maybe coils in this frame
            end

            % show the skeleton
%             for jj = 1:size(skeleton{ii},1)
%                 img(skeleton{ii}(jj,1),skeleton{ii}(jj,2))= gray;
%             end
%             hold on, 
%             imshow(img,[]);
%             hold off,
%             pause(1)
        end
    end
    if exist(fileWrite)
        delete(fileWrite);
    end
        save(fileWrite,'skeleton');

end
