% Input: manually labeled binary tif
% Output: skeleton points in .mat format by using segworm 


clear

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
isNormalized = 1;
skeleton = {};
MSE = [];

%
num_file = 10;

result  = cell(num_file,4);
for nf = 1:  num_file;
    disp([num2str(nf),'/',num2str(num_file)])
    % obtain the filename
    tif_file_ori = file_tif(nf).name(1:end-6);
    
    % read Av
    file_Av =  dir([path_Av,tif_file_ori,'*.hdf5']);
    if abs(size(file_Av,1)-1)>0.1
        error('find 0 or more than 1 file with the similar names')
    end
    Av_info = h5info([path_Av,file_Av.name]);
    Av_skeleton = h5read([path_Av,file_Av.name],'/skeleton'); 
    % fix the head/tail mistake frame by frame
    
    for ii = 2:size(Av_skeleton,3);
        buff = 1;
        if ~isnan(Av_skeleton(1,1,ii))
            error_1 = sum(sum(abs(Av_skeleton(:,:,ii) - Av_skeleton(:,:,ii-buff))))- ...
                    sum(sum(abs(Av_skeleton(:,:,ii) - Av_skeleton(:,end:-1:1,ii-buff))));
            while isnan(error_1)& ii-buff >1.5
                buff = buff+1;
                error_1 = sum(sum(abs(Av_skeleton(:,:,ii) - Av_skeleton(:,:,ii-buff))))- ...
                    sum(sum(abs(Av_skeleton(:,:,ii) - Av_skeleton(:,end:-1:1,ii-buff))));
            end
            if error_1 >0
                Av_skeleton(:,:,ii) = Av_skeleton(:,end:-1:1,ii);
            end
        end
    end
    
    % read Ev
    file_Ev =  dir([path_Ev,tif_file_ori,'*.mat']);
    if abs(size(file_Ev,1)-1)>0.1
        error('find 0 or more than 1 file with the similar names')
    end
    load([path_Ev,file_Ev.name],'skeleton');
    
    Ev_skeleton = NaN(size(Av_skeleton));
    len_ske = size(Av_skeleton,2);
    MSE{nf} = zeros(size(Av_skeleton,3),1);

    % interpolate then calculate the difference between Av and Ev
    for ii = 1: size(Av_skeleton,3)
        t = 1:size(skeleton{ii},1);
        ts = 1:((size(skeleton{ii},1)-1)/(len_ske-1)):size(skeleton{ii},1);
        ske_xy = skeleton{ii};
       
        if abs(ske_xy(1,1)-1)<1e-5 & isnan(Av_skeleton(1,1,ii))
            MSE{nf}(ii) = -1;  % both Ev and Av are NaN
        elseif abs(ske_xy(1,1)-1)<1e-5
            MSE{nf}(ii) = -2;  % means Ev is NaN
        elseif isnan(Av_skeleton(1,1,ii))
            MSE{nf}(ii) = -3;  % means Av is NaN
            Ev_skeleton(1,:,ii) = interp1(t,ske_xy(:,2)',ts,'spline');
            Ev_skeleton(2,:,ii) = interp1(t,ske_xy(:,1)',ts,'spline');
        else
            Ev_skeleton(1,:,ii) = interp1(t,ske_xy(:,2)',ts,'spline');
            Ev_skeleton(2,:,ii) = interp1(t,ske_xy(:,1)',ts,'spline');
            % inverse the skeleton of Ev in each frame if necessary
%            if ii ==1
                if sum(sum(abs(Av_skeleton(:,:,ii) - Ev_skeleton(:,:,ii))))> ...
                        sum(sum(abs(Av_skeleton(:,:,ii) - Ev_skeleton(:,end:-1:1,ii))))
                    Ev_skeleton(:,:,ii) = Ev_skeleton(:,end:-1:1,ii);
                end
%             else
%                 error_2 = sum(sum(abs(Ev_skeleton(:,:,ii) - Ev_skeleton(:,:,ii-1))))- ...
%                         sum(sum(abs(Ev_skeleton(:,:,ii) - Ev_skeleton(:,end:-1:1,ii-1))));
%                 error_3 = sum(sum(abs(Av_skeleton(:,:,ii) - Ev_skeleton(:,:,ii))))- ...
%                         sum(sum(abs(Av_skeleton(:,:,ii) - Ev_skeleton(:,end:-1:1,ii))));    
%                 if error_2 > 0 | error_3 > 0
%                    Ev_skeleton(:,:,ii) = Ev_skeleton(:,end:-1:1,ii);
%                end       
%         end
                MSE{nf}(ii) = sum(sum((Ev_skeleton(:,:,ii) - Av_skeleton(:,:,ii)).^2))/len_ske;
        end
    end
    result{nf,1} = Av_skeleton;
    result{nf,2} = Ev_skeleton;
    result{nf,3} = MSE{nf}(:);
    result{nf,4} = (Av_skeleton+Ev_skeleton)/2;
end

for ii = 1
    for jj = 1: 60
        img = imread([path_from,file_tif(ii).name], jj);
        imshow(double(img))
        hold on,
        plot(result{ii,4}(1,:,jj), result{ii,4}(2,:,jj), 'g*')
        plot(result{ii,1}(1,:,jj), result{ii,1}(2,:,jj), 'r*')
        plot(result{ii,2}(1,:,jj), result{ii,2}(2,:,jj), 'b*')
    end
end