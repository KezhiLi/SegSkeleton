% Input: gray image in .hdf5
% Output:
%

%% add path of segworm algorithms
segworm_path = 'C:\Kezhi\MyCode!!!\ManualVideos\SegSkeleton';
addpath(genpath([segworm_path,'.']));

% current path and folder
Av_folder = 'Ave\';
Ev_folder = 'SegTif\';
path_from = ['X:\Kezhi\DataSet\AllFiles\OutSource_files\All_Label\Tif\'];
path_Av = [path_from(1:end-4),Av_folder];
path_Ev = [path_from(1:end-4),Ev_folder];

root_folder = genpath([path_from,'.']);

file=dir([path_from,'*.tif']);
num_file = size(file,1);
isNormalized = 1;
skeleton = {};
MSE = [];

all_file = subdir('X:\Kezhi\DataSet\AllFiles\MaskedVideos\*.hdf5');
num_all_file = size(all_file,1);

num_file = 10;

num_sec = 60;
diff_Av_stand = zeros(num_file,num_sec);
diff_Ev_stand = zeros(num_file,num_sec);
%%
res  = cell(num_file,5);
for nf = 1:  num_file;
    disp([num2str(nf),'/',num2str(num_file)])
    % obtain the filename
    tif_file_ori = file(nf).name(1:end-6);
    
    try
        %% about Av
        % read Av
        file_Av =  dir([path_Av,tif_file_ori,'*.hdf5']);
        if abs(size(file_Av,1)-1)>0.1
            error('find 0 or more than 1 file with the similar names')
        end
        Av_info = h5info([path_Av,file_Av.name]);
        Av_skeleton_full = h5read([path_Av,file_Av.name],'/skeleton');
        
        candid_ind = [];
        % find the original 1min hdf5 file in another folder
        for ii = 1:num_all_file;
            if strfind(all_file(ii).name, file_Av.name(1:end-15))
                candid_ind = [candid_ind, ii];
            end
        end
        % if find more than 1 or 0 file, there is an error
        if abs(length(candid_ind)-1)> 1e-5
            error('there are more than 1 file with the target file name')
        else
            % create the frames for comparison from full frames
            ori_file = all_file(candid_ind).name;
            %Av_timefile_info = h5info(ori_file);
            Av_time_pos = h5read(ori_file,'/vid_time_pos');
            % if the length of time stamps has problems
            if abs(size(Av_skeleton_full,3)-length(Av_time_pos))>1e-5
                disp([file(nf).name(1:end-6), 'has timestamp problems.'])
                % write to error file txt
                fileID = fopen('timestamp_pro_files.txt','a');
                fprintf(num2str(nf),'%s ',file(nf).name(1:end-6));
                fclose(fileID);
            else
                % write the time stamp to the skeleton hdf5 file
                h5write([path_Av,file_Av.name], '/timestamp/time', Av_time_pos);
                num_frame = floor(Av_time_pos(end))+1;
                Av_skeleton = zeros(size( Av_skeleton_full,1),size( Av_skeleton_full,2),num_frame);
                Av_ske_ind = 2;
                Av_skeleton(:,:,1) = Av_skeleton_full(:,:,1);
                for jj = 2:length(Av_time_pos);
                    if floor(Av_time_pos(jj))~=floor(Av_time_pos(jj-1))
                        Av_skeleton(:,:,Av_ske_ind) = Av_skeleton_full(:,:,jj);
                        Av_ske_ind = Av_ske_ind +1;
                    end
                end
            end
        end
        % fix the head/tail mistake frame by frame
        for ii = 2:size(Av_skeleton,3);
            if sum(sum(abs(Av_skeleton(:,:,ii) - Av_skeleton(:,:,ii-1))))> ...
                    sum(sum(abs(Av_skeleton(:,:,ii) - Av_skeleton(:,end:-1:1,ii-1))))
                Av_skeleton(:,:,ii) = Av_skeleton(:,end:-1:1,ii);
            end
        end
        
        %% about EV
        % read Ev
        file_Ev =  dir([path_Ev,tif_file_ori,'*.mat']);
        if abs(size(file_Ev,1)-1)>0.1
            error('find 0 or more than 1 file with the similar names')
        end
        load([path_Ev,file_Ev.name],'seg_skeleton');
        
        %% do comparison
        Ev_skeleton = zeros(size(Av_skeleton));
        len_ske = size(Av_skeleton,2);
        MSE{nf} = zeros(size(Av_skeleton,3),1);
        
        % interpolate then calculate the difference between Av and Ev
        for ii = 1: size(Av_skeleton,3)
            t = 1:size(seg_skeleton{ii},1);
            ts = 1:((size(seg_skeleton{ii},1)-1)/(len_ske-1)):size(seg_skeleton{ii},1);
            ske_xy = seg_skeleton{ii};
            if abs(ske_xy(1,1)-1)<1e-5 & isnan(Av_skeleton(1,1,ii))
                MSE{nf,ii} = -1;  % both Ev and Av are NaN
            elseif abs(ske_xy(1,1)-1)<1e-5
                MSE{nf,ii} = -2;  % means Ev is NaN
            elseif isnan(Av_skeleton(1,1,ii))
                MSE{nf,ii} = -3;  % means Av is NaN
            else
                % interpolate Ev to the same size of Av
                Ev_skeleton(1,:,ii) = interp1(t,ske_xy(:,1)',ts,'spline');
                Ev_skeleton(2,:,ii) = interp1(t,ske_xy(:,2)',ts,'spline');
                % inverse the skeleton of Ev in each frame if necessary
                if ii ==1
                    if sum(sum(abs(Av_skeleton(:,:,ii) - Ev_skeleton(:,:,ii))))> ...
                            sum(sum(abs(Av_skeleton(:,:,ii) - Ev_skeleton(:,end:-1:1,ii))))
                        Ev_skeleton(:,:,ii) = Ev_skeleton(:,end:-1:1,ii);
                    end
                else
                    error_2 = sum(sum(abs(Ev_skeleton(:,:,ii) - Ev_skeleton(:,:,ii-1))))- ...
                        sum(sum(abs(Ev_skeleton(:,:,ii) - Ev_skeleton(:,end:-1:1,ii-1))));
                    error_3 = sum(sum(abs(Av_skeleton(:,:,ii) - Ev_skeleton(:,:,ii))))- ...
                        sum(sum(abs(Av_skeleton(:,:,ii) - Ev_skeleton(:,end:-1:1,ii))));
                    if error_2 > 0 | error_3 > 0
                        Ev_skeleton(:,:,ii) = Ev_skeleton(:,end:-1:1,ii);
                    end
                end
                MSE{nf,ii} = sum(sum((Ev_skeleton(:,:,ii) - Av_skeleton(:,:,ii)).^2))/len_ske;
            end
        end
        %% save all skeleton to result
        res{nf,1} = result{nf,4};
        res{nf,2} = Av_skeleton;
        res{nf,3} = Ev_skeleton;
        if sum(result{nf,3})>0
            for ii = 1:size(Av_skeleton,3)
                if  MSE{nf,ii}>0 & ~isnan(result{nf,4}(1,1,ii))
                    if sum(sum(abs(Av_skeleton(:,:,ii) - result{nf,4}(:,:,ii))))<sum(sum(abs(Av_skeleton(:,:,ii) - result{nf,4}(:,end:-1:1,ii))))
                        diff_Av_stand(nf,ii) = sum(sum((Av_skeleton(:,:,ii) - result{nf,4}(:,:,ii)).^2))/len_ske;
                        diff_Ev_stand(nf,ii) = sum(sum((Ev_skeleton(:,:,ii) - result{nf,4}(:,:,ii)).^2))/len_ske;
                    else
                        diff_Av_stand(nf,ii) = sum(sum((Av_skeleton(:,:,ii) - result{nf,4}(:,end:-1:1,ii)).^2))/len_ske;
                        diff_Ev_stand(nf,ii) = sum(sum((Ev_skeleton(:,:,ii) - result{nf,4}(:,end:-1:1,ii)).^2))/len_ske;
                    end
                end
            end
            res{nf,4} = diff_Av_stand(nf,:);
            res{nf,5} = diff_Ev_stand(nf,:);
        else
            res{nf,4} = NaN;
            res{nf,5} = NaN;
        end
    catch ME
        
    end
end

comp_mtx= [sum(diff_Av_stand,2),sum(diff_Ev_stand,2)]

for ii = 6

    file_nam = file(ii).name(1:end-6);
    for kk = 1: size(all_file)
        if findstr(all_file(kk).name,file_nam);
            break;
        end
    end

    % read .hdf5
    this_file = all_file(kk).name;
    cur_img_ful_info = h5info(this_file);
    cur_img_full = h5read(this_file,'/mask');
    cur_img_time = h5read(this_file,'/vid_time_pos');

    count = 1;
    jj = 1;
    pp =1;
    pre_timeStamp = -1;
   while (jj<=size(res{ii,1},3))&(pp<=size(cur_img_time,1))

        cur_timeStamp = cur_img_time(pp);
        % when current time stamp increase to anther second, save the frame
        if floor(cur_timeStamp)>floor(pre_timeStamp);
            %img = imread([path_from,file(ii).name], jj);
            img =  cur_img_full(:,:,pp);
            imshow(double(img'),[])
            hold on,
            plot(res{ii,1}(1,:,jj), res{ii,1}(2,:,jj), 'g*')
            plot(res{ii,2}(1,:,jj), res{ii,2}(2,:,jj), 'r*')
            plot(res{ii,3}(1,:,jj), res{ii,3}(2,:,jj), 'b*')
            jj = jj + 1;
        end
        pre_timeStamp = cur_timeStamp;
        pp = pp+1;
    end
end





