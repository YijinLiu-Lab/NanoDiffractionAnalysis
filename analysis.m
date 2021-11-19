clear
clc
scan_list = 131477:3:131537; % modify
theta_list = linspace(-39.7,-37.7,21); %modify
num_scan = size(scan_list,2);

delta = 18.1 * pi / 180; % in-plane rotation % modify
gamma = 28.5 *  pi / 180; % out-of-plane rotation  % modify
detecter_pixel = 1*10^(-4); %detecter pixel
bragg_angle = 28.21/180*pi;  % modify

size_xrf = imfinfo('bare_LCO_heat_raw\xrf\rock_131477_xrf.tif');  % modify
m = size_xrf.Height;
n = size_xrf.Width;

% load data
Img_xrf = zeros(num_scan,m,n);
Img_roi = zeros(num_scan,m,n);
for i = 1 : num_scan
    disp(i);
    Img_xrf(i,:,:) = double(imread(strcat('bare_LCO_heat_raw\xrf\rock_',num2str(scan_list(i)),'_xrf.tif')));
    Img_roi(i,:,:) = double(imread(strcat('bare_LCO_heat_raw\roi\rock_',num2str(scan_list(i)),'_roi.tif')));
end

figure,imshow(squeeze(Img_roi(10,:,:)),[])
title('roi')
figure,imshow(squeeze(Img_xrf(10,:,:)),[])
title('xrf')

% alignment
% offset_x_int = zeros(num_scan,1);
% offset_y_int = zeros(num_scan,1);
% offset_x_tmp = zeros(num_scan,1);
% offset_y_tmp = zeros(num_scan,1);
align_xrf = zeros(size(Img_xrf));
align_xrf(1,:,:) = Img_xrf(1,:,:);
align_roi = zeros(size(Img_roi));
align_roi(1,:,:) = Img_roi(1,:,:);
for s = 2 : num_scan
    disp(s);
    [optimizer, metric] = imregconfig('multimodal');
    ref = log(squeeze(Img_xrf(s-1,:,:)));
    ref(~isfinite(ref)) = 0;
    moving = log(squeeze(Img_xrf(s,:,:)));
    moving(~isfinite(moving)) = 0;
    
    [MOVINGREG] = registerImages(moving,ref);
    align_xrf(s,:,:) = MOVINGREG.RegisteredImage;
    align_roi(s,:,:) = imwarp(squeeze(Img_roi(s,:,:)),MOVINGREG.Transformation,'OutputView',imref2d(size(squeeze(Img_roi(1,:,:)))));
    %
    %     tform = imregtform(moving, ref, 'translation', optimizer, metric);
    %     align_xrf(s,:,:) = imwarp(squeeze(Img_xrf(s,:,:)),tform,'OutputView',imref2d(size(squeeze(Img_xrf(1,:,:)))));
    %     align_roi(s,:,:) = imwarp(squeeze(Img_roi(s,:,:)),tform,'OutputView',imref2d(size(squeeze(Img_roi(1,:,:)))));
    %     offset_x_int(s) = floor(tform.T(3,1));
    %     offset_y_int(s) = floor(tform.T(3,2));
    %     offset_x_tmp(s) = tform.T(3,1);
    %     offset_y_tmp(s) = tform.T(3,2);
end
roi_aligned_sum = squeeze(sum(log(align_roi+0.01),1));
% figure,imshow(squeeze(roi_aligned_sum),[])

Img_xrf = permute(align_xrf,[2 3 1]);
Img_roi = permute(align_roi, [2 3 1]);


%%
mask = (roi_aligned_sum);
mask(mask < -30) = 0; % modify
mask(mask ~= 0) = 1;
figure,imshow(mask,[])
title('mask')
%%

info = imfinfo('bare_LCO_heat_raw\diff\rock_131477_diff_data.tif'); % modify
slice_diff = size(info,1);
width_diff = info.Width;
height_diff = info.Height;
data_masscenter_x = zeros(num_scan,m,n);
data_masscenter_y = zeros(num_scan,m,n);
data_FWHM_x = zeros(num_scan,m,n);
data_FWHM_y = zeros(num_scan,m,n);
data_histcenter_x = zeros(num_scan,m,n);
data_histcenter_y = zeros(num_scan,m,n);

% return

%%
data_slice = zeros(height_diff,width_diff,slice_diff,num_scan);

for ii = 1:num_scan
    %     ii = 1;
    ii
    parfor iii = 1:slice_diff
        disp(ii);
        disp(iii);
        data_slice(:,:,iii,ii) = double(imread(strcat('bare_LCO_heat_raw\diff\rock_',num2str(scan_list(ii)),'_diff_data.tif'),iii));
    end
    
end

save allPixel_heat -v7.3;

%%
% figure;
for a = 1:m
    for b = 1:n
        if mask(a,b) == 1 && ~isfile(['heat_movies/pixel_' num2str(a,'%03d') '_' num2str(b,'%03d') '.mat'])
            disp(['pixel ' num2str(a) '_' num2str(b) '']);
            data_frame = squeeze(data_slice(:,(a-1)*n+b,:,:));
            save(['heat_movies/pixel_' num2str(a,'%03d') '_' num2str(b,'%03d') '.mat'], 'data_frame');
            %aux_stackwrite(uint16(data_frame),['heat_movies/pixel_' num2str(a,'%03d') '_' num2str(b,'%03d') '.tif']);
            
            
            data_frame(data_frame>0.05)=0.05;
            
            selectpixel = [];
            w = 1;
            for i=1:236
                for j=1:211
                    for k=1:21
                        if data_frame(i,j,k)>0.0001
                            selectpixel(w,:) = [i j k*10];
                            values(w,:) = data_frame(i,j,k);
                            w=w+1;
                        end
                    end
                end
            end
            ptCloudB = pcdenoise(pointCloud(selectpixel),'NumNeighbors',100,'Threshold',0.8);
            [pcd_de] = PCD_Filtering_LPAICI_mex(single(ptCloudB.Location),single(5),single(2));
            [pcd_de] = PCD_Filtering_LPAICI_mex(single(pcd_de),single(5),single(0.2));
            ptCloudB = pcdenoise(pointCloud(pcd_de),'NumNeighbors',100,'Threshold',0.8);
            %figure; pcshow(ptCloudB);
            
            write_ply_only_pos(ptCloudB.Location, ['ply_heat/pixel_' num2str(a,'%03d') '_' num2str(b,'%03d') '.ply']);
            
        end
    end
end


return
