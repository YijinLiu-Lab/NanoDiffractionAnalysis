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
    ref = (squeeze(Img_xrf(s-1,:,:))+0.01);
    ref(~isfinite(ref)) = 0;
    moving = (squeeze(Img_xrf(s,:,:))+0.01);
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
roi_aligned_sum = squeeze(sum((align_roi),1));

%%
mask = roi_aligned_sum;
mask(mask < 100) = 0; % modify
mask(mask ~= 0) = 1;
figure,imshow(mask,[])
title('mask')


info = imfinfo('bare_LCO_heat_raw/diff/rock_131477_diff_data.tif'); % modify
slice_diff = size(info,1);
width_diff = info.Width;
height_diff = info.Height;
data_masscenter_x = zeros(num_scan,m,n);
data_masscenter_y = zeros(num_scan,m,n);
data_FWHM_x = zeros(num_scan,m,n);
data_FWHM_y = zeros(num_scan,m,n);
data_histcenter_x = zeros(num_scan,m,n);
data_histcenter_y = zeros(num_scan,m,n);

% figure,imshow(squeeze(data_masscenter_x(1,:,:)),[])
% figure,imshow(squeeze(data_histcenter_x(1,:,:)),[])
% figure,imshow(squeeze(data_FWHM_x(1,:,:)),[])
for ii = 1:num_scan
%     ii = 1;
    data_slice = zeros(slice_diff,height_diff,width_diff);
    parfor iii = 1:slice_diff
        disp(ii);
        disp(iii);
        data_slice(iii,:,:) = double(imread(strcat('bare_LCO_heat_raw/diff/rock_',num2str(scan_list(ii)),'_diff_data.tif'),iii));
    end
    
    for a = 1:m
        for b = 1:n
            if mask(a,b) == 1
                disp(a);
                disp(b);
                data_frame = data_slice(:,:,(a-1)*n+b);
%                 figure,imshow(data_frame,[])
                % frame center and FWHM
                hist_x = sum(data_frame,2);
                hist_x_smooth=smooth(hist_x,20);
                cmx_hist = find(hist_x_smooth == max(hist_x_smooth),1);
                hist_x_FWHM = (hist_x_smooth - 0.5 * max(hist_x_smooth)).^2;
                hist_x_FWHM_l = hist_x_FWHM(1:cmx_hist);
                hist_x_FWHM_r = hist_x_FWHM(cmx_hist:end);
                FWHM_x_l = find(hist_x_FWHM_l == min(hist_x_FWHM_l),1);
                FWHM_x_r = find(hist_x_FWHM_r == min(hist_x_FWHM_r),1,'last');
                FWHM_x = cmx_hist - FWHM_x_l + FWHM_x_r;
                
                hist_y = sum(data_frame,1);
                hist_y_smooth=smooth(hist_y,20);
                cmy_hist = find(hist_y_smooth == max(hist_y_smooth),1);
                hist_y_FWHM = (hist_y_smooth - 0.5 * max(hist_y_smooth)).^2;
                hist_y_FWHM_l = hist_y_FWHM(1:cmy_hist);
                hist_y_FWHM_r = hist_y_FWHM(cmy_hist:end);
                FWHM_y_l = find(hist_y_FWHM_l == min(hist_y_FWHM_l),1);
                FWHM_y_r = find(hist_y_FWHM_r == min(hist_y_FWHM_r),1,'last');
                FWHM_y = cmy_hist - FWHM_y_l + FWHM_y_r;
                
                data_FWHM_x(ii,a,b) = FWHM_x;
                data_FWHM_y(ii,a,b) = FWHM_y;
                data_histcenter_x(ii,a,b) = cmx_hist;
                data_histcenter_y(ii,a,b) = cmy_hist;
                
                [y0,x0] = meshgrid(1:size(data_frame,2),1:size(data_frame,1));
                cmx = sum(x0 .* data_frame) / sum(data_frame);
                cmy = sum(y0 .* data_frame) / sum(data_frame);
                
                data_masscenter_x(ii,a,b) = cmx;
                data_masscenter_y(ii,a,b) = cmy;
            end
        end
    end
end

        
% stop here


rock_img = zeros(m,n);
Cr_masscenter = zeros(m,n);
Ct_masscenter = zeros(m,n);
Cr_histcenter = zeros(m,n);
Ct_histcenter = zeros(m,n);
FWHM_map = zeros(m,n);
sin_alpha = cos(gamma)*sin(delta)/sqrt(1-(cos(gamma)*cos(delta))^2);
cos_alpha = sin(gamma)/sqrt(1-(cos(gamma)*cos(delta))^2);
for a = 1:m
    for b = 1:n
        if mask(a,b) == 1
            disp(a);
            disp(b);
            line = squeeze(align_roi(:,a,b));
            peak_index = find(line==max(line));
            rock_img(a,b) = theta_list(peak_index);
            
            Cr_masscenter(a,b) = (data_masscenter_x(peak_index,a,b)*sin_alpha - data_masscenter_y(peak_index,a,b)*cos_alpha)*detecter_pixel/bragg_angle;
            Ct_masscenter(a,b) = (data_masscenter_x(peak_index,a,b)*cos_alpha + data_masscenter_y(peak_index,a,b)*sin_alpha)*detecter_pixel;
            Cr_histcenter(a,b) = (data_histcenter_x(peak_index,a,b)*sin_alpha - data_histcenter_y(peak_index,a,b)*cos_alpha)*detecter_pixel/bragg_angle;
            Ct_histcenter(a,b) = (data_histcenter_x(peak_index,a,b)*cos_alpha + data_histcenter_y(peak_index,a,b)*sin_alpha)*detecter_pixel;
            FWHM_map(a,b) = data_FWHM_x(peak_index,a,b) * data_FWHM_y(peak_index,a,b);
        end
    end
end

std_Cr = std(Cr_masscenter(Cr_masscenter~=0));
std_Ct = std(Ct_masscenter(Ct_masscenter~=0));

figure,imshow(rock_img,[60 68]) % modify
colormap(gca,'parula')
title('bending y')

Cr_masscenter_nonzero = Cr_masscenter;
Cr_masscenter_nonzero(Cr_masscenter_nonzero == 0) = nan;
figure,imshow(Cr_masscenter_nonzero,[])
colormap(gca,'parula')
title('\Deltad/d')

figure,imshow(Ct_masscenter,[])
colormap(gca,'parula')
title('bending z')

figure,imshow(FWHM_map,[1 10000])
colormap(gca,'jet')
title('FWHM')


% figure,imshow(Cr_histcenter,[])
% colormap(gca,'parula')
% title('\Deltad/d')
% 
% figure,imshow(Ct_histcenter,[])
% colormap(gca,'parula')
% title('bending z')

save conventional_post_heating -v7.3;
return;
out_put_path = 'result_20191106\';

dlmwrite([out_put_path,'NMC_charged_preheating_rock_img.mat'],rock_img);
dlmwrite([out_put_path,'NMC_charged_preheating_rock_img.txt'],rock_img);

dlmwrite([out_put_path,'NMC_charged_preheating_Cr.mat'],Cr_masscenter);
dlmwrite([out_put_path,'NMC_charged_preheating_Cr.txt'],Cr_masscenter);

dlmwrite([out_put_path,'NMC_charged_preheating_Ct.mat'],Ct_masscenter);
dlmwrite([out_put_path,'NMC_charged_preheating_Ct.txt'],Ct_masscenter);
        
