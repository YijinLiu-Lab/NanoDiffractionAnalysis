clear; clc;
filenames = dir('\ply\*.ply');

figure;
parfor i=1:numel(filenames)
    i
    
    a = filenames(i).name;
    
    xind = str2num(a(7:9));
    yind = str2num(a(11:13));
    
     out = load(['\results2\rec_' num2str(i-1) '.mat']);
    
    Feature(:,i) = squeeze(out.feature);
    
end

%%

rng(0); % For reproducibility
idx =kmeans_opt(Feature',15); 
%
numP = max(unique(idx));
CorMap = zeros(140,200);

for i=1:size(Feature,2)
    i
    a = filenames(i).name;
    
    xind = str2num(a(7:9));
    yind = str2num(a(11:13));
    
    CorMap(xind,yind) = idx(i);
end
figure; imshow(CorMap,[]);colormap('jet')

figure;hold on
for j=1:numP
    tmp = zeros(size(CorMap)); tmp(CorMap==j)=1;
    
    [B,L] = bwboundaries(tmp,'noholes');
    for k = 1:length(B)
boundary = B{k};
plot(boundary(:,2), boundary(:,1), 'LineWidth', 2)
end
end
axis off

% optimal cluster number: 13

%%
det_data = pattern05;
coor = 'cryst';
params.energy = 99021; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% unknown
params.delta = 18.1;
params.gamma= 28.5;
params.num_angle =21 ;
% th_rng = params.th_rng;
params.th_step = 0.1;
params.pix = 55; %um
% det_row = params.det_row;
% det_col = params.det_col;
params.det_dist = 0.5*1e6; % 0.5m, to um
params.offset = [0 0];   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% unknown
[xq,yq,zq,vq, h] = dQ_coor_v2(params, det_data, coor);

figure; sliceViewer(vq)
