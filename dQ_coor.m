%%%% det_data contains 2D diffraction data at each rotation angle %%%%
%%%% Number of frames must agree with that of rotation angles %%%%%%%%
%%%% There are three coordinates system, lab, crystal and beam %%%%%%%
%%%% Lab: x, y , z are fixed axes in the lab %%%%%%%%%%%%%%%%%%%%%%%%%
%%%% cryst: z-hkl vector; x-rocking direction %%%%%%%%%%%%%%%%%%%%%%%%
%%%% beam: z-beam direction; y-perpendicular to diffraction plane %%%%
%%%% defined by the beam and hkl vector; x-in the diffraction plane %%
%%%% and is perpendicular to the beam %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xq,yq,zq,vq, h] = dQ_coor_v2(params, det_data, coor)
energy = params.energy;
delta = params.delta;
gamma = params.gamma;
% th = params.th;
num_angle = params.num_angle;
% th_rng = params.th_rng;
th_step = params.th_step;
pix = params.pix;
% det_row = params.det_row;
% det_col = params.det_col;
det_dist = params.det_dist;
offset = params.offset;
k = 1e4/(12.398/energy);

sz = size(det_data);
det_row = sz(1);
det_col = sz(2);

% In right-hand convention, both delta and gamma are negative   
delta = -delta*pi/180;
gamma = -gamma*pi/180;
% th = th*pi/180; 
% th_rng = th_rng*pi/180;
th_step = th_step*pi/180;
% Detector base vectors in lab coordinates system;

Mx = [1, 0, 0; 0, cos(delta), -sin(delta); 0, sin(delta), cos(delta)];
My = [cos(gamma), 0, sin(gamma); 0, 1, 0; -sin(gamma), 0, cos(gamma)];

M_D2L = My*Mx; % convert vectors in det coor into lab coor
M_L2D = inv(M_D2L); % vice versa

kx_lab = M_D2L*[1;0;0]*k*(pix/det_dist);
ky_lab = M_D2L*[0;1;0]*k*(pix/det_dist);

% calculate k_0 direction in detector coordinates
k_0 = k*M_L2D*[0;0;1];
% calculate the direction of the h vector in detector coordinates
h = k*[0;0;1]-k_0;
% rocking direction is perpendicular to the plane defined by h and lab_y in
% detector coordinates
rock_z = cross(M_L2D*[0;1;0],h);
kz = -rock_z*th_step;

kz_lab = M_D2L*kz;

% Matrix convert oblique coordiantes (kx, ky, kz) (in unit of pixel) 
% to lab coordinates

M_O2L = [kx_lab,ky_lab,kz_lab];
M_L2O = inv(M_O2L);

[X, Y, Z] = meshgrid(-round(det_col/2):round(det_col/2)-1,...
    -round(det_row/2):round(det_row/2)-1,-round(num_angle/2):round(num_angle/2)-1);
X = X+offset(2);
Y = Y+offset(1);
% Z = Z*norm(kz);

% unit x, y, z vector of crystal in lab coordinates 
ux_cryst = kz_lab/norm(kz_lab);
uz_cryst = M_D2L*h/norm(M_D2L*h);
uy_cryst = cross(uz_cryst,ux_cryst);

% convertion matrix between lab abd crystal
M_C2L = [ux_cryst,uy_cryst,uz_cryst];
M_C2O = M_L2O*M_C2L;
M_O2C = inv(M_C2O);

uz_beam = [0;0;1];
uy_beam = cross(uz_beam,uz_cryst);
ux_beam = cross(uy_beam, uz_beam); 
        
M_B2L = [ux_beam,uy_beam,uz_beam];
M_O2B = M_L2O*M_B2L;
M_B2O = inv(M_O2B);


if strcmp(coor,'lab')
    [pix_sz, xq, yq, zq] = create_grid(X, Y, Z, M_O2L);
    vq = interp3_oblique(X, Y, Z, det_data, M_L2O, xq, yq, zq);
elseif strcmp(coor,'cryst') 
    [pix_sz, xq, yq, zq] = create_grid(X, Y, Z, M_O2C);
    vq = interp3_oblique(X, Y, Z, det_data, M_C2O, xq, yq, zq);
elseif strcmp(coor,'beam')
    [pix_sz, xq, yq, zq] = create_grid(X, Y, Z, M_O2B);
    vq = interp3_oblique(X, Y, Z, det_data, M_B2O, xq, yq, zq);
else
    print('coor must be lab, cryst or beam');
end
end