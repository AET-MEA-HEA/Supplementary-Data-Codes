%% Main_3_strain_map
% example plot of the strain tensor of MEA 1
% the strain tensor is calculated by lattice fitting method described
% elsewhere

clear; clc; close all;
%%
load('strain_map_MEA1_excl_surf.mat');
disp3D_ch = permute(disp3D,[3,1,2,4]);
err3D_ch = permute(err3D,[3,1,2,4]);
disp3D_ch(:,:,:,2:4) = disp3D_ch(:,:,:,[4,2,3]);
err3D_ch(:,:,:,2:4) = err3D_ch(:,:,:,[4,2,3]);
posDisp_ch = posDisp;
posDisp_ch(:,1:3) = posDisp_ch(:,[3,1,2]);
posDisp_ch(:,8:10) = posDisp_ch(:,[10,8,9]);
posDisp_ch(:,11:13) = posDisp_ch(:,[13,11,12]);
new_arr = linspace(min(posDisp_ch(:,3)),max(posDisp_ch(:,3)),20);
new_arr = [-1000,(new_arr(1:end-1)+new_arr(2:end))/2,1000];

for i = 1:numel(new_arr)-1
    layer_ind = posDisp_ch(:,3)>new_arr(i) & posDisp_ch(:,3)<=new_arr(i+1);
    posDisp_ch(layer_ind,14) = i;
end
strain_map_slice(disp3D_ch,err3D_ch,zv,xv,yv,dxyz,posDisp_ch,10:14);