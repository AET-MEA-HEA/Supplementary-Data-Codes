%% Main_1_calculate_local_displacement
% calculate the local displacements for each NPs
% first fitting by breadth-first algorithm provided elsewhere [Nature 592.7852 (2021): 60-64]
% then use the deviation between the nearest neighbours and fitted template

close all; clear; clc;
%%
HEA_fcc_disp_arr = cell(4,1);
for samp_ind = 1:4
filename_1 = sprintf('ptm_HEA%d_excl_surf.mat',samp_ind);
load(filename_1)

neigh_arr_fcc = cellfun(@numel,neigh);
neigh_arr_fcc(neigh_arr_fcc>12)=12;
%
dist_arr_fcc = zeros(numel(mro_results),12);
for i = 1:numel(mro_results)
mro_fcc_t = mro_results{i}(1);

pos_o = mro_fcc_t.pos-mro_fcc_t.origin;
grid_pt_num = size(mro_fcc_t.grid_pts,1)-1;

[pos_grid,v] = calc_r_grid(mro_fcc_t.param,mro_fcc_t.grid_pts,mro_fcc_t.type);
dist_temp = sqrt(sum((pos_grid - pos_o).^2,2));
dist_temp(1) = [];
dist_arr_fcc(i,1:numel(dist_temp)) = dist_temp;
if grid_pt_num<neigh_arr_fcc(i)
    dist_arr_fcc(i,numel(dist_temp)+1:end) = NaN;
end
end

fcc_op = nanmean(dist_arr_fcc,2)';
HEA_fcc_disp_arr{samp_ind,1} = fcc_op;
end
%%
MEA_fcc_disp_arr = cell(4,1);
for samp_ind = 1:4
filename_1 = sprintf('ptm_MEA%d_excl_surf.mat',samp_ind);
load(filename_1)

neigh_arr_fcc = cellfun(@numel,neigh);
neigh_arr_fcc(neigh_arr_fcc>12)=12;
%
dist_arr_fcc = zeros(numel(mro_results),12);
for i = 1:numel(mro_results)
mro_fcc_t = mro_results{i}(1);

pos_o = mro_fcc_t.pos-mro_fcc_t.origin;
grid_pt_num = size(mro_fcc_t.grid_pts,1)-1;

[pos_grid,v] = calc_r_grid(mro_fcc_t.param,mro_fcc_t.grid_pts,mro_fcc_t.type);
dist_temp = sqrt(sum((pos_grid - pos_o).^2,2));
dist_temp(1) = [];
dist_arr_fcc(i,1:numel(dist_temp)) = dist_temp;
if grid_pt_num<neigh_arr_fcc(i)
    dist_arr_fcc(i,numel(dist_temp)+1:end) = NaN;
end
end

fcc_op = nanmean(dist_arr_fcc,2)';
MEA_fcc_disp_arr{samp_ind,1} = fcc_op;
end
%%
hist_arr = 0:0.05:100;
%%
delta_d_cell_new = MEA_fcc_disp_arr(1);
sigStrain = cell2mat(delta_d_cell_new');
[N1,X] = hist(sigStrain,hist_arr);
N1 = N1/sum(N1);
fprintf('deviation of MEA 1 —— mean: %.02f, std: %.02f\n',mean(sigStrain),std(sigStrain));
delta_d_cell_new = MEA_fcc_disp_arr(2);
sigStrain = cell2mat(delta_d_cell_new');
[N2,~] = hist(sigStrain,hist_arr);
N2 = N2/sum(N2);
fprintf('deviation of MEA 2 —— mean: %.02f, std: %.02f\n',mean(sigStrain),std(sigStrain));
%%
delta_d_cell_new = HEA_fcc_disp_arr(1);
sigStrain = cell2mat(delta_d_cell_new');
[N3,~] = hist(sigStrain,hist_arr);
N3 = N3/sum(N3);
fprintf('deviation of HEA 1 —— mean: %.02f, std: %.02f\n',mean(sigStrain),std(sigStrain));
delta_d_cell_new = HEA_fcc_disp_arr(2);
sigStrain = cell2mat(delta_d_cell_new');
[N4,~] = hist(sigStrain,hist_arr);
N4 = N4/sum(N4);
fprintf('deviation of HEA 2 —— mean: %.02f, std: %.02f\n',mean(sigStrain),std(sigStrain));
%%
figure(12); clf; hold on; set(gcf,'position',[300,250,300,400]);
b1 = bar(X, N1); xlim([0,0.8]); ylim([0,0.3]);
b1.EdgeColor = 'k'; box on; title('MEA 1');

figure(13); clf; hold on; set(gcf,'position',[600,250,300,400]);
b1 = bar(X, N2); xlim([0,0.8]); ylim([0,0.3]);
b1.EdgeColor = 'k'; box on; title('MEA 2');

figure(14); clf; hold on; set(gcf,'position',[900,250,300,400]);
b3 = bar(X, N3); xlim([0,0.8]); ylim([0,0.3]);
b3.EdgeColor = 'k'; box on; title('HEA 1');

figure(15); clf; hold on; set(gcf,'position',[1200,250,300,400]);
b1 = bar(X, N4); xlim([0,0.8]); ylim([0,0.3]);
b1.EdgeColor = 'k'; box on; title('HEA 2');
