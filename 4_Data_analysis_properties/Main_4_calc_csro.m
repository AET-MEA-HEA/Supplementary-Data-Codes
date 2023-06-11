%% Main_4_calc_csro
% plot the twin order parameters for NPs with different structures
% try to rotate the figures to see the zone axis and twin boundaries

clear; clc; close all;
%%
load('MEA_bondlength.mat')
for samp_ind = 1:6

    model_filename = sprintf('model_MEA%d.mat',samp_ind);
    atoms_filename = sprintf('atoms_MEA%d.mat',samp_ind);
    atom_model = importdata(model_filename);
    atom_model = double(atom_model');
    atom_types = importdata(atoms_filename);
    atom_types = double(atom_types');

    Q_bondTh = MEA_bondlength(samp_ind);

    shp = alphaShape(atom_model(:,1),atom_model(:,2),atom_model(:,3),Q_bondTh);
    bf = unique(boundaryFacets(shp));

    atom_model(bf,:) = [];
    atom_types(bf,:) = [];

    neigh = cell(size(atom_model,1),1);

    for i=1:size(atom_model,1)
        dist_arr = pdist2(atom_model,atom_model(i,:));
        neigh{i} = find(dist_arr >0 & dist_arr < Q_bondTh);
    end

    concentration = zeros(3,1);
    for type_ind = 1:3
        concentration(type_ind) = sum(atom_types == type_ind) / numel(atom_types);
    end

    local_alpha_arr = zeros(3,3,length(atom_types));

    for N = 1:length(atom_types)

        temp_model = atom_model(N,:);
        temp_dist = pdist2(temp_model, atom_model);
        atom_ind_radius = temp_dist<Q_bondTh;

        for i = 1:3
            indi = find(atom_types(:)==i & atom_ind_radius(:));
            ind = [];
            for k = 1:numel(indi)
                ind = [ind; neigh{indi(k)}];
            end
            temp_neigh_types_all = atom_types(ind);
            for j = 1:3
                local_alpha_arr(i,j,N) = ( sum(temp_neigh_types_all==j) / length(ind) - concentration(j) ) /...
                    ( (i==j) - concentration(j) );
            end
        end

    end
    local_alpha_arr(isnan(local_alpha_arr))=0;
    Chemsro_alpha = ( local_alpha_arr + permute(local_alpha_arr, [2, 1, 3]) ) / 2; 
    save(sprintf('output/csro/Chemsro_MEA%d_excl_surf.mat',samp_ind),'Chemsro_alpha','concentration','bf');

end
%% interpolate into 3D grid
MEA_model = importdata('MEA_model_info.mat');

mat1 = MatrixQuaternionRot([1,-1,0],acosd(dot([0,0,1],[1,1,1])/sqrt(3)));
mat2 = MatrixQuaternionRot([0,0,1],acosd(dot([1,0,0],[1,1,0])/sqrt(2)));
samp_ind = 1;

load(sprintf('Chemsro_MEA%d_excl_surf.mat',samp_ind))
atom_model = MEA_model(samp_ind).model_100.*MEA_model(samp_ind).pixelsize;
atom_model = double(atom_model-mean(atom_model,2));
atom_model = mat2*mat1*atom_model;

pos_max = ceil(max(atom_model,[],2));
pos_min = floor(min(atom_model,[],2));
bf = bf_finder(MEA_model,samp_ind);
atom_model(:,bf) = [];
ind=[1,1;...
    2,2;...
    3,3;...
    1,2;...
    1,3;...
    2,3];
res = 1;
[xq,yq,zq] = ndgrid(pos_min(1):res:pos_max(1),pos_min(2):res:pos_max(2),pos_min(3):res:pos_max(3));
ngauss = MEA_model(samp_ind).Q_bondTh;

for i = 1:6
fgrid = rbfinterp([xq(:)'; yq(:)'; zq(:)'], rbfcreate(atom_model, squeeze(Chemsro_alpha(ind(i,1),ind(i,2),:))','RBFFunction', 'gaussian', 'RBFConstant', ngauss));
fgrid = reshape(fgrid, size(xq));
fgrid(isnan(fgrid))=0;

pos=(atom_model-pos_min+1)'/res;
shp = alphaShape(atom_model(2,:)',atom_model(1,:)',atom_model(3,:)',4);
in = inShape(shp,yq(:),xq(:),zq(:));
fgrid1=fgrid*0-1;
fgrid1(in)=fgrid(in);
fgrid1(fgrid1==-1) = nan;

if i == 1
    fgrid_arr   = fgrid;
    fgrid1_arr  = fgrid1;
else
    fgrid_arr   = cat(4, fgrid_arr, fgrid);
    fgrid1_arr  = cat(4, fgrid1_arr,fgrid1);
end
end
save(sprintf('output/csro/Chemsro_grid_MEA%d_excl_surf.mat',samp_ind),'fgrid_arr','fgrid1_arr');
%%
cutoff_arr = [0.5,0.25,-0.25,-0.5];

type_arr = [1,5,6];
for ind = 1:3
    type_ind = type_arr(ind);
fgrid1 = fgrid1_arr(:,:,:,type_ind);
fgrid2 = imgaussfilt3(fgrid1,0.5);
% fgrid2 = fgrid1;

mat1 = MatrixQuaternionRot([1,-1,0],acosd(dot([0,0,1],[1,1,1])/sqrt(3)));
mat2 = MatrixQuaternionRot([0,0,1],acosd(dot([1,0,0],[1,1,0])/sqrt(2)));

mean_t = mean(atom_model,2);
pos_min=floor(min(atom_model,[],2));
pos = (atom_model-pos_min+3)';

load('cmap.mat')

atom_model1 = atom_model';
shp = alphaShape(pos(:,2),pos(:,1),pos(:,3),4);

figure(60+type_ind); set(gcf,'Position',[100+500*ind-500,50,500,500]); clf; hold on;
p0 = plot(shp);  p0.FaceColor = 'k';  p0.FaceAlpha = 0.02;  p0.EdgeColor = 'none';

p11 = patch(isosurface(fgrid2,cutoff_arr(3)));
p11.FaceColor = [0.95,0.43,0.44];  p11.FaceAlpha = 0.04;  p11.EdgeColor = 'none';
p12 = patch(isosurface(fgrid2,cutoff_arr(4)));
p12.FaceColor = [0.95,0.43,0.44];  p12.FaceAlpha = 0.12;  p12.EdgeColor = 'none';
p21 = patch(isosurface(fgrid2,cutoff_arr(2)));
p21.FaceColor = [0,0.3,0.7];  p21.FaceAlpha = 0.04;  p21.EdgeColor = 'none';
p22 = patch(isosurface(fgrid2,cutoff_arr(1)));
p22.FaceColor = [0,0.2,0.8];  p22.FaceAlpha = 0.12;  p22.EdgeColor = 'none';

axis image off;  view([1,0.1,0.1]);  lighting gouraud;
line([30,30],[3,13],[3,3],'linewidth',5,'color','k')
end