%% Main_2_plot_twin_op
% plot the twin order parameters for NPs with different structures
% try to rotate the figures to see the zone axis and twin boundaries

clear; clc; close all;
%%
load('HEA_twin_op.mat')
load('MEA_twin_op.mat')
%%
MEA_ind = [1,4,2];
for i = 1:3
    ind = MEA_ind(i);
    figure(20+i); clf; hold on;
    set(gcf,'position',[50+550*i-550,350,550,700]);
    twin_op = MEA_op_arr{ind,3};
    twin_op = min(max(twin_op,-3),3);
    twin_op_size = 80*((twin_op-min(twin_op))/(max(twin_op)-min(twin_op))+0.01);
    load(sprintf('model_MEA%d.mat',ind))

    model_o = double(model_refined);
    shp = alphaShape(model_o(1,:)',model_o(2,:)',model_o(3,:)',3.493);
    bf = boundaryFacets(shp);
    bf_id = true(size(model_o,2),1);
    bf_id(unique(bf)) = 0;

    scatter3(model_o(1,bf_id),model_o(2,bf_id),model_o(3,bf_id),...
        twin_op_size(bf_id),twin_op(bf_id),...
        'filled','marker','o',...
        'markeredgecolor',[0 0 0])

    axis image; caxis([-3,3]);
    colormap('parula'); axis equal off;
%     xlim([-30,30]); ylim([-30,30]); zlim([-30,30]);
end

%%
HEA_ind = [1,2];
for i = 1:2
    ind = HEA_ind(i);
    figure(23+i); clf; hold on;
    set(gcf,'position',[50+550*i-550,50,550,700]);
    twin_op = HEA_op_arr{ind,3};
    twin_op = min(max(twin_op,-3),3);
    twin_op_size = 80*((twin_op-min(twin_op))/(max(twin_op)-min(twin_op))+0.01);
    load(sprintf('model_HEA%d.mat',ind))

    model_o = double(model_refined);
    shp = alphaShape(model_o(1,:)',model_o(2,:)',model_o(3,:)',3.647);
    bf = boundaryFacets(shp);
    bf_id = true(size(model_o,2),1);
    bf_id(unique(bf)) = 0;

    scatter3(model_o(1,bf_id),model_o(2,bf_id),model_o(3,bf_id),...
        twin_op_size(bf_id),twin_op(bf_id),...
        'filled','marker','o',...
        'markeredgecolor',[0 0 0])

    axis image; caxis([-3,3]);
    colormap('parula'); axis equal off;
%     xlim([-30,30]); ylim([-30,30]); zlim([-30,30]);
end