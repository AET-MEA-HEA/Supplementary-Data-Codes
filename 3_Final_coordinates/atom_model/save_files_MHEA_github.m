model_refined = MEA_model(1).model_refined;
save('model_MEA1.mat','model_refined');
model_refined = MEA_model(4).model_refined;
save('model_MEA2.mat','model_refined');
model_refined = MEA_model(2).model_refined;
save('model_MEA3.mat','model_refined');
model_refined = MEA_model(7).model_refined;
save('model_MEA4.mat','model_refined');
model_refined = MEA_model(5).model_refined;
save('model_MEA5.mat','model_refined');
model_refined = MEA_model(6).model_refined;
save('model_MEA6.mat','model_refined');
%%
model_refined = HEA_model_str_nosurface(4).model_refined;
save('model_HEA1.mat','model_refined');
model_refined = HEA_model_str_nosurface(16).model_refined;
save('model_HEA2.mat','model_refined');
model_refined = HEA_model_str_nosurface(1).model_refined;
save('model_HEA3.mat','model_refined');
model_refined = HEA_model_str_nosurface(12).model_refined;
save('model_HEA4.mat','model_refined');
%%
atom_type = MEA_model(1).atom_type;
save('atoms_MEA1.mat','atom_type');
atom_type = MEA_model(4).atom_type;
save('atoms_MEA2.mat','atom_type');
atom_type = MEA_model(2).atom_type;
save('atoms_MEA3.mat','atom_type');
atom_type = MEA_model(7).atom_type;
save('atoms_MEA4.mat','atom_type');
atom_type = MEA_model(5).atom_type;
save('atoms_MEA5.mat','atom_type');
atom_type = MEA_model(6).atom_type;
save('atoms_MEA6.mat','atom_type');
%%
atom_type = HEA_model_str_nosurface(4).atom_type;
save('atoms_HEA1.mat','atom_type');
atom_type = HEA_model_str_nosurface(16).atom_type;
save('atoms_HEA2.mat','atom_type');
atom_type = HEA_model_str_nosurface(1).atom_type;
save('atoms_HEA3.mat','atom_type');
atom_type = HEA_model_str_nosurface(12).atom_type;
save('atoms_HEA4.mat','atom_type');


