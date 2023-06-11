function bf = bf_finder(MEA_model_str,samp_ind)

atom_model  = double(MEA_model_str(samp_ind).model_refined');
Q_bondTh    = MEA_model_str(samp_ind).Q_bondTh;
shp         = alphaShape(atom_model(:,1),atom_model(:,2),atom_model(:,3),Q_bondTh);
bf          = unique(boundaryFacets(shp));

end