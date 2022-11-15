 function EPI2T2W_mlb=MY_get_default_coreg_batch_struct(ref, source, other)
    EPI2T2W_mlb{1}.spm.spatial.coreg.estwrite.ref =  ref;
    EPI2T2W_mlb{1}.spm.spatial.coreg.estwrite.source = source;
    EPI2T2W_mlb{1}.spm.spatial.coreg.estwrite.other = other;
    EPI2T2W_mlb{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
    EPI2T2W_mlb{1}.spm.spatial.coreg.estwrite.eoptions.sep = [2 1];
    EPI2T2W_mlb{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    EPI2T2W_mlb{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [1 1];    %% 1 1
    EPI2T2W_mlb{1}.spm.spatial.coreg.estwrite.roptions.interp = 2;      %% 2
    EPI2T2W_mlb{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
    EPI2T2W_mlb{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
    EPI2T2W_mlb{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'c';
end