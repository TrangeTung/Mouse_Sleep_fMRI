function EPI2Template_mlb=MY_get_default_oldnormalize_batch_struct(ref, source, other)

     V = spm_vol(ref{:});
     M = V.mat;
     re_voxel = sqrt(sum(M(1:3,1:3).^2));

    EPI2Template_mlb{1}.spm.tools.oldnorm.estwrite.subj.source = source;
    EPI2Template_mlb{1}.spm.tools.oldnorm.estwrite.subj.wtsrc = '';
    EPI2Template_mlb{1}.spm.tools.oldnorm.estwrite.subj.resample = other;
    EPI2Template_mlb{1}.spm.tools.oldnorm.estwrite.eoptions.weight = '';
    EPI2Template_mlb{1}.spm.tools.oldnorm.estwrite.eoptions.smosrc = 8;     %  1
    EPI2Template_mlb{1}.spm.tools.oldnorm.estwrite.eoptions.smoref = 8;     %
    EPI2Template_mlb{1}.spm.tools.oldnorm.estwrite.eoptions.template = ref;
    EPI2Template_mlb{1}.spm.tools.oldnorm.estwrite.eoptions.regtype = 'mni';
    EPI2Template_mlb{1}.spm.tools.oldnorm.estwrite.eoptions.cutoff = 25;
    EPI2Template_mlb{1}.spm.tools.oldnorm.estwrite.eoptions.nits = 16;          %
    EPI2Template_mlb{1}.spm.tools.oldnorm.estwrite.eoptions.reg = 10;
    EPI2Template_mlb{1}.spm.tools.oldnorm.estwrite.roptions.bb = [NaN,NaN,NaN;NaN,NaN,NaN];
    EPI2Template_mlb{1}.spm.tools.oldnorm.estwrite.roptions.vox = re_voxel(:);
    EPI2Template_mlb{1}.spm.tools.oldnorm.estwrite.roptions.interp = 2;     %
    EPI2Template_mlb{1}.spm.tools.oldnorm.estwrite.roptions.wrap = [1 1 0];   %
    EPI2Template_mlb{1}.spm.tools.oldnorm.estwrite.roptions.prefix = 'n';
    EPI2Template_mlb{1}.spm.tools.oldnorm.estwrite.roptions.preserve = 0;
end