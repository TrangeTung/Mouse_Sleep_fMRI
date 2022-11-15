function calculateVDM_mlb=MY_get_default_calculateVDM_batch_struct(phasemap,magmap,epi_first_frame,EffectiveTE,tert)

    calculateVDM_mlb{1,1}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.phase = phasemap;
    calculateVDM_mlb{1,1}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.magnitude = magmap;
    calculateVDM_mlb{1,1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.et = EffectiveTE ;
    calculateVDM_mlb{1,1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.maskbrain = 0;
    calculateVDM_mlb{1,1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.blipdir = 1;
    calculateVDM_mlb{1,1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.tert = tert;
    calculateVDM_mlb{1,1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.epifm = 0;
    calculateVDM_mlb{1,1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.ajm = 0;
    calculateVDM_mlb{1,1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.method = 'Mark3D';
    calculateVDM_mlb{1,1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.fwhm = 10;
    calculateVDM_mlb{1,1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.pad = 0;
    calculateVDM_mlb{1,1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.ws = 1;
    calculateVDM_mlb{1,1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.template = epi_first_frame;
    calculateVDM_mlb{1,1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.fwhm = 5;
    calculateVDM_mlb{1,1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.nerode = 2;
    calculateVDM_mlb{1,1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.ndilate = 4;
    calculateVDM_mlb{1,1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.thresh = 0.5;
    calculateVDM_mlb{1,1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.reg = 0.02;
    calculateVDM_mlb{1,1}.spm.tools.fieldmap.calculatevdm.subj.session.epi = epi_first_frame;
    calculateVDM_mlb{1,1}.spm.tools.fieldmap.calculatevdm.subj.matchvdm = 1;
    calculateVDM_mlb{1,1}.spm.tools.fieldmap.calculatevdm.subj.sessname = 'session';
    calculateVDM_mlb{1,1}.spm.tools.fieldmap.calculatevdm.subj.writeunwarped = 1;
    calculateVDM_mlb{1,1}.spm.tools.fieldmap.calculatevdm.subj.anat = '';
    calculateVDM_mlb{1,1}.spm.tools.fieldmap.calculatevdm.subj.matchanat = 0;
end