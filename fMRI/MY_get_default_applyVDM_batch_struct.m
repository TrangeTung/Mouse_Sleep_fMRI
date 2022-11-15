function applyVDM_mlb=MY_get_default_applyVDM_batch_struct(all_epi,vdmfile,PEdirection)

    applyVDM_mlb{1,1}.spm.tools.fieldmap.applyvdm.data.scans = all_epi;
    applyVDM_mlb{1,1}.spm.tools.fieldmap.applyvdm.data.vdmfile = vdmfile;
    applyVDM_mlb{1,1}.spm.tools.fieldmap.applyvdm.roption.pedir = PEdirection;
    applyVDM_mlb{1,1}.spm.tools.fieldmap.applyvdm.roption.which = [2 1];
    applyVDM_mlb{1,1}.spm.tools.fieldmap.applyvdm.roption.rinterp = 4;
    applyVDM_mlb{1,1}.spm.tools.fieldmap.applyvdm.roption.wrap = [0 0 0];
    applyVDM_mlb{1,1}.spm.tools.fieldmap.applyvdm.roption.mask = 1;
    applyVDM_mlb{1,1}.spm.tools.fieldmap.applyvdm.roptions.prefix = 'b';
end