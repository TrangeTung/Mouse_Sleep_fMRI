function first_level_analysis_mlb = MY_1st_level_analysis_1rodentNscan_get_default_batch_struct(output_dir,design_unit,EPI_TR,sess)
first_level_analysis_mlb{1}.spm.stats.fmri_spec.dir = output_dir; 
first_level_analysis_mlb{1}.spm.stats.fmri_spec.timing.units = design_unit; 
first_level_analysis_mlb{1}.spm.stats.fmri_spec.timing.RT = EPI_TR;
first_level_analysis_mlb{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
first_level_analysis_mlb{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;

first_level_analysis_mlb{1}.spm.stats.fmri_spec.sess = sess;

first_level_analysis_mlb{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
first_level_analysis_mlb{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
first_level_analysis_mlb{1}.spm.stats.fmri_spec.volt = 1;
first_level_analysis_mlb{1}.spm.stats.fmri_spec.global = 'None';
first_level_analysis_mlb{1}.spm.stats.fmri_spec.mthresh = 0.8;
first_level_analysis_mlb{1}.spm.stats.fmri_spec.mask = {''};
first_level_analysis_mlb{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
first_level_analysis_mlb{1}.spm.stats.fmri_spec.bases.hrf.order = 'mouse';
end