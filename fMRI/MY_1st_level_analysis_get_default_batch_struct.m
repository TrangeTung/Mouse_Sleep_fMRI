function first_level_analysis_mlb = MY_1st_level_analysis_get_default_batch_struct(output_dir,design_unit,all_epi,name,onset,duration,EPI_TR,multi_regress)

first_level_analysis_mlb{1}.spm.stats.fmri_spec.dir = output_dir; 
first_level_analysis_mlb{1}.spm.stats.fmri_spec.timing.units = design_unit; 
first_level_analysis_mlb{1}.spm.stats.fmri_spec.timing.RT = EPI_TR;
first_level_analysis_mlb{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
first_level_analysis_mlb{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
first_level_analysis_mlb{1}.spm.stats.fmri_spec.sess.scans = all_epi; 
first_level_analysis_mlb{1}.spm.stats.fmri_spec.sess.cond.name = name;   
first_level_analysis_mlb{1}.spm.stats.fmri_spec.sess.cond.onset = onset; 
first_level_analysis_mlb{1}.spm.stats.fmri_spec.sess.cond.duration = duration;  
first_level_analysis_mlb{1}.spm.stats.fmri_spec.sess.cond.tmod = 0;
first_level_analysis_mlb{1}.spm.stats.fmri_spec.sess.cond.pmod = struct('name', {}, 'param', {}, 'poly', {});
first_level_analysis_mlb{1}.spm.stats.fmri_spec.sess.cond.orth = 1;
first_level_analysis_mlb{1}.spm.stats.fmri_spec.sess.multi = {''};
first_level_analysis_mlb{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
first_level_analysis_mlb{1}.spm.stats.fmri_spec.sess.multi_reg = {multi_regress}; 
first_level_analysis_mlb{1}.spm.stats.fmri_spec.sess.hpf = 128;
first_level_analysis_mlb{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
first_level_analysis_mlb{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
first_level_analysis_mlb{1}.spm.stats.fmri_spec.volt = 1;
first_level_analysis_mlb{1}.spm.stats.fmri_spec.global = 'None';
first_level_analysis_mlb{1}.spm.stats.fmri_spec.mthresh = 0.8;
first_level_analysis_mlb{1}.spm.stats.fmri_spec.mask = {''};
first_level_analysis_mlb{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
end