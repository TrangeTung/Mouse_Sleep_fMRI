function estimate_mlb = MY_get_default_estimate_batch_struct(SPM_file)
estimate_mlb{1}.spm.stats.fmri_est.spmmat = SPM_file;
estimate_mlb{1}.spm.stats.fmri_est.write_residuals = 0;
estimate_mlb{1}.spm.stats.fmri_est.method.Classical = 1;
end