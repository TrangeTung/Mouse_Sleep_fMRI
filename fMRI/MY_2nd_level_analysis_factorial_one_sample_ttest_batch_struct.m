function one_sample_ttest_mlb = MY_2nd_level_analysis_factorial_one_sample_ttest_batch_struct(dir,scans)
%% dir,scans cell

one_sample_ttest_mlb{1}.spm.stats.factorial_design.dir = dir;
one_sample_ttest_mlb{1}.spm.stats.factorial_design.des.t1.scans = scans;
one_sample_ttest_mlb{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
one_sample_ttest_mlb{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
one_sample_ttest_mlb{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
one_sample_ttest_mlb{1}.spm.stats.factorial_design.masking.im = 1;
one_sample_ttest_mlb{1}.spm.stats.factorial_design.masking.em = {''};
one_sample_ttest_mlb{1}.spm.stats.factorial_design.globalc.g_omit = 1;
one_sample_ttest_mlb{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
one_sample_ttest_mlb{1}.spm.stats.factorial_design.globalm.glonorm = 1;

end
