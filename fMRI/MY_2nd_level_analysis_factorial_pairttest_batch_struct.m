function pair_t_test_mlb = MY_2nd_level_analysis_factorial_pairttest_batch_struct(dir,scans)
%% dir,scans cell

pair_t_test_mlb{1}.spm.stats.factorial_design.dir = dir;
pair_t_test_mlb{1}.spm.stats.factorial_design.des.pt.pair = scans;
pair_t_test_mlb{1}.spm.stats.factorial_design.des.pt.gmsca = 0;
pair_t_test_mlb{1}.spm.stats.factorial_design.des.pt.ancova = 0;
pair_t_test_mlb{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
pair_t_test_mlb{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
pair_t_test_mlb{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
pair_t_test_mlb{1}.spm.stats.factorial_design.masking.im = 1;
pair_t_test_mlb{1}.spm.stats.factorial_design.masking.em = {''};
pair_t_test_mlb{1}.spm.stats.factorial_design.globalc.g_omit = 1;
pair_t_test_mlb{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
pair_t_test_mlb{1}.spm.stats.factorial_design.globalm.glonorm = 1;

end