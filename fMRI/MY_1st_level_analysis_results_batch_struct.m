function results_mlb = MY_1st_level_analysis_results_batch_struct(SPM_filepath,name,weights,sessrep)

results_mlb{1}.spm.stats.con.spmmat = {SPM_filepath};
results_mlb{1}.spm.stats.con.consess{1}.tcon.name = name;
results_mlb{1}.spm.stats.con.consess{1}.tcon.weights = weights{1};
results_mlb{1}.spm.stats.con.consess{1}.tcon.sessrep = sessrep;
results_mlb{1}.spm.stats.con.consess{2}.tcon.name = name;
results_mlb{1}.spm.stats.con.consess{2}.tcon.weights = weights{2};
results_mlb{1}.spm.stats.con.consess{2}.tcon.sessrep = sessrep;
results_mlb{1}.spm.stats.con.consess{3}.tcon.name = name;
results_mlb{1}.spm.stats.con.consess{3}.tcon.weights = weights{3};
results_mlb{1}.spm.stats.con.consess{3}.tcon.sessrep = sessrep;
results_mlb{1}.spm.stats.con.delete = 0;
results_mlb{1}.spm.job.bases = 'mouse';
end