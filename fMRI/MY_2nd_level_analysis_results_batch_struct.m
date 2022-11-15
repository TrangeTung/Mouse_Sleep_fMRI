function results_mlb = MY_2nd_level_analysis_results_batch_struct(SPM_filepath,name,weights)

results_mlb{1}.spm.stats.con.spmmat = SPM_filepath;
results_mlb{1}.spm.stats.con.consess{1}.tcon.name = name;
results_mlb{1}.spm.stats.con.consess{1}.tcon.weights = weights;
results_mlb{1}.spm.stats.con.consess{1}.tcon.sessrep = 'replsc';
results_mlb{1}.spm.stats.con.delete = 0;

end