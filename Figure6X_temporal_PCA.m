clc;clear
codepath = 'H:\ECoG\code_Trange\';
addpath(genpath(codepath));
cd(codepath)

FilePath = 'I:\raw_data_yyl\';
ihdr = spm_vol(fullfile('H:\ECoG\pca\','mask_RAM.nii'));
lmask = spm_read_vols(ihdr);

WholePath = 'H:\ECoG\pca\';
SCORE = [];
for idx=1:46
    load(fullfile(WholePath,num2str(idx,'%02d'),'pca_t.mat'));
    SCORE = cat(1,SCORE,score);
end
[Gcoeff, Gscore, ~, ~, explained, ~] = pca(SCORE','NumComponents',100);
save(fullfile(WholePath,'GroupWeight.mat'),'Gcoeff','Gscore','-v7.3');



Vall = [];
for idx=1:46
    idx
    filename = fullfile(FilePath,num2str(idx),'sLnrrsm2dseq.nii');
    fMRI_4D = spm_read_vols_4D(spm_vol(filename));
    V = fmask(fMRI_4D,lmask);
    V = (V-mean(V,2))./std(V,0,2);
    V(isnan(V))=0;
    clear fMRI_4D
    
    [coeff1, score1] = pca(V','NumComponents',1000);
    cd(fullfile('H:\ECoG\pca\',num2str(idx,'%02d')));
    save('coeff_1000.mat','coeff');
    save('score_1000.mat','score');

end
