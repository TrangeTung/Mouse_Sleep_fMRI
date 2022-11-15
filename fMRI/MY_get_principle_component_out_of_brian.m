function coeff = MY_get_principle_component_out_of_brian(fMRI_data,lmask,pc_num)
%% fMRI_data: The raw data after Bruker2nifiti or 4D data include the signal out of brain.
% lmask: the mask for fMRI_data.
% ps_num: the number of principle components we need.

% stage 01 : signal intensity normalization (out of the brain)
lmask = abs(1-lmask);
SE = strel('square',3);
mask_RAM = reshape(lmask,[size(lmask,1) size(lmask,2)*size(lmask,3)]);
mask_RAM =  imerode(mask_RAM,SE);
mask_RAM = reshape(mask_RAM,size(lmask));
fMRI_noise = fmask(fMRI_data,mask_RAM);
fMRI_noise(isnan(fMRI_noise)) = 0;
fMRI_intensity = mean(fMRI_noise,2);
fMRI_intensity_descend = sort(fMRI_intensity,'descend');
lmask = zeros(size(lmask));
lmask(fMRI_intensity>=fMRI_intensity_descend(floor(numel(fMRI_intensity)*0.005))) = 1;
fMRI_noise = fmask(fMRI_data,lmask);
fMRI_noise_ready_PCA = (fMRI_noise-mean(fMRI_noise,2))./std(fMRI_noise,0,2);
% stage 02 : principle component analysis
[coeff,score,latent,~]= pca(fMRI_noise_ready_PCA','algorithm','svd');

% latent_test = sort(latent,'descend');
% if sum(latent_test-latent) ~= 0
%     error('Something was wrong during PCA');
% else
%    % PC = score;%(:,1:pc_num);
% end

end