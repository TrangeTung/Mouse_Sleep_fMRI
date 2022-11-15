function smooth_mlb=MY_get_default_smooth_batch_struct(img,smooth_para)
    smooth_mlb{1}.spm.spatial.smooth.data = img;
    smooth_mlb{1}.spm.spatial.smooth.fwhm = smooth_para;
    smooth_mlb{1}.spm.spatial.smooth.dtype = 0;
    smooth_mlb{1}.spm.spatial.smooth.im = 0;
    smooth_mlb{1}.spm.spatial.smooth.prefix = 's';
end