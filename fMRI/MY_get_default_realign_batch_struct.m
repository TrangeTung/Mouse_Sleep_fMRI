function realign_mlb=MY_get_default_realign_batch_struct(all_func)
    realign_mlb{1}.spm.spatial.realign.estwrite.data = all_func;
    realign_mlb{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
    realign_mlb{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
    realign_mlb{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;     %
    realign_mlb{1}.spm.spatial.realign.estwrite.eoptions.rtm = 0;
    realign_mlb{1}.spm.spatial.realign.estwrite.eoptions.interp = 4;
    realign_mlb{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
    realign_mlb{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
    realign_mlb{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
    realign_mlb{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
    realign_mlb{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
    realign_mlb{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
    realign_mlb{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
end