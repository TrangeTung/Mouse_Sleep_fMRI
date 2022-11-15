function SliceTiming_mlb = MY_get_default_slicetiming_batch_struct(all_func,Nslice,TR,sliceorder)
    SliceTiming_mlb{1}.spm.temporal.st.scans = all_func;
    SliceTiming_mlb{1}.spm.temporal.st.nslices = Nslice;
    SliceTiming_mlb{1}.spm.temporal.st.tr = TR;
    SliceTiming_mlb{1}.spm.temporal.st.ta = TR-TR/Nslice;
    switch sliceorder
        case 'IA'
          SliceTiming_mlb{1}.spm.temporal.st.so = [1:2:Nslice 2:2:Nslice];
        case 'IB'
          SliceTiming_mlb{1}.spm.temporal.st.so = [Nslice:-2:1, Nslice-1:-2:1];
        case 'A'
          SliceTiming_mlb{1}.spm.temporal.st.so = [1:1:Nslice];
        case 'D'
          SliceTiming_mlb{1}.spm.temporal.st.so = [Nslice:-1:1];
    end
    SliceTiming_mlb{1}.spm.temporal.st.refslice = 1;
    SliceTiming_mlb{1}.spm.temporal.st.prefix = 's';
end