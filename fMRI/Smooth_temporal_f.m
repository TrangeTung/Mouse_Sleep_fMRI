function raw_fil=Smooth_temporal_f(raw,raw_mask, TR)
%h = fspecial('gaussian',[4,4],0.5);
    %wbhi = waitbar(0, 'Temporally Smoothing Session 1');
    img_fmask = fmask(raw,raw_mask);
    clear raw;
    [point,T] = size(img_fmask);
    raw_filtered = zeros(point,T);
    for j = 1:point
        temp = squeeze(img_fmask(j,:));
        raw_filtered(j,:) =  Rec_Filter(temp, TR,0.5,0.0005);% 0.1,1/T/3);
    end
    clear img_fmask;
    raw_fil = funmask(raw_filtered,raw_mask);

end

function New = Rec_Filter(Original, TR, High, Low)
    scan_num = length(Original);
    Freq_res = 1/TR/scan_num;
    Fre_up_lim = floor(High/Freq_res);
    Fre_low_lim = ceil(Low/Freq_res);
    temp = fftshift(fft(Original));
    temp(round(scan_num/2)-(Fre_low_lim-2):round(scan_num/2)+(Fre_low_lim-1)) = 0+0i;
    temp(1:round(scan_num/2)-Fre_up_lim) = 0+0i;
    temp(round(scan_num/2)+Fre_up_lim+1:scan_num) = 0+0i;
    New = abs(ifft(temp));

end

function out=fmask(fdata,mask)
temp=squeeze(fdata(:,:,:,1));
if ~isequal(size(temp), size(mask))
    error('mask size is not equal to data size');
end
idx=find(mask>0);
fdata_reshape=reshape(fdata,[],size(fdata,4));
out=fdata_reshape(idx,:);
end


function out=funmask(fdata,mask)
out=zeros(size(mask,1),size(mask,2),size(mask,3),size(fdata,2));
out_reshape=reshape(out,[],size(fdata,2));
out_reshape(mask>0,:)=fdata;
out=reshape(out_reshape,size(mask,1),size(mask,2),size(mask,3),size(fdata,2));
end
