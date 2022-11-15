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