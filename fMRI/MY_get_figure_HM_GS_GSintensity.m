function F = MY_get_figure_HM_GS_GSintensity(rp,G_Signal)
%% rp : head motion parameters (6)
% G_Signal: the whole brain signal intensity after ''fmask''
% F: Figure
rp(:,4:6) = rp(:,4:6).*500;
F = figure;
%% rp
positionVector1=[0 0.8 1 0.2];
subplot('Position',positionVector1)
plot(rp);
xlim([1 length(rp)]);
ylim([min(rp(:)) max(rp(:))]);
%% Global Signal
positionVector1=[0 0.55 1 0.2];
subplot('Position',positionVector1)
Global_Signal = mean(G_Signal);
plot(Global_Signal);
xlim([1 length(rp)]);
ylim([min(Global_Signal(:)) max(Global_Signal(:))*1.2]);
%% whole brain signal intensity
positionVector1=[0 0 1 0.5];
subplot('Position',positionVector1)
imagesc(G_Signal), caxis([mean(G_Signal(:))-2*std(G_Signal(:)) mean(G_Signal(:))+2*std(G_Signal(:))]);

end