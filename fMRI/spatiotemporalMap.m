function F = spatiotemporalMap(data, spatial_interval, temporal_interval)
F = figure;
set(F,'position',[50 50 1400 900]);
data(data==0)=nan;
slice = size(data,3);
timepoint = size(data,4);
slicePicked = 1:spatial_interval:slice;
timePicked = 1:temporal_interval:timepoint;
width = 1/length(slicePicked);
height = 1/length(timePicked);
bottom = 1-height;
for yNum = 1:length(timePicked)
    left = 0;
    for xNum = 1:length(slicePicked)
        subplot('Position',[left bottom width height])
        imagesc(rot90(data(:,:,slicePicked(xNum),timePicked(yNum))));
        left = left+width;
        axis off
         caxis([-1 1])
    end
    bottom = bottom - height;
end
    
end