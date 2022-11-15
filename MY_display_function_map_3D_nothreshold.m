
function Img_RGB = MY_display_function_map_3D_nothreshold(Func_Img_3D,bar_value,slice)

map_nothre_reshape = reshape(permute(flip(Func_Img_3D(:,:,slice),1),[2 1 3]),[size(Func_Img_3D,2) size(Func_Img_3D,1)*numel(slice)]);
% defaultMap = MY_viridis_colormap(64)/255;
% defaultMap = jet(64);
gdmap = [(0:31)'/31,(0:31)'/31,ones(32,1)];
drmap = [ones(32,1),(0:31)'/31,(0:31)'/31];
defaultMap = [gdmap;flipud(drmap);];
%defaultMap = jet(64);
%close(fig);
nothrebar = [min(bar_value(:)) max(bar_value(:))];
tmap = ones([size(map_nothre_reshape),3]);
map_normalize = (map_nothre_reshape-nothrebar(1))/(nothrebar(2)-nothrebar(1));
map_normalize(map_normalize<=0) = 0.000000001;
map_normalize(map_normalize>=1) = 1;
nan_mask = isnan(map_normalize);
map_normalize(isnan(map_normalize)) = 0.000000001;
tmap(:,:,1) = reshape(defaultMap(ceil(map_normalize*64),1)*255,size(map_normalize));
tmap(:,:,2) = reshape(defaultMap(ceil(map_normalize*64),2)*255,size(map_normalize));
tmap(:,:,3) = reshape(defaultMap(ceil(map_normalize*64),3)*255,size(map_normalize));

mask_index = find(double(nan_mask)==1);
tmap(mask_index+numel(map_nothre_reshape)*0) = 0;
tmap(mask_index+numel(map_nothre_reshape)*1) = 0;
tmap(mask_index+numel(map_nothre_reshape)*2) = 0;


tmap = flip(tmap,1);
Img_RGB = tmap;
end