function out=fmask(fdata,mask)
temp=squeeze(fdata(:,:,:,1));
if ~isequal(size(temp), size(mask))
    error('mask size is not equal to data size');
end
idx=find(mask>0);
fdata_reshape=reshape(fdata,[],size(fdata,4));
out=fdata_reshape(idx,:);
end