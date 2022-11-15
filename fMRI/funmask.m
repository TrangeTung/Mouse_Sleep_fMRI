function out=funmask(fdata,mask)
out=zeros(size(mask,1),size(mask,2),size(mask,3),size(fdata,2));
out_reshape=reshape(out,[],size(fdata,2));
out_reshape(mask>0,:)=fdata;
out=reshape(out_reshape,size(mask,1),size(mask,2),size(mask,3),size(fdata,2));
end
