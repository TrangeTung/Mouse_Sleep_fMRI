function [Y,XYZ] = spm_read_vols_4D(V,mask)
% Read in entire image volumes
% FORMAT [Y,XYZ] = spm_read_vols(V,mask)
% V    - vector of mapped image volumes to read in (from spm_vol)
% mask - implicit zero mask?
%
% Y    - 4D matrix of image data, fourth dimension indexes images
% XYZ  - 3xn matrix of XYZ locations returned (in mm)
%__________________________________________________________________________
%
% For image data types without a representation of NaN (see spm_type),
% implicit zero masking can be used. If mask is set, then zeros are
% treated as masked, and returned as NaN.
%__________________________________________________________________________
% Copyright (C) 1999-2013 Wellcome Trust Centre for Neuroimaging

% Andrew Holmes
% $Id: spm_read_vols.m 5731 2013-11-04 18:11:44Z guillaume $

close all
%-Argument checks
%--------------------------------------------------------------------------
if nargin<2, mask = 0; end
if nargin<1, error('insufficient arguments'), end

spm_check_orientations(V);

%      datatype__________________________________________________________________
%		2 - uint8,  4 - int16,  8 - int32,  16 - float32,
%		32 - complex64,  64 - float64,  128 - RGB24,
%		256 - int8,  511 - RGB96,  512 - uint16,
%		768 - uint32,  1792 - complex128

p=V(1).dt(1);
switch p(1)
   case  2;  precision = 'uint8';
   case  4;  precision = 'int16';
   case  8;  precision = 'int32';
   case  16;  precision = 'float32';
   case  32;  precision = 'complex64';
   case  64;  precision = 'float64';
   case  128;  precision = 'RGB24';
   case  256;  precision = 'int8';
   case  512;  precision = 'uint16';
   case  511;  precision = 'RGB96';
   case  768;  precision = 'uint32';
   case  1792;  precision = 'complex128';
end


%-Read in image data
%--------------------------------------------------------------------------
n = numel(V);                       %-#images
l = ones(V(1).dim(1:3));            %-wholebrain mask
s = prod([V(1).dim(1:3),n]);

fid = fopen(V(1).fname,'r');
offset = V(1).pinfo(end);
fseek(fid,offset, 'bof');
k = fread(fid, s,sprintf('*%s', precision));
fclose(fid); 

r=reshape(k,[s/n,n]); clear k

Y=zeros([V(1).dim(1:3),n]);
Yr=reshape(Y,[],size(r,2));
Yr(l>0,:)=r;
Y=reshape(Yr,[V(1).dim(1:3),n]);

clear Yr r l s
%-Apply implicit zero mask for image datatypes without a NaNrep
%--------------------------------------------------------------------------
if mask
    %-Work out images without NaNrep
    im = false(n,1);
    for i=1:n, im(i)=~spm_type(V(i).dt(1),'NaNrep'); end
    %-Mask
    Y(Y(:,:,:,im)==0) = NaN;
end

%-Return as 3D matrix if single image
%--------------------------------------------------------------------------
if n==1, Y=Y(:,:,:,1); end

%-Compute XYZ co-ordinates (if required)
%--------------------------------------------------------------------------
if nargout>1
    [R,C,P]  = ndgrid(1:V(1).dim(1),1:V(1).dim(2),1:V(1).dim(3));
    RCP      = [R(:)';C(:)';P(:)'];
    clear R C P
    RCP(4,:) = 1;
    XYZ      = V(1).mat(1:3,:)*RCP;
end


end
