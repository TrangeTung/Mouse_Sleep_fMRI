function V = spm_write_vol_4D(V,Y)
% Write an image volume to disk, setting scales and offsets as appropriate
% FORMAT V = spm_write_vol(V,Y)
% V (input)  - a structure containing image volume information (see spm_vol)
% Y          - a one, two or three dimensional matrix containing the image voxels
% V (output) - data structure after modification for writing.
%
% Note that if there is no 'pinfo' field, then SPM will figure out the
% max and min values from the data and use these to automatically determine
% scalefactors.  If 'pinfo' exists, then the scalefactor in this is used.
%__________________________________________________________________________
% Copyright (C) 1999-2013 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_write_vol.m 5731 2013-11-04 18:11:44Z guillaume $


use_offset = false;

if ndims(Y) == 3
    
    [V,Y] = MY_design_hdr_img(V,Y);
    
    %-Create and write image
    %--------------------------------------------------------------------------
    V = spm_create_vol(V);
    V = spm_write_plane_4D(V,Y,':');
else
    
%     [V1,~] = MY_design_hdr_img(V(1),Y(:,:,:,1));
%     for n = 1:numel(V)
%         V(n).mat = V1.mat;        
%     end
    %-Create and write image
    %--------------------------------------------------------------------------
    V = spm_create_vol(V);
    V = spm_write_plane_4D(V,Y,':');
end

end

function [V,Y] = MY_design_hdr_img(V,Y)
    
    dim = [size(Y) 1 1 1];
    if ~all(dim(1:3) == V.dim(1:3))
        error('Incompatible dimensions.');
    end
    
    if ~isfield(V,'pinfo')
        V.pinfo = [1;0;0];
        rescal  = true;
    elseif ~all(isfinite(V.pinfo(1:2))) || V.pinfo(1) == 0
        V.pinfo(1:2) = [1;0];
        rescal  = true;
    else
        rescal  = false;
    end
    
    if rescal
        % Set scalefactors and offsets
        %----------------------------------------------------------------------
        dt           = V.dt(1);
        s            = find(dt == [2 4 8 256 512 768]);
        if isempty(s)
            V.pinfo(1:2) = [1;0];
        else
            dmnmx        = [0 -2^15 -2^31 -2^7 0 0 ; 2^8-1 2^15-1 2^31-1 2^7-1 2^16-1 2^32-1];
            dmnmx        = dmnmx(:,s);
            mxs          = zeros(dim(3),1)+NaN;
            mns          = zeros(dim(3),1)+NaN;
            
            for p=1:dim(3)
                tmp      = double(Y(:,:,p));
                tmp      = tmp(isfinite(tmp));
                if ~isempty(tmp)
                    mxs(p) = max(tmp);
                    mns(p) = min(tmp);
                end
            end
            
            mx = max(mxs(isfinite(mxs)));
            mn = min(mns(isfinite(mns)));
            if isempty(mx), mx = 0; end
            if isempty(mn), mn = 0; end
            if mx ~= mn
                if use_offset
                    V.pinfo(1,1) = (mx-mn)/(dmnmx(2)-dmnmx(1));
                    V.pinfo(2,1) = (dmnmx(2)*mn-dmnmx(1)*mx)/(dmnmx(2)-dmnmx(1));
                else
                    if dmnmx(1) < 0
                        V.pinfo(1) = max(mx/dmnmx(2),mn/dmnmx(1));
                    else
                        V.pinfo(1) = mx/dmnmx(2);
                    end
                    V.pinfo(2) = 0;
                end
            else
                V.pinfo(1,1) = mx/dmnmx(2);
                V.pinfo(2,1) = 0;
            end
        end
    end

end

function V = spm_write_plane_4D(V,dat,n)
% Write transverse plane(s) of image data
% FORMAT V = spm_write_plane(V,dat,n)
% V   - data structure containing image information (see spm_vol)
% dat - the two/three dimensional image to write
% n   - the plane number(s) (beginning from 1). If an entire volume
%       should be written, n should contain the single character ':'
%       instead of plane numbers.
%
% V   - (possibly) modified data structure containing image information.
%       It is possible that future versions of spm_write_plane may
%       modify scalefactors (for example).
%_______________________________________________________________________
% Copyright (C) 1999-2014 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_write_plane.m 6079 2014-06-30 18:25:37Z spm $

% For performance reasons, on network filesystems one should write
% out as large contiguous blocks data as possible at once. Therefore,
% multiple planes or even entire volumes should be handled here.
% Dimension checking is left to mat2file.
if numel(size(dat)) == 3
    if isfield(V,'n')
        n1 = num2cell(V(1).n);
        n  = {n n1{:}};
    else
        n  = {n};
    end
    S      = struct('type','()','subs',{{':',':',n{:}}});
else
    for k = 1:size(dat,4)
        n1 = num2cell(V(k).n);
        N  = {n n1{:}};
        S(k,:)      = struct('type','()','subs',{{':',':',N{:}}});
    end
end
for i = 1:size(dat,4)
    V(i).private.dat = subsasgn(squeeze(V(i).private.dat),S(i),squeeze(dat(:,:,:,i)));
end
end
