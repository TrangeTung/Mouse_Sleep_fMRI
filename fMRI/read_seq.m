function [Img]=read_seq(path,pars)

% FUNCTION read_seq.m
% Reads and Reorients 4D Bruker volume series

        fid     =   fopen(deblank(path),'r');
        Img     =   fread(fid,pars.tp,pars.endian);
        fclose(fid);
       
        dims    =   pars.dims(1:3);
        dims    =   squeeze([dims(pars.vect); pars.dims(4)]);
        dims(3) =   numel(Img)/prod(dims)*dims(3);
        switch lower(pars.method)
            case 'fieldmap'
                Img     =   Img(1:prod(dims));
            case {'epi'}
%                 slope   =   repmat(pars.frame.slope ,[prod(dims(1:2)),1]);
%                 offset  =   repmat(pars.frame.offset,[prod(dims(1:2)),1]);
%                 Img     =   Img(:).*slope(:) + offset(:);
                % divide 20 for compresing the data to int 16
%                 Img     =   Img/20;
        end
        Img     =   reshape(Img,dims');   
        
end
