function Bruker2nifti_multislice(Main_path,folder,species)

% FUNCTION Bruker2nifti.m
% Extracts scanner parameters and builds the Nifti volume


for kk = folder

    pars   =   get_pars([Main_path,num2str(kk),'\pdata\1\2dseq']);
    switch lower(pars.method)
        case {'fieldmap'}
            subfolder = 1:2;
        otherwise
            subfolder = 1;
    end
    for jj = subfolder
        try
            path = [Main_path,num2str(kk),'\pdata\',num2str(jj),'\2dseq1'];
            pars            =   get_pars(path);
            Img             =   read_seq(path,pars); 
        catch
            path = [Main_path,num2str(kk),'\pdata\',num2str(jj),'\2dseq'];
            pars            =   get_pars(path);
            Img             =   read_seq(path,pars); 
        end

        switch lower(pars.method)
            case {'rare'}
                dt = 4;
                outpath = [Main_path,'Results\T2'];
            case {'fieldmap'}
                dt = 16;
                outpath = [Main_path,'Results\FieldMap',num2str(jj)];
            otherwise
                dt = 4;
                outpath = [Main_path,'Results\',num2str(kk)];
        end
        
        
        switch lower(species)
            case {'mouse'};      Voxel_Multiple = 20;
            case {'rat'};        Voxel_Multiple = 10;
            case {'marmoset'};   Voxel_Multiple = 6;
            case {'monkey'};     Voxel_Multiple = 2;
            case {'human'};      Voxel_Multiple = 1;
        end
        pars.resol = pars.resol * Voxel_Multiple;
        pars.pos0  = pars.pos0  * Voxel_Multiple;
        
        % Converts Vol to Nifti format with the help of all the image parameters
        % (or=orientation, r_out=readout, idist=scanner isodist, m_or=orientation
        % matrix, dims= dimensions, FOV=field of view, resol=resolution,
        % offset=scanner offset, tp= data type)
        
        or      =   pars.orient;
        orient  =   pars.m_or;
        dims    =   pars.dims;
        resol   =   pars.resol;   %in mm
        tp      =   pars.tp;
        pos0    =   pars.pos0;
        vect    =   pars.vect;
        
        
        
        %      datatype__________________________________________________________________
        %		2 - uint8,  4 - int16,  8 - int32,  16 - float32,
        %		32 - complex64,  64 - float64,  128 - RGB24,
        %		256 - int8,  511 - RGB96,  512 - uint16,
        %		768 - uint32,  1792 - complex128
        
        
        dims        =   dims(vect);
        
        %MATRIX________________________________________________________________
        orig = pos0./resol;
        off  = double(orig(1:3)).*double(resol(1:3));
        mat  = [ resol(1)     0         0      off(1)
            0      resol(2)     0      off(2)
            0         0      resol(3)  off(3)
            0         0         0        1   ];
        
        if strcmp(or,'axial')
            mat             =   spm_matrix([0 0 0 pi/2 0 0 1 1 1])*mat;
        elseif strcmp(or,'sagittal')
            mat             =   spm_matrix([0 0 0 0 -pi/2 pi 1 1 1])*mat;
        elseif strcmp (or,'coronal')
            mat             =   spm_matrix([0 0 0 0 0 pi 1 1 1])*mat;
        end
        
        
        % Write
        if ~exist(outpath,'dir'), mkdir(outpath);end
        path_out=[outpath '\' '2dseq.nii'];
        for k = 1:size(Img,4)
            Vol(k,1) = struct(  'fname',    path_out,...
                                'dim',      size(Img(:,:,:,k)),...
                                'mat',      mat,...
                                'n',        [k,1],...
                                'pinfo',    [1;0;0],...
                                'descrip',  ['fm_',species,' image'],...
                                'dt',       [dt 0]);
        end
        spm_write_vol_4D(Vol,Img);
        
        
        cd(outpath)
        ihdr = Vol(1);
        ihdr.fname = [outpath '\' 'mean2dseq.nii'];
        img = mean(Img,4);
        spm_write_vol(ihdr,img);
        
    end
end
end


