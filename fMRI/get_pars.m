function pars = get_pars(path)

% FUNCTION get_pars.m
% Extracts Bruker aquisition parameters from acqp, method and reco files


    pars.acq_dim     =   0;
%     pars.scale       =   ones(4,1);    %%%%20171010 seems to have only 3
%     values
    pars.scale       =   ones(2,1);
    pars.cmpx        =   0;
    pars.n_coils     =   1;
    pars.orient      =   '';
    pars.r_out       =   '';
    pars.idist       =   0;
    pars.m_or        =   zeros(3,2);
    pars.dims        =   zeros(4,1);
    pars.FOV         =   zeros(3,1);      %in mm
    pars.resol       =   zeros(3,1);   %in mm
    pars.offset      =   zeros(3,1);   %in mm
    pars.tp          =   '';
    pars.endian      =   'l';
    pars.n_acq       =   '';
    pars.TR          =   0;
    
    
    matching    =   ['.*Di?a?([0-9]+).*\' filesep '.*'];
    if isunix 
        [a b c d e f]   =   regexpi(path,matching,'match','tokens');
    else
        [a b, c d e f]   =   regexpi(path,matching,'match','split');
    end
%     day         =   str2num(char(path(e{1,:})));
    path_acqp   =   [fileparts(fileparts(fileparts(path))) filesep 'acqp'];
%     path_acqp   =   [path, '\acqp'];
    fid         =   fopen(deblank(path_acqp), 'rt');
    
    
    %Reading method
    path_method     =   [fileparts(fileparts(fileparts(path))) filesep 'method'];
%     path_method     =   [path, '\method'];
    fid2            =   fopen(deblank(path_method), 'rt');
    while feof(fid2) == 0
        line    =   fgetl(fid2);
        tag     =   strread(line,'%s','delimiter','=');
        switch tag{1}
            case '##$PVM_SPackArrNSlices'
                pars.dims(3)     =   eval(fgetl(fid2));
            case '##$PVM_NRepetitions' 
                pars.dims(4)     =   eval(tag{2}); 
            case '##$PVM_ScanTime'
                scan_time   =   eval(tag{2}); 
            case '##$PVM_SPackArrSliceOrient'
                pars.orient  =   fgetl(fid2);
            case '##$PVM_SPackArrReadOrient'
                pars.r_out   =   fgetl(fid2);
            case '##$Method'
                tag     =	strread(line,'%s','delimiter',':');
                pars.method   =   tag{2}(1:end-1);
            case '##$PVM_SPackArrGradOrient'
                m       =   0;
                while true
                    read    =   fgetl(fid2);
                    header  =   strread(read,'%c','delimiter','#');
                    if header(1)=='#' break; end
                    m       =   [m;strread(deblank(read))'];
                end
                pars.m_or   =   m(2:10);
                pars.m_or   =   reshape(pars.m_or,[3,3]);

        end
    end
    fclose(fid2);
    if exist('scan_time','var')
        pars.TR          =   scan_time/pars.dims(4) ;
    end
    
    
    %Reading acqp
    while feof(fid) == 0
        line    =   fgetl(fid);
        tag     =   strread(line,'%s','delimiter','=');
        switch tag{1}
            case '##$ACQ_dim'
                pars.acq_dim     =   str2num(tag{2});
            case '##$ACQ_slice_sepn'
                if pars.acq_dim==2
                    pars.idist       =   eval(fgetl(fid));
                end
            case '##$ACQ_slice_thick' 
                if pars.acq_dim==2              
                    pars.resol(3)    =   eval(tag{2});
                end
            case '##$ACQ_slice_offset'
                          slices        =   [];
                          while true
                            read        =   fgetl(fid);
                            header      =   strread(read,'%c','delimiter','#');
                            if strcmp(header(1),'#') ||strcmp(header(1),'$') break; end
                            slices      =   [slices;strread(read,'%f','delimiter','\\ ')];
                            pars.offset(3)   =   mean(slices);
                          end
            case '##$ACQ_read_offset'
                x = fgetl(fid);
               try 
                   read         =   strread(x,'%f','delimiter','\\ '); 
                   pars.offset(1)    =   read(1);
               catch
                   pars.offset(1)    =   eval(x(2:end))/eval(tag{2});
               end
               
            case '##$ACQ_phase2_offset' 
                x = fgetl(fid);
               try 
                   read         =   strread(x,'%f','delimiter','\\ '); 
                   pars.offset(2)    =   read(1);  
               catch
                   pars.offset(2)    =   eval(x(2:end))/eval(tag{2});
               end
               
        end
    end
    fclose(fid);

   
    
    %Reading reco 
    [pathstr,nam,ext]   =   fileparts(fileparts(fileparts(fileparts(path))));
    pars.n_acq          =   nam;
    fid                 =   fopen(deblank([fileparts(path) filesep 'reco']),'rt');
%     fid                 =   fopen(deblank([path, '\pdata\1\reco']),'r');
    while feof(fid) == 0        
            line    =   fgetl(fid);
            tag     =	strread(line,'%s','delimiter','=');
            switch tag{1}
                case '##$RECO_fov'
                    if strcmp(tag{2},'( 2 )')
                        [a b]               =   strread(fgetl(fid),'%f %f');
                        pars.FOV(1)              =   10*a; 
                        pars.FOV(2)              =   10*b; %from cm to mm                        
                    elseif strcmp(tag{2},'( 3 )')
                        [a b c]               =   strread(fgetl(fid),'%f %f %f');
                        pars.FOV(1)              =   10*a; 
                        pars.FOV(2)              =   10*b; %from cm to mm
                        pars.FOV(3)              =   10*c;
                    end
                case '##$RECO_size'
                    if strcmp(tag{2},'( 2 )')
                        [pars.dims(1) pars.dims(2)]   =   strread(fgetl(fid),'%d %d');    
                    elseif strcmp(tag{2},'( 3 )')
                        [pars.dims(1) pars.dims(2) pars.dims(3)]   =   strread(fgetl(fid),'%d %d %d');
                        pars.idist               =   pars.resol(3);
                    end                    
                case '##$RECO_wordtype'
                    pars.tp  =   tag{2};
                    switch pars.tp
                        case '_8BIT_USGN_INT'
                            pars.tp  =   'uint8';
                        case '_16BIT_SGN_INT' 
                            pars.tp  =   'int16'; 
                        case '_32BIT_SGN_INT' 
                            pars.tp  =   'int32'; 
                        case '_32BIT_FLT' 
                            pars.tp  =   'float32';
                        case '_64BIT_FLT' 
                            pars.tp  =   'float64'; 
                        case '_8BIT_SGN_INT' 
                            pars.tp  =   'int8'; 
                        case '_16BIT_USGN_INT' 
                            pars.tp  =   'uint16'; 
                        case '_32BIT_USGN_INT' 
                            pars.tp  =   'uint32';              
                    end
                case '##$RECO_byte_order'
                    if strcmp (tag{2},'littleEndian')
                        pars.endian  =   'l';
                    elseif strcmp (tag{2},'bigEndian')
                        pars.endian  =   'b';                        
                    end
                case '##$RecoNumInputChan'
                    pars.n_coils     =   str2num(tag{2});
                case '##$RECO_image_type'
                    if strmatch(tag{2},'COMPLEX_IMAGE')
                        pars.cmpx    =   1;
                    end
                case '##$RECO_map_offset'
                      offset        =   [];
                      while true
                        read        =   fgetl(fid);
                        header      =   strread(read,'%c','delimiter','#');
                        if strcmp(header(1),'#') ||strcmp(header(1),'$');break; end
                        try
                            offset      =   [offset;strread(read,'%f','delimiter','\\ ')];
                        catch
                            offset      =   [offset;eval(read(2:end))];
                        end
                      end
                      pars.frame.offset = offset;
                case '##$RecoScaleChan'
                    [a b c d]   =   strread(fgetl(fid),'%f %f %f %f');
                    if nnz([a b c d]) == 1
                        pars.scale   =   1;
                    else
                        pars.scale(1)=   a(1); 
                        pars.scale(2)=   b(1); 
%                         pars.scale(3)=   c(1); 
%                         pars.scale(1)=   a; 
%                         pars.scale(2)=   b; 
%                         pars.scale(3)=   c; 
%                         pars.scale(4)=   d;
                    end
            end  
    end 
    fclose(fid);
    
    %Reading reco 
    [pathstr,nam,ext]   =   fileparts(fileparts(fileparts(fileparts(path))));
    pars.n_acq          =   nam;
    fid                 =   fopen(deblank([fileparts(path) filesep 'reco']),'rt');
%     fid                 =   fopen(deblank([path, '\pdata\1\reco']),'r');
    while feof(fid) == 0        
            line    =   fgetl(fid);
            tag     =	strread(line,'%s','delimiter','=');
            switch tag{1}
    
                    case '##$RECO_map_slope'
                      slope        =   [];
                      while true
                        read        =   fgetl(fid);
                        header      =   strread(read,'%c','delimiter','#');
                        if strcmp(header(1),'#') ||strcmp(header(1),'$');break; end
                        
                        try
                            slope      =   [slope;strread(read,'%f','delimiter','\\ ')];
                        catch
                            slope      =   [slope;eval(read(2:end))];
                        end
                        
                      end
                      pars.frame.slope = slope;

            end
    end
    fclose(fid);
    
    %resol
    pars.resol(1:2)     =   [pars.FOV(1)/pars.dims(1); pars.FOV(2)/pars.dims(2)];
    if pars.acq_dim==2
        pars.FOV(3)     =   (pars.dims(3)-1)*pars.idist+pars.resol(3);
    end
    pars.dims           =   cast(pars.dims,'int16');    
    
    % Position calculation. FOV, resol order
    dim             =   pars.dims;
    pars.resol(3)   =   pars.FOV(3)/double(pars.dims(3)); 
    switch pars.orient
        case 'axial' 
            if strcmp(pars.r_out,'L_R')                     
                pars.vect      =   [1,2,3];
                pars.FOV    =   pars.FOV(pars.vect);                
                half_vx     =   pars.resol(pars.vect)/2;  
                pars.offset =   [pars.offset(1);pars.offset(2);-pars.offset(3)];                 
                half_vx     =   pars.resol(pars.vect)/2;                  
                m_or2       =   pars.m_or;                
            elseif   strcmp(pars.r_out,'A_P')
                pars.vect      =   [2,1,3];
                pars.FOV    =   pars.FOV(pars.vect);                  
                half_vx     =   pars.resol(pars.vect)/2;  
                pars.offset =   [pars.offset(2);pars.offset(1);-pars.offset(3)]; 
                m_or2       =   [pars.m_or(:,2) pars.m_or(:,1) pars.m_or(:,3)];                
            end
        case 'sagittal'
            if strcmp(pars.r_out,'H_F')                 
                pars.vect      =   [2,1,3];
                pars.FOV    =   pars.FOV(pars.vect);  
                half_vx     =   pars.resol(pars.vect)/2;                 
                pars.offset =   [pars.offset(2);pars.offset(1);-pars.offset(3)];                 
                m_or2       =   [pars.m_or(:,2) pars.m_or(:,1) pars.m_or(:,3)];                
            elseif strcmp(pars.r_out,'A_P')
                pars.vect      =   [1,2,3];
                pars.FOV    =   pars.FOV(pars.vect); 
                half_vx     =   pars.resol(pars.vect)/2;                  
                pars.offset =   [pars.offset(1);pars.offset(2);-pars.offset(3)]; 
                m_or2       =   [pars.m_or(:,1) pars.m_or(:,2) pars.m_or(:,3)];                
            end
        case 'coronal'
            if strcmp(pars.r_out,'H_F')                 
                pars.vect      =   [2,1,3];
                pars.FOV    =   pars.FOV(pars.vect);        
                half_vx     =   pars.resol(pars.vect)/2;                  
                pars.offset =   [pars.offset(2);pars.offset(1);-pars.offset(3)];                 
                m_or2       =   [pars.m_or(:,2) pars.m_or(:,1) pars.m_or(:,3)];                
            elseif strcmp(pars.r_out,'L_R')
                pars.vect      =   [1,2,3];
                pars.FOV    =   pars.FOV(pars.vect);  
                half_vx     =   pars.resol(pars.vect)/2;                  
                pars.offset =   [pars.offset(1);pars.offset(2);-pars.offset(3)]; 
                m_or2       =   [pars.m_or(:,1) pars.m_or(:,2) pars.m_or(:,3)];                
            end        
    end
    shift       =   -pars.FOV/2-half_vx-pars.offset;
    pars.pos0   =   m_or2*shift;
    pars.m_or   =   m_or2;
    pars.dims(pars.dims==0) = 1;

end

