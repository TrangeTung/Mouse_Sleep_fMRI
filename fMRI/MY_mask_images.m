function MY_mask_images(path,folder,file_name,mask,output_file_name,file_type)


for kk = folder
    switch file_type
        case 'EPI'
            file_path = [path,'Results\',num2str(kk),'\'];
            head = spm_vol([file_path,file_name]);
            img_RAM = fmask(spm_read_vols(head),mask);
            %% detrend
%             RegBasFuc = [ones(size(img_RAM,2),1) (1:size(img_RAM,2))'];
%             img_RAM(isnan(img_RAM)) = 0;
%             for iii = 1:size(img_RAM,1)
%                 [Beta,~,Residual] = regress(squeeze(img_RAM(iii,:))',RegBasFuc);
%                 img_RAM(iii,:) = Residual + Beta(1);
%             end
            img = funmask(img_RAM,mask);
            for i = 1:numel(head)
               head(i).fname = [file_path,output_file_name]; 
            end
            spm_write_vol_4D(head,img);
 
        case 'Fieldmap'
            for jj = 2:-1:1
                file_path = [path,'Results\FieldMap',num2str(jj),'\'];
                head = spm_vol([file_path,file_name]);
                if jj == 1; mask=ones(size(mask)); end
                img = spm_read_vols(head).*mask;
                head.fname = [file_path,output_file_name]; 
                spm_write_vol_4D(head,img);
            end
        case 'T2'
            file_path = [path,'Results\T2\'];
            head = spm_vol([file_path,file_name]);
            if numel(head)>1;head(2:end) = [];end
            img = spm_read_vols(head).*mask;
            head.fname = [file_path,output_file_name]; 
            spm_write_vol_4D(head,img);
    end
    
end
end

