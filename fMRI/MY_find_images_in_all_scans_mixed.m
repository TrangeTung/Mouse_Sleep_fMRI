function all_img_list_cell=MY_find_images_in_all_scans_mixed(animal_root_path, prefix, suffix, frames, concatenating_type, experimenttype)
%     This function is to select images with specifc prefix and
%     sufix in all scans of a animal. This structure of the directory must
%     be like: under the "animal_root_path", there are scan folders named
%     by their scan ids with pure number in range from 1 to 99. And under
%     each scan folders, there are the images to be selected, whose names matches
%     regular expressions "prefix scan_id sufix". If the "concatenating_type"is all mixed,
%     all the scan images will be put in a cell, while if it is 'separate_cells', images of 
%     each scan will be put into one cell and there is one outer cell wrapping up all these cells.
%     The result will be returned in the form compatible with spm preprocessing functions.
    all_img_list_cell={};
    for scan_inx=1:99
        scan_id=num2str(scan_inx);
        if (experimenttype == 1)
         working_folder = fullfile(animal_root_path, scan_id,'\NoFilter');
        else
         working_folder = fullfile(animal_root_path, scan_id);    
        end
        
        if exist(working_folder,'file')
            cell_func= MY_select_file_for_SPM(working_folder, [prefix suffix], frames);
            switch concatenating_type
                case 'all_mixed',
                    all_img_list_cell=[all_img_list_cell; cell_func];
                case 'separate_cells',
                    all_img_list_cell=[all_img_list_cell, {cell_func}];
                otherwise,
                    error('concatenating_type is not correct')
            end
        end
    end
end