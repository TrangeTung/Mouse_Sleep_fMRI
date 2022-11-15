function info = MY_get_basic_information_from_hdr(head,param)


vox_info = head(1).mat(1:3,1:3);
re_voxel = (abs(vox_info(vox_info~=0)));
re_voxel(re_voxel<0.01) = [];

switch param
    case 'origin'
    info = abs(head(1).mat(1:3,4))./re_voxel;
    case 'voxel_size'
    info = re_voxel;
end
    
end