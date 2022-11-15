%clear all;
% [EPIfilename, EPIpathname, filterindex] = uigetfile( ...
%        {'*.spr','spr-files (*.spr)'; ...
%         '*.sdt','sdt-files (*.sdt)'; ...
%         '*.*',  'All Files (*.*)'}, ...
%         'Select mc data for the subject', ...
%         'MultiSelect', 'on');
function raw_fil=Smooth_temporal_f(raw,raw_mask, TR)
%h = fspecial('gaussian',[4,4],0.5);
%TR = 1;
% if ~iscell(EPIfilename)
    %[raw interx intery interz] = loadSdt([EPIpathname,EPIfilename]);

    %wbhi = waitbar(0, 'Temporally Smoothing Session 1');
    img_fmask = fmask(raw,raw_mask);
    clear raw;
    [point,T] = size(img_fmask);
    raw_filtered = zeros(point,T);
    for j = 1:point
        temp = squeeze(img_fmask(j,:));
        raw_filtered(j,:) =  Rec_Filter(temp, TR, 0.1,1/T/3);
    end
    clear img_fmask;
    raw_fil = funmask(raw_filtered,raw_mask);
    
%     for x = 1:ImgX
%         for y = 1:ImgY
%             for z = 1:ImgZ
%                 temp = squeeze(raw(x,y,z,:));
%                 raw_fil(x,y,z,:) =  Rec_Filter(temp, TR, 0.1,1/T/3);
%             end
%         end
%      %   waitbar(x/ImgX);
%     end
    %close(wbhi);
    %saveSdt(raw_fil, [EPIpathname, strrep(EPIfilename,'.spr',''),'_tem'],interx, intery, interz);
% else
%     runnum = size(EPIfilename, 2);
%     for r = 1:runnum
%         [raw interx intery interz] = loadSdt([EPIpathname,EPIfilename{:,r}]);
%         [ImgX ImgY ImgZ T] = size(raw);
%         wbhi = waitbar(0, ['Temporally Smoothing Session ',num2str(r),'...']);
%         for x = 1:ImgX
%             for y = 1:ImgY
%                 for z = 1:ImgZ
%                     temp = squeeze(raw(x,y,z,:));
%                     raw_fil(x,y,z,:) =  Rec_Filter(temp, TR, 0.1, 1/T/3);
%                 end
%             end
%             waitbar(x/ImgX);
%         end
%         close(wbhi);
%         saveSdt(raw_fil, [EPIpathname, strrep(EPIfilename{:,r},'.spr',''),'_tem'],interx, intery, interz);
%         clear raw; clear raw_fil;
%     end
end