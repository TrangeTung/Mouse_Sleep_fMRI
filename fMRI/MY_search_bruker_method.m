function data = MY_search_bruker_method(param,folder,path)
% cd(path);
% %% search parameter from bruker method file
% cd(num2str(folder));
fid = fopen([path,num2str(folder),'\method'],'rb');
temp = fread(fid, inf, 'uchar');
string = char(temp);
string = string(:)';
switch param
    case 'total_EPI_readout_time'
        keyword = '$PVM_EpiNEchoes=';
        loc = strfind(string,keyword)+numel(keyword);
        [token,~] = strtok(string(loc:loc+10),'#');
        NEchoes = str2double(token);
        keyword = '$PVM_EpiEchoSpacing=';
        loc = strfind(string,keyword)+numel(keyword);
        [token,~] = strtok(string(loc:loc+10),'#');
        EchoSpacing = str2double(token);
        data = NEchoes*EchoSpacing;
        
    case 'Effective_TE'
        keyword = '$EffectiveTE=';
        loc = strfind(string,keyword)+numel(keyword);
        data(1) = str2double(string(loc+6:loc+11));
        data(2) = str2double(string(loc+13:loc+18));
        
    case 'EPI_TR'
        keyword = '$PVM_RepetitionTime=';
        loc = strfind(string,keyword)+numel(keyword);
        [token,~] = strtok(string(loc:loc+5),'#');
        token(strfind(token,newline):end) = [];
        data = str2double(token);
        
    case 'Segments'
        keyword = '$NSegments=';
        loc = strfind(string,keyword)+numel(keyword);
        [token,~] = strtok(string(loc:loc+5),'#');
        token(strfind(token,newline):end) = [];
        data = str2double(token);
        
    case 'Nslice'    
        keyword = '$PVM_SPackArrNSlices=( 1 )';
        loc = strfind(string,keyword)+numel(keyword);
        [token,~] = strtok(string(loc:loc+5),'#');
        data = str2double(token);
        
    case 'Repetitions'
        keyword = '$PVM_NRepetitions=';
        loc = strfind(string,keyword)+numel(keyword);
        [token,~] = strtok(string(loc:loc+5),'#');
        token(strfind(token,newline):end) = [];
        data = str2double(token);
        
    case 'DummyScansDuration'
        keyword = '$PVM_DummyScansDur=';
        loc = strfind(string,keyword)+numel(keyword);
        [token,~] = strtok(string(loc:loc+5),'#');
        token(strfind(token,newline):end) = [];
        data = str2double(token);
        
        
        
end