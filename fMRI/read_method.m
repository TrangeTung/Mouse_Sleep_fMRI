function Method = read_method(path)
%% search parameter from bruker method file
fid = fopen(path,'rb');
temp = fread(fid, inf, 'uchar');
string = char(temp);
string = string(:)';
keyword = '$Method=<Bruker:';
loc = strfind(string,keyword)+numel(keyword);
[token,~] = strtok(string(loc:loc+10),'>');
Method = token;
    
end