function f_file = get_sampling_frequency(filename)
    [fid, ~] = fopen(filename, 'rt');
    fmt = '%s';
    frewind(fid)
    C = textscan(fid, fmt, 3);
    fclose(fid);
    
    C = C{1};
    f_file = str2double(C{2});
    if isnan(f_file)
        f_file = str2double(C{3});
    end
end