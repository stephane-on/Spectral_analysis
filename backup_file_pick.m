function backup_file_pick(f_name)
if nargin() ~= 1
    disp('Function backup_file_pick must be called with an argument (file name)')
    return
end
if exist(f_name,'file') == 2
    f_name_save=[f_name '_' date '.sav'];
    i=0;
    while exist(f_name_save,'file') == 2
       i=i+1;
       f_name_save=[f_name '_' date '.sav.' num2str(i)];
    end
    system(['cp ' f_name ' ' f_name_save]);
end
endfunction