function [Mw,sigMw,fc,sigfc,list_dir_event]=get_source_param_from_inversion(dir_inversion)

%dir_inversion='/home/stephane/WORK/ON/Spectral_analysis/inversion_juin2014/';

file_source_param_inversion=[dir_inversion 'source_params.out'];
event_names_codes=[dir_inversion 'event_names_codes.txt'];

if (~exist(file_source_param_inversion,'file') || ~exist(event_names_codes,'file'))
    disp(['inversion source param file or file with event-code association does not exist'])
end

fid=fopen(file_source_param_inversion,'r');
info_source=textscan(fid,'%d %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
fclose(fid);

fid=fopen(event_names_codes,'r');
codes_names=textscan(fid,'%d %s');
fclose(fid);

%check if info_source and codes_names have the same size
if size(info_source{1},1) ~= size(codes_names{1},1)
    disp('Inversion source_param file and events codes_names association do not have the same size. Paused...')
    pause
end

list_dir_event={};
Mw=[]; sigMw=[];
fc=[]; sigfc=[];
for i=1:size(codes_names{2},1)
    
    dir1=fileparts(char(codes_names{2}(i)));
    dir2=fileparts(dir1);
    [dir_base, dir_event]=fileparts(dir2);
    
    check='n';
    for j=1:size(info_source{1},1)
        tmp=char(info_source{2}(j));
        ev_code2=str2num(tmp(7:10));
        if codes_names{1}(i) == ev_code2
            %disp(['Event found' num2str(i) num2str(j)]);
            check='y';
            list_dir_event=strvcat(list_dir_event,dir_event);
            Mw=[Mw;info_source{5}(j)];
            sigMw=[sigMw;info_source{6}(j)];
            fc=[fc;info_source{7}(j)];
            sigfc=[sigfc;info_source{8}(j)];
        end
    end
    
    if strcmp(check,'n') == 1
       disp(['event ' dir_event ' number ' num2str(codes_names{1}(i)) ' not found in source_param file'])
    end
    
end

endfunction