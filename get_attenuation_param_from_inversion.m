function [gamma,siggamma,Qo,sigQo,alpha,sigalpha]=get_attenuation_param_from_inversion(dir_inversion)

attenuation_param_file=[dir_inversion 'attenuation.out'];
fid=fopen(attenuation_param_file,'r');
info_attenuation=textscan(fid,'%s %f %f');
fclose(fid);

gamma=info_attenuation{2}(1);
siggamma=info_attenuation{3}(1);
Qo=info_attenuation{2}(2);
sigQo=info_attenuation{3}(2);
alpha=info_attenuation{2}(3);
sigalpha=info_attenuation{3}(3);

endfunction