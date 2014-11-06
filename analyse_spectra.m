function analyse_spectra
format long
clear
%  close all
addpath(genpath('~/octave'));

figure(2); clf; figure(3); clf

pname=uigetdir('/home/stephane/DATA/ON/Earthquakes/','Select an event directory (which includes spectra/ sub-directory):');
pname=strrep(strcat(pname,'/'),'//','/');
if exist(strcat(pname,'spectra'),'dir') == 0
    disp('Selected directory does not include spectra/ sub-directory. Stopping.')
    return
end

answer=inputdlg({'Enter gamma:','Enter Qo:','Enter alpha:','Enter vs (m/s):'},'Give attenuation parameters',1,{'1.0','2500.0','0.0','3500.0'});
gamma=str2double(answer{1});
Qo=str2double(answer{2});
alpha=str2double(answer{3});
vs=str2double(answer{4});

fich_info=[pname 'spectra/info_att_par.txt'];
if exist(fich_info,'file') == 2
    backup_file(fich_info);
end
fid=fopen(fich_info,'w');
fprintf(fid,'%s %5.2f\n','Gamma_used: ',gamma);
fprintf(fid,'%s %8.2f\n','Qo_used: ',Qo);
fprintf(fid,'%s %5.2f\n','alpha_used: ',alpha);
fprintf(fid,'%s %8.2f\n','vs_used: ',vs);
fclose(fid);
disp('---------------------------------------------------------------------------')
fprintf(1,'%s %5.2f %s %8.2f %s %5.2f\n','Gamma used: ',gamma,' - Qo used: ',Qo,' - alpha used: ',alpha);
disp('---------------------------------------------------------------------------')

fid_picks_f=fopen([pname 'spectra/list_picks_f.txt'],'r');
info_peaks_f=textscan(fid_picks_f,'%s %f %f %f %s');
fclose(fid_picks_f);

fid_out=fopen('data_corrected.in','w');
fid_out2=fopen('data_corrected_smooth_for_gmt.txt','w');
i_rec_good=0;
for i=1:length(info_peaks_f{1,1})
    spcfile=deblank(char(info_peaks_f{1,1}(i)));
    sacfile=strrep(spcfile,'.spc','');
%      A=load([pname 'spectra/' spcfile]);
%      % Convert velocity to acceleration
%      A(:,2)=A(:,2).*(2*pi*A(:,1));
%      A(:,3)=A(:,3).*(2*pi*A(:,1));
    [frequency,vel,velnoise,durS,durN]=read_spectra_file([pname 'spectra/' spcfile]);
    % Convert velocity to acceleration
    acc=vel.*(2*pi*frequency);
%      % Normalise noise by sqrt(durS/durNoise)
%      noise=sqrt(durS/durN)*(velnoise.*(2*pi*frequency));
    
    % Need the distance, read sac file
    S=readsac([pname sacfile]);    
    fmin=info_peaks_f{1,2}(i);
    fmax=info_peaks_f{1,3}(i);
    if fmax <= fmin
        fprintf(1,'%s %s\n','Data excluded because fmax lower than fmin',spcfile);
    else
        % Data between fmin and fmax
%          index=find(A(:,1) >= fmin & A(:,1) < fmax);
        index=find(frequency >= fmin & frequency < fmax);
        clear f_SN_good amp_SN_good
        f_SN_good=frequency(index);
        amp_SN_good=acc(index);
        if isempty(index) ~= 1
            i_rec_good=i_rec_good+1;
            att_geo=1/(S.DIST*1000)^gamma;
            att_ane=exp(-pi*S.DIST*1000*f_SN_good.^(1-alpha)/(Qo*vs));
            acc_spec_source=(1/att_geo)*amp_SN_good./att_ane;
            for j=1:length(acc_spec_source)
                fprintf(fid_out,'%f %e %f %d %s\n',f_SN_good(j),acc_spec_source(j),0.3,i_rec_good,sacfile);
            end
            % Smoothed data with > indicating new record (for use in gmt)
            fprintf(fid_out2,'%s\n','>');
            fprintf(fid_out2,'%f %e\n',[f_SN_good smooth(acc_spec_source,length(acc_spec_source)/10)]');
        else
            fprintf(1,'%s\n','No data between fmin and fmax',spcfile);
        end
    end
    
end
fclose(fid_out);
fclose(fid_out2);

%%
fid_in=fopen('data_corrected.in','r');
A=textscan(fid_in,'%f %f %f %d %s');
fclose(fid_in);
f_all=double(A{1,1});
amp_all=double(A{1,2});
uncertainty_all=double(A{1,3});
irec_all=double(A{1,4});

fid=fopen('Mw_fc_per_record.res','w');
Nrec=max(A{1,4});
irec_ok=0;
Mw_start=3.0;
for i=1:Nrec
    index=find(irec_all == i);
    freq=f_all(index);
    f_name=deblank(char(A{1,5}(index(1))));
    % For some records there may be no data with good signal-to-noise ratio
    if ! isempty(freq)
        log10_acc_obs=log10(amp_all(index));
        sig_log10_acc_obs=uncertainty_all(index);
        [Mw,fc,std_Mw_fc,mod,convergence]=inversion_Mo_fc_lsqr(freq,log10_acc_obs,sig_log10_acc_obs,Mw_start);
        if convergence == 1 && fc < 100
            % REMEMBER THAT fc > 100 Hz is an excluded case, it should not happen for good quality data
            irec_ok=irec_ok+1;
            Mw_tmp(irec_ok)=Mw;
            fc_tmp(irec_ok)=fc;
            fprintf(fid,'%s %3.1f %6.2f\n',f_name,Mw,fc);
        else
           if convergence ~= 1
              fprintf(1,'%s\n',['Inversion did not converge' f_name])
           end
           if fc >= 100
              fprintf(1,'%s\n',['Inverted fc >= 100 Hz' f_name])
           end
        end
    end
end
fclose(fid);

figure(2)
subplot(2,1,1)
plot(Mw_tmp,'*')
ylim([1 5])
xlabel('record number')
ylabel('M_w')
subplot(2,1,2)
semilogy(fc_tmp,'*')
ylim([0.01 50])
xlabel('record number')
ylabel('f_c (Hz)')

% Now invert using all the data from all the records simultaneously
freq=f_all;
log10_acc_obs=log10(amp_all);
sig_log10_acc_obs=uncertainty_all;
[Mw,fc,std_Mw_fc,mod,convergence]=inversion_Mo_fc_lsqr(freq,log10_acc_obs,sig_log10_acc_obs,Mw_start);

if convergence == 1 && fc < 100
    figure(3)
    subplot(2,1,1)
    [NN,X]=hist(Mw_tmp,20);
    bar(X,NN)
    hold on
    xlabel('M_w')
    ylabel('count')
    plot(Mw,0.8*max(NN),'*r')
    hold off
    subplot(2,1,2)
    [NN,X]=hist(log10(fc_tmp),20);
    bar(X,NN)
    hold on
    xlabel('log_{10}(f_c)')
    ylabel('count')
    plot(log10(fc),0.8*max(NN),'*r')
    hold off

    fid_out1=fopen('Mw_fc.out','w');
    fprintf(fid_out1,'%s %3.1f %3.1f\n','Mw=',Mw,std(Mw_tmp));
    %% UNCERTAINTY ON FC IS ESTIMATED FROM std(log10(fc))
    fprintf(fid_out1,'%s %6.2f %6.2f\n','fc=',fc,std(log10(fc_tmp))*fc);
    fclose(fid_out1);
    fid_out2=fopen('model_source_acc.out','w');
    fprintf(fid_out2,'%f %e \n',[freq mod]');
    fclose(fid_out2);

    copyfile('data_corrected.in',[pname 'spectra/data_corrected.in']);
    copyfile('data_corrected_smooth_for_gmt.txt',[pname 'spectra/data_corrected_smooth_for_gmt.txt']);
    copyfile('model_source_acc.out',[pname 'spectra/model_source_acc.out']);
    copyfile('Mw_fc.out',[pname 'spectra/Mw_fc.out']);
    copyfile('Mw_fc_per_record.res',[pname 'spectra/Mw_fc_per_record.res']);
end
