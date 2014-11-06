function compute_fourier_spectra
clear
format long
close all
addpath(genpath('~/octave'));

pname=uigetdir('/home/stephane/DATA/ON/Earthquakes','Select an event directory (which includes spectra/ sub-directory):');
pname=strrep(strcat(pname,'/'),'//','/');

fich_t_picks=[pname 'spectra/list_picks.txt'];
fid_picks_t=fopen(fich_t_picks,'r');
info_peaks_t=textscan(fid_picks_t,'%s %s %s %f %f %f %f %f %f');
fclose(fid_picks_t);
%figure(1)
%clf
%plot_all_time_series(pname,info_peaks_t);
%pause(1)

% Data are assumed to be in m/s
answer=inputdlg('Conversion factor from data unit to m/s (if data are in m/s, conv=1, if data are in nm/s, conv=1e-9).','Enter unit converion factor',1,{'1'});
conv=str2double(answer{1});

for i=1:length(info_peaks_t{1,1})
    sacfile=char(info_peaks_t{1,3}(i));
    fprintf(1,'%s\n',[pname sacfile]);
    %% Read SAC data
    S=readsac([pname sacfile]);
    Fnyquist=1/(2*S.DELTA);
    clear time data_tmp data
    time=[0:S.DELTA:S.DELTA*(length(S.DATA1)-1)]';
    %% demean, detrend
    data_tmp=detrend(S.DATA1*conv,'constant');
    data=detrend(data_tmp,'linear');
    clear data_tmp
    data=data.*tukeywin(length(data),0.05);
    %% Filter for FFT
    w1=0.1/Fnyquist;
    w2=0.9999;
    clear B A data_filter
    [B,A] = butter(4,[w1 w2]);
    data_filter=filter(B,A,data);
    tnbeg=info_peaks_t{1,4}(i);
    tnend=info_peaks_t{1,5}(i);
    tpbeg=info_peaks_t{1,6}(i);
    tpend=info_peaks_t{1,7}(i);
    tsbeg=info_peaks_t{1,8}(i);
    tsend=info_peaks_t{1,9}(i);

    %% Compute Fourier spectra
    % Fourier transform of velocity
    NdatNoise=length(data(round(tnbeg/S.DELTA):round(tnend/S.DELTA)));
    NdatS=length(data(round(tsbeg/S.DELTA):round(tsend/S.DELTA)));
    if NdatS >= NdatNoise
        NFFT=2^nextpow2(NdatS);
    else
        NFFT=2^nextpow2(NdatNoise);
    end
    clear data_for_fft
    data_for_fft=data_filter(round(tsbeg/S.DELTA):round(tsend/S.DELTA));
    data_for_fft=data_for_fft.*tukeywin(length(data_for_fft),0.05);
    [f2,spec2,df,norm]=fourier_transform(data_for_fft,S.DELTA,NFFT,1,length(data(round(tsbeg/S.DELTA):round(tsend/S.DELTA))));
    % Take out f=0 Hz
    f=f2(2:length(f2));
    vel_spec_S=abs(spec2(2:NFFT/2+1));
    clear data_for_fft
    data_for_fft=data_filter(round(tnbeg/S.DELTA):round(tnend/S.DELTA));
    data_for_fft=data_for_fft.*tukeywin(length(data_for_fft),0.05);
    [fn2,specn2,df,norm]=fourier_transform(data_for_fft,S.DELTA,NFFT,1,length(data(round(tnbeg/S.DELTA):round(tnend/S.DELTA))));
    vel_spec_N=abs(specn2(2:NFFT/2+1));
    if strcmp(S.IDEP,'IACC') == 1
        fprintf(1,'%s\n','Data are not velocity data');
        vel_spec_S=abs(spec2(2:NFFT/2+1))./(2*pi*f)';
        vel_spec_N=abs(specn2(2:NFFT/2+1))./(2*pi*f)';
    end
    fid_spectra=fopen([pname '/spectra/' sacfile '.spc'],'w');
    fprintf(fid_spectra,'%f %f\n',NdatS*S.DELTA,NdatNoise*S.DELTA);
    fprintf(fid_spectra,'%e %e %e\n',[f' vel_spec_S vel_spec_N]');
    fclose(fid_spectra);
end

endfunction
