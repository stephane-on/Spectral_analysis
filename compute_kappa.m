function compute_kappa(pname_tmp)
format long
%  clear
addpath(genpath('~/octave'));

if nargin == 0
%     disp('no input arg')
   pname=uigetdir('/home/stephane/DATA/ON/Earthquakes/','Select an event directory (which includes spectra/ sub-directory):');
elseif nargin == 1
   pname=pname_tmp;
   disp(pname)
else
   disp('Number of input arguments in compute_kappa not correct')
   return
end

%  pname=uigetdir('/home/stephane/DATA/ON/Earthquakes/','Select an event directory (which includes spectra/ sub-directory):');
pname=strrep(strcat(pname,'/'),'//','/');

if exist(strcat(pname,'spectra'),'dir') == 0
    disp('Selected directory does not include spectra/ sub-directory. Stopping.')
    return
end

%------------------- Plot time series
%  fid_picks_t=fopen([pname 'spectra/list_picks.txt'],'r');
%  info_peaks_t=textscan(fid_picks_t,'%s %s %s %f %f %f %f %f %f');
%  fclose(fid_picks_t);
%  figure(1)
%  clf
%  plot_all_time_series(pname,info_peaks_t)
%-------------------

if exist([pname 'spectra/list_picks_f.txt'],'file') == 2
   fid_picks_f=fopen([pname 'spectra/list_picks_f.txt'],'r');
   info_peaks_f=textscan(fid_picks_f,'%s %f %f %f %s');
   fclose(fid_picks_f);
else
   disp([pname 'spectra/list_picks_f.txt does not exist. Quitting.']);
   return
end

% Retrieve source parameters is exists
if exist([pname 'spectra/Mw_fc.out'],'file') == 2
   fid_info_source=fopen([pname 'spectra/Mw_fc.out'],'r');
   info_source=textscan(fid_info_source,'%s %f %f');
   fclose(fid_info_source);
   Mw=info_source{2}(1);
   fc=info_source{2}(2);
   Mo=10^(1.5*Mw+9.1);
else 
   if exist([pname 'mR.txt'],'file') == 2
      fid_info_source=fopen([pname 'mR.txt'],'r');
      info_source=textscan(fid_info_source,'%f %f %f','headerlines',1);
      fclose(fid_info_source);
      Mw=info_source{1}(1);
   else
      answer=inputdlg('No Mw or mR found','Enter Mw',1,{'3'});
      Mw=str2double(answer{1});
   end
   % Compute fc for a stress drop of 1MPa
   Mo=10^(1.5*Mw+9.1);
   fc=0.37*3500*(16*1e6/(7*Mo))^(1/3);
end
% p=[Mw, fc, gamma, Q, R]
rad=0.55;
vs=3500;
rho=2800;
Q=2500;
gamma=1;
model_FAS=@(p,f) ((2*rad/(4*pi*rho*vs^3))*(10^(1.5*p(1)+9.1))*(((2*pi*f).^2)./(1+(f./p(2)).^2)).*exp(-pi*p(5)*f./(p(4)*vs)))/p(5)^p(3);

fich_kappa=[pname 'spectra/kappa.out'];
if exist(fich_kappa,'file') == 2
    backup_file(fich_kappa);
end
fid_kappa=fopen(fich_kappa,'w');

for i=1:length(info_peaks_f{1,1})
    spcfile=deblank(char(info_peaks_f{1,1}(i)));
    sacfile=strrep(spcfile,'.spc','');
    fmin=info_peaks_f{1,2}(i);
    fmax=info_peaks_f{1,3}(i);
    fe=info_peaks_f{1,4}(i);
    % Need the distance, read sac file and sampling rate
    S=readsac([deblank(pname) sacfile]);
    %% Read the spectra - Convert velocity to acceleration
%      A=load([pname 'spectra/' sacfile '.spc']);
%      freq=A(:,1);
%      acc=A(:,2).*(2*pi*freq);
%      noise=A(:,3).*(2*pi*freq);
    [freq,vel,velnoise,durS,durN]=read_spectra_file([pname 'spectra/' sacfile '.spc']);
    % Convert velocity to acceleration
    acc=vel.*(2*pi*freq);
    % Normalise noise by sqrt(durS/durNoise)
    noise=sqrt(durS/durN)*(velnoise.*(2*pi*freq));
    
    figure(2)
    clf
    %% plot raw signals
    % log-log plot
    subplot(1,2,1)
    loglog(freq,noise,'Color',[211/255 211/255 211/255])
    hold on
    loglog(freq,acc,'Color',[0.68 0.85 0.90])
    %-------------------------------------------------------
    % multi-taper spectra
    if (exist([pname 'spectra/' sacfile '.spc_n'],'file') == 2 && exist([pname 'spectra/' sacfile '.spc_s'],'file') == 2)
       [freq2,vel2,velnoise2,durS2,durN2]=read_spectra_file_mt([pname 'spectra/' sacfile '.spc_s']);
       % Convert velocity to acceleration
       acc2=vel2.*(2*pi*freq2);
       % Normalise noise by sqrt(durS/durNoise)
       noise2=sqrt(durS2/durN2)*(velnoise2.*(2*pi*freq2));
       loglog(freq2,noise2,'k')
       loglog(freq2,acc2,'b')
    end
    % Plot model_FAS
    loglog(freq,model_FAS([Mw fc gamma Q S.DIST*1000],freq),'c')
    %-------------------------------------------------------
    xlim([min(freq) max(freq)])
    ylim([min(acc) max(acc)])
    title([sacfile ' - ' num2str(S.DIST) ' km'],'Interpreter', 'none')
    xlabel('frequency (Hz)')
    ylabel('FAS')
    acc_smooth=smooth(acc,20);
    
    % lin-log plot
    subplot(1,2,2)
    semilogy(freq,noise,'Color',[211/255 211/255 211/255])
    hold on
    semilogy(freq,acc,'Color',[0.68 0.85 0.90])
    %-------------------------------------------------------
    % multi-taper spectra
    if (exist([pname 'spectra/' sacfile '.spc_n'],'file') == 2 && exist([pname 'spectra/' sacfile '.spc_s'],'file') == 2)
       semilogy(freq2,noise2,'k')
       semilogy(freq2,acc2,'b')
    end
    % Plot model_FAS
    semilogy(freq,model_FAS([Mw fc gamma Q S.DIST*1000],freq),'c')
    semilogy([fc fc],[min(acc) max(acc)],'k','LineWidth',2)
    %-------------------------------------------------------
    xlim([min(freq) max(freq)])
    ylim([min(acc) max(acc)])
    xlabel('frequency (Hz)')
    ylabel('FAS')

    % Data between fe and fmax
    index_kappa=find(freq >= fe & freq < fmax);
    if isempty(index_kappa) == 0
        subplot(1,2,1)
        loglog(freq(index_kappa),acc(index_kappa),'r')
        subplot(1,2,2)
        semilogy(freq(index_kappa),acc(index_kappa),'r')
    end
    pause(1)

    clear ButtonName
    if isempty(index_kappa) == 0
%          figure(2)
%          subplot(1,2,2)
%          hold on
%          semilogy(freq(index_kappa),acc_smooth(index_kappa),'g')
        
        ButtonName = questdlg('Use this data to compute kappa?','Yes','No');
        if strcmp(ButtonName,'Yes') == 1
            % Kappa computed from smoothed data
            [kappa_ls, int_ls, kappa_rob, int_rob] = kappa1file([freq(index_kappa) acc_smooth(index_kappa)]);
            % Kappa computed from raw data
%              [kappa_ls, int_ls, kappa_rob, int_rob] = kappa1file([freq(index_kappa) acc(index_kappa)])
            fprintf(fid_kappa,'%f %f %f %f %f %s\n',S.DIST,kappa_ls,int_ls,kappa_rob,int_rob,[sacfile '.spc']);
%              fprintf(1,'%f %f %f %f %f %s\n',S.DIST,kappa_ls,int_ls,kappa_rob,int_rob,[sacfile '.spc']);
%              pause(1)
            
            subplot(1,2,1)
            loglog([min(freq(index_kappa)) max(freq(index_kappa))],[exp(int_rob-kappa_rob*pi*min(freq(index_kappa))) exp(int_rob-kappa_rob*pi*max(freq(index_kappa)))],'r')
            subplot(1,2,2)
            semilogy([min(freq(index_kappa)) max(freq(index_kappa))],[exp(int_rob-kappa_rob*pi*min(freq(index_kappa))) exp(int_rob-kappa_rob*pi*max(freq(index_kappa)))],'r')

            figure(3)
            plot(S.DIST,kappa_ls,'r+')
            plot(S.DIST,kappa_rob,'bo')
            xlim([0 1500])
            ylim([-0.05 0.3])
            hold on
        end
        figure(2)
        subplot(1,2,1)
        hold off
        subplot(1,2,2)
        hold off
    else
        % If no data between fe and fmax, no kappa computed
        ButtonName='No';
    end
end
fclose(fid_kappa);
