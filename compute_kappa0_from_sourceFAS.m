function compute_kappa0_from_sourceFAS(sta_name)

% Define the attenuation parameters
%  answer=inputdlg({'Enter gamma','Enter Qo','Enter alpha'},'Attenuation parameters',1,{'1.0','2000','0.3'});
%  gamma=answer{1};
%  Qo=answer{2};
%  alpha=answer{3};
vs=3.5; % S-wave velocity in km/s
dir_inversion='/home/stephane/WORK/ON/Spectral_analysis/inversion_dec2014/';
[Mw,sigMw,fc,sigfc,list_dir_event]=get_source_param_from_inversion(dir_inversion);
[gamma,siggamma,Qo,sigQo,alpha,sigalpha]=get_attenuation_param_from_inversion(dir_inversion);

% Search all spectra files which include sta_name in the filename
list_files={};
list_files{1}=glob(['/home/stephane/DATA/ON/Earthquakes/Ev_from_SC3/201[1-2]/*/spectra/*' sta_name '*.spc']);
list_files{2}=glob(['/home/stephane/DATA/ON/Earthquakes/Ev_from_SC3/201[1-2]/*/spectra/*' tolower(sta_name) '*.spc']);
list_files{3}=glob(['/home/stephane/DATA/ON/IPT/sac/*/spectra/*' sta_name '*.spc']);
list_files{4}=glob(['/home/stephane/DATA/ON/IPT/sac/*/spectra/*' tolower(sta_name) '*.spc']);
list_files{5}=glob(['/home/stephane/DATA/ON/Earthquakes/Ev_from_Marcelo/*/spectra/*' sta_name '*.spc']);
list_files{6}=glob(['/home/stephane/DATA/ON/Earthquakes/Ev_from_Marcelo/*/spectra/*' tolower(sta_name) '*.spc']);
% Create a list of unique directories from the list of spectra files
k=0;
list_dirs={};
for i=1:size(list_files,2)
    if length(list_files{i}) ~= 0
        for j=1:length(list_files{i})
            [dir, name, ext, ver] = fileparts (list_files{i}{j});
            k=k+1;
            list_dirs{k}=strrep(list_files{i}{j},[name ext],"");
        end
    end
end
list_dirs_final=unique(list_dirs);

% Loop over directories
kappa0=[];
for j=1:length(list_dirs_final)
    % Get information about fmin, fmax, fE
    if exist([list_dirs_final{j} 'list_picks_f.txt'],'file')
    
    fid_picks_f=fopen([list_dirs_final{j} 'list_picks_f.txt'],'r');
    info_peaks_f=textscan(fid_picks_f,'%s %f %f %f %s');
    fclose(fid_picks_f);

    % Check if Mw, fc exist for that event
    dir1=fileparts(list_dirs_final{j});
    dir2=fileparts(dir1);
    [dir3,dir_ev]=fileparts(dir2);
    check_source_params='n';
    for i=1:size(list_dir_event,1)
        if strcmp(list_dir_event(i,:),dir_ev)
            check_source_params='y';
            Mw_ev=Mw(i);
            fc_ev=fc(i);
            sigMw_ev=sigMw(i);
            sigfc_ev=sigfc(i);
        end
    end
    if strcmp(check_source_params,'y') == 1
        disp(['Event: ' dir_ev ' - Mw=' num2str(Mw_ev) ' - fc=' num2str(fc_ev)]);
    else
        disp(['Event: ' dir_ev ' - No Mw, fc found in the inversion results.']);
    end
    
    % list all spectra files with filenames including sta_name
    list_files_spc=glob([list_dirs_final{j} '*' sta_name '*.spc']);
    % Search for files with lower case station name
    if length(list_files_spc) == 0
        list_files_spc=glob([list_dirs_final{j} '*' tolower(sta_name) '*.spc']);
    end
    for k=1:length(list_files_spc)
        [dir_spc, sac_name, ext] = fileparts (list_files_spc{k});
        dir_sac=strrep(dir_spc,'spectra','');
        S=readsac([dir_sac sac_name]);
        comp=get_comp_name(deblank(S.KCMPNM),sac_name);
        
        % Use only horizontal components
        if strcmp(comp,'z') == 0
            fprintf(1,'%s %s %s %s\n','filename: ',list_files_spc{k},' component code: ',comp);
            if comp == -99; disp('Component not recognized. Paused...'); pause; end
            fid_spc=fopen(list_files_spc{k},'r');
            dur_tmp=fgetl(fid_spc);
            dur=textscan(dur_tmp,'%f %f');
            spc=textscan(fid_spc,'%f %f %f');
            fclose(fid_spc);
            
            % Get fmin, fmax, fE
            for i=1:size(info_peaks_f{1,1},1)
                if strcmp(info_peaks_f{1,1}(i),[sac_name ext]) == 1
                    fmin=info_peaks_f{1,2}(i);
                    fmax=info_peaks_f{1,3}(i);
                    fE=info_peaks_f{1,4}(i);
                end
            end
            
            % Exclude cases with fmax <= fmin bad data
            if (~isempty(fmin) && ~isempty(fmax) && ~isempty(fE) && fmax > fmin)
                freq=spc{1,1};
                spc_acc_signal=2*pi*freq.*spc{1,2};
                spc_acc_noise=2*pi*freq.*spc{1,3};
                index=find(freq >= fmin & freq <= fmax);
                
                % Correct acceleration spectra from attenuation
                Rh=sqrt(S.DIST^2+S.EVDP^2);
                cor_att=exp(-(pi*Rh*freq.^(1-alpha))./(Qo*vs))./Rh^gamma;
                spc_acc_signal_corr_att=spc_acc_signal./cor_att;
                
                figure(1)
                loglog(freq,spc_acc_noise,'c')
                hold on
                loglog(freq,spc_acc_signal)
                loglog(freq(index),spc_acc_signal(index),'r')
                if strcmp(check_source_params,'y') == 1
                    % Compute the full FAS acclelration model to plot
                    site=ones(length(freq),1);
                    rad=0.53;
                    rho=2800;
                    [acc_mod]=FAS_acc_Mo_fc_gamma_Qo_alpha_site(freq,Mw_ev,fc_ev,Rh,gamma,Qo,alpha,site,rad,vs,rho);
                    loglog(freq,acc_mod,'k','linewidth',2)
                end
                hold off
                
                figure(2)
                semilogy(freq,spc_acc_signal_corr_att)
                hold on
                semilogy(freq(index),spc_acc_signal_corr_att(index),'r')
                if strcmp(check_source_params,'y') == 1
                
                    semilogy([fc_ev fc_ev],[min(spc_acc_signal_corr_att) max(spc_acc_signal_corr_att)],'g','linewidth',2)
                    semilogy([fE fE],[min(spc_acc_signal_corr_att) max(spc_acc_signal_corr_att)],'m','linewidth',2)
                    semilogy(freq,acc_mod./cor_att,'k','linewidth',2)
                    % If fE > fmin and fmax > 10
                    if (fE > fc_ev && fmax > 10)
                        % Frequencies between fE and fmax
                        index2=find(freq >= fE & freq <= fmax);
                        [P,S]=polyfit(freq(index2),log10(spc_acc_signal_corr_att(index2)),1);
                        kappa0=[kappa0; -P(1)];
                        semilogy(freq(index2),power(10,P(1)*freq(index2)+P(2)),'c','linewidth',2)
                    end
                end
                hold off
                
            end
            
%              pause
            
        end
    end
    end
end

if ~isempty(kappa0)
    hist(kappa0)
    xlim([-0.03 0.07])
    xlabel('kappa0 (s)','fontsize',14)
    ylabel('counts','fontsize',14)
    title(sta_name,'fontsize',14)
    if ~exist('../plot_kappa','dir'), system([mkdir ' ../plot_kappa']);, end
    print(['../plot_kappa/hist_kappa0_from_source_FAS_' sta_name '.png'],'-dpng')
end

endfunction

function [acc_mod]=FAS_acc_Mo_fc_gamma_Qo_alpha_site(freq,Mw,fc,Rh,gam,Qo,alpha,site,rad,vs,rho)
% -------------------------------------------------------------------------
% Input:
% freq= vector of frequencies
% Rh= hypocentral distance in km
% vs= S-wave velocity (km/s)
%--------------------------------------------------------------------------

%acc_mod=zeros(length(freq),1);
fac=pi*Rh*1000/(vs*1000);
cste=2*rad/(4*pi*rho*(vs*1000)^3);
Mo=10^(1.5*Mw+9.1);
source=(cste*Mo*(2*pi*freq).^2)./(1+(freq./fc).^2);
att_geo=1/(Rh*1000)^gam;
att_ane=exp(-(fac*freq.^(1-alpha))./Qo);
acc_mod=((source.*att_geo).*att_ane).*site;

endfunction