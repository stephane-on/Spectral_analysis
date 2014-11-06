function compute_mR(out_dir_base)
format long
%  clear
close all
addpath(genpath('~/octave'));

[fname,pname]=uigetfile('list_*.txt','Select a list of filenames for one event','/home/stephane/DATA/ON/Earthquakes/');
% locate final dir name within the path
pos=strfind(pname,'/');
if length(pos) == 1
   disp('Definition of output directory will not work')
   return
else
   ev_name_dir=pname(pos(length(pos)-1)+1:length(pname));
end
%  out_dir_base=uigetdir('/home/stephane/WORK/ON/Spectral_analysis','Select directory where you want to put the results:');
%  outdir=strrep(strcat(out_dir_base,'/',ev_name_dir),'//','/');
%  if exist(outdir,'dir') == 0
%      mkdir(outdir)
%  end
outdir=pname

answer=menu("------------------------------------------------------\nWhat is the input data unit (the will convert to micro-m/s):\n------------------------------------------------------","m/s (conversion *1e6)","nm/s (conversion *1e-3)");
if answer == 1
   % Data in m/s, conv will convert in micro-meter/s
   conv=1e6;
else
   % Data in nm/s, conv will convert in micro-meter/s
   conv=1e-3;
end

fid_out=fopen([outdir '/mR_per_record.txt'],'w');
fprintf(fid_out,'%s\n','file_name mR_vel_max mr_vel_peak2peak period_max_vel mr_dis_max mr_dis_peak2peak period_max_dis)')
list_sacfiles=read_file_list_string([pname fname]);

imr=0;
for i=1:size(list_sacfiles,1)
    %% Read SAC data
    if exist('S') == 1
        clear S
    end
    S=readsac([pname list_sacfiles(i,:)]);

    %% Test if components is recognised from S.KCMPNM or from filename
    if exist('comp') == 1
        clear comp
    end
    comp=get_comp_name(deblank(S.KCMPNM),list_sacfiles(i,:));
    if comp == -99
       disp('component not recognized'),disp([pname list_sacfiles(i,:)])
       % If component not recognized, jump to next file
       continue
    end
    
    % Test is necessary data are not empty: S exists means reading sac is ok, and header variables
    % S.DIST and S.DELTA are defined
    if (comp == 'z' && ~isempty(S) && isnan(S.DIST) ~=1 && isnan(S.DELTA) ~=1 )
        clear time Fnyquist w1 w2 B A filtered_data Np Ns
        time=[0:S.DELTA:S.DELTA*(length(S.DATA1)-1)]';
        Fnyquist=1/(2*S.DELTA);
        % w for the filter is a fraction of f/fNyquist
        w1=1.0/Fnyquist;
        %w2=min(0.99999999999,10.0/Fnyquist);
        w2=min(0.9999,10.0/Fnyquist);
        % 4th order ButterWorth filter between 1 and 10 Hz
        [B,A] = butter(4,[w1 w2]);
        % Filter the data and converts to micro-m/s
        %------------Note modify the initial conversion test in order to get micro-meters/s
        filtered_data=filter(B,A,S.DATA1)*conv;

         % Test on the distance mR is computed only if 200<=distance<=1500 km
         if (S.DIST >= 200 && S.DIST <= 1500)
             disp(list_sacfiles(i,:))
             Np=floor((S.A-S.B)/S.DELTA);
             Ns=ceil((S.T0-S.B)/S.DELTA);
%               length(filtered_data)
%               length(S.DATA1)
                          
             [Amax,iAmax]=max(abs(filtered_data(Np:Ns)));
             % Need to have amplitude at the max > 0 for the functions find_next_min and find_previous_min
             % Consequently use signal or -signal depending on the sign of filtered_data(Np+iAmax)
             if filtered_data(Np+iAmax) < 0
                 data_tmp=-filtered_data(Np:Ns);
             else
                 data_tmp=filtered_data(Np:Ns);
             end
             [i1,i2,i3,i4,i5,amp1,amp2]=find_extremas(iAmax,data_tmp);

             % Note: period is not used to compute mR from velocity, it is just for comparison
             % with the period of the max in displacement
             period_max_vel=((i3-i1)*S.DELTA+(i5-i3)*S.DELTA)/2;
             % Warning, Amax always >0, Aprev, Anext always < 0
             amp_max_from_velocity=Amax;
             half_peak2peak_max_from_velocity=(Amax+max(abs(amp1),abs(amp2)))/2;
             
             figure(1)
             % Plot velocity data
             subplot(2,2,1)
             plot(time,filtered_data,'k')
             hold on
             xlabel('time (s)')
             ylabel('velocity (m/s)')
             title(list_sacfiles(i,:),'interpreter','none')
             plot([(S.A-S.B) (S.A-S.B)],[-max(filtered_data) max(filtered_data)],'r:')
             plot([(S.T0-S.B) (S.T0-S.B)],[-max(filtered_data) max(filtered_data)],'r:')
             hold off
             % Zoom to the p-s time window
             subplot(2,2,2)
             plot(time,filtered_data,'k')
             hold on
             xlabel('time (s)')
             ylabel('velocity (m/s)')
             plot(time(Np+iAmax),Amax,'ro')
             plot(time(Np+i1:Np+i5),filtered_data(Np+i1:Np+i5),'b')
             plot(time(Np+i2),amp1,'g+')
             plot(time(Np+i4),amp2,'g+')
             xlim([(Np+iAmax)*S.DELTA-3 (Np+iAmax)*S.DELTA+3])
             ylim([-1.1*Amax 1.1*Amax])
             hold off
             
             % Plot integral of velocity (i.e. displacement)
             % QUESTION DO WE INTEGRATE RAW DATA AND FILTER AFTERWARDS OR DO WE INTEGRATE
             % FILTERED DATA? THE RESULTA ARE VERY DIFFERENT!
%               displacement_data=cumtrapz(time,filtered_data);
%               displacement_data=filter(B,A,cumtrapz(time,filtered_data));
             displacement_data=filter(B,A,cumtrapz(time,S.DATA1*conv));
 
             [Amax,iAmax]=max(abs(displacement_data(Np:Ns)));
             % Need to have amplitude at the max > 0 for the functions find_next_min and find_previous_min
             % Consequently use signal or -signal depending on the sign of filtered_data(Np+iAmax)
             if displacement_data(Np+iAmax) < 0
                 data_tmp=-displacement_data(Np:Ns);
             else
                 data_tmp=displacement_data(Np:Ns);
             end
             [i1,i2,i3,i4,i5,amp1,amp2]=find_extremas(iAmax,data_tmp);

             % Compute period of the max using the average period between second previous
             % and second next maximax
             period_max_dis=((i3-i1)*S.DELTA+(i5-i3)*S.DELTA)/2;
             % Warning, Amax always >0, Aprev, Anext always < 0
             amp_max_from_displacement=Amax;
             half_peak2peak_max_from_displacement=(Amax+max(abs(amp1),abs(amp2)))/2;
             
             % Plot displacement data
             subplot(2,2,3)
             plot(time,displacement_data,'k')
             hold on
             xlabel('time (s)')
             ylabel('displacement (m)')
             plot([(S.A-S.B) (S.A-S.B)],[-max(displacement_data) max(displacement_data)],'r:')
             plot([(S.T0-S.B) (S.T0-S.B)],[-max(displacement_data) max(displacement_data)],'r:')
             hold off
             % Zoom to the p-s time window
             subplot(2,2,4)
             plot(time,displacement_data,'k')
             hold on
             xlabel('time (s)')
             ylabel('displacement (m)')
             plot(time(Np+iAmax),Amax,'ro')
             plot(time(Np+i1:Np+i5),displacement_data(Np+i1:Np+i5),'b')
             plot(time(Np+i2),amp1,'g+')
             plot(time(Np+i4),amp2,'g+')
             xlim([(Np+iAmax)*S.DELTA-3 (Np+iAmax)*S.DELTA+3])
             ylim([-1.1*Amax 1.1*Amax])
             hold off
             
             % Save the graphic
             print('test.ps','-dpsc')
             copyfile('test.ps',[deblank(outdir) deblank(list_sacfiles(i,:)) '.ps']);
             
             imr=imr+1;
             mr_vel_max(imr)=log10(amp_max_from_velocity)+2.3*log10(S.DIST)-2.28;
             mr_vel_peak2peak(imr)=log10(half_peak2peak_max_from_velocity)+2.3*log10(S.DIST)-2.28;
             mr_dis_max(imr)=log10(amp_max_from_displacement/period_max_dis)+2.3*log10(S.DIST)-1.48;
             mr_dis_peak2peak(imr)=log10(half_peak2peak_max_from_displacement/period_max_dis)+2.3*log10(S.DIST)-1.48;
             % Write output
             fprintf(fid_out,'%s %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f\n',list_sacfiles(i,:),mr_vel_max(imr),mr_vel_peak2peak(imr),period_max_vel,mr_dis_max(imr),mr_dis_peak2peak(imr),period_max_dis)
%               fprintf(1,'%s %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f\n',list_sacfiles(i,:),mr_vel_max(imr),mr_vel_peak2peak(imr),period_max_vel,mr_dis_max(imr),mr_dis_peak2peak(imr),period_max_dis)
%               pause(1)
         
         end

    end
end

%%%%%%%%%%mr_moy=0.0;
%%%%%%%%%%if imr >= 1
%%%%%%%%%%    mr_moy=sum(mr)/imr;
%%%%%%%%%%end
%%%%%%%%%%sig_mr_moy=0.0
%%%%%%%%%%if imr >= 2
%%%%%%%%%%    sig_mr_moy=sqrt(sum((mr-mr_moy).^2)/(imr-1));
%%%%%%%%%%end
%%%%%%%%%%fprintf(fid_out,'%s %4.2f %s %4.2f\n','mR=',mr_moy,'+-',sig_mr_moy)
fclose(fid_out);

mR=mean(mr_vel_peak2peak);
sig_mR=std(mr_vel_peak2peak);
Ndata=length(mr_vel_peak2peak);
fid_out=fopen([outdir 'mR.txt'],'w');
fprintf(fid_out,'%s\n','mean_mR std_mR Nsta_mR');
fprintf(fid_out,'%4.2f %4.2f %d\n',mR,sig_mR,Ndata);
fclose(fid_out);

save([outdir 'mR.mat'],'mr_vel_max','mr_vel_peak2peak','mr_dis_max','mr_dis_peak2peak')

figure(2)
plot(mr_vel_max,'*r')
hold on
xlabel('record index')
ylabel('mR')
plot(mr_vel_peak2peak,'or')
plot(mr_dis_max,'*b')
plot(mr_dis_peak2peak,'ob')
legend('vel max','vel peak-to-peak','dis max','dis peak-to-peak')
hold off

figure(3)
[N,X]=hist([mr_vel_max mr_vel_peak2peak mr_dis_max mr_dis_peak2peak]);
bar(X,N)
hold on
xlabel('mR')
ylabel('count')
title(list_sacfiles(i,:),'interpreter','none')
hold off

figure(4)
subplot(2,2,1)
hist(mr_vel_max)
hold on
xlabel('mR')
ylabel('count')
title(list_sacfiles(i,:),'interpreter','none')
subplot(2,2,2)
hist(mr_vel_peak2peak)
subplot(2,2,3)
hist(mr_dis_max)
subplot(2,2,4)
hist(mr_dis_peak2peak)
hold off
             

endfunction

function [i1,i2,i3,i4,i5,Aused1,Aused2]=find_extremas(iAmax,data_tmp)
             % Hypothesis if the edge of the time window is found in one direction
             % then it is long in the other

             if iAmax == 1
                 exit_flag_prev=0;
             else
                 [iprev,Aprev,exit_flag_prev]=find_previous_min(iAmax,data_tmp);
             end
             if iAmax == length(data_tmp) % Case the max is already at the upper edge of the window
                 exit_flag_next=0;
             else
                 [inext,Anext,exit_flag_next]=find_next_min(iAmax,data_tmp);
             end
             if (exit_flag_prev == 0) % Search all next extremas after iAmax
                 i1=iAmax;
                 [i2,Anext1,exit_flag_next]=find_next_min(i1,-data_tmp);
                 Aused1=Anext1;
                 Aused2=Aused1;
                 [i3,Anext3,exit_flag_next]=find_next_min(i2,data_tmp);
                 [i4,Anext4,exit_flag_next]=find_next_min(i3,-data_tmp);
                 [i5,Anext5,exit_flag_next]=find_next_min(i4,data_tmp);
             elseif (exit_flag_next == 0)
                 i5=iAmax;
                 [i4,Aprev4,exit_flag_prev]=find_previous_min(i5,-data_tmp);
                 Aused1=Aprev4;
                 Aused2=Aused1;
                 [i3,Aprev3,exit_flag_prev]=find_previous_min(i4,data_tmp);
                 [i2,Aprev4,exit_flag_prev]=find_previous_min(i3,-data_tmp);
                 [i1,Aprev5,exit_flag_prev]=find_previous_min(i2,data_tmp);
             elseif (exit_flag_prev ~= 0 && exit_flag_next ~= 0)
                 [iprev2,Aprev2,exit_flag_prev2]=find_previous_min(iprev,-data_tmp);
                 [inext2,Anext2,exit_flag_next2]=find_next_min(inext,-data_tmp);
                 if exit_flag_prev2 == 0
                     i1=iprev;
                     Aused1=Aprev;
                     i2=iAmax;
                     i3=inext2;
                     Aused2=Anext2;
                     [i4,Anext4,exit_flag_next]=find_next_min(i3,data_tmp);
                     [i5,Anext5,exit_flag_next]=find_next_min(i4,-data_tmp);
                 elseif exit_flag_next2 == 0
                     i5=inext;
                     Aused2=Anext;
                     i4=iAmax;
                     i3=iprev2;
                     Aused1=Aprev2;
                     [i2,Aprev2,exit_flag_prev]=find_previous_min(i3,data_tmp);
                     [i1,Aprev1,exit_flag_prev]=find_previous_min(i2,-data_tmp);
                 elseif (exit_flag_prev2 ~= 0 && exit_flag_next2 ~= 0)
                     i1=iprev2;
                     i2=iprev;
                     Aused1=Aprev;
                     i3=iAmax;
                     i4=inext;
                     Aused2=Anext;
                     i5=inext2;
                 else
                     disp('Case not taken into account')
                     pause
                 end
             else
                 disp('Case not taken into account')
                 pause
             end
endfunction

function [iprev,Aprev,exit_flag]=find_previous_min(iref,X)
% Find the previous minimum centered on X(iref)
Xtest=X(iref);
iprev=0;
for i=iref-1:-1:1
    if X(i) < Xtest
       Xtest=X(i);
    else
       iprev=i;
       Aprev=X(i);
       exit_flag=1;
       break
    end
    if iprev == 0
        Aprev=0.0;
        exit_flag=0;
    end
end
endfunction

function [inext,Anext,exit_flag]=find_next_min(iref,X)
% Find next minimum centered on X(iref)
Xtest=X(iref);
inext=0;
for i=iref+1:length(X)
   if X(i) < Xtest
       Xtest=X(i);
    else
       inext=i;
       Anext=X(i);
       exit_flag=1;
       break
    end
    % Case no minimum is found in the window (i.e. if the max is very close to the
    % edge of the window)
    if inext == 0
        Anext=0.0;
        exit_flag=0;
    end
end
endfunction
