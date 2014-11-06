function plot_ts(pname)
%  clear
format long
close all
addpath(genpath('~/octave'));
if nargin() ~= 1
    [fname,pname]=uigetfile('list_*.txt','Select a list of filenames for one event','/home/stephane/DATA/ON/Earthquakes/Ev_from_SC3');
else
    fname='list_sac_files.txt'
end
% locate final dir name within the path
pos=strfind(pname,'/');
if length(pos) == 1
   disp('Definition of output directory will not work')
   return
else
   ev_name_dir=pname(pos(length(pos)-1)+1:length(pname));
end

[dist,dur,sta,fich]=textread(strcat(pname,'distance_duration.txt'),'%f %f %s %s');
sta_unique=unique(sta);
for i=1:length(sta_unique)
   figure(i)
   sta_unique(i)
   index=find(strcmp(sta_unique(i),sta) == 1);
   max_time=0;
   for j=1:length(index)
      fich(index(j))
      S=readsac([pname char(fich(index(j)))]);
      time=[0:S.DELTA:S.DELTA*(length(S.DATA1)-1)]';
      data=detrend((S.DATA1-mean(S.DATA1)),'linear');
      Fnyquist=1/(2*S.DELTA);
      w1=1.0/Fnyquist;
      w2=10.0/Fnyquist;
      [B,A] = butter(4,[w1 w2]);
      data_filter1=filter(B,A,data);
      subplot(3,1,j)
      plot(time,data_filter1)
      hold on
      plot([time(ceil(S.T0/S.DELTA)) time(ceil(S.T0/S.DELTA))],[max(data_filter1) min(data_filter1)],'k')
      plot([time(ceil((S.T0+dur(index(j)))/S.DELTA)) time(ceil((S.T0+dur(index(j)))/S.DELTA))],[max(data_filter1) min(data_filter1)],'k')
      hold off
      if S.NPTS*S.DELTA > max_time, max_time=S.NPTS*S.DELTA,end
   end
   for j=1:length(index)
      subplot(3,1,j)
      xlim([0 max_time])
   end
end

