function plot_all_time_series(pname,TMP)

figure(1)
clf
Tmax=0;
for i=1:length(TMP{1,1})
    subplot(length(TMP{1,1}),1,i)
    %% Read the time series
    sacfile=char(TMP{1,3}(i));
    tnbeg=TMP{1,4}(i);
    tnend=TMP{1,5}(i);
    tpbeg=TMP{1,6}(i);
    tpend=TMP{1,7}(i);
    tsbeg=TMP{1,8}(i);
    tsend=TMP{1,9}(i);
    pname
    sacfile
    S=readsac(strcat(deblank(pname),sacfile));
    time=[0:S.DELTA:S.DELTA*(length(S.DATA1)-1)]';
    if max(time) > Tmax, Tmax=max(time);, end
    Fnyquist=1/(2*S.DELTA);
    w1=1.0/Fnyquist;
    w2=min(0.9999,10/Fnyquist);
    [B,A] = butter(4,[w1 w2]);
    figure(1)
    titre=[sacfile ' - ' deblank(S.KSTNM) ' - ' deblank(S.KCMPNM) ' - ' num2str(S.DIST) ' km'];
    plot(time,filter(B,A,S.DATA1),["-;" titre ";"],'color','r')
    hold on
    plot([time(ceil(tpbeg/S.DELTA)) time(ceil(tpbeg/S.DELTA))],[max(filter(B,A,S.DATA1)) min(filter(B,A,S.DATA1))],'k')
    plot([time(ceil(tpend/S.DELTA)) time(ceil(tpend/S.DELTA))],[min(filter(B,A,S.DATA1)) max(filter(B,A,S.DATA1))],'k')
    plot([time(ceil(tsbeg/S.DELTA)) time(ceil(tsbeg/S.DELTA))],[max(filter(B,A,S.DATA1)) min(filter(B,A,S.DATA1))],'b')
    plot([time(ceil(tsend/S.DELTA)) time(ceil(tsend/S.DELTA))],[min(filter(B,A,S.DATA1)) max(filter(B,A,S.DATA1))],'b')
    plot([time(ceil(tnbeg/S.DELTA)) time(ceil(tnbeg/S.DELTA))],[max(filter(B,A,S.DATA1)) min(filter(B,A,S.DATA1))],'g')
    plot([time(floor(tnend/S.DELTA)) time(floor(tnend/S.DELTA))],[min(filter(B,A,S.DATA1)) max(filter(B,A,S.DATA1))],'g')
end
for i=1:length(TMP{1,1})
    subplot(length(TMP{1,1}),1,i)
    xlim([0 Tmax])
    hold off
end

end