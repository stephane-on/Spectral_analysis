function durationS2=compute_5_95_nrj(t,a,ts,dt);

%njr_from_ts=sum(a(find(t >= ts)).^2);
njr_from_ts=0;
for i=min(find(t >= ts)):length(t)
    njr_from_ts=njr_from_ts+a(i)^2;
    husid(i)=njr_from_ts;
end

%  figure(5)
%  plot(t,a./max(a),'b')
%  hold on
%  plot(t,husid./max(husid),'r')

sum=0;
testistart=0;
testiend=0;
for i=min(find(t >= ts)):length(t)
    sum=sum+a(i)^2;
    if (sum >= 0.05*njr_from_ts && testistart == 0)
        istart=i;
        testistart=1;
    end
    if (sum >= 0.95*njr_from_ts && testiend == 0)
        iend=i;
        testiend=1;
    end
end

durationS2=(iend-istart)*dt;
