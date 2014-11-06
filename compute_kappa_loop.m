pname_base='/home/stephane/DATA/ON/Earthquakes/Ev_from_SC3/2012'
list_dirs=glob([pname_base '/*']);
for i=1:length(list_dirs)
    compute_kappa(list_dirs{i})
end
