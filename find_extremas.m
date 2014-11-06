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
