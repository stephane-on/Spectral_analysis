%%%% 
%%%% Reads files of x-y pairs of frequency and acceleration FFT
%%%%   (these are already cut in sac to fE, fX - see Douglas et al., 2010)
%%%% Discards last lines if zero
%%%% Computes kappa based on the downward linear trend of data 
%%%%   when plotted on lin-ln axis (x-lny)
%%%% Uses both standard and robust linear least square regression
%%%%
%%%% olga, irsn, mar 2011
%%%%

function [kappa_ls, int_ls, kappa_rob, int_rob] = kappa1file(xy)  

N=length(xy);

x=zeros(N,1);
y=zeros(N,1);

for i = 1:N;
    x(i)=xy(i,1);
    y(i)=xy(i,2); 
end

% x=x';
% y=y';
%size(x)
yy=log(y);     %%%%%%   !!!!!! this is loge, not log10 !!!!!!

%  ls=regress(yy,[ones(N,1) x]);
ls=polyfit(x,yy,1);
%  rob=robustfit(x,yy);
rob=ls;

%%% get values from vectors (dont want to produce vectors as results)
%%% the line is of the form: yy=lny=ax+b, 
%%% the vector (e.g. rob) gives rob(1)=b and rob(2)=a
%%% we name a=lamda (slope) and b=int (intercept)
%  lamda_ls = ls(2);
%  int_ls = ls(1);
%  lamda_rob = rob(2);
%  int_rob = rob(1);
lamda_ls = ls(1);
int_ls = ls(2);
lamda_rob = rob(1);
int_rob = rob(2);

%%% get kappa from lamda
kappa_ls=-lamda_ls/pi;
kappa_rob=-lamda_rob/pi;

