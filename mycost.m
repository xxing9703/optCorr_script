%y_exp=[63.95  4.80 18.77 10.03  0.35];
%ext = [70,0,10,10];  % input x vector trial
function [err,y_simulated]=mycost(x,y_exp,mz,atoms,type,abd,imp,fwhm,dmz,option)
H=1.00728;
%abd=para.abundance;
%imp=para.impurity;
%ab=[0.0110694,0.003663,0.00015];
%impurity=0.01;
%formula = 'C6H13O9P';
%charge = 1;
%type = 1;
% [mz,atoms]=str2mass(formula);
% mz=mz+H*charge;
%ppm=5;

dim=atoms(type);
[M,Mp]=get_isomatrix(dim,abd(type),imp(type));
ext=propagate(x,dim,M,Mp);
[out,tb]=mid_simulation(mz,atoms,dim,abd,ext,type);
maxM=min(length(out)-1,dim);
if type==4
    maxM=maxM*2;  % updated 7/22/2024
end
mzwindow=fwhm*15;
stepsize=mzwindow*0.0005;
S = simul(out,tb,fwhm,maxM,mzwindow,stepsize);
y_simulated=mid_reading(type,maxM,dmz,S,fwhm,option)*100; %output simulated reading
%padding zeros to make y_exp and y_simulated the same length

n=max(length(y_exp),length(y_simulated));

y_exp(end+1:n)=deal(0);
y_simulated(end+1:n)=deal(0);

[~,topM]=max(y_exp);
y_exp=y_exp/y_exp(topM);
y_simulated=y_simulated/(y_simulated(topM)+1e-9);

err= log10(sum((y_exp-y_simulated).^2));



