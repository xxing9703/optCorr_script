function meta=elmaven_optcor(meta,para)
if nargin==1 %default parameter
    para.abd=[0.0110694,0.003663,0.00015]; %13C 15N 2H, natural abundance
    para.imp=[0.01,0.01,0.01]; %tracer impurity
    para.purity=0.99;
    para.ppm=5;
    para.resolution=140000;
    para.option=1; %signal reading method 1. highest peak within ppm 2. nearest peak, 3 sum
    para.parallet=false;
end
sz=length(meta);
for i=1:sz
   items=size(meta(i).original_pct,2);
   fprintf(['working on ',num2str(i),'/',num2str(sz)]);
 formula=meta(i).formula;
 mz=meta(i).mz;
 type=meta(i).type;
 [~,atoms]=str2mass(formula);
 dmz=mz*para.ppm*1e-6;
 fwhm=mz/para.resolution*sqrt(mz/200);
 corr_abs=[];corr_pct=[];
 for j=1:items
    fprintf('.');
    distin=meta(i).original_pct(:,j);
    para.imp(meta(i).type)=1-para.purity; % update impurity value from main input
    [x,~,ysim]=optcorr(mz,atoms,type,para.abd,para.imp,fwhm,dmz,distin,para.option);
    x=x/sum(x);
    ysim=ysim/sum(ysim);
    corr_pct(:,j)=x;  %normalized
    factor=meta(i).original_abs(1,j)/ysim(1);
    corr_abs(:,j)=x*factor;
 end
 fprintf('\n');
    %idx=meta(i).reduced_idx;  
    meta(i).corr_abs=corr_abs;
    meta(i).corr_pct=corr_pct;
    meta(i).corr_tic=sum(meta(i).corr_abs,1);  
 end
end