function meta=elmaven_isocor(meta,para)
if nargin==1 %default parameter
    para.abd=[0.0110694,0.003663,0.00015]; %13C 15N 2H, natural abundance
    para.imp=[0.01,0.01,0.01]; %13C 15N 2H, tracer impurity
end
for i=1:length(meta)
 atom1.num=meta(i).atom_num(meta(i).type);
 atom1.abd=para.abd(meta(i).type);
 %atom1.imp=para.imp(meta(i).type);
 atom1.imp=1-para.purity; 
 [corr_abs,~,corr_pct]=isocorr_due(meta(i).original_abs,atom1);
 %idx=meta(i).reduced_idx;  
 meta(i).corr_abs=corr_abs;
 meta(i).corr_pct=corr_pct;
 meta(i).corr_tic=sum(meta(i).corr_abs,1);  

end