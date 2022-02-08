function elmaven_cor_out(fname,fout,A,meta,start_col,para)
cat_abs=[];cat_pct=[];cat_tic=[];Compound=[];Label=[];
sample_name=A.Properties.VariableNames(start_col:end)';
for i=1:length(meta)
    cat_abs=[cat_abs;meta(i).corr_abs];  %concatenate for csv output
    cat_pct=[cat_pct;meta(i).corr_pct]; 
    cat_tic=[cat_tic;meta(i).tic];
    sz=size(meta(i).original_abs,1);
    Compound=[Compound;repmat({meta(i).name},sz,1)];

    Label=[Label;[0:sz-1]'];
end
cat_pct=floor(cat_pct*1e5)/1e5;
 T_abs=array2table(cat_abs);
 T_abs.Properties.VariableNames=sample_name;
 T_abs=[cell2table(Compound),array2table(Label),T_abs];

 T_pct=array2table(cat_pct);
 T_pct.Properties.VariableNames=sample_name;
 T_pct=[cell2table(Compound),array2table(Label),T_pct];

 T_tic=array2table(cat_tic);
 T_tic.Properties.VariableNames=sample_name;
 T_tic=[cell2table({meta.name}'),T_tic];
 T_tic.Properties.VariableNames{1}='Compound';

 C={'solver',para.solver;...     
    'purity',para.purity;...
    'resolution',para.resolution;...
    'ppm', para.ppm;...
    'runtime',para.toc};
 T_para=cell2table(C);

 [filepath,filename,~] = fileparts(fname);
 fname_S=fullfile(filepath,[filename,fout,'.xlsx']); %fname_S: filename for save
 writetable(A,fname_S,'sheet','Original');
 writetable(T_abs,fname_S,'sheet','Corrected');
 writetable(T_pct,fname_S,'sheet','Normalized');
 writetable(T_tic,fname_S,'sheet','PoolAfterDF');
 writetable(T_para,fname_S,'sheet','logs','WriteVariableNames',0);