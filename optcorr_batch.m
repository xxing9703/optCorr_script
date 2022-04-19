%example usage: it takes Maven/elMaven output .csv file, exports the corrected in .xlxs
%optcorr_batch();  use openDialog to load input file, using all default settings  
%optcorr_batch('resolution',100000,'purity',0.98)
%optcorr_batch('solver','isocorr','resolution',480000,'ppm',10,'purity',0.995)

%optional input settings (pairwise)
% 'fname': [] %input path/file name.  default is blank (use openDialog)
% 'solver': 'optcorr' (default) or 'isocorr'
% 'purity': 0.99(default) purity of tracer
% 'abundance': [0.0110694,0.003663,0.00015] (default) %natural abundance for 13N, 15N, 2D
% 'resolution': 140000 (default) % instrument resolution @mz=200 (optcorr only) 
% 'ppm':  5(default)  ppm settings in elMaven signal reading setting.(optcorr only) 
% 'option': 1(default) signal reading option: 1)heighest, 2)nearest,3.sum (optcorr only) 
% 'parallel': false(default)  This is to be implemented in the future
function optcorr_batch(varargin)
  p = inputParser;
  addParameter(p,'fname',[]);
  addParameter(p,'solver','optcorr');
  addParameter(p,'resolution',140000);
  %addParameter(p,'abundance',[0.0110694,0.003663,0.00015]);
  addParameter(p,'abundance',[0.0107,0.00364,0.000115]);%use the same as accucorr
  addParameter(p,'purity',0.99);
  addParameter(p,'ppm',5);
  addParameter(p,'option',1);
  addParameter(p,'parallel',false);
  parse(p,varargin{:});   
  fname=p.Results.fname;
  solver=p.Results.solver;
  para.resolution=p.Results.resolution; 
  para.abd=p.Results.abundance;  
  para.purity=p.Results.purity; %scaler, depending on tracer type
  para.imp=[0.01,0.01,0.01]; % default array
  para.ppm=p.Results.ppm;
  para.option=p.Results.option;
  para.parallel=p.Results.parallel;
  para.solver=solver;
  
  

if isempty(fname)
    [filename,filepath]=uigetfile('*.csv');   
    fname=fullfile(filepath,filename);
end
fname
para.fname=fname;
fprintf('Please verify the following settings are correct!\n');
fprintf(['purity = ',num2str(para.purity),'\n']);
fprintf(['solver = ',solver,'\n']);
warning('off');
tic
[meta,A,start_col]=elmaven_cor_in(fname); 
if strcmp(solver,'optcorr')
    fprintf(['--resolution = ',num2str(para.resolution),'\n']);
    fprintf(['--ppm = ',num2str(para.ppm),'\n']);
    meta=elmaven_optcor(meta,para);
    fout='_optcor';
elseif strcmp(solver,'isocorr')
    meta=elmaven_isocor(meta,para);
    fout='_isocor';
elseif strcmp(solver,'accucor')
    para.option=3;
    meta=elmaven_optcor(meta,para);
    fout='_accucor';
else
    fprintf('select an algorithm: optcorr, isocorr or accucor');
    return
end
para.toc=toc;
elmaven_cor_out(fname,fout,A,meta,start_col,para)
fprintf('Done!\n');
