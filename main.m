% usage: 
% optcorr_batch()  use FileOpenDialog to load elmaven output file(.csv). Use all default settings,
% optcorr_batch(name,value) use name-value pair to change the default parameter settings: 
% example:

fname ='example1.csv'; % default is []
solver = 'optcorr';  % 'optcorr' or 'isocorr'  default is optcorr
resolution = 140000;  % instrument resolution, default=140000
ppm = 5;  %ppm setting used in ElMaven, default = 5 
purity = 0.99;  %purity of tracer, default = 0.99
%----------------------------------------------------------------RUN
optcorr_batch('fname',fname,'solver',solver,'resolution',resolution,'ppm',ppm,'purity',purity)

%results will be saved as output file: "example1_optcorr.xlsx"