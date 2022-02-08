% input: out, tb from MID simulation, 
%   res-- MS resolution
function S=simul(out,tb,fwhm,maxM,mzWindow,stepsize)

n=length(out);
%mzWindow=0.04; % in Da, for each integer m/z, centered at 0,  watch [-0.04, 0.04]
%stepsize=0.00001;% delta m/z window on each side of integer
%maxM=5;  % M start from 0, up to maxM
maxM=min(maxM,n-1); % maxM cannot exceed n-1 from out.

xx0=-mzWindow:stepsize:mzWindow;
sigma = fwhm;
S(1).sigma=sigma;
S(1).stickX=0;
S(1).stickY=out(1).pct/100;
S(1).spectX=xx0;
S(1).spectY=out(1).pct/100*exp(-xx0.^2/sigma^2*4*log(2));
S(1).str='P';
% S(1).pct={num2str([S(1).stickY*100]','%.2f%%')};
% S(1).pctrelative={'100%'};
S(1).pct=S(1).stickY;
S(1).pctrelative=1;
S(1).top=S(1).stickY;
S(1).toprelative=1;
S(1).sum=S(1).top;
S(1).sumrelative=1;
S(1).toploc=0;
S(1).pk=S(1).top; %pk signal and location from findpeaks
S(1).loc=0;
mz_defect=[tb([tb.pctGrp]==100).dmz];
for i=2:maxM+1
   M=mz_defect(i);
   %sigma = out(i).mz/res; %mass resolution in delta mz  at current mz
   xx = M-mzWindow:stepsize:M+mzWindow;  %x array
   tb_select=tb(abs([tb.dmz]-M)<0.5);
   S(i).stickX=[tb_select.dmz];
   S(i).stickY=[tb_select.ab];  
      y=zeros(length(tb_select),length(xx)); %stores yy      
     for k=1:length(tb_select)
      y(k,:)=tb_select(k).ab*exp(-(xx-tb_select(k).dmz).^2/sigma^2*4*log(2));
     end
     if k==1
        yy=y;
     else
        yy=sum(y);        
     end
   S(i).spectX=xx; 
   S(i).spectY=yy;
   S(i).sigma=sigma;
   S(i).str={tb_select.str};
   S(i).pct=S(i).stickY;
   S(i).pctrelative=S(i).stickY/S(1).stickY;
   [top,topindex]=max(S(i).spectY);
   S(i).top=top;
   S(i).toprelative=top/S(1).stickY;   
   S(i).toploc=S(i).spectX(topindex);
   S(i).sum=sum(S(i).stickY);
   S(i).sumrelative=S(i).sum/S(1).stickY; 
   try
   [pk,loc]=findpeaks(S(i).spectY,S(i).spectX);
    S(i).pk=pk;
    S(i).loc=loc;
   catch
    S(i).pk=0;
    S(i).loc=mean(S(i).spectX);
   end
end

