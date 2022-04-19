% charge = -1,0,1
% type = 1, 2, 3,  (C, N, N)
% dmz = m/z window, tolerance for mid_reading
% distin is the input MID. 
% option: signal reading option
function [x,err,ysim]=optcorr(mz,atoms,type,abd,imp,fwhm,dmz,distin,option)
    if isnan(sum(distin)) || sum(distin)==0  %input all zeros or nan, bypass        
       x=distin;
       err=0;
       ysim=distin;
       return
    end
    %dim=app.currentDim;
    %distin=str2num(app.EditFieldObserved.Value);
    distin=distin(:)';
    ind=find(distin>0); % removing trailing zeros;
    distin=distin(1:ind(end));
    %[mz,atoms]=str2mass(formula,charge);
    dim=atoms(type); % atom counts of type (C, N, or H)
     if length(distin)>dim+1
           distin=distin(1:dim+1);
     end
     y_exp=distin/sum(distin)*100; %normalize
     n=length(y_exp);
     x0=y_exp; % initialized x
            [~,maxM]=max(x0);
            lb=zeros(n,1);         
            ub=ones(n,1)*100;
            lb(maxM)=x0(maxM); %fix the value of the maxM in distin 
            ub(maxM)=x0(maxM);   
            options = optimoptions('fmincon','Display', 'off');
            x=fmincon(@(x)mycost(x,y_exp,mz,atoms,type,abd,imp,fwhm,dmz,option),x0,[],[],[],[],lb,ub,[],options);
            % % for GA optimization, very slow
            %options = optimoptions('ga','MaxGenerations',100,'PopulationSize',50,'PlotFcn', @gaplotbestf);
            %x=ga(@(x)mycost(x,y_exp,formula,charge,type,dmz),length(distin),[],[],[],[],lb,ub,[],options);            
            x=x/sum(x)*100;
            x(end+1:dim+1)=deal(0);
            x=x(:)'; %output x:  MID after correction
            [err, ysim]=mycost(x,y_exp,mz,atoms,type,abd,imp,fwhm,dmz,option); % output err            

end

