function [out,tb,M]=mid_simulation(mass,atoms,maxM,abd,exts,types)
%##input
% atoms, a vector of atom_number C,N,H,O,S,[P],[F],Cl,Br
% maxM is the maximum M to calculate up to. default is 8(if not specified)
% abd: a vector of natural abudance consts for C,N,H
% ext and type are optional inputs, the combination of these two inputs are for labeled MID simulation. 
% type =1,2 or 3 representing (C, N, or H) 
% exts are the corresponding MID replacement from labeled 

%##output
% out is a structure containing the m/z and abundance of M+0, M+1,....
% tb is a structure containing the detailed info for each individual monoisotopic
% peaks (absolute abundance, name in str such as '13C2 15N 18O', m/z, dm/z, type:a flag used internally; pct: normalized to each M+i)  
% M is abundance for singly labeled only, but not combined.

%##isotopes considered in this simulation 
% 13C; 15N; 2H; 17O,18O, 33S,34S,36S; 37Cl; 81Br
if nargin==2
    maxM=8;
end
%----------------------------------------------------
n=7;
ty=logical([zeros(1,n*3)]); % 7 elements, 3 possible mass shift (+1, +2, +4)
str={'13C','15N','2H','17O','33S',    '',    '',...
       '' ,''   ,''  ,'18O','34S','37Cl','81Br',...
       '' ,''   ,''  ,''   ,'36S',''    ,    '',...
    };
atoms=atoms([1,2,3,4,5,8,9]); %select C,N,H,O,S,[P],[F],Cl,Br
%default abuandances
A0=[0.98893,0.996337,0.99985,0.9975904,0.9502,0.7576,0.5069]; %abundance for M0
A1=[0.0110694,0.003663,0.00015,3.7387E-04,0.0074962,0,0]; %abundance type A for +1
B1=[0,0,0,0.0020358,0.042099,0.2524,0.4931]; %abundance type B for +2
C1=[0,0,0,0,0.00020529,0,0]; %abundance type C for +4

if nargin>3  % use the abudances from input
    A1(1:3)=abd(1:3);    
    A0(1:3)=1-A1(1:3);
   
end

mzA=[1.00335,0.99703,1.00630,1.00422,0.99940,0,0]; %mz diff for +1
mzB=[0,0,0,2.00424,1.99580,1.997,1.998]; %mz diff for +2
mzC=[0,0,0,0,3.9950,0,0]; %mz diff for +4

%------------------------------------------

ab_0=prod(A0.^atoms(1:n)); %abundance for M0
% calculate abundance up to m
m=maxM;
for i=1:n
   for j=1:m
    ab_A(i,j)=ab_0*(A1(i)/A0(i))^j*nck(atoms(i),j);
    ab_B(i,j)=ab_0*(B1(i)/A0(i))^j*nck(atoms(i),j);
    ab_C(i,j)=ab_0*(C1(i)/A0(i))^j*nck(atoms(i),j);
    dmz_A(i,j)=mzA(i)*j;
    dmz_B(i,j)=mzB(i)*j;  
    dmz_C(i,j)=mzC(i)*j; 
   end
 end
M=[ab_A;ab_B;ab_C]; %MID of singly labeled(e.g, "33S" is single, but "13C 33S" is not single)
ct=0; % 
%------------ revise the MID of singly labeled based on: exts, types
if nargin>4
    if ~iscell(exts)
        tp{1}=exts;
        exts=tp;
    end
  for extNum=1:length(exts)
    ext=exts{extNum};
    type=types(extNum);
    ab_0_new=ab_0/prod(A0(type)^atoms(type))*(ext(1)/sum(ext)); %revise ab_0
    for k=1:size(M,1)
     for i=1:size(M,2)    
          M(k,i)=M(k,i)/ab_0*ab_0_new;  %revise M table by rescaling
     end
    end
    for i=1:length(ext)-1 
        if type<=3
          M(type,i)=ext(i+1)/ext(1)*ab_0_new; %modify type
        elseif type==4
          M(type+size(ab_A,1),i)=ext(i+1)/ext(1)*ab_0_new;  %updated 7/19/2024
        end
    end
    ab_0=ab_0_new;      
  end 
end

% ----------- write into tb
N=[dmz_A;dmz_B;dmz_C];
for i=1:size(M,1)
  for j=1:size(M,2)
      if M(i,j)>0
       ct=ct+1;
       tb(ct).ab=M(i,j);
       tb(ct).str=[str{i},num2str(j)];
       tb(ct).dmz=N(i,j); 
       tb(ct).mz=N(i,j)+mass; 
       %tb(ct).type=i;
       tb(ct).type=ty; 
       tb(ct).type(i)=true;
      end
  end
end

[~,ind]=sort([tb.ab],'descend'); %sorting
tb=tb(ind);
cutoff=1e-6;%apply cutoff
ind=find([tb.ab]>cutoff);  
tb=tb(ind);
%tb_save1=tb; % for debug
%% ----tb append to consider all combinations
ct=length(tb);
for i=1:length(tb)
    for j=i:length(tb)
        if ~sum(tb(i).type & tb(j).type)
            tp=tb(i).ab*tb(j).ab/ab_0;
          if tp>cutoff
              ct=ct+1;
              tb(ct).ab=tp;
               str_in=[tb(i).str,' ',tb(j).str];
              tb(ct).str=strjoin( sort(strsplit(str_in,' ')));
              tb(ct).dmz=tb(i).dmz+tb(j).dmz; 
              tb(ct).mz=tb(i).dmz+tb(j).dmz+mass; 
              tb(ct).type=tb(i).type|tb(j).type;
          end              
        end
    end
end
 [~,ind]=unique({tb.str});  %remove duplicates
 tb=tb(ind);

 tb=decouple(atoms,tb); %decouple, fix coeffients of exclusive items (e.g, "17O 18O")
 
 [~,ind]=find([tb.ab]>cutoff); %apply cutoff, remove ~zero items
 tb=tb(ind);
%% ---------- sorting option1 (not used)
 [~,ind]=sort([tb.ab],'descend'); %sorting by abudance
tb=tb(ind);
%tb_save2=tb; %for debug
%% -------------sorting option2 (used)
 [~,ind]=sort([tb.dmz],'ascend'); %sorting by dmz
tb=tb(ind);
%tb_save3=tb; %for debug

tb0.ab=ab_0;
tb0.str='Parent';
tb0.dmz=0;
tb0.mz=mass;
tb0.type=ty;
tb=[tb0,tb];
tb(1).pctGrp=100;
%% -------------- sum up for M+i and normalize
out(1).IsotopeNumber=0;
out(1).mz=mass;
out(1).pct=ab_0*100;
out(1).pctMax=100;
for i=1:round(max([tb.dmz]))
  out(i+1).IsotopeNumber=i;
  [~,ind]=find(round([tb.dmz])==i);
  out(i+1).mz=sum([tb(ind).mz].*[tb(ind).ab])/sum([tb(ind).ab]);
  out(i+1).pct=sum([tb(ind).ab])*100;
  for j=1:length(ind)
    tb(ind(j)).pctGrp=tb(ind(j)).ab/max([tb(ind).ab])*100;
  end  
end
for i=1:round(max([tb.dmz]))+1
   out(i).pctMax=out(i).pct/max([out.pct])*100;
end
end

function flag = hasNoRepeats(A)
    flag = true;
    A = sort(A);
    for i=1:(numel(A)-1)
        if A(i)==A(i+1)
            flag = false;
            break;
        end
    end
end
%------------------------------------------------------------------------
%if a combination contains items from the same element, e.g., "17O 18O2"
%the factor needs to be revised e.g.,num_O=6, nck(6,1)*nck(6,2) --> nck(6,1)*nck(5,2)
%this applies to [17O,18O]  and [33S 34S 36S]
function tb=decouple(atoms,tb)
for i=1:length(tb)
  item=tb(i);
  tp=strfind(item.str,'O');
  if length(tp)==2
    n1=str2num(item.str(tp(1)+1));
    n2=str2num(item.str(tp(2)+1));
    nS=atoms(4);
    factor=nck(nS,n1)*nck(nS-n1,n2)/(nck(nS,n1)*nck(nS,n2));
    item.ab=item.ab*factor;    
  end
  
  tp=strfind(item.str,'S');
  if length(tp)==2
    n1=str2num(item.str(tp(1)+1));
    n2=str2num(item.str(tp(2)+1));
    nS=atoms(5);
    factor=nck(nS,n1)*nck(nS-n1,n2)/(nck(nS,n1)*nck(nS,n2));
    item.ab=item.ab*factor;  
  elseif length(tp)==3
    n1=str2num(item.str(tp(1)+1));
    n2=str2num(item.str(tp(2)+1));
    n3=str2num(item.str(tp(3)+1));
    nS=atoms(5);
    factor=nck(nS,n1)*nck(nS-n1,n2)*nck(nS-n1-n2,n3)/(nck(nS,n1)*nck(nS,n2)*nck(nS,n3));
    item.ab=item.ab*factor;      
  end
  tb(i)=item;
end
end
% ----------------------------------------
% modified function of nchoosek, accepting n<k which returns 0
function out=nck(n,k)
if n<k 
    out=0;
else
    out=nchoosek(n,k);
end
end
