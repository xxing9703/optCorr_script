%distin  input isotopmer distribution matrix, one sample per column
%atom1/atom2 is a structure containing:
% .num  %atom numbers in the formula
% .abd  %its natural isotope abudance
% .imp  %its impurity for labeled 
% takes 13C, 15N, 2H
function [distout,mm,distout_pct]=isocorr_due(distin,atom1,atom2)
Cnum=atom1.num; %atom number, abundance and impurity of atom1
Ca=atom1.abd;
Cb=atom1.imp;

if nargin==3
Nnum=atom2.num; %atom number, abundance and impurity of atom2
Na=atom2.abd;
Nb=atom2.imp;

elseif nargin==2
   Nnum=0;  
   Na=0.01;
   Nb=0.01;
end

%C matrix
dim=(Cnum+1)*(Nnum+1);
M=zeros(dim,dim);
mm=zeros(dim,dim,4);
mmp=zeros(dim,dim,4);
for i=1:dim
   C_num_row=floor((i-1)/(Nnum+1));
   N_num_row=mod((i-1),(Nnum+1));
    for j=1:i
         C_num_col=floor((j-1)/(Nnum+1));
         N_num_col=mod((j-1),(Nnum+1));
         
         a1=Cnum-C_num_col;
         b1=C_num_row-C_num_col;    
        
      
         a2=Nnum-N_num_col;
         b2=N_num_row-N_num_col;      
      mm(i,j,:)=[a1,b1,a2,b2];
      if b1>=0 && b2>=0     
        M(i,j)=dbinom(a1,b1,Ca)*dbinom(a2,b2,Na);
      end
    end
end

Mp=zeros(dim,dim);
for i=1:dim
   C_num_row=floor((i-1)/(Nnum+1));
   N_num_row=mod((i-1),(Nnum+1));
    for j=i:dim
         C_num_col=floor((j-1)/(Nnum+1));
         N_num_col=mod((j-1),(Nnum+1));
         
         a1=C_num_col;
         b1=C_num_col-C_num_row;    
        
      
         a2=N_num_col;
         b2=N_num_col-N_num_row;      
      mmp(i,j,:)=[a1,b1,a2,b2];
      if b1>=0 && b2>=0     
        Mp(i,j)=dbinom(a1,b1,Cb)*dbinom(a2,b2,Nb);
      end        

    end
end

N=M*Mp;
distout=distin;
 for i=1:size(distin,2)
     distout(:,i)=lsqnonneg(N,distin(:,i));%non-negative least square fitting
 end
 distout_pct=distout./(sum(distout,1)+1e-10);

% distout=max(distout,0); %nonzero
% distout=distout./sum(distout,2); 


function out=dbinom(a,b,r) %binormial distribution
out=nchoosek(a,b)*(1-r)^(a-b)*r^b;
