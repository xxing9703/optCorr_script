% isocorrection matrix for single or double tracer
% use get_isomatrix(n,ab1,impu1)  for single tracer
% n, m : atom number
% ab1, ab2:  natural abundance (i.e., 13C = 0.0110694, 15N =0.003663, 2D= 0.00015) 
% impu1, impu2: impurity (i,e, 0.01)
function [M,Mp]=get_isomatrix(n,ab1,impu1,m,ab2,impu2)
Ca=ab1;
Cb=impu1;
if nargin>3
Na=ab2;
Nb=impu2;
else
Na=0;
Nb=0;
m=0;
end

dim=(n+1)*(m+1);
M=zeros(dim,dim);
mm=zeros(dim,dim,4);
mmp=zeros(dim,dim,4);
for i=1:dim
   C_num_row=floor((i-1)/(m+1));
   N_num_row=mod((i-1),(m+1));
    for j=1:i
         C_num_col=floor((j-1)/(m+1));
         N_num_col=mod((j-1),(m+1));
         
         a1=n-C_num_col;
         b1=C_num_row-C_num_col;    
        
      
         a2=m-N_num_col;
         b2=N_num_row-N_num_col;      
      mm(i,j,:)=[a1,b1,a2,b2];
      if b1>=0 && b2>=0     
        M(i,j)=dbinom(a1,b1,Ca)*dbinom(a2,b2,Na);
      end
    end
end

Mp=zeros(dim,dim);
for i=1:dim
   C_num_row=floor((i-1)/(m+1));
   N_num_row=mod((i-1),(m+1));
    for j=i:dim
         C_num_col=floor((j-1)/(m+1));
         N_num_col=mod((j-1),(m+1));
         
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

function out=dbinom(a,b,r)
out=nchoosek(a,b)*(1-r)^(a-b)*r^b;