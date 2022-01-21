function [mass,atoms,flag]=str2mass(str,charge)
if nargin==1
    charge=0;
end
H=1.00728;
try
[mass,~,atoms]=formula2mass(str); %ordering: C,N,H,O,S
mass=mass+H*charge;
flag=1;
catch
    try 
      [mass,atoms]=pep2mass(str);
      flag=2;
    catch
       flag=0;
       mass=0;
       atoms=0;
       fprint('error') 
    end
end