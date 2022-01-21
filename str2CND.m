%This function parses string of isotopeLabel(from Maven output) and
%finds out C N D labels. Deliminator '-' is used to make the judgement.
%examples:  str='C13-label-1', [Cnum,Nnum,Dnum]=[1,0,0]
%str='N15-label-1',[Cnum,Nnum,Dnum]=[0,1,0]
%str='C13N15-label-2-1',[Cnum,Nnum,Dnum]=[2,1,0]
%str='C12 PARENT', [Cnum,Num,Dnum]=[0,0,0]

function [Cnum,Nnum,Dnum,errmsg]=str2CND(str)
      Cnum=0; Nnum=0; Dnum=0;
      errmsg=0;
      sub_str=split(str,'-');
      len=length(sub_str);
      if len==1
      elseif len==3
          num=str2num(sub_str{end}); 
          if strcmp(str(1),'C')              
              Cnum=num;
          elseif strcmp(str(1),'N')
              Nnum=num;
          elseif strcmp(str(1),'D')
              Dnum=num;
          end
      elseif len==4          
            num1=str2num(sub_str{end-1});
            num2=str2num(sub_str{end});  
          if strcmp(str(1),'C')              
            Cnum=num1;
          elseif strcmp(str(1),'N')
            Nnum=num1;
          elseif strcmp(str(1),'D')
            Dnum=num1;
          end
          if contains(str(2:4),'C')
            Cnum=num2;
          elseif contains(str(2:4),'N')
            Nnum=num2;
          elseif contains(str(2:4),'D')
            Dnum=num2;
          end
      else
          fprintf('something is wrong with the string');
          errmsg=2;          
      end

