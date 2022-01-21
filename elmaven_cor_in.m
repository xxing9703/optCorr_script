function [meta,A,start_col]=elmaven_cor_in(fname)
          A=readtable(fname,'readvariablename',true);
          A=A(1:length(find([A.medMz]>0)),:); %cut empty rows in case user edited csv. 
          start_col=find(strcmp(A.Properties.VariableNames,'parent'))+1;
          sample_names=A.Properties.VariableNames(start_col:end)';
          grpHead=find(strcmp(A.isotopeLabel,'C12 PARENT'));
          grpHead(end+1)=size(A,1)+1;
            % read formula, extract C/N/D number
         for i=1:length(grpHead)-1          
              ID(i)=i;          
              ids= grpHead(i):grpHead(i+1)-1;
              A.metaGroupId(ids)=ones(1,length(ids))*i;
              A_sub=A(ids,:); %A_sub: data sheet for the selected metabolite ID 
              try  %8/16/2020  added warning message for incorrect formulas
                [~,~,tp]=formula2mass(A_sub.formula{1});
              catch
                msgbox(['check row#',num2str(ids(1)+1),' for errors in the formula name: ',A_sub.formula{1}],'Error detected!');
                return
              end
              atom_num=tp(1:3);              
              
              lb=A_sub.isotopeLabel; %string analysis to get number of C,N,D
              counts=zeros(length(lb),2);
              for j=1:length(lb)
                 str=lb{j};
                 [Cnum,Nnum,Dnum,errmsg]=str2CND(str);
                 if errmsg==0
                   counts(j,1)=Cnum;
                   counts(j,2)=Nnum;
                   counts(j,3)=Dnum;               
                 else
                   msgbox('erros in isotopeLabel detected');
                   return
                 end
              end
              type=find(sum(counts)>0); % type can be 1,2,3 (C,N,H), single only
              distin = A_sub{:,start_col:end};
              distin_full=zeros(atom_num(type)+1,size(distin,2));
              for j=1:size(distin,1)
                  tp=counts(j,type)+1;
                  distin_full(tp,:)=distin(j,:);
              end
                meta(i).ID=ID(i);      
                meta(i).name=A_sub.compound{1};
                meta(i).formula=A_sub.formula{1};
                meta(i).mz=A_sub.medMz(1);
                meta(i).rt=A_sub.medRt(1);
                meta(i).type=type;
                meta(i).atom_num=atom_num;
                meta(i).original_abs=distin_full;
                meta(i).original_pct=distin_full./sum(distin_full,1);
                meta(i).tic=sum(distin_full,1);              
         end