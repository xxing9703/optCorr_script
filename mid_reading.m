% this function reads signal at target locations up to M=maxM
% dmz defines mz tolerance window [-dmz, +dmz] 
% S.loc and S.pk are m/z & intensity pairs, can be either exp or simul
% if multiple peaks are detected within the tolerance window, 
% option 1:maximum; option 2: nearest 3. sum up

function [sig,mzerr]=mid_reading(type,maxM,dmz,S,option)
mzShift=[1.00335,0.99703,1.00630]; %mzshift for 13C, 15N, 2H
%type=1;
%maxM=6;
%ppm=5;
%mz=out(1).mz;
%dmz=mz*ppm/1e6;
locs=[0:maxM]*mzShift(type); % target m/z locations to extract signals from
locs_ob=[S.loc]; % observed m/z locations
pk_ob=[S.pk];

for i=1:length(locs)
   ind=find(abs(locs_ob(:)-locs(i))<dmz);
   if ~isempty(ind)
      switch option
          case 1
            [~,q]=max(pk_ob(ind));
            sig(i)=pk_ob(ind(q));
            mzerr(i)=locs_ob(ind(q))-locs(i); %mz err in da
          case 2
            [~,q]=min(abs(locs_ob(ind)-locs(i)));
            sig(i)=pk_ob(ind(q));
            mzerr(i)=locs_ob(ind(q))-locs(i);
          case 3
            sig(i)=sum(pk_ob(ind));
            mzerr(i)=nan;
      end
   else
       sig(i)=0;
       mzerr(i)=nan;
   end
end