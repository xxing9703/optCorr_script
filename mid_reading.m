% this function reads signal at target locations up to M=maxM
% dmz defines mz tolerance window [-dmz, +dmz] 
% S.loc and S.pk are m/z & intensity pairs, can be either exp or simul
% if multiple peaks are detected within the tolerance window, 
% option 1:maximum (observed); option 2: nearest(observed) 3. sum up (stick based)

function [sig,mzerr]=mid_reading(type,maxM,dmz,S,fwhm,option)
mzShift=[1.00335,0.99703,1.00630,2.004245]; %mzshift for 13C, 15N, 2H,18O
%type=1;
%maxM=6;
%ppm=5;
%mz=out(1).mz;
%dmz=mz*ppm/1e6;
locs_target=[0:maxM]*mzShift(type); % target m/z locations to extract signals from
if option==1 || option==2
 locs_obs=[S.loc]; % observed m/z locations from peakpicking of simulated profile
 pk_obs=[S.pk];
else
 locs_obs=[S.stickX]; % observed m/z locations from theoretical sticks
 pk_obs=[S.stickY];   
end

for i=1:length(locs_target)
   if option==1 || option==2
     ind=find(abs(locs_obs(:)-locs_target(i))<dmz);
   else
     ind=find(abs(locs_obs(:)-locs_target(i))<1.66*fwhm);    
   end
   if ~isempty(ind)
      switch option
          case 1
            [~,q]=max(pk_obs(ind));
            sig(i)=pk_obs(ind(q));
            mzerr(i)=locs_obs(ind(q))-locs_target(i); %mz err in da
          case 2
            [~,q]=min(abs(locs_obs(ind)-locs_target(i)));
            sig(i)=pk_obs(ind(q));
            mzerr(i)=locs_obs(ind(q))-locs_target(i);
          case 3
            sig(i)=sum(pk_obs(ind));
            mzerr(i)=nan;
      end
   else
       sig(i)=0;
       mzerr(i)=nan;
   end
end