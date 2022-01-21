function out=propagate(in,dim,M,Mp)
   in=in(:)';  % make it a row vector
   if length(in)>dim+1  %truncate
      in=in(1:dim+1);
   elseif length(in)<dim+1  % padding 0
      in(end+1:dim+1)=deal(0);
   end
   out=M*Mp*in';  % propagate
   %out=out(1:min(8,length(out)));
  