function A = D8DmultiSimil(FDsingle,FDinf)
%A = D8DmultiSimil(FDsingle,FDinf)
% Computes the similarity betwwen the D8 and multiple flow algorithm
% (developed by Qin et al, 2017)
%INPUT:
%- FDsingle : FLOWobj computed the 'D8' algorithm
%- FDinf : FLOWobj computed using the 'multi' algorihtm 
%OUTPUT:
%A : vector containing the same number of elements as the flow direction
%objects, indicating for each pixel the fraction of the flow that is
%distributed to the lowest neighbouring pixel in the mi=ultiple flow
%algorithm

ix       = FDinf.ix;
fraction = FDinf.fraction;
A        = zeros(FDsingle.size);
[X,I]    = sort(ix);
fraction = fraction(I);
inode    = 1;
imax     = inode;
while  inode<numel(X)
    while inode<numel(X) && X(inode+1)== X(inode) 
        if fraction(inode)>fraction(imax)
            imax  = inode+1;
            inode = inode+1;
        else
            inode = inode+1;
        end
        
    end
    A(X(inode)) = fraction(imax);
    inode       = inode+1;
    imax        = inode;
end
            
end
