function [head] = CO2CHAIN(DEM)
% CO2CHAIN retrieves channel heads in a basin following the CO2CHAIN algorithm  developped in Lurin et al, 2023 
%
% input : DEM : Topotoolbox GRIDobj for the DEM on which channel heads will
% be extracted. The resolution should be 3m or 10m
% output: head: list of indexes of the channel heads found by CO2CHAIN
%
% The heads will be extracted from the largest basin in the DEM. If you
% want to extract all heads regardless of the basin, replace lines 29-30
% by:
% DB = DEM;
% DB.Z(:) = 1;
%
% Parametrization : 
% The default best parameters are 
% -U_c_thresh = 12
% -L_c_thresh = 0.45
% They can be changed on line 53-54




%% Preliminary computations

res             = DEM.cellsize;

%flow direction and accumulation
FD              = FLOWobj(DEM);
A               = flowacc(FD);
A.Z             = A.Z*A.cellsize^2;
FDinf           = FLOWobj(DEM,'dinf');
Ainf            = flowacc(FDinf);
Ainf.Z          = Ainf.Z*Ainf.cellsize^2;
FDmulti         = FLOWobj(DEM,'multi');
Amulti          = flowacc(FDmulti);
Amulti.Z        = Amulti.Z*Amulti.cellsize^2;
FDsimil_multi   = D8DmultiSimil(FD,FDmulti);
FD.fastindexing = 1;
G               = gradient8(DEM);

%extract main basin                            
outlet  = find(A.Z==max(A.Z(:)));
DB      = drainagebasins(FD,outlet);

%extract crests
DEM2    = filter(DEM,'mean',floor(100/res));
tpi     = DEM-DEM2;
tpimask = tpi.Z>0.5;
Amask   = (A.Z<=1*A.cellsize^2) .* double(DB.Z).*tpimask;
crests  = find(Amask==1);
clear('DEM2','tpi','tpimask','Amask','Amulti','FDmulti')

%% Identification of hillslpes/channel transition

%prepare loop
L_c_thresh = 0.45;
U_c_thresh = 12;
heads      = zeros(size(crests));
gthresh    = 0.45.*median(G.Z(DB.Z==1),'omitnan');
zthresh    = prctile(DEM.Z(:),20);

%Ignore one out of two crest pixels for efficiency on three-meter-resolution DEMs
if res == 3
    step = 2;
else
    step = 1;
end

%Find hillslope/channel transition for each crest pixel
for ic = 1:step:numel(crests)
    
%extract flow path and compute distance and area for flow paths above 400m
[path,dist]     = flowpathextract(FD,crests(ic));
if numel(path) >= floor(400/res) 
dist(1)         = 0.1;
path            =  path(1:floor(400/res));
dist            =  dist(1:floor(400/res));
a               =  Ainf.Z(path);

%upstream convergence differential
a_L              = a./dist;
L_c              = FDsimil_multi(path);
[dist_bin,a_bin] = logbin_sumred(dist(2:end),a(2:end),3,0,30); 
[dist_bin,U_c]   = logbin_sumred(dist(2:end),a_L(2:end),3,0,30); 
dU_c             = diff(U_c);
L_c_env          = localmaxeenveloppe(L_c,5); 

%Find transition
idx                 = 4;
[trans_Uc,found_Uc] = findjump(dU_c,U_c_thresh,idx); %all significative Uc step
trans_Lc            = find(L_c_env >= L_c_thresh);

%find linear indexes
i_Lc = trans_Lc;
i_Uc = zeros(1,numel(trans_Uc));
if found_Uc
   for j = 1:numel(trans_Uc) 
       [~,i_Uc(j)] = min(abs(dist-dist_bin(trans_Uc(j))));
   end
end

%Select transition
ifound = 0;
if ~isempty(i_Lc) && found_Uc
    ifound = 0;
    j      = 0;
    while ifound ==0 && j<min(numel(i_Uc),numel(i_Lc))
      j          = j+1;
      [ecart,jj] = min(abs(i_Uc(j)-i_Lc));
      ifound     = ecart<=floor(15/res);
      i_head     = max(i_Uc(j),i_Lc(jj));
    end

    %Store channel head localization
    if ifound
    dist_lims    = [(dist_bin(trans_Uc(j))+dist_bin(trans_Uc(j)-1))/2,dist_bin(trans_Uc(j)+1)];
    [~,minlin]   = min(abs(dist-dist_lims(1)));
    [~,maxlin]   = min(abs(dist-dist_lims(2)));
    jumps        = diff(a_L(minlin:maxlin));
    [~,ilin]     = max(jumps);
    i_head       = minlin+ilin;
    heads(ic)    = path(i_head);
    else
        heads(ic) = 0;
    end
end
end
end


%% Drainage network correction
 
heads = heads(heads>0);

% raw drainage network

% 1st correct drainage network
S1 = STREAMobj(FD,'minarea',1);
S1 = klargestconncomps(S1,1);
S  = modify(S1,'downstreamto',heads);
%find heads
segs    = networksegment(S,FD,DEM,A,0.39);
order1  = find(segs.strahler==1);
head    = segs.IX(order1,2);
g       = segs.fslope(order1);
z       = DEM.Z(head);
suppr   = (g<gthresh).*(z<zthresh);
head    = head(suppr==0);

%Remove small branches
S2      = modify(S1,'downstreamto',head);
segs    = networksegment(S2,FD,DEM,A,0.39);
rm      = ((segs.strahler==1).*(segs.length<=28.3));
suppr   = segs.ix(rm==1);
eal     = zeros(size(S2.ix));
for iS = 1:numel(suppr)
    del      = find(S2.ix == suppr(iS));
    eal(del) =1;
end
S3      = rmedge(S2,eal==1);
S3      = klargestconncomps(S3,1);

%remove parallels
maskS            = zeros(size(DEM.Z));
maskS(S3.IXgrid) = 1;
segs             = networksegment(S3,FD,DEM,A,0.39);
order1           = find(segs.strahler==1);
eal              = zeros(size(S3.ix));
for jj=1:numel(order1)
    phead = segs.IX(order1(jj),2);
    [x,y] = ind2sub(size(DEM.Z),phead);
    try
    Spoints = sum(sum(maskS(x-2:x+2,y-2:y+2)));
    if Spoints>5
        seg      = segs.ix(order1(jj));
        seg      = find(S3.ix==seg);
        eal(seg) = 1;
    end
    catch
    end
    
end

S4     = rmedge(S3,eal==1);
S4     = klargestconncomps(S4,1);

%one more round of small segments deletion
segs   = networksegment(S4,FD,DEM,A,0.39);
order1 = find(segs.strahler==1);
rm     = ((segs.strahler==1).*(segs.length<=28.3));
suppr  = segs.ix(rm==1);
eal = zeros(size(S4.ix));
for iS = 1:numel(suppr)
    del = find(S4.ix == suppr(iS));
    eal(del)=1;
end
S = rmedge(S4,eal==1);
S = klargestconncomps(S,1);
segs = networksegment(S,FD,DEM,A,0.39);
order1 = find(segs.strahler==1);
head = segs.IX(order1,2);

end
