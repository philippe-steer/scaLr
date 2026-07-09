function [LandslideMap,LandslideStats,LandslideStatsCC]=IdentifyLandslide(DEM,FD,Sriver,phi,c,rho,g,iopt)
% Identify landslides using a 1D hillslope safety factor criterion.
% Based on Steer et al. (2026) and on Jeandet et al. (2019) but
% extrapolated to a DEM (not only 1D). Landslides are identified as
% unstable lines that connect to the same outlet point. 
% This algorithm checks stability for each cell of the DEM and extract the
% most unstable rupture line between this point and a downstream point (in
% an hydrological sense). This is associated to an optimal depth.
%
% --- Inputs:
% DEM    - DEM GRIDobj as produced by Topotoolbox
% FD     - Flow direction FLOWobj as produced by Topotoolbox
% Sriver - River STREAMobj computed using for instance a minimum drainage area
% phi    - Hillslope friction angle (degrees)
% c      - Hillslope cohesion (Pa)  
% rho    - Hillslope density (kg/m3)
% g      - Standard gravity (m/s2)
% iopt   - if opt==1 all landslides are kept (this is useful to compute
%          unstable areas and depths)
%        - if opt==2 only 2D landslides are kept (landslides with more than
%          one hydrologic path)
%
% --- Outputs:
% LandslideMap   - Structure that contains information on landslide in a map ready format
% LandslideStats - Structure that contains statistical information on landslides% IDENTIFYLANDSLIDE_FAST  Faster, corrected replacement for IdentifyLandslide.
%
% Same inputs and outputs as the original IdentifyLandslide.
%
% Speed:
%   1. flowpathextract(FD2,i) is called once per DEM cell in the original
%      (millions of TopoToolbox calls). Here the single-flow-direction
%      receiver map is built ONCE from FD2.ix / FD2.ixc, and each downstream
%      (hydrologic) path is rebuilt by following that map with plain array
%      indexing. This reproduces flowpathextract(FD2,i) exactly (source
%      dropped, path to the outlet).
%   2. tand(phi) is hoisted out of the loop.
%   3. The independent per-source work is run with parfor (falls back to a
%      serial loop if no parallel pool is available).
%
% Correction (rupture-plane selection):
%   The optimal-depth condition dF/dh = 0 is a quadratic in h, so there are
%   TWO roots (the -/+ square-root branches). On the admissible interval
%   (0, Dz), F(h) -> +Inf at both ends, so it has a single interior minimum:
%   exactly one of the two roots lies in (0, Dz) and is the physical
%   minimiser (not necessarily the '-' root). The radicand becomes negative
%   only when a tested downstream node lies above the source (Dz < 0), a
%   sink-filling artifact with no physical failure plane; both roots are then
%   complex conjugates and the node is discarded.
%
%   For each downstream node this function therefore evaluates BOTH roots,
%   keeps the real one in (0, Dz), and selects, over the whole path, the
%   admissible plane with the minimum safety factor. All outputs are real.
%
% See IdentifyLandslide.m for the full input/output documentation.

% --- Crop the flowdir obj to hillslopes only (remove the river network)
GS  = STREAMobj2GRIDobj(Sriver); GS = ~GS;
FD2 = crop(FD,GS);                          % Flow dir without the fluvial network
FD1 = FD;  FD1.fastindexing = true;         % Flow dir with the fluvial network (used by iopt==2)

% --- Extract DEM coordinates
[X,Y] = getcoordinates(DEM,'matrix'); siz = size(DEM.Z); N = numel(DEM.Z);
x = X(:); y = Y(:); z = DEM.Z(:);

% --- Precompute the single-flow-direction receiver map for the hillslope network.
%     R(p) = downstream neighbour of p in FD2 (0 if p drains out of the hillslope
%     network, e.g. into the river). Following R reproduces flowpathextract(FD2,p).
R = zeros(N,1);
R(FD2.ix) = FD2.ixc;

% --- Precompute the number of downstream hillslope nodes of each cell
%     (= length of flowpathextract(FD2,p) minus 1). L(p) = L(R(p)) + 1.
%     Filled downstream-first, i.e. reverse of the topologically sorted edge list.
L   = zeros(N,1);
ixh = FD2.ix;  ixch = FD2.ixc;
for r = numel(ixh):-1:1
    L(ixh(r)) = L(ixch(r)) + 1;
end

% --- Hoist constant out of the loop
tanphi = tand(phi);

% --- Source cells (non-NaN). Results are stored per-source, then scattered onto the grid.
ind = find(~isnan(DEM.Z));
ns  = numel(ind);

Failure_s = zeros(ns,1);
h_Fmin_s  = zeros(ns,1);
Fmin_s    = zeros(ns,1);
iout_s    = zeros(ns,1);
Length_s  = zeros(ns,1);
Stopo_s   = zeros(ns,1);
Srup_s    = zeros(ns,1);

parfor j = 1:ns
    i = ind(j);
    n = L(i);                                  % number of downstream nodes
    if n > 0
        % --- Rebuild the downstream (hydrologic) path without flowpathextract.
        %     idown = [R(i); R(R(i)); ... ; outlet], i.e. flowpathextract(FD2,i)(2:end).
        idown = zeros(n,1);
        p = R(i);
        for t = 1:n
            idown(t) = p;
            p = R(p);
        end

        xi = x(i); yi = y(i); zi = z(i);
        Dx = sqrt( (x(idown)-xi).^2 + (y(idown)-yi).^2 );   % horizontal distances
        Dz = zi - z(idown);                                 % vertical distances

        % --- Both roots of dF/dh = 0 (optimal rupture depth).
        sq   = (c.*(c + Dz.*g.*tanphi.*rho).*(Dx.^2 + Dz.^2)).^(1./2);  % complex iff Dz<0
        den  = c.*Dz - g.*tanphi.*rho.*Dx.^2;
        num0 = c.*Dx.^2 + c.*Dz.^2;
        hm   = (num0 - Dx.*sq)./den;           % '-' root (original h_Fmin_vec)
        hp   = (num0 + Dx.*sq)./den;           % '+' root

        % --- Safety factor at each root (real part; invalid roots masked to Inf below)
        hmr = real(hm);
        Fm  = c./(rho.*g.*hmr).*(Dx.^2+(Dz-hmr).^2)./(Dx.*(Dz-hmr)) + tanphi./(Dz./Dx - hmr./Dx);
        hpr = real(hp);
        Fp  = c./(rho.*g.*hpr).*(Dx.^2+(Dz-hpr).^2)./(Dx.*(Dz-hpr)) + tanphi./(Dz./Dx - hpr./Dx);

        % --- Keep only physically admissible planes: node below the source,
        %     real depth, 0 < h < Dz.
        okm = (Dz > 0) & (imag(hm) == 0) & (hmr > 0) & (hmr < Dz);
        okp = (Dz > 0) & (imag(hp) == 0) & (hpr > 0) & (hpr < Dz);
        Fm(~okm) = Inf;
        Fp(~okp) = Inf;

        % --- Per node, pick the admissible root; over the path, pick the min safety factor
        useplus = Fp < Fm;
        Fnode   = min(Fm,Fp);
        hnode   = hmr;  hnode(useplus) = hpr(useplus);

        [Fsel,k] = min(Fnode);
        if isfinite(Fsel)
            Failure_s(j) = 1;                             % this point can be a failure
            h_Fmin_s(j)  = hnode(k);                      % potential rupture depth
            Fmin_s(j)    = Fsel;                          % safety factor (real, > 0)
            Length_s(j)  = Dx(k);                         % horizontal length of rupture plane
            Stopo_s(j)   = Dz(k)./Dx(k);                  % slope of topographic profile
            Srup_s(j)    = (Dz(k)-hnode(k))./Dx(k);       % slope of rupture plane
            iout_s(j)    = idown(k);                      % node where the rupture daylights
        end
    end
end

% --- Scatter per-source results back onto the full grid (column vectors, as in the original)
Failure = zeros(N,1); Failure(ind) = Failure_s;
h_Fmin  = zeros(N,1); h_Fmin(ind)  = h_Fmin_s;
Fmin    = zeros(N,1); Fmin(ind)    = Fmin_s;
iout    = zeros(N,1); iout(ind)    = iout_s;
Length  = zeros(N,1); Length(ind)  = Length_s;
Stopo   = zeros(N,1); Stopo(ind)   = Stopo_s;
Srup    = zeros(N,1); Srup(ind)    = Srup_s;

% Add outlets
ind = find(iout>0); iout(iout(ind)) = iout(ind);

if iopt==2
    % Identify and remove outlets that are associated to a unique unstable
    % hydrologic path (motivation: landslides should be 2D)
    iout_map=reshape(iout,siz);[iout_unique,~,ic]=unique(iout);n=zeros(size(iout_unique));
    [II,JJ] = ind2sub(siz,iout_unique);
    % Extract flow direction
    ixcix = FD1.ixcix;ixcix(ixcix==0)=min(setdiff(1:numel(DEM.Z),iout_unique)); % Remove some spurious zeros that appear on the side of ixcix
    ix   = FD1.ix;
    for i=2:numel(iout_unique)
        % Extract iout direction of the direct neighbous to the central point (iout outlet)
        temp1=iout_map(max(1,II(i)-1):min(II(i)+1,siz(1)),max(1,JJ(i)-1):min(JJ(i)+1,siz(2)));
        % Extract flow givers of the direct neighbors
        temp2=zeros(size(temp1));temp2(temp1~=0)=ix(ixcix(temp1(temp1~=0)));
        % Remove the central node (the outlet iout)
        temp1(2,2)=0; temp2(2,2)=0;
        % Count the umbers of direct neighbours that flow towards the
        % outlet and being unstable towards him;
        % n(i)=numel(find(temp1==iout_unique(i) & temp2==iout_unique(i)));
        n(i)=numel(find(temp1==iout_unique(i)));
    end
    % Remove iout nodes that does not respect the condition of being
    % associated to at least two unstable hydrologic paths (different from
    % each other)
    nn=ones(size(iout));nn(n<2)=0;
    iout=iout.*nn(ic);
end

% Generate landslide map information
ind=find(Failure==1 & Fmin<=1 & Fmin>0 & iout>0);
iout_cleaned=zeros(size(iout));iout_cleaned(ind)=iout(ind);[iout_cleaned,~,ic] = unique(iout_cleaned);nl=1:numel(iout_cleaned);iout_cleaned=nl(ic);
ind=find(Failure==0 | Fmin>1 | Fmin<=0 | iout==0 | imag(Fmin)>0);
LandslideMap.Label=reshape(iout_cleaned,siz);   LandslideMap.Label(ind)=NaN;
LandslideMap.IdxOut=reshape(iout,siz);          LandslideMap.IdxOut(ind)=NaN;
LandslideMap.Failure=reshape(Failure,siz);      LandslideMap.Failure(ind)=0;
LandslideMap.Depth=reshape(h_Fmin,siz);         LandslideMap.Depth(ind)=NaN;
LandslideMap.Fmin=reshape(Fmin,siz);            LandslideMap.Fmin(ind)=NaN;
LandslideMap.Srup=reshape(Srup,siz);            LandslideMap.Srup(ind)=NaN;
LandslideMap.Stopo=reshape(Stopo,siz);          LandslideMap.Stopo(ind)=NaN;
LandslideMap.Length=reshape(Length,siz);        LandslideMap.Length(ind)=NaN;
LandslideMap.Label_shuffled = shufflelabel_simple(LandslideMap.Label);
LandslideMap.IdxOut_shuffled = shufflelabel_simple(LandslideMap.IdxOut);

% Generate a connected compnent catalog (to assess amalgamation compared to observations)
temp=LandslideMap.Label;temp(temp>0)=1;temp(isnan(temp))=0;LandslideMap.LabelCC = bwlabel(temp);LandslideMap.LabelCC(LandslideMap.LabelCC==0)=NaN;
LandslideMap.LabelCC_shuffled = shufflelabel_simple(LandslideMap.LabelCC);

% Identify landslides and compute their geometrical descriptors (original
% catalog)
tempLabel=LandslideMap.Label; tempLabel(ind)=-1;
stats = regionprops(tempLabel,'PixelList','PixelIdxList','centroid','Area','Perimeter','Orientation','MinorAxisLength','MajorAxisLength','Eccentricity','Circularity');
LandslideStats.n=size(stats,1);
for i=1:LandslideStats.n;    LandslideStats.Idx{i}=stats(i).PixelIdxList; end
LandslideStats.Area=[stats.Area].*DEM.cellsize.^2;
for i=1:LandslideStats.n;    LandslideStats.Depth(1,i)=mean(LandslideMap.Depth(stats(i).PixelIdxList));end
LandslideStats.Volume=LandslideStats.Area.*LandslideStats.Depth;
LandslideStats.Perimeter=[stats.Perimeter].*DEM.cellsize;
LandslideStats.MinorAxisLength=[stats.MinorAxisLength];
LandslideStats.MajorAxisLength=[stats.MajorAxisLength];
LandslideStats.LWratio=[stats.MajorAxisLength]./[stats.MinorAxisLength];
temp=reshape([stats.Centroid],2,LandslideStats.n);temp=sub2ind(siz, temp(2,:), temp(1,:));
LandslideStats.CentroidIdx=temp;

% Identify landslides and compute their geometrical descriptors (connected
% component catalog)
tempLabel=LandslideMap.LabelCC; tempLabel(ind)=-1;
stats = regionprops(tempLabel,'PixelList','PixelIdxList','centroid','Area','Perimeter','Orientation','MinorAxisLength','MajorAxisLength','Eccentricity','Circularity');
LandslideStatsCC.n=size(stats,1);
for i=1:LandslideStatsCC.n;    LandslideStatsCC.Idx{i}=stats(i).PixelIdxList; end
LandslideStatsCC.Area=[stats.Area].*DEM.cellsize.^2;
for i=1:LandslideStatsCC.n;    LandslideStatsCC.Depth(1,i)=mean(LandslideMap.Depth(stats(i).PixelIdxList));end
LandslideStatsCC.Volume=LandslideStatsCC.Area.*LandslideStatsCC.Depth;
LandslideStatsCC.Perimeter=[stats.Perimeter].*DEM.cellsize;
LandslideStatsCC.MinorAxisLength=[stats.MinorAxisLength];
LandslideStatsCC.MajorAxisLength=[stats.MajorAxisLength];
LandslideStatsCC.LWratio=[stats.MajorAxisLength]./[stats.MinorAxisLength];
temp=reshape([stats.Centroid],2,LandslideStatsCC.n);temp=sub2ind(siz, temp(2,:), temp(1,:));
LandslideStatsCC.CentroidIdx=temp;
