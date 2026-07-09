clear all; close all;
addpath(genpath('topotoolbox-master/'),'-end');

%% Load DEM
DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');

%% Compute drainage network with Topotoolbox
G   = gradient8(DEM);Gdeg=gradient8(DEM,'degree');
DEMf = fillsinks(DEM);
FD = FLOWobj(DEMf);
A  = flowacc(FD);
Ariver=2e5; W = A>Ariver/(DEM.cellsize^2); 
Sriver = STREAMobj(FD,W);

%% Identify Landslides
phi  = 32;                                                           % Hillslope friction angle (degrees)
c    = 1e3;                                                          % Hillslope cohesion (Pa)  
rho  = 2600;                                                         % Hillslope density (kg/m3)
g    = 9.81;                                                         % Standard gravity (m/s2)
iopt = 2;                                                            % Keep only "2D landslides"
[LandslideMap,LandslideStats,LandslideStatsCC]=IdentifyLandslide(DEM,FD,Sriver,phi,c,rho,g,iopt);

%% Plot some classical landslide metrics
% Landslide label
figure; imageschs(DEM,LandslideMap.LabelCC_shuffled);title('Landslide label')
% Landslide depth
figure; imageschs(DEM,LandslideMap.Depth,'caxis',[0 prctile(LandslideMap.Depth,99.9,'all')]);title('Landslide depth (m)')
% Landslide volume-area distribution
figure; loglog(LandslideStats.Area,LandslideStats.Volume,'.');xlabel('landslide area (m^2)');ylabel('landslide volume (m^3)');title('Landslide volume-area')
% Landslide depth-area distribution
figure; loglog(LandslideStats.Area,LandslideStats.Depth,'.');xlabel('landslide area (m^2)');ylabel('landslide depth (m)');title('Landslide depth-area')
% PDF of landslide area
figure; histogram(LandslideStats.Area,logspace(log10(DEM.cellsize.^2),log10(max(LandslideStats.Area)),10),'Normalization','pdf');xscale log; yscale log;xlabel('landslide area (m^2)');ylabel('pdf (m^{-2})');title('Probability density function of landdslide area')
