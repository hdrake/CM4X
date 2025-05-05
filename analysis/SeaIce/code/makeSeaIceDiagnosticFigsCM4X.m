%%this function creates a standard set of sea ice diagnostic figures

function makeSeaIceDiagnosticFigsCM4X(model,run,year_start,year_end,year_shift,short_name,SIC_clim_years,SIT_clim_years,vel_clim_years,isPI,loadFigData) 

%note: year_shift is typically zero (for historical runs, etc), but can be used to line up a control run (which starts from year 1) with an appropriate obs period.

%addpath('~/scripts/export_fig/')
addpath('~/scripts/export_fig_2021/')
addpath('~/scripts/utils/')
addpath('~/scripts/m_map/')
addpath('~/scripts/')

close all

N=length(model)

compare_string=short_name{1};
for k=2:N
compare_string=[compare_string,'_',short_name{k}];
end

compare_string

label_str={'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)'};

fig_dir=['/home/mib/Figures/compare_sea_ice/',model{1},'/',compare_string];
if(~isdir(fig_dir))
mkdir(fig_dir)
end
work_dir=['/work/mib/processed_data/compare_sea_ice/',model{1},'/',compare_string];
if(~isdir(work_dir))
mkdir(work_dir)
end

modelDir=['/work/mib/raw_data/',model{1}];
load([modelDir,'/ice_grid_regions.mat'],'ifXY_reg','reg_list')

%%%%%%load monthly model data and compute sea ice diagnostics

%%define minimum year
plot_year_min=1978;
plot_year_max=2020;
for k=1:N
time_years{k}=year_start{k}+year_shift{k}:1:year_end{k}+year_shift{k};

tmp=time_years{k};
plot_year_min=min([plot_year_min min(tmp)])
plot_year_max=max([plot_year_max max(tmp)])

end

obs_year_end=2023;

time_obs=1979+1/24:1/12:obs_year_end+1-1/24;
obs_years=1979:obs_year_end;
PIOMAS_years=1979:obs_year_end;

if(loadFigData==0)

for k=1:N

MODEL{k}=struct;

%time{k}=year_start{k}+1/24:1/12:year_end{k}+1-1/24;
time_years{k}=year_start{k}+year_shift{k}:1:year_end{k}+year_shift{k};

dataDir=['/work/mib/raw_data/',model{k},'/',run{k}];

%%NH data
%var=['siconc_monthly'];
var=['siconc'];
filename=[dataDir,'/',var,'.mat']
load(filename,'data','time','x','y','x_c','y_c','ifXY','w')

size(data)

if(exist('x_c')==1)
MODEL{k}.xC_nh=x_c;
MODEL{k}.yC_nh=y_c;
else
MODEL{k}.xC_nh=x;
MODEL{k}.yC_nh=y;
end
MODEL{k}.x_nh=x;
MODEL{k}.y_nh=y;
MODEL{k}.w_nh=w;
MODEL{k}.ifXY_nh=ifXY;
if(max(max(data))>=99)
MODEL{k}.SIC_nh=data/100;
else
MODEL{k}.SIC_nh=data;
end


%var=['sithick_monthly'];
var=['sithick'];
filename=[dataDir,'/',var,'.mat']
load(filename,'data')

size(data)

MODEL{k}.SIT_nh=data;

[MODEL{k}.SIC_clim_nh MODEL{k}.SIC_std_nh MODEL{k}.SIT_clim_nh MODEL{k}.SIT_std_nh MODEL{k}.extent_nh MODEL{k}.extent_clim_nh MODEL{k}.extent_anom_nh MODEL{k}.volume_nh MODEL{k}.volume_clim_nh MODEL{k}.volume_anom_nh]=computeSeaIceDiagnostics(MODEL{k}.SIC_nh,MODEL{k}.SIT_nh,MODEL{k}.w_nh);

%compute clim over SIC satellite era
ind1=(SIC_clim_years(1)-year_start{k}-year_shift{k})*12+1
ind2=(SIC_clim_years(2)+1-year_start{k}-year_shift{k})*12

[MODEL{k}.SIC_clim_nh_common MODEL{k}.SIC_std_nh_common MODEL{k}.SIT_clim_nh_common MODEL{k}.SIT_std_nh_common MODEL{k}.extent_nh_common MODEL{k}.extent_clim_nh_common MODEL{k}.extent_anom_nh_common MODEL{k}.volume_nh_common MODEL{k}.volume_clim_nh_common MODEL{k}.volume_anom_nh_common]=computeSeaIceDiagnostics(MODEL{k}.SIC_nh(:,ind1:ind2),MODEL{k}.SIT_nh(:,ind1:ind2),MODEL{k}.w_nh);

%compute SIV/SIA regional timeseries and correlations over SIC satellite era

modelDir=['/work/mib/raw_data/',model{k}];
load([modelDir,'/ice_grid_regions.mat'],'ifXY_reg','reg_list')

MODEL{k}.ifXY_reg_nh=ifXY_reg;

[MODEL{k}.regionalSIA MODEL{k}.regionalSIAclim]=computeRegionalSIA(MODEL{k}.SIC_nh(:,ind1:ind2),MODEL{k}.w_nh,MODEL{k}.ifXY_nh,MODEL{k}.ifXY_reg_nh);
[MODEL{k}.regionalSIV MODEL{k}.regionalSIVclim]=computeRegionalSIV(MODEL{k}.SIC_nh(:,ind1:ind2),MODEL{k}.SIT_nh(:,ind1:ind2),MODEL{k}.w_nh,MODEL{k}.ifXY_nh,MODEL{k}.ifXY_reg_nh);


%compute clim over SIT satellite era
ind1=(SIT_clim_years(1)-year_start{k}-year_shift{k})*12+1;
ind2=(SIT_clim_years(2)+1-year_start{k}-year_shift{k})*12;
[MODEL{k}.SIT_clim_nh_common SIT_anom]=computeclim(MODEL{k}.SIT_nh(:,ind1:ind2));


%var=['siu_monthly'];
var=['siu_trueE'];
filename=[dataDir,'/',var,'.mat']
load(filename,'data')

size(data)

MODEL{k}.siu_nh=100*data;

%var=['siv_monthly'];
var=['siv_trueN'];
filename=[dataDir,'/',var,'.mat']
load(filename,'data')

MODEL{k}.siv_nh=100*data;

size(data)

%compute clim over OSISAF satellite era
ind1=(vel_clim_years(1)-year_start{k}-year_shift{k})*12+1;
ind2=(vel_clim_years(2)+1-year_start{k}-year_shift{k})*12;
[MODEL{k}.siu_clim_nh_common siu_anom]=computeclim(MODEL{k}.siu_nh(:,ind1:ind2));
[MODEL{k}.siv_clim_nh_common siv_anom]=computeclim(MODEL{k}.siv_nh(:,ind1:ind2));


%%SH data
%var=['siconc_monthly_sh'];
var=['siconc_sh'];
filename=[dataDir,'/',var,'.mat']
load(filename,'data','time','x','y','x_c','y_c','ifXY','w')

if(exist('x_c')==1)
MODEL{k}.xC_sh=x_c;
MODEL{k}.yC_sh=y_c;
else
MODEL{k}.xC_sh=x;
MODEL{k}.yC_sh=y;
end
MODEL{k}.x_sh=x;
MODEL{k}.y_sh=y;
MODEL{k}.w_sh=w;
MODEL{k}.ifXY_sh=ifXY;
if(max(max(data))>=99)
MODEL{k}.SIC_sh=data/100;
else
MODEL{k}.SIC_sh=data;
end

%var=['sithick_monthly_sh'];
var=['sithick_sh'];
filename=[dataDir,'/',var,'.mat']
load(filename,'data')

MODEL{k}.SIT_sh=data;

[MODEL{k}.SIC_clim_sh MODEL{k}.SIC_std_sh MODEL{k}.SIT_clim_sh MODEL{k}.SIT_std_sh MODEL{k}.extent_sh MODEL{k}.extent_clim_sh MODEL{k}.extent_anom_sh MODEL{k}.volume_sh MODEL{k}.volume_clim_sh MODEL{k}.volume_anom_sh]=computeSeaIceDiagnostics(MODEL{k}.SIC_sh,MODEL{k}.SIT_sh,MODEL{k}.w_sh);

%compute clim over SIC satellite era
ind1=(SIC_clim_years(1)-year_start{k}-year_shift{k})*12+1;
ind2=(SIC_clim_years(2)+1-year_start{k}-year_shift{k})*12;
[MODEL{k}.SIC_clim_sh_common MODEL{k}.SIC_std_sh_common MODEL{k}.SIT_clim_sh_common MODEL{k}.SIT_std_sh_common MODEL{k}.extent_sh_common MODEL{k}.extent_clim_sh_common MODEL{k}.extent_anom_sh_common MODEL{k}.volume_sh_common MODEL{k}.volume_clim_sh_common MODEL{k}.volume_anom_sh_common]=computeSeaIceDiagnostics(MODEL{k}.SIC_sh(:,ind1:ind2),MODEL{k}.SIT_sh(:,ind1:ind2),MODEL{k}.w_sh);

%compute clim over SIT satellite era
ind1=(SIT_clim_years(1)-year_start{k}-year_shift{k})*12+1;
ind2=(SIT_clim_years(2)+1-year_start{k}-year_shift{k})*12;
[MODEL{k}.SIT_clim_sh_common SIT_anom]=computeclim(MODEL{k}.SIT_sh(:,ind1:ind2));

%var=['siu_monthly_sh'];
%var=['siu_monthly_sh'];
var=['siu_trueE_sh'];
filename=[dataDir,'/',var,'.mat']
load(filename,'data')

MODEL{k}.siu_sh=100*data;

%var=['siv_monthly_sh'];
var=['siv_trueN_sh'];
filename=[dataDir,'/',var,'.mat']
load(filename,'data')

MODEL{k}.siv_sh=100*data;

%compute clim over OSISAF satellite era
ind1=(vel_clim_years(1)-year_start{k}-year_shift{k})*12+1;
ind2=(vel_clim_years(2)+1-year_start{k}-year_shift{k})*12;
[MODEL{k}.siu_clim_sh_common siu_anom]=computeclim(MODEL{k}.siu_sh(:,ind1:ind2));
[MODEL{k}.siv_clim_sh_common siv_anom]=computeclim(MODEL{k}.siv_sh(:,ind1:ind2));

%%%%load monthly obs data

%%%%NSIDC SIC on MODEL grid


%obsdir=['/work/mib/raw_data/NSIDC/NASATEAM_merged_v1.1_nrt_2018'];
%obsdir=['/work/mib/raw_data/NSIDC/NASATEAM_v1.1_',num2str(obs_year_end)];
%obsdir=['/work/mib/raw_data/NSIDC/CDR_',num2str(obs_year_end),'_monthly'];
obsdir=['/work/mib/raw_data/NSIDC/CDR_',num2str(obs_year_end),'_monthly_fromdaily'];
load([obsdir,'/sic_regrid_GFDL_',model{k},'.mat'],'data','x','y','ifXY','w')
OBS_SIC{k}=struct;
OBS_SIC{k}.SIC_nh=data/100;
OBS_SIC{k}.x_nh=x;
OBS_SIC{k}.y_nh=y;
OBS_SIC{k}.w_nh=w;
OBS_SIC{k}.ifXY_nh=ifXY;

OBS_SIC{k}.SIT_nh=zeros(size(OBS_SIC{k}.SIC_nh));

[OBS_SIC{k}.SIC_clim_nh OBS_SIC{k}.SIC_std_nh OBS_SIC{k}.SIT_clim_nh OBS_SIC{k}.SIT_std_nh OBS_SIC{k}.extent_nh OBS_SIC{k}.extent_clim_nh OBS_SIC{k}.extent_anom_nh OBS_SIC{k}.volume_nh OBS_SIC{k}.volume_clim_nh OBS_SIC{k}.volume_anom_nh]=computeSeaIceDiagnostics(OBS_SIC{k}.SIC_nh,OBS_SIC{k}.SIT_nh,OBS_SIC{k}.w_nh);

%compute clim over common period
if(isPI)
ind1=(1979-1979)*12+1;
ind2=(1990+1-1979)*12;
else
ind1=(SIC_clim_years(1)-1979)*12+1;
ind2=(SIC_clim_years(2)+1-1979)*12;
end

[OBS_SIC{k}.SIC_clim_nh_common OBS_SIC{k}.SIC_std_nh_common OBS_SIC{k}.SIT_clim_nh_common OBS_SIC{k}.SIT_std_nh_common OBS_SIC{k}.extent_nh_common OBS_SIC{k}.extent_clim_nh_common OBS_SIC{k}.extent_anom_nh_common OBS_SIC{k}.volume_nh_common OBS_SIC{k}.volume_clim_nh_common OBS_SIC{k}.volume_anom_nh_common]=computeSeaIceDiagnostics(OBS_SIC{k}.SIC_nh(:,ind1:ind2),OBS_SIC{k}.SIT_nh(:,ind1:ind2),OBS_SIC{k}.w_nh);

%compute regional SIE over common period
%load regional mask on SPEAR grid
modelDir=['/work/mib/raw_data/',model{k}];
load([modelDir,'/ice_grid_regions.mat'],'ifXY_reg')
[OBS_SIC{k}.regionalSIE OBS_SIC{k}.regionalSIEclim]=computeRegionalSIE(OBS_SIC{k}.SIC_nh(:,ind1:ind2),OBS_SIC{k}.w_nh,OBS_SIC{k}.ifXY_nh,ifXY_reg);

%obsdir=['/work/mib/raw_data/NSIDC/NASATEAM_merged_v1.1_nrt_2018_sh'];
%obsdir=['/work/mib/raw_data/NSIDC/NASATEAM_v1.1_',num2str(obs_year_end),'_sh'];
%obsdir=['/work/mib/raw_data/NSIDC/CDR_',num2str(obs_year_end),'_monthly_sh'];
obsdir=['/work/mib/raw_data/NSIDC/CDR_',num2str(obs_year_end),'_monthly_fromdaily_sh'];
load([obsdir,'/sic_regrid_GFDL_',model{k},'_45s.mat'],'data','x','y','ifXY','w')
OBS_SIC{k}.SIC_sh=data/100;
OBS_SIC{k}.x_sh=x;
OBS_SIC{k}.y_sh=y;
OBS_SIC{k}.w_sh=w;
OBS_SIC{k}.ifXY_sh=ifXY;

OBS_SIC{k}.SIT_sh=zeros(size(OBS_SIC{k}.SIC_sh));

[OBS_SIC{k}.SIC_clim_sh OBS_SIC{k}.SIC_std_sh OBS_SIC{k}.SIT_clim_sh OBS_SIC{k}.SIT_std_sh OBS_SIC{k}.extent_sh OBS_SIC{k}.extent_clim_sh OBS_SIC{k}.extent_anom_sh OBS_SIC{k}.volume_sh OBS_SIC{k}.volume_clim_sh OBS_SIC{k}.volume_anom_sh]=computeSeaIceDiagnostics(OBS_SIC{k}.SIC_sh,OBS_SIC{k}.SIT_sh,OBS_SIC{k}.w_sh);

%compute clim over common period
if(isPI)
ind1=(1979-1979)*12+1;
ind2=(1990+1-1979)*12;
else
ind1=(SIC_clim_years(1)-1979)*12+1;
ind2=(SIC_clim_years(2)+1-1979)*12;
end

[OBS_SIC{k}.SIC_clim_sh_common OBS_SIC{k}.SIC_std_sh_common OBS_SIC{k}.SIT_clim_sh_common OBS_SIC{k}.SIT_std_sh_common OBS_SIC{k}.extent_sh_common OBS_SIC{k}.extent_clim_sh_common OBS_SIC{k}.extent_anom_sh_common OBS_SIC{k}.volume_sh_common OBS_SIC{k}.volume_clim_sh_common OBS_SIC{k}.volume_anom_sh_common]=computeSeaIceDiagnostics(OBS_SIC{k}.SIC_sh(:,ind1:ind2),OBS_SIC{k}.SIT_sh(:,ind1:ind2),OBS_SIC{k}.w_sh);

end

%%%AWI SIT (spans 2011-obs_year_end, inclusive)

version='v2p6';

dataDir=['/work/mib/raw_data/CRYOSAT_AWI/',version];
filename=[dataDir,'/sit_',num2str(obs_year_end),'.mat'];

load(filename,'SIT','SIT_clim','SIC','volume','volume_clim','time','x','y','ifXY','w')

OBS_SIT.SIT_nh=SIT;
OBS_SIT.SIT_clim_nh=SIT_clim;
OBS_SIT.x_nh=x;
OBS_SIT.y_nh=y;
OBS_SIT.w_nh=w;
OBS_SIT.ifXY_nh=ifXY;
OBS_SIT.volume=volume;

%%%%%PIOMAS SIT (spans 1979-obs_year_end, inclusive)

dataDir=['/work/mib/raw_data/PIOMAS/v2.1_',num2str(obs_year_end)];
%load([dataDir,'/sit.mat'],'data','x','y','ifXY','w')
load([dataDir,'/heff.mat'],'data','x','y','ifXY','w')

nT=size(data,2);
vol_PIOMAS=zeros(nT,1);

for i=1:nT
vol_PIOMAS(i)=sum(data(:,i).*w);
end

for month=1:12
vol_PIOMAS_clim(month)=mean(vol_PIOMAS(month:12:end));
end

PIOMAS.volume_nh=vol_PIOMAS;
PIOMAS.volume_clim_nh=vol_PIOMAS_clim;

modelDir=['/work/mib/raw_data/PIOMAS'];
load([modelDir,'/ice_grid_regions.mat'],'ifXY_reg')

%compute clim over common period
if(isPI)
ind1=(1979-1979)*12+1;
ind2=(1990+1-1979)*12;
else
ind1=(SIC_clim_years(1)-1979)*12+1;
ind2=(SIC_clim_years(2)+1-1979)*12;
end

[PIOMAS.regionalSIV PIOMAS.regionalSIVclim]=computeRegionalSIV(ones(size(data(:,ind1:ind2))),data(:,ind1:ind2),w,ifXY,ifXY_reg);

[PIOMAS.SIC_clim_nh_common PIOMAS.SIC_std_nh_common PIOMAS.SIT_clim_nh_common PIOMAS.SIT_std_nh_common PIOMAS.extent_nh_common PIOMAS.extent_clim_nh_common PIOMAS.extent_anom_nh_common PIOMAS.volume_nh_common PIOMAS.volume_clim_nh_common PIOMAS.volume_anom_nh_common]=computeSeaIceDiagnostics(ones(size(data(:,ind1:ind2))),data(:,ind1:ind2),w);

%%OSI-SAF sivel
yrLim=[2010 obs_year_end];
dataDir=['/work/mib/raw_data/OSISAF/drift_lr_',num2str(yrLim(1)),'-',num2str(yrLim(2))];

filename=[dataDir,'/siu_trueE.mat'];
load(filename,'data','x','y','ifXY','w')

OBS_vel.siu_nh=100*data;
OBS_vel.x_nh=x;
OBS_vel.y_nh=y;
OBS_vel.w_nh=w;
OBS_vel.ifXY_nh=ifXY;

filename=[dataDir,'/siv_trueN.mat'];
load(filename,'data','time','x','y','ifXY','w')
OBS_vel.siv_nh=100*data;

yrLim=[2013 obs_year_end];
dataDir=['/work/mib/raw_data/OSISAF/drift_lr_',num2str(yrLim(1)),'-',num2str(yrLim(2))];

filename=[dataDir,'/siu_trueE_sh.mat'];
load(filename,'data','x','y','ifXY','w')

OBS_vel.siu_sh=100*data;
OBS_vel.x_sh=x;
OBS_vel.y_sh=y;
OBS_vel.w_sh=w;
OBS_vel.ifXY_sh=ifXY;

filename=[dataDir,'/siv_trueN_sh.mat'];
load(filename,'data','time','x','y','ifXY','w')
OBS_vel.siv_sh=100*data;

if(isPI)
ind1=(2010-2010)*12+1;
ind2=(obs_year_end+1-2010)*12;
else
%ind1=(vel_clim_years(1)-2010)*12+1;
%ind2=(vel_clim_years(2)+1-2010)*12;
%compute velocity clim over 2010-obs_year_end period, inclusive
ind1=(2010-2010)*12+1;
ind2=(obs_year_end+1-2010)*12;
end

[OBS_vel.siu_clim_nh siu_anom]=computeclim(OBS_vel.siu_nh(:,ind1:ind2));
[OBS_vel.siv_clim_nh siv_anom]=computeclim(OBS_vel.siv_nh(:,ind1:ind2));

if(isPI)
ind1=(2013-2013)*12+1;
ind2=(obs_year_end+1-2013)*12;
else
%ind1=(vel_clim_years(1)-2010)*12+1;
%ind2=(vel_clim_years(2)+1-2010)*12;
%compute velocity clim over 2010-obs_year_end period, inclusive
ind1=(2013-2013)*12+1;
ind2=(obs_year_end+1-2013)*12;
end

[OBS_vel.siu_clim_sh siu_anom]=computeclim(OBS_vel.siu_sh(:,ind1:ind2));
[OBS_vel.siv_clim_sh siv_anom]=computeclim(OBS_vel.siv_sh(:,ind1:ind2));

%%%%SAVE ALL DATA NEEDED TO MAKE FIGURES

%prior to saving, remove all time-resolved spatial data from MODEL and OBS structs

for k=1:N
MODEL{k}=rmfield(MODEL{k},'SIC_nh');
MODEL{k}=rmfield(MODEL{k},'SIT_nh');
MODEL{k}=rmfield(MODEL{k},'siu_nh');
MODEL{k}=rmfield(MODEL{k},'siv_nh');
MODEL{k}=rmfield(MODEL{k},'SIC_sh');
MODEL{k}=rmfield(MODEL{k},'SIT_sh');
MODEL{k}=rmfield(MODEL{k},'siu_sh');
MODEL{k}=rmfield(MODEL{k},'siv_sh');

OBS_SIC{k}=rmfield(OBS_SIC{k},'SIC_nh');
OBS_SIC{k}=rmfield(OBS_SIC{k},'SIT_nh');
OBS_SIC{k}=rmfield(OBS_SIC{k},'SIC_sh');
OBS_SIC{k}=rmfield(OBS_SIC{k},'SIT_sh');
end

OBS_vel=rmfield(OBS_vel,'siu_nh');
OBS_vel=rmfield(OBS_vel,'siv_nh');
OBS_vel=rmfield(OBS_vel,'siu_sh');
OBS_vel=rmfield(OBS_vel,'siv_sh');

work_dir=['/work/mib/processed_data/compare_sea_ice/',model{1},'/',compare_string];
filename=[work_dir,'/DiagnosticFigData.mat'];

save(filename,'MODEL','OBS_SIC','OBS_SIT','OBS_vel','PIOMAS','reg_list','-v7.3')

elseif(loadFigData==1)

work_dir=['/work/mib/processed_data/compare_sea_ice/',model{1},'/',compare_string];
filename=[work_dir,'/DiagnosticFigData.mat'];

load(filename,'MODEL','OBS_SIC','OBS_SIT','OBS_vel','PIOMAS','reg_list')

end

%%Arctic ice free date
%for month=1:12
%month
%for k=1:N
%k
%tmp=MODEL{k}.extent_nh(month:12:end);
%tmp_time=time_years{k};
%
%ind=find(tmp<1e12,1,'first')
%tmp(ind)
%tmp_time(ind)
%end
%end
%
%%Antarctic ice free date
%for month=1:12
%month
%for k=1:N
%k
%tmp=MODEL{k}.extent_sh(month:12:end);
%tmp_time=time_years{k};
%
%ind=find(tmp<1e12,1,'first')
%tmp(ind)
%tmp_time(ind)
%end
%end

%%compute trends
%[SIE_trend_obs]=computeMonthlyTrends(OBS_SIC{1}.extent_nh_common)
%
%for k=1:N
%k
%[SIE_trend]=computeMonthlyTrends(MODEL{k}.extent_nh_common)
%end
%
%
%%compute trends
%[SIE_trend_obs]=computeMonthlyTrends(OBS_SIC{1}.extent_sh_common)
%
%for k=1:N
%k
%[SIE_trend]=computeMonthlyTrends(MODEL{k}.extent_sh_common)
%end

%%%%%%MAKE FIGURES

colorlist={rgb('DodgerBlue'),rgb('Crimson'),rgb('LimeGreen'),'c','m',rgb('Gray'),rgb('Indigo')};

%%%%%%%FIGURE 1: NH/SH SIE

    nTileX=2;
    nTileY = 1;
        
    figWidth = 9.5; %6.5 inches = 39 picas, for a two column figure
    deltaX   = 0.6;
    deltaX2  = 0.2;
    deltaY   = 0.5;
    deltaY2  = 0.4;
    gapX     = 0.4;
    gapY     = 0.6;
    panel_scale=0.8;

ax=make_axes(nTileX,nTileY,figWidth,deltaX,deltaX2,deltaY,deltaY2,gapX,gapY,panel_scale);

axes(ax(1,1))
setFigProp(gca)

fact=1e-12;
L=1;

%albdiff=MODEL{3}.extent_clim_nh_common-MODEL{2}.extent_clim_nh_common;

hold on
plot(1:12,fact*OBS_SIC{1}.extent_clim_nh_common,'k','LineWidth',2)
for k=1:N
plot(1:12,fact*MODEL{k}.extent_clim_nh_common,'Color',colorlist{k},'LineWidth',L)
end
%plot(1:12,fact*(MODEL{1}.extent_clim_nh_common+2*albdiff),'Color',colorlist{1},'LineWidth',L,'LineStyle','--')
xlabel('Month')
ylabel('SIE (Million km$^2$)')
xlim([0.5 12.5])
grid on
if(isPI)
title(['Arctic SIE Climatology (',num2str(1979),'-',num2str(1990),')'])
else
title(['Arctic SIE Climatology (',num2str(SIC_clim_years(1)),'-',num2str(SIC_clim_years(2)),')'])
end
set(gca,'xTick',1:12,'xTickLabel',{'J','F','M','A','M','J','J','A','S','O','N','D'} )
legend(['NSIDC',short_name],'location','SouthWest')
yl=ylim;
text(-0.8,yl(1)+1.05*(yl(2)-yl(1)),'(a)','FontWeight','bold')

text(11.01,yl(1)+(0.3-0.05*0)*(yl(2)-yl(1)),'RMSE','Color','k')
for k=1:N
rms_tmp=sqrt(sum((fact*MODEL{k}.extent_clim_nh_common-fact*OBS_SIC{1}.extent_clim_nh_common).^2)/12);
text(11.1,yl(1)+(0.3-0.05*k)*(yl(2)-yl(1)),num2str(rms_tmp,'%.2f'),'Color',colorlist{k})
end

axes(ax(2,1))
setFigProp(gca)

%albdiff=MODEL{3}.extent_clim_sh_common-MODEL{2}.extent_clim_sh_common;

hold on
plot(1:12,fact*OBS_SIC{1}.extent_clim_sh_common,'k','LineWidth',2)
for k=1:N
plot(1:12,fact*MODEL{k}.extent_clim_sh_common,'Color',colorlist{k},'LineWidth',L)
end
%plot(1:12,fact*(MODEL{1}.extent_clim_sh_common+2*albdiff),'Color',colorlist{1},'LineWidth',L,'LineStyle','--')
xlabel('Month')
xlim([0.5 12.5])
grid on
if(isPI)
title(['Antarctic SIE Climatology (',num2str(1979),'-',num2str(1990),')'])
else
title(['Antarctic SIE Climatology (',num2str(SIC_clim_years(1)),'-',num2str(SIC_clim_years(2)),')'])
end
set(gca,'xTick',1:12,'xTickLabel',{'J','F','M','A','M','J','J','A','S','O','N','D'} )
yl=ylim;
text(-0.2,yl(1)+1.05*(yl(2)-yl(1)),'(b)','FontWeight','bold')
text(10.01,yl(1)+(0.3-0.05*0)*(yl(2)-yl(1)),'RMSE','Color','k')
for k=1:N
rms_tmp=sqrt(sum((fact*MODEL{k}.extent_clim_sh_common-fact*OBS_SIC{1}.extent_clim_sh_common).^2)/12);
text(10.1,yl(1)+(0.3-0.05*k)*(yl(2)-yl(1)),num2str(rms_tmp,'%.2f'),'Color',colorlist{k})
end

set(gcf,'Color','w')

filename=[fig_dir,'/SIE_clim'];
export_fig(gcf,filename,'-a1','-pdf');
export_fig(gcf,filename,'-a1','-png','-r400','-opengl');

%%%%%%FIGURE 2: NH/SH SIE timeseries

    nTileX=2;
    nTileY = 2;
        
    figWidth = 9.5; %6.5 inches = 39 picas, for a two column figure
    deltaX   = 0.7;
    deltaX2  = 0.2;
    deltaY   = 0.6;
    deltaY2  = 0.4;
    gapX     = 0.4;
    gapY     = 1.2;
    panel_scale=0.8;

ax=make_axes(nTileX,nTileY,figWidth,deltaX,deltaX2,deltaY,deltaY2,gapX,gapY,panel_scale);


fact=1e-12;
L=1;

axes(ax(1,1))
setFigProp(gca)
month=3;
hold on
plot(obs_years,fact*OBS_SIC{1}.extent_nh(month:12:end),'k','LineWidth',2)
for k=1:N
plot(time_years{k},fact*MODEL{k}.extent_nh(month:12:end),'Color',colorlist{k},'LineWidth',L)
end

grid on
l=set_legend(['NSIDC',short_name],'horizontal','SouthWest',0.14,-0.1,1,1);

%l=legend(['NSIDC',short_name],'location','SouthWest','orientation','horizontal');
%p=get(l,'position');
%set(l,'position',[p(1) p(2)-0.12 p(3) p(4)]);



xlim([plot_year_min plot_year_max])
xlabel('Time (years)')
ylabel('SIE (Million km$^2$)')
title('March Arctic SIE')
set_text_label('(a)','k',-0.1,1.07,12)

axes(ax(2,1))
setFigProp(gca)
month=9;
hold on
plot(obs_years,fact*OBS_SIC{1}.extent_nh(month:12:end),'k','LineWidth',2)
for k=1:N
plot(time_years{k},fact*MODEL{k}.extent_nh(month:12:end),'Color',colorlist{k},'LineWidth',L)
end
grid on
xlim([plot_year_min plot_year_max])
xlabel('Time (years)')
title('September Arctic SIE')
set_text_label('(b)','k',-0.05,1.07,12)


axes(ax(1,2))
setFigProp(gca)
month=3;
hold on
plot(obs_years,fact*OBS_SIC{1}.extent_sh(month:12:end),'k','LineWidth',2)
for k=1:N
plot(time_years{k},fact*MODEL{k}.extent_sh(month:12:end),'Color',colorlist{k},'LineWidth',L)
end
grid on
xlim([plot_year_min plot_year_max])
xlabel('Time (years)')
ylabel('SIE (Million km$^2$)')
title('March Antarctic SIE')
set_text_label('(c)','k',-0.1,1.07,12)

axes(ax(2,2))
setFigProp(gca)
month=9;
hold on
plot(obs_years,fact*OBS_SIC{1}.extent_sh(month:12:end),'k','LineWidth',2)
for k=1:N
plot(time_years{k},fact*MODEL{k}.extent_sh(month:12:end),'Color',colorlist{k},'LineWidth',L)
end
grid on
xlim([plot_year_min plot_year_max])
xlabel('Time (years)')
title('September Antarctic SIE')
set_text_label('(d)','k',-0.05,1.07,12)

set(gcf,'Color','w')
filename=[fig_dir,'/SIE_interannual'];
export_fig(gcf,filename,'-a1','-pdf');
export_fig(gcf,filename,'-a1','-png','-r400','-opengl');

%%%%%Figure 3: volume climatology

    nTileX=2;
    nTileY = 1;
        
    figWidth = 9.5; %6.5 inches = 39 picas, for a two column figure
    deltaX   = 0.6;
    deltaX2  = 0.2;
    deltaY   = 0.5;
    deltaY2  = 0.4;
    gapX     = 0.4;
    gapY     = 0.6;
    panel_scale=0.8;

ax=make_axes(nTileX,nTileY,figWidth,deltaX,deltaX2,deltaY,deltaY2,gapX,gapY,panel_scale);

axes(ax(1,1))
setFigProp(gca)

fact=1e-12;
L=1;

hold on
plot(1:12,fact*PIOMAS.volume_clim_nh_common,'k','LineWidth',2)
for k=1:N
plot(1:12,fact*MODEL{k}.volume_clim_nh_common,'Color',colorlist{k},'LineWidth',L)
end
%plot(1:12,fact*(MODEL{1}.volume_clim_nh_common+2*albdiff),'Color',colorlist{1},'LineWidth',L,'LineStyle','--')
xlabel('Month')
ylabel('SIV (10$^{12}$ m$^3$)')
xlim([0.5 12.5])
grid on
legend(['PIOMAS',short_name],'location','NorthEast')
if(isPI)
title(['Arctic SIV Climatology (',num2str(1979),'-',num2str(1990),')'])
else
title(['Arctic SIV Climatology (',num2str(SIC_clim_years(1)),'-',num2str(SIC_clim_years(2)),')'])
end
set(gca,'xTick',1:12,'xTickLabel',{'J','F','M','A','M','J','J','A','S','O','N','D'} )
yl=ylim;
set_text_label('(a)','k',-0.1,1.07,12)

axes(ax(2,1))
setFigProp(gca)

%albdiff=MODEL{3}.volume_clim_sh_common-MODEL{2}.volume_clim_sh_common;

hold on
%plot(1:12,fact*OBS_SIC.extent_clim_sh_common,'k','LineWidth',2)
for k=1:N
plot(1:12,fact*MODEL{k}.volume_clim_sh_common,'Color',colorlist{k},'LineWidth',L)
end
%plot(1:12,fact*(MODEL{1}.volume_clim_sh_common+2*albdiff),'Color',colorlist{1},'LineWidth',L,'LineStyle','--')
xlabel('Month')
xlim([0.5 12.5])
grid on
if(isPI)
title(['Antarctic SIV Climatology (',num2str(1979),'-',num2str(1990),')'])
else
title(['Antarctic SIV Climatology (',num2str(SIC_clim_years(1)),'-',num2str(SIC_clim_years(2)),')'])
end
set(gca,'xTick',1:12,'xTickLabel',{'J','F','M','A','M','J','J','A','S','O','N','D'} )
yl=ylim;
set_text_label('(b)','k',-0.05,1.07,12)

set(gcf,'Color','w')
filename=[fig_dir,'/SIV_clim'];
export_fig(gcf,filename,'-a1','-pdf');
export_fig(gcf,filename,'-a1','-png','-r400','-opengl');

%%%%%%FIGURE 4: NH/SH SIV timeseries

    nTileX=2;
    nTileY = 2;
        
    figWidth = 9.5; %6.5 inches = 39 picas, for a two column figure
    deltaX   = 0.7;
    deltaX2  = 0.2;
    deltaY   = 0.5;
    deltaY2  = 0.6;
    gapX     = 0.4;
    gapY     = 1.2;
    panel_scale=0.8;

ax=make_axes(nTileX,nTileY,figWidth,deltaX,deltaX2,deltaY,deltaY2,gapX,gapY,panel_scale);


fact=1e-12;
L=1;

axes(ax(1,1))
setFigProp(gca)
month=3;
hold on
plot(PIOMAS_years,fact*PIOMAS.volume_nh(month:12:end),'k','LineWidth',2)
for k=1:N
plot(time_years{k},fact*MODEL{k}.volume_nh(month:12:end),'Color',colorlist{k},'LineWidth',L)
end

grid on
l=legend(['PIOMAS',short_name],'location','SouthWest','orientation','horizontal');
p=get(l,'position');
set(l,'position',[p(1)+0.14 p(2)-0.1 p(3) p(4)]);

xlim([plot_year_min plot_year_max])
xlabel('Time (years)')
ylabel('SIV (10$^{12}$ m$^3$)')
title('March Arctic SIV')
set_text_label('(a)','k',-0.05,1.07,12)

axes(ax(2,1))
setFigProp(gca)
month=9;
hold on
plot(PIOMAS_years,fact*PIOMAS.volume_nh(month:12:end),'k','LineWidth',2)
for k=1:N
plot(time_years{k},fact*MODEL{k}.volume_nh(month:12:end),'Color',colorlist{k},'LineWidth',L)
end
grid on
xlim([plot_year_min plot_year_max])
xlabel('Time (years)')
title('September Arctic SIV')
set_text_label('(b)','k',-0.05,1.07,12)


axes(ax(1,2))
setFigProp(gca)
month=3;
hold on
%plot(obs_years,fact*OBS_SIC.extent_sh(month:12:end),'k','LineWidth',2)
for k=1:N
plot(time_years{k},fact*MODEL{k}.volume_sh(month:12:end),'Color',colorlist{k},'LineWidth',L)
end
grid on
xlim([plot_year_min plot_year_max])
xlabel('Time (years)')
ylabel('SIV (10$^{12}$ m$^3$)')
title('March Antarctic SIV')
set_text_label('(c)','k',-0.05,1.07,12)

axes(ax(2,2))
setFigProp(gca)
month=9;
hold on
%plot(obs_years,fact*OBS_SIC.extent_sh(month:12:end),'k','LineWidth',2)
for k=1:N
plot(time_years{k},fact*MODEL{k}.volume_sh(month:12:end),'Color',colorlist{k},'LineWidth',L)
end
grid on
xlim([plot_year_min plot_year_max])
xlabel('Time (years)')
title('September Antarctic SIV')
set_text_label('(d)','k',-0.05,1.07,12)

set(gcf,'Color','w')
filename=[fig_dir,'/SIV_interannual'];
export_fig(gcf,filename,'-a1','-pdf');
export_fig(gcf,filename,'-a1','-png','-r400','-opengl');

%%%FIGURE 5: SIC biases

    nTileX=N;
    nTileY = 4;
        
    figWidth = 6.5; %6.5 inches = 39 picas, for a two column figure
    deltaX   = 0.4;
    deltaX2  = 0.02;
    deltaY   = 0.7;
    deltaY2  = 0.3;
    gapX     = 0.01;
    gapY     = 0.01;
    panel_scale=1;

ax=make_axes(nTileX,nTileY,figWidth,deltaX,deltaX2,deltaY,deltaY2,gapX,gapY,panel_scale);

month=3;
std_threshold=0.05;

for k=1:N

axes(ax(k,1))
setFigProp(gca)
tmp_obs=OBS_SIC{k}.SIC_clim_nh_common(:,month);
tmp=MODEL{k}.SIC_clim_nh_common(:,month);

make_polar_plot_pcolor_modelmask(MODEL{k}.xC_nh,MODEL{k}.yC_nh,MODEL{k}.ifXY_nh,tmp-tmp_obs,45,'',0,-1,1);
make_polar_contour_plot_modelmask(MODEL{k}.x_nh,MODEL{k}.y_nh,MODEL{k}.ifXY_nh,tmp_obs,45,'','k',1)
%text(-0.75,0.9,'A','FontWeight','bold')
if(k==1)
text(-1,-0.5,'March SIC Bias','FontWeight','bold','rotation',90)
end
%text(-0.7,0.94,short_name{k},'FontWeight','bold','FontSize',12)
title(short_name{k})

obs_var_mask=(OBS_SIC{k}.SIC_std_nh_common(:,month)>=std_threshold);
model_var_mask=(MODEL{k}.SIC_std_nh_common(:,month)>=std_threshold);
var_mask=(obs_var_mask+model_var_mask)>0;
tmp_sic_error=sqrt(nansum((tmp-tmp_obs).^2.*var_mask.*MODEL{k}.w_nh)/nansum(var_mask.*MODEL{k}.w_nh));
set_text_label(num2str(tmp_sic_error,'%.2f'),'k',0.7,0.5,10)

end

month=9;

for k=1:N

axes(ax(k,2))
setFigProp(gca)

tmp_obs=OBS_SIC{k}.SIC_clim_nh_common(:,month);
tmp=MODEL{k}.SIC_clim_nh_common(:,month);

make_polar_plot_pcolor_modelmask(MODEL{k}.xC_nh,MODEL{k}.yC_nh,MODEL{k}.ifXY_nh,tmp-tmp_obs,30,'',0,-1,1);
make_polar_contour_plot_modelmask(MODEL{k}.x_nh,MODEL{k}.y_nh,MODEL{k}.ifXY_nh,tmp_obs,30,'','k',1)

if(k==1)
text(-0.65,-0.3,'Sept SIC Bias','FontWeight','bold','rotation',90)
end

obs_var_mask=(OBS_SIC{k}.SIC_std_nh_common(:,month)>=std_threshold);
model_var_mask=(MODEL{k}.SIC_std_nh_common(:,month)>=std_threshold);
var_mask=(obs_var_mask+model_var_mask)>0;
tmp_sic_error=sqrt(nansum((tmp-tmp_obs).^2.*var_mask.*MODEL{k}.w_nh)/nansum(var_mask.*MODEL{k}.w_nh));
set_text_label(num2str(tmp_sic_error,'%.2f'),'k',0.8,0.5,10)

end


%SH SIC

month=9;

for k=1:N

axes(ax(k,3))
setFigProp(gca)

tmp_obs=OBS_SIC{k}.SIC_clim_sh_common(:,month);
tmp=MODEL{k}.SIC_clim_sh_common(:,month);

make_polar_plot_pcolor_modelmask_sh(MODEL{k}.x_sh,MODEL{k}.y_sh,MODEL{k}.ifXY_sh,tmp-tmp_obs,40,'',0,-1,1);
make_polar_contour_plot_modelmask_sh(MODEL{k}.x_sh,MODEL{k}.y_sh,MODEL{k}.ifXY_sh,tmp_obs,40,'','k',1)

if(k==1)
text(-0.85,-0.3,'Sept SIC Bias','FontWeight','bold','rotation',90)
end

obs_var_mask=(OBS_SIC{k}.SIC_std_sh_common(:,month)>=std_threshold);
model_var_mask=(MODEL{k}.SIC_std_sh_common(:,month)>=std_threshold);
var_mask=(obs_var_mask+model_var_mask)>0;
tmp_sic_error=sqrt(nansum((tmp-tmp_obs).^2.*var_mask.*MODEL{k}.w_sh)/nansum(var_mask.*MODEL{k}.w_sh));
set_text_label(num2str(tmp_sic_error,'%.2f'),'k',0.52,0.5,10)

end
%text(-0.7,0.8,'G','FontWeight','bold')

month=3;

for k=1:N

axes(ax(k,4))
setFigProp(gca)

tmp_obs=OBS_SIC{k}.SIC_clim_sh_common(:,month);
tmp=MODEL{k}.SIC_clim_sh_common(:,month);

make_polar_plot_pcolor_modelmask_sh(MODEL{k}.x_sh,MODEL{k}.y_sh,MODEL{k}.ifXY_sh,tmp-tmp_obs,30,'',0,-1,1);
make_polar_contour_plot_modelmask_sh(MODEL{k}.x_sh,MODEL{k}.y_sh,MODEL{k}.ifXY_sh,tmp_obs,30,'','k',1)

obs_var_mask=(OBS_SIC{k}.SIC_std_sh_common(:,month)>=std_threshold);
model_var_mask=(MODEL{k}.SIC_std_sh_common(:,month)>=std_threshold);
var_mask=(obs_var_mask+model_var_mask)>0;
tmp_sic_error=sqrt(nansum((tmp-tmp_obs).^2.*var_mask.*MODEL{k}.w_sh)/nansum(var_mask.*MODEL{k}.w_sh));
set_text_label(num2str(tmp_sic_error,'%.2f'),'k',0.55,0.5,10)

if(k==1)
text(-0.6,-0.3,'March SIC Bias','FontWeight','bold','rotation',90)


                 pos=get(gca,'position');

                 c=colorbar('horizontal');

                 set(gca,'position',pos)

                 x_pos=get(c, 'Position');
                 y_pos=[x_pos(1)+0.05 x_pos(2)+0.05 2.5*x_pos(3) 0.7*x_pos(4)];
                 set(c,'Position',y_pos)

end

end

set(gcf,'Color','w')
filename=[fig_dir,'/SIC_bias.png'];
%export_fig(gcf,filename,'-a1','-png','-r400','-painters');
export_fig(gcf,filename,'-a1','-png','-r400','-opengl');

%%%FIGURE 6: SIT clim

    nTileX=N+1;
    nTileY = 4;
        
    figWidth = 6.5; %6.5 inches = 39 picas, for a two column figure
    deltaX   = 0.4;
    deltaX2  = 0.02;
    deltaY   = 0.5;
    deltaY2  = 0.3;
    gapX     = 0.02;
    gapY     = 0.1;
    panel_scale=1;

ax=make_axes(nTileX,nTileY,figWidth,deltaX,deltaX2,deltaY,deltaY2,gapX,gapY,panel_scale);


%load('/home/Mitchell.Bushuk/scripts/utils/cmap_red_seq.mat')
load('/home/Mitchell.Bushuk/scripts/utils/cmap_red_seq2.mat')
%cmap=cmocean('amp',20);
colormap(cmap)

month=3;

SIT_max=3.5;
SIT_max_summer=2.5;
SIT_max_sh=2;
%SIT_max=20;
%SIT_max_summer=20;
%SIT_max_sh=20;

for k=1:N

axes(ax(k,1))
setFigProp(gca)
tmp=MODEL{k}.SIT_clim_nh_common(:,month);

make_polar_plot_pcolor_modelmask(MODEL{k}.xC_nh,MODEL{k}.yC_nh,MODEL{k}.ifXY_nh,tmp,45,'',0,0,SIT_max);
%text(-0.75,0.9,'A','FontWeight','bold')
if(k==1)
text(-1,-0.3,'March SIT','FontWeight','bold','rotation',90)

end
text(-0.7,0.94,short_name{k},'FontWeight','bold','FontSize',12)

end

axes(ax(N+1,1))
setFigProp(gca)
tmp_obs=OBS_SIT.SIT_clim_nh(:,5);%March SIT
make_polar_plot_pcolor_modelmask(OBS_SIT.x_nh,OBS_SIT.y_nh,OBS_SIT.ifXY_nh,tmp_obs,45,[''],0,0,SIT_max)
hold on
%text(-0.75,0.9,'D','FontWeight','bold')
text(-0.4,0.94,'CryoSat-2','FontWeight','bold','FontSize',12)

       pos=get(gca,'position');

       c=colorbar('vertical');

       set(gca,'position',pos)

       x_pos=get(c, 'Position');
       y_pos=[x_pos(1) x_pos(2) 1*x_pos(3) 1*x_pos(4)];
       set(c,'Position',y_pos)

month=9;

for k=1:N

axes(ax(k,2))
setFigProp(gca)

tmp=MODEL{k}.SIT_clim_nh_common(:,month);

make_polar_plot_pcolor_modelmask(MODEL{k}.xC_nh,MODEL{k}.yC_nh,MODEL{k}.ifXY_nh,tmp,30,'',0,0,SIT_max_summer);

if(k==1)
text(-0.6,-0.3,'Sept SIT','FontWeight','bold','rotation',90)
end

end

axes(ax(N,2))
setFigProp(gca)


       pos=get(gca,'position');

       c=colorbar('vertical');

       set(gca,'position',pos)

       x_pos=get(c, 'Position');
       y_pos=[x_pos(1)+0.02 x_pos(2) 1*x_pos(3) 1*x_pos(4)];
       set(c,'Position',y_pos)

%SH SIC

month=9;

for k=1:N

axes(ax(k,3))
setFigProp(gca)

tmp=MODEL{k}.SIT_clim_sh_common(:,month);

make_polar_plot_pcolor_modelmask_sh(MODEL{k}.x_sh,MODEL{k}.y_sh,MODEL{k}.ifXY_sh,tmp,40,'',0,0,SIT_max_sh);

if(k==1)
text(-0.85,-0.3,'Sept SIT','FontWeight','bold','rotation',90)
end

end
%text(-0.7,0.8,'G','FontWeight','bold')

axes(ax(N,3))
setFigProp(gca)


       pos=get(gca,'position');

       c=colorbar('vertical');

       set(gca,'position',pos)

       x_pos=get(c, 'Position');
       y_pos=[x_pos(1)+0.05 x_pos(2) 1*x_pos(3) 1*x_pos(4)];
       set(c,'Position',y_pos)

month=3;

for k=1:N

axes(ax(k,4))
setFigProp(gca)

tmp=MODEL{k}.SIT_clim_sh_common(:,month);

make_polar_plot_pcolor_modelmask_sh(MODEL{k}.x_sh,MODEL{k}.y_sh,MODEL{k}.ifXY_sh,tmp,30,'',0,0,SIT_max_sh);

if(k==1)
text(-0.6,-0.3,'March SIT','FontWeight','bold','rotation',90)
end

end

axes(ax(N,4))
setFigProp(gca)


       pos=get(gca,'position');

       c=colorbar('vertical');

       set(gca,'position',pos)

       x_pos=get(c, 'Position');
       y_pos=[x_pos(1)+0.05 x_pos(2) 1*x_pos(3) 1*x_pos(4)];
       set(c,'Position',y_pos)

colormap(cmap)

set(gcf,'Color','w')
filename=[fig_dir,'/SIT.png'];
%export_fig(gcf,filename,'-a1','-png','-r400','-painters');
export_fig(gcf,filename,'-a1','-png','-r400','-opengl');

%%%FIGURE 7: SIT.*SIC clim

    nTileX=N+1;
    nTileY = 4;
        
    figWidth = 6.5; %6.5 inches = 39 picas, for a two column figure
    deltaX   = 0.4;
    deltaX2  = 0.02;
    deltaY   = 0.5;
    deltaY2  = 0.3;
    gapX     = 0.02;
    gapY     = 0.1;
    panel_scale=1;

ax=make_axes(nTileX,nTileY,figWidth,deltaX,deltaX2,deltaY,deltaY2,gapX,gapY,panel_scale);


%load('/home/Mitchell.Bushuk/scripts/utils/cmap_red_seq.mat')
load('/home/Mitchell.Bushuk/scripts/utils/cmap_red_seq2.mat')
colormap(cmap)

month=3;

for k=1:N

axes(ax(k,1))
setFigProp(gca)
tmp=MODEL{k}.SIT_clim_nh_common(:,month).*MODEL{k}.SIC_clim_nh_common(:,month);

make_polar_plot_pcolor_modelmask(MODEL{k}.xC_nh,MODEL{k}.yC_nh,MODEL{k}.ifXY_nh,tmp,45,'',0,0,SIT_max);
%text(-0.75,0.9,'A','FontWeight','bold')
if(k==1)
text(-1,-0.3,'March SIT mean','FontWeight','bold','rotation',90)

end
text(-0.7,0.94,short_name{k},'FontWeight','bold','FontSize',12)

end

axes(ax(N+1,1))
setFigProp(gca)
tmp_obs=OBS_SIT.SIT_clim_nh(:,5);%March SIT
make_polar_plot_pcolor_modelmask(OBS_SIT.x_nh,OBS_SIT.y_nh,OBS_SIT.ifXY_nh,tmp_obs,45,[''],0,0,SIT_max)
hold on
%text(-0.75,0.9,'D','FontWeight','bold')
text(-0.4,0.94,'CryoSat-2','FontWeight','bold','FontSize',12)

       pos=get(gca,'position');

       c=colorbar('vertical');

       set(gca,'position',pos)

       x_pos=get(c, 'Position');
       y_pos=[x_pos(1) x_pos(2) 1*x_pos(3) 1*x_pos(4)];
       set(c,'Position',y_pos)

month=9;

for k=1:N

axes(ax(k,2))
setFigProp(gca)

tmp=MODEL{k}.SIT_clim_nh_common(:,month).*MODEL{k}.SIC_clim_nh_common(:,month);

make_polar_plot_pcolor_modelmask(MODEL{k}.xC_nh,MODEL{k}.yC_nh,MODEL{k}.ifXY_nh,tmp,30,'',0,0,SIT_max_summer);

if(k==1)
text(-0.6,-0.3,'Sept SIT mean','FontWeight','bold','rotation',90)
end

end

axes(ax(N,2))
setFigProp(gca)


       pos=get(gca,'position');

       c=colorbar('vertical');

       set(gca,'position',pos)

       x_pos=get(c, 'Position');
       y_pos=[x_pos(1)+0.02 x_pos(2) 1*x_pos(3) 1*x_pos(4)];
       set(c,'Position',y_pos)

%SH SIC

month=9;

for k=1:N

axes(ax(k,3))
setFigProp(gca)

tmp=MODEL{k}.SIT_clim_sh_common(:,month).*MODEL{k}.SIC_clim_sh_common(:,month);

make_polar_plot_pcolor_modelmask_sh(MODEL{k}.x_sh,MODEL{k}.y_sh,MODEL{k}.ifXY_sh,tmp,40,'',0,0,SIT_max_sh);

if(k==1)
text(-0.85,-0.3,'Sept SIT mean','FontWeight','bold','rotation',90)
end

end
%text(-0.7,0.8,'G','FontWeight','bold')

axes(ax(N,3))
setFigProp(gca)


       pos=get(gca,'position');

       c=colorbar('vertical');

       set(gca,'position',pos)

       x_pos=get(c, 'Position');
       y_pos=[x_pos(1)+0.05 x_pos(2) 1*x_pos(3) 1*x_pos(4)];
       set(c,'Position',y_pos)

month=3;

for k=1:N

axes(ax(k,4))
setFigProp(gca)

tmp=MODEL{k}.SIT_clim_sh_common(:,month).*MODEL{k}.SIC_clim_sh_common(:,month);

make_polar_plot_pcolor_modelmask_sh(MODEL{k}.x_sh,MODEL{k}.y_sh,MODEL{k}.ifXY_sh,tmp,30,'',0,0,SIT_max_sh);

if(k==1)
text(-0.6,-0.3,'March SIT mean','FontWeight','bold','rotation',90)
end

end

axes(ax(N,4))
setFigProp(gca)


       pos=get(gca,'position');

       c=colorbar('vertical');

       set(gca,'position',pos)

       x_pos=get(c, 'Position');
       y_pos=[x_pos(1)+0.05 x_pos(2) 1*x_pos(3) 1*x_pos(4)];
       set(c,'Position',y_pos)

colormap(cmap)

set(gcf,'Color','w')
filename=[fig_dir,'/SIT_mean.png'];
%export_fig(gcf,filename,'-a1','-png','-r400','-painters');
export_fig(gcf,filename,'-a1','-png','-r400','-opengl');

%%%%%Figure 8: Ice thickness and motion
%this figure plots obs and the first three models available

    nTileX=2;
    nTileY = 2;
        
    figWidth = 9; %6.5 inches = 39 picas, for a two column figure
    deltaX   = 0.05;
    deltaX2  = 0.6;
    deltaY   = 0.5;
    deltaY2  = 0.5;
    gapX     = 0.02;
    gapY     = 0.1;
    panel_scale=1;

ax=make_axes(nTileX,nTileY,figWidth,deltaX,deltaX2,deltaY,deltaY2,gapX,gapY,panel_scale);

%vecplot parameters
K=3;
head_length=2.5;
shaft_width=0.8;
radius=34;
scale=50;

SIT_max=3.5;


%load('/home/Mitchell.Bushuk/scripts/utils/cmap_red_seq.mat')
load('/home/Mitchell.Bushuk/scripts/utils/cmap_red_seq2.mat')
colormap(cmap)

month=3;

for k=1:min(N,3)

col=mod(k,2);
if(col==0)
col=2;
end

row=1+(k-col)/2;

axes(ax(col,row))
setFigProp(gca)
tmp=MODEL{k}.SIT_clim_nh_common(:,month);

make_polar_plot_pcolor_modelmask(MODEL{k}.xC_nh,MODEL{k}.yC_nh,MODEL{k}.ifXY_nh,tmp,radius,'',0,0,SIT_max);

tmp1=mean(MODEL{k}.siu_clim_nh_common(:,[1 2 3 12]),2);
tmp2=mean(MODEL{k}.siv_clim_nh_common(:,[1 2 3 12]),2);
SIC=mean(MODEL{k}.SIC_clim_nh_common(:,[1 2 3 12]),2);
SIC_mask=SIC>=0.9;
tmp1=tmp1.*SIC_mask;
tmp2=tmp2.*SIC_mask;

make_polar_plot_vec_modelmask(MODEL{k}.x_nh,MODEL{k}.y_nh,MODEL{k}.ifXY_nh,tmp1,tmp2,radius,'',K,head_length,shaft_width,scale)

set_text_label(label_str{k},'k',0.05,0.9,12)
if(k==1)
text(-1,-0.3,'March SIT','FontWeight','bold','rotation',90)

end
text(-0.4,0.35,short_name{k},'FontWeight','bold','FontSize',12)

%%%%compute mean SIT over mass budget region from Bushuk and Polvani 2023

%%load SIMIP sea ice mask

%use for ODS data
load('/work/mib/raw_data/LE_data/CESM_ODS/ice_grid_regions.mat','x','y','ifXY','w','ifXY_reg')
inds=[1 4 5 6 7 10 11];

x_SIMIP=x;
y_SIMIP=y;
w_SIMIP=w;
ifXY_SIMIP=ifXY;
ifXY_reg_SIMIP=squeeze(sum(ifXY_reg(:,:,inds),3));

%%%%regrid SIMIP to GFDL grid

GFDL_gridpoints=ones(size(MODEL{k}.x_nh));
GFDL_gridpoints=GFDL_gridpoints>0;

x_SIMIP=mod(x_SIMIP+360,360);
SIMIP_gridpoints=ones(size(x_SIMIP));
SIMIP_gridpoints=SIMIP_gridpoints>0;
SIMIP_land=~ifXY_SIMIP;

nD=nnz(MODEL{k}.ifXY_nh);
regrid_data=zeros(size(MODEL{k}.x_nh,1),size(MODEL{k}.x_nh,2));

A=ifXY_reg_SIMIP;

dataTmp=A(:);
x=x_SIMIP(:);
y=y_SIMIP(:);

xq=MODEL{k}.x_nh(GFDL_gridpoints);
yq=MODEL{k}.y_nh(GFDL_gridpoints);

disp('Computing interpolant...')
F = scatteredInterpolant(x,y,dataTmp,'nearest','nearest');
vq=F(xq,yq);

disp('Done!')

A=zeros(size(MODEL{k}.x_nh));
A(GFDL_gridpoints)=vq;

regrid_data=A;

MODEL{k}.ifXY_reg_nh=regrid_data;

[mean_SIT regionalSumclim regionalSumsig]=computeRegionalMean(tmp,MODEL{k}.w_nh,MODEL{k}.ifXY_nh,MODEL{k}.ifXY_reg_nh);

%%%%%%%%%%%%%

text(0.125,0.35,['Mean SIT: ',num2str(mean_SIT,'%.1f'),' m'],'FontWeight','bold','FontSize',12)
end

tmp1=mean(OBS_vel.siu_clim_nh(:,[1 2 3 12]),2);
tmp2=mean(OBS_vel.siv_clim_nh(:,[1 2 3 12]),2);

axes(ax(2,2))
setFigProp(gca)
tmp_obs=OBS_SIT.SIT_clim_nh(:,5);%March SIT
make_polar_plot_pcolor_modelmask(OBS_SIT.x_nh,OBS_SIT.y_nh,OBS_SIT.ifXY_nh,tmp_obs,radius,[''],0,0,SIT_max)
hold on
make_polar_plot_vec_modelmask(OBS_vel.x_nh,OBS_vel.y_nh,OBS_vel.ifXY_nh,tmp1,tmp2,radius,'',K,head_length,shaft_width,scale)
set_text_label('(d)','k',0.05,0.9,12)
%text(-0.4,0.35,'CryoSat-2','FontWeight','bold','FontSize',12)
text(-0.41,0.35,'Observations','FontWeight','bold','FontSize',12)

%%%%compute mean SIT over mass budget region from Bushuk and Polvani 2023

%%load SIMIP sea ice mask

%use for ODS data
load('/work/mib/raw_data/LE_data/CESM_ODS/ice_grid_regions.mat','x','y','ifXY','w','ifXY_reg')
inds=[1 4 5 6 7 10 11];

x_SIMIP=x;
y_SIMIP=y;
w_SIMIP=w;
ifXY_SIMIP=ifXY;
ifXY_reg_SIMIP=squeeze(sum(ifXY_reg(:,:,inds),3));

%%%%regrid SIMIP to GFDL grid

GFDL_gridpoints=ones(size(OBS_SIT.x_nh));
GFDL_gridpoints=GFDL_gridpoints>0;

x_SIMIP=mod(x_SIMIP+360,360);
SIMIP_gridpoints=ones(size(x_SIMIP));
SIMIP_gridpoints=SIMIP_gridpoints>0;
SIMIP_land=~ifXY_SIMIP;

nD=nnz(OBS_SIT.ifXY_nh);
regrid_data=zeros(size(OBS_SIT.x_nh,1),size(OBS_SIT.x_nh,2));

A=ifXY_reg_SIMIP;

dataTmp=A(:);
x=x_SIMIP(:);
y=y_SIMIP(:);

xq=OBS_SIT.x_nh(GFDL_gridpoints);
yq=OBS_SIT.y_nh(GFDL_gridpoints);

disp('Computing interpolant...')
F = scatteredInterpolant(x,y,dataTmp,'nearest','nearest');
vq=F(xq,yq);

disp('Done!')

A=zeros(size(OBS_SIT.x_nh));
A(GFDL_gridpoints)=vq;

regrid_data=A;

OBS_SIT.ifXY_reg_nh=regrid_data;

[mean_SIT regionalSumclim regionalSumsig]=computeRegionalMean(tmp_obs,OBS_SIT.w_nh,OBS_SIT.ifXY_nh,OBS_SIT.ifXY_reg_nh);

%%%%%%%%%%%%%

text(0.125,0.35,['Mean SIT: ',num2str(mean_SIT,'%.1f'),' m'],'FontWeight','bold','FontSize',12)

       set_colorbar(0.07,0.05,1,2)

colormap(cmap)

set(gcf,'Color','w')
filename=[fig_dir,'/SIT_vel.png'];
%export_fig(gcf,filename,'-a1','-png','-r400','-painters');
export_fig(gcf,filename,'-a1','-png','-r400','-opengl');

%%%%%%Figure 9: Ice thickness and motion
%this figure plots obs and the first model and models 4,5 (if available)

    nTileX=2;
    nTileY = 2;
        
    figWidth = 9; %6.5 inches = 39 picas, for a two column figure
    deltaX   = 0.4;
    deltaX2  = 0.02;
    deltaY   = 0.5;
    deltaY2  = 0.5;
    gapX     = 0.02;
    gapY     = 0.1;
    panel_scale=1;

ax=make_axes(nTileX,nTileY,figWidth,deltaX,deltaX2,deltaY,deltaY2,gapX,gapY,panel_scale);

%vecplot parameters
K=3;
head_length=2.5;
shaft_width=0.8;
radius=34;
scale=50;

%load('/home/Mitchell.Bushuk/scripts/utils/cmap_red_seq.mat')
load('/home/Mitchell.Bushuk/scripts/utils/cmap_red_seq2.mat')
colormap(cmap)

month=3;

if(N>=5)
tmp_list=[1 4 5];
elseif(N==4)
tmp_list=[1 3 4];
elseif(N<=3)
tmp_list=[1:min(N,3)];
end

for k=1:min(N,3)

col=mod(k,2);
if(col==0)
col=2;
end

row=1+(k-col)/2;

axes(ax(col,row))
setFigProp(gca)
tmp=MODEL{tmp_list(k)}.SIT_clim_nh_common(:,month);

make_polar_plot_pcolor_modelmask(MODEL{tmp_list(k)}.xC_nh,MODEL{tmp_list(k)}.yC_nh,MODEL{tmp_list(k)}.ifXY_nh,tmp,radius,'',0,0,SIT_max);

tmp1=mean(MODEL{tmp_list(k)}.siu_clim_nh_common(:,[1 2 3 12]),2);
tmp2=mean(MODEL{tmp_list(k)}.siv_clim_nh_common(:,[1 2 3 12]),2);
SIC=mean(MODEL{tmp_list(k)}.SIC_clim_nh_common(:,[1 2 3 12]),2);
SIC_mask=SIC>=0.8;
tmp1=tmp1.*SIC_mask;
tmp2=tmp2.*SIC_mask;

make_polar_plot_vec_modelmask(MODEL{tmp_list(k)}.x_nh,MODEL{tmp_list(k)}.y_nh,MODEL{tmp_list(k)}.ifXY_nh,tmp1,tmp2,radius,'',K,head_length,shaft_width,scale)

%text(-0.75,0.9,'A','FontWeight','bold')
if(k==1)
text(-0.8,-0.3,'March SIT','FontWeight','bold','rotation',90)

end
text(-0.4,0.35,short_name{tmp_list(k)},'FontWeight','bold','FontSize',12)

end

tmp1=mean(OBS_vel.siu_clim_nh(:,[1 2 3 12]),2);
tmp2=mean(OBS_vel.siv_clim_nh(:,[1 2 3 12]),2);

axes(ax(2,2))
setFigProp(gca)
tmp_obs=OBS_SIT.SIT_clim_nh(:,5);%March SIT
make_polar_plot_pcolor_modelmask(OBS_SIT.x_nh,OBS_SIT.y_nh,OBS_SIT.ifXY_nh,tmp_obs,radius,[''],0,0,SIT_max)
hold on
make_polar_plot_vec_modelmask(OBS_vel.x_nh,OBS_vel.y_nh,OBS_vel.ifXY_nh,tmp1,tmp2,radius,'',K,head_length,shaft_width,scale)
%text(-0.75,0.9,'D','FontWeight','bold')
text(-0.4,0.35,'CryoSat-2','FontWeight','bold','FontSize',12)

       pos=get(gca,'position');

       c=colorbar('vertical');

       set(gca,'position',pos)

       x_pos=get(c, 'Position');
       y_pos=[x_pos(1) x_pos(2) 1*x_pos(3) 1*x_pos(4)];
       set(c,'Position',y_pos)

colormap(cmap)

set(gcf,'Color','w')
filename=[fig_dir,'/SIT_vel2.png'];
%export_fig(gcf,filename,'-a1','-png','-r400','-painters');
export_fig(gcf,filename,'-a1','-png','-r400','-opengl');

%%%%%%Figure 10: Ice concentration and motion
%this figure plots obs and the first three models available

    nTileX=2;
    nTileY = 2;
        
    figWidth = 9; %6.5 inches = 39 picas, for a two column figure
    deltaX   = 0.05;
    deltaX2  = 0.7;
    deltaY   = 0.5;
    deltaY2  = 0.5;
    gapX     = 0.02;
    gapY     = 0.1;
    panel_scale=1;

ax=make_axes(nTileX,nTileY,figWidth,deltaX,deltaX2,deltaY,deltaY2,gapX,gapY,panel_scale);

%vecplot parameters
%K=3;
%head_length=2.5;
%shaft_width=0.8;
%radius=34;
%scale=50;

%K=2;
%head_length=1.5;
%shaft_width=0.3;
%radius=40;
%scale=110;

K=3;
head_length=2.0;
shaft_width=0.5;
radius=40;
scale=110;


%load('/home/Mitchell.Bushuk/scripts/utils/cmap_red_seq.mat')
%load('/home/Mitchell.Bushuk/scripts/utils/cmap_red_seq2.mat')
cmap = cmocean('ice');
colormap(cmap)

month=9;

for k=1:min(N,3)

col=mod(k,2);
if(col==0)
col=2;
end

row=1+(k-col)/2;

axes(ax(col,row))
setFigProp(gca)
tmp=MODEL{k}.SIC_clim_sh_common(:,month);

max(tmp)

make_polar_plot_pcolor_modelmask_sh(MODEL{k}.x_sh,MODEL{k}.y_sh,MODEL{k}.ifXY_sh,tmp,radius,'',0,0,1);

tmp1=mean(MODEL{k}.siu_clim_sh_common(:,[7:9]),2);
tmp2=mean(MODEL{k}.siv_clim_sh_common(:,[7:9]),2);
SIC=mean(MODEL{k}.SIC_clim_sh_common(:,[7:9]),2);
SIC_mask=SIC>=0.3;
tmp1=tmp1.*SIC_mask;
tmp2=tmp2.*SIC_mask;

max(tmp1)
max(tmp2)

make_polar_plot_vec_modelmask_sh(MODEL{k}.x_sh,MODEL{k}.y_sh,MODEL{k}.ifXY_sh,tmp1,tmp2,radius,'',K,head_length,shaft_width,scale)

set_text_label(label_str{k},'k',0.05,0.9,12)
if(k==1)
text(-1,-0.3,'Sept SIC','FontWeight','bold','rotation',90)

end
text(-0.1,0.05,short_name{k},'FontWeight','bold','FontSize',12)

end

tmp1=mean(OBS_vel.siu_clim_sh(:,[7:9]),2);
tmp2=mean(OBS_vel.siv_clim_sh(:,[7:9]),2);

max(tmp1)
max(tmp2)

axes(ax(2,2))
setFigProp(gca)
tmp_obs=OBS_SIC{1}.SIC_clim_sh_common(:,month);
max(tmp_obs)
make_polar_plot_pcolor_modelmask_sh(OBS_SIC{1}.x_sh,OBS_SIC{1}.y_sh,OBS_SIC{1}.ifXY_sh,tmp_obs,radius,[''],0,0,1)
hold on
make_polar_plot_vec_modelmask_sh(OBS_vel.x_sh,OBS_vel.y_sh,OBS_vel.ifXY_sh,tmp1,tmp2,radius,'',K,head_length,shaft_width,scale)
set_text_label('(d)','k',0.05,0.9,12)
text(-0.15,0.05,'Observations','FontWeight','bold','FontSize',12)

       set_colorbar(0.07,0.05,1,2)

colormap(cmap)

set(gcf,'Color','w')
filename=[fig_dir,'/SIC_vel_sh.png'];
export_fig(gcf,filename,'-a1','-png','-r400','-opengl');


