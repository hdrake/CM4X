%%%this function creates a standard set of sea ice mass budget diagnostic figures
%This follows the mass budget approach used by Keen et al. (2021).

function makeSeaIceDiagnosticFigs_massbudget(model,run,year_start,year_end,year_shift,short_name,SIC_clim_years,SIT_clim_years,vel_clim_years,loadFigData) 

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

fig_dir=['/home/mib/Figures/compare_sea_ice/',model{1},'/',compare_string];
if(~isdir(fig_dir))
mkdir(fig_dir)
end

%%%%%%load monthly model data and compute sea ice diagnostics

%%define minimum year
plot_year_min=1978;
for k=1:N
time_years{k}=year_start{k}+year_shift{k}:1:year_end{k}+year_shift{k};

tmp=time_years{k};
plot_year_min=min([plot_year_min min(tmp)])

end

if(loadFigData==0)

for k=1:N

MODEL{k}=struct;

%time{k}=year_start{k}+1/24:1/12:year_end{k}+1-1/24;
time_years{k}=year_start{k}+year_shift{k}:1:year_end{k}+year_shift{k};

dataDir=['/work/mib/raw_data/',model{k},'/',run{k}];

%%load SIMIP sea ice mask

load('/work/mib/raw_data/NSIDC/SIMIP_mass_region.mat','x','y','ifXY','w','ifXY_reg','reg_list')

x_SIMIP=x;
y_SIMIP=y;
w_SIMIP=w;
ifXY_SIMIP=ifXY;
ifXY_reg_SIMIP=ifXY_reg;

%%NH data
var=['simass'];
%var=['MI'];
filename=[dataDir,'/',var,'.mat']
load(filename,'data','time','x','y','ifXY','w')

size(data)

MODEL{k}.x_nh=x;
MODEL{k}.y_nh=y;
MODEL{k}.w_nh=w;
MODEL{k}.ifXY_nh=ifXY;

%%%%regrid SIMIP to GFDL grid

GFDL_gridpoints=ones(size(MODEL{k}.x_nh));
GFDL_gridpoints=GFDL_gridpoints>0;

x_SIMIP=mod(x_SIMIP+360,360);
SIMIP_gridpoints=ones(size(x_SIMIP));
SIMIP_gridpoints=SIMIP_gridpoints>0;
SIMIP_land=~ifXY_SIMIP;

nReg=size(ifXY_reg_SIMIP,3);

nD=nnz(MODEL{k}.ifXY_nh);
regrid_data=zeros(size(MODEL{k}.x_nh,1),size(MODEL{k}.x_nh,2),nReg);

for i=1:nReg

A=ifXY_reg_SIMIP(:,:,i);

dataTmp=A(SIMIP_gridpoints);

xq=MODEL{k}.x_nh(GFDL_gridpoints);
yq=MODEL{k}.y_nh(GFDL_gridpoints);

x=x_SIMIP(SIMIP_gridpoints);
y=y_SIMIP(SIMIP_gridpoints);

disp('Computing interpolant...')
F = scatteredInterpolant(x,y,dataTmp,'nearest','nearest');
vq=F(xq,yq);

disp('Done!')

A=zeros(size(MODEL{k}.x_nh));
A(GFDL_gridpoints)=vq;

regrid_data(:,:,i)=A;
i
end

MODEL{k}.ifXY_reg_nh=regrid_data;

%%%%done regridding

[MODEL{k}.MI_nh regionalSumclim regionalSumsig]=computeRegionalSum(data,MODEL{k}.w_nh,MODEL{k}.ifXY_nh,MODEL{k}.ifXY_reg_nh);

var=['sisnmass'];
filename=[dataDir,'/',var,'.mat']
load(filename,'data')

size(data)

[MODEL{k}.sisnmass_nh regionalSumclim regionalSumsig]=computeRegionalSum(data,MODEL{k}.w_nh,MODEL{k}.ifXY_nh,MODEL{k}.ifXY_reg_nh);

%note: TSNK shows errors when computed as LSNK-BSNK. There appears to be an inconsistency in the way that the BSNK diagnostic is computed (sometimes TSNK>0, which should never happen). BSNK_RES computes BSNK as a residual: BSNK_RES=LSNK-TMELT_MASS
%var=['BSNK'];
var=['BSNK_RES'];
filename=[dataDir,'/',var,'.mat']
load(filename,'data')

size(data)

[MODEL{k}.BSNK_nh regionalSumclim regionalSumsig]=computeRegionalSum(data,MODEL{k}.w_nh,MODEL{k}.ifXY_nh,MODEL{k}.ifXY_reg_nh);

%note: TSNK shows errors when computed as LSNK-BSNK. There appears to be an inconsistency in the way that the BSNK diagnostic is computed (sometimes TSNK>0, which should never happen). Instead use TMELT_MASS, which computes TSNK using the TMELT diagnostic.
%var=['TSNK'];
var=['TMELT_MASS'];
filename=[dataDir,'/',var,'.mat']
load(filename,'data')

size(data)

[MODEL{k}.TSNK_nh regionalSumclim regionalSumsig]=computeRegionalSum(data,MODEL{k}.w_nh,MODEL{k}.ifXY_nh,MODEL{k}.ifXY_reg_nh);

var=['LSRC'];
filename=[dataDir,'/',var,'.mat']
load(filename,'data')

size(data)

[MODEL{k}.LSRC_nh regionalSumclim regionalSumsig]=computeRegionalSum(data,MODEL{k}.w_nh,MODEL{k}.ifXY_nh,MODEL{k}.ifXY_reg_nh);

var=['XPRT'];
filename=[dataDir,'/',var,'.mat']
load(filename,'data')

size(data)

[MODEL{k}.XPRT_nh regionalSumclim regionalSumsig]=computeRegionalSum(data,MODEL{k}.w_nh,MODEL{k}.ifXY_nh,MODEL{k}.ifXY_reg_nh);

[MODEL{k}.SIMIP_area regionalSumclim regionalSumsig]=computeRegionalSum(ones(size(MODEL{k}.w_nh)),MODEL{k}.w_nh,MODEL{k}.ifXY_nh,MODEL{k}.ifXY_reg_nh);

%compute clim over SIC satellite era
ind1=(SIC_clim_years(1)-year_start{k}-year_shift{k})*12+1
ind2=(SIC_clim_years(2)+1-year_start{k}-year_shift{k})*12

[MODEL{k}.LSRC_nh_clim MODEL{k}.LSRC_nh_anom data_sig]=computeclim1D(MODEL{k}.LSRC_nh(ind1:ind2));
[MODEL{k}.XPRT_nh_clim MODEL{k}.XPRT_nh_anom data_sig]=computeclim1D(MODEL{k}.XPRT_nh(ind1:ind2));
[MODEL{k}.BSNK_nh_clim MODEL{k}.BSNK_nh_anom data_sig]=computeclim1D(MODEL{k}.BSNK_nh(ind1:ind2));
[MODEL{k}.TSNK_nh_clim MODEL{k}.TSNK_nh_anom data_sig]=computeclim1D(MODEL{k}.TSNK_nh(ind1:ind2));
[MODEL{k}.MI_nh_clim MODEL{k}.MI_nh_anom data_sig]=computeclim1D(MODEL{k}.MI_nh(ind1:ind2));
[MODEL{k}.sisnmass_nh_clim MODEL{k}.sisnmass_nh_anom data_sig]=computeclim1D(MODEL{k}.sisnmass_nh(ind1:ind2));

[MODEL{k}.LSRC_nh_annualMean]=computeAnnualMean1D(MODEL{k}.LSRC_nh);
[MODEL{k}.XPRT_nh_annualMean]=computeAnnualMean1D(MODEL{k}.XPRT_nh);
[MODEL{k}.BSNK_nh_annualMean]=computeAnnualMean1D(MODEL{k}.BSNK_nh);
[MODEL{k}.TSNK_nh_annualMean]=computeAnnualMean1D(MODEL{k}.TSNK_nh);
[MODEL{k}.MI_nh_annualMean]=computeAnnualMean1D(MODEL{k}.MI_nh);
[MODEL{k}.sisnmass_nh_annualMean]=computeAnnualMean1D(MODEL{k}.sisnmass_nh);

%%SH data

var=['simass_sh'];
%var=['MI_sh'];
filename=[dataDir,'/',var,'.mat']
load(filename,'data','time','x','y','ifXY','w')

size(data)

MODEL{k}.x_sh=x;
MODEL{k}.y_sh=y;
MODEL{k}.w_sh=w;
MODEL{k}.ifXY_sh=ifXY;

%%%%define mass budget region

lat_lim=-63;

MODEL{k}.ifXY_reg_sh = (ifXY & (y<lat_lim));

[MODEL{k}.MI_sh regionalSumclim regionalSumsig]=computeRegionalSum(data,MODEL{k}.w_sh,MODEL{k}.ifXY_sh,MODEL{k}.ifXY_reg_sh);

var=['sisnmass_sh'];
filename=[dataDir,'/',var,'.mat']
load(filename,'data')

size(data)

[MODEL{k}.sisnmass_sh regionalSumclim regionalSumsig]=computeRegionalSum(data,MODEL{k}.w_sh,MODEL{k}.ifXY_sh,MODEL{k}.ifXY_reg_sh);

%var=['BSNK_sh'];
var=['BSNK_RES_sh'];
filename=[dataDir,'/',var,'.mat']
load(filename,'data')

size(data)

[MODEL{k}.BSNK_sh regionalSumclim regionalSumsig]=computeRegionalSum(data,MODEL{k}.w_sh,MODEL{k}.ifXY_sh,MODEL{k}.ifXY_reg_sh);

%var=['TSNK_sh'];
var=['TMELT_MASS_sh'];
filename=[dataDir,'/',var,'.mat']
load(filename,'data')

size(data)

[MODEL{k}.TSNK_sh regionalSumclim regionalSumsig]=computeRegionalSum(data,MODEL{k}.w_sh,MODEL{k}.ifXY_sh,MODEL{k}.ifXY_reg_sh);

var=['LSRC_sh'];
filename=[dataDir,'/',var,'.mat']
load(filename,'data')

size(data)

[MODEL{k}.LSRC_sh regionalSumclim regionalSumsig]=computeRegionalSum(data,MODEL{k}.w_sh,MODEL{k}.ifXY_sh,MODEL{k}.ifXY_reg_sh);

var=['XPRT_sh'];
filename=[dataDir,'/',var,'.mat']
load(filename,'data')

size(data)

[MODEL{k}.XPRT_sh regionalSumclim regionalSumsig]=computeRegionalSum(data,MODEL{k}.w_sh,MODEL{k}.ifXY_sh,MODEL{k}.ifXY_reg_sh);

[MODEL{k}.sh_area regionalSumclim regionalSumsig]=computeRegionalSum(ones(size(MODEL{k}.w_sh)),MODEL{k}.w_sh,MODEL{k}.ifXY_sh,MODEL{k}.ifXY_reg_sh);

%compute clim over SIC satellite era
ind1=(SIC_clim_years(1)-year_start{k}-year_shift{k})*12+1
ind2=(SIC_clim_years(2)+1-year_start{k}-year_shift{k})*12

[MODEL{k}.LSRC_sh_clim MODEL{k}.LSRC_sh_anom data_sig]=computeclim1D(MODEL{k}.LSRC_sh(ind1:ind2));
[MODEL{k}.XPRT_sh_clim MODEL{k}.XPRT_sh_anom data_sig]=computeclim1D(MODEL{k}.XPRT_sh(ind1:ind2));
[MODEL{k}.BSNK_sh_clim MODEL{k}.BSNK_sh_anom data_sig]=computeclim1D(MODEL{k}.BSNK_sh(ind1:ind2));
[MODEL{k}.TSNK_sh_clim MODEL{k}.TSNK_sh_anom data_sig]=computeclim1D(MODEL{k}.TSNK_sh(ind1:ind2));
[MODEL{k}.MI_sh_clim MODEL{k}.MI_sh_anom data_sig]=computeclim1D(MODEL{k}.MI_sh(ind1:ind2));
[MODEL{k}.sisnmass_sh_clim MODEL{k}.sisnmass_sh_anom data_sig]=computeclim1D(MODEL{k}.sisnmass_sh(ind1:ind2));

[MODEL{k}.LSRC_sh_annualMean]=computeAnnualMean1D(MODEL{k}.LSRC_sh);
[MODEL{k}.XPRT_sh_annualMean]=computeAnnualMean1D(MODEL{k}.XPRT_sh);
[MODEL{k}.BSNK_sh_annualMean]=computeAnnualMean1D(MODEL{k}.BSNK_sh);
[MODEL{k}.TSNK_sh_annualMean]=computeAnnualMean1D(MODEL{k}.TSNK_sh);
[MODEL{k}.MI_sh_annualMean]=computeAnnualMean1D(MODEL{k}.MI_sh);
[MODEL{k}.sisnmass_sh_annualMean]=computeAnnualMean1D(MODEL{k}.sisnmass_sh);


end

%%%%SAVE ALL DATA NEEDED TO MAKE FIGURES

work_dir=['/work/mib/processed_data/compare_sea_ice/',model{1},'/',compare_string];
filename=[work_dir,'/DiagnosticFigData_massbudget.mat'];

save(filename,'MODEL')

elseif(loadFigData==1)

work_dir=['/work/mib/processed_data/compare_sea_ice/',model{1},'/',compare_string];
filename=[work_dir,'/DiagnosticFigData_massbudget.mat'];

load(filename,'MODEL')

end

%%%%%%MAKE FIGURES

colorlist={'b','r',rgb('LimeGreen'),'c','m',rgb('Gray'),rgb('Indigo')};

%%%%%%%FIGURE 1: NH mass budget

    nTileX=2;
    nTileY = 2;
        
    figWidth = 7; %6.5 inches = 39 picas, for a two column figure
    deltaX   = 0.6;
    deltaX2  = 0.2;
    deltaY   = 0.5;
    deltaY2  = 0.4;
    gapX     = 0.7;
    gapY     = 0.9;
    panel_scale=0.8;

ax=make_axes(nTileX,nTileY,figWidth,deltaX,deltaX2,deltaY,deltaY2,gapX,gapY,panel_scale);


fact=1/(905*MODEL{k}.SIMIP_area);
L=2;

axes(ax(1,1))
setFigProp(gca)

hold on
for k=1:N
plot(1:12,MODEL{k}.LSRC_nh_clim*fact,'Color',colorlist{k},'LineWidth',L)
end
xlabel('Month')
ylabel('Thickness tendency (m/yr)')
xlim([0.5 12.5])
ylim([0 4])
grid on
title(['Ice growth climatology'])
set_text_label('(a)','k',-0.1,1.1,12)
set(gca,'xTick',1:12,'xTickLabel',{'J','F','M','A','M','J','J','A','S','O','N','D'},'XTickLabelRotation',0)

axes(ax(2,1))
setFigProp(gca)

hold on
for k=1:N
plot(1:12,MODEL{k}.XPRT_nh_clim*fact,'Color',colorlist{k},'LineWidth',L)
end
xlabel('Month')
ylabel('Thickness tendency (m/yr)')
xlim([0.5 12.5])
ylim([-2 2])
grid on
title(['Ice export climatology'])
set_text_label('(b)','k',-0.1,1.1,12)
set(gca,'xTick',1:12,'xTickLabel',{'J','F','M','A','M','J','J','A','S','O','N','D'},'XTickLabelRotation',0)
legend([short_name],'location','NorthWest')

axes(ax(1,2))
setFigProp(gca)

hold on
for k=1:N
plot(1:12,MODEL{k}.BSNK_nh_clim*fact,'Color',colorlist{k},'LineWidth',L)
end
xlabel('Month')
ylabel('Thickness tendency (m/yr)')
xlim([0.5 12.5])
ylim([-3.8 0.2])
grid on
title(['Basal melt climatology'])
set_text_label('(c)','k',-0.1,1.1,12)
set(gca,'xTick',1:12,'xTickLabel',{'J','F','M','A','M','J','J','A','S','O','N','D'},'XTickLabelRotation',0)

axes(ax(2,2))
setFigProp(gca)

hold on
for k=1:N
plot(1:12,MODEL{k}.TSNK_nh_clim*fact,'Color',colorlist{k},'LineWidth',L)
end
xlabel('Month')
ylabel('Thickness tendency (m/yr)')
xlim([0.5 12.5])
ylim([-3.8 0.2])
grid on
title(['Top melt climatology'])
set_text_label('(d)','k',-0.1,1.1,12)
set(gca,'xTick',1:12,'xTickLabel',{'J','F','M','A','M','J','J','A','S','O','N','D'},'XTickLabelRotation',0)

tmp=MODEL{1}.ifXY_reg_nh;
hold on
%axes('Position',[0.075 0.09 0.17 0.15])
axes('Position',[0.57 0.09 0.17 0.15])
box on
fact=1e6;
make_polar_plot_pcolor(MODEL{1}.x_nh,MODEL{1}.y_nh,MODEL{1}.ifXY_nh,fact*tmp(MODEL{1}.ifXY_nh),35,[''],0,-1,1);
box on


set(gcf,'Color','w')
filename=[fig_dir,'/massbudget_clim_nh'];
%export_fig(gcf,filename,'-a1','-pdf');
%export_fig(gcf,filename,'-a1','-png','-r400','-opengl');
exportgraphics(gcf,[filename,'.png'],'Resolution',400)
exportgraphics(gcf,[filename,'.pdf'],'Resolution',400)

%%%%%%%FIGURE 2: SH mass budget

    nTileX=2;
    nTileY = 2;
        
    figWidth = 7; %6.5 inches = 39 picas, for a two column figure
    deltaX   = 0.6;
    deltaX2  = 0.2;
    deltaY   = 0.5;
    deltaY2  = 0.4;
    gapX     = 0.7;
    gapY     = 0.9;
    panel_scale=0.8;

ax=make_axes(nTileX,nTileY,figWidth,deltaX,deltaX2,deltaY,deltaY2,gapX,gapY,panel_scale);

fact=1/(905*MODEL{k}.sh_area);

axes(ax(1,1))
setFigProp(gca)

hold on
for k=1:N
plot(1:12,MODEL{k}.LSRC_sh_clim*fact,'Color',colorlist{k},'LineWidth',L)
end
xlabel('Month')
ylabel('Thickness tendency (m/yr)')
xlim([0.5 12.5])
ylim([0 4.0])
grid on
title(['Ice growth climatology'])
set_text_label('(a)','k',-0.1,1.1,12)
set(gca,'xTick',1:12,'xTickLabel',{'J','F','M','A','M','J','J','A','S','O','N','D'},'XTickLabelRotation',0)

axes(ax(2,1))
setFigProp(gca)

hold on
for k=1:N
plot(1:12,MODEL{k}.XPRT_sh_clim*fact,'Color',colorlist{k},'LineWidth',L)
end
xlabel('Month')
ylabel('Thickness tendency (m/yr)')
xlim([0.5 12.5])
ylim([-2 2])
grid on
title(['Ice export climatology'])
set_text_label('(b)','k',-0.1,1.1,12)
legend([short_name],'location','NorthWest')
set(gca,'xTick',1:12,'xTickLabel',{'J','F','M','A','M','J','J','A','S','O','N','D'},'XTickLabelRotation',0)

axes(ax(1,2))
setFigProp(gca)

hold on
for k=1:N
plot(1:12,MODEL{k}.BSNK_sh_clim*fact,'Color',colorlist{k},'LineWidth',L)
end
xlabel('Month')
ylabel('Thickness tendency (m/yr)')
xlim([0.5 12.5])
ylim([-3.8 0.2])
grid on
title(['Basal melt climatology'])
set_text_label('(c)','k',-0.1,1.1,12)
set(gca,'xTick',1:12,'xTickLabel',{'J','F','M','A','M','J','J','A','S','O','N','D'},'XTickLabelRotation',0)

axes(ax(2,2))
setFigProp(gca)

hold on
for k=1:N
plot(1:12,MODEL{k}.TSNK_sh_clim*fact,'Color',colorlist{k},'LineWidth',L)
end
xlabel('Month')
ylabel('Thickness tendency (m/yr)')
xlim([0.5 12.5])
ylim([-3.8 0.2])
grid on
title(['Top melt climatology'])
set_text_label('(d)','k',-0.1,1.1,12)
set(gca,'xTick',1:12,'xTickLabel',{'J','F','M','A','M','J','J','A','S','O','N','D'},'XTickLabelRotation',0)
tmp=MODEL{1}.ifXY_reg_sh;
hold on
%axes('Position',[0.075 0.09 0.17 0.15])
axes('Position',[0.57 0.09 0.17 0.15])
box on
fact=1e6;
make_polar_plot_pcolor_sh(MODEL{1}.x_sh,MODEL{1}.y_sh,MODEL{1}.ifXY_sh,fact*tmp(MODEL{1}.ifXY_sh),35,[''],0,-1,1);
box on

set(gcf,'Color','w')
filename=[fig_dir,'/massbudget_clim_sh'];
%export_fig(gcf,filename,'-a1','-pdf');
%export_fig(gcf,filename,'-a1','-png','-r400','-opengl');
exportgraphics(gcf,[filename,'.png'],'Resolution',400)
exportgraphics(gcf,[filename,'.pdf'],'Resolution',400)

%%%%%%FIGURE 2: NH Mass Budget timeseries

    nTileX=3;
    nTileY = 2;
        
    figWidth = 9.5; %6.5 inches = 39 picas, for a two column figure
    deltaX   = 0.7;
    deltaX2  = 0.2;
    deltaY   = 0.5;
    deltaY2  = 0.3;
    gapX     = 0.6;
    gapY     = 1.4;
    panel_scale=0.8;

ax=make_axes(nTileX,nTileY,figWidth,deltaX,deltaX2,deltaY,deltaY2,gapX,gapY,panel_scale);

fact=1/(905*MODEL{k}.SIMIP_area);
fact_snow=1/(330*MODEL{k}.SIMIP_area);
L=2;

axes(ax(1,1))
setFigProp(gca)
for k=1:N
plot(time_years{k},fact*MODEL{k}.LSRC_nh_annualMean,'Color',colorlist{k},'LineWidth',L)
end

grid on
l=legend([short_name],'location','SouthWest','orientation','horizontal');
p=get(l,'position');
set(l,'position',[p(1) p(2)-0.15 p(3) p(4)]);

xlim([plot_year_min 2020])
xlabel('Time (years)')
ylabel('Thickness tendency (m/yr)')
title(['Ice Growth'])

axes(ax(2,1))
setFigProp(gca)
for k=1:N
plot(time_years{k},fact*MODEL{k}.XPRT_nh_annualMean,'Color',colorlist{k},'LineWidth',L)
end

grid on

xlim([plot_year_min 2020])
xlabel('Time (years)')
title(['Ice Export'])

axes(ax(1,2))
setFigProp(gca)
for k=1:N
plot(time_years{k},fact*MODEL{k}.BSNK_nh_annualMean,'Color',colorlist{k},'LineWidth',L)
end

grid on

xlim([plot_year_min 2020])
xlabel('Time (years)')
ylabel('Thickness tendency (m/yr)')
title(['Basal Melt'])

axes(ax(2,2))
setFigProp(gca)
for k=1:N
plot(time_years{k},fact*MODEL{k}.TSNK_nh_annualMean,'Color',colorlist{k},'LineWidth',L)
end

grid on

xlim([plot_year_min 2020])
xlabel('Time (years)')
title(['Top Melt'])

axes(ax(3,1))
setFigProp(gca)
for k=1:N
plot(time_years{k},fact*MODEL{k}.MI_nh_annualMean,'Color',colorlist{k},'LineWidth',L)
end

grid on

xlim([plot_year_min 2020])
xlabel('Time (years)')
ylabel('Thickness (m)')
title(['Ice Thickness'])

axes(ax(3,2))
setFigProp(gca)
for k=1:N
plot(time_years{k},fact_snow*MODEL{k}.sisnmass_nh_annualMean,'Color',colorlist{k},'LineWidth',L)
end

grid on

xlim([plot_year_min 2020])
xlabel('Time (years)')
ylabel('Thickness (m)')
title(['Snow Thickness'])


set(gcf,'Color','w')
filename=[fig_dir,'/massbudget_interannual_nh'];
%export_fig(gcf,filename,'-a1','-pdf');
%export_fig(gcf,filename,'-a1','-png','-r400','-opengl');
exportgraphics(gcf,[filename,'.png'],'Resolution',400)
exportgraphics(gcf,[filename,'.pdf'],'Resolution',400)

%%%%%%FIGURE : SH Mass Budget timeseries

    nTileX=3;
    nTileY = 2;
        
    figWidth = 9.5; %6.5 inches = 39 picas, for a two column figure
    deltaX   = 0.7;
    deltaX2  = 0.2;
    deltaY   = 0.5;
    deltaY2  = 0.3;
    gapX     = 0.6;
    gapY     = 1.4;
    panel_scale=0.8;

ax=make_axes(nTileX,nTileY,figWidth,deltaX,deltaX2,deltaY,deltaY2,gapX,gapY,panel_scale);

fact=1/(905*MODEL{k}.sh_area);
fact_snow=1/(330*MODEL{k}.sh_area);

axes(ax(1,1))
setFigProp(gca)
for k=1:N
plot(time_years{k},fact*MODEL{k}.LSRC_sh_annualMean,'Color',colorlist{k},'LineWidth',L)
end

grid on
l=legend([short_name],'location','SouthWest','orientation','horizontal');
p=get(l,'position');
set(l,'position',[p(1) p(2)-0.15 p(3) p(4)]);

xlim([plot_year_min 2020])
xlabel('Time (years)')
ylabel('Thickness tendency (m/yr)')
title(['Ice Growth'])

axes(ax(2,1))
setFigProp(gca)
for k=1:N
plot(time_years{k},fact*MODEL{k}.XPRT_sh_annualMean,'Color',colorlist{k},'LineWidth',L)
end

grid on

xlim([plot_year_min 2020])
xlabel('Time (years)')
title(['Ice Export'])

axes(ax(1,2))
setFigProp(gca)
for k=1:N
plot(time_years{k},fact*MODEL{k}.BSNK_sh_annualMean,'Color',colorlist{k},'LineWidth',L)
end

grid on

xlim([plot_year_min 2020])
xlabel('Time (years)')
ylabel('Thickness tendency (m/yr)')
title(['Basal Melt'])

axes(ax(2,2))
setFigProp(gca)
for k=1:N
plot(time_years{k},fact*MODEL{k}.TSNK_sh_annualMean,'Color',colorlist{k},'LineWidth',L)
end

grid on

xlim([plot_year_min 2020])
xlabel('Time (years)')
title(['Top Melt'])

axes(ax(3,1))
setFigProp(gca)
for k=1:N
plot(time_years{k},fact*MODEL{k}.MI_sh_annualMean,'Color',colorlist{k},'LineWidth',L)
end

grid on

xlim([plot_year_min 2020])
xlabel('Time (years)')
ylabel('Thickness (m)')
title(['Ice Thickness'])

axes(ax(3,2))
setFigProp(gca)
for k=1:N
plot(time_years{k},fact_snow*MODEL{k}.sisnmass_sh_annualMean,'Color',colorlist{k},'LineWidth',L)
end

grid on

xlim([plot_year_min 2020])
xlabel('Time (years)')
ylabel('Thickness (m)')
title(['Snow Thickness'])

set(gcf,'Color','w')
filename=[fig_dir,'/massbudget_interannual_sh'];
%export_fig(gcf,filename,'-a1','-pdf');
%export_fig(gcf,filename,'-a1','-png','-r400','-opengl');
exportgraphics(gcf,[filename,'.png'],'Resolution',400)
exportgraphics(gcf,[filename,'.pdf'],'Resolution',400)






