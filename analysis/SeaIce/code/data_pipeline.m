
model{1}='CM4p125_new';
run{1}='CM4_c192_OM4p125_historical_Scenario_ssp585';
year_start{1}=1850;
year_end{1}=2099;
year_shift{1}=0;%finish year is 2015
short_name{1}='CM4X-p125';

model{2}='CM4p25';
run{2}='CM4_c192_OM4p25_historical_Scenario_ssp585';
year_start{2}=1850;
year_end{2}=2099;
year_shift{2}=0;%finish year is 2015
short_name{2}='CM4X-p25';

model{3}='CM4p25';
run{3}='CM4_historical_Scenario_ssp585';
year_start{3}=1850;
year_end{3}=2099;
year_shift{3}=0;%finish year is 2015
short_name{3}='CM4.0';

SIC_clim_years=[1979 2014];
SIT_clim_years=[2005 2014];
vel_clim_years=[2005 2014];

isPI=0;
loadFigData=1;

makeSeaIceDiagnosticFigsCM4X(model,run,year_start,year_end,year_shift,short_name,SIC_clim_years,SIT_clim_years,vel_clim_years,isPI,loadFigData)
makeSeaIceDiagnosticFigs_massbudget(model,run,year_start,year_end,year_shift,short_name,SIC_clim_years,SIT_clim_years,vel_clim_years,loadFigData)


