% Load surface tempeature data for all of the CESM simulations
%%%%%%%%%

inpath='/net/meso/climphys/cesm104/b.e104.B_1850-2000_CN.f19_g16.'; % the path to the data

% All of the different simulations
simulations=['im128.500';...
    'im128.860';...
    'im64.510 ';...
    'im128.380';...
    'im128.520';...
    'im64.390 ';...
    'im64.530 ';...
    'im128.400';...
    'im128.540';...
    'im64.410 ';...
    'im64.550 '];

% loop through the simulations and load the data
for sim=1:length(simulations(:,1))
    
    % loop through all of the files and load them
    data=nan(1872,96,144); % size of the final dataset
    time=nan(1872,2);
    i=1;
    for yr=1850:2005
        for mth=1:12
            if mth<10
                indata=[inpath deblank(simulations(sim,:)) '/archive/atm/hist/b.e104.B_1850-2000_CN.f19_g16.' deblank(simulations(sim,:)) '.cam2.h0.' num2str(yr) '-0' num2str(mth) '.nc'];
            else
                indata=[inpath deblank(simulations(sim,:)) '/archive/atm/hist/b.e104.B_1850-2000_CN.f19_g16.' deblank(simulations(sim,:)) '.cam2.h0.' num2str(yr) '-' num2str(mth) '.nc'];
            end
            data(i,:,:)=getnc(indata,'TREFHT');
            time(i,1)=yr; time(i,2)=mth;
            i=i+1;
        end
    end
    
    % convert file from monthly to annual
    data_ann=my_mnd2ann3D(data);
    time_ann=time(1:12:end,1);
    % fill in the path where you want to save the data
    save(['/net/h2o/climphys/ryserp/cesm/trefht_ann_' deblank(simulations(sim,:))],'data_ann','time_ann');
    clear data data_ann time time_ann i
end
