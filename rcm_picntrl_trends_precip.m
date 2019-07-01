% RCM control simulation
% fra Silje sine data: '/net/o3/hymet/ssilje/EXAR/merge/'
% tas_day=getnc('/net/o3/hymet/ssilje/EXAR/merge/tas_day_EUR-44_CLMcom-CCLM5-0-6_EXAR_millennium_2501-3500.nc');

indata='/net/h2o/climphys/medhaugi/data/ch2018/cntrl/';
outfig='/home/medhaugi/ch2018/cntrl/';

period=['mam';'jja';'son';'djf';'ann'];
regs=['chne';'chnw';'chse';'chea';'chsw';'chal';'euro';'chna'];

mask2=getnc([indata '../regions_CH2018_all_EUR44.nc'],'regions');
lat =getnc([indata 'tas_mon_EUR-44_CLMcom-CCLM5-0-6_EXAR_millennium_2501-3500.nc'],'lat');
lon=getnc([indata 'tas_mon_EUR-44_CLMcom-CCLM5-0-6_EXAR_millennium_2501-3500.nc'],'lon');
mask_ls=getnc('/net/o3/hymet/ssilje/EXAR/merge/CCLM5-0-6_FR_LAND.nc','FR_LAND');
mask_ls(mask_ls>0)=1; % land-sea mask europe

% different regions
ch_ne=find(mask2==1); map_chne=lat*NaN; map_chne(ch_ne)=1;
ch_nw=find(mask2==2); map_chnw=lat*NaN; map_chnw(ch_nw)=1;
ch_se=find(mask2==3); map_chse=lat*NaN; map_chse(ch_se)=1;
ch_ea=find(mask2==4); map_chea=lat*NaN; map_chea(ch_ea)=1;
ch_sw=find(mask2==5); map_chsw=lat*NaN; map_chsw(ch_sw)=1;
ch_all=find(mask2>0); map_chal=lat*NaN;map_chal(ch_all)=1;
eu_all=find(mask_ls>0); map_euro=lat*NaN;map_euro(eu_all)=1;
map_chna=lat*NaN; map_chna(ch_ne)=1; map_chna(ch_nw)=1; % switzerland north of the alps (NE+NW-mine, CHNE+CHW-offisiell CH2018)

%%%%%%%%%
% Limit to only the first 600 years!!!!!!
%%%%%%%%%
% temperature
tas_mon=getnc([indata 'pr_mon_EUR-44_CLMcom-CCLM5-0-6_EXAR_millennium_2501-3500.nc'],'pr');
time =timenc([indata 'pr_mon_EUR-44_CLMcom-CCLM5-0-6_EXAR_millennium_2501-3500.nc'],'time',-1,-1,'noleap'); % NB: remove first time step
tas_mon=tas_mon(2:7189,:,:); time=time(2:7189,:);

% seasonal; after removal -> starts with spring
%mam jja son djf
tas_sea=getnc([indata 'pr_seasonal_EUR-44_CLMcom-CCLM5-0-6_EXAR_millennium_2501-3500.nc'],'pr');
time_sea=timenc([indata 'pr_seasonal_EUR-44_CLMcom-CCLM5-0-6_EXAR_millennium_2501-3500.nc'],'time',-1,-1,'noleap');
tas_sea=tas_sea(2:2396,:,:); time_sea=time_sea(2:2396,:);

% make annual values from monthly:
tas_ann=getnc([indata 'pr_ann_EUR-44_CLMcom-CCLM5-0-6_EXAR_millennium_2501-3500.nc'],'pr');
time_ann=timenc([indata 'pr_ann_EUR-44_CLMcom-CCLM5-0-6_EXAR_millennium_2501-3500.nc'],'time',-1,-1,'noleap');
tas_ann=tas_ann(2:600,:,:); time_ann=time_ann(2:600,:);
% tas_ann=my_mnd2ann3D(tas_mon);
% time_ann=time(1:12:end,:);


%%

[Nt,Nx,Ny]=size(tas_ann); % (Nt,Nx,Ny)

for regions=1:length(regs(:,1));
    eval(['test=reshape(map_' regs(regions,:) ',[],1);']);
    for seas=1:length(period(:,1));
        if seas <5;
            data=tas_sea(seas:4:end,:,:);
        elseif seas==5;
            data=tas_ann;
        end
        nt=length(data(:,1,1));
        map_temp=meshgrid(test,1:nt);
        map_temp=reshape(map_temp,nt,Nx,Ny);
        
        temp=data.*map_temp;
        eval(['tas_' regs(regions,:) '_' period(seas,:) '=nanmean(reshape(temp,nt,[]),2);']);
        eval(['pr_' regs(regions,:) '_' period(seas,:) '=nanmean(reshape(temp,nt,[]),2);']);
        clear temp data map_temp nt
    end
    clear test pref
end

% calculate trends
Mt_mam=nan(length(regs(:,1)),50,Nt);
Mt_jja=Mt_mam;
Mt_son=Mt_mam;
Mt_djf=Mt_mam;
Mt_ann=Mt_mam;
for seas=1:length(period(:,1));
    eval(['Mt=Mt_' period(seas,:) ';']);
    for regions=1:length(regs(:,1));
        display(['Region ' regs(regions,:) ': Period: ' period(seas,:) ]);
        eval(['ts=tas_' regs(regions,:) '_' period(seas,:) ';']);
        Nt=length(ts);
        t=1:Nt;
        % calculating trends for trend length 2-50 years
        for j=2:50;
            peri=nan(Nt-j,2);
            for i=1:Nt-j;
                peri(i,:)=polyfit(t(i:i+(j-1)),ts(i:i+(j-1))',1);
            end
            Mt(regions,j,1:1+j-1)=NaN;
            Mt(regions,j,1+j:Nt)=peri(:,1)'*j;
            clear peri
        end
        clear t ts Nt
    end
    Mt(Mt==0)=NaN;
    eval(['Mt_' period(seas,:) '=Mt;']);
end


pct=0:1:100;
for seas=1:length(period(:,1));
    eval(['Mt=Mt_' period(seas,:) ';']);
    
    for regions=1:length(regs(:,1));
%         figure(seas);
%         subplot(2,4,regions); 
%         plot(1:50,squeeze(Mt(regions,:,:)),'.k');
%         set(gca,'ylim',[-6 6]);
%         title([period(seas,:) ' ' regs(regions,:) ]);
        
        % percentiles
        prct=nan(101,50);
        for i=2:50;
            test=squeeze(Mt(regions,i,:));
            prct(:,i)=prctile(test(isnan(test)==0),pct);
            clear test
        end
        figure(seas+90);
        set(gcf,'PaperType','A4','paperOrientation','landscape','paperunits','CENTIMETERS','PaperPosition',[.63, .63, 28.41, 19.72]);
     
        subplot(2,4,regions);hold on; box on;
        %         plot(1:50,prct(pct==0,:),'.k','linewidth',2); % 0th percentile
        plot(1:50,prct(pct==5,:),'--k','linewidth',2); % 5th percentile
        plot(1:50,prct(pct==25,:),'-.k','linewidth',2,'color',[.5 .5 .5]); % 25th percentile
        plot(1:50,prct(pct==50,:),'k','linewidth',2); % 50th percentile
        plot(1:50,prct(pct==75,:),'-.k','linewidth',2,'color',[.5 .5 .5]); % 75th percentile
        plot(1:50,prct(pct==95,:),'--k','linewidth',2); % 95th percentile
        %         plot(1:50,prct(pct==100,:),'.k','linewidth',2); % 100th percentile
        
        if seas==1; 
            set(gca,'ylim',[-5.5 5.5],'xlim',[2 50]);
        elseif seas==2;
            set(gca,'ylim',[-3.5 3.5],'xlim',[2 50]);
        elseif seas==3; 
            set(gca,'ylim',[-4 4],'xlim',[2 50]);
        elseif seas==4;
            set(gca,'ylim',[-6 6],'xlim',[2 50]);
        elseif seas==5; 
            set(gca,'ylim',[-2.2 2.2],'xlim',[2 50]);
        end
        
        title([regs(regions,:) ]);
        eval(['prct_' period(seas,:) '_' regs(regions,:) '=prct;']);
        
    end
    figure(seas+90);
    subplot(2,4,8); hold on; set(gca,'Xcolor','none','Ycolor','none');
    plot(1:50,prct(pct==5,:),'--k','linewidth',2); % 5th percentile
    plot(1:50,prct(pct==25,:),'-.k','linewidth',2,'color',[.5 .5 .5]); % 25th percentile
    plot(1:50,prct(pct==50,:),'k','linewidth',2); % 50th percentile
    set(gca,'ylim',[6 8],'xlim',[2 50],'ytick',{},'xtick',{});
    legend('5th-95th prctile','25th-75th prctile','50th prctile')
    
    print('-depsc',[outfig 'regional_trend_precip_prctl_' period(seas,:) '_1_600y']);
end


peri=1:50;
% save([indata 'trend_prec_matrices_region_season_ch_eur_1_600'],'Mt_ann','Mt_djf','Mt_jja','Mt_mam','Mt_son','peri','regs');

pr_chea_ann=tas_chea_ann; pr_chea_djf=tas_chea_djf; pr_chea_jja=tas_chea_jja; pr_chea_mam=tas_chea_mam; pr_chea_son=tas_chea_son;
pr_chne_ann=tas_chne_ann; pr_chne_djf=tas_chne_djf; pr_chne_jja=tas_chne_jja; pr_chne_mam=tas_chne_mam; pr_chne_son=tas_chne_son;
pr_chnw_ann=tas_chnw_ann; pr_chnw_djf=tas_chnw_djf; pr_chnw_jja=tas_chnw_jja; pr_chnw_mam=tas_chnw_mam; pr_chnw_son=tas_chnw_son;
pr_chse_ann=tas_chse_ann;pr_chse_djf=tas_chse_djf;pr_chse_jja=tas_chse_jja;pr_chse_mam=tas_chse_mam;pr_chse_son=tas_chse_son;
pr_chsw_ann=tas_chsw_ann;pr_chsw_djf=tas_chsw_djf;pr_chsw_jja=tas_chsw_jja;pr_chsw_mam=tas_chsw_mam;pr_chsw_son=tas_chsw_son;
pr_euro_ann=tas_euro_ann;pr_euro_djf=tas_euro_djf;pr_euro_jja=tas_euro_jja;pr_euro_mam=tas_euro_mam;pr_euro_son=tas_euro_son;
pr_chal_ann=tas_chal_ann;pr_chal_djf=tas_chal_djf;pr_chal_jja=tas_chal_jja;pr_chal_mam=tas_chal_mam;pr_chal_son=tas_chal_son;
pr_chna_ann=tas_chna_ann;pr_chna_djf=tas_chna_djf;pr_chna_jja=tas_chna_jja;pr_chna_mam=tas_chna_mam;pr_chna_son=tas_chna_son;

% % save([indata 'pr_chea_1_600'],'pr_chea_ann','pr_chea_djf','pr_chea_jja','pr_chea_mam','pr_chea_son','time_ann','time_sea');
% % save([indata 'pr_chne_1_600'],'pr_chne_ann','pr_chne_djf','pr_chne_jja','pr_chne_mam','pr_chne_son','time_ann','time_sea');
% % save([indata 'pr_chnw_1_600'],'pr_chnw_ann','pr_chnw_djf','pr_chnw_jja','pr_chnw_mam','pr_chnw_son','time_ann','time_sea');
% % save([indata 'pr_chse_1_600'],'pr_chse_ann','pr_chse_djf','pr_chse_jja','pr_chse_mam','pr_chse_son','time_ann','time_sea');
% % save([indata 'pr_chsw_1_600'],'pr_chsw_ann','pr_chsw_djf','pr_chsw_jja','pr_chsw_mam','pr_chsw_son','time_ann','time_sea');
% % save([indata 'pr_euro_1_600'],'pr_euro_ann','pr_euro_djf','pr_euro_jja','pr_euro_mam','pr_euro_son','time_ann','time_sea');
% % save([indata 'pr_chal_1_600'],'pr_chal_ann','pr_chal_djf','pr_chal_jja','pr_chal_mam','pr_chal_son','time_ann','time_sea');
% % save([indata 'pr_chna_1_600'],'pr_chna_ann','pr_chna_djf','pr_chna_jja','pr_chna_mam','pr_chna_son','time_ann','time_sea');


figure;
hold on;
plot(1:50,prct(pct==5,:),'k--','linewidth',2)
plot(1:50,prct(pct==25,:),'-.','linewidth',2,'color',[.5 .5 .5])
plot(1:50,prct(pct==50,:),'k','linewidth',2)
plot(1:50,prct(pct==75,:),'-.','linewidth',2,'color',[.5 .5 .5])
plot(1:50,prct(pct==95,:),'k--','linewidth',2)
% plot(1:50,zeros(50,1),'k.','linewidth',2)
set(gca,'xlim',[1 50],'ylim',[-.4 .4],'fontsiz',16);
box on;
xlabel('Trend length (Years)','fontsize',16)
ylabel('Trend (mm/day)');
set(gcf, 'PaperPositionMode', 'manual', 'PaperType','A4','renderer','painters');
print('-depsc','-r300','/home/medhaugi/hiatus_review/trend_length_precip_5_50_95_pctile_all_years_mods_picntrl_1_600');

% save('/net/h2o/climphys/medhaugi/data/CMIP5/CMIP5_variability/CMIP5_variability_pr_trend_piControl_1_600','prct','pct');


%%%%%%%%%%
load([indata 'trend_prec_matrices_region_season_ch_eur_1_600'],'Mt_ann','Mt_djf','Mt_jja','Mt_mam','Mt_son','peri','regs');

load([indata 'pr_chea_1_600'],'pr_chea_ann','pr_chea_djf','pr_chea_jja','pr_chea_mam','pr_chea_son','time_ann','time_sea');
load([indata 'pr_chne_1_600'],'pr_chne_ann','pr_chne_djf','pr_chne_jja','pr_chne_mam','pr_chne_son','time_ann','time_sea');
load([indata 'pr_chnw_1_600'],'pr_chnw_ann','pr_chnw_djf','pr_chnw_jja','pr_chnw_mam','pr_chnw_son','time_ann','time_sea');
load([indata 'pr_chse_1_600'],'pr_chse_ann','pr_chse_djf','pr_chse_jja','pr_chse_mam','pr_chse_son','time_ann','time_sea');
load([indata 'pr_chsw_1_600'],'pr_chsw_ann','pr_chsw_djf','pr_chsw_jja','pr_chsw_mam','pr_chsw_son','time_ann','time_sea');
load([indata 'pr_euro_1_600'],'pr_euro_ann','pr_euro_djf','pr_euro_jja','pr_euro_mam','pr_euro_son','time_ann','time_sea');
load([indata 'pr_chal_1_600'],'pr_chal_ann','pr_chal_djf','pr_chal_jja','pr_chal_mam','pr_chal_son','time_ann','time_sea');
load([indata 'pr_chna_1_600'],'pr_chna_ann','pr_chna_djf','pr_chna_jja','pr_chna_mam','pr_chna_son','time_ann','time_sea');

for i=1:length(Mt_ann(:,1,1));
    eval(['Mt_ann(i,:,:)=100*Mt_ann(i,:,:)/mean(pr_' regs(i,:) '_ann);']);
    eval(['Mt_djf(i,:,:)=100*Mt_djf(i,:,:)/mean(pr_' regs(i,:) '_djf);']);
    eval(['Mt_jja(i,:,:)=100*Mt_jja(i,:,:)/mean(pr_' regs(i,:) '_jja);']);
    eval(['Mt_mam(i,:,:)=100*Mt_mam(i,:,:)/mean(pr_' regs(i,:) '_mam);']);
    eval(['Mt_son(i,:,:)=100*Mt_son(i,:,:)/mean(pr_' regs(i,:) '_son);']);
end

%% save([indata 'trend_prec_matrices_region_season_ch_eur_percentage_1_600'],'Mt_ann','Mt_djf','Mt_jja','Mt_mam','Mt_son','peri','regs');

