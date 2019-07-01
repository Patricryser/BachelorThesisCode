function [data]=my_mnd2ann3D(ts);
s=size(ts);
data=nan(int64(s(1)/12),s(2),s(3));
for i=1:length(ts(:,1,1))/12;
    data(i,:,:)=nanmean(ts(i*12-11:i*12,:,:));
end