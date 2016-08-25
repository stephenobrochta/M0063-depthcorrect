function [data]=depthcorrect(data,hole,varargin)

% specify depths (y values may be attached), 'hole' as string. add a 1  to plot
%
%[ambsf]=depthcorrect(data,'C',0)
%returns ambsf for data from hole c and plots a figure showing the adjustment of each data point.
%obviously plotting with large datasets is a bad idea.
%
%returns ambsf, 1 sigma, and core top and bottom depths. These are needed as bound on 1 sigma.

if nargin>2
	plotme=1;
else
	plotme=0;
end

%load private/data.mat
s=load(['private/63' upper(hole) 'summary.txt']);
n=1000;

core=data(:,1);
x=data(:,2);

%build all expansion profiles which is simpler to index when coverting discreet samples
for i=1:length(s)
	mbsf_top(i,1)=s(i,2);
	mbsf_bot_corr(i,1)=s(i,2)+s(i,3);
	mbsf_bot(i,1)=s(i,2)+s(i,4);
	
	%use a logspace vector instead of the above polyfit because can precisely control the start and end values
	d_orig(:,i)=logspace(0,1,n);
	d_orig(:,i)=mbsf_top(i)+((mbsf_bot(i)-mbsf_top(i))*(d_orig(:,i)-min(d_orig(:,i))))/(max(d_orig(:,i))-min(d_orig(:,i)));
	d_corr(:,i)=linspace(mbsf_top(i),mbsf_bot_corr(i),n);
end

x_corr=[];
for i=1:length(x)
	x_corr=[x_corr; interp1(d_orig(:,core(i)),d_corr(:,core(i)),x(i)) (s(core(i),3)*.1) mbsf_top(core(i)) mbsf_bot_corr(core(i))];
	if plotme==1;
		plusonesig=d_corr(:,core(i))+(s(core(i),3)*.1);
		plusonesig(plusonesig>mbsf_bot_corr(core(i)))=mbsf_bot_corr(core(i));
		minusonesig=d_corr(:,core(i))-(s(core(i),3)*.1);
		minusonesig(minusonesig<mbsf_top(core(i)))=mbsf_top(core(i));
		figure
		hold(gca,'on')
		plot(d_orig(:,core(i)),d_corr(:,core(i)))
		plot(d_orig(:,core(i)),plusonesig)
		plot(d_orig(:,core(i)),minusonesig)
		plot(x(i),x_corr(i,1),'or')
		axis tight
		xlabel('Depth (mbsf)')
		ylabel('Adjusted Depth (ambsf)')
		set(gca,'tickdir','out','box','on','xtick',0:.25:max(mbsf_bot(end,end)),'ytick',0:.25:max(mbsf_bot_corr(end,end)))
		title({['M0063' upper(hole) '-' num2str(core(i))] []},'fontsize',12)
		line([x(i) x(i)],[d_corr(1,core(i)) d_corr(end,core(i))],'color','r')
		line([d_orig(1,core(i)) d_orig(end,core(i))],[x_corr(i) x_corr(i)],'color','r')
		grid on
	end
end

data=[x_corr data(:,3:end)];
