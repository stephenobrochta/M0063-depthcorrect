function [data,varargout]=depthcorrect(core,data,hole,varargin)

% Input: core number, mbsf depth, and hole. Add a 1 to plot a figure for each core.
% Output: adjusted mbsf (ambsf) with 1 sigma range.
% Optional output: adjusted mcd (amcd) with 1 sigma range.
%
% No adjustment is performed for cores with <=100% recovery.
%
% Note that no composite splice exists Site M0063.
% mcd determined in Correlator using magnetic susceptibility and is an estimate
% Only use mcd for the purposes of transferring age models between holes.
% mcd uncertainty includes an additional 10% that must be considered when assigning age
%
% Convert 2.5 mbsf in Core M0063C-1H and plot a figure showing the adjustment:
% ambsf=depthcorrect(1, 2.5,'C',1)
%
% Core number required because 2.5 mbsf occurs in Core M0063C-1H and Core M0063C-2H
%
% Calcuate composite depth based on drilling offsets
% [ambsf,amcd]=depthcorrect(1, 2.5,'C',1)
%
% Convert a range of depths from one core. Example does not plot
% ambsf=depthcorrect(1, [2.5; 2.51; 2.52],'C')
%
% Convert a range of depths from more than one core. Example does not plot
% ambsf=depthcorrect([1; 1; 1; 2; 2; 2], [2.5; 2.51; 2.52; 2.5; 2.51; 2.52],'C')
%
% Obviously plotting with large datasets is a bad idea.
%
% S.P. Obrochta, Sept. 15, 2016

s=load(['private/63' upper(hole) 'summary.txt']);

for i=1:length(core)
	if isempty(find(core(i)==s(:,1),1))==1
		error(['Core M0063' upper(hole) '-' num2str(core(i)) ' is not included in curated database. Are you sure it exists?'])
	end
end

if length(core)==1
	core=core*ones(length(data),1);
elseif length(core)~=length(data)
	error(['There are ' num2str(length(data)) ' depths input but only ' num2str(length(core)) ' cores specified. If each depth is from the same core, input the core number as a scalar. If the depths are from multiple cores, enter one core number per depth.'])
end

x=data(:,1);
mbsf_top=s(:,2);
pen=s(:,3);
mbsf_bot_corr=s(:,2)+s(:,3);
rec=s(:,4);
mbsf_bot=mbsf_top+rec;

%fit parameters
a_mean=0.3330;
b_mean=0.6437;
c_mean=11.6602;
d_mean=0.3631;

a_minus=0.1203;
b_minus=0.7501;
c_minus=11.8222;
d_minus=0.7065;

a_plus=0.4857;
b_plus=0.5016;
c_plus=14.5016;
d_plus=0.1844;

fun=@(x,a,b,c,d,n) (a+b*log10((1/n)*c*x+d));
ambsf=NaN(length(x),3);
skip=[];
for i=1:length(x)
	ind=s(:,1)==core(i);
	if x(i)>=mbsf_top(ind) && x(i)<=mbsf_bot(ind)
		if rec(ind)>pen(ind)
			% Xscale is the scaling factor of the X-axis
			Xscale=rec(ind);
			ambsf(i,1)=mbsf_top(ind)+pen(ind)*integral(@(x)fun(x,a_mean,b_mean,c_mean,d_mean,Xscale),0,x(i)-mbsf_top(ind))/integral(@(x)fun(x,a_mean,b_mean,c_mean,d_mean,Xscale),0,rec(ind));
			ambsf(i,2)=mbsf_top(ind)+pen(ind)*integral(@(x)fun(x,a_minus,b_minus,c_minus,d_minus,Xscale),0,x(i)-mbsf_top(ind))/integral(@(x)fun(x,a_minus,b_minus,c_minus,d_minus,Xscale),0,rec(ind));
			ambsf(i,3)=mbsf_top(ind)+pen(ind)*integral(@(x)fun(x,a_plus,b_plus,c_plus,d_plus,Xscale),0,x(i)-mbsf_top(ind))/integral(@(x)fun(x,a_plus,b_plus,c_plus,d_plus,Xscale),0,rec(ind));
		else
			ambsf(i,1)=x(i);
			ambsf(i,2:3)=NaN;
			skip=[skip core(i)];
		end
	else
		ambsf(i,1)=NaN;
	end
	%check for 1 sigma exceeding core top or bottom
	if ambsf(i,2)<mbsf_top(core(i))
		ambsf(i,2)=mbsf_top(core(i));
	end
	if ambsf(i,3)>mbsf_bot_corr(core(i))
		ambsf(i,3)=mbsf_bot_corr(core(i));
	end
end

data=[ambsf data(:,3:end)];

%calculate mcd based on offsets
if nargout==2
	switch strcmpi(hole,'A')==1
		case 1
			warning('No published offsets for Hole A. Calculating ambsf only')
			varargout={NaN(size(data,1),1)};
		case 0
			o=load(['private/63' upper(hole) 'offset.txt']);
			amcd=ambsf;
			for i=1:size(amcd,1)
				amcd(i,:)=amcd(i,:)+o(core(i)==o(:,1),2);
			end
			varargout={amcd};
	end
else
	varargout={[]};
end

if isempty(skip)==0
skip=unique(skip);
	for i=1:length(skip)
		skipstr{i}=['M0063' upper(hole) '-' num2str(skip(i)) '; '];
	end
	skipstr{end}=strrep(skipstr{end},'; ','.');
	warning(['The following cores were not adjusted because recovery is <=100% :' skipstr{:}])
end

if nargin>3
	plotcore=unique(core);
	for i=1:length(plotcore)
		ind=s(:,1)==plotcore(i);
		xline=[mbsf_top(ind):.01:mbsf_bot(ind)];
		Xscale=rec(ind);
		clear linemid linelo lineup
		for j=1:length(xline)
			interval(j,1)=xline(j)-mbsf_top(ind);
			linemid(j,1)=mbsf_top(ind)+pen(ind)*integral(@(x)fun(x,a_mean,b_mean,c_mean,d_mean,Xscale),0,xline(j)-mbsf_top(ind))/integral(@(x)fun(x,a_mean,b_mean,c_mean,d_mean,Xscale),0,rec(ind));
			linelo(j,1)=mbsf_top(ind)+pen(ind)*integral(@(x)fun(x,a_minus,b_minus,c_minus,d_minus,Xscale),0,xline(j)-mbsf_top(ind))/integral(@(x)fun(x,a_minus,b_minus,c_minus,d_minus,Xscale),0,rec(ind));
			lineup(j,1)=mbsf_top(ind)+pen(ind)*integral(@(x)fun(x,a_plus,b_plus,c_plus,d_plus,Xscale),0,xline(j)-mbsf_top(ind))/integral(@(x)fun(x,a_plus,b_plus,c_plus,d_plus,Xscale),0,rec(ind));
		end
		plotind=find(core==s(ind,1));
		figure
		axes
		hold(gca,'on')
		plot(xline,linemid,'k')
		plot(xline,[linelo lineup],'linestyle','--','color',[.5 .5 .5])
		plot(x(plotind),ambsf(plotind,1),'or')
		for j=1:length(plotind)
			line([x(plotind(j)) x(plotind(j))],[ambsf(plotind(j),2) ambsf(plotind(j),3)],'color','r','linestyle','-')
		end
		xlabel('Depth (mbsf)')
		ylabel('Adjusted Depth (ambsf)')
		set(gca,'tickdir','out','box','on','xlim',[mbsf_top(ind) mbsf_bot(ind)])
		title({['M0063' upper(hole) '-' num2str(plotcore(i))] []},'fontsize',12)
		grid on
	end
end
