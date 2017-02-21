function y=normal(x)

%checks to see if the Statistics Toolbox is installed.
%If so, then nanmean and nanstd are used.
%If not, mean and std are used.
%Do not pass NaN values if you do not have the Statistics Toolbox installed

v=ver;
if isempty(strfind([v.Name],'Statistics'))==0
	if length(x)>1
		ms=[nanmean(x) nanstd(x)];
		y=(x-ms(1))./ms(2);
	else
		y=x;
	end
else
	if length(x)>1
		ms=[mean(x) std(x)];
		y=(x-ms(1))./ms(2);
	else
		y=x;
	end
end
