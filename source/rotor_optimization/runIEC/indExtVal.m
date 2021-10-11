function [i_peak,y_peak] = indExtVal(y)
th = mean(y)+1.4.*std(y);
zerIdx=[];
y_peak = [];
i_peak = [];

if y(1)>th
    y_peak(1) = y(1);
    i_peak(1) = 1;
end

zerIdy=[];
for i=1:length(y)-1
    if (y(i)<=th && y(i+1)>=th)
        zerIdy(end+1)=i; % save index of threshold crossing
    end
end

for i=1:length(zerIdy)-1
    [y_peak(i),i_peak(i)] =  max(y(zerIdy(i):zerIdy(i+1)));
    i_peak(i) =  i_peak(i) + zerIdy(i) -1;
end

[y_peak(end+1),i_peak(end+1)] = max(y(zerIdy(end):end));
 i_peak(end) =  i_peak(end) + zerIdy(end) -1;
 
 
  plot(y), hold on
plot(th.*ones(1,length(y)))
plot(i_peak,y_peak,'ro')
 
 
end