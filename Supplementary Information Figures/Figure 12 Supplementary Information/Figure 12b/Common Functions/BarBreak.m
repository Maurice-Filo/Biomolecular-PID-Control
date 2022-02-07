function Handle_Bar = BarBreak(Handle_Axis, X, Y, Y1, Y2, Type, Scale, BarWidth)
Y_Up = Y;
Y_Up(Y_Up>=Y2) = Y_Up(Y_Up>=Y2) - Y_Up(Y_Up>=Y2) * Scale;

%find the max and min and cut max to 1.5 times the min
Handle_Bar = bar(X, Y_Up, BarWidth);
XLim = Handle_Axis.XLim;
XTick = Handle_Axis.XTick;
YTick = Handle_Axis.YTick;
[~,i] = min(YTick<=Y1);
y = (YTick(i) + YTick(i-1)) / 2;
dy = (YTick(2) - YTick(1)) / 10;
x = XTick(1);
dx = (XTick(2) - XTick(1)) / 2;
switch Type
    case 'Patch'
		% this can be vectorized
        dx = (XLim(2) - XLim(1)) / 10;
        yy = repmat([y - 2*dy, y - dy], 1, 6);
        xx = XLim(1) + dx*(0:11);
		patch([xx(:); flipud(xx(:))], ...
            [yy(:); flipud(yy(:) - 2*dy)], ...
            [.8 .8 .8])
    case 'RPatch'
		% this can be vectorized
        dx=(XLim(2)-XLim(1))./100;
        yy=y+rand(101,1).*2.*dy;
        xx=XLim(1)+dx.*(0:100);
		patch([xx(:);flipud(xx(:))], ...
            [yy(:);flipud(yy(:)-2.*dy)], ...
            [.8 .8 .8])
    case 'Line'
		line([x-dx x   ],[y-2.*dy y-dy   ]);
		line([x    x+dx],[y+dy    y+2.*dy]);
		line([x-dx x   ],[y-3.*dy y-2.*dy]);
		line([x    x+dx],[y+2.*dy y+3.*dy]);
end

%ytick(ytick>y_break_start)=ytick(ytick>y_break_start)+y_break_mid;

YTick(YTick>Y1)=YTick(YTick>Y1)+(Y(Y>=Y2)*Scale);

for i=1:length(YTick)
   yticklabel{i}=sprintf('%d',YTick(i));
end
set(gca,'yticklabel',yticklabel);