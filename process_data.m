clear all
clc
dbstop if error
warning off
%choose the track 
map = 'Shanghai';%Shanghai Silverstone
str = [map '.csv'];
track = readtable(str);
data = table2array(track);

%Frenet coordinate（SL coordinate） 
x = data(:,1);
y = data(:,2);
twr = data(:,3);
twl = data(:,4);

% Number of interpolation points
nseg = 1500;
pathXY = [x y];
stepLengths = sqrt(sum(diff(pathXY,[],1).^2,2));
stepLengths = [0; stepLengths]; % Add 0 at the start
cumulativeLen = cumsum(stepLengths);% Cumulative distance along the track

% Interpolation
finalStepLocs = linspace(0,cumulativeLen(end), nseg);
finalPathXY = interp1(cumulativeLen, pathXY, finalStepLocs);
xt = finalPathXY(:,1);
yt = finalPathXY(:,2);
twrt = interp1(cumulativeLen, twr, finalStepLocs,'spline')';
twlt = interp1(cumulativeLen, twl, finalStepLocs,'spline')';

% Gradients along center line
dx = gradient(xt);
dy = gradient(yt);
dL = hypot(dx,dy);%dL = sqrt(dx.^2+dy.^2)

% Generate offset mapping function using the normal vector
xoff = @(a) -a*dy./dL + xt;
yoff = @(a)  a*dx./dL + yt;



% Calculate inner and outer edges
offset = [-twrt twlt];
for i = 1:nseg
    %inner
    xinner = xoff(offset(i,1));      
    yinner = yoff(offset(i,1));
    %outer
    xouter  = xoff(offset(i,2));      
    youter  = yoff(offset(i,2));
end

% Delta between inner and outer edges
delx = xouter - xinner;
dely = youter - yinner;

trackData = [xt yt xinner yinner xouter youter];
amin = -1;
amax = 1;
vmax = 10;
vmin = 0;
save('data.mat','trackData','nseg','vmin','vmax','amax','amin','map')