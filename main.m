clear all
clc
dbstop if error
warning off

%choose the track 
map = 'Shanghai';% Shanghai or Silverstone
str = [map '.csv'];
track = readtable(str);
data = table2array(track);

importTrackToScenarioDesigner('Shanghai.csv')


% pick = 1 => shortest path, 2 => minimum-time approximation
pick = 1;

%Frenet coordinate（SL coordinate）
x = data(:,1);% Extract raw data
y = data(:,2);
twr = data(:,3);
twl = data(:,4);

% Number of interpolation points
nseg = 1500;
pathXY = [x y];

% Compute segment distances
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
dL = hypot(dx,dy);% same as dL = sqrt(dx.^2+dy.^2)

% Generate offset mapping function using the normal vector
xoff = @(a) -a*dy./dL + xt;
yoff = @(a)  a*dx./dL + yt;


% Calculate inner and outer edges
offset = [-twrt twlt];
for i = 1:nseg
    % Inner edge
    xinner = xoff(offset(i,1));      
    yinner = yoff(offset(i,1));
    % outer edge
    xouter  = xoff(offset(i,2));      
    youter  = yoff(offset(i,2));
end

% Delta between inner and outer edges
delx = xouter - xinner;
dely = youter - yinner;

trackData = [xt yt xinner yinner xouter youter];

% Construct H and B matrices for the quadratic program
H1 = zeros(nseg);
B1 = zeros(size(delx)).';

% --- Shortest Path Formulation ---
for i=1:nseg-1
    H1(i,i)     = H1(i,i)     + delx(i)^2          + dely(i)^2;
    H1(i+1,i)   = H1(i+1,i)   - delx(i)*delx(i+1)  - dely(i)*dely(i+1);
    H1(i,i+1)   = H1(i,i+1)   - delx(i)*delx(i+1)  - dely(i)*dely(i+1);
    H1(i+1,i+1) = H1(i+1,i+1) + delx(i+1)^2        + dely(i+1)^2;
end

for i=1:nseg-1
    B1(1,i)   = B1(1,i)   - 2*(xinner(i+1)-xinner(i))*delx(i)   - 2*(yinner(i+1)-yinner(i))*dely(i);
    B1(1,i+1) = B1(1,i+1) + 2*(xinner(i+1)-xinner(i))*delx(i+1) + 2*(yinner(i+1)-yinner(i))*dely(i+1);
end

H2 = zeros(nseg);
B2 = zeros(size(delx)).';

% --- Minimum-Curvature Formulation ---
for i=2:nseg-1
    
    % first row
    H2(i-1,i-1) = H2(i-1,i-1) + delx(i-1)^2         + dely(i-1)^2;
    H2(i-1,i)   = H2(i-1,i)   - 2*delx(i-1)*delx(i) - 2*dely(i-1)*dely(i);
    H2(i-1,i+1) = H2(i-1,i+1) + delx(i-1)*delx(i+1) + dely(i-1)*dely(i+1);
    
    %second row
    H2(i,i-1)   = H2(i,i-1)   - 2*delx(i-1)*delx(i) - 2*dely(i-1)*dely(i);
    H2(i,i)     = H2(i,i )    + 4*delx(i)^2         + 4*dely(i)^2;
    H2(i,i+1)   = H2(i,i+1)   - 2*delx(i)*delx(i+1) - 2*dely(i)*dely(i+1);
    
    % third row
    H2(i+1,i-1) = H2(i+1,i-1) + delx(i-1)*delx(i+1) + dely(i-1)*dely(i+1);
    H2(i+1,i)   = H2(i+1,i)   - 2*delx(i)*delx(i+1) - 2*dely(i)*dely(i+1);
    H2(i+1,i+1) = H2(i+1,i+1) + delx(i+1)^2         + dely(i+1)^2;
    
end

% formation of B matrix (1xn)
for i=2:nseg-1
    B2(1,i-1) = B2(1,i-1) + 2*(xinner(i+1)+xinner(i-1)-2*xinner(i))*delx(i-1) + 2*(yinner(i+1)+yinner(i-1)-2*yinner(i))*dely(i-1);
    B2(1,i)   = B2(1,i)   - 4*(xinner(i+1)+xinner(i-1)-2*xinner(i))*delx(i)   - 4*(yinner(i+1)+yinner(i-1)-2*yinner(i))*dely(i);
    B2(1,i+1) = B2(1,i+1) + 2*(xinner(i+1)+xinner(i-1)-2*xinner(i))*delx(i+1) + 2*(yinner(i+1)+yinner(i-1)-2*yinner(i))*dely(i+1);
    
end

% Boundaries
lb = zeros(nseg,1);
ub = ones(nseg,1);

% Equality constraint: offset(1) = offset(end)
Aeq      =   zeros(1,nseg);
Aeq(1)   =   1;
Aeq(end) =   -1;
beq      =   0;
    
% QP options
options = optimoptions('quadprog','Display','iter');

% Choose which objective to solve
if pick==1
    [resSP,fval,exitflag,output] = quadprog(2*H1,B1,[],[],Aeq,beq,lb,ub,[],options);
elseif pick==2
    [resSP,fval,exitflag,output] = quadprog(2*H2,B2,[],[],Aeq,beq,lb,ub,[],options);
end

% Visualization
% Map the solution offsets back to real coordinates
xresSP = zeros(size(xt));
yresSP = zeros(size(xt));
for i = 1:numel(xt)
    xresSP(i) = xinner(i)+resSP(i)*delx(i);
    yresSP(i) = yinner(i)+resSP(i)*dely(i);
end

trajSP = [xresSP yresSP];

% Simple kinematic check
amin = -1;
amax = 1;
vmax = 10;
vmin = 0;
nseg = nseg-1;
dist = zeros(nseg,1);

% Distance for each pair of adjacent points
for t = 1:nseg
    dist(t) = sqrt((yresSP(t+1)-yresSP(t))^2+(xresSP(t+1)-xresSP(t))^2);
end

% Total distance
dist_total = sum(dist);

% Estimate curvature
dx = gradient(xresSP);%speed
d2x = gradient(dx);%acceleration
dy = gradient(yresSP);
d2y = gradient(dy);
k = sum(abs(dx.*d2y-dy.*d2x)./((dx.^2+dy.^2).^(3/2)));%curvature

% Estimate speed profile and time
v_origin = sqrt(dx.^2+dy.^2);
a_origin = gradient(v_origin);
a_origin_min = min(a_origin);
a_origin_max = max(a_origin);

scale_v = vmax/max(v_origin);
scale_a = min(amax/a_origin_max,amin/a_origin_min);
scale = min(scale_a,scale_v);
v = v_origin*scale;
t_toal = sum(dist./v(2:end));

% Plot final path
figure()
plot(xresSP,yresSP,'color','r','linew',2)
hold on
plot(xt,yt,'--')
plot(xinner,yinner,'color','g')
plot(xouter,youter,'color','k')
hold off
axis equal
xlabel('x(m)')
ylabel('y(m)')
if pick==1
    legend('shortest path','Reference Line','Inner Track','Outer Track',"FontSize",24)
elseif pick==2
    legend('fastest path','Reference Line','Inner Track','Outer Track',"FontSize",24)
end
title([sprintf(map) ' ' 'distance:' num2str(dist_total) ' ' 'cost time:' num2str(t_toal)],"FontSize",24)

% Plot speed heat map
figure()
scatter(xresSP,yresSP,5,v,'filled')
colorbar
title([sprintf(map) ' ' 'speed-heatmap'],"FontSize",24)

% Build a driving scenario for track visualization
%s = drivingScenario;
%roadcenter = [x y];
% roadcenter = fliplr(roadcenter);
%road(s,roadcenter,twr+twl,'Lanes',lanespec([4,4]));

% figure()
%plot(s)
%%%%% Weighted Sum of Shortest Path & Min Curvature %%%%%
% We'll create an extra figure that shows how combining
% H1,B1 (shortest path) and H2,B2 (min curvature)
% yields a "trade-off" between total distance & curvature.

% 1) Define weight values from 0 to 1
w_values = 0:0.1:1;

% 2) Prepare arrays to store results
numW = length(w_values);
distVals = zeros(numW,1);
curvVals = zeros(numW,1);

% We'll use the original 'nC' = numel(xt) for coordinate arrays
nC = numel(xt);

for iw = 1:numW
    w = w_values(iw);
    
    % 3) Combine the QP matrices 
    Hcomb = w*H1 + (1-w)*H2;
    Bcomb = w*B1 + (1-w)*B2;
    
    % Solve using quadprog (we multiply H by 2 for MATLAB's QP form)
    [alphaCmb, ~, exitflag] = quadprog(2*Hcomb, Bcomb, ...
                                       [],[], Aeq, beq, lb, ub, [], options);
    if exitflag <= 0
        warning('Quadprog did not converge at w=%.2f', w);
    end

    % 4) Map alphaCmb -> (xCmb, yCmb)
    xCmb = zeros(nC,1);
    yCmb = zeros(nC,1);
    for i = 1:nC
        xCmb(i) = xinner(i) + alphaCmb(i)*delx(i);
        yCmb(i) = yinner(i) + alphaCmb(i)*dely(i);
    end

    % 5) Compute total distance for this combined solution
    distVal = 0;
    for t = 1:(nC-1)
        distVal = distVal + sqrt( (xCmb(t+1)-xCmb(t))^2 + (yCmb(t+1)-yCmb(t))^2 );
    end
    distVals(iw) = distVal;

    % 6) Compute total curvature (similar to your "Estimate curvature" code)
    dxC = gradient(xCmb);
    dyC = gradient(yCmb);
    d2xC = gradient(dxC);
    d2yC = gradient(dyC);
    kC = sum( abs(dxC.*d2yC - dyC.*d2xC) ./ ((dxC.^2 + dyC.^2).^(3/2)) , 'omitnan');
    curvVals(iw) = kC;
end

% 7) Plot the trade-off curve
figure('Name','Distance vs. Curvature (Weighted Sum)','Color','white');
plot(distVals, curvVals, 'o-','LineWidth',1.5); 
xlabel('Total Distance',"FontSize",24); 
ylabel('Total Curvature',"FontSize",24);
title('Trade-off: Combining Shortest Path & Min Curvature',"FontSize",24);
grid on;

% Optionally label each point with w
for iw=1:numW
    text(distVals(iw), curvVals(iw), ...
         sprintf('w=%.1f', w_values(iw)), ...
         'HorizontalAlignment','left','FontSize',8);
end
