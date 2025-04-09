clear all  
clc
dbstop if error
load('data.mat')% Load pre-processed track data

% pick = 1 => shortest path, 2 => minimum-time approximation
pick = 2;

% Lower and upper bounds for the decision variable alpha
lb = 0;
ub = 1;

% Dimension of the decision space = number of segments
dim = nseg;

% Grey Wolf Optimizer parameters
pop_num=100;  % Number of search agents 
Max_iter=20000;    % Maximum number of iterations 

% Objective function handle
if pick==1
    fun = @(x) evaluate(x);
else
    fun = @(x) evaluate2(x);
end

% Call the GWO optimizer
[Alpha_score,Alpha_pos,Convergence_curve]=GWO(pop_num,Max_iter,lb,ub,dim,fun);

% Plot the convergence     
figure()
plot(Convergence_curve)
hold on 
hold off
xlabel('Iterations')
if pick==1
    ylabel('Total Distance(m)')
else
    ylabel('Time Cost(s)')
end
title('GWO')

% Convert the best alpha solution into actual x,y coordinates
% trackData(:,3), trackData(:,4) => xInner, yInner
% trackData(:,5), trackData(:,6) => xOuter, yOuter
Alpha_pos = Alpha_pos';
Alpha_pos = [(trackData(:,5)-trackData(:,3)).*Alpha_pos+trackData(:,3),(trackData(:,6)-trackData(:,4)).*Alpha_pos+trackData(:,4)];

% Close the loop by appending the first point at the end
Alpha_pos = [Alpha_pos;Alpha_pos(1,:)];

% Compute total distance
dist = zeros(nseg,1);
for i = 1:nseg
    dist(i) = sqrt(sum((Alpha_pos(i,:)-Alpha_pos(i+1,:)).^2));
end
dist_total = sum(dist);


% Extract final X and Y results
xresSP = Alpha_pos(:,1);
yresSP = Alpha_pos(:,2);

% Numerical gradients for velocity/acceleration estimates
dx = gradient(xresSP);%speed
d2x = gradient(dx);%acceleration
dy = gradient(yresSP);
d2y = gradient(dy);
% Calculate curvature (discrete approximation)
k = sum(abs(dx.*d2y-dy.*d2x)./((dx.^2+dy.^2).^(3/2)));

% Approximate velocity profile and time
v_origin = sqrt(dx.^2+dy.^2);
a_origin = gradient(v_origin);
a_origin_min = min(a_origin);
a_origin_max = max(a_origin);


% Scale velocity due to max speed (vmax) and max acceleration (amin, amax)
scale_v = vmax/max(v_origin);
scale_a = min(amax/a_origin_max,amin/a_origin_min);
scale = min(scale_a,scale_v);
v = v_origin*scale;
t_toal = sum(dist./v(2:end));

% Plot speed heatmap over the resulting path
figure()
scatter(xresSP,yresSP,5,v,'filled')
colorbar
title([sprintf(map) ' ' 'GWO' ' ' 'speed-heatmap'])

% Finally, plot the resulting path
figure()
plot(Alpha_pos(:,1),Alpha_pos(:,2),'color','r','linew',2)
if pick==1
    legend('shortest path')
elseif pick==2
    legend('fastest path')
end
title([sprintf(map) ' ' 'distance:' num2str(dist_total) ' ' 'cost time:' num2str(t_toal)])