function combinedCostDemo()
% COMBINEDCOSTDEMO  Demonstrate combining "shortest path" and "min curvature"
% objectives in one script. It uses a weighted sum of the two QP forms 
% (H1,B1) and (H2,B2). 
%
% Steps:
%  1) Build track geometry + delx, dely, etc.
%  2) Construct H1,B1 for shortest path and H2,B2 for min curvature
%  3) For w in [0,1], form Hcomb,Bcomb => quadprog => alpha
%  4) Compute distance, curvature => plot trade-off

    %% 1) Load / Generate track data & interpolation
    % ------------------------------------------------
    % TODO: Insert your own code that:
    %  - Reads the CSV (map = 'Silverstone', etc.)
    %  - Interpolates nseg points -> xt, yt
    %  - Calculates xinner, yinner, xouter, youter
    %  - Defines delx = xouter - xinner, dely = youter - yinner
    %  - We assume you have nseg, trackData, etc. at the end.
    
    map = 'Silverstone';
    str = [map '.csv'];
    track = readtable(str);
    data  = table2array(track);
    
    % Example placeholder: do your actual interpolation 
    % Here we just assume you end up with:
    %   xt, yt, xinner, yinner, xouter, youter, delx, dely, nseg
    % ...
    % (Fill in your existing steps from "main" or "process_data".)
    
    % Once done, store trackData if you like
    trackData = [xt, yt, xinner, yinner, xouter, youter];
    
    %% 2) Construct H1,B1 (shortest path) and H2,B2 (min curvature)
    % ------------------------------------------------
    H1 = zeros(nseg);
    B1 = zeros(size(delx)).';
    
    % --- Shortest Path Formulation ---
    for i = 1:nseg-1
        H1(i,i)     = H1(i,i)     + delx(i)^2 + dely(i)^2;
        H1(i+1,i)   = H1(i+1,i)   - delx(i)*delx(i+1) - dely(i)*dely(i+1);
        H1(i,i+1)   = H1(i,i+1)   - delx(i)*delx(i+1) - dely(i)*dely(i+1);
        H1(i+1,i+1) = H1(i+1,i+1) + delx(i+1)^2 + dely(i+1)^2;
    end
    
    for i = 1:nseg-1
        B1(1,i)   = B1(1,i)   - 2*(xinner(i+1)-xinner(i))*delx(i)   - 2*(yinner(i+1)-yinner(i))*dely(i);
        B1(1,i+1) = B1(1,i+1) + 2*(xinner(i+1)-xinner(i))*delx(i+1) + 2*(yinner(i+1)-yinner(i))*dely(i+1);
    end
    
    % --- Min-Curvature Formulation ---
    H2 = zeros(nseg);
    B2 = zeros(size(delx)).';
    
    for i = 2:nseg-1
        % first row
        H2(i-1,i-1) = H2(i-1,i-1) + delx(i-1)^2 + dely(i-1)^2;
        H2(i-1,i)   = H2(i-1,i)   - 2*delx(i-1)*delx(i) - 2*dely(i-1)*dely(i);
        H2(i-1,i+1) = H2(i-1,i+1) + delx(i-1)*delx(i+1) + dely(i-1)*dely(i+1);
        
        % second row
        H2(i,i-1)   = H2(i,i-1)   - 2*delx(i-1)*delx(i) - 2*dely(i-1)*dely(i);
        H2(i,i)     = H2(i,i )    + 4*delx(i)^2 + 4*dely(i)^2;
        H2(i,i+1)   = H2(i,i+1)   - 2*delx(i)*delx(i+1) - 2*dely(i)*dely(i+1);
        
        % third row
        H2(i+1,i-1) = H2(i+1,i-1) + delx(i-1)*delx(i+1) + dely(i-1)*dely(i+1);
        H2(i+1,i)   = H2(i+1,i)   - 2*delx(i)*delx(i+1) - 2*dely(i)*dely(i+1);
        H2(i+1,i+1) = H2(i+1,i+1) + delx(i+1)^2 + dely(i+1)^2;
    end
    
    % formation of B2 (1xn)
    for i = 2:nseg-1
        B2(1,i-1) = B2(1,i-1) + 2*(xinner(i+1)+xinner(i-1)-2*xinner(i))*delx(i-1) ...
                                + 2*(yinner(i+1)+yinner(i-1)-2*yinner(i))*dely(i-1);
        B2(1,i)   = B2(1,i)   - 4*(xinner(i+1)+xinner(i-1)-2*xinner(i))*delx(i) ...
                                - 4*(yinner(i+1)+yinner(i-1)-2*yinner(i))*dely(i);
        B2(1,i+1) = B2(1,i+1) + 2*(xinner(i+1)+xinner(i-1)-2*xinner(i))*delx(i+1) ...
                                + 2*(yinner(i+1)+yinner(i-1)-2*yinner(i))*dely(i+1);
    end
    
    %% 3) Define constraints (lb, ub, Aeq, beq)
    lb = zeros(nseg,1);
    ub = ones(nseg,1);
    Aeq = zeros(1,nseg);
    Aeq(1,1)   = 1;
    Aeq(1,nseg)= -1;
    beq = 0;

    %% 4) Weighted-sum loop for w in [0,1]
    w_values = linspace(0,1,11);
    distArray = zeros(size(w_values));
    curvArray = zeros(size(w_values));
    alphaSolutions = cell(size(w_values));
    
    options = optimoptions('quadprog','Display','none');

    for idx = 1:length(w_values)
        w = w_values(idx);
        
        % Combined cost
        Hcomb = w * H1 + (1 - w)*H2;
        Bcomb = w * B1 + (1 - w)*B2;
        
        % In MATLAB's quadprog: objective = 0.5 * x^T * (2*H) * x + ...
        Hquad = 2 * Hcomb;
        fquad = Bcomb;
        
        [alphaOpt, fval, exitflag] = quadprog(Hquad, fquad, ...
                                              [],[], Aeq, beq, lb, ub, [], options);
        if exitflag <= 0
            warning('Quadprog did not converge at w=%.2f', w);
        end
        
        alphaSolutions{idx} = alphaOpt;
        
        % Compute distance & curvature
        distVal = computeDistance(alphaOpt, xinner, yinner, delx, dely);
        curvVal = computeCurvature(alphaOpt, xinner, yinner, delx, dely);
        
        distArray(idx) = distVal;
        curvArray(idx) = curvVal;
    end

    %% 5) Plot the distance-curvature trade-off
    figure('Name','Distance vs. Curvature','Color','white');
    plot(distArray, curvArray, 'o-','LineWidth',1.5);
    xlabel('Total Distance');
    ylabel('Total Curvature');
    title('Trade-off: Weighted Sum of ShortestPath & MinCurvature');
    grid on;
    
    % Optionally label points
    for i=1:length(w_values)
        text(distArray(i), curvArray(i), ...
            sprintf('w=%.1f', w_values(i)), ...
            'HorizontalAlignment','left','FontSize',8);
    end
    
    %% Show sample best result (if you want to pick by some criterion)
    % [bestScore, bestIdx] = min(distArray + curvArray);
    % fprintf('Best combined at w=%.1f (dist=%.2f, curv=%.2f)\n',...
    %          w_values(bestIdx), distArray(bestIdx), curvArray(bestIdx));
end

%% Example placeholders for computeDistance & computeCurvature
function d = computeDistance(alpha, xinner, yinner, delx, dely)
    % This is a trivial example. You should do your 
    % actual mapping alpha->(x,y) and sum up the pairwise distances.
    
    % map alpha -> x,y
    x = xinner + alpha.*delx;
    y = yinner + alpha.*dely;
    n = length(x);
    dsum = 0;
    for i=1:n-1
        dsum = dsum + sqrt((x(i+1)-x(i))^2 + (y(i+1)-y(i))^2);
    end
    d = dsum;
end

function c = computeCurvature(alpha, xinner, yinner, delx, dely)
    % Another simplified approach: 
    % you might replicate the numeric gradient logic 
    % from your main code to approximate sum of kappa^2.
    x = xinner + alpha.*delx;
    y = yinner + alpha.*dely;
    dx = gradient(x);
    dy = gradient(y);
    d2x = gradient(dx);
    d2y = gradient(dy);
    kArray = abs(dx.*d2y - dy.*d2x)./((dx.^2 + dy.^2).^(3/2));
    c = sum(kArray, 'omitnan');  % sum of curvature (or sum of curvature^2 if you prefer)
end
