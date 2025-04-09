function importTrackToScenarioDesigner(csvFile)
    % IMPORTTRACKTOSCENARIODESIGNER  Create a drivingScenario from a CSV track file
    %
    % csvFile: A string specifying the CSV file path
    %
    % CSV format assumption:
    %    Col1: x_m        (track center X)
    %    Col2: y_m        (track center Y)
    %    Col3: w_tr_right (right boundary distance)
    %    Col4: w_tr_left  (left boundary distance)

    % 1) Read CSV
    data = readmatrix(csvFile);
    x   = data(:,1);
    y   = data(:,2);
    wr  = data(:,3);  % right offset
    wl  = data(:,4);  % left offset

    % 2) Create drivingScenario object
    scenario = drivingScenario;

    % 3) Build the road from track center points
    roadCenters = [x, y];
    totalWidth  = wr + wl;  
    % Lanespec(2) => 2-lane road (1 lane each direction), adjust if needed
    road(scenario, roadCenters, totalWidth, 'Lanes', lanespec(2));

    % 4) Directly open Driving Scenario Designer with the newly created scenario
    drivingScenarioDesigner(scenario);

    % figure; plot(scenario);
end
