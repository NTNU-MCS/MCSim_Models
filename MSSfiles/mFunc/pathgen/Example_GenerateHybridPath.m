% Function that calculates the hybrid continuous path based on waypoints.

%% Initialization
% Waypoints:
% X_wp = [0 5 15];
% Y_wp = [0 9 10];
% X_wp = [0 5 10 15 20 25 30 35 40];
% Y_wp = [0 9 10 11 20 10 20 10 20];
X_wp = [10  5  0  5 10  5  0  5 10  5  0  5 10];
Y_wp = [10 15 10  5 10 15 10  5 10 15 10  5 10];

s           = 0.5;  % Path parameter in [0, NumSubPaths]
lambda      = 0.3; % Curvature constant.
r           = 3;    % Differentiability order.
PlotHandle  = 1;    % 0: No plotting, 1: Plot the path only, 2: Plot all derivative curves

%% Run scripts
WP = [X_wp' Y_wp'];
Path = GenerateHybridPath(WP,r,lambda,PlotHandle);
PathSignals = GetHybridPathSignals(Path,s);


%% Extra plotting
if PlotHandle > 1
    figure(PlotHandle); hold on;
    plot(PathSignals.pd(1),PathSignals.pd(2),'sm','LineWidth',1.5);
    plot([PathSignals.pd(1) PathSignals.pd(1)+PathSignals.pd_der{1}(1)],[PathSignals.pd(2) PathSignals.pd(2)+PathSignals.pd_der{1}(2)],'m','LineWidth',1.25);

    s = [0:.01:Path.NumSubpaths];
    for k=1:Path.Order
        pd_der = zeros(length(s),2);
        for j=1:length(s)
            PS = GetHybridPathSignals(Path,s(j));
            pd_der(j,:) = PS.pd_der{k};
        end
    
        % Plotting derivatives
        figure(PlotHandle+k); clf; 
        plot(s,pd_der(:,1),s,pd_der(:,2),'LineWidth',1.25); 
        text_x = ['x_d^{s^',num2str(k),'}'];
        text_y = ['y_d^{s^',num2str(k),'}'];
        legend(text_x,text_y); grid on;
    end
end

