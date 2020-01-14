% Function that tests the online continuous path based on waypoints that are produced on-the-fly.

%% Initialization
% Waypoints that we will test flying through:
% X_wp = [0 20 40];
% Y_wp = [0 15 30];
X_wp = [0 5 15 20];
Y_wp = [0 9 10 15];
% X_wp = [0 5 10 15 20 25 30 35 40];
% Y_wp = [0 9 10 11 20 10 20 10 20];
% X_wp = [10  5  0  5 10  5  0  5 10  5  0  5 10];
% Y_wp = [10 15 10  5 10 15 10  5 10 15 10  5 10];
% X_wp = [10   40  80 120 120 80 40 40 80 120 120 80 40 10  8  8 10  40  80];
% Y_wp = [98  100 100  98  72 70 68 42 40  38  12 10 10 12 40 70 98 100 100];

s           = 1.2;  % Path parameter in [0, NumSubPaths]
lambda      = 2; % Tuning constant.
r           = 3;    % Differentiability order.
PlotHandle  = 2;    % 1: Plot the path only, 2: Plot all derivative curves

%% Run scripts
WP = [X_wp' Y_wp'];
N  = length(WP)-1;

t0 = WP(2,:)-WP(1,:);
Path = {};
figure(1); clf;

for j=1:N
    p0 = WP(j,:);
    pt = WP(j+1,:);
    Path{j} = GenerateOnlinePathSegment(p0,pt,t0,r,lambda,1);
    t0 = WP(j+1,:)-WP(j,:);
    
    figure(1); hold on;
end
axis equal; axis([min(WP(:,2))-1.0 max(WP(:,2))+1.0 min(WP(:,1))-1.0 max(WP(:,1))+1.0]); 

ii    = floor(s) + 1;
theta = s - ii + 1;
PathSignals = GetHybridPathSignals(Path{ii},theta);

%% Extra plotting
if PlotHandle > 1
    figure(PlotHandle); hold on;
    plot(PathSignals.pd(2),PathSignals.pd(1),'sm','LineWidth',1.5);
    plot([PathSignals.pd(2) PathSignals.pd(2)+PathSignals.pd_der{1}(2)],[PathSignals.pd(1) PathSignals.pd(1)+PathSignals.pd_der{1}(1)],'m','LineWidth',1.25);

    s = [0:.01:N-0.01];
    for k=1:Path{1}.Order
        pd_der = zeros(length(s),2);
        for j=1:length(s)
            ii    = floor(s(j)) + 1;
            theta = s(j) - ii + 1
            PS = GetHybridPathSignals(Path{ii},theta);
            pd_der(j,:) = PS.pd_der{k};
        end

        % Plotting derivatives
        figure(PlotHandle+k-1); clf; 
        plot(s,pd_der(:,1),s,pd_der(:,2),'LineWidth',1.25); 
        text_x = ['x_d^{s^',num2str(k),'}'];
        text_y = ['y_d^{s^',num2str(k),'}'];
        legend(text_x,text_y); grid on;
    end
end

