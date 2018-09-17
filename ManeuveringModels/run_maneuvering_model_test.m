clear all;
close all;
clc;

%% Collecting Data For R/V Gunnerus Vessel Model
InitGunnerusMan;

%% Run Simulation
model='maneuvering';
load_system(model);
simulation_time = 200;
%Running the simulation
tic;
warning off
simOut=sim(model,'SaveOutput','on',...
        'SaveFormat','StructureWithTime',...
        'ReturnWorkspaceOutputs','on',...
        'SaveOutput','on','OutputSaveName','res_gun',...
        'StartTime','0','StopTime',num2str(Tsim) );
    
%% Extracting Values from Simulation
%True Values
CSEI_eta = simOut.find('eta_true');
CSEI_nu=simOut.find('nu_true');

RVG_eta = simOut.find('eta_true_gunnerus');
RVG_nu=simOut.find('nu_true_gunnerus');

%Unscaled Models
nomoto=simOut.find('nomoto');
nomoto2nd=simOut.find('nomoto_2nd');
sway_yaw=simOut.find('sway_yaw');
surge_speed=simOut.find('surge_speed');
inputs=simOut.find('inputs');

%Scaled Models
nomoto_scaled=simOut.find('nomoto_scaled');

surge_speed_scaled=simOut.find('surge_speed_scaled');
scaled_inputs=simOut.find('scaled_inputs');


%% Plots

%Maneuvering Models for Model Vessel
figure(1)
plot(CSEI_nu.time(:,1),CSEI_nu.signals.values(:,1),'o','MarkerSize',2);
hold on;
plot(surge_speed.time(:,1),surge_speed.signals.values(:,2),'-','linewidth',1);
legend('u_{true}','u_{spm}')
ylabel('u[m/s]')
xlabel('Time[s]')
grid on

figure(2)
subplot(2,2,1)
plot(CSEI_eta.time(:,1),CSEI_eta.signals.values(:,2),'o','MarkerSize',2);
hold on;
plot(sway_yaw.time(:,1),sway_yaw.signals.values(:,3),'-','linewidth',1);
legend('y_{true}','y_{s-y}')
ylabel('y[m]')
xlabel('Time[s]')
grid on

subplot(2,2,2)
plot(CSEI_eta.time(:,1),CSEI_eta.signals.values(:,3),'o','MarkerSize',2);
hold on;
plot(sway_yaw.time(:,1),sway_yaw.signals.values(:,1),'-','linewidth',1);
legend('\psi_{true}','\psi_{s-y}')
ylabel('\psi[rad]')
xlabel('Time[s]')
grid on

subplot(2,2,3)
plot(CSEI_nu.time(:,1),CSEI_nu.signals.values(:,2),'o','MarkerSize',2);
hold on;
plot(sway_yaw.time(:,1),sway_yaw.signals.values(:,4),'-','linewidth',1);
legend('v_{true}','v_{s-y}')
ylabel('v[m/s]')
xlabel('Time[s]')
grid on

subplot(2,2,4)
plot(CSEI_nu.time(:,1),CSEI_nu.signals.values(:,3),'o','MarkerSize',2);
hold on;
plot(sway_yaw.time(:,1),sway_yaw.signals.values(:,2),'-','linewidth',1);
legend('r_{true}','r_{s-y}')
ylabel('r[rad/s]')
xlabel('Time[s]')
grid on 

figure(3)
subplot(2,1,1)
plot(CSEI_eta.time(:,1),CSEI_eta.signals.values(:,3),'o','MarkerSize',2);
hold on;
plot(CSEI_nu.time(:,1),CSEI_nu.signals.values(:,3),'o','MarkerSize',2);
hold on;
plot(nomoto.time(:,1),nomoto.signals.values(:,1),'-','linewidth',1);
hold on;
plot(nomoto.time(:,1),nomoto.signals.values(:,2),'-','linewidth',1);
legend('\psi_{true}','r_{true}','\psi_{1st}','r_{1st}')
ylabel('\psi[rad],r[rad/s]')
xlabel('Time[s]')
grid on

subplot(2,1,2)
plot(CSEI_eta.time(:,1),CSEI_eta.signals.values(:,3),'o','MarkerSize',2);
hold on;
plot(CSEI_nu.time(:,1),CSEI_nu.signals.values(:,3),'o','MarkerSize',2);
hold on;
plot(nomoto2nd.time(:,1),nomoto2nd.signals.values(:,1),'-','linewidth',1);
hold on;
plot(nomoto2nd.time(:,1),nomoto2nd.signals.values(:,2),'-','linewidth',1);
legend('\psi_{true}','r_{true}','\psi_{2nd}','r_{2nd}')
ylabel('\psi[rad],r[rad/s]')
xlabel('Time[s]')
grid on

figure(4)
plot(inputs.time(:,1),inputs.signals.values,'-','linewidth',1);
axis([0 200 -1 2])
legend('\tau_{X}','\tau_{Y}','\tau_{N}')
ylabel('\tau_X[N],\tau_Y[N],\tau_N[Nm]')
xlabel('Time[s]')
grid on

%Plots for Scaled Models

figure(5)
plot(RVG_nu.time(:,1),RVG_nu.signals.values(:,1),'o','MarkerSize',2);
hold on;
plot(surge_speed_scaled.time(:,1),surge_speed_scaled.signals.values(:,2),'-','linewidth',1);
legend('u_{true}','u_{spm}')
ylabel('u[m/s]')
xlabel('Time[s]')
legend('u_{true}','u_{spm}')
ylabel('u[m/s]')
xlabel('Time[s]')
grid on

figure(6)
plot(RVG_eta.time(:,1),RVG_eta.signals.values(:,6),'o','MarkerSize',2);
hold on;
plot(RVG_nu.time(:,1),RVG_nu.signals.values(:,6),'o','MarkerSize',2);
hold on;
plot(nomoto_scaled.time(:,1),nomoto_scaled.signals.values(:,1),'-','linewidth',1);
hold on;
plot(nomoto_scaled.time(:,1),nomoto_scaled.signals.values(:,2),'-','linewidth',1);
legend('\psi_{true}','r_{true}','\psi_{1st}','r_{1st}')
ylabel('\psi[rad/s]')
xlabel('Time[s]')
grid on

figure(7)
plot(scaled_inputs.time(:,1),scaled_inputs.signals.values(:,1),'-','linewidth',1);
hold on;
plot(scaled_inputs.time(:,1),scaled_inputs.signals.values(:,2),'-','linewidth',1);
hold on;
plot(scaled_inputs.time(:,1),scaled_inputs.signals.values(:,3),'-','linewidth',1);
legend('\tau_{X}[N]','\tau_{Y}[N]','\tau_{N}[Nm]')
ylabel('Scaled Input')
xlabel('Time[s]')
grid on