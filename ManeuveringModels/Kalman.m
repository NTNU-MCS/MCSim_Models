function [sys,x0,str,ts] = Kalman(t,x,u,struct_kalman)


    size = simsizes; % do not modify

    size.NumContStates  = 0; % Number of continuous states in the system, do not modify
    size.NumDiscStates  = 35; % Number of discrete states in the system, modify. 
    size.NumOutputs     = 1; % Number of outputs, the hint states 2
    size.NumInputs      = 2; % Number of inputs, the hint states 2
    size.DirFeedthrough = 1; % 1 if the input is needed directly in the
    % update part
    size.NumSampleTimes = 1; % Do not modify  

    sys = simsizes(size); % Do not modify  

    x0  = [struct_kalman.x0_hat_apri', 0 0 0 0 0, mat2vec(struct_kalman.P0_apri)']'; % Initial values for the discrete states.

    str = []; % Do not modify

    ts  = [-1 0]; 
    
    P_= vec2mat(x(11:35));                                                      %Extract P from x
    K = P_*struct_kalman.C_d'/(struct_kalman.C_d*P_*struct_kalman.C_d'+struct_kalman.R);            %New Kalman gain
    xh_= x(1:5);                                                                %Get Apriori state estimate
    xh = xh_ + K*(u(1)-struct_kalman.C_d*xh_);                                       %Calculate state estimate
    P=(eye(5)-K*struct_kalman.C_d)*P_*(eye(5)-K*struct_kalman.C_d)'+K*struct_kalman.R*K';      %Error covariance qe (4.2.11) Brown&Hwang
    xh_=struct_kalman.A_d*xh + struct_kalman.B_d*u(2);                                    %New Apriori state estimate
    P_= struct_kalman.A_d*P*struct_kalman.A_d' + struct_kalman.E_d*struct_kalman.Q*struct_kalman.E_d';   %New Apriori covariance

    sys=[xh_', xh', mat2vec(P_)']';

    sys=[x(8)];