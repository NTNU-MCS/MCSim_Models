
%Makes struct for Discrete Kalman Filter
samples = 10;    %10Hz
T_sample = 1/samples;
omega_0=0.7823;
T=10;
lambda=0.1;
m=vessel.MRB(6,6);


A=[0 1 0 0 0;-(omega_0^2) -(2*lambda*omega_0) 0 0 0;0 0 0 1 0;0 0 0 -(1/T) 1;0 0 0 0 0];
B=[0;0;0;1/m;0];
E=[0 0;1 0;0 0;0 0;0 1];
C=[0 1 1 0 0];

[Adisc,Bdisc,Cdisc,Ddisc] = c2dm(A,B,C,1,T_sample,'zoh');
D=[1 1];
[Adisc,Edisc,Cdisc,Ddisc] = c2dm(A,E,C,D,T_sample,'zoh');

R = 6.0954*10^(-7)/(T_sample);
Q=[30 0;0 1e-6];

P0_=[1 0 0 0 0;0 1 0 0 0;0 0 1 0 0;0 0 0 1 0;0 0 0 0 1];
xh_=[0 0 0 0 0]';

%struct_kalman=struct('Adisc',Adisc,'B_d',Bdisc,'Cdisc',Cdisc,'Edisc',Edisc,'Q',Q,'R',R,'P0_apri',P0_,'x0_hat_apri',x0hat_);

