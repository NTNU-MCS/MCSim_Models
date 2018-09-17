clear all;
close all;
clc;

set(0,'defaulttextinterpreter','none');
x = [0:0.1:(2*pi)];
fig1 = figure(1);
plot(x,sin(x),x,cos(x));
legend('sin($\theta$)','cos($\theta$)');
xlabel('$0 \leq \theta \leq 2\pi$')
legend2latex(fig1);
laprint(1,'figure1');