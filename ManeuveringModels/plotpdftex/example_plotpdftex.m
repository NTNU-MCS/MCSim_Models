close all
clear all


x = 0:0.01:1;
y_1 = sin(2*pi*x);
y_2 = cos(4*pi*x);

F1 = figure;
plot(x,y_1,'r--')
hold on
plot(x,y_2,'b-')
legend('$\sin(2\pi x)$','$\cos(4\pi x)$','Location','Best')
xlabel('$x$ [\si{\metre}]')
ylabel('$y$ [\si{\metre}]')

legend2latex(F1);   
plotpdftex(F1,'filename',[1 1],'siunitx')