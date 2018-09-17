%function alp=alph(x,r,p,k,plot_flag)
plot_flag = 1;
x=0:.01:5;

for count=1:length(x)
    b(count)=exp(-.08*(.45*(x(count)))^15);
end

if plot_flag == 1
    figure
    plot(x,b,'b')
    grid on
end
