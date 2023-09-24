rmp=@(t,a,b,t1,t2) a+(b-a)/(t2-t1)*(t-t1)*double(t>t1)-(b-a)/(t2-t1)*(t-t2)*double(t>t2);
Xm=@(t) rmp(t,-7,-0.1,3,5);
Q=@(t)  rmp(t,0,100,1,2);
color=@(t) rmp(t,1,0.01,22.5,25);
 close all
% for t=0:0.1:5
%     plot(t,Xm(t))
%     hold on
% end
t=5
a=-3
b=2
t1=1
t2=4
fplot(color,[22,26])