%% initialize
clear
close all
format compact
R=3.3;
stmp = char(datetime("now"));
%stmp=stmp(length(stmp)-7:length(stmp));
stmp = strrep(stmp, ':', '_');
txt='my_movie_';
name=[txt stmp];
name=char(name);
fps=60;
T=30;
boundsx=18;
boundsy=12;
simwinx=18;
simwiny=10;
noppr=55;
parspace=0.6;
upplus=8;




% Populate x and y positions with sequential values
%[xPositions, yPositions] = meshgrid(linspace(-boundsx,boundsx,noppr), linspace(-boundsy,boundsy,noppr));
[xPositions, yPositions] = meshgrid(-boundsx:parspace:-boundsx,-boundsy:parspace:boundsy+upplus);
N=numel(xPositions);
particle = struct('x',0,'y',0);

particles(N) = particle;
for i = 1:N
    particles(i).x = xPositions(i);
    particles(i).y = yPositions(i);
   
end
figure('Position', [-3000, +500, 1920, 1080]);
axis equal
xlim([-simwinx,simwinx]);
ylim([-simwiny,simwiny]);
wheel=linspace(0,2*pi(),9);
%% velocity function

eps=0;
uinf=12;
Q=50;
K=uinf*R^2;
omega=@(t) 7*(t-5)/10*double(t>5)-7*(t-15)/10*double(t>15);
B=@(t) 2*pi()*R^2*omega(t);
VB=@(x,y,t) B(t)/2/pi()/sqrt(x^2+y^2+eps^2);
Vr=@(x,y,t) -K*cos(atan2(y,x))/(x^2+y^2+eps^2);
Vtheta=@(x,y,t) -K*sin(atan2(y,x))/(x^2+y^2+eps^2)+VB(x,y,t);

V=@(x,y,t) [
    Vr(x,y,t)*cos(atan2(y,x))-Vtheta(x,y,t)*sin(atan2(y,x))+uinf...
    Vr(x,y,t)*sin(atan2(y,x))+Vtheta(x,y,t)*cos(atan2(y,x))];
% V=@(x,y) [
%     Vr(x,y)*cos(atan2(y,x))-Vtheta(x,y)*sin(atan2(y,x))...
%     Vr(x,y)*sin(atan2(y,x))+Vtheta(x,y)*cos(atan2(y,x))];

ttr=parspace/uinf;
noftr=round(ttr*fps);
%% time loop
stpperframe=1;
h=1/stpperframe/fps;
thres=1;
f=0;
for t=0:1/fps:T
    f=f+1;
    t
    % size(particles)
    
    % if t>5
    % % omega=1/15*t-5/15;
    % % B=2*pi()*R^2*omega;
    % % VB=@(x,y) B/2/pi()/sqrt(x^2+y^2+eps^2);
    % end
%% particle creation
    if mod(f,noftr)<0.1
        for p=-boundsy:parspace:boundsy+upplus
            particle.x=-boundsx; particle.y=p;
            particles=[particles particle];
        end
        % for p=-boundsx:parspace:0
        %     particle.x=p; particle.y=boundsy+upplus;
        %     particles=[particles particle];
        % end
    end

%% pre plot
    clf;
    axis equal
    xlim([-simwinx,simwinx]);
    ylim([-simwiny,simwiny]);
    hold on;
%% bool checks
    good=[particles.x] < simwinx;
    particles = particles(good);

    bad=[particles.x] > 0 & [particles.y] > boundsy;
    good = ~bad;
    particles = particles(good);
    N=length(particles);
%% RK4 update
    for i=1:N
        x=particles(i).x;
        y=particles(i).y;
    for j=1:stpperframe
        if x^2+y^2<thres^2&&x>0
            r=0.1+0.3*rand;
            theta=2*rand*pi();
            x=r*cos(theta);
            y=r*sin(theta);
        end
            
        U=V(x,y,t);
        u=U(1);v=U(2);
        k1=V(x,y,t)*h;
        k2=V(x+0.5*k1(1),y+0.5*k1(2),t)*h;
        k3=V(x+0.5*k2(1),y+0.5*k2(2),t)*h;
        k4=V(x+k3(1),y+k3(2),t)*h;
        x=x+(k1(1)+2*k2(1)+2*k3(1)+k4(1))/6;
        y=y+(k1(2)+2*k2(2)+2*k3(2)+k4(2))/6;
        
    end

        particles(i).x=x;
        particles(i).y=y;
    end
%% for source
    % newparticle=struct('x',01*(rand()-0.5),'y',01*(rand()-0.5));
    % particles=[particles,newparticle];
    % N=N+1;
%% plot
    viscircles([0,0],R);
    plot(0,0,'r.')
    txt=['\omega = ' sprintf('%.2f', abs(omega(t)))];
    
    text(8,-8,txt,"FontSize",40)
    set(gca,'FontSize',20)
    
    wheel=wheel+omega(t)/fps;
    plot(R*cos(wheel),R*sin(wheel),'.','LineWidth',3,'MarkerSize',23,'Color',[0.8500 0.3250 0.0980])
    scatter([particles.x],[particles.y],16,'blue','filled')
%% frame

    frame = getframe(gcf);
    % Write the frame to the video file
    if t == 0
        % Create a new video file
        writerObj = VideoWriter(name,'MPEG-4');
        writerObj.FrameRate = 30;
        open(writerObj);
    else
        writeVideo(writerObj, frame);
    end

end
close(writerObj);
%% 
% for t=0:0.1:150
% plot(t,omega(t))
% hold on
% end