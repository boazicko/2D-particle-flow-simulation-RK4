%% initialize
clear
close all
format compact

stmp = char(datetime("now"));
%stmp=stmp(length(stmp)-7:length(stmp));
stmp = strrep(stmp, ':', '_');
txt='my_movie_';
name=[txt stmp];
name=char(name);
fps=60;

boundsx=18;
boundsy=12;
simwinx=18;
simwiny=10;
noppr=55;
parspace=0.4;
upplus=0;

%%

T=45;

R=3.5;
[xPositions, yPositions] = meshgrid(-boundsx:parspace:+boundsx,-boundsy:parspace:boundsy+upplus);
N=numel(xPositions);
particle = struct('x',0,'y',0);

particles(N) = particle;
particles=[];
for i = 1:N
    particles(i).x = xPositions(i);
    particles(i).y = yPositions(i);

end
%figure('Position', [-3000, +500, 1920, 1080]);
%figure('Position', [-2400, 500, 1920, 1080]);
figure('Position', [0, -50, 1920, 1080]);
axis equal
hold on
xlim([-simwinx,simwinx]);
ylim([-simwiny,simwiny]);
wheel=linspace(0,2*pi(),9);
%% velocity function
eps=10^-1;
rmp=@(t,a,b,t1,t2) a+(b-a)/(t2-t1)*(t-t1)*double(t>t1)-(b-a)/(t2-t1)*(t-t2)*double(t>t2);

%makor
Q=100;
boost=100;
rm=0.6;
tz1=15;
tz2=18;
d1=7;
d2=0.8;
Xm=@(t) 0;
Qm=@(t) rmp(t,0,Q,0,2);
Vrm=@(x,y,t) Qm(t)/2/pi()/sqrt(x^2+y^2+eps^2);
Vm=@(x,y,t) [Vrm(x,y,t)*cos(atan2(y,x))...
             Vrm(x,y,t)*sin(atan2(y,x))];
%bor
rb=0.8;
Xb=@(t) rmp(t,d1,d2,tz1,tz2);
Qb=@(t) 0;
Vrb=@(x,y,t) -Qb(t)/2/pi()/sqrt(x^2+y^2+eps^2);
Vb=@(x,y,t) [Vrb(x,y,t)*cos(atan2(y,x))...
             Vrb(x,y,t)*sin(atan2(y,x))];

%metsifa
uinfty=4;
%uinf=@(t) rmp(t,0,uinfty,5,6)+rmp(t,0,uinfty,15,16)+rmp(t,0,uinfty,25,26)+rmp(t,0,uinfty,30,31)+rmp(t,0,uinfty,35,36)+rmp(t,0,uinfty,40,41);
uinf=@(t) rmp(t,0,25,0,40);
%zugan
%K=@(t) uinf(t)*R^2;
K=@(t) 0;

%arbol
w=0;
omega=@(t) rmp(t,0,w,15,20);
A=@(t) 2*pi()*R^2*omega(t);
VA=@(x,y,t) A(t)/2/pi()/sqrt(x^2+y^2+eps^2);

%polar
Vr=@(x,y,t) -K(t)*cos(atan2(y,x))/(x^2+y^2+eps^2);
Vtheta=@(x,y,t) -K(t)*sin(atan2(y,x))/(x^2+y^2+eps^2)+VA(x,y,t);

%kartezi
V=@(x,y,t) [
    Vr(x,y,t)*cos(atan2(y,x))-Vtheta(x,y,t)*sin(atan2(y,x))+uinf(t)...
    Vr(x,y,t)*sin(atan2(y,x))+Vtheta(x,y,t)*cos(atan2(y,x))]...
    +Vm(x-Xm(t),y,t)+Vb(x-Xb(t),y,t);


% ttr=parspace/uinf;
% noftr=round(ttr*fps);
vtest=@(x) V(x,0,22);
% figure (2)
% fplot(vtest,[-5,-1.5])
%% time loop
stpperframe=3;
h=1/stpperframe/fps;
thres=1;
f=0;
dis=0;
disturb=0.02*parspace;
for t=0:1/fps:T
    
    dis=dis+uinf(t)/fps;
    f=f+1;
    t
    %% pre plot
    clf;
    axis equal
    xlim([-simwinx,simwinx]);
    ylim([-simwiny,simwiny]);
    hold on;

%% bool checks
    good=[particles.x] < simwinx+5;
    particles = particles(good);
    good=[particles.x] > -simwinx-0.25;
    particles = particles(good);
    % bad=[particles.x] > 0 & [particles.y] > boundsy;
    % good = ~bad;
    % particles = particles(good);
    N=length(particles);

    %bor
    if Qb(t)>0
    good=([particles.x]-Xb(t)).^2+[particles.y].^2>rb^2;
    particles = particles(good);
    N=length(particles);
    end
    
%% RK4 update
    N=length(particles);
    for i=1:N
        x=particles(i).x;
        y=particles(i).y;
    for j=1:stpperframe
    
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
    %% particle creation
    if Qm(t)>0
    plot(Xm(t),0,'r*','MarkerSize',12)
    for p=1:round(Qm(t)/parspace/25)
    r=0.01+rm*rand();
    theta=2*pi()*rand;
    particle.x=r*cos(theta)+Xm(t); particle.y=r*sin(theta);
    particles=[particles particle];
    N=length(particles);
    end
    if Qb(t)>0
    plot(Xb(t),0,'r*','MarkerSize',12)
    end
    end

   if dis>parspace
       dis=0;
       disturb=-disturb;
        for p=-boundsy:parspace:boundsy+upplus
            particle.x=-boundsx; particle.y=p+disturb;
            particles=[particles particle];
        end
        N=length(particles);
        % for p=-boundsx:parspace:0
        %     particle.x=p; particle.y=boundsy+upplus;
        %     particles=[particles particle];
        % end
   end
%% plot
    color=@(t) [1 rmp(t,1,0.01,22.5,25) rmp(t,1,0.01,22.5,25)];
    if t>400
    viscircles([0,0],R,'Color',color(t),'LineWidth',5);
    wheel=wheel+omega(t)/fps;
    plot(R*cos(wheel),R*sin(wheel),'.','LineWidth',3,'MarkerSize',30,'Color',color(t))
    end
    plot(0,0,'r.')
    %txt
    txt=['U_{\infty} = ' sprintf('%.2f', abs(uinf(t)))];
    text(10.8,8,txt,"FontSize",40)
    txt=['\omega = ' sprintf('%.2f', abs(omega(t)))];
    text(12,5,txt,"FontSize",40)
    txt=['Q_{sink} = ' sprintf('%.2f', Qb(t))];
    text(8.3,-5,txt,"FontSize",40)
    txt=['Q_{source} = ' sprintf('%.2f', Qm(t))];
    text(7,-8,txt,"FontSize",40)
    set(gca,'FontSize',20)
    txt=['t = ' sprintf('%.2f', abs(t))];
    text(-17,-8,txt,"FontSize",40)    
    wheel=wheel+omega(t)/fps;
    %plot(R*cos(wheel),R*sin(wheel),'.','LineWidth',3,'MarkerSize',23,'Color',[0.8500 0.3250 0.0980])
    scatter([particles.x],[particles.y],16,'blue','filled')
%% frame

    frame = getframe(gcf);
    % Write the frame to the video file
    if f == 1
        % Create a new video file
        writerObj = VideoWriter(name,'MPEG-4');
        writerObj.FrameRate = fps;
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