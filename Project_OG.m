% Project for MAT 5800

clear;
close all;
clc;
% We can play with parameters that are (P+A)/2 = sqrt(PA)
par1 = [12,9,pi/6];
par2 = [5,2,-2*pi/3];
par = [par1;par2];
[m,n] = size(par);
syms t

for i = 1:m
   alpha = (par(i,2) - par(i,1)) / 2;
   beta = (par(i,2) + par(i,1)) / 2;
   gamma = sqrt(par(i,2) * par(i,1));
   orbit_vec(:,i) = [cos(par(i,3)), sin(par(i,3));...
       -sin(par(i,3)), cos(par(i,3))]*[alpha + beta*cos(t); ...
       gamma * sin(t)];
end

[t1,t2] = meshgrid(-8:.1:7);
x1 = subs(orbit_vec(1,1),t,t1);
y1 = subs(orbit_vec(2,1),t,t1);
x2 = subs(orbit_vec(1,2),t,t2);
y2 = subs(orbit_vec(2,2),t,t2);

dist_cont = .5 .* ((x1-x2).^2 + (y1-y2).^2);

figure(1)
contourf(t1,t2,dist_cont,15);
xlabel('t_1 points')
ylabel('t_2 points')
title('The Optimal t_1 and t_2 Points to Find The Minimum Distance Between The Two Orbits')
hold on
figure(2)
contour3(t1,t2,dist_cont,200);
xlabel('t_1 points')
ylabel('t_2 points')
title('The Contour Plot For t_1 and t_2 in 2D')


syms t t1 t2

x1 = subs(orbit_vec(1,1),t,t1);
y1 = subs(orbit_vec(2,1),t,t1);
x2 = subs(orbit_vec(1,2),t,t2);
y2 = subs(orbit_vec(2,2),t,t2);

figure(3)
fplot(orbit_vec(1,1),orbit_vec(2,1),'linewidth',2)
hold on
fplot(orbit_vec(1,2),orbit_vec(2,2),'linewidth',2)
hold on
axis([-15, 15, -15, 15])
grid on
xlabel('x-points')
ylabel('y-points')
title('The Elliptical Orbits')
legend('Orbit 1', 'Orbit 2','Location','NorthEast')


x_1=@(t1) subs(orbit_vec(1,1),t,t1);
y_1=@(t1) subs(orbit_vec(2,1),t,t1);
x_2=@(t2) subs(orbit_vec(1,2),t,t2);
y_2=@(t2) subs(orbit_vec(2,2),t,t2);

dist=@(x1,x2,y1,y2) .5 * ((x1-x2)^2 + (y1-y2)^2);

    % Creating quiver plot
    [x,y] = meshgrid(-8:.1:7);
    dist1=@(t1,t2) dist(x_1(t1),x_2(t2),y_1(t1),y_2(t2));

    for i = 1:length(x)
        for j = 1:length(y)
        dist2(i,j) = double(dist1(x(i,j),y(i,j)));
        end
    end

    dist2 = double(dist2);

    [Dt1,Dt2] = gradient(dist2,.5);
    
    figure(1)
    quiver(x,y,Dt1,Dt2)
    hold on

% I want to perform steepest decent 10 times.
num_point=4;
t1_0 = [-2, -3.2, 1.8, -1.5]; t2_0 = [2.5, -1, -.5, -2];
figure(1)
plot(t1_0(num_point),t2_0(num_point),'mh','MarkerSize',20,'LineWidth',2)
hold on

% My initial search
n = 20;
[x,y] = Project_Steepest_Decent_MAT_5800(dist, x1, x2, y1, y2,...
    x_1, x_2, y_1, y_2, t1_0(num_point), t2_0(num_point), n);
t1 = x;
t2 = y;


Opt_t1 = t1(n+1);
Opt_t2 = t2(n+1);

figure(1)
plot(t1,t2,'rd-','LineWidth',1.3)
hold on
plot(Opt_t1,Opt_t2,'gs','MarkerSize',20,'LineWidth',2)
hold on
hold off


x1 = double(subs(orbit_vec(1,1),t,Opt_t1))
y1 = double(subs(orbit_vec(2,1),t,Opt_t1))
x2 = double(subs(orbit_vec(1,2),t,Opt_t2))
y2 = double(subs(orbit_vec(2,2),t,Opt_t2))

distance = dist(x1,x2,y1,y2);

figure(3)
plot(x1,y1,'g*','LineWidth',2)
hold on
plot(x2,y2,'m*','LineWidth',2)
hold on
plot([x1,x2],[y1,y2],'k-','LineWidth',1.2)
hold on
legend('(x_1,y_1)','(x_2,y_2)','Optimal (x_1,y_1)','Optimal (x_2,y_2)',...
    'Location','NorthEast')
hold on

