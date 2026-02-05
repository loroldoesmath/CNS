%Hodgkin-Huxley ODE solver. See functions file for the stimulus and
%parameters

%the vector of variables is Y=[v,m,h,n]'. Make sure the order is the same
%in the equations file 

%initial condition
in_cond=[-65 .05 0.6 .317];

%define simulation tims, in msec
t_total=200;

[t,Y]=ode15s(@hh_functions,[0 t_total],in_cond);
v=Y(:,1);
m=Y(:,2);
h=Y(:,3);
n=Y(:,4);

figure(1)
clf
plot(t,v)

figure(2)
clf
plot(t,m,'-',t,n,'--',t,h,'.')
