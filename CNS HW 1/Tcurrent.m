%% Script for HW1 Prob 4a
clear;
close all;
%% first lets plot m_inf and h_inf as a function of V
V = -100:.01:0;
% define the default parameter values
V_half_m=-70;
beta_m=1;
V_half_h=-80;
beta_h=1;
% calculate m_inf and h_inf

% if you want to do it using a function, then save the following lines in 
%the file called m_inf_fun.m in the same folder as this file
%% the m_inf function
%function [m_inf] = m_inf_fun(V_half_m, beta_m, V)
%m_inf = 1./(1+exp((V_half_m-V)/beta_m));
%end
%and then uncomment the following line here
%m_inf = m_inf_fun(V_half_m,beta_m,V);
%you can do the same thing for h as well.
%h_inf = h_inf_fun(V_half_h,beta_h,V);

%If you want to compute the values of the infnity curves directly, 
%then do this:
m_inf=1./(1+exp((V_half_m-V)/beta_m));
h_inf=1./(1+exp(-(V_half_h-V)/beta_h));


% plot the curves
plot(V,m_inf,'red','linewidth',1.5)
hold on;
plot(V,h_inf,'blue','linewidth',1.5)
% adjust plot features
legend('m_{inf}','h_{inf}')
xlabel('V (mV)')
ylabel('Activation')

%% Now lets vary V_half_m and beta_m to see the effects
figure(2)
% first, replot the default case from above
plot(V,m_inf,'r')
hold on;
% changes V_half_m to be -60

%if you used a function, you can simply recompute m_inf withe the 
%new parameter value like this
%m_inf = m_inf_fun(-60,beta_m,V);

%othewise, change V_half_m first (remember that you changed it!)
V_half_m=-60;
%recompute m_inf
m_inf=1./(1+exp((V_half_m-V)/beta_m));

%plot
plot(V,m_inf,'blue')


%Now you can, for example, vary other things such as beta_m