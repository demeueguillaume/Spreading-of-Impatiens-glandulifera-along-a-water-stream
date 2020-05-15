function [] = Impatiens_glandulifera()

clear; close all; clc;

% Physical parameters
L = 0.2*1e3;                       % Domain size (m)                         
Smax = 20;                         % Maximal population density (plants/m)
   
alpha = 10;                        % Diffusivity in absence of 
                                   % other plants (m^2/an)
beta = alpha/Smax;                 % Diffusivity contribution due to
                                   % other plants (m^3/(ind*an))

mu = 0.4;                          % rate of entry in the stream (1/an)
sigma = 1;                         % rate of return to banks (1/an)
u = 30;                            % stream advection (m/an)
k = 20;                            % stream diffusivity (m^2/an)

sPhi = 10;                         % Spread of potential for meander 
                                   % attraction (m)
umax = 5*L;                        % Max velocity in direction of 
                                   % meander (m/s)
xPhi = 0.65*L;                     % meander position (m)
Phimax = sPhi*umax;                 % Max of potential for environmental attraction (m^2/s)


% Definition of a structure containing 
% all the problem parameters
my_para = struct('L',L,'Smax',Smax,'alpha',alpha,'beta',beta,...
    'mu',mu,'sigma',sigma,'u',u,'k',k,...
    'sPhi',sPhi,'umax',umax,'xPhi',xPhi,'Phimax',Phimax);

% Numerical parameters
m = 0;                              
x = linspace(0, L ,501);       % 300 pas d'espace          
t = linspace(0,5,6);              % 10 pas de temps


% PDE solver
options = odeset;
sol = pdepe(m,@pde,@ic,@bc,x,t,...
    options,my_para);

R = sol(:,:,1); B = sol(:,:,2);

figure; set(gcf,'Color',[1 1 1])    % 2D snapshots of the solution
% -> plot with R
ax1 = gca;
set(ax1,'Xlim',[0 my_para.L],'YLim',[0 1.001*my_para.Smax],'fontsize',...
    15,'Ycolor','b','xtick',[]);
set(get(ax1,'YLabel'),'String','R(x,t)','fontsize',15)
for i=1:length(t)
    hl1 = line(x,R(i,:),'Color','b');
    if i==1
        set(hl1,'linewidth',2)
    end
end
% -> plot with B
ax2 = axes('Position',get(ax1,'Position'),'YAxisLocation','right',...
    'Color','none');
set(ax2,'Xlim',[0 my_para.L],'YLim',[0 my_para.Smax],'fontsize',15,...
    'Ycolor','r');
set(get(ax2,'YLabel'),'String','B(x,t)','fontsize',15)
set(get(ax2,'XLabel'),'String','x [miles]','fontsize',15)
for i=1:length(t)
    hl2 = line(x,B(i,:),'Color','r','Parent',ax2);
    if i==1
        set(hl2,'linewidth',2)
    end
end


% -> 2D snapshots of the solution B
figure
set(gcf,'Color',[1 1 1]); set(gca,'fontsize',15); 
plot(x,sol(1,:,2),'-b'); hold on      
for i=2:length(t)-1
    plot(x,sol(i,:,2),'--'); hold on
end
plot(x,sol(length(t),:,2),'-r');
box off; xlabel('x [m]'); ylabel('B(x,t) [ind./m]');

% -> 2D snapshots of the solution R
figure
set(gcf,'Color',[1 1 1]); set(gca,'fontsize',15); 
plot(x,sol(1,:,1),'-b'); hold on      
for i=2:length(t)-1
    plot(x,sol(i,:,1),'--'); hold on
end
plot(x,sol(length(t),:,1),'-r');
box off; xlabel('x [m]'); ylabel('R(x,t) [ind./m]');

% -> 2D snapshots of the solution B+R
figure
set(gcf,'Color',[1 1 1]); set(gca,'fontsize',15); 
plot(x,sol(1,:,2)+sol(1,:,1),'-b'); hold on      
for i=2:length(t)-1
    plot(x,sol(i,:,2)+sol(i,:,1),'--'); hold on
end
plot(x,sol(length(t),:,2)+sol(length(t),:,1),'-r');
box off; xlabel('x [m]'); ylabel('(B+R)(x,t) [ind./m]');

% -> Potential for meander attraction
figure
set(gcf,'Color',[1 1 1]); set(gca,'fontsize',15);     
plot(x,get_Phi(x,my_para));
box off; xlabel('x'); ylabel('\Phi(x)');
return

function [c,f,s] = pde(x,t,u,dudx,my_para)
% Function to define the PDE
c = [1; 1];   
K = [my_para.k; my_para.alpha + my_para.beta*u(2)]; 
V = [my_para.u + get_advvel(x,my_para); 0];
f = K .* dudx - V .* u;    
s = [my_para.mu*u(2)-my_para.sigma*u(1);...
    -my_para.mu*u(2)+my_para.sigma*u(1) + ...
    get_r(x,t,my_para)*(u(2)+u(1)) .* (1-((u(2)+u(1))/my_para.Smax))];
return

function u0 = ic(x,my_para)
% Function to define the initial condition
R0 = zeros(size(x));
B0 = 0.25 .* (x>my_para.L*0.25-5 & x<my_para.L*0.25+5) + ...
    0 .* not(x>my_para.L*0.25-5 & x<my_para.L*0.25+5);
u0 = [R0; B0];

function [pl,ql,pr,qr] = bc(xl,ul,xr,ur,t,my_para)
% Function to define the left and right bc: p+q*f = 0
pl = [0; 0]; ql = [1; 1];
pr = [0; 0]; qr = [1; 1];

function Phi = get_Phi(x,my_para)
% Function to define the potential of env. attraction
%
Phi = (-my_para.Phimax * exp(-(x-my_para.xPhi)...
    .*(x-my_para.xPhi)/(my_para.sPhi*my_para.sPhi))).* (x<=my_para.xPhi) ...
    + 0 .* not(x<=my_para.xPhi);
return

function V = get_advvel(x,my_para)
% Function to define the advection vel. V = -dPhi/dx
V = (x-my_para.xPhi) .* get_Phi(x,my_para)...
    / (my_para.sPhi*my_para.sPhi);
return

function r = get_r(x,t,my_para)
% Function to define r
r = 2;                      % no intervention
% intervention year 2
%r = -1 .*(t<=3 & t>=2) + 2 .* not(t<=3 & t>=2);
% intervention year > 2
%r = 2 .* (t<=2) + -1 .* not(t<=2);
% intervention year 4, x > 120m
%r = -1 .* (x>120 & t>4 & t<5) + 2.* not(x>120 & t>4 & t<5);

return