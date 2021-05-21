Code
Part 1
clear all
close all
clc
 
%Problem 5 
disp('Blade Element Momentum Iterator for a wind turbine')
disp('to be used to build up a whole blade BEM model')
disp('By Diego Ruiz & Ignacio Losada')
 
 
%Given: 
 
B = 3;
TSR = 8;
alpha = 6.1;
STA = 0;
Cl_max = 1.23;
R = 67.11;
 
%Solution. PROBLEM 4 Homework #5
 
r_R = (0.05:0.1:0.95);
LSR = TSR*r_R;
 
% if the twist angle is the sum of the pitch angle and the initial pitch
%angle where twist angle is zero, the initial pitch angle is -1.6.
 
alpha1=zeros(1,length(r_R));
alpha1(1,:) = 6.1;
 
%Solution. PROBLEM 1 Homework #6
 
phi_6 = (2/3)*atand(1./LSR);
r_6 = r_R*R;
c_6 = ((8*pi*r_6)/(B*Cl_max)).*(1-cosd(phi_6));
SP_6 = phi_6-alpha;
ST_6 = SP_6-SP_6(1,10)
 
%figure 6
Solidity = (B*c_6)./(2*pi*r_6);
a_6 =(4*sind(phi_6).^2);
a_6_b = (Solidity.*cosd(phi_6));
a_6 = 1./(1+(a_6)./a_6_b);
a(1:length(r_R)) = 0.3;
 
fig1 = figure;
c_6_R = c_6./R;
xlabel('r/R');ylabel('c/R');hold on;grid on
plot(r_R,c_6_R,'LineWidth',2);grid on
legend('with wake')
 
fig2 = figure;
grid on;hold on;
plot(r_R,ST_6,'LineWidth',2);grid on;
xlabel('r/R')
ylabel('Angle^\circ')
legend('Twist angle^\circ with wake')
 
fig3 = figure;
grid on;hold on;
plot(r_R,phi_6,'LineWidth',2);grid on;
xlabel('r/R')
ylabel('Angle^\circ')
legend('Angle of relative wind^\circ with wake');
 
 
fig4 = figure;
grid on;hold on;
plot(r_R,SP_6,'LineWidth',2);
xlabel('r/R')
ylabel('Angle^\circ')
legend('Section pitch angle^\circ with wake')
 
 
 
fig5 = figure;
 
grid on; hold on
plot(r_R,alpha1,'LineWidth',2);grid on
xlabel('r/R');
ylabel('angle of attack')
legend('Angle of Attack^\circ with wake')
 
fig6 = figure;
 
 
grid on;hold on;
plot(r_R,a_6,'LineWidth',2);grid on
legend('a with wake')
xlabel('r/R');
ylabel('Induction factor')
 
 
fig7 = figure;
 
a_prime(1:length(r_R)) = 0.0;
a_prime_6 = (1-3*a_6)./((4*a_6)-1);
 
grid on;hold on;
plot(r_R,a_prime_6,'LineWidth',2);grid on
xlabel('r/R');
ylabel('Induction factor')
legend('a_(_p_r_i_m_e_) with wake')
 
 
table = [r_6;c_6_R;ST_6;phi_6;SP_6;alpha1;a_6;a_prime_6]'
 
filename = 'designoftheblade.xlsx';
xlswrite(filename,table)

Part 2
clear all
clc
close all
 
 
%Inicialitation. % Blade Shape Optimum Rotor with Wake Rotation.
 
% DATA for 1 N load.
filename ='designoftheblade.xlsx';
num1 = xlsread(filename);
 
r_R_311 = (0.05:0.1:0.95)';
 
B = 3; %# number of blades
 
TSR_311 = 8; % Tip Speed Ratio
 
LSR_311 = TSR_311*r_R_311;
 
 
R = 67.11;
r_311 = r_R_311*R;
Chord_311_R = num1(1:end,2);
Chord_311 = Chord_311_R.*R;
STA_311_deg = num1(1:end,3);
STA_311_rad=STA_311_deg*pi/180;
SP_311_deg = num1(1:end,5);
SP_311_rad = SP_311_deg*pi/180;
Cl_max_311 = 1.23;
Solidity_311 = B*Chord_311./(2*pi*r_311);
 
 

 
 
 
 
relax=0.1;
converge = 0.0001;
m=0;
ite = 1;
for TSR=1:14
 
    m=m+1;
    TSR_cp(m)=TSR;
    
    r_R = 0.05:0.1:0.95;
    LSR= (TSR_cp(m)*r_R)';
    
%INITITALIZATION.
%     a(i,1) = 1./(1+((4*sin(phi_6_rad).^2)./(Cl_max_6*Solidity_6.*cos(phi_6_rad))))
%     a_prime(i,1) = (1-3*a(i,1))./((4.*a(i,1))-1)
 
    error_axial = 1;
    error_rotation = 1;
        a(1:10,1)=1/3;
        a(1:10,2)=0;
        a_prime(1:10,1)=0;
        a_prime(1:10,2)=0;
%Iterations
 
    for i=1:10
 
        n = 1;
 
 
        while error_axial(n)>converge && error_rotation(n)>converge 
 
            
            
        phi_ite(i,ite) = atan(((1-a(i,ite))./(LSR(i,ite).*(1+a_prime(i,ite)))));
 
        F(i,ite) = (2/pi)*acos(exp(-(((B/2)*(1-r_R(ite,i)))'./(r_R(ite,i)'.*sin(phi_ite(i,ite))))));
 
        alpha(i,ite) = phi_ite(i,ite)-SP_311_rad(i,ite);
 
        alpha_deg(i,ite) = (alpha(i,ite)*180)/pi;

if alpha_deg(i,ite) >= -90 && alpha_deg(i,ite) < -20.2
              Cl_ite(i,ite) = 0*(alpha_deg(i,ite)^6)- 5.1692355828E-09*(alpha_deg(i,ite)^5)...
                  - 1.6355119057E-06*(alpha_deg(i,ite)^4)- 1.9133970760E-04*(alpha_deg(i,ite)^3)...
                  - 9.5705220382E-03*(alpha_deg(i,ite)^2)- 1.8165835654E-01*(alpha_deg(i,ite))...
                  - 1.6588452104;
              Cd_ite(i,ite) = 0*(alpha_deg(i,ite)^6)- 1.1152917967E-08*(alpha_deg(i,ite)^5)...
                  - 2.8195257909E-06*(alpha_deg(i,ite)^4)- 2.5518305432E-04*(alpha_deg(i,ite)^3)...
                  - 1.0185092709E-02*(alpha_deg(i,ite)^2)- 2.1598008964E-01*(alpha_deg(i,ite))...
                  - 1.6032899751;
            elseif alpha_deg(i,ite) >= -20.2 && alpha_deg(i,ite) < 15.2 
              Cl_ite(i,ite) = + 5.8680602716E-08*(alpha_deg(i,ite)^6)+ 1.1307018778E-06*(alpha_deg(i,ite)^5)...
                  - 2.2279581308E-05*(alpha_deg(i,ite)^4)- 5.5537885629E-04*(alpha_deg(i,ite)^3)...
                  + 2.0710034350E-03*(alpha_deg(i,ite)^2)+ 1.2284570320E-01*(alpha_deg(i,ite))...
                  + 6.4951257976E-02;
              Cd_ite(i,ite) = - 1.2376869750E-08*(alpha_deg(i,ite)^6)- 1.0239769291E-07*(alpha_deg(i,ite)^5)...
                  + 6.4410496667E-06*(alpha_deg(i,ite)^4)+ 1.4537948353E-05*(alpha_deg(i,ite)^3)...
                  - 3.9178699176E-04*(alpha_deg(i,ite)^2)- 1.3992764977E-04*(alpha_deg(i,ite))...
                 + 1.4782368968E-02;
            elseif alpha_deg(i,ite) >= 15.2 && alpha_deg(i,ite) < 30.2
              Cl_ite(i,ite) = - 4.906789560853660E-10*(alpha_deg(i,ite)^6)+ 1.442768813298790E-07*(alpha_deg(i,ite)^5)...
                  - 1.608493976865220E-05*(alpha_deg(i,ite)^4)+ 8.360188208672010E-04*(alpha_deg(i,ite)^3)...
                  - 2.022769271535290E-02*(alpha_deg(i,ite)^2)+ 2.062537651900840E-01*(alpha_deg(i,ite))...
                  + 0.105;
              Cd_ite(i,ite) = + 1.466647557343170E-10*(alpha_deg(i,ite)^6)- 3.605590522578590E-08*(alpha_deg(i,ite)^5)...
                  + 3.234069226504270E-06*(alpha_deg(i,ite)^4)- 1.384701212572280E-04*(alpha_deg(i,ite)^3)...
                  + 3.432617357852050E-03*(alpha_deg(i,ite)^2)- 2.068911716673940E-02*(alpha_deg(i,ite))...
                  + 1.17e-02;
            elseif alpha_deg(i,ite) >= 30.2 && alpha_deg(i,ite) <= 90
              Cl_ite(i,ite) = 0*(alpha_deg(i,ite)^6)- 1.351495434687250E-08*(alpha_deg(i,ite)^5)...
                  + 4.323360471115850E-06*(alpha_deg(i,ite)^4)- 5.289720600966720E-04*(alpha_deg(i,ite)^3)...
                  + 3.018583104282480E-02*(alpha_deg(i,ite)^2)- 7.916138518339260E-01*(alpha_deg(i,ite))...
                  + 8.636443181225;
              Cd_ite(i,ite) = 0*(alpha_deg(i,ite)^6)+ 1.029729429071550E-08*(alpha_deg(i,ite)^5)...
                  - 2.480789504877290E-06*(alpha_deg(i,ite)^4)+ 2.042689103922690E-04*(alpha_deg(i,ite)^3)...
                  - 6.541889658905120E-03*(alpha_deg(i,ite)^2)+ 9.175984462513040E-02*(alpha_deg(i,ite))...
                  + 1.17e-02;
            elseif alpha_deg(i,ite) <-90 || alpha_deg(i,ite) > 90
                disp 'alpha = '; disp (alpha_deg(i,ite));
                disp ' i = '; disp (i);
            end
 
        Ct(i,ite) = Solidity_311(i,ite).*((1-a(i,ite)).^2).*(Cl_ite(i,ite).*cos(phi_ite(i,ite))+Cd_ite(i,ite).*sin(phi_ite(i,ite)))...
        ./(sin(phi_ite(i,ite)).^2);
 
    
            if Ct(i,ite)<0.96
        
                 a(i,ite+1) = 1./(1+((4.*F(i,ite).*sin(phi_ite(i,ite)).^2)./(Cl_ite(i,ite).*Solidity_311(i,ite).*cos(phi_ite(i,ite)))));
    
            elseif  Ct(i,ite)>0.96
        
                 a(i,ite+1) = (1./F(i,ite)).*(0.143+sqrt(0.0203-0.6427*(0.889-Ct(i,ite))));
            end
 
                 a_prime(i,ite+1) = 1./(((4.*F(i,ite).*cos(phi_ite(i,ite)))./(Cl_ite(i,ite).*Solidity_311(i,ite)))-1);
 
            n = n+1;
 
            error_axial(n) = abs((a(i,ite)-a(i,ite+1)));
            error_rotation(n) = abs((a_prime(i,ite)-a_prime(i,ite+1)));
 
            a_diff = a(i,ite+1)-(a(i,ite));
            a_prime_diff = ((a_prime(i,ite+1))-(a_prime(i,ite)));
 
 
a(i,ite) = a(i,ite)+relax*a_diff;
a_prime(i,ite) = a_prime(i,ite)+relax*a_prime_diff;
 
 
 
 
 
 
end
 
phi_cp(i,m)=phi_ite(i,ite);
cd_cp(i,m) = Cd_ite(i,ite);
cl_cp(i,m) = Cl_ite(i,ite);
LSR_cp(i,m) = LSR(i,ite);
F_cp(i,m) = F(i,ite);
 
 
end
end
 
 
 
Cp= (F_cp.*(sin(phi_cp).^2))...
.*(cos(phi_cp)-LSR_cp.*sin(phi_cp))...
.*(sin(phi_cp)+LSR_cp.*cos(phi_cp))...
.*((1-(cd_cp./cl_cp).*cot(phi_cp)).*LSR_cp.^2);
 
for m =1:length(TSR_cp)
  Cp_final(m) = (8/(TSR_cp(m)*10))*sum(Cp(:,m));
end 
 
table=[r_R_311 Chord_311_R STA_311_deg phi_ite*360/(2*pi) SP_311_deg alpha_deg a(:,2) a_prime(:,2) F r_311 Chord_311 Cl_ite Cd_ite Cp*(8/(TSR_cp(m)*10)) Ct];
 
z=[Cp_final(1) Cp_final(6)]
x2=[TSR_cp(1) TSR_cp(6)]
p = polyfit(x2,z,3);
y1 = polyval(p,TSR_cp(1:6));
Cp_final(1:6)=y1
 
plot(TSR_cp,Cp_final,'o','LineWidth',2);hold on
plot(TSR_cp,Cp_final)
title('Power Curve vs Tip Speed Ratio')
ylabel('Cp')
xlabel('Tip speed ratio');
 
% %%%%%Initialization%%%%%
% 
nu = 0.8;   %drivetrain losses.
R = 67.11;     %rotor radio in meters.
rho = 1.225; %density at sea level in kg/m^3.
 
 
 
j=0;
i=0;
    
    
for U = 3.5:0.5:25
%for U = 2:2:16
 
    j=j+1;
  U_matrix(j)=U;
for  TSR=1:14
    i=i+1 
    
    Omega(i,j) = (TSR*U)/R;
    Omega_RPM(i,j)=(Omega(i,j)*60)/(2*pi);
    P(i,j) = Cp_final(1,TSR)*0.5*rho*(pi*R^2)*U^3;
 
end
i=0;
end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%g)
 
for i = 1:length(P(1,:))
    
    P_max(i) = max(P(:,i))
    n = find(P(:,i) == P_max(i))
    Omega_RPM_max(i)=Omega_RPM(n,i);
     
    Omega_max(i)=Omega(n,i);
 
end
TSR_max = (Omega_max*R)/U
 
u_plot=[0,3.5,U_matrix(1:18),25,25.001]
p_plot=[0,0,P_max(1:18),7*10^6,0]
figure(3)
plot(u_plot,p_plot,'Linewidth',2)
xlabel('Wind Speed (m/s)')
ylabel('Power (Watts)')
 
fig2 = figure
 
hold on
 for j= 1:length(Omega_RPM(1,:))
 
    plot(Omega_RPM(:,j),P(:,j),'LineWidth',1.2);
    plot(Omega_RPM_max(j),P_max(j),'ro','LineWidth',2);
 end
  hold off
  
title('Power vs. Angular Velocity for different velocities from 3.5 to 25 m/s in incremets of 0.5 m/s ')
ylabel('Power [Watts]')
xlabel('\omega [RPM]');
 
 

table222=[P_max'/10^6 Omega_RPM_max']
 
Part 3

clear all
clc
close all
 
u=3.5:1:25
c=7.6
k=2.5
p=[];
for i=1:length(u)
    p=[p;(k/c)*(u(i)/c)^(k-1)*exp(-(u(i)/c)^k)];
end
 
figure(1)
bar(u,p)
hold on
plot(u,p','Linewidth',2,'Color','red')
xlim([3 25])
xlabel('Wind Speed (m/s)')
ylabel('Probability')
 
 
Power=0.5*1.225*pi*67.11^2*0.8*(16/27).*u.^3.*p'
CF=sum(Power)*100/(7*10^6)
