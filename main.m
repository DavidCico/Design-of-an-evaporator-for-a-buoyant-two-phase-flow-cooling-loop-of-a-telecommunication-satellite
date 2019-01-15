%% Main program

% This program contains all other sub programs.
% It calculates the void fraction, vapor mass fraction, pressure, 
% and temperature inside the tubes (calculation for one tube)

% We calculate these values, with and without gravity


close all
clear all


global Tsat Psat rhol_R245FA rhog_R245FA Hl_R245FA Hg_R245FA Cpl_R245FA Cpg_R245FA mul_R245FA mug_R245FA lambdal_R245FA lambdag_R245FA sigma_R245FA
global R_contact
global rhol rhog mul mug sigma Hg Hl hlg lambda Cpl Cpg g q m D G Z1 Z2 Tf P0 k fpl R lmt
global x Zz P ggeom qgeom N Jll Jgg Q_geom L_geom l_geom Np_geom

sim=input('What type of simulation do you want? (2 = with and without gravity / 1 = with gravity / 0 = without gravity)  ');




m=0.01;             % mass flux : from 5 to 30 g/s (sum of the 3 tubes)

D=12e-3;            % tubes diameter

x0=0.3;             % initial vapor mass fraction between 0 and 0.3

Tf=45;              % initial fluid temperature

G=4*m/(pi*D^2);     % surface flux

P0=0.3e6;           % saturation vapor pressure at 45 degsC

lmt=0.99;           % validity limit of the models

%% loading fluid properties
R245FA;



%% loading thermal properties
R_contact=5000;     % contact resistance between the imitator and the fin

%% geometry definition
geom=importdata('elbow_geometry.txt');

Zgeom1=geom(:,1)/1000;  % geometry position
Zgeom2=geom(:,2)/1000;

k_coude=geom(:,3);      % elbow's presence

qgeom=geom(:,4);        % surface flux
ggeom=geom(:,5);        % gravity value

Q_geom=geom(:,6);       % component power
L_geom=geom(:,7);       % section's length
l_geom=geom(:,8);       % width of the component
Np_geom=geom(:,9);      % number of passages at the component level


    % geometry initialisation
k=k_coude(1);
q=qgeom(1);
g=ggeom(1);

Z1=Zgeom1(1);
Z2=Zgeom2(1);

R=(Zgeom2(1)-Zgeom2(1))/pi;     % radius of curvature (if existent)
            

if (sim==1) || (sim==2)

for i=2:length(Tsat)
    if (Tf<=Tsat(i)) && (Tf>Tsat(i-1))
        rhol=(rhol_R245FA(i)-rhol_R245FA(i-1))/(Tsat(i)-Tsat(i-1))*(Tf-Tsat(i))+rhol_R245FA(i);
        rhog=(rhog_R245FA(i)-rhog_R245FA(i-1))/(Tsat(i)-Tsat(i-1))*(Tf-Tsat(i))+rhog_R245FA(i);
        mul=(mul_R245FA(i)-mul_R245FA(i-1))/(Tsat(i)-Tsat(i-1))*(Tf-Tsat(i))+mul_R245FA(i);
        mug=(mug_R245FA(i)-mug_R245FA(i-1))/(Tsat(i)-Tsat(i-1))*(Tf-Tsat(i))+mug_R245FA(i);
        Hl=(Hl_R245FA(i)-Hl_R245FA(i-1))/(Tsat(i)-Tsat(i-1))*(Tf-Tsat(i))+Hl_R245FA(i);
        Hg=(Hg_R245FA(i)-Hg_R245FA(i-1))/(Tsat(i)-Tsat(i-1))*(Tf-Tsat(i))+Hg_R245FA(i);
        hlg=Hg-Hl;
    end
end

%% initial Rg calculation
Rg0=Rg_initial(x0,g,rhol,rhog,mul,mug,m,D);


%% solving equations with gravity
[Z Rg]=ode45(@calcul_Rg_comp, [Zgeom1(1) Zgeom2(1)], [Rg0 x0 P0]);

Zz=Z;
Rgg=Rg(:,1);
x=Rg(:,2);
P=Rg(:,3);

Rg0=Rgg(length(Zz));
x0=x(length(Zz));
P0=P(length(Zz));

N(1)=length(Z);



    % calculate values for all geometry
for j=2:length(Zgeom1)
    q=qgeom(j);
    g=ggeom(j);
    k=k_coude(j);
    Z1=Zgeom1(j);
    Z2=Zgeom2(j);
    
    R=(Zgeom2(j)-Zgeom2(j))/pi;
    
    % condition on change of sign of  g => recalculate initial Rg
    if (ggeom(j)*ggeom(j-1)==0) && (Rg0<lmt)
        Rg0=Rg_initial(x0,g,rhol,rhog,mul,mug,m,D);
    end
    
    [Z Rg]=ode45(@calcul_Rg_comp, [Zgeom1(j) Zgeom2(j)], [Rg0 x0 P0]);

    
    n1=length(Zz);
    n2=length(Z);
    
    Rg0=Rg(n2,1);
    x0=Rg(n2,2);
    P0=Rg(n2,3);
    
    
    N(j)=length(Z);     % number of iterations in each g
    
    for ij=n1+1:n1+n2
        Zz(ij)=Z(ij-n1);       
        Rgg(ij)=Rg(ij-n1,1);
        x(ij)=Rg(ij-n1,2);
        P(ij)=Rg(ij-n1,3);
    end


end

Z_g=Zz;
x_g=x;       % to leave x free
Rg_g=Rgg;
P_g=P;
N_g=N;

    %% velocity calculation

Jgg=G.*x./rhog;
Jll=G.*(1-x)./rhol;
Ugg=Jgg./Rgg;
Ull=Jll./(1.-Rgg);

    %% temperature computation
[Tp_S_G Tc_S_G h_S_G]=temp_S_G;
[Tp_K Tc_K h_K]=temp_K;
[Tp_G_W Tc_G_W h_G_W]=temp_G_W;


%% compute pressure and temperature loss between 
%% inlet and outlet


Delta_P=P_g(1)-P_g(length(P_g));

for i=2:length(Psat)
    if (P(length(P))<=Psat(i)) && (P(length(P))>Psat(i-1))
        Tf_min=(Tsat(i)-Tsat(i-1))/(Psat(i)-Psat(i-1))*(P(length(P))-Psat(i))+Tsat(i);
    end
end

Delta_T=45-Tf_min;

%% Tmax curves
T_max=zeros(length(x_g),1);
for i=1:length(x_g)
    if Z_g(i)<45
        T_max(i)=65;
    else
        T_max(i)=80;
    end
end

end


%% solving equations without gravity

if (sim==2) || (sim==0)
    


ggeom=zeros(length(ggeom),1);       % null gravity here

Tf=45;      % initial temperature

for i=2:length(Tsat)
    if (Tf<=Tsat(i)) && (Tf>Tsat(i-1))
        rhol=(rhol_R245FA(i)-rhol_R245FA(i-1))/(Tsat(i)-Tsat(i-1))*(Tf-Tsat(i))+rhol_R245FA(i);
        rhog=(rhog_R245FA(i)-rhog_R245FA(i-1))/(Tsat(i)-Tsat(i-1))*(Tf-Tsat(i))+rhog_R245FA(i);
        mul=(mul_R245FA(i)-mul_R245FA(i-1))/(Tsat(i)-Tsat(i-1))*(Tf-Tsat(i))+mul_R245FA(i);
        mug=(mug_R245FA(i)-mug_R245FA(i-1))/(Tsat(i)-Tsat(i-1))*(Tf-Tsat(i))+mug_R245FA(i);
        Hl=(Hl_R245FA(i)-Hl_R245FA(i-1))/(Tsat(i)-Tsat(i-1))*(Tf-Tsat(i))+Hl_R245FA(i);
        Hg=(Hg_R245FA(i)-Hg_R245FA(i-1))/(Tsat(i)-Tsat(i-1))*(Tf-Tsat(i))+Hg_R245FA(i);
        hlg=Hg-Hl;
    end
end

g=ggeom(1);            
x0=0.3;         % initial title between 0 and 0.3
P0=0.3e6;       % saturation vapor pressure at 45 degsC
Rg0=Rg_initial(x0,g,rhol,rhog,mul,mug,m,D);


[Z Rg]=ode45(@calcul_Rg_comp, [Zgeom1(1) Zgeom2(1)], [Rg0 x0 P0]);

Zz=Z;
Rgg=Rg(:,1);
x=Rg(:,2);
P=Rg(:,3);

Rg0=Rgg(length(Zz));
x0=x(length(Zz));
P0=P(length(Zz));

N(1)=length(Z);
    
for j=2:length(Zgeom1)
    q=qgeom(j);
    g=ggeom(j);
    k=k_coude(j);
    Z1=Zgeom1(j);
    Z2=Zgeom2(j);
    
    R=(Zgeom2(j)-Zgeom2(j))/pi;
    
    [Z Rg]=ode45(@calcul_Rg_comp, [Zgeom1(j) Zgeom2(j)], [Rg0 x0 P0]);

    n1=length(Zz);
    n2=length(Z);
    
    Rg0=Rg(n2,1);
    x0=Rg(n2,2);
    P0=Rg(n2,3);
    
    
    N(j)=length(Z);     % number of iterations in each g
    
    for ij=n1+1:n1+n2
        Zz(ij)=Z(ij-n1);       
        Rgg(ij)=Rg(ij-n1,1);
        x(ij)=Rg(ij-n1,2);
        P(ij)=Rg(ij-n1,3);
    end


end

Z_no_g=Zz;
x_no_g=x;
Rg_no_g=Rgg;
P_no_g=P;

    %% velocity calculation

Jgg=G.*x./rhog;
Jll=G.*(1-x)./rhol;
Ugg_no_g=Jgg./Rg_no_g;
Ull_no_g=Jll./(1.-Rg_no_g);

    %% temperature calculation
[Tp_S_G_no_g Tc_S_G_no_g h_S_G_no_g]=temp_S_G;
[Tp_K_no_g Tc_K_no_g h_K_no_g]=temp_K;
[Tp_G_W_no_g Tc_G_W_no_g h_G_W_no_g]=temp_G_W;

                                         
%% compute pressure and temperature loss between 
%% inlet and outlet


Delta_P_no_g=P_no_g(1)-P_no_g(length(P_no_g));


for i=2:length(Psat)
    if (P_no_g(length(P_no_g))<=Psat(i)) && (P_no_g(length(P_no_g))>Psat(i-1))
        Tf_min_no_g=(Tsat(i)-Tsat(i-1))/(Psat(i)-Psat(i-1))*(P_no_g(length(P_no_g))-Psat(i))+Tsat(i);
    end
end


Delta_T_no_g=45-Tf_min_no_g;
       

%% Tmax curves

T_max_no_g=zeros(length(x_no_g),1);
for i=1:length(x_no_g)
    if Z_no_g(i)<45
        T_max_no_g(i)=65;
    else
        T_max_no_g(i)=80;
    end
end

end

%% figures

if sim==1
figure(1)
plot(Z_g,Rg_g)
xlabel('Length of tubing')
ylabel('Vapour fraction Rg')
title('Evolution of vapour fraction in tubing, with gravity')


figure(2)
plot(Z_g,x_g)
xlabel('Length of tubing')
ylabel('Quality x')
title('Evolution of quality in tubing, with gravity')


figure(3)
plot(Z_g,Ugg)
hold on
plot(Z_g,Ull,'cyan')
hold off
xlabel('Length of tubing')
ylabel('Phase velocity')
title('Evolution of velocity in tubing, with gravity')
legend('vapour velocity, with gravity','liquid velocity, with gravity','location','Northwest')


figure(4)
plot(Z_g,P_g)
xlabel('Length of tubing')
ylabel('Pressure P')
title('Evolution of pressure in tubing, with gravity')


figure(5)
plot(Z_g,Tp_S_G)
hold on
plot(Z_g,Tp_G_W,'g')
plot(Z_g,Tp_K,'r')
hold off
xlabel('Length of tubing')
ylabel('Wall temperature')
title('Evolution of wall temperature in tubing, with gravity')
legend('model of Schrock & Grossman','model of Gunger & Winterton','model of Kandlikar','location','Northwest')


figure(7)
plot(Z_g,Tc_S_G)
hold on
plot(Z_g,Tc_G_W,'g')
plot(Z_g,Tc_K,'r')
plot(Z_g,T_max,'color',[1 1/2 0])
hold off
str1(1)={'Maximum temperature'};
text(20,67,str1,'color',[1 1/2 0])
xlabel('Length of tubing')
ylabel('Component temperature')
title('Evolution of component temperature in tubing, with gravity')
legend('model of Schrock & Grossman','model of Gunger & Winterton','model of Kandlikar','location','Northwest')



%% write results in a txt file

write_results;

fclose('all');
end


if sim==0
figure(1)
plot(Z_no_g,Rg_no_g,'r')
xlabel('Length of tubing')
ylabel('Vapour fraction Rg')
title('Evolution of vapour fraction in tubing, without gravity')


figure(2)
plot(Z_no_g,x_no_g,'r')
xlabel('Length of tubing')
ylabel('Quality x')
title('Evolution of quality in tubing, without gravity')


figure(3)
hold on
plot(Z_no_g,Ugg_no_g,'r')
plot(Z_no_g,Ull_no_g,'color',[1 0.5 0])
hold off
xlabel('Length of tubing')
ylabel('Phase velocity')
title('Evolution of velocity in tubing, without gravity')
legend('vapour velocity, without gravity','liquid velocity, without gravity','location','Northwest')


figure(4)
plot(Z_no_g,P_no_g,'r')
xlabel('Length of tubing')
ylabel('Pressure P')
title('Evolution of pressure in tubing, without gravity')


figure(6)
plot(Z_no_g,Tp_S_G_no_g)
hold on
plot(Z_no_g,Tp_G_W_no_g,'g')
plot(Z_no_g,Tp_K_no_g,'r')
hold off
xlabel('Length of tubing')
ylabel('Wall temperature')
title('Evolution of wall temperature in tubing, without gravity')
legend('model of Schrock & Grossman','model of Gunger & Winterton','model of Kandlikar','location','Northwest')


figure(8)
plot(Z_no_g,Tc_S_G_no_g)
hold on
plot(Z_no_g,Tc_G_W_no_g,'g')
plot(Z_no_g,Tc_K_no_g,'r')
plot(Z_no_g,T_max_no_g,'color',[1 1/2 0])
hold off
str1(1)={'Maximum temperature'};
text(20,67,str1,'color',[1 1/2 0])
xlabel('Length of tubing')
ylabel('Component temperature')
title('Evolution of component temperature in tubing, without gravity')
legend('model of Schrock & Grossman','model of Gunger & Winterton','model of Kandlikar','location','Northwest')


%% write results in a txt file

write_results_no_g;

fclose('all');
end


if sim==2
figure(1)
plot(Z_g,Rg_g)
hold on
plot(Z_no_g,Rg_no_g,'r')
hold off
xlabel('Length of tubing')
ylabel('Vapour fraction Rg')
title('Evolution of vapour fraction in tubing')
legend('with gravity','without gravity','location','Southeast')

figure(2)
plot(Z_g,x_g)
hold on
plot(Z_no_g,x_no_g,'r')
hold off
xlabel('Length of tubing')
ylabel('Quality x')
title('Evolution of quality in tubing')
legend('with gravity','without gravity','location','Southeast')

figure(3)
plot(Z_g,Ugg)
hold on
plot(Z_g,Ull,'cyan')
plot(Z_no_g,Ugg_no_g,'r')
plot(Z_no_g,Ull_no_g,'color',[1 0.5 0])
hold off
xlabel('Length of tubing')
ylabel('Phase velocity')
title('Evolution of velocity in tubing')
legend('vapour velocity, with gravity','liquid velocity, with gravity','vapour velocity, without gravity','liquid velocity, without gravity','location','Northwest')



figure(4)
plot(Z_g,P_g)
hold on
plot(Z_no_g,P_no_g,'r')
hold off
xlabel('Length of tubing')
ylabel('Pressure P')
title('Evolution of pressure in tubing')
legend('with gravity','without gravity')

figure(5)
plot(Z_g,Tp_S_G)
hold on
plot(Z_g,Tp_G_W,'g')
plot(Z_g,Tp_K,'r')
hold off
xlabel('Length of tubing')
ylabel('Wall temperature')
title('Evolution of wall temperature in tubing, with gravity')
legend('model of Schrock & Grossman','model of Gunger & Winterton','model of Kandlikar','location','Northwest')


figure(6)
plot(Z_no_g,Tp_S_G_no_g)
hold on
plot(Z_no_g,Tp_G_W_no_g,'g')
plot(Z_no_g,Tp_K_no_g,'r')
hold off
xlabel('Length of tubing')
ylabel('Wall temperature')
title('Evolution of wall temperature in tubing, without gravity')
legend('model of Schrock & Grossman','model of Gunger & Winterton','model of Kandlikar','location','Northwest')


figure(7)
plot(Z_g,Tc_S_G)
hold on
plot(Z_g,Tc_G_W,'g')
plot(Z_g,Tc_K,'r')
plot(Z_g,T_max,'color',[1 1/2 0])
hold off
str1(1)={'Maximum temperature'};
text(20,67,str1,'color',[1 1/2 0])
xlabel('Length of tubing')
ylabel('Component temperature')
title('Evolution of component temperature in tubing, with gravity')
legend('model of Schrock & Grossman','model of Gunger & Winterton','model of Kandlikar','location','Northwest')


figure(8)
plot(Z_no_g,Tc_S_G_no_g)
hold on
plot(Z_no_g,Tc_G_W_no_g,'g')
plot(Z_no_g,Tc_K_no_g,'r')
plot(Z_no_g,T_max_no_g,'color',[1 1/2 0])
hold off
str1(1)={'Maximum temperature'};
text(20,67,str1,'color',[1 1/2 0])
xlabel('Length of tubing')
ylabel('Component temperature')
title('Evolution of component temperature in tubing, without gravity')
legend('model of Schrock & Grossman','model of Gunger & Winterton','model of Kandlikar','location','Northwest')


%% write results in a txt file

write_results;
write_results_no_g;

fclose('all');
end


