%% Function to solve the equation of vertical annular model with two fluids


function dRg=calcul_Rg_comp(z,Rg)


global Psat rhol_R245FA rhog_R245FA Hl_R245FA Hg_R245FA mul_R245FA mug_R245FA
global rhol rhog mul mug  Hl Hg hlg g q D G Z1 Z2 k fpl R lmt


dRg=zeros(3,1);

%% Parameters description

% Rg(1)=Rg       | void fraction |
% Rg(2)=x        | vapor mass fraction |
% Rg(3)=P        | pressure |


% dRg(1)=dRg/dz        | void fraction derivative in the tuibe |
% dRg(2)=dx/dz         | vapor mass fraction derivative in the tube |
% dRg(3)=dP/dz         | pressure derivative in the tube |



%% Re-calculate fluid properties in function of pressure
for i=2:length(Psat)
    if (Rg(3)<=Psat(i)) && (Rg(3)>Psat(i-1))
        rhol=(rhol_R245FA(i)-rhol_R245FA(i-1))/(Psat(i)-Psat(i-1))*(Rg(3)-Psat(i))+rhol_R245FA(i);
        rhog=(rhog_R245FA(i)-rhog_R245FA(i-1))/(Psat(i)-Psat(i-1))*(Rg(3)-Psat(i))+rhog_R245FA(i);
        mul=(mul_R245FA(i)-mul_R245FA(i-1))/(Psat(i)-Psat(i-1))*(Rg(3)-Psat(i))+mul_R245FA(i);
        mug=(mug_R245FA(i)-mug_R245FA(i-1))/(Psat(i)-Psat(i-1))*(Rg(3)-Psat(i))+mug_R245FA(i);
        Hl=(Hl_R245FA(i)-Hl_R245FA(i-1))/(Psat(i)-Psat(i-1))*(Rg(3)-Psat(i))+Hl_R245FA(i);
        Hg=(Hg_R245FA(i)-Hg_R245FA(i-1))/(Psat(i)-Psat(i-1))*(Rg(3)-Psat(i))+Hg_R245FA(i);
        hlg=Hg-Hl;
    end
end


%% Velocity calculation
Jg=G*Rg(2)/rhog;
Jl=G*(1-Rg(2))/rhol;

Ug=Jg/Rg(1);
Ul=Jl/(1-Rg(1));


%% interface drag correlation of Wallis
fi=0.005*(1+150*(1-sqrt(Rg(1))));
taui=-1/2*fi*rhog*abs(Ug-Ul)*(Ug-Ul);   


%% skin friction drag

Rel=Ul*D*rhol/mul;
Reg=Ug*D*rhog/mug;

    % condition on C and n, function of the flow

    if (Rel<=2000) && (Reg<=2000)
        n_l=1;
        K_l=16;
    elseif (Rel>3000) && (Reg<=2000)
        n_l=0.25;
        K_l=0.079;
    elseif (Rel<=2000) && (Reg>3000)
        n_l=1;
        K_l=16;
    elseif (Rel>3000) && (Reg>3000)
        n_l=0.25;
        K_l=0.079;
    elseif (Rel<=2000) && ((Reg>2000) && (Reg<=3000))
        n_l=1;
        K_l=16;
    elseif (Rel>3000) && ((Reg>2000) && (Reg<=3000))
        n_l=0.25;
        K_l=0.079;
    end
    
    if (Rel>2000) && (Rel<=3000)
        log10f = ( log10( (0.079*3000^(-0.25))/(16*2000^(-1)) ) / log10(3000/2000) ) * log10(Rel/2000) + log10(16/2000);
        fpl = 10^log10f;
    else
        fpl=K_l.*(Rel).^(-n_l);
    end

taup=-1/2*fpl*rhol*Ul^2;              


%% drop loss in the elbow pipes
Ksp=fpl*(Z2-Z1)/D+0.294*(R/D)^0.5;
Dpsp=Ksp*G^2/(2*rhol);
b=1+2.2/(Ksp*(2+R/D));
phi=1+(rhol/rhog-1)*Rg(2)*(b*(1-Rg(2))+Rg(2));
Dp_coude=phi*Dpsp*(z-Z1)/(Z2-Z1);


%% computation of equations on void fraction Rg, mass fraction x, and
% pressure loss dP/dz

    % simplification of terms for Rg
A=G^2*(((1-Rg(1))*Rg(2)^2)/(rhog*(Rg(1)^2))+(Rg(1)*(1-Rg(2))^2)/(rhol*((1-Rg(1))^2)));
B=taui*4*sqrt(Rg(1))/D;
C=Rg(1)*taup*4/D;
Dd=(rhol-rhog)*Rg(1)*(1-Rg(1))*(-g);
E=G^2*dRg(2)*(2*Rg(2)*(1-Rg(1))/(rhog*Rg(1))+(1-Rg(2))*(2*Rg(1)-1)/(rhol*(1-Rg(1))));

    % different conditions on the parameters
if Rg(1)<=lmt      % condition on Rg
    if Rg(2)<=lmt      % condition on x
        if Rg(3)>=0         % condition on P
            if q~=0             % if a component is present at the tube level
                dRg(1)=(-B+C-Dd+E)/(A);
                dRg(2)=4*q/(G*D*hlg);
                dRg(3)=-G^2/rhog*(2*Rg(2)/Rg(1)*dRg(2)-Rg(2)^2/Rg(1)^2*dRg(1))...
                    -G^2/rhol*(-2*(1-Rg(2))/(1-Rg(1))*dRg(2)-(1-Rg(2))^2/(1-Rg(1))^2*dRg(1))...
                    +4*taup/D-(rhog*Rg(1)+rhol*(1-Rg(1)))*(-g)+k*Dp_coude;
            else                % absence of component
                dRg(1)=0;
                dRg(2)=0;
                dRg(3)=-G^2/rhog*(2*Rg(2)/Rg(1)*dRg(2)-Rg(2)^2/Rg(1)^2*dRg(1))...
                    -G^2/rhol*(-2*(1-Rg(2))/(1-Rg(1))*dRg(2)-(1-Rg(2))^2/(1-Rg(1))^2*dRg(1))...
                    +4*taup/D-(rhog*Rg(1)+rhol*(1-Rg(1)))*(-g)+k*Dp_coude;
            end
        else
            if q~=0
                dRg(1)=(-B+C-Dd+E)/(A);
                dRg(2)=4*q/(G*D*hlg);
                Rg(3)=0;
            else
                dRg(1)=0;
                dRg(2)=0;
                Rg(3)=0;
            end
        end
    else
        if q~=0              % if a component is present at the tube level
            dRg(1)=(-B+C-Dd+E)/(A);
            dRg(2)=0;
            dRg(3)=0;
        else
            dRg(1)=0;
            dRg(2)=0;
            dRg(3)=0;
        end
    end
else
    if Rg(2)<=lmt      % condition on x
        Rg(1)=lmt;
        dRg(1)=0;
        dRg(2)=4*q/(G*D*hlg);
        dRg(3)=0;
    else
        Rg(1)=lmt;
        dRg(1)=0;
        dRg(2)=0;
        dRg(3)=0;
    end
end
