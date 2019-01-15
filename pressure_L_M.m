 %% Calculate pressure gradient with Lockhart & Martinelli model
% Pressure gradient L&M must be done once the main program
% has been run once, to get x and Z values.

% If a comparison between the different approximations is wanted,
% first launch L&M, then the others, to have the legend and all
% figures on the same graph.

close all


% to no slow the compuational process, fluid properties are not
% recalculated.
rhol=1282.5;
rhog=16.345;
mul=0.000317;
mug=1.0928e-5;
sigma=0.01148;
hlg=(438.75-260)*1000;
lambda=0.0839;
Cp=1385.5;
P0=3e5;


m=0.0092;                   % mass debit range : from 5 to 30 g/s
D=12e-3;

G=4*m/(pi*D^2);

dPdz_frl=ones(length(x)-1,1);
dPdz_frg=ones(length(x)-1,1);

P_LMl=ones(length(x),1);
P_LMl(1)=P0;
P_LMg=ones(length(x),1);
P_LMg(1)=P0;

% skin friction L&M

for i=2:length(x)
    x_LM=x(i);

    Jl_LM=G.*(1-x_LM)/rhol;
    Jg_LM=G.*x_LM/rhog;

    Rel=Jl_LM(i)*D*rhol/mul;
    Reg=Jg_LM(i)*D*rhog/mug;


    if (Rel<=2000) && (Reg<=2000)
        C=5;
        n_l=1;
        K_l=16;
        n_g=1;
        K_g=16;
    elseif (Rel>3000) && (Reg<=2000)
        C=10;
        n_l=0.25;
        K_l=0.079;
        n_g=1;
        K_g=16;
    elseif (Rel<=2000) && (Reg>3000)
        C=12;
        n_l=1;
        K_l=16;
        n_g=0.25;
        K_g=0.079;
    elseif (Rel>3000) && (Reg>3000)
        C=20;
        n_l=0.25;
        K_l=0.079;
        n_g=0.25;
        K_g=0.079;
    elseif (Rel<=2000) && ((Reg>2000) && (Reg<=3000))
        C=((12-5)/log10(3000/2000))*log10(Reg/2000)+5;
        n_l=1;
        K_l=16;
    elseif (Rel>3000) && ((Reg>2000) && (Reg<=3000))
        C=((20-10)/log10(3000/2000))*log10(Reg/2000)+10;
        n_l=0.25;
        K_l=0.079;
    elseif (Reg<=2000) && ((Rel>2000) && (Rel<=3000))
        C=((10-5)/log10(3000/2000))*log10(Rel/2000)+5;
        n_g=1;
        K_g=16;
    elseif (Reg>3000) && ((Rel>2000) && (Rel<=3000))
        C=((20-12)/log10(3000/2000))*log10(Rel/2000)+12;
        n_g=0.25;
        K_g=0.079;
    elseif ((Reg>2000) && (Reg<=3000)) && ((Rel>2000) && (Rel<=3000))
        if Reg>=Rel
            C=((20-12)/log10(3000/2000))*log10(Rel(i)/2000)+((12-5)/log10(3000/2000))*log10(Reg/3000)+12;
        else
            C=((10-5)/log10(3000/2000))*log10(Rel(i)/3000)+((20-10)/log10(3000/2000))*log10(Reg/2000)+10;
        end
        
    end

    if (Reg>2000) && (Reg<=3000)
        log10f = ( log10( (0.079*3000^(-0.25))/(16*2000^(-1)) ) / log10(3000/2000) ) * log10(Reg(i)/2000) + log10(16/2000);
        fpg_LM = 10^log10f;
    else
        fpg_LM=K_g.*(Reg).^(-n_g);
    end
    
    if (Rel>2000) && (Rel<=3000)
        log10f = ( log10( (0.079*3000^(-0.25))/(16*2000^(-1)) ) / log10(3000/2000) ) * log10(Rel(i)/2000) + log10(16/2000);
        fpl_LM = 10^log10f;
    else
        fpl_LM=K_l.*(Rel).^(-n_l);
    end
   
    

    X=Jl_LM(i)./Jg_LM.*sqrt(rhol*fpl_LM/(rhog*fpg_LM));

    Phi2_l=(1+C/X+1/X^2);
    Phi2_g=(1+C*X+X^2);

    dPdz_l=-4/D*fpl_LM*rhol*Jl_LM^2/2;
    dPdz_g=-4/D*fpg_LM*rhog*Jg_LM^2/2;

    dPdz_frl(i)=Phi2_l*dPdz_l;
    dPdz_frg(i)=Phi2_g*dPdz_g;

    P_LMl(i)=P_LMl(i-1)+dPdz_frl(i)*(Z(i)-Z(i-1));
    P_LMg(i)=P_LMg(i-1)+dPdz_frg(i)*(Z(i)-Z(i-1));
    
    DP_LM=P_LMl(i)-P_LMg(i);
end



figure (1)
plot(Zz,P_LMl)
xlabel('Length of tubing')
ylabel('Pressure P')
title('Evolution of pressure in tubing with approximation')


