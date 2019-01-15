% Calculate pressure gradient with Awad model


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

P_Awad_l=ones(length(x),1);
P_Awad_l(1)=P0;
P_Awad_g=ones(length(x),1);
P_Awad_g(1)=P0;

% skin friction Awad

%  p parameter of Awad definition

p = 2/7;

for i=2:length(x)
    x_Awad=X;

    Jl_Awad=G.*(1-x_Awad)/rhol;
    Jg_Awad=G.*x_Awad/rhog;

    Rel=Jl_Awad*D*rhol/mul;
    Reg=Jg_Awad*D*rhog/mug;
    
    if (Rel<=2000) && (Reg<=2000)
        n_l=1;
        K_l=16;
        n_g=1;
        K_g=16;
    elseif (Rel>3000) && (Reg<=2000)
        n_l=0.25;
        K_l=0.079;
        n_g=1;
        K_g=16;
    elseif (Rel<=2000) && (Reg>3000)
        n_l=1;
        K_l=16;
        n_g=0.25;
        K_g=0.079;
    elseif (Rel>3000) && (Reg>3000)
        n_l=0.25;
        K_l=0.079;
        n_g=0.25;
        K_g=0.079;
    elseif (Rel<=2000) && ((Reg>2000) && (Reg<=3000))
        n_l=1;
        K_l=16;
    elseif (Rel>3000) && ((Reg>2000) && (Reg<=3000))
        n_l=0.25;
        K_l=0.079;
    elseif (Reg<=2000) && ((Rel>2000) && (Rel<=3000))
        n_g=1;
        K_g=16;
    elseif (Reg>3000) && ((Rel>2000) && (Rel<=3000))
        n_g=0.25;
        K_g=0.079;      
    end
    

    if (Reg>2000) && (Reg<=3000)
        log10f = ( log10( (0.079*3000^(-0.25))/(16*2000^(-1)) ) / log10(3000/2000) ) * log10(Reg/2000) + log10(16/2000);
        fpg_Awad = 10^log10f;
     else
        fpg_Awad=K_g.*(Reg).^(-n_g);
    end
    
    if (Rel>2000) && (Rel<=3000)
        log10f = ( log10( (0.079*3000^(-0.25))/(16*2000^(-1)) ) / log10(3000/2000) ) * log10(Rel/2000) + log10(16/2000);
        fpl_Awad = 10^log10f;
     else
        fpl_Awad=K_l.*(Rel).^(-n_l);
    end

    X=Jl_Awad./Jg_Awad.*sqrt(rhol*fpl_Awad/(rhog*fpg_Awad));

    Phi2_l=(1+(1/(X^2))^p)^(1/p);
    Phi2_g=(1+(X^2)^p)^(1/p);

    dPdz_l=-4/D*fpl_Awad*rhol*Jl_Awad^2/2;
    dPdz_g=-4/D*fpg_Awad*rhog*Jg_Awad^2/2;

    dPdz_frl(i)=Phi2_l*dPdz_l;
    dPdz_frg(i)=Phi2_g*dPdz_g;

    P_Awad_l(i)=P_Awad_l(i-1)+dPdz_frl(i)*(Zz(i)-Zz(i-1));
    P_Awad_g(i)=P_Awad_g(i-1)+dPdz_frg(i)*(Zz(i)-Zz(i-1));
    
end




hold  on
plot(Zz,P_Awad_l,'r')
hold off
