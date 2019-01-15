% Calculate pressure gradient with Baroczy model



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

P_B_l=ones(length(x),1);
P_B_l(1)=P0;
P_B_g=ones(length(x),1);
P_B_g(1)=P0;

% skin friction Baroczy

for i=2:length(x)
    x_B=x(i);

    Jl_B=G.*(1-x_B)/rhol;
    Jg_B=G.*x_B/rhog;
    
    vl=G/rhol;
    vv=G/rhog;

    Rel=vl*D*rhol/mul;
    Reg=vv*D*rhog/mug;
    
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
    

    if Reg>2000 && Reg<=3000
        log10f = ( log10( (0.079*3000^(-0.25))/(16*2000^(-1)) ) / log10(3000/2000) ) * log10(Reg/2000) + log10(16/2000);
        fpg_B = 10^log10f;
     else
        fpg_B=K_g.*(Reg).^(-n_g);
    end
    
    if (Rel>2000) && (Rel<=3000)
        log10f = ( log10( (0.079*3000^(-0.25))/(16*2000^(-1)) ) / log10(3000/2000) ) * log10(Rel/2000) + log10(16/2000);
        fpl_B = 10^log10f;
        n_l=0.25;
    else
        fpl_B=K_l.*(Rel).^(-n_l);
    end

    dPdz_l=-4/D*fpl_B*rhol*Jl_B^2/2;
    
    % Baroczy specific
    
    Y=sqrt(rhol*fpg_B/(rhog*fpl_B));
    
    if (Y>0) && (Y<9.5)
        B=55/(G^0.5);
    elseif (Y>=9.5) && (Y<28)
        B=520/(Y*G^0.5);
    elseif (Y>=28)
        B=15000/(Y^2*G^0.5);
    end
          
    phi2_l0 = 1+(Y.^2-1)*(B*x(i)^((2-n_l)/n_l)*(1-x(i))^((2-n_l)/2)+x(i)^(2-n_l));

    dPdz_frl(i)=phi2_l0*dPdz_l;
    

    P_B_l(i)=P_B_l(i-1)+dPdz_frl(i)*(Zz(i)-Zz(i-1));
    
end




hold on
plot(Zz,P_B_l,'g')
hold off