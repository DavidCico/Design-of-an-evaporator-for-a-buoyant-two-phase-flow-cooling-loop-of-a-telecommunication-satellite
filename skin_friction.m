%% Parameter X calculus
 
function X=skin_friction(Jl_LM,Jg_LM,rhol,rhog,mul,mug,D)



%% skin friction coefficients, fpl et fpg calculation
    Rel=Jl_LM*D*rhol/mul;
    Reg=Jg_LM*D*rhog/mug;


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
        fpg = 10^log10f;
    else
        fpg=K_g.*(Reg).^(-n_g);
    end
    
    if Rel>2000 && Rel<=3000
        log10f = ( log10( (0.079*3000^(-0.25))/(16*2000^(-1)) ) / log10(3000/2000) ) * log10(Rel/2000) + log10(16/2000);
        fpl = 10^log10f;
    else
        fpl=K_l.*(Rel).^(-n_l);
    end
   
    

    X=Jl_LM./Jg_LM.*sqrt(rhol*fpl/(rhog*fpg));