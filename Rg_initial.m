%% Program to calculate the initial Rg (depends on g)



function Rg2=Rg_initial(x0,g,rhol,rhog,mul,mug,m,D)

%% initialise the program
Rg1=-1;
Rg0=-1;

i=0;

      
lmt=0.99;  % model's limit of validity

G=4*m/(pi*D^2);     % surface flux

%% look for the value for cancels in the equation of Rg
while Rg1*Rg0>0
    i=i+1;
    Rg2=(i)/10000;
    
    Rg0=Rg1;

    Ug1=G*x0/(rhog*Rg2);
    Ul1=G*(1-x0)/(rhol*(1-Rg2));
    fi1=0.005*(1+150*(1-sqrt(Rg2)));
    taui1=-1/2*fi1*rhog*abs(Ug1-Ul1)*(Ug1-Ul1);
    

    Rel1=abs(Ul1)*D*rhol/mul;
    Reg1=abs(Ug1)*D*rhog/mug;

    if (Rel1<=2000) && (Reg1<=2000)
        n_l1=1;
        K_l1=16;
    elseif (Rel1>3000) && (Reg1<=2000)
        n_l1=0.25;
        K_l1=0.079;
    elseif (Rel1<=2000) && (Reg1>3000)
        n_l1=1;
        K_l1=16;
    elseif (Rel1>3000) && (Reg1>3000)
        n_l1=0.25;
        K_l1=0.079;
    elseif (Rel1<=2000) && ((Reg1>2000) && (Reg1<=3000))
        n_l1=1;
        K_l1=16;
    elseif (Rel1>3000) && ((Reg1>2000) && (Reg1<=3000))
        n_l1=0.25;
        K_l1=0.079;
    end
    
    if (Rel1>2000) && (Rel1<=3000)
        log10f = ( log10( (0.079*3000^(-0.25))/(16*2000^(-1)) ) / log10(3000/2000) ) * log10(Rel1/2000) + log10(16/2000);
        fpl1 = 10^log10f;
    else
        fpl1=K_l1.*(Rel1).^(-n_l1);
    end

    taup1=-1/2*fpl1*rhol*Ul1^2;  

    Rg1=taui1*4/D*sqrt(Rg2)-taup1*Rg2*4/D-(rhol-rhog)*Rg2*(1-Rg2)*(g);
    
   
    
end

if Rg2>lmt  % if Rg is superior to our limit
    Rg2=lmt;
end



