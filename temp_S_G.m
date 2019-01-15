% Schrock and Grossman temperature calculation

function [Tp Tc h]=temp_S_G

global Tsat Psat rhol_R245FA rhog_R245FA Hl_R245FA Hg_R245FA Cpl_R245FA Cpg_R245FA mul_R245FA mug_R245FA lambdal_R245FA lambdag_R245FA sigma_R245FA
global R_contact
global rhol rhog mul mug hlg lambda Cpl Cpg D G Tf lmt
global x Zz P ggeom qgeom N Jll Jgg Q_geom L_geom l_geom Np_geom

Tp=ones(length(x),1);
Tc=ones(length(x),1);
qq=qgeom(1);


for j=1:N(1)
    for i=2:length(Psat)
        if (P(j)<=Psat(i)) && (P(j)>Psat(i-1))
            rhol=(rhol_R245FA(i)-rhol_R245FA(i-1))/(Psat(i)-Psat(i-1))*(P(j)-Psat(i))+rhol_R245FA(i);
            rhog=(rhog_R245FA(i)-rhog_R245FA(i-1))/(Psat(i)-Psat(i-1))*(P(j)-Psat(i))+rhog_R245FA(i);
            mul=(mul_R245FA(i)-mul_R245FA(i-1))/(Psat(i)-Psat(i-1))*(P(j)-Psat(i))+mul_R245FA(i);
            mug=(mug_R245FA(i)-mug_R245FA(i-1))/(Psat(i)-Psat(i-1))*(P(j)-Psat(i))+mug_R245FA(i);
            Hl=(Hl_R245FA(i)-Hl_R245FA(i-1))/(Psat(i)-Psat(i-1))*(P(j)-Psat(i))+Hl_R245FA(i);
            Hg=(Hg_R245FA(i)-Hg_R245FA(i-1))/(Psat(i)-Psat(i-1))*(P(j)-Psat(i))+Hg_R245FA(i);
            lambda=(lambdal_R245FA(i)-lambdal_R245FA(i-1))/(Psat(i)-Psat(i-1))*(P(j)-Psat(i))+lambdal_R245FA(i);
            Cpl=(Cpl_R245FA(i)-Cpl_R245FA(i-1))/(Psat(i)-Psat(i-1))*(P(j)-Psat(i))+Cpl_R245FA(i);
            Cpg=(Cpg_R245FA(i)-Cpg_R245FA(i-1))/(Psat(i)-Psat(i-1))*(P(j)-Psat(i))+Cpg_R245FA(i);
            Tf=(Tsat(i)-Tsat(i-1))/(Psat(i)-Psat(i-1))*(P(j)-Psat(i))+Tsat(i);
            hlg=Hg-Hl;
        end
    end
    if x(j)<lmt
        
        Pr=Cpl*mul/lambda;

        X=frottement_par(Jll(j),Jgg(j),rhol,rhog,mul,mug,D);

        hl=lambda/D*0.023.*(G.*(1.-x(j))*D/mul).^0.8*Pr^(1/3);

        h(j)=7.39e3.*hl.*(qq./(G*hlg)+0.00015.*(1./X).^0.66);
        
        Tp(j)=Tf+qq./(h(j));        
        Tc(j)=Tp(j)+1/Conductance(h(j))*Q_geom(1)/(L_geom(1)*l_geom(1))*pi*D*3*Np_geom(1)/l_geom(1)+1/R_contact*Q_geom(1)/(L_geom(1)*l_geom(1));
        if qq==0
            h(j)=0;
            Tc(j)=Tf;
        end
    else
        Pr=Cpl*mul/lambda;
        hl=0.023*lambda/D*(G*D/mul)^0.8*Pr^(1/3);
        h(j)=hl;
        Tp(j)=Tf+qq*(4/(D*G*Cpg)*(Zz(j)-Zz(1))+1/hl);   
        Tc(j)=Tp(j)+1/Conductance(h(j))*Q_geom(1)/(L_geom(1)*l_geom(1))*pi*D*3*Np_geom(1)/l_geom(1)+1/R_contact*Q_geom(1)/(L_geom(1)*l_geom(1));
        if qq==0
            h(j)=0;
            Tc(j)=Tf;
        end
    end
end


for k=2:length(N)
    g=ggeom(k);
    qq=qgeom(k);
    for j=sum(N(1:k-1))+1:sum(N(1:k))
        for i=2:length(Psat)
            if (P(j)<=Psat(i)) && (P(j)>Psat(i-1))
                rhol=(rhol_R245FA(i)-rhol_R245FA(i-1))/(Psat(i)-Psat(i-1))*(P(j)-Psat(i))+rhol_R245FA(i);
                rhog=(rhog_R245FA(i)-rhog_R245FA(i-1))/(Psat(i)-Psat(i-1))*(P(j)-Psat(i))+rhog_R245FA(i);
                mul=(mul_R245FA(i)-mul_R245FA(i-1))/(Psat(i)-Psat(i-1))*(P(j)-Psat(i))+mul_R245FA(i);
                mug=(mug_R245FA(i)-mug_R245FA(i-1))/(Psat(i)-Psat(i-1))*(P(j)-Psat(i))+mug_R245FA(i);
                Hl=(Hl_R245FA(i)-Hl_R245FA(i-1))/(Psat(i)-Psat(i-1))*(P(j)-Psat(i))+Hl_R245FA(i);
                Hg=(Hg_R245FA(i)-Hg_R245FA(i-1))/(Psat(i)-Psat(i-1))*(P(j)-Psat(i))+Hg_R245FA(i);
                lambda=(lambdal_R245FA(i)-lambdal_R245FA(i-1))/(Psat(i)-Psat(i-1))*(P(j)-Psat(i))+lambdal_R245FA(i);
                Cpl=(Cpl_R245FA(i)-Cpl_R245FA(i-1))/(Psat(i)-Psat(i-1))*(P(j)-Psat(i))+Cpl_R245FA(i);
                Cpg=(Cpg_R245FA(i)-Cpg_R245FA(i-1))/(Psat(i)-Psat(i-1))*(P(j)-Psat(i))+Cpg_R245FA(i);
                Tf=(Tsat(i)-Tsat(i-1))/(Psat(i)-Psat(i-1))*(P(j)-Psat(i))+Tsat(i);
                hlg=Hg-Hl;
            end
        end
        if x(j)<lmt

            Pr=Cpl*mul/lambda;

            X=frottement_par(Jll(j),Jgg(j),rhol,rhog,mul,mug,D);

            hl=lambda/D*0.023.*(G.*(1.-x(j))*D/mul).^0.8*Pr^(1/3);

            h(j)=7.39e3.*hl.*(qq./(G*hlg)+0.00015.*(1./X).^0.66);

            Tp(j)=Tf+qq./(h(j));
            Tc(j)=Tp(j)+1/Conductance(h(j))*Q_geom(1)/(L_geom(k)*l_geom(k))*pi*D*3*Np_geom(k)/l_geom(k)+1/R_contact*Q_geom(k)/(L_geom(k)*l_geom(k));
            if qq==0
                h(j)=0;
                Tc(j)=Tf;
            end
        else
            Pr=Cpl*mul/lambda;
            hl=0.023*lambda/D*(G*D/mul)^0.8*Pr^(1/3);
            h(j)=hl;
            Tp(j)=Tf+qq*(4/(D*G*Cpg)*(Zz(j)-Zz(sum(N(1:k-1))))+1/hl);
            Tc(j)=Tp(j)+1/Conductance(h(j))*Q_geom(1)/(L_geom(k)*l_geom(k))*pi*D*3*Np_geom(k)/l_geom(k)+1/R_contact*Q_geom(k)/(L_geom(k)*l_geom(k));
            if qq==0
                h(j)=0;
                Tc(j)=Tf;
            end
        end
    end
end