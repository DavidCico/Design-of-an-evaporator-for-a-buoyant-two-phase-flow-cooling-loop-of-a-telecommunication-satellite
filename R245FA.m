global Tsat Psat rhol_R245FA rhog_R245FA Hl_R245FA Hg_R245FA Cpl_R245FA Cpg_R245FA mul_R245FA mug_R245FA lambdal_R245FA lambdag_R245FA sigma_R245FA

data=importdata('R123-R245FA.xls');

Tsat=data.data.LiquidR245FA(75:180,1);
Psat=data.data.LiquidR245FA(75:180,2)*1e6;

% liquid
rhol_R245FA=data.data.LiquidR245FA(75:134,3);
Hl_R245FA=data.data.LiquidR245FA(75:134,6)*1e3;
Cpl_R245FA=data.data.LiquidR245FA(75:134,9)*1e3;
mul_R245FA=data.data.LiquidR245FA(75:134,12);
lambdal_R245FA=data.data.LiquidR245FA(75:134,13);
sigma_R245FA=data.data.LiquidR245FA(75:134,14);

% vapor
rhog_R245FA=data.data.VaporR245FA(75:134,3);
Hg_R245FA=data.data.VaporR245FA(75:134,6)*1e3;
Cpg_R245FA=data.data.VaporR245FA(75:134,9)*1e3;
mug_R245FA=data.data.VaporR245FA(75:134,12);
lambdag_R245FA=data.data.VaporR245FA(75:134,13);

% plot(Tsat,Psat)