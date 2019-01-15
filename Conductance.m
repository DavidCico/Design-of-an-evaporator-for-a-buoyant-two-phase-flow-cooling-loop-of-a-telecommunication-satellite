%% Function to calculate the fin thermal conductance

function R_conduct=Conductance(h)

R2=4756.30728;
R1=4163.084707;
h2=1300;
h1=500;

R_conduct=(R2-R1)/(h2-h1)*(h-h1)+R1;