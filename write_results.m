%% write results in a txt file

% with gravity

fID1=fopen('results_with_gravity.txt','w+');

fprintf(fID1,'%c','position z');            fprintf(fID1, '%c','    ');
fprintf(fID1,'%c','title x');           fprintf(fID1, '%c','    ');
fprintf(fID1,'%c','Rg');                fprintf(fID1, '%c','      ');
fprintf(fID1,'%c','pressure');          fprintf(fID1, '%c','    ');
fprintf(fID1,'%c','Tc G&W');            fprintf(fID1, '%c','      ');
fprintf(fID1,'%c','Tc K');              fprintf(fID1, '%c','      ');
fprintf(fID1,'%c','Tc S&G');            fprintf(fID1, '%c','     ');

fprintf(fID1, '%c','     ');            fprintf(fID1,'%c','gravity not null');


fprintf(fID1,'%c\n','');                fprintf(fID1, '%c','                                                                                '); 
fprintf(fID1,'%c','Delta T = ');        fprintf(fID1,'%2.3f',Delta_T);      fprintf(fID1,'%c','°C');
fprintf(fID1,'%c\n','');

for i=1:length(x_g)
    fprintf(fID1,'%2.3f',Z_g(i));       fprintf(fID1, '%c','     ');
    fprintf(fID1,'%2.3f',x_g(i));       fprintf(fID1, '%c','     ');
    fprintf(fID1,'%2.3f',Rg_g(i));      fprintf(fID1, '%c','     ');
    fprintf(fID1,'%2.0f',P_g(i));       fprintf(fID1, '%c','     ');
    fprintf(fID1,'%2.3f',Tc_G_W(i));    fprintf(fID1, '%c','     ');
    fprintf(fID1,'%2.3f',Tc_K(i));      fprintf(fID1, '%c','     ');
    fprintf(fID1,'%2.3f',Tc_S_G(i));    fprintf(fID1, '%c','     ');
    
    
    fprintf(fID1,'%c\n','');
end

fclose(fID1);