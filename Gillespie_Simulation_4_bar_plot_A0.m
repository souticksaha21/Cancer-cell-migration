%Plot of mean and median of CI for Delta TGF for different values of A0
clear all
hold off

trials=400;
lg = [.75 1 .75];
gr = [0 .75 0];
gy = [.5 .5 .5];
or = [1 0.6 0.5];
%A0=1000
fCI1 = fopen('../data/Gillespie_Simulation_4_CI_T_A0_1000.dat');
t1 = textscan(fCI1, '%f%f%f%f%f');
fclose(fCI1);
CI1=t1{3};
CIT=median(CI1); %median of CI
CIT_avg=mean(CI1); %mean CI
CIT_ste=std(CI1)/sqrt(trials);%standard error
CIT_stme=median(abs(CI1-median(CI1)*ones(400,1)))/sqrt(trials); %standard error of median

%A0=3000
fCI2 = fopen('../data/Gillespie_Simulation_4_CI_T_A0_3000.dat');
t2 = textscan(fCI2, '%f%f%f%f%f');
fclose(fCI2);
CI2=t2{3};
CIE=median(CI2);
CIE_avg=mean(CI2);
CIE_ste=std(CI2)/sqrt(trials);
CIE_stme=median(abs(CI2-median(CI2)*ones(400,1)))/sqrt(trials);

%A0=5000
fCI3 = fopen('../data/Gillespie_Simulation_4_CI_T_A0_5000.dat');
t3 = textscan(fCI3, '%f%f%f%f%f');
fclose(fCI3);
CI3=t3{3};
CITE=median(CI3);
CITE_avg=mean(CI3);
CITE_ste=std(CI3)/sqrt(trials);
CITE_stme=median(abs(CI3-median(CI3)*ones(400,1)))/sqrt(trials);

%A0=7000
fCI4 = fopen('../data/Gillespie_Simulation_4_CI_T_A0_7000.dat');
t4 = textscan(fCI4, '%f%f%f%f%f');
fclose(fCI4);
CI4=t4{3};
CITE0=median(CI4);
CITE0_avg=mean(CI4);
CITE0_ste=std(CI4)/sqrt(trials);
CITE0_stme=median(abs(CI4-median(CI4)*ones(400,1)))/sqrt(trials);

%A0=9000
fCI5 = fopen('../data/Gillespie_Simulation_4_CI_T_A0_9000.dat');
t5 = textscan(fCI5, '%f%f%f%f%f');
fclose(fCI5);
CI5=t5{3};
CITT0=median(CI5);
CITT0_avg=mean(CI5);
CITT0_ste=std(CI5)/sqrt(trials);
CITT0_stme=median(abs(CI5-median(CI5)*ones(400,1)))/sqrt(trials);

%A0=25000
fCI6 = fopen('../data/Gillespie_Simulation_4_CI_T_A0_25000.dat');
t6 = textscan(fCI6, '%f%f%f%f%f');
fclose(fCI6);
CI6=t6{3};
CIM6=median(CI6);
CIA6=mean(CI6);
CIA6_ste=std(CI6)/sqrt(trials);
CIM6_stme=median(abs(CI6-median(CI6)*ones(size(CI6,1),1)))/sqrt(size(CI6,1));

%A0=50000
fCI7 = fopen('../data/Gillespie_Simulation_4_CI_T_A0_50000.dat');
t7 = textscan(fCI7, '%f%f%f%f%f');
fclose(fCI7);
CI7=t7{3};
CIM7=median(CI7);
CIA7=mean(CI7);
CIA7_ste=std(CI7)/sqrt(trials);
CIM7_stme=median(abs(CI7-median(CI7)*ones(size(CI7,1),1)))/sqrt(size(CI7,1));

%A0=100000
fCI11 = fopen('../data/Gillespie_Simulation_4_CI_T_A0_100000.dat');
t11 = textscan(fCI11, '%f%f%f%f%f');
fclose(fCI11);
CI11=t11{3};
CIT11=median(CI11);
CIT11_avg=mean(CI11);
CIT11_ste=std(CI11)/sqrt(trials);
CIT11_stme=median(abs(CI11-median(CI11)*ones(400,1)))/sqrt(trials);

%A0=300000
fCI12 = fopen('../data/Gillespie_Simulation_4_CI_T_A0_300000.dat');
t12 = textscan(fCI12, '%f%f%f%f%f');
fclose(fCI12);
CI12=t12{3};
CIT12=median(CI12);
CIT12_avg=mean(CI12);
CIT12_ste=std(CI12)/sqrt(trials);
CIT12_stme=median(abs(CI12-median(CI12)*ones(400,1)))/sqrt(trials);

%A0=500000
fCI13 = fopen('../data/Gillespie_Simulation_4_CI_T_A0_500000.dat');
t13 = textscan(fCI13, '%f%f%f%f%f');
fclose(fCI13);
CI13=t13{3};
CIT13=median(CI13);
CIT13_avg=mean(CI13);
CIT13_ste=std(CI13)/sqrt(trials);
CIT13_stme=median(abs(CI13-median(CI13)*ones(400,1)))/sqrt(trials);

%A0=700000
fCI14 = fopen('../data/Gillespie_Simulation_4_CI_T_A0_700000.dat');
t14 = textscan(fCI14, '%f%f%f%f%f');
fclose(fCI14);
CI14=t14{3};
CIT14=median(CI14);
CIT14_avg=mean(CI14);
CIT14_ste=std(CI14)/sqrt(trials);
CIT14_stme=median(abs(CI14-median(CI14)*ones(400,1)))/sqrt(trials);

k22=1.48;


hold off
figure(4); clf;
x = [1000,3000,5000,7000,9000,25000,50000,100000,300000,500000,700000];
prices = [CIT CIE CITE CITE0 CITT0 CIM6 CIM7 CIT11 CIT12 CIT13 CIT14];
stdme = 0.0008*sqrt(x);
fig4 = loglog(x,prices,'-ko',...
    'LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',or,...
    'MarkerSize',10); hold on;
%h2=loglog(x,stdme,'-b','LineWidth',2);
%errorbar(x,prices,stde,'.k','LineWidth',1.5);
ylim([0.1 1.2])
xticks([1000,3000,5000,7000,9000,25000,50000,100000,300000,500000,700000])
yticks([0.1,0.3,0.5,0.7,0.9,1.1])
%legend([fig4 h2],'CI median','fit','Location','northwest');
ylabel('CI median');
xlabel('M_0');
title('CI median')
set(gca,'FontSize',13)
xticklabels({'1000','3000','5000','7000','9000','25000','50000','100000','300000','500000','700000'})
xtickangle(45);
saveas(fig4,'CI_median_fig_1_final.png');

hold off
figure(5); clf;
x = [1000,3000,5000,7000,9000,25000,50000,100000,300000,500000,700000];
prices = [CIT_avg CIE_avg CITE_avg CITE0_avg CITT0_avg CIA6 CIA7 CIT11_avg CIT12_avg CIT13_avg CIT14_avg];
stde = 0.0008*sqrt(x);
fig5 = loglog(x,prices,'-ko',...
    'LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',lg,...
    'MarkerSize',10); hold on;
%h2=loglog(x,stde,'-b','LineWidth',2);
%errorbar(x,prices,stde,'.k','LineWidth',1.5);
ylim([0.1 1.2])
xticks([1000,3000,5000,7000,9000,25000,50000,100000,300000,500000,700000])
yticks([0.1,0.3,0.5,0.7,0.9,1.1])
%legend([fig5 h2],'CI mean','fit','Location','northwest');
ylabel('CI average');
xlabel('M_0');
title('CI mean')
set(gca,'FontSize',13)
xticklabels({'1000','3000','5000','7000','9000','25000','50000','100000','300000','500000','700000'})
xtickangle(45);
saveas(fig5,'CI_avg_fig_1_final.png');

hold off
figure(5); clf;
x = [1000,3000,5000,7000,9000,25000,50000,100000,300000,500000,700000];
prices = [CIT CIE CITE CITE0 CITT0 CIM6 CIM7 CIT11 CIT12 CIT13 CIT14];
stdme = 0.0008*sqrt(x);
fig4 = loglog(x,prices,'-ko',...
    'LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',or,...
    'MarkerSize',10); hold on;
%h2=loglog(x,stdme,'-b','LineWidth',2);
%errorbar(x,prices,stde,'.k','LineWidth',1.5);
ylim([0.1 1.2])
xticks([1000,5000,50000,500000])
yticks([0.1,0.3,0.5,0.7,0.9,1.1])
%legend([fig4 h2],'CI median','fit','Location','northwest');
ylabel('CI median');
xlabel('M_0');
%title('CI median')
set(gca,'FontSize',15)
xticklabels({'1000','5000','50000','500000'})
xtickangle(45);
saveas(fig4,'CI_median_fig_nolte_final.png');

