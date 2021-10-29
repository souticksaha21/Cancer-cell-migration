%This code compares the CI between 2D and 3D migration
clear all;
clc;
%hold off;
pi=3.14159; %defining pi
trials=2500; %number of trajectories for each value of bias (epsilon)
for i=1:5 %loop over number of bias
    bias=i*0.02; %deifining bias for each loop
for tr=1:trials %loop over trajectories
%initialising position of a cell in 2D
x2=0; 
y2=0;
for n = 1:1000 %loop over time points    
p=0;
x=1;
%choosing run direction
while p < x
    x=rand;
    t=2*pi*rand;
    p=(1+bias*cos(t))/(2*pi);
end
%updating position
x2=x2+cos(t);
y2=y2+sin(t);
end
%defining CI for a trajectory
CI2=x2/(x2^2+y2^2)^0.5;
%storing CI into an array
A2(tr)=CI2;

%code for 3D migration
x3=0;
y3=0;
z3=0;
for n = 1:1000
p=0;
x=1;
while p < x
    x=rand;
    t=2*pi*rand;
    p=(1+bias*cos(t))/(2*pi);
end
k5=rand;
x3=x3+sin(pi*k5)*cos(t);
y3=y3+sin(pi*k5)*sin(t);
z3=z3+cos(pi*k5);
end
CI3=x3/(x3^2+y3^2)^0.5;
A3(tr)=CI3;
end

%defininf moments of CI
mean2=mean(A2);
std2=std(A2)/sqrt(trials);
CImean2(i)=mean2;
Xmean2(i)=bias;
CImeanstd2(i)=std2;

median2=median(A2);
stme2=median(abs(A2-median(A2)*ones(1,trials)))/sqrt(trials);
CImedian2(i)=median2;
Xmedian2(i)=bias;
CImedianstd2(i)=stme2;

mean3=mean(A3);
std3=std(A3)/sqrt(trials);
CImean3(i)=mean3;
Xmean3(i)=bias+0.006;
CImeanstd3(i)=std3;

median3=median(A3);
stme3=median(abs(A3-median(A3)*ones(1,trials)))/sqrt(trials);
CImedian3(i)=median3;
Xmedian3(i)=bias+0.006;
CImedianstd3(i)=stme3;

%hold off

end

lg = [.75 1 .75];
gr = [0 .75 0];
gy = [.5 .5 .5];
or = [1 0.6 0.5];

hold off
figure(2); clf;
%x = [2 3 4 5];
%prices = [mean2 mean3 median2 median3];
%stde = [std2 std3 stme2 stme3];
fig1 = plot(Xmean2,CImean2,'-ko',...
    'LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',or,...
    'MarkerSize',10); hold on;
h1=plot(Xmean2,CImean3,'-ko',...
    'LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',lg,...
    'MarkerSize',10); hold on;
ylim([-0.09 1.2])
%xlim([0 10])
xticks([0.03 0.05 0.07 0.09 0.11])
ylabel('CI');
xlabel('bias');
legend([fig1 h1],'2D','3D');
set(gca,'FontSize',13)
%xticklabels({'mean 2D','mean 3D','median 2D','median 3D'})
%xtickangle(45);
xlabel('\epsilon');
title('CI Mean comparison')
saveas(fig1,'CI_2D_vs_3D_mean.png');

hold off
figure(3); clf;
%x = [2 3 4 5];
%prices = [mean2 mean3 median2 median3];
%stde = [std2 std3 stme2 stme3];
fig3 = plot(Xmedian2,CImedian2,'-ko',...
    'LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',or,...
    'MarkerSize',10); hold on;
h3=plot(Xmedian2,CImedian3,'-ko',...
    'LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',lg,...
    'MarkerSize',10); hold on;
ylim([-0.09 1.2])
%xlim([0 10])
xticks([0.03 0.05 0.07 0.09 0.11])
ylabel('CI');
xlabel('\epsilon');
legend([fig3 h3],'2D','3D');
set(gca,'FontSize',13)
title('CI Median comparison')
%xticklabels({'mean 2D','mean 3D','median 2D','median 3D'})
%xtickangle(45);
saveas(fig3,'CI_2D_vs_3D_median.png');