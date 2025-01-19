%% demo of the experiment 3 (efficacy of the L0BPG step in Algorithm 2)
clear all;clc
lambda = 1.5;
load('data/Experiment3/A6.mat');
load('data/Experiment3/b6.mat');
load('data/Experiment3/x_true.mat');

% keeping the 12th element
out = sparse_l0_effect_first12(A,lambda, b,1);
norm(A*out.x-b)^2/2

% keeping the 14th element
out = sparse_l0_effect(A,lambda, b,1);
norm(A*out.x-b)^2/2

indt = find(x_true>0);
indo = find(out.x>0);
length(intersect(indt,indo))
track = out.track;

% plot(track(3,:));hold on;
plot(track(4,:));hold on;
plot(track(5,:));hold on;
plot(track(6,:));hold on;
plot(track(7,:));hold on;
legend('1','2','3','4')
%% plot
trackn = track;
for i = 4:7
    for j = 1:811
        if trackn(i,j) == 0
            trackn(i,j) = 0.0145;
        end
    end
end
len = length(trackn(4,:));
set(gcf, 'PaperSize', [25 25])
plot(300:len,trackn(4,300:end),'-o','MarkerSize',10,'Linewidth', 3,'MarkerIndices',1:50:800);hold on; 
plot(300:len,trackn(5,300:end),'-*','MarkerSize',10,'Linewidth', 3,'MarkerIndices',1:50:800);hold on;
plot(300:len,trackn(6,300:end),'-s','MarkerSize',10,'Linewidth', 3,'MarkerIndices',1:50:800);hold on;
plot(300:len,trackn(7,300:end),'-x','MarkerSize',10,'Linewidth', 3,'MarkerIndices',1:50:800);hold on;
line([354 354],[0.0145 0.02],'LineStyle','--','LineWidth',3,'Color','black'); hold on;
line([300 850], [0.015224 0.015224],'LineStyle','--','LineWidth',3,'Color','black'); hold on;
grid on
set(gca,'FontSize',18);
xlabel('Iteration','Fontname','Times New Roman','FontSize',25)
ylabel("Value",'Fontname','Times New Roman','FontSize',25)
legend({'12th largest','13th largest','14th largest','15th largest'}, 'NumColumns',1,...
    'location','northeast' ,'Fontname','Times New Roman','FontSize',20)

ylim([0.0145,0.0195])
yticks([0.0145 0.015,0.016,0.017,0.018,0.019,0.02])
yticklabels({'0','0.015','0.016','0.017','0.018','0.019','0.02'})

xlim([300,850])
title('Numerical change','FontSize',25)
