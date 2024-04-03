clear all;clc;
numerator = 0;
denominator = 0;
numerator1 = 0;
denominator1 = 0;
denominator2 = 0;
SNR = [];
RSNR = [];
time = [];
num = 100;
for k = 1:num
    load('SNR_new/50/A.mat'); % 
    load(['SNR_new/50/b',num2str(k),'.mat'],'b'); % 
    load(['SNR_new/50/x',num2str(k),'.mat'],'x_true'); % 
    SNR = [SNR;10*log10(norm(A*x_true,2)^2/norm(b-A*x_true,2)^2)];
    numerator = numerator + norm(A*x_true,2)^2;
    denominator = denominator + norm(b - A*x_true,2)^2;
    
    lambda = 2;
    out = sparse_l0(A,lambda,b,1); %1.8-4

    indt = find(x_true>0);
    indx = find(out.x>0);
    length(intersect(indt,indx))
    RSNR = [RSNR;10*log10(norm(x_true)^2/norm(out.x-x_true)^2)];
    time = [time;out.time];
    denominator1 = denominator1 + norm(out.x-x_true)^2;
end

SNR_mean = 10*log10((numerator)/(denominator));
RSNR_mean = 10*log10((numerator1)/(denominator1));


% rowNames = {'1','2','3','4','5','6','7','8','9','10','mean'};
rowNames = cell(1, num+1);  
for i = 1:num
    rowNames{i} = num2str(i); 
end
rowNames{num+1} = 'mean';
result = [SNR SNR_mean*ones(num,1) RSNR RSNR_mean*ones(num,1) time 2*ones(num,1)];
result = [result;mean(result)];
colNames = {'SNR','SNR_mean','RSNR','RSNR_mean','time','lambda'};
result = array2table(result,'RowNames',rowNames,'VariableNames',colNames)
% writetable(result,'SNR_new/SNR50.txt') % 
%%
clear all;clc
path = 'SNR_new/';
namelist = dir([path,'*.txt']);
RSNR = []; time = []; RSNR_his = []; time_his = [];
for i = 1:length(namelist)
    filename = [path,namelist(i).name];
    result = readtable(filename,VariableNamingRule='preserve');
    RSNR = [RSNR result{1:100,'RSNR(our)'}];
    time = [time result{1:100,"time(our)"}];
    % RSNR_his = [RSNR_his result{1:100,"RSNR(his)"}];
    % time_his = [time_his result{1:100,"time(his)"}];
end
mean_data_o = mean(RSNR);
std_data_o = std(RSNR);
% mean_data_h = mean(RSNR_his);
% std_data_h = std(RSNR_his);

mean_time_o = mean(time);
std_time_o = std(time);
% mean_time_h = mean(time_his);
% std_time_h = std(time_his);
% SNR15 = readtable('SNR_new/SNR15.txt')
%%
set(gcf, 'PaperSize', [25 25])
% set(fig,'defaultAxesColorOrder',['black'; 'black']);

e1 = errorbar(15:5:50,mean_data_o, std_data_o,'-s',"MarkerSize",10, ...
    "MarkerEdgeColor","black","MarkerFaceColor",'#219ebc','CapSize',12,'LineWidth',1); hold on;
p1 = plot(15:5:50,mean_data_o,'-s',"MarkerSize",10, ...
    "MarkerEdgeColor","black","MarkerFaceColor",'#219ebc','LineWidth',1, 'Color','black');

% e2 = errorbar(15:5:50,mean_data_h, std_data_h,'-o',"MarkerSize",10, ...
%     "MarkerEdgeColor","black","MarkerFaceColor",'#ffb703','CapSize',12,'LineWidth',1);
% p2 = plot(15:5:50,mean_data_h,'-o',"MarkerSize",10, ...
%     "MarkerEdgeColor","black","MarkerFaceColor",'#ffb703','LineWidth',1, 'Color','black');
set(gca,'FontSize',15);
xlim([13,52]);
grid on;
xlabel('SNR (dB)','Fontname','Times New Roman')
ylabel("RSNR (dB)",'Fontname','Times New Roman')


yyaxis right
set(gca,'YColor','red');
p3 = plot(15:5:50,mean_time_o,'-s',"MarkerSize",10, ...
    "MarkerEdgeColor","black","MarkerFaceColor",'#219ebc','LineWidth',1.5,'Color','red'); hold on;
% p4 = plot(15:5:50,mean_time_h,'--o',"MarkerSize",10, ...
%     "MarkerEdgeColor","black","MarkerFaceColor",'#ffb703','LineWidth',1.5,'Color','red');
ylim([-0.5,25])
ylabel('Time (s)','Fontname','Times New Roman')
% legend('time(our)','time(GPG)')
legend({'RSNR (Alg.2)','','RSNR (GPG)','','time (Alg.2)','time (GPG)'},'Location','northwest','NumColumns',2,'Fontname','Times New Roman')
hold off;