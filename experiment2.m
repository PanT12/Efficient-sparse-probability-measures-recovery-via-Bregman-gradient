%% demo of the experiment 2 (support accuracy)
clear all;clc;

% save result
mkdir("result/Experiment2")

num = 100;

% load data
A_sizes_300 = [30,33,35,37,40,43,45,47,50,60,70];
A_sizes_900 = [90,100,110,120,130,140,150,170,180,210];
n_size = [300, 900];

for ni = 1:length(n_size)
    n = n_size(ni)
    path = ['result/Experiment2/', num2str(n), 'result'];
    mkdir(path)
    if ni == 1
        A_sizes = A_sizes_300;
        lambdas = [14,13,12,12,10,9,9,7,5,1.5,0.5]./10;
    elseif ni == 2
        A_sizes = A_sizes_900;
        lambdas = [11,12,11,9,6,4,2.5,1,1,0.5]./10;
    end
    for si = 1:length(A_sizes)
        m = A_sizes(si);
        time = zeros(num,1);
        norm1 = zeros(num,1);
        actualtrue = zeros(num,1);
        predtrue = zeros(num,1);
        TP = zeros(num,1);
        TN = zeros(num,1);
        FN = zeros(num,1);
        FP = zeros(num,1);
        accuracy = zeros(num,1);
        precision = zeros(num,1);
        recall = zeros(num,1);
        F1 = zeros(num,1);
        load(['data/Experiment2/',num2str(m),'*', num2str(n),'/A.mat'], 'A'); % 
        lambda = lambdas(si);
        for k = 1:num
            load(['data/Experiment2/',num2str(m),'*', num2str(n),'/b',num2str(k),'.mat'],'b'); % 
            load(['data/Experiment2/',num2str(m),'*', num2str(n),'/x',num2str(k),'.mat'],'x_true'); % 
            
            indt = find(x_true>0);
            IZ1 = find(x_true == 0);
            actualtrue(k)=length(indt);
        
            out = sparse_l0(A,lambda,b,1,1e-7,0); 
            indx = find(out.x>0);
            IZ = find(out.x == 0);
            predtrue(k) = length(indx);
            TP(k) = length(intersect(indt,indx));
            time(k) = out.time;
            norm1(k) = norm(A*out.x-b,2)^2/2;   
            TN(k) = length(intersect(IZ,IZ1));
            FN(k) = length(intersect(indt,IZ));
            FP(k) = n-TP(k)-TN(k)-FN(k);
            accuracy(k) = (TP(k)+TN(k))/n;
            precision(k) = TP(k)/(TP(k)+FP(k));
            recall(k) = TP(k)/(TP(k)+FN(k));
            F1(k) = 2*precision(k)*recall(k)/(precision(k)+recall(k));
        
        end
        % rowNames = {'1','2','3','4','5','6','7','8','9','10','mean'};
        rowNames = cell(1, num);  
        for i = 1:num
            rowNames{i} = num2str(i); 
        end
        result = [time norm1 actualtrue predtrue TP FP TN FN accuracy precision recall F1 lambda*ones(num,1)];
        colNames = {'time','norm','actualtrue','predtrue','TP','FP','TN','FN','accuracy','precision','recall','F1','lambda'};
        result = array2table(result,'RowNames',rowNames,'VariableNames',colNames)
        mean(result)
        writetable(result,[path,'/',num2str(m)]) % 
        % writetable(result,'confusion matrix_new/900result/210900.txt') % 
    end
end
%% 300 size
% showing result
clear all;clc
path = 'result/Experiment2/300result/';
namelist = dir([path,'*.txt']);
[~,ind] = sort([namelist(:).datenum],'ascend');
l = length(namelist);
namelist = namelist(ind);
num = 100;
F1 = zeros(num,l);
accuracy = zeros(num,l);
precision = zeros(num,l);
recall = zeros(num,l);
time = zeros(num,l);
diff = zeros(num,l);

for i = 1:l
    filename = [path,namelist(i).name]
    result = readtable(filename);
    accuracy(:,i) = result{1:num,'accuracy'};
    precision(:,i) = result{1:num,'precision'};
    recall(:,i) = result{1:num,'recall'};
    F1(:,i) = result{1:num,'F1'};
    time(:,i) = result{1:num,"time"};
    diff(:,i) = result{1:num, "norm"};
end
%%
set(gcf, 'PaperSize', [25 25])
x = [30,33,35,37,40,43,45,47,50,60,70];
% errorbar(x,mean(accuracy), std(accuracy),'-s',"MarkerSize",10, ...
%     "MarkerEdgeColor","black","MarkerFaceColor",'#219ebc','CapSize',12,'LineWidth',1,'Color','#219ebc'); hold on;
% errorbar(x,mean(precision), std(precision),'-s',"MarkerSize",10, ...
%     "MarkerEdgeColor","black","MarkerFaceColor",'#ffb703','CapSize',12,'LineWidth',1,'Color','#ffb703'); hold on;
% errorbar(x,mean(recall), std(recall),'-s',"MarkerSize",10, ...
%     "MarkerEdgeColor","black","MarkerFaceColor",'red','CapSize',12,'LineWidth',1,'Color','red'); hold on;
% errorbar(x,mean(F1), std(F1),'-s',"MarkerSize",10, ...
%     "MarkerEdgeColor","black","MarkerFaceColor",'#9c3587','CapSize',12,'LineWidth',1,'Color','#9c3587'); hold on;

plot(x,mean(accuracy),'-s',"MarkerSize",10, ...
    "MarkerEdgeColor","black","MarkerFaceColor",'#219ebc','LineWidth',1,'Color','#219ebc'); hold on;
plot(x,mean(precision),'-s',"MarkerSize",10, ...
    "MarkerEdgeColor","black","MarkerFaceColor",'#ffb703','LineWidth',1,'Color','#ffb703'); hold on;
plot(x,mean(recall),'-s',"MarkerSize",10, ...
    "MarkerEdgeColor","black","MarkerFaceColor",'red','LineWidth',1,'Color','red'); hold on;
plot(x,mean(F1),'-s',"MarkerSize",10, ...
    "MarkerEdgeColor","black","MarkerFaceColor",'#9c3587','LineWidth',1,'Color','#9c3587'); hold on;

set(gca,'FontSize',18);
xlim([28,72]);
set(gca,'xtick',30:5:70) 
grid on;
xlabel('The number of row in A','Fontname','Times New Roman','FontSize',25)
ylabel("Value",'Fontname','Times New Roman','FontSize',25)

% legend('time(our)','time(GPG)')
legend({'accuracy','precision','recall','F1'},'Location','southeast','NumColumns',1,'Fontname','Times New Roman','FontSize',18)
hold off;

%% 900 size
% showing result
clear all;clc
path = 'result/Experiment2/900result/';
namelist = dir([path,'*.txt']);
[~,ind] = sort([namelist(:).datenum],'ascend');
l = length(namelist);
namelist = namelist(ind);
num = 100;
F1 = zeros(num,l);
accuracy = zeros(num,l);
precision = zeros(num,l);
recall = zeros(num,l);
time = zeros(num,l);
diff = zeros(num,l);

for i = 1:l
    filename = [path,namelist(i).name]
    result = readtable(filename);
    accuracy(:,i) = result{1:num,'accuracy'};
    precision(:,i) = result{1:num,'precision'};
    recall(:,i) = result{1:num,'recall'};
    F1(:,i) = result{1:num,'F1'};
    time(:,i) = result{1:num,"time"};
    diff(:,i) = result{1:num, "norm"};
end
%%
set(gcf, 'PaperSize', [25 25])
x = [90,100,110,120,130,140,150,170,180,210];
% errorbar(x,mean(accuracy), std(accuracy),'-s',"MarkerSize",10, ...
%     "MarkerEdgeColor","black","MarkerFaceColor",'#219ebc','CapSize',12,'LineWidth',1,'Color','#219ebc'); hold on;
% errorbar(x,mean(precision), std(precision),'-s',"MarkerSize",10, ...
%     "MarkerEdgeColor","black","MarkerFaceColor",'#ffb703','CapSize',12,'LineWidth',1,'Color','#ffb703'); hold on;
% errorbar(x,mean(recall), std(recall),'-s',"MarkerSize",10, ...
%     "MarkerEdgeColor","black","MarkerFaceColor",'red','CapSize',12,'LineWidth',1,'Color','red'); hold on;
% errorbar(x,mean(F1), std(F1),'-s',"MarkerSize",10, ...
%     "MarkerEdgeColor","black","MarkerFaceColor",'#9c3587','CapSize',12,'LineWidth',1,'Color','#9c3587'); hold on;

plot(x,mean(accuracy),'-s',"MarkerSize",10, ...
    "MarkerEdgeColor","black","MarkerFaceColor",'#219ebc','LineWidth',1,'Color','#219ebc'); hold on;
plot(x,mean(precision),'-s',"MarkerSize",10, ...
    "MarkerEdgeColor","black","MarkerFaceColor",'#ffb703','LineWidth',1,'Color','#ffb703'); hold on;
plot(x,mean(recall),'-s',"MarkerSize",10, ...
    "MarkerEdgeColor","black","MarkerFaceColor",'red','LineWidth',1,'Color','red'); hold on;
plot(x,mean(F1),'-s',"MarkerSize",10, ...
    "MarkerEdgeColor","black","MarkerFaceColor",'#9c3587','LineWidth',1,'Color','#9c3587'); hold on;

set(gca,'FontSize',18);
xlim([85,215]);
grid on;
xlabel('The number of row in A','Fontname','Times New Roman','FontSize',25)
ylabel("Value",'Fontname','Times New Roman','FontSize',25)

% legend('time(our)','time(GPG)')
legend({'accuracy','precision','recall','F1'},'Location','southeast','NumColumns',1,'Fontname','Times New Roman','FontSize',18)
hold off;
