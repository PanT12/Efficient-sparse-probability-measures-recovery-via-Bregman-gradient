%% demo of the experiment 6 (portfolio optimization)
% You can download the dataset in http://people.brunel.ac.uk/~mastjjb/jeb/orlib/portinfo.html
%% standard efficient frontier
clear all;clc;

% save data
save_path = 'result/Experiment6';
mkdir(save_path)

% load data
path = 'data/Experiment6/';


data = [[1,31];[2,85];[3,89];[4,98];[5,225]];
len = 2000;

for j = 1:size(data,1)
    k = data(j,1);
    num = data(j,2);
    filename = [path, 'port',num2str(k),'.txt'];
    [number,mean,deviance,correlation] = importfile(filename,[1,(1+num)*num/2+1+num]); %
    cov = diag(deviance)*correlation*diag(deviance);
    A = cov;
    % A = cov + eye(number)/(100/sqrt(number));
    
    time = zeros(1,len);
    x = zeros(number,len);
    variance = zeros(len,1);
    mean_val = zeros(len,1);
    lambda = zeros(len,1);
    for i = 1:len
        mu = (i)/len
        lambda(i) = mu;
        out = sparse_portfolio(mu*A,0,(1-mu)*mean,1,number);
        variance(i) = out.x'*A*out.x;
        mean_val(i) = out.x'*mean;
        x(:,i) = out.x;
        time(i) = out.time;
    end
    % set(gcf, 'PaperSize', [25 25])
    % plot(variance,mean_val,'LineWidth',1.5);hold on
    % xlabel('variance of return','FontSize',15)
    % ylabel('mean of return','FontSize',15)
    % title('DAX 100','FontSize',15)
    
    % rowNames = {'1','2','3','4','5','6','7','8','9','10','mean'};
    result = [lambda,variance,mean_val];
    colNames = {'lambda','variance of return','mean of return'};
    result = array2table(result,'VariableNames',colNames)
    name = [save_path,'/port',num2str(k),'mean-variance.txt'];
    writetable(result,name)
    
    
    
    colNames = cell(1, 2000);  
    for i = 1:2000
        colNames{i} = num2str(i); 
    end
    
    rowNames = cell(1, number+1);  
    for i = 1:number
        rowNames{i} = num2str(i); 
    end
    rowNames{number+1} = 'time';
    result2 = [x;time];
    
    
    
    % result2 = [x;time];
    result2 = array2table(result2,'RowNames',rowNames,'VariableNames',colNames)
    name = [save_path, '/port',num2str(k),'x-time.txt'];
    writetable(result2,name)
end

%% efficient frontier
% sparsity
n = 10;

data = [[1,31];[2,85];[3,89];[4,98];[5,225]];
len=50;
for j = 1:size(data,1)
    k = data(j,1)
    num = data(j,2);
    filename = [path, 'port',num2str(k),'.txt'];
    [number,mean,deviance,correlation] = importfile(filename,[1,(1+num)*num/2+1+num]); %
    cov = diag(deviance)*correlation*diag(deviance);
    A = cov;
    lambda1 = max(max(abs(A)))*(log(n/(n-1)));

    time = zeros(1,len);
    x = zeros(number,len);
    variance = zeros(len,1);
    mean_val = zeros(len,1);
    lambda = zeros(len,1);
    for i = 1:len
        mu = (i)/len
        lambda(i) = mu;
        out = sparse_portfolio(mu*A,lambda1,(1-mu)*mean,1,n);
        variance(i) = out.x'*cov*out.x;
        mean_val(i) = out.x'*mean;
        x(:,i) = out.x;
        time(i) = out.time;
    end
    
    % plot(variance,mean_val,'o','MarkerFaceColor','r')
    
    
    result = [lambda,variance,mean_val];
    colNames = {'lambda','variance of return','mean of return'};
    result = array2table(result,'VariableNames',colNames);
    name = [save_path,'/port',num2str(k),'mean-variance_efficient.txt'];
    writetable(result,name) % 
    
    colNames = cell(1, 50);  
    for i = 1:50
        colNames{i} = num2str(i); 
    end
    
    rowNames = cell(1, number+1);  
    for i = 1:number
        rowNames{i} = num2str(i); 
    end
    rowNames{number+1} = 'time';
    result2 = [x;time];
    
    result2 = array2table(result2,'RowNames',rowNames,'VariableNames',colNames);
    name = [save_path, '/port',num2str(k),'x-time_efficient.txt'];
    writetable(result2,name)
    
    
    % legend('no sparsity','sparsity 10','FontSize',15)
end
%% calculate metrics
clear all;clc;
save_path = 'result/Experiment6';
for k = 1:5
    standard = readtable([save_path,'/port',num2str(k),'mean-variance.txt'],VariableNamingRule = 'preserve');
    efficient = readtable([save_path,'/port',num2str(k),'mean-variance_efficient.txt'],VariableNamingRule = 'preserve');
    standard = table2array(standard);
    efficient = table2array(efficient);

    standard = standard(:,2:3);
    Euclidean_distance = 0;
    variance_return = 0;
    mean_return = 0;
    
    for i = 1:50
        variance_mean_efficient = efficient(i,2:3);
        diff = variance_mean_efficient - standard;
        distance = vecnorm(diff,2,2);
        index = find(distance == min(distance));
        index = index(1);
        Euclidean_distance = Euclidean_distance + min(distance);
        variance_return = variance_return + abs(efficient(i,2) - standard(index,1))/efficient(i,2);
        mean_return = mean_return + abs(efficient(i,3) - standard(index,2))/efficient(i,3);
    end
    
    Euclidean_distance = Euclidean_distance/50;
    variance_return = variance_return*100/50;
    mean_return = mean_return*100/50;
    
    time1 = table2array(readtable([save_path, '/port',num2str(k),'x-time.txt'],VariableNamingRule = 'preserve')); %
    time2 = table2array(readtable([save_path, '/port',num2str(k),'x-time_efficient.txt'],VariableNamingRule = 'preserve')); %
    time = sum(time1(end,:))+sum(time2(end,:));
    
    result3 = [Euclidean_distance,variance_return,mean_return,time];
    colNames = {'Euclidean distance','variance of return error','mean of return error','time'};
    result3 = array2table(result3,'VariableNames',colNames)
    writetable(result3,[save_path, '/port',num2str(k),'_performance.txt']) % 
end
%% plot
titlenames = ["Hang Seng", "DAX 100", "FTSE 100", "S&P 100", "Nikkei"];
for k = 1:5
    portfolio_plot(save_path,k,titlenames(k))
end