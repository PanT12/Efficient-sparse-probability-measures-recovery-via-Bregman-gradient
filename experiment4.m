%% demo of the experiment 4 (Huber loss)

% random seed
rng("default")

% save result
mkdir('result/Experiment4')

% fix SNR=20
snr = 20;

% Salt and Pepper Impulse noise range
impulse_noise_densities = 0.1:0.01:0.2;

num = 100;
RSNR = zeros(length(impulse_noise_densities),num,4);

for i = 1:length(impulse_noise_densities)
    impulse_noise_density = impulse_noise_densities(i);
    fprintf('impulse noise density\t%e\n',impulse_noise_density)
    A = randn(200, 400);
    for k = 1:num
        fprintf('Iter\t%d\n',k)
        % ground truth x
        x_true = sprandn(400,1,0.02);
        x_true = abs(x_true)/norm(x_true,1);
    
        % gaussian noise
        n = randn(200,1);
        len = exp(log(norm(A*x_true)) - snr*log(10)/20); 
        m = n/norm(n)*len;
        
        % impulse noise
        % generate impluse noise
        noise_impulse = imnoise(zeros(size(m)), 'salt & pepper', impulse_noise_density); 
        % the maximum value of the impluse noise
        noise_impulse(noise_impulse == 1) = max(abs(m))*20; 
        

        b = A*x_true + m + noise_impulse;

        out = sparse_l0(A,10.0,b,1,1e-6,false); %1.8-4
        RSNR(i,k,1) = 10*log10(norm(x_true)^2/norm(out.x-x_true)^2);        

        out = huber_sparse_l0(A,10.0,b,1,1.0,false); %1.8-4
        RSNR(i,k,2) = 10*log10(norm(x_true)^2/norm(out.x-x_true)^2);

        out = huber_sparse_l0(A,5.0,b,1,1.0,false); %1.8-4
        RSNR(i,k,3) = 10*log10(norm(x_true)^2/norm(out.x-x_true)^2);

        out = huber_sparse_l0(A,1.0,b,1,1.0,false); %1.8-4
        RSNR(i,k,4) = 10*log10(norm(x_true)^2/norm(out.x-x_true)^2);
    end
end
path = 'result/Experiment4';
writematrix(RSNR,'result/Experiment4/RSNR.txt') % 


%%
clear all;clc
num=100;
RSNR = readmatrix("result/Experiment4/RSNR.txt");
RSNR = reshape(RSNR,11,num,4);
[m,~,~] = size(RSNR);
mean_RSNR = zeros(4,m);
std_RSNR = zeros(4,m);
for i = 1:m
    result = reshape(RSNR(i,:,:),num,4);
    mean_RSNR(:,i) = mean(result)';
    std_RSNR(:,i) = std(result)';
end
%%
set(gcf, 'PaperSize', [25 25])
p1 = plot(0.1:0.01:0.2,mean_RSNR(1,:),'-s',"MarkerSize",10, ...
    "MarkerEdgeColor","black","MarkerFaceColor",'#219ebc','LineWidth',1, 'Color','black');

hold on;
p2 = plot(0.1:0.01:0.2,mean_RSNR(2,:),'-s',"MarkerSize",10, ...
    "MarkerEdgeColor","black","MarkerFaceColor",'#ffb703','LineWidth',1, 'Color','black');
hold on;
p3 = plot(0.1:0.01:0.2,mean_RSNR(3,:),'-s',"MarkerSize",10, ...
    "MarkerEdgeColor","black","MarkerFaceColor",'red','LineWidth',1, 'Color','black');
hold on;
p4 = plot(0.1:0.01:0.2,mean_RSNR(4,:),'-s',"MarkerSize",10, ...
    "MarkerEdgeColor","black","MarkerFaceColor",'#9c3587','LineWidth',1, 'Color','black');

set(gca,'FontSize',15);
% xlim([13,52]);
ylim([0,20])
grid on;
xlabel('impulse noise density','Fontname','Times New Roman')
ylabel("RSNR (dB)",'Fontname','Times New Roman')


legend({'Quadratic loss (\lambda=10)','Huber loss (\lambda=10)','Huber loss (\lambda=5)','Huber loss (\lambda=1)'},'Location','northeast','NumColumns', 2, 'Fontname','Times New Roman')
hold off;