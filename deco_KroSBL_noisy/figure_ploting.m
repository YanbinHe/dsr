%% some pre-defined values
% this file aims to reproduce the figures in the paper
run('Configuration.m');


% [all_themes, all_colors] = GetColors();
% all_colors = all_themes{1,1};
all_colors = [0, 70/255, 222/255;
              0.6350, 0.0780, 0.1840;
              0, 0.4470, 0.7410
              0.9290, 0.6940, 0.1250
              255/255,0,1;
              0.49 0.18 0.56;];

all_colors_noise = [0, 70/255, 222/255;
              0.6350, 0.0780, 0.1840;
              0, 0.4470, 0.7410
              0.9290, 0.6940, 0.1250
              255/255,0,1;];

line_type_set{1} = '->';
line_type_set{2} = '-<';
line_type_set{3} = '-o';
line_type_set{4} = '-x';
line_type_set{5} = '-^';
line_type_set{6} = '-+';
line_type_set{7} = '-.>';
line_type_set{8} = '-.<';
line_type_set{9} = '-.o';
line_type_set{10} = '-.x';
line_type_set{11} = '-.^';
line_type_set{12} = '-.+';

line_type_noise{1} = '->';
line_type_noise{2} = '-<';
line_type_noise{3} = '-o';
line_type_noise{4} = '-x';
line_type_noise{5} = '-+';

legend_type_set{1} = '>';
legend_type_set{2} = '<';
legend_type_set{3} = 'o';
legend_type_set{4} = 'x';
legend_type_set{5} = '^';
legend_type_set{6} = '+';

legend_type_set_noise{1} = '>';
legend_type_set_noise{2} = '<';
legend_type_set_noise{3} = 'o';
legend_type_set_noise{4} = 'x';
legend_type_set_noise{5} = '+';

algo_name{1} = 'cSBL';
algo_name{2} = 'AM-KroSBL';
algo_name{3} = 'SVD-KroSBL';
algo_name{4} = 'dSBL';
algo_name{5} = 'OMP';
algo_name{6} = 'dOMP';

algo_name_noise{1} = 'cSBL';
algo_name_noise{2} = 'AM-KroSBL';
algo_name_noise{3} = 'SVD-KroSBL';
algo_name_noise{4} = 'dSBL';
algo_name_noise{5} = 'Truth';

num_alg = 6;
ratio = M1*M2*M3/N^3;

fontsizeman = 20;
opengl hardware
%%
% data process and figures plot
% take results.mat as input
% preprocessing
trials = 100;
% first index: algorithm; second index: 3 values (nmse/srr/time); third
% index: change with different conditions
resultSaggre = zeros(7,3,length(SNR));
resultMaggre = zeros(7,3,length(M1));
resultKaggre = zeros(7,3,length(K));

% fifth row is the true one
resultNoise = zeros(5,length(SNR));


for t = 1:trials
    filename = ['./results/noisy_compare2_', num2str(t),'.mat'];
    load(filename)
    for algo = 1:num_alg % for each algorithm
        metric = 1;
        for s = 1:length(SNR)
            resultSaggre(algo,metric,s) = resultSaggre(algo,metric,s) + resultS{s}{algo,1}{metric,2};
        end

        for m = 1:length(M1)
            resultMaggre(algo,metric,m) = resultMaggre(algo,metric,m) + resultM{m}{algo,1}{metric,2};
        end

        for k = 1:length(K)
            resultKaggre(algo,metric,k) = resultKaggre(algo,metric,k) + resultK{k}{algo,1}{metric,2};
        end

        metric = 3;
        for s = 1:length(SNR)
            resultSaggre(algo,metric-1,s) = resultSaggre(algo,metric-1,s) + resultS{s}{algo,1}{metric,2};
        end

        for m = 1:length(M1)
            resultMaggre(algo,metric-1,m) = resultMaggre(algo,metric-1,m) + resultM{m}{algo,1}{metric,2};
        end

        for k = 1:length(K)
            resultKaggre(algo,metric-1,k) = resultKaggre(algo,metric-1,k) + resultK{k}{algo,1}{metric,2};
        end

        metric = 4;
        for s = 1:length(SNR)
            resultSaggre(algo,metric-1,s) = resultSaggre(algo,metric-1,s) + resultS{s}{algo,1}{metric,2};
        end

        for m = 1:length(M1)
            resultMaggre(algo,metric-1,m) = resultMaggre(algo,metric-1,m) + resultM{m}{algo,1}{metric,2};
        end

        for k = 1:length(K)
            resultKaggre(algo,metric-1,k) = resultKaggre(algo,metric-1,k) + resultK{k}{algo,1}{metric,2};
        end
    end

    for algo = 1:3
        for s = 1:length(SNR)
            resultNoise(algo,s) = resultNoise(algo,s) + resultS{s}{algo,1}{5,2};
        end
    end
    algo = 4;
    for s = 1:length(SNR)
        resultNoise(algo,s) = resultNoise(algo,s) + resultS{s}{algo,1}{5,2}{1}+resultS{s}{algo,1}{5,2}{2}+resultS{s}{algo,1}{5,2}{3};
    end
    algo = 5;
    for s = 1:length(SNR)
        resultNoise(algo,s) = resultNoise(algo,s) + resultS{s}{algo+2,1};
    end
end

resultSaggre = resultSaggre/trials;
resultMaggre = resultMaggre/trials;
resultKaggre = resultKaggre/trials;
resultNoise = resultNoise/trials;
%% fig 1/2 (a) NMSE/SRR vs SNR
figure
metric = 1;
for algo_index = 1:num_alg
    plot(SNR,reshape(resultSaggre(algo_index,metric,:),[6,1]),line_type_set{algo_index},'Color',all_colors(algo_index, :),'Display',algo_name{algo_index});
    hold on
end

grid on
set(gca, 'yscale', 'log');
legend('boxoff')
set(0,'DefaultLineLineWidth',3)
set(0,'DefaultAxesFontSize',fontsizeman)
set(0,'DefaultLineMarkerSize',14)
set(0,'DefaultAxesFontWeight','bold')
set(gca,'FontSize',fontsizeman)
set(get(gca,'Xlabel'),'FontSize',fontsizeman)
set(get(gca,'Ylabel'),'FontSize',fontsizeman)
set(get(gca,'Title'),'FontSize',fontsizeman)
set(get(gca,'Xlabel'),'FontWeight','bold')
set(get(gca,'Ylabel'),'FontWeight','bold')
set(get(gca,'Title'),'FontWeight','bold')
xlabel('SNR (dB)')
ylabel('NMSE')
% axis square

h1 = plot(SNR,reshape(resultSaggre(1,metric,:),[6,1]),legend_type_set{1},'Color',all_colors(1, :),'Display',algo_name{1});
hold on
h2 = plot(SNR,reshape(resultSaggre(2,metric,:),[6,1]),legend_type_set{2},'Color',all_colors(2, :),'Display',algo_name{2});
hold on
h3 = plot(SNR,reshape(resultSaggre(3,metric,:),[6,1]),legend_type_set{3},'Color',all_colors(3, :),'Display',algo_name{3});
hold on
h4 = plot(SNR,reshape(resultSaggre(4,metric,:),[6,1]),legend_type_set{4},'Color',all_colors(4, :),'Display',algo_name{4});
hold on
h5 = plot(SNR,reshape(resultSaggre(5,metric,:),[6,1]),legend_type_set{5},'Color',all_colors(5, :),'Display',algo_name{5});
hold on
h6 = plot(SNR,reshape(resultSaggre(6,metric,:),[6,1]),legend_type_set{6},'Color',all_colors(6, :),'Display',algo_name{6});
hold on
legend([h1 h2 h3 h4 h5 h6],{algo_name{1},algo_name{2},algo_name{3},algo_name{4},algo_name{5},algo_name{6}},'Location','northeast','Interpreter','LaTex')


% axes('Position',[.55 .18 .3 .3])
% 
% for algo_index = 1:5
%     plot(SNR,resultNoise(algo_index,:),line_type_noise{algo_index},'Color',all_colors_noise(algo_index, :),'Display',algo_name_noise{algo_index});
%     hold on
% end
% set(gca, 'yscale', 'log');
% fontsizemanlocal = 20;
% set(0,'DefaultLineLineWidth',3)
% set(0,'DefaultAxesFontSize',fontsizemanlocal)
% set(0,'DefaultLineMarkerSize',14)
% set(0,'DefaultAxesFontWeight','bold')
% set(gca,'FontSize',fontsizemanlocal)
% set(get(gca,'Xlabel'),'FontSize',fontsizemanlocal)
% set(get(gca,'Ylabel'),'FontSize',fontsizemanlocal)
% set(get(gca,'Title'),'FontSize',fontsizemanlocal)
% set(get(gca,'Xlabel'),'FontWeight','bold')
% set(get(gca,'Ylabel'),'FontWeight','bold')
% set(get(gca,'Title'),'FontWeight','bold')
% box on
% grid on
% legend off



figure
metric = 2;
for algo_index = 1:num_alg
    plot(SNR,reshape(resultSaggre(algo_index,metric,:),[6,1]),line_type_set{algo_index},'Color',all_colors(algo_index, :),'Display',algo_name{algo_index});
    hold on
end

grid on
% set(gca, 'yscale', 'log');
legend('boxoff')
set(0,'DefaultLineLineWidth',3)
set(0,'DefaultAxesFontSize',fontsizeman)
set(0,'DefaultLineMarkerSize',14)
set(0,'DefaultAxesFontWeight','bold')
set(gca,'FontSize',fontsizeman)
set(get(gca,'Xlabel'),'FontSize',fontsizeman)
set(get(gca,'Ylabel'),'FontSize',fontsizeman)
set(get(gca,'Title'),'FontSize',fontsizeman)
set(get(gca,'Xlabel'),'FontWeight','bold')
set(get(gca,'Ylabel'),'FontWeight','bold')
set(get(gca,'Title'),'FontWeight','bold')
xlabel('SNR (dB)')
ylabel('SRR')
% axis square

h1 = plot(SNR,reshape(resultSaggre(1,metric,:),[6,1]),legend_type_set{1},'Color',all_colors(1, :),'Display',algo_name{1});
hold on
h2 = plot(SNR,reshape(resultSaggre(2,metric,:),[6,1]),legend_type_set{2},'Color',all_colors(2, :),'Display',algo_name{2});
hold on
h3 = plot(SNR,reshape(resultSaggre(3,metric,:),[6,1]),legend_type_set{3},'Color',all_colors(3, :),'Display',algo_name{3});
hold on
h4 = plot(SNR,reshape(resultSaggre(4,metric,:),[6,1]),legend_type_set{4},'Color',all_colors(4, :),'Display',algo_name{4});
hold on
h5 = plot(SNR,reshape(resultSaggre(5,metric,:),[6,1]),legend_type_set{5},'Color',all_colors(5, :),'Display',algo_name{5});
hold on
h6 = plot(SNR,reshape(resultSaggre(6,metric,:),[6,1]),legend_type_set{6},'Color',all_colors(6, :),'Display',algo_name{6});
hold on
legend([h1 h2 h3 h4 h5 h6],{algo_name{1},algo_name{2},algo_name{3},algo_name{4},algo_name{5},algo_name{6}},'Location','northeast','Interpreter','LaTex')
%% figure 1/2 (b) NMSE/SRR vs under-sampling
figure
metric = 1;
for algo_index = 1:num_alg
    plot(ratio,reshape(resultMaggre(algo_index,metric,:),[7,1]),line_type_set{algo_index},'Color',all_colors(algo_index, :),'Display',algo_name{algo_index});
    hold on
end

grid on
set(gca, 'yscale', 'log');
legend('boxoff')
set(0,'DefaultLineLineWidth',3)
set(0,'DefaultAxesFontSize',fontsizeman)
set(0,'DefaultLineMarkerSize',14)
set(0,'DefaultAxesFontWeight','bold')
set(gca,'FontSize',fontsizeman)
set(get(gca,'Xlabel'),'FontSize',fontsizeman)
set(get(gca,'Ylabel'),'FontSize',fontsizeman)
set(get(gca,'Title'),'FontSize',fontsizeman)
set(get(gca,'Xlabel'),'FontWeight','bold')
set(get(gca,'Ylabel'),'FontWeight','bold')
set(get(gca,'Title'),'FontWeight','bold')
xlabel('Under-sampling ratio')
ylabel('NMSE')
% axis square

h1 = plot(ratio,reshape(resultMaggre(1,metric,:),[7,1]),legend_type_set{1},'Color',all_colors(1, :),'Display',algo_name{1});
hold on
h2 = plot(ratio,reshape(resultMaggre(2,metric,:),[7,1]),legend_type_set{2},'Color',all_colors(2, :),'Display',algo_name{2});
hold on
h3 = plot(ratio,reshape(resultMaggre(3,metric,:),[7,1]),legend_type_set{3},'Color',all_colors(3, :),'Display',algo_name{3});
hold on
h4 = plot(ratio,reshape(resultMaggre(4,metric,:),[7,1]),legend_type_set{4},'Color',all_colors(4, :),'Display',algo_name{4});
hold on
h5 = plot(ratio,reshape(resultMaggre(5,metric,:),[7,1]),legend_type_set{5},'Color',all_colors(5, :),'Display',algo_name{5});
hold on
h6 = plot(ratio,reshape(resultMaggre(6,metric,:),[7,1]),legend_type_set{6},'Color',all_colors(6, :),'Display',algo_name{6});
hold on
legend([h1 h2 h3 h4 h5 h6],{algo_name{1},algo_name{2},algo_name{3},algo_name{4},algo_name{5},algo_name{6}},'Location','southwest','Interpreter','LaTex')
% saveas(gcf,'deco2','epsc')
figure
metric = 2;
% axis square
h1 = plot(ratio,reshape(resultMaggre(1,metric,:),[7,1]),legend_type_set{1},'Color',all_colors(1, :),'Display',algo_name{1});
hold on
h2 = plot(ratio,reshape(resultMaggre(2,metric,:),[7,1]),legend_type_set{2},'Color',all_colors(2, :),'Display',algo_name{2});
hold on
h3 = plot(ratio,reshape(resultMaggre(3,metric,:),[7,1]),legend_type_set{3},'Color',all_colors(3, :),'Display',algo_name{3});
hold on
h4 = plot(ratio,reshape(resultMaggre(4,metric,:),[7,1]),legend_type_set{4},'Color',all_colors(4, :),'Display',algo_name{4});
hold on
h5 = plot(ratio,reshape(resultMaggre(5,metric,:),[7,1]),legend_type_set{5},'Color',all_colors(5, :),'Display',algo_name{5});
hold on
h6 = plot(ratio,reshape(resultMaggre(6,metric,:),[7,1]),legend_type_set{6},'Color',all_colors(6, :),'Display',algo_name{6});
hold on
for algo_index = 1:num_alg
    plot(ratio,reshape(resultMaggre(algo_index,metric,:),[7,1]),line_type_set{algo_index},'Color',all_colors(algo_index, :),'Display',algo_name{algo_index});
    hold on
end
legend([h1 h2 h3 h4 h5 h6],{algo_name{1},algo_name{2},algo_name{3},algo_name{4},algo_name{5},algo_name{6}},'Location','northwest','Interpreter','LaTex')
grid on
legend('boxoff')
set(0,'DefaultLineLineWidth',3)
set(0,'DefaultAxesFontSize',fontsizeman)
set(0,'DefaultLineMarkerSize',14)
set(0,'DefaultAxesFontWeight','bold')
set(gca,'FontSize',fontsizeman)
set(get(gca,'Xlabel'),'FontSize',fontsizeman)
set(get(gca,'Ylabel'),'FontSize',fontsizeman)
set(get(gca,'Title'),'FontSize',fontsizeman)
set(get(gca,'Xlabel'),'FontWeight','bold')
set(get(gca,'Ylabel'),'FontWeight','bold')
set(get(gca,'Title'),'FontWeight','bold')
xlabel('Under-sampling ratio')
ylabel('SRR')

% saveas(gcf,'deco3','epsc')
%% figure 1/2 (c) NMSE/SRR vs sparsity level
figure
metric = 1;

for algo_index = 1:num_alg
    plot(K,reshape(resultKaggre(algo_index,metric,:),[5,1]),line_type_set{algo_index},'Color',all_colors(algo_index, :),'Display',algo_name{algo_index});
    hold on
end
set(gca, 'yscale', 'log');
legend('boxoff')
set(0,'DefaultLineLineWidth',3)
set(0,'DefaultAxesFontSize',fontsizeman)
set(0,'DefaultLineMarkerSize',14)
set(0,'DefaultAxesFontWeight','bold')
set(gca,'FontSize',fontsizeman)
set(get(gca,'Xlabel'),'FontSize',fontsizeman)
set(get(gca,'Ylabel'),'FontSize',fontsizeman)
set(get(gca,'Title'),'FontSize',fontsizeman)
set(get(gca,'Xlabel'),'FontWeight','bold')
set(get(gca,'Ylabel'),'FontWeight','bold')
set(get(gca,'Title'),'FontWeight','bold')
grid on
% ylim([1e-3,5]);
% axis square
h1 = plot(K,reshape(resultKaggre(1,metric,:),[5,1]),legend_type_set{1},'Color',all_colors(1, :),'Display',algo_name{1});
hold on
h2 = plot(K,reshape(resultKaggre(2,metric,:),[5,1]),legend_type_set{2},'Color',all_colors(2, :),'Display',algo_name{2});
hold on
h3 = plot(K,reshape(resultKaggre(3,metric,:),[5,1]),legend_type_set{3},'Color',all_colors(3, :),'Display',algo_name{3});
hold on
h4 = plot(K,reshape(resultKaggre(4,metric,:),[5,1]),legend_type_set{4},'Color',all_colors(4, :),'Display',algo_name{4});
hold on
h5 = plot(K,reshape(resultKaggre(5,metric,:),[5,1]),legend_type_set{5},'Color',all_colors(5, :),'Display',algo_name{5});
hold on
legend([h1 h2 h3 h4 h5],{algo_name{1},algo_name{2},algo_name{3},algo_name{4},algo_name{5}},'Location','southeast','Interpreter','LaTex')

xlabel('Sparsity level')
ylabel('NMSE')
% saveas(gcf,'deco4','epsc')
% % create a new pair of axes inside current figure
% axes('position',[.55 .2 .35 .4])
% box on % put box around new pair of axes
% indexOfInterest = 3:5; % range of t near perturbation
% for algo_index = [1,2,3,5]
%     plot(K(indexOfInterest),reshape(resultKaggre(algo_index,metric,indexOfInterest),[3,1]),line_type_set{algo_index},'Color',all_colors(algo_index, :));
%     hold on
% end
% axis tight
% grid on
% set(gca, 'yscale', 'log');
% set(0,'DefaultLineLineWidth',3)
% set(0,'DefaultAxesFontSize',fontsizeman)
% set(0,'DefaultLineMarkerSize',14)
% set(0,'DefaultAxesFontWeight','bold')
% set(gca,'FontSize',fontsizeman)
% set(gca,'YTick',[])
% set(get(gca,'Xlabel'),'FontSize',fontsizeman)
% set(get(gca,'Ylabel'),'FontSize',fontsizeman)
% set(get(gca,'Title'),'FontSize',fontsizeman)
% set(get(gca,'Xlabel'),'FontWeight','bold')
% set(get(gca,'Ylabel'),'FontWeight','bold')
% set(get(gca,'Title'),'FontWeight','bold')
%
figure
metric = 2;
for algo_index = 1:num_alg
    plot(K,reshape(resultKaggre(algo_index,metric,:),[5,1]),line_type_set{algo_index},'Color',all_colors(algo_index, :),'Display',algo_name{algo_index});
    hold on
end

grid on
% set(gca, 'yscale', 'log');
legend('boxoff')
set(0,'DefaultLineLineWidth',3)
set(0,'DefaultAxesFontSize',fontsizeman)
set(0,'DefaultLineMarkerSize',14)
set(0,'DefaultAxesFontWeight','bold')
set(gca,'FontSize',fontsizeman)
set(get(gca,'Xlabel'),'FontSize',fontsizeman)
set(get(gca,'Ylabel'),'FontSize',fontsizeman)
set(get(gca,'Title'),'FontSize',fontsizeman)
set(get(gca,'Xlabel'),'FontWeight','bold')
set(get(gca,'Ylabel'),'FontWeight','bold')
set(get(gca,'Title'),'FontWeight','bold')
xlabel('Sparsity level')
ylabel('SRR')
% axis square

h1 = plot(K,reshape(resultKaggre(1,metric,:),[5,1]),legend_type_set{1},'Color',all_colors(1, :),'Display',algo_name{1});
hold on
h2 = plot(K,reshape(resultKaggre(2,metric,:),[5,1]),legend_type_set{2},'Color',all_colors(2, :),'Display',algo_name{2});
hold on
h3 = plot(K,reshape(resultKaggre(3,metric,:),[5,1]),legend_type_set{3},'Color',all_colors(3, :),'Display',algo_name{3});
hold on
h4 = plot(K,reshape(resultKaggre(4,metric,:),[5,1]),legend_type_set{4},'Color',all_colors(4, :),'Display',algo_name{4});
hold on
h5 = plot(K,reshape(resultKaggre(5,metric,:),[5,1]),legend_type_set{5},'Color',all_colors(5, :),'Display',algo_name{5});
hold on
legend([h1 h2 h3 h4 h5],{algo_name{1},algo_name{2},algo_name{3},algo_name{4},algo_name{5}},'Location','southwest','Interpreter','LaTex')
% saveas(gcf,'deco5','epsc')
%% plot computation time

timeAggre = zeros(6,7);

for t = 1:trials % for different independent runs
    filename = ['./results/noisy_compare2_', num2str(t),'.mat'];
    load(filename)
    for algo = 1:num_alg % for each algorithm
        for m = 1:length(M1)
            timeAggre(algo,m) = timeAggre(algo,m) + resultM{m}{algo,1}{4,2};
        end
    end
end

timeAggre = timeAggre/trials;

figure
for algo_index = 1:num_alg
    plot(ratio,timeAggre(algo_index,:),line_type_set{algo_index},'Color',all_colors(algo_index, :),'Display',algo_name{algo_index});
    hold on
end
legend('boxoff')
% title('(a)')
grid on
set(gca, 'yscale', 'log');
set(0,'DefaultLineLineWidth',3)
set(0,'DefaultAxesFontSize',fontsizeman)
set(0,'DefaultLineMarkerSize',14)
set(0,'DefaultAxesFontWeight','bold')
set(gca,'FontSize',fontsizeman)
set(get(gca,'Xlabel'),'FontSize',fontsizeman)
set(get(gca,'Ylabel'),'FontSize',fontsizeman)
set(get(gca,'Title'),'FontSize',fontsizeman)
set(get(gca,'Xlabel'),'FontWeight','bold')
set(get(gca,'Ylabel'),'FontWeight','bold')
set(get(gca,'Title'),'FontWeight','bold')
xlabel('Under-sampling ratio')
ylabel('Seconds (s)')

%%
h1 = plot(ratio,timeAggre(1,:),legend_type_set{1},'Color',all_colors(1, :),'Display',algo_name{1});
hold on
h2 = plot(ratio,timeAggre(2,:),legend_type_set{2},'Color',all_colors(2, :),'Display',algo_name{2});
hold on
h3 = plot(ratio,timeAggre(3,:),legend_type_set{3},'Color',all_colors(3, :),'Display',algo_name{3});
hold on
h4 = plot(ratio,timeAggre(4,:),legend_type_set{4},'Color',all_colors(4, :),'Display',algo_name{4});
hold on
h5 = plot(ratio,timeAggre(5,:),legend_type_set{5},'Color',all_colors(5, :),'Display',algo_name{5});
hold on
legend([h1 h2 h3 h4 h5],{algo_name{1},algo_name{2},algo_name{3},algo_name{4},algo_name{5}},'Location','southwest','Interpreter','LaTex')
