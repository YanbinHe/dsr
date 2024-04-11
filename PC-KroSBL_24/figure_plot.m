%%
run('Configuration.m');

trial = 50;

lengt = 6;
num_alg = 6;
errorAggr = zeros(num_alg,lengt);
srrAggr = zeros(num_alg,lengt);
timeAggr = zeros(num_alg,lengt);
serAggr = zeros(num_alg,lengt);
noiseAggr = zeros(num_alg,lengt);

for t = 1:trial
    filename = ['./Results/ce_', num2str(t),'.mat'];
    load(filename)
    for algo_index = 1:num_alg
        errorAggr(algo_index,:) = errorAggr(algo_index,:) + error(algo_index,:);
        timeAggr(algo_index,:) = timeAggr(algo_index,:) + time(algo_index,:);
        serAggr(algo_index,:) = serAggr(algo_index,:) + ser(algo_index,:);
    end
end
errorAggr = errorAggr/trial;
timeAggr = timeAggr/trial;
serAggr = serAggr/trial;

% for algo_index = 1:num_alg
%     errorAggr(algo_index,:) = sum(error{algo_index,2}(:,:),2)/trial;
%     srrAggr(algo_index,:) = sum(supprecovery{algo_index,2}(:,:),2)/trial;
%     timeAggr(algo_index,:) = sum(time{algo_index,2}(:,:),2)/trial;
%     serAggr(algo_index,:) = sum(ser{algo_index,2}(:,:),2)/trial;
% end
% for algo_index = [1,2,4]
% noiseAggr(algo_index,:) = sum(noise_var_error{algo_index,2}(:,:),2)/trial;
% end
%
fontsizeman = 20;
          
all_colors = [0      , 0.4470, 0.7410 ;
              0.6350 , 0.0780, 0.1840 ;
              0.9290 , 0.6940, 0.1250 ;
              0      , 70/255, 222/255;
              1      , 0     , 1      ;
              0.49   , 0.18  , 0.56   ;
              ];


line_type_set{1} = '-o';
line_type_set{2} = '-<';
line_type_set{3} = '-x';
line_type_set{4} = '->';
line_type_set{5} = '-^';
line_type_set{6} = '-+';
line_type_set{7} = '-.o';
line_type_set{8} = '-.<';
line_type_set{9} = '-.x';
line_type_set{10} = '-.>';
line_type_set{11} = '-.^';
line_type_set{12} = '-.+';

legend_type_set{1} = 'o';
legend_type_set{2} = '<';
legend_type_set{3} = 'x';
legend_type_set{4} = '>';
legend_type_set{5} = '^';
legend_type_set{6} = '+';


algo_name{1} = 'SVD-KroSBL';
algo_name{2} = 'AM-KroSBL';
algo_name{3} = 'dSBL';
algo_name{4} = 'cSBL';
algo_name{5} = 'OMP';
algo_name{6} = 'dOMP';
%
figure
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
% yyaxis left
for algo_index = 1:num_alg
    plot(SNRl,errorAggr(algo_index,:),line_type_set{algo_index},'Display',algo_name{algo_index},'Color',all_colors(algo_index, :));
    hold on
end

h1 = plot(SNRl,errorAggr(1,:),legend_type_set{1},'Display',algo_name{1},'Color',all_colors(1, :));
hold on
h2 = plot(SNRl,errorAggr(2,:),legend_type_set{2},'Display',algo_name{2},'Color',all_colors(2, :));
hold on
h3 = plot(SNRl,errorAggr(3,:),legend_type_set{3},'Display',algo_name{3},'Color',all_colors(3, :));
hold on
h4 = plot(SNRl,errorAggr(4,:),legend_type_set{4},'Display',algo_name{4},'Color',all_colors(4, :));
hold on
h5 = plot(SNRl,errorAggr(5,:),legend_type_set{5},'Display',algo_name{5},'Color',all_colors(5, :));
hold on
h6 = plot(SNRl,errorAggr(6,:),legend_type_set{6},'Display',algo_name{6},'Color',all_colors(6, :));
hold on
legend([h4 h5 h2 h1 h3 h6],{algo_name{4},algo_name{5},algo_name{2},algo_name{1},algo_name{3},algo_name{6}},'Location','northeast','Interpreter','LaTex')

grid on
set(gca, 'yscale', 'log');
legend('boxoff')

xticks(SNRl)
xlabel('SNR (dB)')
ylabel('NMSE')

% yyaxis right
% fontsizemanlocal = 10;
% ylim([1e-3 1e0])
% axes('Position',[.18 .18 .3 .3])
% box on
% for algo_index = [1,2,4]
%     plot(SNRl,noiseAggr(algo_index,:),line_type_set{algo_index+5},'Color',all_colors(algo_index, :));
%     hold on
% end
% legend off
% grid on
% set(gca, 'yscale', 'log');
% set(0,'DefaultLineLineWidth',3)
% set(0,'DefaultAxesFontSize',fontsizemanlocal)
% set(0,'DefaultLineMarkerSize',7)
% set(0,'DefaultAxesFontWeight','bold')
% set(gca,'FontSize',fontsizemanlocal)
% set(get(gca,'Xlabel'),'FontSize',fontsizemanlocal)
% set(get(gca,'Ylabel'),'FontSize',fontsizemanlocal)
% set(get(gca,'Title'),'FontSize',fontsizemanlocal)
% set(get(gca,'Xlabel'),'FontWeight','bold')
% set(get(gca,'Ylabel'),'FontWeight','bold')
% set(get(gca,'Title'),'FontWeight','bold')
% xticks(SNRl)
% xlabel('SNR (dB)','Interpreter','Latex')
% ylabel('NMSE','Interpreter','Latex')

%%
% figure
% grid on
% % set(gca, 'yscale', 'log');
% legend('boxoff')
% set(0,'DefaultLineLineWidth',3)
% set(0,'DefaultAxesFontSize',fontsizeman)
% set(0,'DefaultLineMarkerSize',14)
% set(0,'DefaultAxesFontWeight','bold')
% set(gca,'FontSize',fontsizeman)
% set(get(gca,'Xlabel'),'FontSize',fontsizeman)
% set(get(gca,'Ylabel'),'FontSize',fontsizeman)
% set(get(gca,'Title'),'FontSize',fontsizeman)
% set(get(gca,'Xlabel'),'FontWeight','bold')
% set(get(gca,'Ylabel'),'FontWeight','bold')
% set(get(gca,'Title'),'FontWeight','bold')
% xticks(SNRl)
% xlabel('SNR (dB)')
% ylabel('SRR')
% for algo_index = 1:num_alg
%     plot(SNRl,srrAggr(algo_index,:),line_type_set{algo_index},'Display',algo_name{algo_index},'Color',all_colors(algo_index, :));
%     hold on
% end
% 
% h1 = plot(SNRl,srrAggr(1,:),legend_type_set{1},'Display',algo_name{1},'Color',all_colors(1, :));
% hold on
% h2 = plot(SNRl,srrAggr(2,:),legend_type_set{2},'Display',algo_name{2},'Color',all_colors(2, :));
% hold on
% h3 = plot(SNRl,srrAggr(3,:),legend_type_set{3},'Display',algo_name{3},'Color',all_colors(3, :));
% hold on
% h4 = plot(SNRl,srrAggr(4,:),legend_type_set{4},'Display',algo_name{4},'Color',all_colors(4, :));
% hold on
% h5 = plot(SNRl,srrAggr(5,:),legend_type_set{5},'Display',algo_name{5},'Color',all_colors(5, :));
% hold on
% legend([h4 h5 h2 h1 h3],{algo_name{4},algo_name{5},algo_name{2},algo_name{1},algo_name{3}},'Location','northeast','Interpreter','LaTex')



%
figure
for algo_index = 1:num_alg
    plot(SNRl,serAggr(algo_index,:),line_type_set{algo_index},'Display',algo_name{algo_index},'Color',all_colors(algo_index, :));
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
xticks(SNRl)
xlabel('SNR (dB)')
ylabel('SER')

figure
for algo_index = 1:num_alg
    plot(SNRl,timeAggr(algo_index,:),line_type_set{algo_index},'Display',algo_name{algo_index},'Color',all_colors(algo_index, :));
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
xticks(SNRl)
xlabel('SNR (dB)','Interpreter','Latex')
ylabel('Second (s)','Interpreter','Latex')