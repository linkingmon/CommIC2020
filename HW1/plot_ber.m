clc, clear
Nt = 4;
K = 4;
EbN0 = 0:2:20;
sim_algorithm = ["ZF", "MMSE", "Kbest", "sorted_Kbest", "SD", "sorted_SD"];
sim_algorithm_name = ["ZF", "MMSE", "Kbest", "sorted Kbest", "SD", "sorted SD"];
color = ['r', 'b', 'm', 'c', 'k', 'g'];
line_type = ["-s", "-x", "-o", "-^"];
M = [16 64];
algorithm_cnt = numel(sim_algorithm);

for iM = 1 : numel(M)
    ber = [];
    for i_algorithm=1:algorithm_cnt
        filename = [int2str(Nt), 'x', int2str(K),'_', int2str(M(iM)), '-QAM_', char(sim_algorithm(i_algorithm))];
        ber1 = load(['data/', filename, '.mat']);
        BER = [ber1.cur_ber.' ; zeros(length(EbN0)-length(ber1.cur_ber),1)];
        ber = [ber BER];
    end

    for i_algorithm=1:algorithm_cnt
        legend_name = strcat(sim_algorithm_name(i_algorithm), " ", int2str(M(iM)), "-QAM");
        leg{(i_algorithm-1)*numel(M)+iM} = sprintf(legend_name);
        plot(EbN0,ber(:,i_algorithm), line_type(iM), 'LineWidth', 2, 'Color', color(i_algorithm));
        hold on;
    end
end

set(gca, 'YScale', 'log')
ylim([1e-5 1]);
xlabel('EbN0 (dB)');
ylabel('Bit error rate');
title(['BER - EbN0 (Nt=',int2str(Nt),',K=',int2str(K),',',int2str(M),'QAM)' ]);
legend(leg);
grid on;