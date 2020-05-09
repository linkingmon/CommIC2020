% Set Nt to K (MIMO channel)
% Apply Rayleigh channel

% initialize channel and algorithm parameters
Nt = 4;
K = 4;
EbN0 = 0:2:20; % in dB
sim_algorithm = ["ZF", "MMSE", "Kbest", "SD", "sorted_Kbest", "sorted_SD"];
min_error = 30;
max_sim_time = min_error/(K*log2(M))/min_ber;
total_trans_num = round(linspace(500, max_sim_time, length(EbN0)));

% initialize statistic paramters
algorithm_cnt = numel(sim_algorithm);
Nsnr = numel(EbN0);
BER = zeros(algorithm_cnt, Nsnr,1);
rng(0);

% initialize optimizer & simulator
my_simulator = Simulator(Nt, K, M);
my_detector = Detector(Nt, K, M);

fprintf("Running under system {Nt = %d ; K = %d ; M = %d} ...\n",Nt, K, M);
fprintf("EbN0 from %3d to %3d dB\n", EbN0(1), EbN0(end))
for m = 1:Nsnr
    fprintf("Running EbN0 %d dB ...\n", EbN0(m));
    snr = db2pow(EbN0(m)) * log2(M);
    err = zeros(algorithm_cnt, 1);
    tic;
    for n = 1:total_trans_num(m)
        % generate channels
        H = my_simulator.generate_channel();
        
        % generate transmit symbols
        [b_vec, transmit_symbol] = my_simulator.generate_TX_bit_and_symbol();

        % generate noise
        noise = my_simulator.generate_noise(transmit_symbol, snr);

        % transmit through channels
        receive_symbol = H * transmit_symbol + noise;

        % running different detection algorithm
        for i_algorithm = 1 : algorithm_cnt

            % MIMO detection
            detect_symbol = eval(strcat("my_detector.", sim_algorithm(i_algorithm),"_detect(receive_symbol, H, snr)"));

            % decode Rx symbols
            b_vec_r = my_simulator.decode_RX_symbol(detect_symbol);

            % record number of errors
            err(i_algorithm) = err(i_algorithm) + my_simulator.error_bits(b_vec, b_vec_r);
        end
    end
    fprintf("Consuming %3d min , %2d sec in EbN0 %2d dB\n",floor(ceil(toc)/60),mod(ceil(toc),60),EbN0(m))
    for i_algorithm = 1 : algorithm_cnt
        BER(i_algorithm, m,1) = err(i_algorithm)/total_trans_num(m)/K/log2(M);
    end
end

% save BER result
for i_algorithm = 1 : algorithm_cnt
    filename = [int2str(Nt), 'x', int2str(K),'_', int2str(M), '-QAM_', char(sim_algorithm(i_algorithm))];
    fprintf("Save result to \' %s\' ... \n\n", ['data/', filename]);
    cur_ber = BER(i_algorithm, :);
    save(['data/', filename], 'cur_ber');
end