% Simulator
classdef Simulator
    properties
        Nt;             % number of transmit antennas
        K;              % number of users, current supporting values (equals to Nt)
        M;              % QAM order, current supporting values (4, 16, 64, 256)
    end
    methods
        % constructor
        function obj = Simulator(Nt, K, M)
            % global initialize parameters
            obj.Nt = Nt;
            obj.K = K;
            obj.M = M;
        end
        % generate rayleigh channels
        function H = generate_channel(obj)
            H = ( randn(obj.K,obj.Nt) + 1i*randn(obj.K,obj.Nt) ) / sqrt(2);
        end
        % transmit a random generated bit stream
        function [b_vec, transmit_symbol] = generate_TX_bit_and_symbol(obj)
            b_vec = randi([0 1],obj.K,log2(obj.M));
            id_vec = bi2de(b_vec);
            transmit_symbol = qammod(id_vec, obj.M, 'gray');
        end
        % receive symbol
        function b_vec_r = decode_RX_symbol(obj, receive_symbol)
            id_vec_r = qamdemod(receive_symbol, obj.M, 'gray');
            b_vec_r = de2bi(id_vec_r, log2(obj.M));
        end
        % calculate number of error bits
        function err = error_bits(obj, b_vec, b_vec_r)
            err_ary = (b_vec ~= b_vec_r);
            err  = sum(err_ary(:));
        end
        % generate noise based in transmit symbol energy and SNR
        function noise = generate_noise(obj, transmit_symbol, snr)            
            Es = transmit_symbol' * transmit_symbol;
            sigma = Es/snr;
            noise = sqrt(sigma/2).*(randn(obj.K,1)+1i*randn(obj.K,1));
        end
    end
end