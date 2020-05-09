% CI-optimizer
classdef Detector
    properties
        Nt;             % number of transmit antennas
        K;              % number of users, current supporting values (equals to Nt)
        M;              % QAM order, current supporting values (4, 16, 64, 256)
        symbol_set;     % the allowable symbol for M-QAM modulation
        s_out;          % for recording the best solution found in SD methods
        cur_min;        % for recording the best calue found in the SD methods
        max_K;          % the K-best parameter
    end
    methods
        % constructor
        function obj = Detector(Nt, K, M)
            % global initialize parameters
            obj.Nt = Nt;
            obj.K = K;
            obj.M = M;
            obj.symbol_set = [];
            for ii = 0 : M-1
                obj.symbol_set = [obj.symbol_set qammod(ii, M, 'gray')];
            end
        end

        % ZF decoding
        function s_out = ZF_detect(obj, s_in, H, snr)
            s_out = pinv(H)*s_in;
        end
        
        % MMSE decoding
        function s_out = MMSE_detect(obj, s_in, H, snr)
            s_out = inv(H'*H+eye(obj.K)/snr)*H'*s_in;
        end

        % Kbest decoding
        function s_out = Kbest_detect(obj, s_in, H, snr)
            function traverse(lev, s_in, cur_sols, cur_costs)
                next_sols = {};
                next_costs = [];
                % calculate max_K * M solutions
                for i_kbest = 1 : size(cur_costs, 2)
                    for node_idx = 1 : obj.M
                        next_sol = cur_sols{i_kbest};
                        next_cost = cur_costs(i_kbest);
                        next_sol(lev) = obj.symbol_set(node_idx);
                        cur_partial = abs(( R(lev,lev:obj.K)*next_sol(lev:obj.K, 1) ) - s_in(lev))^2;
                        next_cost = next_cost + cur_partial;
                        next_sols{ (i_kbest-1)*obj.M + node_idx } = next_sol;
                        next_costs( (i_kbest-1)*obj.M + node_idx ) = next_cost;
                    end
                end
                % traverse next level if not leaf
                if lev == 1
                    [min_val, min_idx] = min(next_costs);
                    obj.s_out = next_sols{min_idx};
                    return
                end
                % only save max_K solution
                if size(next_costs,2) > obj.max_K
                    [next_costs, mink_indexs] = mink(next_costs, obj.max_K);
                    next_sols = next_sols(mink_indexs);
                end
                traverse(lev-1, s_in, next_sols, next_costs);
            end
            obj.max_K = 16;
            [Q, R] = qr(H);                                 % Q' * s_in = R * s_out
            init_sol = {obj.symbol_set(1) .* ones(obj.K,1)};  % initialize a solution
            init_costs = 0;                                 % initial cost 0
            traverse(obj.K, Q'*s_in, init_sol, init_costs);
            s_out = obj.s_out;
        end
        
        % SD decoding
        function s_out = SD_detect(obj, s_in, H, snr)
            function traverse(lev, idx, cur_sum, cur_sol, s_in)
                % strcat("traverse ", int2str(lev), " ", int2str(idx), " ", num2str(cur_sum), " ", num2str(obj.cur_min), "\n")
                % cur_sol.'
                next_sol = cur_sol;
                next_sol(lev) = obj.symbol_set(idx);
                cur_partial = abs(( R(lev,lev:obj.K)*next_sol(lev:obj.K, 1) ) - s_in(lev))^2;
                if (lev == 1)
                    if (cur_sum + cur_partial < obj.cur_min)
                        obj.cur_min = cur_sum + cur_partial;
                        obj.s_out = next_sol;
                    end
                    return
                end
                if (cur_sum + cur_partial > obj.cur_min)
                    return
                end
                for ii = 1 : obj.M
                    traverse(lev-1, ii, cur_sum + cur_partial, next_sol, s_in);
                end
            end
            [Q, R] = qr(H);                                % Q' * s_in = R * s_out
            obj.cur_min = inf;                             % initialize min with inf
            cur_sol = obj.symbol_set(1) .* ones(obj.K,1);  % initialize a solution
            for node_idx = 1 : obj.M
                traverse(obj.K, node_idx, 0, cur_sol, Q'*s_in);
            end
            s_out = obj.s_out;
        end

        % optimal sorting and suboptimal sorting (refer to p.9)
        function idx = suboptimal_sorting(obj, H)
            col_norm = vecnorm(H,2,1);
            [val, idx] = sort(col_norm);
        end

        % sorted_Kbest decoding
        function s_out = sorted_Kbest_detect(obj, s_in, H, snr)
            idx = obj.suboptimal_sorting(H);             % find sorting index
            s_out = obj.Kbest_detect(s_in, H(:,idx), snr);  % fun algorithm in sorting H
            r_idx(idx) = 1:obj.K;                        % calculate reverse permute vector
            s_out = s_out(r_idx);                        % reverse permute s_out
        end
        
        % sorted_SD decoding
        function s_out = sorted_SD_detect(obj, s_in, H, snr)
            idx = obj.suboptimal_sorting(H);             % find sorting index
            s_out = obj.SD_detect(s_in, H(:,idx), snr);  % fun algorithm in sorting H
            r_idx(idx) = 1:obj.K;                        % calculate reverse permute vector
            s_out = s_out(r_idx);                        % reverse permute s_out
            s_out_org = obj.SD_detect(s_in, H, snr);
            err = (s_out ~= s_out_org);
            if sum(err(:) > 0)
                stop
            end
        end
    end
end