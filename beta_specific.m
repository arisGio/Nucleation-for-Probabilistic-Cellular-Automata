%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MARKOV CHAIN MONTE CARLO (MCMC) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% FUNCTION DEFINITIONS %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    function [ magn, st_magn, spins ] = beta_specific(L,h,kappa,n,T,beta)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% SIMULATION PARAMETERS %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % L is the lattice-size
        % h is the magnetic field
        % kappa is the self-interaction
        % n is the assumed spin value
        % T is the total time

        rng('shuffle');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% INITIAL CONFIGURATION, PREALLOCATION %%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        spins = -ones(L); % -1 spin configuration
        e_spins = zeros(L+2,L+2); % extended spin matrix
        magn = zeros(1,T+1); % magnetization
        st_magn = zeros(1,T+1); % staggered magnetization
        z = zeros(L,L); % random variables matrix
        probab = zeros(L,L); % probabilities matrix

        for t = 1:T+1

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%% EXTENDED SPIN CONFIGURATION (PBCs) %%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            e_spins(1,2:L+1) = spins(L,:); % up
            e_spins(L+2,2:L+1) = spins(1,:); % down
            e_spins(2:L+1,1) = spins(:,L); % left
            e_spins(2:L+1,L+2) = spins(:,1); % right
            e_spins(2:L+1,2:L+1) = spins(:,:); % filling the inside

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%% CALC OF TRANSITION PROBABILITY FOR ALL SITES %%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            for i = 2:L+1
                for j = 2:L+1
                    probab(i-1,j-1) = t_prob(i, j, e_spins, n, beta, kappa, h);
                end
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% CREATING RANDOM VARIABLES MATRIX %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            z = rand(L);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%% TIME EVOLUTION - MAGNETIZATION %%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            for x = 1:L
                for y = 1:L
                    magn(t) = magn(t) + spins(x,y)/L^2;
                end
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%% TIME EVOLUTION - STAGGERED MAGNETIZATION %%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            for x = 1:L
                for y = 1:L
                    st_magn(t) = st_magn(t) + (-1)^(x+y)*spins(x,y)/L^2;
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%% UPDATE SPIN CONFIGURATION %%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if t<T+1
                spins = sign(probab-z);
            end

        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% TRANSITION PROBABILITY %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function f = t_prob(i, j, e_spins, n, beta, kappa, h)
        % VON - NEUMANN NEIGHBORS SUM
        n_sum = e_spins(i-1,j) + e_spins(i+1,j) + e_spins(i,j+1) + e_spins(i,j-1);
        % TRANSITION PROBABILITY
        f = (1+n*tanh(beta*(kappa*e_spins(i,j)+n_sum+h)))/2;
    end

