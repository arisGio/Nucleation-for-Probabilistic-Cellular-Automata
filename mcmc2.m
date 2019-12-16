%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% NO SELF-INTERACTION WITH STAGGERED MAGNETIZATION %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    tic % inital time

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% SIMULATION PARAMETERS %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    L=350; % lattice-size
    h=0.2;
    kappa=0; % self-interaction
    n=1;
    T=8000; %total time
    
    beta=0.5;  %inverse temperature
    
    rng('shuffle');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% INITIAL CONFIGURATION, PREALLOCATION %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    spins = -ones(L);
    e_spins = zeros(L+2,L+2); % extended spin matrix including PBCs
    magn = zeros(1,T+1); % magnetization initial value
    z = zeros(L,L); % random variables matrix
    probab = zeros(L,L); % probabilities matrix

    for t = 1:T+1

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% EXTENDED SPIN CONFIGURATION (PBCs) %%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        e_spins(1,2:L+1) = spins(L-1,:); %up
        e_spins(L,2:L+1) = spins(2,:); %down
        e_spins(2:L+1,1) = spins(:,L-1); %left
        e_spins(2:L+1,L) = spins(:,2); %right
        e_spins(2:L+1,2:L+1) = spins(:,:); %filling the inside
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%% CALC OF TRANSITION PROBABILITY FOR ALL SITES %%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        for i = 2:L+1
            for j = 2:L+1
                probab(i-1,j-1) = t_prob(i, j, e_spins, n, beta, kappa, h);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% PERIODIC BOUNDARY CONDITIONS ON RANDOM VARIABLE Z %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        z = rand(L);
        z(1,:) = z(L,:);
        z(:,1) = z(:,L);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%% TIME EVOLUTION - STAGGERED MAGNETIZATION %%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        for x = 1:L
            for y = 1:L
                magn(t) = magn(t) + (-1)^(x+y)*spins(x,y)/L^2;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%% UPDATE SPIN CONFIGURATION %%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        spins(:,:) = sign(probab(:,:)-z(:,:));

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% PLOTTING SPIN CONFIGURATIONS %%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %imagesc(spins(:,:))
    
    %plot(0:T,magn,'g',0:T,magn,'r',0:T,magn,'b')
    
    plot(0:T,magn,'g')

    toc % final time

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% FUNCTION DEFINITIONS %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% TRANSITION PROBABILITY %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function f = t_prob(i, j, e_spins, n, beta, kappa, h)
        % VON - NEUMANN NEIGHBORS SUM
        n_sum = e_spins(i-1,j) + e_spins(i+1,j) + e_spins(i,j+1) + e_spins(i,j-1);
        % TRANSITION PROBABILITY
        f = (1+n*tanh(beta*(kappa*e_spins(i,j)+n_sum+h)))/2;
    end

