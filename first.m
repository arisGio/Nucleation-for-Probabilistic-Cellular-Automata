%6  24 85
plot(0:8000,beta_spec(0.5),'g',0:8000,beta_spec(0.55),'r',0:8000,beta_spec(0.57),'b')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% FUNCTION DEFINITIONS %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MCMC FOR SPECIFIC TEMPERATURE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function magn = beta_spec(beta)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% SIMULATION PARAMETERS %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    L=250;
    h=0.4;
    kappa=1;
    n=1;
    T=8000; %total time
    %beta=0.5;  %inverse temperature
    spins = zeros(L,L,T+1); %all the spin configurations
    z = zeros(L,L,T); %random variables matrix
    probab = zeros(L,L,T); %probabilities matrix
    rng('shuffle');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% INITIAL CONFIGURATION %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    spins(:,:,1) = -ones(L);
    e_spins = zeros(L+2,L+2,T+1); %extended spin matrix including PBCs

    for t = 1:T

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% EXTENDED SPIN CONFIGURATION %%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        e_spins(1,2:L+1,t) = spins(1,:,t); %up
        e_spins(L,2:L+1,t) = spins(L,:,t); %down
        e_spins(2:L+1,1,t) = spins(:,1,t); %left
        e_spins(2:L+1,L,t) = spins(:,L,t); %right
        e_spins(2:L+1,2:L+1,t) = spins(:,:,t); %filling the inside
        %e_spins
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%% CALC OF TRANSITION PROBABILITY FOR ALL SITES %%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for i = 2:L+1
            for j = 2:L+1
                probab(i-1,j-1,t) = t_prob(i, j, t, e_spins, n, beta, kappa, h);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% PERIODIC BOUNDARY CONDITIONS ON RANDOM VARIABLE Z %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&
        z(:,:,t)=rand(L);
        z(1,:,t) = z(L,:,t);
        z(:,1,t) = z(:,L,t);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%% UPDATE SPIN CONFIGURATION %%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        spins(:,:,t+1) = sign(probab(:,:,t)-z(:,:,t));

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% PLOTTING SPIN CONFIGURATIONS %%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %imagesc(spins(:,:,5));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% TIME EVOLUTION - MAGNETIZATION %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for y = 1:(T+1)
        magn(y) = sum(spins(:,:,y),'all')/L^2;
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% TRANSITION PROBABILITY %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = t_prob(i, j, t, e_spins, n, beta, kappa, h)
    % VON - NEUMANN NEIGHBORS SUM
    n_sum = e_spins(i-1,j,t) + e_spins(i+1,j,t) + e_spins(i,j+1,t) + e_spins(i,j-1,t);
    % TRANSITION PROBABILITY
    f = (1+n*tanh(beta*(kappa*e_spins(i,j,t)+n_sum+h)))/2;
end

