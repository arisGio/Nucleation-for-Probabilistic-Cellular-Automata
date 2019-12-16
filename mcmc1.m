%%%%%%%%%%%%%%% MCMC WITH SELF-INTERACTION %%%%%%%%%%%%%%%%%%

tic

%rng('shuffle');


%[ magn, st_magn, spins ] = beta_specific(L,h,kappa,n=1,T,beta)


[a1,b1,c1]=beta_specific(250,0.4,1,1,8000,0.5); % temperature #1

[a2,b2,c2]=beta_specific(250,0.4,1,1,8000,0.55); % temperature #2

[a3,b3,c3]=beta_specific(250,0.4,1,1,8000,0.57); % temperature #3

%%%%%%%%%%%%%%%%%%%%% MAGNETIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
plot(1:length(a1),a1,'g',1:length(a2),a2,'r',1:length(a3),a3,'b');

%%%%%%%%%%%%%%%%%%%%% STAGGERED - MAGNETIZATION %%%%%%%%%%%%%%%%%%%

%plot(1:length(b1),b1,'g',1:length(b2),b2,'r',1:length(b3),b3,'b');

%%%%%%%%%%%%%%%%%%%%% SPINS CONFIGURATION %%%%%%%%%%%%%%%%%%%%%%%%%

figure(2)
imagesc(c1)
figure(3) 
imagesc(c2)
figure(4)
imagesc(c3)

toc