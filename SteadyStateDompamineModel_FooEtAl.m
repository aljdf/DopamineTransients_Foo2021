% This MATLAB script accompanies Foo et al., (2021) Current Biology, 
% "Reinforcement learning links spontaneous cortical dopamine impulses to reward"
% Written Monday June 21st 2021 by Johnatan Aljadeff (email: aljadeff@ucsd.edu)

% The code produces a fit to the steady state model of learning dopamine 
% impluse rate and delay (relative to reward). The model is fit to days 2
% and 3 of the experiment.

clear ;

load('C:\Users\Admin\Documents\MATLAB\DopamineTransients\Data\d2_animals_lags_ntransients_May21.mat') ;
r_data = [0.00942857142857143, 0.0162857142857143; ...
          0.00571428571428571, 0.0100000000000000; ...
          0.0108571428571429 , 0.0171428571428571; ...
          0.0291428571428571 , 0.0188571428571429; ...
          0.00914285714285714, 0.0105714285714286; ...
          0.0137142857142857 , 0.0222857142857143; ...
          0.00400000000000000, 0.0100000000000000; ...
          0.0422857142857143 ,0.0134285714285714 ] ; 
% Average dopamine impluse rates in Hz for 8 animals over 2 days.
% One animal was excluded from this analysis, see manuscript for details. 

d_data = [ 1.88300000000000  ,-6.11800000000000  ; ...
           2.35300000000000  ,-9.64600000000000  ; ...
          -4.23600000000000  ,-0.470000000000000 ; ...
           0.236000000000000 ,-2.35300000000000  ; ...
           3.29400000000000  ,-10.1141000000000  ; ...
           0.236000000000000 ,-1.64680000000000  ; ...
           2.11700000000000  ,-6.58800000000000  ; ...
           0.234999999999999 ,-1.88200000000000] ;
% Average delay in seconds of dopamine impluses relative to reward delivery
       
[n_animal,n_day] = size(r_data) ;
n_data = n_animal * n_day ; 


% fixed parameters: 
% (see manuscript for explanation on how parameters were chosen) 
tau     = 30   ; % decay timescale of dopamine transients (seconds)
Delta0  = 5    ; % baseline value of delay
Delta1  = -21  ; % change in delay due to threshold crossing (seconds)
lambda0 = 0.01 ; % baseline value of impulse rate 

n_lambda1 = 100 ;                      % number of values of lambda1 where model is evaluated
lambda1x = linspace(0,0.2,n_lambda1) ; % values of lambda1 where model is evaluated

n_Deltaw = 100 ;                       % number of values of Deltaw where model is evaluated
Deltawx = linspace(0.1,10,n_Deltaw) ;  % values of Deltaw where model is evaluated

n_theta   = 1000 ;                      % number of values of theta where model is evaluated     
theta    = logspace(-7,3,n_theta) ;    % values of theta where model is evaluated

f_Delta    = @(r,th) Delta0 + Delta1*gammainc(th,tau*r,'upper') ;
                                       % steady state value of Delta (see Eq. 7 in manuscript)

SqErr_r = zeros(n_lambda1,n_Deltaw,n_data) ;
                                       % squared error for each data point and each combination of parameters

for jr1 = 1:n_lambda1                  % loop over lambda1 values
    for jdwin = 1:n_Deltaw             % loop over Deltaw  values
        lambda_i = zeros(n_theta,1) ;  % lambda as a function of theta 
        Delta_i  = zeros(n_theta,1) ;  % Delta as a function of theta
        f_lambda = @(r,th) lambda0 + lambda1x(jr1)*gammainc(th,tau*r,'upper').*exp(-abs(f_Delta(r,th)/Deltawx(jdwin))) ;
        % implicit equation for lambda as a function of itself and theta
        % (see Eq. 8 in manuscript)
        
        for jth = 1:n_theta
            fz_lambda = @(lambda) lambda-f_lambda(lambda,theta(jth)) ; 
            lambda_i(jth) = fzero(fz_lambda,[lambda0,1]) ;
            % solving implicit equation for lambda as a function of itself and theta (Eq. 8)
            Delta_i(jth) = f_Delta(lambda_i(jth),theta(jth)) ;
            % solving equation for Delta as a function of lambda (Eq. 7)
        end
        [Delta_i,i] = unique(Delta_i) ;
        lambda_i = lambda_i(i) ;
        SqErr_r(jr1,jdwin,:) = (interp1(Delta_i,lambda_i,d_data(:))-r_data(:)).^2 ;
        % evaluation of squared error of each data point.
        % relies on interpolation because of gap in solution (see text for details)
    end
end

%%
[minSqErr,jfit] = min(sum(SqErr_r,3),[],'all','linear') ; 
[jfit_r1,jfit_dwin] = ind2sub([n_lambda1,n_Deltaw],jfit) ;
% finding indices of model parameters corresponding to minimal squared error

lambda_fit = zeros(n_theta,1) ;
Delta_fit  = zeros(n_theta,1) ;

f_lambda = @(r,th) lambda0 + lambda1x(jfit_r1)*gammainc(th,tau*r,'upper').*exp(-abs(f_Delta(r,th)/Deltawx(jfit_dwin))) ;

for jth = 1:n_theta
    fz_lambda        = @(r) r-f_lambda(r,theta(jth)) ;
    lambda_fit(jth) = fzero(fz_lambda,[lambda0 1]) ;
    Delta_fit(jth)  = f_Delta(lambda_fit(jth),theta(jth)) ;
    % solving for lambda as a function of Delta
end
[~,jj] = max(abs(diff(Delta_fit))) ; % finding gap in solutions (see text for details)

% plotting 
h(3) = plot(Delta_fit,lambda_fit,'-','Color',[0.5 0 0.5],'LineWidth',10) ; hold on ;
h(4) = plot(Delta_fit(jj:(jj+1)),lambda_fit(jj:(jj+1)),'-','Color',[1 0 1],'LineWidth',10) ; hold on ;
h(1) = plot(d_data(:,1),r_data(:,1),'ow','MarkerFaceColor','r','MarkerSize',16) ; hold on ;
h(2) = plot(d_data(:,2),r_data(:,2),'ow','MarkerFaceColor','b','MarkerSize',16) ; hold on ;
legend(h,{'Day 2','Day 3','Best fitting model','Gap in model solutions'}) ;
xlabel('lag (s)') ;
ylabel('pulse rate (1/s)') ;
title(['\lambda_1 = ' num2str(lambda1x(jfit_r1),3) 'Hz, \Delta_w = ' num2str(Deltawx(jfit_dwin),3) 's']) ;