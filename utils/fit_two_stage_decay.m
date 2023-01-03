function [taus, gm, ir_fit] = fit_two_stage_decay(ir_env,t_env,t60_lims)

%% Fits two stage decay to an energy envelope of the form
% gm(1) + gm(2) * exp(-t/tau(1)) + gm(3) * exp(-t/tau(2))
% INPUTS : 
% ir_env - energy envelope of RIR
% t_env - time stamps corresponding to energy envelope
% t60_lims - 2x2 matrix of limits of two reverberation times (seconds)
% OUTPUTS :
% taus - 2x1 vector of decay slopes
% gm - 3x1 vector of amplitude and noise floor
% ir_fit - the fit energy decay curve 

tau_lims = t60_lims./log(1000);
%give more weight to the tail of the exponential
w = 1./ir_env;
lb = tau_lims(:,1);
ub = tau_lims(:,2);
tau0 = [lb(1);lb(2)];   %initial time constants in seconds
A = [1 -1];
b = 0;

function[err] = costfn(tau)
    
    E = [ones(length(t_env),1), exp(-t_env./tau(1)), exp(-t_env./tau(2))];
    gamma = (w.*E)\(w.*ir_env); %weighted least squares
    ir_hat = E*gamma;
    err = 0.5*sum((w.*(ir_env - ir_hat)).^2);
    
%      figure(3);grid on;
%      plot(t_env, 20*log10([ir_env, ir_hat]), 'Linewidth',2); hold off;
%      xlabel('time (s)'); ylabel('Magnitude (dB)');
%      xlim([0, 5]); ylim([-250 0]);
%      drawnow;
      
end

options = optimoptions('fmincon','Display','iter');
taus = fmincon(@costfn, tau0, A, b,[],[], lb, ub,[],options);
E = [ones(length(t_env),1), exp(-t_env./taus(1)), exp(-t_env./taus(2))];
gm = real((w.*E)\(w.*ir_env));
ir_fit = E*gm;

end


