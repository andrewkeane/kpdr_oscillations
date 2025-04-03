function f = ocean_box_rhs(xx, par)

% par: F1, F2, sigma, tau, Ts

S1 = xx(1,1,:); % S1(t)
S1d = xx(1,2,:); % S1(t-\tau)
S2 = xx(2,1,:); % S2(t)

Ts = par(5); %T2 - T1;

S0 = 35; % psu
St = 3*S0;
S3 = St - S1 - S2;

V = 3.5e17; % m^3

Gamma = 365*24*60*60; % Use Gamma to rescale everything in terms of years
F1 = par(1)*1e6*Gamma; % Sv * 10^6 * Gamma = m^3 per year
F2 = par(2)*1e6*Gamma;
sigma = par(3)*1e6*Gamma;
k = 23*1e17; % m^3 per year
alpha_p = 1.7e-4; % per K
beta_p = 0.8e-3; % per psu

m = k/Gamma*(beta_p*(S2 - S1) - alpha_p*Ts)/1e6;

m_abs = abs(m); 
% This isn't smooth at psi=0, but really we're only concerned with psi>0

m_plus = Gamma*0.5*(m + m_abs)*1e6;
m_minus = Gamma*0.5*(m - m_abs)*1e6;

t_scale = 1/V;

%d/dt S1
f(1,1,:)= t_scale*(S0*F1 + m_plus.*(S2 - S1) + m_minus.*(S1 - S3)  + sigma*(S1d - S1));
%d/dt S2
f(2,1,:)= t_scale*(-S0*F2 + m_plus.*(S3 - S2) + m_minus.*(S2 - S1));
end

