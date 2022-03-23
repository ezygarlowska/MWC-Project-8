function Dr = dissRoller(beta,Er,c,g)
%DISSROLLER Compute roller dissipation.
%   DR = DISSROLLER(BETA, ER, C) computes the roller dissipation DR
%   using the slope of the wave front BETA, the roller energy ER, the
%   wave celerity C, and the gravitational acceleration G. DISSROLLER
%   is based on the equations: (1) Dr = c * tau, where tau is the
%   shear stress between roller and wave interface; (2) tau = (rho * g
%   * A * sin(beta)) / L, where rho is the density of water, g the
%   gravitational acceleration, L the wavelength, and A the roller
%   area; (3) A = (Er * 2 * L) / (rho * c^2).

Dr = 2*sin(beta)*g*Er/c;
