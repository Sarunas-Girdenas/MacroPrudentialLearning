

% instatiate SteadyState solution class

SteadySLN = SteadyStateSolution;

% define parameters

SteadySLN.epsilon   = 11;
SteadySLN.v         = 0.03;
SteadySLN.gamma     = 0.98;
SteadySLN.beta      = 0.99;
SteadySLN.m_star_st = 0.92;
SteadySLN.j         = 0.1;
SteadySLN.eta       = 1.01;
SteadySLN.theta     = 0.75;
SteadySLN.m         = 0.89;
SteadySLN.r_pi      = 0.27;
SteadySLN.chi       = 1; %problematic

% solve for steady state and return struct of values.
% to access the steady state values, use 
% steadyStateValues.VariableName

steadyStateValues = SteadySLN.solveSteadyState();