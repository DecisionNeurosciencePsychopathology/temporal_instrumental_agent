function a = initialize_agents_struct()
%This defines the simulation parameters for all agents to be tested for optimality.
%These are stored in a struct where the field name is the model naem
clear a;

% NB: a struct array, despite its elegance, is causing many headaches in parfor loops since agents(a) indexing is invalid
% unless a is the parfor loop counter. This is not the best way to parallelize (since there are only 4-6 agents).
% revert to single struct and use fieldnames.

%MODEL 1: KALMAN SKEPTIC
% Parameters: epsilon (proportion reduction in uncertainty), sig_spread (interval over which to generalize PE)
skeptic.init_params = [.6 .08]; %proportion reduction in uncertainty, proportion spread for temporal generalization
skeptic.lower_bounds = [.1 .01];
skeptic.upper_bounds = [.99 .25];
skeptic.k = length(skeptic.init_params); %number of free parameters
skeptic.uvsum = 0;
skeptic.name = 'skeptic'; %add explicit name parameter since each field in struct gets pulled out separately
a(1) = skeptic;

%MODEL 2: Q LEARNING
% qlearning.init_params = [0.9 0.2 0.08 0.99];
% qlearning.lower_bounds = [0.8 0.01 0.01 0.90];
% qlearning.upper_bounds = [0.99 0.35 0.30 0.999];
% qlearning.k = length(qlearning.init_params); %number of free parameters
% qlearning.name = 'qlearning';
% a(2) = qlearning;
% 
% %MODEL 3: SARSA
% sarsa.init_params = [0.9 0.2 0.08 0.99];
% sarsa.lower_bounds = [0.8 0.01 0.01 0.90];
% sarsa.upper_bounds = [0.99 0.35 0.30 0.999];
% sarsa.k = length(sarsa.init_params); %number of free parameters
% sarsa.name = 'sarsa';
% a(3) = sarsa;
% 
% %MODEL 4: FRANK TC
% franktc.init_params = [ 0.2, 3000, 0.3, 0.3, 1000, 0.1, 300 ];
% franktc.lower_bounds = [ 0, 0, 0.01, 0.01, 1, 0, 0 ];
% franktc.upper_bounds = [1, 100000, 5, 5, 5000, 5000, 10000 ];
% franktc.k = length(franktc.init_params); %number of free parameters
% franktc.name = 'franktc';
% a(4) = franktc;
% 
% skeptic_uvsum.init_params = [.6 .08]; %proportion reduction in uncertainty, proportion spread for temporal generalization
% skeptic_uvsum.lower_bounds = [.1 .01];
% skeptic_uvsum.upper_bounds = [.99 .25];
% skeptic_uvsum.k = length(skeptic.init_params); %number of free parameters
% skeptic_uvsum.uvsum = 1;
% skeptic_uvsum.name = 'skeptic_uvsum';
% a(5) = skeptic_uvsum;

end