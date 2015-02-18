function [cost] = fitnessfcn(params)
    cost = multirun_clock_logistic_operator(20, 888, {params, 0, 'IEV', 100});
end
