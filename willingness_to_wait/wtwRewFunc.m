function prr=wtwRewFunc(trial)

% ct = 1:trial;

pars = {{0, 1},{.4, .575, 0},{0.25, 0.25}, {5, 1}};
conds = {'unif', 'gp', 'beta', 'beta'};
output_names = {'unif', 'gp', 'beta', 'late_beta'};

prr=[];
%pev=[];



for c = 1:length(conds)
    v = (random(conds{c}, pars{c}{:}, 1, trial));
    
    prr = setfield(prr, output_names{c}, v);
    %pev = setfield(pev, output_names{c}, e_v);
end


z = camel_distrib(trial);
prr = setfield(prr, 'camel', z);
save prr
end



function z = camel_distrib(trial)

sig = .05;
x = random('normal',1/4,sig,1,trial/2);
y = random('normal',3/4,sig,1,trial/2);

z = 1:trial;
rp = randperm(trial);
z(1:length(x))=x;
z(length(x)+1:end)=y;
z = z(rp);


return 
end





% 
%             prr.gp_t10(i) = random(conds{cond},pars{cond}{:});
% 
%             prr.unif_t10(i) = random(cond{cond},pars{cond}{:});
% 
%             prr.late_beta_t10(i) = random(conds(cond),pars{cond}{:});
