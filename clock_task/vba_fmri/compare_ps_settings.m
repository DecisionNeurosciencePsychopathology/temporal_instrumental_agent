%basedir='/Users/michael/ics/temporal_instrumental_agent/clock_task/vba_fmri/vba_out';
basedir='/storage/group/mnh5174_collab/temporal_instrumental_agent/clock_task/vba_fmri/vba_out';
rundirs=dir(basedir);
%first 2 elements are always . and .. -- remove these
rundirs(1:2) = [];

models=cell(1,length(rundirs));
energy=NaN(76, length(rundirs));
for d = 1:length(rundirs)
   res = extract_sceptic_fmri([basedir, '/', rundirs(d).name]);
   energy(:,d) = res.pars(:,2);
   models{d} = rundirs(d).name;
end

%notes: fmriold, replicate_prefix, and sigmamax1p0, breaknew are all identical
%these use the old h_sceptic code with a .08 PS starting point, freely estimated prop_spread, but
%actually use prop_spread as sigma of Gaussian, leading to a max SD of 1.0. Also, refspread is not re-estimated
%but instead the refspread of .08 Gaussian is used.
%breaknew uses the h_sceptic code, but tweaks the top piece to match 

addpath(genpath('/storage/home/mnh5174/MATLAB/VBA-toolbox'));

%multisession constrain 0125 is now in the second position
%try it compared to other good contstrained models
%VBA_groupBMC(energy(:,[1:5])') %.0125, .0125 multisession, .025
VBA_groupBMC(energy(:,[1 2 10 11 12])') %.0125, .0125 multisession, .025
VBA_groupBMC(energy(:,[2 10])') %.0125, .0125 multisession


VBA_groupBMC(energy')
VBA_groupBMC(energy(:,[1 2])') %.0125 versus .025
VBA_groupBMC(energy(:,[2 3])') %.025 versus .03
VBA_groupBMC(energy(:,[2 4])') %.025 versus .04
VBA_groupBMC(energy(:,[2 5])') %.025 versus .05
VBA_groupBMC(energy(:,[2:4])') %all constrained

%compare fixed models: .0125 fixed wins
VBA_groupBMC(energy(:,[9:11])') %all constrained
VBA_groupBMC(energy(:,[2 9])') %.025 free versus .0125 fixed
VBA_groupBMC(energy(:,[1:3 9 10])') %.025 free versus .0125 fixed



VBA_groupBMC(energy(:,[1:4, 6])')
VBA_groupBMC(energy(:,1:5)')
VBA_groupBMC(energy(:,5:6)')
VBA_groupBMC(energy(:,[2 5])')

VBA_groupBMC(energy(:,[1 7])')

[a, r] = VBA_groupBMC(energy(:,[2 3 8])') %.025 versus old de facto sigma=1.0, ps=.025 version (but the latter has the refspread defect)


ediff = energy(:,2) - energy(:,3);
figure; hist(ediff)

find(ediff > 40)