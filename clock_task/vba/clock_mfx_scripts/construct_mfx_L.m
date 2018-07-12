function L=construct_mfx_L

%Find and load all o_sub mat files
file_list = glob('*o_sub*');

L=[]; %Initialize output variable
for i = 1:length(file_list)
   load(file_list{i}) %Load in the data
   L(i,:)=cellfun(@(x) x.F,o_sub); %Pull log model evidence to L struct
   stop=1;
end

%Save the data
save('model_read_in_order.mat','file_list')
save('mfx_L_struct', 'L')
