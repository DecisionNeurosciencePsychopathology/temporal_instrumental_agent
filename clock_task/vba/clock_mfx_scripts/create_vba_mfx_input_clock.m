function vba_df=create_vba_mfx_input_clock()

specc_flag=0; %1 for fMRI sample 0 for behavioral

models = {'fixed', 'fixed_decay', 'fixed_uv', 'kalman_logistic', 'kalman_softmax', 'kalman_uv_sum'};

for model = models
    
    %Initialize the tbale
    vba_df = table();
    
    %Use glob to pull all the vba files
    if specc_flag
        data_path = 'E:/data/sceptic/vba_out/specc_clock/';
        vba_files = glob([data_path '*' model{:} '.mat']);
        study = 'specc';
    else
        data_path = 'E:/data/sceptic/vba_out/explore_clock/';
        vba_files = glob([data_path '*' model{:} '.mat']);
        study = 'explore';
    end
    
    %File loop
    for vba_file = vba_files'
        load(vba_file{:}, 'out') %Load in the file should contain b,out,posterior
        id = regexp(vba_file{:}, '[0-9]+', 'match');
        id = str2double(id{:});
        
        %Initialize temporay dataframe
        tmp_table = table();
        
        %Grab id
        tmp_table.ID = id;
        
        %Grab y & u
        %tmp_table.y = {out.y};
        %tmp_table.u = {out.u};
        
        %Grab the options used or perhaps create a sub function to create this
        tmp_table.options = {out.options};
        
        vba_df = [vba_df; tmp_table];
    end
    
    %Save the data
    save([data_path 'vba_mfx_input/' sprintf('vba_mfx_input_%s_%s',model{:},study)], 'vba_df');
end