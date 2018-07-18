function vba_df=create_vba_mfx_input_clock()

sample=1; %1 for n = 76, 2 for SPECC, 3 for EXPLORE

% models = {'fixed', 'fixed_decay', 'fixed_uv', 'kalman_logistic', 'kalman_softmax', 'kalman_uv_sum'};
models = {'fixed_decay'};

for model = models
    
    %Initialize the table
    vba_df = table();
    
    %Use glob to pull all the vba files
    % AD: placed all the vba output in
    % '/Volumes/bek/datafolderfromjonEdrive/sceptic/vba_out/specc_clock/'
    % to preserve directory structure
    if sample==1
        data_path = '/Volumes/bek/datafolderfromjonEdrive/sceptic/vba_out/n76/';
        vba_files = glob([data_path '*' model{:} '_multinomial1*']);
        study = 'n76';
    elseif sample==2
        data_path = '/Volumes/bek/datafolderfromjonEdrive/sceptic/vba_out/specc_clock/';
        vba_files = glob([data_path '*' model{:} '.mat']);
        study = 'specc';
    elseif sample==3
        data_path = '/Volumes/bek/datafolderfromjonEdrive/sceptic/vba_out/explore_clock/';
        vba_files = glob([data_path '*' model{:} '.mat']);
        study = 'explore';
    end
    
    %File loop
    for vba_file = vba_files'
        load(vba_file{:}, 'out') %Load in the file should contain b,out,posterior
        id = regexp(vba_file{:}, '[0-9]+', 'match');
        if sample==1
        id = str2double(id{2});
        else
        id = str2double(id{:});
        end
        %Initialize temporay dataframe
        tmp_table = table();
        
        %Grab id
        tmp_table.ID = id;
        
        %Grab y & u
        %tmp_table.y = {out.y};
        %tmp_table.u = {out.u};
        
        %Grab the options used or perhaps create a sub function to create this
        tmp_table.options = {out.options};
        tmp_table.y = {out.y};
        tmp_table.u = {out.u};
        vba_df = [vba_df; tmp_table];
    end
    
    %Save the data
    name = sprintf('%s/vba_mfx_input/',data_path );
    if 7~=exist(name,'dir')
        mkdir(name)
    end
    save([data_path '/vba_mfx_input/' sprintf('vba_mfx_input_%s_%s',model{:},study)], 'vba_df');
end