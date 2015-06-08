for i = 1:5
    for o = 1:100
    s.(agents{i}).evAll(1,o) = sum(s.(agents{i}).eV(1,o).ev_i); %IEV
    s.(agents{i}).evAll(2,o) = sum(s.(agents{i}).eV(2,o).ev_i); %DEV
    s.(agents{i}).evAll(3,o) = sum(s.(agents{i}).eV(3,o).ev_i); %QUADUP
    end
end



for i = 1:5
    s.(agents{i}).mean_costEV(1,1) = mean(s.(agents{i}).evAll(1,:));
    s.(agents{i}).mean_costEV(2,1) = mean(s.(agents{i}).evAll(2,:));
    s.(agents{i}).mean_costEV(3,1) = mean(s.(agents{i}).evAll(3,:));
end


%Reversals....

s_tmp = s;
save s_tmp s_tmp

for i = 1:5
    for o = 1:100
    s.(agents{i}).evIEVtoDEV(o) = sum(s.(agents{i}).eV(4,o).ev_i);
    end
end

for i = 1:5
    s.(agents{i}).mean_costIevtoDev_EV = mean(s.(agents{i}).evIEVtoDEV);
end


s_tmp = s;
save s_tmp s_tmp

for i = 1:5
    for o = 1:100
    s.(agents{i}).evDEVtoIEV(o) = sum(s.(agents{i}).eV(5,o).ev_i);
    end
end

for i = 1:5
    s.(agents{i}).mean_costDevtoIev_EV = mean(s.(agents{i}).evDEVtoIEV);
end



%Debug this is to double check to make run every agent has 200 data points
%per run, as it seems this is notthe case...
for i = 1:5
    for k = 1:5
        idx=1;
        for j = 1:100
            if length(s.(agents{i}).eV(k,j).ev_i) <200
                bad_run_index{k,idx} = {i k j}; %The agent, the row, the column
                idx = idx+1;
            end
        end
    end
end