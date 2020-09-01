
function [E_repaired] = close_meshpolygon(E_initial, P)
%   Given an unclosed polygon defined by vertex list P and 
%   connection list E_initial, close it by matching singly-referenced
%   vertices based on distance between them.

%   Copyright WAW/SNM 2020

    E1 = E_initial(:,1);
    E2 = E_initial(:,2);

    E1E2 = sort([E1; E2]);

    E1_sorted = E1E2(mod(1:length(E1E2), 2)==1);
    E2_sorted = E1E2(mod(1:length(E1E2), 2)==0);

    %Find the first vertex that is not repeated
    nonrepeated_idx = find(E1_sorted ~= E2_sorted, 1);
    v = zeros(0, 1);
    while( ~isempty(nonrepeated_idx) && nonrepeated_idx <= length(E1_sorted) - length(v))

        if(E1_sorted(nonrepeated_idx) < E2_sorted(nonrepeated_idx)) %If E1 contains the unrepeated vertex
            v(end+1) = E1_sorted(nonrepeated_idx);
            %Clean out the first nonrepeated vertex and place it at the end of its respective list
            E1_sorted(nonrepeated_idx) = [];
            E1_sorted(end+1) = v(end);

        else
            v(end+1, 1) = E2_sorted(nonrepeated_idx);
            E2_sorted(nonrepeated_idx) = [];
            E2_sorted(end+1) = v(end);
        end
        nonrepeated_idx = find(E1_sorted ~= E2_sorted, 1);
    end

    if(isempty(v))
        E_repaired = E_initial;
        return;
    end
    
    %Find the vertices corresponding to entries of v that are nearest each other and pair them
    [v_unique, ~, ic] = unique(v);
    required_connections = accumarray(ic, 1); %Number of times each element of v_unique occurs
    
    c = [];
    while(~isempty(v_unique))
        [~, start] = min(required_connections);
        v_temp = v_unique(start);
        nearest_idx = knnsearch(P(v_unique,:), P(v_temp,:), 'K', 2);
        nearest_idx = nearest_idx(2);
        
        %Connect v_temp and v(nearest_idx)
        c = [c; v_temp, v_unique(nearest_idx)];
        %Clear out one possible connection for each of the involved elements
        required_connections(start) = required_connections(start) - 1;
        required_connections(nearest_idx) = required_connections(nearest_idx) - 1;
        v_unique(required_connections == 0) = [];
        required_connections(required_connections == 0) = [];
    end
    %idx = knnsearch(P(v,:), P(v,:), 'K', 2);
    %idx = idx(:,2);
    %v = [v v(idx)];
    %Delete duplicate rows
%     v = sort(v, 2);
%     for j = length(v):-1:2
%         if(any(v(j-1,1) == v(j,1)))
%             v(j,:) = [];
%         end
%     end

    %Append v to EofXY
    %E_repaired = [E_initial; v];
    E_repaired = [E_initial; c];

end