classdef flow_node
    properties
        downstream_connections(:, 1) flow_link
        local_leaf int % leaf index
        head double = 0 % Head feet at node (ft)
        flow double  = 0 % flow through the node
    end
    methods
        function obj = flow_node(dc, leaf)
            obj = flow_node;
            obj.downstream_connections = dc;
            obj.local_leaf = leaf;
        end
    end
end