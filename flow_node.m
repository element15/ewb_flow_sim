classdef flow_node
    properties
        downstream_connections(:, 1) flow_link
        local_leaf int8 % leaf index
        head double = 0 % Head feet at node (ft)
        flow double  = 0 % flow through the node
        head_limit double = -1 % pressure regulator installed if not -1
    end
    methods
        function obj = flow_node(dc, leaf)
            if nargin < 2
                return;
            end
            obj = flow_node;
            if ~isempty(dc)
                obj.downstream_connections = dc;
            end
            obj.local_leaf = leaf;
        end
    end
end