classdef flow_node
    properties
        downstream_connections(:, 1) flow_link
        local_leaf int % leaf index
        head double % Head feet at node (ft)
        flow double % flow through the node
    end
end