classdef flow_leaf
    properties
        head double = 0; % Head of the connected node
        flow_function % Function relating pressure to flow
    end
    methods
        function obj = flow_leaf(func)
            if nargin < 1
                return;
            end
            obj = flow_leaf;
            obj.flow_function = func;
        end
    end
end