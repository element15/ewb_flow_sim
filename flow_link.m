classdef flow_link
    properties
        % Intrensic parameters
        name                        string
        downstream_node             flow_node   % links connected to lower end
        diameter                    double      % pipe diameter, ft
        length                      double      % pipe length, ft
        delta_z                     double      % height differential end-over-end, ft
        all_downstream_leaves(:, 1) int8        % list of leaf indices
        % Primary state parameters
        pressure_iteration          uint8 = 0;  % last time pressure values were recalculated
        flow_iteration              uint8 = 0;  % last time flow values were recalculated
        velocity                    double = 0; % fluid velocity, ft/s
        head_in                     double = 0; % heat feet at inlet, ft
        % Secondary state parameters
        head_out                    double = 0; % heat feet at outlet, ft
        % Minor head loss at upstream end
        K                           double      % dimensionless
    end
    methods
        function obj = flow_link(name_in, dn, dia, l, dz, leaves, k)
            if nargin < 7
                return;
            end
            obj = flow_link;
            obj.name = name_in;
            obj.downstream_node = dn;
            obj.diameter = dia;
            obj.length = l;
            obj.delta_z = dz;
            obj.all_downstream_leaves = leaves;
            obj.K = k;
        end
    end
end