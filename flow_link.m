classdef flow_link
    properties
        % Intrensic parameters
        downstream_node             flow_node   % links connected to lower end
        diameter                    double      % pipe diameter, ft
        length                      double      % pipe length, ft
        delta_z                     double      % height differential end-over-end, ft
        all_downstream_leaves(:, 1) int         % list of leaf indices
        % Primary state parameters
        pressure_iteration          int = 0;    % last time pressure values were recalculated
        flow_iteration              int = 0;    % last time flow values were recalculated
        velocity                    double = 0; % fluid velocity, ft/s
        head_in                     double = 0; % heat feet at inlet, ft
        % Secondary state parameters
        head_out                    double = 0; % heat feet at outlet, ft
        % Minor head loss at upstream end (3K Method)
        K_1                         double      % dimensionless
        K_inf                       double      % dimensionless
        K_d                         double      % ft^0.3
    end
    methods
        function obj = flow_link(dn, dia, l, dz, leaves, k1, kinf, kd)
            obj = flow_link;
            obj.downstream_node = dn;
            obj.diameter = dia;
            obj.length = l;
            obj.delta_z = dz;
            obj.all_downstream_leaves = leaves;
            obj.K_1 = k1;
            obj.K_inf = kinf;
            obj.K_d = kd;
        end
    end
end