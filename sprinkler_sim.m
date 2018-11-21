% Defined values
global mu;
global rho;
global gamma;
global rough;
global g;
global a;
% Properties of water at 40 deg F
mu    = 32.34e-6; % viscosity, lbf*s/ft^2
rho   =  1.94   ; % density, slug/ft^3
gamma = 62.43   ; % spc. gravity, lbf/ft^3
rough = 20e-6   ; % upper bound roughness for PVC, ft
g     = 32.2    ; % acceleration of gravity, ft/s^2


global leaf_list;
leaf_list = [];



% Data for sprinkler 1
p_data = 5:2.5:35; % psi
v_dot_data = [6.83, 7.77, 8.75, 9.53, 10.30, 10.82, 11.87, 12.29, ...
    12.84, 13.18, 13.53, 14.70, 15.53] ./ 12^3; % ft^3/s
[a, R2] = square_root_fit(p_data, v_dot_data);
v_dot = @(p) a .* p.^0.5; % Function for the curve fit

global leaf_list;
[root, leaf_list] = init();




% Given a node in the system, update the inbound head of every connected
% link, and call update_link_pressure() on each of those links, which will
% recursively solve all downstream parts of the system. Note that minor
% head losses due to node geometry are handled by the affected downstream
% link.
function node_out = update_node_pressure(node_in)
global leaf_list;
% Exit condition
if isempty(node_in.downstream_connections)
    node_out = node_in;
    return;
end

for i = 1:length(node_in.downstream_connections)
    node_in.downstream_connections(i).head_in = node_in.head;
    node_in.downsteram_connections(i) = ...
        update_link_pressure(node_in.downstream_connections(i));
end
i = node_in.local_leaf;
if i ~= -1
    leaf_list(i).head = node_in.head;
end

node_out = node_in;
end

% Given a link in the system, calculate the head loss in that link,
% including minor losses incurred at the upstream connection
function link_out = update_link_pressure(link_in)
% Note that there is no formal exit condition, because every link should be
% connected to nodes at both the upstream and downstream ends
global g;
link_in.head_out = link_in.head_in - link_in.delta_z - ...
    link_in.velocity^2 / 2 / g * head_loss_coefficient(link_in);
link_in.pressure_iteration = link_in.pressure_iteration + 1;
link_in.downstream_node = update_node_pressure(link_in.downstream_node);

link_out = link_in;
end

function node_out = update_node_flow(node_in)
flow = 0;
ind = node_in.local_leaf;
if ind ~= -1
    flow = flow + get_leaf_flow(ind);
end
if ~isempty(node_in.downstream_connections) % Exit condition
    for i = 1:length(node_in.downstream_connections)
        for j = node_in.downstream_connections(i).all_downstream_leaves(:)
            flow = flow + get_leaf_flow(j);
        end
        update_link_flow(node_in.downstream_connections(i));
    end
end
node_in.flow = flow;
node_out = node_in;
end

function link_out = update_link_flow(link_in)
link_in.velocity = link_in.downstream_node.flow / ...
    (pi/4 * node_in.diameter^2);
link_in.flow_iteration = link_in.flow_iteration + 1;
update_node_flow(link_in.downstream_node);
link_out = link_in;
end

function flow = get_leaf_flow(index)
global leaf_list;
flow = leaf_list(index).flow_function(leaf_list(index).head);
end

% Given a flow link, calculate the coefficient of head loss, or the
% coefficient of v^2/2g for the modified Bernoulli equation. This
% coefficient is given by f*L/D + K, and includes both major and 3K-method
% minor head losses
function c = head_loss_coefficient(link)
f = friction_factor(link.diameter, link.velocity);
K = link.K_1 / reynolds(link.diameter, link.velocity) + ...
    link.K_inf * (1 + link.K_d / link.diameter^0.3);
c = f * link.length / link.diameter + K;
end

% Given a nominal pipe size, return the inner diameter of the corresponding
% Schedule 40 PVC pipe. Returned diameter is given in *feet*
function d_inner = pvc_diameter(d_nominal)
d = 0;
switch d_nominal
    case '1/2'
        d = 0.622;
    case '3/4'
        d = 0.824;
    case '1'
        d = 1.049;
    case '1 1/4'
        d = 1.380;
    case '1 1/2'
        d = 1.610;
    case '2'
        d = 2.067;
    case '2 1/2'
        d = 2.469;
    case '3'
        d = 3.068;
    case '4'
        d = 4.026;
end
if d == 0
    excep = MException('SprinklerSim:PVCLookupError', ...
        ['An error occured looking up the inner diameter for the ' ...
        'specified pipe size `%s`.'], d_nominal);
    throw(excep);
end
d_inner = d / 12; % Convert from inches to feet
end

% Given a pipe diameter and fluid velocity, calculate the friction factor
% for major head loss. This function uses values for roughness, fluid
% density, and viscosity defined in the file header, as these are values
% which are considered to be constant throughout the system.
% NOTE: Diameter and velocity are given in ft and ft/s, respectively
function f = friction_factor(diameter, velocity)
global rough;
re = reynolds(diameter, velocity);
% Haaland approximation
f = (-1.8 * log10( ...
    (rough / diameter / 3.7)^1.11 + (6.9 / re) ...
    ))^-2;
end

% Calculate the Reynolds number for a given velocity and diameter, and
% using global values of density and viscosity
function re = reynolds(diameter, velocity)
global rho;
global mu;
re = rho * velocity * diameter / mu;
end

% Given a dataset (x_v, y_v), find the least squares curve fit of the form
% y = a * sqrt(x), and also return the R^2 value for the curve fit. This
% curve fit assumes that the coefficient `a` is strictly positive.
function [a_out, R2_out] = square_root_fit(x_v, y_v)
a_lo = 1e-3; % Initial lower bound
a_hi = 1e+3; % Initial upper bound

y  = @(a, x) a .* x.^0.5; % This equation is the form desired
y_bar = sum(y_v) / length(y_v); % Average of the input data
R2 = @(a) 1 - ... % Computes statistical R2 value for a given `a`
    sum((y_v - y(a, x_v)).^2) / ...
    sum((y(a, x_v) - y_bar).^2);
mid = @(u, v) (u + v) / 2; % midpoint

for ind = 1:100
    a_mid = mid(a_lo, a_hi);
    if R2(a_lo) > R2(a_hi)
        a_hi = a_mid;
    else
        a_lo = a_mid;
    end
end

a_val = mid(a_lo, a_hi);
R2_val = R2(a_val);

% Ensure that the values for `a` and `R2` are valid
if R2_val < 0 || R2_val > 1 || a_val < 0
    excep = MException('SprinklerSim:CurveFitError', ...
        ['An unidentified error has occured in `square_root_fit()`, ' ...
        'resulting in values for `a` or `R2` that are invalid: ' ...
        'a = %f, R2 = %f'], a_val, R2_val);
    throw(excep);
end

% Return results
a_out = a_val;
R2_out = R2_val;
end





















































