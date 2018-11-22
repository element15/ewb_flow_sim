clear;
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

% Data for sprinkler 1
p_data = (5:2.5:35) .* 12^2; % psf
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
if ~isempty(node_in.downstream_connections)
    for i = 1:length(node_in.downstream_connections)
        node_in.downstream_connections(i).head_in = node_in.head;
        node_in.downstream_connections(i) = ...
            update_link_pressure(node_in.downstream_connections(i));
    end
end
i = node_in.local_leaf;
if i ~= -1
    leaf_list{i}.head = node_in.head;
end

node_out = node_in;
end

% Given a link in the system, calculate the head loss in that link,
% including minor losses incurred at the upstream connection
function link_out = update_link_pressure(link_in)
% Note that there is no formal exit condition, because every link should be
% connected to nodes at both the upstream and downstream ends
global g;
link_in.head_out = max(0, link_in.head_in - link_in.delta_z - ...
    link_in.velocity^2 / 2 / g * head_loss_coefficient(link_in));
link_in.pressure_iteration = link_in.pressure_iteration + 1;
link_in.downstream_node.head = link_in.head_out;
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
        node_in.downstream_connections(i) = ...
            update_link_flow(node_in.downstream_connections(i));
    end
end
node_in.flow = flow;
node_out = node_in;
end

function link_out = update_link_flow(link_in)
link_in.velocity = link_in.downstream_node.flow / ...
    (pi/4 * link_in.diameter^2);
link_in.flow_iteration = link_in.flow_iteration + 1;
link_in.downstream_node = update_node_flow(link_in.downstream_node);
link_out = link_in;
end

function flow = get_leaf_flow(index)
global leaf_list;
l = leaf_list{index};
flow = l.flow_function(l.head);
end

% Given a node, generate a list including the local leaf (if any), as well
% as leaves connected to all downstream nodes.
function list = countl(start_node)
list = [];
if isempty(start_node)
    return;
end
if start_node.local_leaf ~= -1
    list = [start_node.local_leaf];
end
if ~isempty(start_node.downstream_connections)
    for link = start_node.downstream_connections
        list = cat(1, list, link.all_downstream_leaves);
    end
end
end

% Given a flow link, calculate the coefficient of head loss, or the
% coefficient of v^2/2g for the modified Bernoulli equation. This
% coefficient is given by f*L/D + K, and includes both major and 3K-method
% minor head losses
function c = head_loss_coefficient(link)
f = friction_factor(link.diameter, link.velocity);
if link.velocity == 0
    K = 0;
else
    K = link.K_1 / reynolds(link.diameter, link.velocity) + ...
        link.K_inf * (1 + link.K_d / link.diameter^0.3);
end
c = f * link.length / link.diameter + K;
end

% Given a nominal pipe size, return the inner diameter of the corresponding
% Schedule 40 PVC pipe. Returned diameter is given in *feet*
function d_inner = pvc(d_nominal)
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

% This function builds the full system structure, and returns (1) the root
% node of the system, `root`, as well as (2) a list of sprinkler heads
% (leaves) which are referenced by nodes in the system, `leaves`.
function [root, leaves] = init()
% Define the leaves first
leaf_count = 32;
leaves = cell(leaf_count, 1);
global a;
flow_func = @(v) a * v.^0.5;
for i = 1:leaf_count
    leaves{i} = flow_leaf();
    leaves{i}.flow_function = flow_func;
end

% Start from the bottom and work towards the root.
% NOTE: 4 in^0.3 = 1.9 ft^0.3 (units of k_d)
% All nodes modeled as tees, assuming flanged joins. Data from Neutrium 
% (https://neutrium.net/fluid_flow/pressure-loss-from-fittings-3k-method)
k_thru = [150, 0.05, 1.9];
k_turn = [800, 0.28, 1.9];

AN = flow_node([], 32);
AMAN = flow_link(AN, pvc('1/2'), 330, -160.08, countl(AN), k_thru);
AM = flow_node(AMAN, 31);
ALAM = flow_link(AM, pvc('1'), 100, 0, countl(AM), k_thru);
AL = flow_node(ALAM, 30);
AKAL = flow_link(AL, pvc('1'), 100, 0, countl(AL), k_thru);
AK = flow_node(AKAL, 29);
AHAK = flow_link(AK, pvc('1'), 100, -49.33, countl(AK), k_thru);

AU = flow_node([], 27);
AHAU = flow_link(AU, pvc('1/2'), 108.9, 0, countl(AU), k_turn);
AI = flow_node([], 28);
AHAI = flow_link(AI, pvc('1/2'), 108.9, 0, countl(AI), k_turn);

AH = flow_node([AHAU, AHAK, AHAI], -1);
AGAH = flow_link(AH, pvc('1'), 50, 0, countl(AH), k_thru);
AG = flow_node(AGAH, -1);
ADAG = flow_link(AG, pvc('1'), 50, -25.88, countl(AG), k_thru);

AF = flow_node([], 25);
ADAF = flow_link(AF, pvc('1/2'), 99, 0, countl(AF), k_turn);
AE = flow_node([], 26);
ADAE = flow_link(AE, pvc('1/2'), 99, 0, countl(AE), k_turn);

AD = flow_node([ADAF, ADAE, ADAG], -1);
AAAD = flow_link(AD, pvc('1'), 99, 0, countl(AD), k_thru);

AC = flow_node([], 23);
AAAC = flow_link(AC, pvc('1/2'), 89.1, 0, countl(AC), k_turn);
AB = flow_node([], 24);
AAAB = flow_link(AB, pvc('1/2'), 89.1, 0, countl(AB), k_turn);

AA = flow_node([AAAB, AAAC, AAAD], -1);
ZAA = flow_link(AA, pvc('1'), 49.5, -117.79, countl(AA), k_thru);
Z = flow_node(ZAA, 22);
WZ = flow_link(Z, pvc('1'), 49.5, 0, countl(Z), k_thru);

Y = flow_node([], 20);
WY = flow_link(Y, pvc('1/2'), 99, 0, countl(Y), k_turn);
X = flow_node([], 21);
WX = flow_link(X, pvc('1/2'), 99, 0, countl(X), k_turn);

W = flow_node([WX, WY, WZ], -1);
SW = flow_link(W, pvc('1'), 100, -20.88, countl(W), k_thru);
S = flow_node(SW, -1);
TS = flow_link(S, pvc('2'), 99, 0, countl(S), k_thru);
T = flow_node(TS, 19);
RT = flow_link(T, pvc('2'), 66, -18.29, countl(T), k_turn);

U = flow_node([], 18);
RU = flow_link(U, pvc('1/2'), 136, -43.21, countl(U), k_turn);
R = flow_node([RU, RT], -1);

PR = flow_link(R, pvc('2'), 36, -15.5, countl(R), k_thru);
Q = flow_node([], 17);
PQ = flow_link(Q, pvc('1/2'), 52, -15.5, countl(Q), k_turn);
P = flow_node([PQ, PR], -1);

NP = flow_link(P, pvc('2'), 124, -30.66, countl(P), k_thru);
O = flow_node([], 16);
NO = flow_link(O, pvc('1/2'), 61, 0, countl(O), k_turn);
N = flow_node([NO, NP], 15);
LN = flow_link(N, pvc('2'), 16, 0, countl(N), k_thru);

M = flow_node([], 14);
LM = flow_link(M, pvc('1/2'), 51, -1.25, countl(M), k_turn);
L = flow_node([LM, LN], -1);

KL = flow_link(L, pvc('2'), 84, 0, countl(L), k_turn);
n13 = flow_node([], 13);

% According to the layout diagram, this node flows uphill, as indicated by
% the positive delta-z.
PFB5 = flow_link(n13, pvc('1/2'), 55.8, 30.5, countl(n13), k_turn);

K = flow_node([PFB5, KL], -1);
GK = flow_link(K, pvc('2'), 14, -23.13, countl(K), k_turn);

J = flow_node([], 12);
IJ = flow_link(J, pvc('1'), 155, -24.88, countl(J), k_turn);
I = flow_node(IJ, 11);
HI = flow_link(I, pvc('1'), 95, 0, countl(I), k_turn);
H = flow_node(HI, 10);
GH = flow_link(H, pvc('1'), 61, -35.08, countl(H), k_turn);

G = flow_node([GH, GK], 9);
FG = flow_link(G, pvc('2'), 148, 0, countl(G), k_turn);
F = flow_node(FG, -1);
EF = flow_link(F, pvc('2'), 91, -38.62, countl(F), k_thru);

% Many of the lengths and delta-z values in this set of unsurveyed nodes
% and links may be extremely rough estimates. Unsurveyed points are
% designated by sprinkler number, rather than letter. Nodes given by a `Z`
% designation are non-sprinkler, non-survey-point nodes.
n8 = flow_node([], 8);
n7n8 = flow_link(n8, pvc('1'), 200, 0, countl(n8), k_thru);
n7 = flow_node(n7n8, 7);
% One half of PFB3
Z1n7 = flow_link(n7, pvc('1'), 269/2, -18.08/2, countl(n7), k_thru);
n6 = flow_node([], 6);
Z1n6 = flow_link(n6, pvc('1'), 20, 0, countl(n6), k_turn);
Z1 = flow_node([Z1n6, Z1n7], -1);
% Other half of PFB3
n5Z1 = flow_link(Z1, pvc('1'), 269/2, -18.08/2, countl(Z1), k_turn);
n5 = flow_node(n5Z1, 5);
En5 = flow_link(n5, pvc('1'), 20, 0, countl(n5), k_turn);

E = flow_node([En5, EF], -1);
DE = flow_link(E, pvc('2'), 89, -16.63, countl(E), k_thru);

% More non-survey sprinkler nodes
n4 = flow_node([], 4);
PFB2 = flow_link(n4, pvc('1'), 119, 0, countl(n4), k_thru);
n3 = flow_node(PFB2, 3);
Dn3 = flow_link(n3, pvc('1'), 20, 0, countl(n3), k_turn);

D = flow_node([Dn3, DE], -1);
CD = flow_link(D, pvc('2'), 35, -31.06, countl(D), k_thru);
C = flow_node(CD, -1);
BC = flow_link(C, pvc('2'), 34, 0, countl(C), k_thru);

% More non-survey sprinkler nodes
n2 = flow_node([], 2);
PFB1 = flow_link(n2, pvc('2'), 137, 0, countl(n2), k_thru);
n1 = flow_node(PFB1, 1);
Bn1 = flow_link(n1, pvc('2'), 20, 0, countl(n1), k_turn);

B = flow_node([BC, Bn1], -1);
AB = flow_link(B, pvc('3'), 124.4, 0, countl(B), k_thru);
A = flow_node(AB, -1);

% Connection to A from the upper tank `ut`
utA = flow_link(A, pvc('3'), 203.6, -50, countl(A), k_thru);
ut = flow_node(utA, -1);


root = ut;
end
