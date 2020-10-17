% A) Define material and model parameters

%Model
R0     = 8;   %Radius, 8m
w      = 42;  %Distance from wall to drill
h      = 42;  %Distance from roof to drill
t      = 1;   %Thickness 
ptype  = 2;   % planes strain = 2 (plan stress = 1)
ep     = [ptype t];
%Material
G      = 300; %Young modulus, GPa
nu     = 0.3; %Poisson ratio
gamma  = 30*pi/180;  %How to intrepet? Which angle? Need to change to radians
A      = 0.0067;
B      = 48.2;

mp     = [G,nu,gamma,A,B];

% B) Define inital elastic state

%{
Displacement boundary conditions  bc
External force (pressure)         f
Initial displacements             a
Element degrees of freedom        edof
Element nodes                     enod
Element x-y-coordinates           coord
Nodal degrees of freedom          dof
%}

load initial_state.mat

% C) Define unloading parameters
tol       = 1e-3;
NbrSteps  = 100; %Number of steps
unload    = false;

% D) Define various iteration quanteties
%Allocatin memory
a         = zeros(length(dof),1);               %Displacement vector
K         = zeros(length(dof),1);                 %Stiffness matrix
f         = zeros(length(dof),1);               %External force vector
f_int     = zeros(length(dof),1);               %Internal force vector
%eps_his   = zeros(nbr_elem,3*NbrSteps);     %Strain history
eps       = zeros(length(enod),3);           %Current strain

% E) Build initial tangent and internal force

% for i=1:length(enod)
%     %Element coordinates
%     ex = coord(enod(i,:),1)'; 
%     ey = coord(enod(i,:),2)';
%     
%     % Insert my elastic D, create function
%     Ke = plante(ex,ey,ep,Dt);
%     indx = edof(i,2:end);
%     K(indx,indx) = K(indx,indx)+Ke;
% end
