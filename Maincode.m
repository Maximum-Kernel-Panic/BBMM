% A) Define material and model parameters
clear all
clc
%Model
R0     = 8;   %Radius, 8m
w      = 42;  %Distance from wall to drill
h      = 42;  %Distance from roof to drill
t      = 1;   %Thickness 
ptype  = 2;   % planes strain = 2 (plan stress = 1)
ep     = [ptype t];
%Material
G      = 3e8; %Young modulus, GPa
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
NbrSteps  = 75; %Number of steps
unload    = false;

% D) Define various iteration quanteties
%Allocatin memory
K         = zeros(2*length(dof));                 %Stiffness matrix
f_int     = zeros(2*length(dof),1);               %Internal force vector
eps_his   = zeros(2*length(dof),4);     %Strain history
eps       = zeros(length(enod),4);           %Current strain
sigma_old = zeros(4,length(enod));
%sigma     = zeros(length(enod),4);



ep_eff_old = 0;
ep_eff     = 0;
dlambda    = 0;
% E) Build initial tangent and internal force
Dstar = elastic_tan_stiff(mp);

for i=1:length(enod)
    %Element coordinates
    ex = coord(enod(i,:),1)'; 
    ey = coord(enod(i,:),2)';
    
    % Insert my elastic D, create function
    Ke = plante(ex,ey,ep,Dstar);
    indx = edof(i,2:end);
    K(indx,indx) = K(indx,indx)+Ke;
end

f_int = f;
%% F) Unloading loop (force controlled)

Dats = Dstar;
for load_step=1:NbrSteps

    disp(' ')
    disp(['Load step number: ', num2str(load_step)])
    
    % G) Unload and compute residual
    f = 0.99*f;
    res = f - f_int;           %check that out of balance force is zero
    res(bc(:,1)) = 0;
    res = max(abs(res));
    
    % H) Newton loop
    while res>tol
        disp(['Residual: ', num2str(res)])
        
        % I) Solve for displacement increment and update
        da = solveq(K, f-f_int, bc );
        a  = a+da;
        epshistory = eps;
        ed = extract(edof,a);
        % J) Update tangent and internal force
        K       = zeros(2*length(dof));           %Reset stiffness matrix
        f_int   = zeros(2*length(dof),1);     %Reset internal force vector
        for el=1:length(enod)
            
            % K) Compute current total strain and strain increment
            %   (requires to have access to previous equlibrium state)
            % Element coordinates
            ex = coord(enod(el,:),1)';
            ey = coord(enod(el,:),2)';
 
            % Update strains
            [~,eps_e] = plants(ex, ey, ep, Dats, ed(el,:));
            eps(el,:) = eps_e;
            delta_eps = eps(el,:) - epshistory(el,:);
            
            % L) Update plastic variables, (check for plasticity)
            [sigma,dlambda,ep_eff] = ...
            update_variables(sigma_old(:,el),ep_eff_old,delta_eps',Dstar,mp);
            sigma_old(:,el) = sigma;
            % M) Compute element algorithmic tangent, D_ats           
            Dats = alg_tan_stiff(sigma,dlambda,ep_eff,Dstar,mp);
            
            % N) Compute element internal forces and stiffness matrix
            Ke      = plante(ex,ey,ep,Dats);
            f_int_e = plantf(ex, ey, ep, sigma');
            
            % O) Assemble element matrices to global matrices
            indx         = edof(el,2:end);
            K(indx,indx) = K(indx,indx)+Ke;
            f_int(indx)  = f_int(indx)+f_int_e;
            
        end
        
        % P) Update residual
        res = f - f_int;           %check that out of balance force is zero
        res(bc(:,1)) = 0;
        res = max(abs(res));
        disp(['new Residual: ', num2str(res)])
    end
    
    % Q) Accept and save quantities as an equilibrium state
    
end
