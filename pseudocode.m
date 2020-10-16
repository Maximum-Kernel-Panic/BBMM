%Pseudocode

%%
close all, clear all;

% A) Define material and model parameters

%{
Displacement boundary conditions  bc
External force (pressure)         f
Initial displacements             a
Element degrees of freedom        edof
Element nodes                     enod
Element x-y-coordinates           coord
Nodal degrees of freedom          dof
%}

% B) Define inital elastic state
load initial_state.mat

% C) Define unloading parameters

% D) Define various iteration quanteties

% E) Build initial tangent and internal force

%% F) Unloading loop (force controlled)
for load_step=1:number_of_steps
    disp(' ')
    disp(['Load step number: ', num2str(load_step)])
    
    % G) Unload and compute residual
    
    % H) Newton loop
    while max(abs(res))>1e-4
        disp(['Residual: ', num2str(load_step)])
        
        % I) Solve for displacement increment and update
        
        % J) Update tangent and internal force
        for el=1:nelm
            
            % K) Compute current total strain and strain increment
            %   (requires to have access to previous equlibrium state)
            
            % L) Update plastic variables, (check for plasticity)
            [sigma,dlambda,ep_eff] = ...
            update_variables(sigma_old,ep_eff_old,delta_eps,Dstar,mp)
            
            % M) Compute element algorithmic tangent, D_ats
            Dats = alg_tan_stiff(sigma,dlambda,ep_eff,Dstar,mp);
            
            % N) Compute element internal forces and stiffness matrix

            % O) Assemble element matrices to global matrices
            
        end
        
        % P) Update residual
        
    end
    
    % Q) Accept and save quantities as an equilibrium state

end

%% R) Plot results
    
    % Note that some aditional information from the simulation need be
    % extracted and saved during the unloading loop to be able to create
    % all of these plots.
    
    % 1) Plot percentual radius decrease to wall pressure

    % 2) Plot load path of element closest to cone tip

    % 3) Plot tanB vs plastic effective strain of element closest to cone tip

    % 4) Plot Original and deformed mesh

    % 5) Plot von Mises effective stress in the body

    % 6) Plot I_1 stress in the body

    % 7) Plot plastic elements