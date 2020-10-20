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

Initial_load_percent = 100; %Initial load in percent of initial load
End_load_percent  = 17; %Final load fraction of initial load
step_size_big = 1; %Big step size percentage of initial load
step_size_small = 0.1; %Small step size percentage of inital load
break_percent = 20; %Breakpoint percentage for big step size
unload    = false;

% D) Define various iteration quanteties
%Allocatin memory
K         = zeros(2*length(dof));                 %Stiffness matrix
f_int     = zeros(2*length(dof),1);               %Internal force vector
eps_his   = zeros(2*length(dof),4);     %Strain history
eps       = zeros(length(enod),4);           %Current strain
sigma_old = zeros(4,length(enod));
plasticitycheck = zeros(1,length(enod));
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
current_load_percent = Initial_load_percent;
step_nbr = 1;


while current_load_percent > End_load_percent
    tic
    if current_load_percent <= break_percent
        step_size = step_size_small/100;
    else
        step_size = step_size_big/100;
    end
        
    disp(' ')
    disp(['Force %: ', num2str(current_load_percent)])
    current_load_percent = current_load_percent - step_size*100;
    
    % G) Unload and compute residual
    f = (1-step_size)*f;
    res = f - f_int;           %check that out of balance force is zero
    res(bc(:,1)) = 0;
    res = max(abs(res));

    
    % H) Newton loop
    while res>tol
%         if mod(iteration,1) == 0
%             disp(['Residual: ', num2str(res)])
%         end
        
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
            
            if dlambda ~= 0
                plasticitycheck(el) = 1;
            end
            
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
%         if mod(iteration,1) == 0
%             disp(['new Residual: ', num2str(res)])
%         end
    end
    step_time(step_nbr) = toc;
    step_nbr = step_nbr + 1;
    % Q) Accept and save quantities as an equilibrium state
end

%% ------------------Plot step time-----------------------------------------
n = linspace(1,length(step_time),length(step_time));
hold on;
grid on;
plot(n,step_time,'*')

%% ------------------ PLOT -------------------------------------------------
    
%
figure('Renderer', 'painters', 'Position', [400 100 800 600])
[ex,ey] = coordxtr(edof,coord,dof,3);
eldraw2(ex,ey,[1 2 0]); %green
%eldraw2(ex,-ey,[1 2 0]);
% eldraw2(-ex+0.29,ey, [1 2 0]);
% eldraw2(-ex+0.29,-ey, [1 2 0]);

eldisp2(ex,ey,ed,[1 4 0], 1); %red
% eldisp2(ex,-ey,ed,[1 4 0], 1); %red
% eldisp2(-ex+0.29,ey,-ed,[1 4 0], 1); %red
% eldisp2(-ex+0.29,-ey,-ed,[1 4 0], 1); %red
legend('Reference configuration');
axis equal
title('Displacement field (m), disp controlled')

%% Plot von Mises
vMises = zeros(length(enod),1);   %Von Mises stress
for el=1:length(enod)
    vMises(el) = stress_invariant_J2(sigma_old(:,el));
end
figure('Renderer', 'painters', 'Position', [400 100 800 600])
[ex,ey] = coordxtr(edof,coord,dof,3);
nnod = length(coord(:,1));
eff_node = zeros(nnod,1);
for node = 1:nnod
    [c0,c1] = find(enod==node);
    eff_node(node,1) = sum(vMises(c0)/size(c0,1));
end
enodtemp = [(1:length(enod))',enod];
eff_field = extract(enodtemp,eff_node(:));
fill(ex', ey', eff_field');
title('Von Mises effective stressfield [N/m^2]')
colorbar;

%% Plot volumetric stress
for el=1:length(enod)
    I1plot(el) = stress_invariant_I1(sigma_old(:,el));
end
figure('Renderer', 'painters', 'Position', [400 100 800 600])
[ex,ey] = coordxtr(edof,coord,dof,3);
nnod = length(coord(:,1));
eff_node = zeros(nnod,1);
for node = 1:nnod
    [c0,c1] = find(enod==node);
    eff_node(node,1) = sum(I1plot(c0)/size(c0,1));
end
enodtemp = [(1:length(enod))',enod];
eff_field = extract(enodtemp,eff_node(:));
fill(ex', ey', eff_field');
title('Volumetric stress field [N/m^2]')
colorbar;

%% Plot plasticity check

figure('Renderer', 'painters', 'Position', [400 100 800 600])
[ex,ey] = coordxtr(edof,coord,dof,3);
ex1=[];
ex2=[];
ey1=[];
ey2=[];
for el=1:length(enod)
    if plasticitycheck(el) == 1
        ex1 = [ex1; ex(el,:)];
        ey1 = [ey1; ey(el,:)];
    else
        ex2 = [ex2; ex(el,:)];
        ey2 = [ey2; ey(el,:)];
    end
end

eldraw2(ex1,ey1,[1 2 0]); %green

eldraw2(ex2,ey2,[1 4 0]); %red

%legend('Reference configuration');
%axis equal
%title('Displacement field (m), disp controlled')

%%



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
% end

