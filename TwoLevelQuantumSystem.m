
clear, clc
close all
guardar = 1; % Cambiar a 1 para guardar figuras


%% 1. Constantes
hbar = 1.054571817e-34;
e_carga = 1.602176634e-19;
ao = 0.0529177e-9;
w = 1;
wo = 1; 
deg = 2*e_carga*ao;
Eo = 10e-6;
phi = 0;

% Adimensionalización

wbar = w/wo;
% lambda = (deg*Eo)/(hbar*wo);
lambda = 0.01;


% CONDICIONES INICIALES

ce0 = 0 + 0i;
cg0 = 1 + 0i;
psi0 = [ce0; cg0];
tauspan = [0 30]; 
dtau = 1e-5;

Omega = lambda * wo;  
Delta = w - wo;

%% 2. Representaciones gráficas

% RÉGIMEN IDEAL

lambda = 0.01;
wbar = 1; 
Omega = lambda * wo;  
Delta = (wbar*wo) - wo; 

tauspan = [0 (16*pi)/sqrt(Omega.^2 + Delta.^2)]; 
[tau, ~, ce2, cg2] = exact_sol(tauspan,dtau,psi0,lambda,wbar,phi);
[Pe_RWA, Pg_RWA] = RWA_sol(psi0, tau, lambda, Delta, wo);

fig = Represent_sol(tau, ce2, cg2, Pe_RWA, Pg_RWA, wo, wbar, phi, Omega, Delta, "R\'egimen ideal de la RWA", 1, 'regimenideal', guardar);


% AUMENTO PROGRESIVO DE INTENSIDAD

lambda_values = [0.05, 0.1, 0.2, 0.3];
wbar = 1;
Delta = (wbar*wo) - wo;

for nvalue = 1:length(lambda_values)
    lval = lambda_values(nvalue);
    Omega = lval * wo;
    
    tauspan = [0 (2*pi)/sqrt(Omega.^2 + Delta.^2)]; 
    [tau, ~, ce2, cg2] = exact_sol(tauspan,dtau,psi0,lval,wbar,phi);
    [Pe_RWA, Pg_RWA] = RWA_sol(psi0, tau, lval, Delta, wo);
    
    fig = Represent_sol(tau, ce2, cg2, Pe_RWA, Pg_RWA, wo, wbar, phi, Omega, Delta, 'Aumento progresivo de la intensidad', 1, sprintf("intensidad%d",nvalue), guardar);
end


% ESTUDIO DEL CORRIMIENTO BLOCH-SIEGERT
lambda = 0.1;
Omega = lambda * wo;  
wbar_values = 0.99:0.0001:1;
maxtransf = 0;
wbarmaxtransf = 0;

for nvalue = 1:length(wbar_values)
    wbval = wbar_values(nvalue);
    Delta = (wbval * wo) - wo;
    
    tauspan = [0 (2*pi)/sqrt(Omega.^2 + Delta.^2)]; 
    [~, ~, ce2, ~] = exact_sol(tauspan,1e-4,psi0,lambda,wbval,phi);

    if max(ce2) > maxtransf
        maxtransf = max(ce2);
        wbarmaxtransf = wbval;
    end
    
end

% Obtenemos los valores para la máxima transferencia
Delta = (wbarmaxtransf * wo) - wo;

tauspan = [0 (2*pi)/sqrt(Omega.^2 + Delta.^2)]; 
[tau, ~, ce2, cg2] = exact_sol(tauspan,dtau,psi0,lambda,wbarmaxtransf,phi);


fprintf("\n\nLa máxima transferencia se obtiene en Δ = %.4f y es %f\n\n", wbarmaxtransf,maxtransf)

fig = Represent_sol(tau, ce2, cg2, Pe_RWA, Pg_RWA, wo, wbarmaxtransf, phi, Omega, Delta, "Gráfica de máxima transferencia", 2, sprintf("maximatransferencia",nvalue), guardar);


% FUERA DE RESONANCIA

lambda = 0.05;
Omega = lambda * wo;  
wbar_values = [0.8, 1.2, 1.5];

for nvalue = 1:length(wbar_values)
    wbval = wbar_values(nvalue);
    Delta = (wbval * wo) - wo;
    
    tauspan = [0 (2*pi)/sqrt(Omega.^2 + Delta.^2)]; 
    [tau, ~, ce2, cg2] = exact_sol(tauspan,dtau,psi0,lambda,wbval,phi);
    [Pe_RWA, Pg_RWA] = RWA_sol(psi0, tau, lambda, Delta, wo);
    
    fig = Represent_sol(tau, ce2, cg2, Pe_RWA, Pg_RWA, wo, wbval, phi, Omega, Delta, 'Fuera de resonancia', 1, sprintf("fueraderesonancia%d",nvalue), guardar);
end


% RÉGIMEN ULTRA-FUERTE

lambda = 0.6;
wbar = 1;
Delta = (wbar*wo) - wo;
Omega = lambda * wo;  

tauspan = [0 (2*pi)/sqrt(Omega.^2 + Delta.^2)]; 
[tau, psi, ce2, cg2] = exact_sol(tauspan,dtau,psi0,lambda,wbar,phi);
[Pe_RWA, Pg_RWA] = RWA_sol(psi0, tau, lambda, Delta, wo);

fig = Represent_sol(tau, ce2, cg2, Pe_RWA, Pg_RWA, wo, wbar, phi, Omega, Delta, 'R\''egimen ultra-fuerte', 1, 'ultrafuerte', guardar);

disp("Fin.")


%% --- FUNCIONES ---

function [tau, psi, ce2, cg2] = exact_sol(tauspan,dtau, psi0, lambda, wbar, phi)

options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8, MaxStep = dtau*abs(tauspan(2) - tauspan(1)) );

fprintf("\nCalculando solución numérica (wbar = %.2e, lambda = %.2e)...\n",wbar,lambda)

[tau, psi] = ode45(@(t,y) edocc(t,y,lambda,wbar,phi), tauspan, psi0, options);
ce2 = abs(psi(:,1)).^2;
cg2 = abs(psi(:,2)).^2;

end

function [Pe_RWA, Pg_RWA] = RWA_sol(psi0, tau, lambda, Delta, wo)

    Omega = lambda * wo;  

    Omega_R = sqrt(Omega.^2 + Delta.^2);

    arg = Omega_R .* (tau ./ wo) ./ 2;
    
    fprintf("Calculando RWA (Omega=%.2e, Delta=%.2e)\n", Omega, Delta);
    
    
    Pe_RWA = abs(1i .* (Omega ./ Omega_R) .* psi0(2) .* sin(arg) + psi0(1) .* ( cos(arg) - 1i .* (Delta ./ Omega_R) .* sin(arg) )).^2;
    
    Pg_RWA = abs( psi0(2) .* ( cos(arg) + 1i .* (Delta ./ Omega_R) .* sin(arg) ) + 1i .* (Omega ./ Omega_R) .* psi0(1) .* sin(arg) ).^2;
end

function representation = Represent_sol(tau, ce2, cg2, Pe_RWA, Pg_RWA, wo, wbar, phi, Omega, Delta, titl, type, savetitle, save)
%
% type (int): parámetro que determina si se representa solo la parte de
%   1: solución numérica y RWA
%   2: solución numérica
%   3: solución RWA
%

    representation = figure('Units','centimeters','Position',[5 5 18 14]); 
    
    % Configuración global de fuentes para esta figura
    set(representation, 'DefaultAxesFontSize', 12);
    set(representation, 'DefaultTextInterpreter', 'latex');
    set(representation, 'DefaultLegendInterpreter', 'latex');
    set(representation, 'DefaultAxesTickLabelInterpreter', 'latex');
    
    sgtitle(titl, 'Interpreter','latex', 'FontSize', 16, 'FontWeight', 'bold');

    subplot(2,1,1)
    hold on
    if type == 1 || type == 2
        plot(tau, ce2, '-', 'LineWidth', 2, 'Color',[0 0.5 0]) 
        plot(tau, cg2, '--', 'LineWidth', 2, 'Color',[0.61 0.31 0.02]) 
    end
    if type == 1 || type == 3
        plot(tau, Pe_RWA, ':k', 'LineWidth', 1.5); 
        plot(tau, Pg_RWA, ':r', 'LineWidth', 1.5); 
    end
    hold off
    ylabel('Probabilidad', 'FontSize', 14)
    grid on
    box on
    xlim([min(tau) max(tau)])
    ylim([-0.05 1.05])
    
    if type == 1
        legend("$|c_e|^2$ Num","$|c_g|^2$ Num","$P_{e}$ RWA","$P_{g}$ RWA", ...
            'Location', 'eastoutside', 'FontSize', 11)
    elseif type == 2
        legend("$|c_e|^2$","$|c_g|^2$", 'Location', 'eastoutside', 'FontSize', 11)
    elseif type == 3
        legend("$P_{e}$ RWA","$P_{g}$ RWA", 'Location', 'eastoutside', 'FontSize', 11)
    end
    
    lambda_val = Omega/wo;
    str = sprintf('$\\lambda=%.2f, \\delta=%.2f, \\Delta=%.2f, \\Omega=%.2f$', lambda_val, Delta, wbar,Omega);
    title(str, 'FontSize', 12, 'FontWeight', 'normal')

    subplot(2,1,2)
    hold on
    if type == 1 || type == 2
        plot(tau, ce2 - cg2, '-', 'LineWidth', 2, 'Color',[0.5 0 0]) 
    end
    if type == 1 || type == 3
        plot(tau, Pe_RWA - Pg_RWA, '--b', 'LineWidth', 1.5);
    end
    hold off
    xlabel('$\tau$', 'FontSize', 14)
    ylabel('Inversi\''on ($P_e - P_g$)', 'FontSize', 14)
    grid on
    box on
    xlim([min(tau) max(tau)])
    ylim([-1.1 1.1])
    
    if type == 1
        legend("Num\'erica", "RWA", 'Location', 'eastoutside', 'FontSize', 11)
    elseif type == 2
        legend("Num\'erica", 'Location', 'eastoutside', 'FontSize', 11)
    elseif type == 3
        legend("RWA", 'Location', 'eastoutside', 'FontSize', 11)
    end

    if save == 1
        exportgraphics(representation, strcat(savetitle, '.png'), 'Resolution', 350);
    end
end

% Sistema de EDOs

function dcdt = edocc(tau,psi,lambda,wbar,phi)
dcdt = zeros(2,1);
dcdt(1) = (-1i) * ( 0.5*psi(1) - lambda * cos(wbar*tau + phi)*psi(2) );
dcdt(2) = (-1i) * (-lambda * cos(wbar*tau + phi) * psi(1) - 0.5 * psi(2) );
end