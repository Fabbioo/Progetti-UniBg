
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Autore: Fabio Filippo Mandalari.                      %
%                       Matricola: 1047426.                               %
%       Corso: CAM (Controllo Avanzato Multivariabile), UniBg.            %
%           Docenti: Prof. Antonio Ferramosca, Ing. Marco Polver.         %
%   Progetto relativo al controllo di un servomeccanismo mediante MPC.    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Funzione.

function x_dot = MyServomechanism(t, x, u, K_el, D_el, J_m, J_l, n, tauf_segnato)

x_dot = zeros(4,1);

x_dot(1) = x(3);

x_dot(2) = x(4);

x_dot(3) = ...
    - (K_el / J_m) * x(1) ...
    + ((n * K_el) / J_m) * x(2) ...
    - (D_el / J_m) * x(3) - (tauf_segnato / J_m) * tanh(x(3)) ...
    + ((n * D_el) / J_m * x(4)) ...
    + (1 / J_m) * u(1);

x_dot(4) = ...
    + ((n * K_el) / J_l) * x(1) ...
    - ((n^2 * K_el) / J_l) * x(2) ...
    + ((n * D_el) / J_l) * x(3) ...
    - ((n^2 * D_el) / J_l) * x(4) ...
    - (1 / J_l) * u(2);

end