
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Autore: Fabio Filippo Mandalari.                      %
%                       Matricola: 1047426.                               %
%       Corso: CAM (Controllo Avanzato Multivariabile), UniBg.            %
%           Docenti: Prof. Antonio Ferramosca, Ing. Marco Polver.         %
%   Progetto relativo al controllo di un servomeccanismo mediante MPC.    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Funzione.

% Input:
% _ sys_d     -> sistema linearizzato e discretizzato;
% _ Hx, hx    -> Matrice e vettore relativi ai vincoli sugli stati;
% _ Hu, hu    -> Matrice e vettore relativi ai vincoli sugli ingressi;
% _ x_ref_lin -> Riferimenti degli stati nelle coordinate del sistema linearizzato;
% _ u_ref_lin -> Riferimenti degli ingressi nelle coordinate del sistema linearizzato;
% _ Q         -> Matrice per il peso degli stati;
% _ R         -> Matrice per il peso degli ingressi;
% _ N         -> Orizzonte di predizione.

% Output:
% _ mpc       -> Struttura dati contenente tutte le informazioni che servono
%                per poter effettuare un controllo mediante MPC.

function mpc = designMPC(sys_d, Hx, hx, Hu, hu, x_ref_lin, u_ref_lin, Q, R, N)

    % Obiettivo: calcolo di tutti gli elementi utili per la corretta
    % formulazione del problema di programmazione quadratica.

    %% 1. Definizione della dinamica del sistema.
    
    % Calcolo i movimenti libero e forzato dello stato (pag. 12 exe. 4)
    Acalligrafica = calcola_MovimentoLibero(sys_d.A, N); % Movimento libero
    Bcalligrafica = calcola_MovimentoForzato(sys_d.A, sys_d.B, N); % Movimento forzato
    
    %% 2. Trasformazione dei vincoli.
    
    % Calcolo le matrici Htilde_u e htilde_u (pag. 13 exe. 4)
    Htilde_u = calcola_Hxutilde(Hu, N-1);
    htilde_u = calcola_hxutildee(hu, N-1);

    % Calcolo le matrici Htilde_x e htilde_x (pag. 14 exe. 4)
    Htilde_x = calcola_Hxutilde(Hx, N);
    htilde_x = calcola_hxutildee(hx, N);

    % Definizione dei vincoli (pag. 17 exe. 4)
    A_qp = [Htilde_u; Htilde_x * Bcalligrafica];
    b_cost = Htilde_x * Acalligrafica;
        
    %% 3. Trasformazione del costo.
    
    % Calcolo il peso dello stato finale (soluzione dell'equazione di Riccati)
    [~, P] = dlqr(sys_d.A, sys_d.B, Q, R);
    
    % Calcolo le matrici Q_tilde e R_tilde (pag. 18 exe. 4)
    Qtilde = calcola_Qtilde(Q, P, N);
    Rtilde = calcola_Rtilde(R, N);
    
    H = Rtilde + Bcalligrafica' * Qtilde * Bcalligrafica;
    H_qp = 2 * H;

    xR = calcola_xuRefLin(x_ref_lin, N);
    uR = calcola_xuRefLin(u_ref_lin, N-1);

    f1_cost = 2 * Acalligrafica' * Qtilde * Bcalligrafica;
    f2 = - 2 * xR' * Qtilde * Bcalligrafica - 2 * uR' * Rtilde;

    %% Costruzione MPC

    mpc = struct( ...
        'Acalligrafica', Acalligrafica, ...
        'Bcalligrafica', Bcalligrafica, ...
        'Htilde_u', Htilde_u, 'htilde_u', htilde_u, ...
        'Htilde_x', Htilde_x, 'htilde_x', htilde_x, ...
        'Qtilde', Qtilde, ...
        'Rtilde', Rtilde, ...
        'H_qp', H_qp, ...
        'f1_cost', f1_cost, ...
        'f2', f2, ...
        'A_qp', A_qp, ...
        'b_cost', b_cost ...
    );

end
