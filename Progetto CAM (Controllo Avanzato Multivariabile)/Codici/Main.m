
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Autore: Fabio Filippo Mandalari.                      %
%                       Matricola: 1047426.                               %
%       Corso: CAM (Controllo Avanzato Multivariabile), UniBg.            %
%           Docenti: Prof. Antonio Ferramosca, Ing. Marco Polver.         %
%   Progetto relativo al controllo di un servomeccanismo mediante MPC.    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 0. Setup ambiente di lavoro.

clear;
clc;
close all;
format rational;

set(0,'DefaultLineLineWidth', 1.5);
set(0,'defaultAxesFontSize', 14);
set(0,'DefaultFigureWindowStyle', 'docked');

addpath Libreria Simulazioni;

%% 1. Setup delle strutture.

% Obiettivo della sezione: definizione di tutte le variabili utili per
% l'impostazione del problema da risolvere.

syms x1 x2 x3 x4 u1 u2 u3 u4 J_m J_l K_el D_el tauf_segnato n;

% Vettore contenente tutti gli stati che compongono il sistema
stati = [x1; x2; x3; x4];
ordineStati = size(stati,1);

% Vettore contenente tutti gli ingressi del sistema
ingressi = [u1; u2; u3; u4];
ordineIngressi = size(ingressi,1);

% Vettore contenente tutte le uscite da analizzare
uscite = stati; % Suppongo lo stato accessibile

%% 2. Definizione del sistema di riferimento.

% Obiettivo della sezione: definizione delle equazioni che descrivono il
% sistema dinamico a tempo continuo.

x1_dot = x3;

x2_dot = x4;

x3_dot = ...
    - (K_el / J_m) * x1 ...
    + ((n * K_el) / J_m) * x2 ...
    - (D_el / J_m) * x3 - (tauf_segnato / J_m) * tanh(x3) ...
    + ((n * D_el) / J_m) * x4 ...
    + (1 / J_m) * u1;

x4_dot = ...
    + ((n * K_el) / J_l) * x1 ...
    - ((n^2 * K_el) / J_l) * x2 ...
    + ((n * D_el) / J_l) * x3 ...
    - ((n^2 * D_el) / J_l) * x4 ...
    - (1 / J_l) * u2;

% Creazione di una variabile contenente tutte le equazioni sopra definite
sys = [x1_dot; x2_dot; x3_dot; x4_dot]; % Sistema non lineare a tempo continuo

%% 3. Definizione dei parametri delle equazioni.

J_m = 0.001;        % [kg*m^2]
J_l = 9;            % [kg*m^2]
K_el = 0.7;         % [Nm/rad]
D_el = 0.03;        % [Nms/rad]
tauf_segnato = 0.2; % [Nm]
n = 250;            % Adimensionale

% Valutazione del sistema di riferimento nei parametri forniti dal testo
sys = subs(sys);

%% 4. Estrazione delle matrici del sistema linearizzato.

% Obiettivo della sezione: calcolo della matrice Jacobiana (matrice di
% derivate parziali prime) del sistema e delle uscite, entrambi sia 
% rispetto gli stati sia rispetto gli ingressi.

A_lin = jacobian(sys, stati);
B_lin = jacobian(sys, ingressi);
C_lin = jacobian(uscite, stati);
D_lin = jacobian(uscite, ingressi);

%% 5. Valutazione all'equilibrio.

% Obiettivi della sezione:
% _ estrazione delle matrici che descrivono la dinamica del sistema
%   linearizzato nell'intorno del punto di equilibrio;
% _ calcolo degli ingressi di equilibrio.

% Punto di equilibrio
x1 = 425.0286;
x2 = 1.7;
x3 = 0;
x4 = 0;
x_eq = [x1; x2; x3; x4];

% Valutazione delle matrici del sistema di riferimento all'equilibrio
A_eq = double(subs(A_lin));
B_eq = double(subs(B_lin));
C_eq = double(subs(C_lin));
D_eq = double(subs(D_lin));

if rank(ctrb(A_eq, B_eq)) ~= ordineStati
    fprintf("Sistema non completamente raggiungibile: rank(A_eq, B_eq) = %d != %d.\n", rank(ctrb(A_eq, B_eq)), ordineStati);
    return
end

% Calcolo degli ingressi di equilibrio
u_eq = pinv(B_eq) * ((eye(ordineStati) - A_eq) * x_eq);

%% 6. Modello del sistema nello spazio degli stati.

sys_c = ss(A_eq, B_eq, C_eq, D_eq); % Modello a tempo continuo

Ts = 0.1; % Tempo di campionamento

sys_d = c2d(sys_c, Ts); % Modello a tempo discreto (Eulero in avanti)

%% Sezione di cancellazione.

% Obiettivo della sezione: pulizia della workspace dalle variabili inutili
% per il proseguo.

clear x1 x2 x3 x4
clear u1 u2 u3 u4
clear ingressi stati uscite
clear x1_dot x2_dot x3_dot x4_dot
clear sys sys_c
clear A_lin B_lin C_lin D_lin
clear A_eq B_eq C_eq D_eq

%% 7. Setup parametri per MPC.

% Obiettivo della sezione: definizione di tutte le variabili utili per la
% risoluzione del problema di programmazione quadratica.

% Parametri dei vincoli
coppiaMotoreMin = -0.25;            % Lower bound del vincolo su u1 [Nm] 
coppiaMotoreMax = 0.25;             % Upper bound del vincolo su u1 [Nm]
tau_l = 5.0050;                     % Vincolo su u2 [Nm]
posizioneAngolareCaricoMin = -0.07; % Lower bound del vincolo su x2 [rad]
posizioneAngolareCaricoMax = 1.85;  % Upper bound del vincolo su x2 [rad]
ub = 10^6;                          % Upper bound per le variabili non vincolate (scelta arbitraria)
lb = -ub;                           % Lower bound per le variabili non vincolate (scelta arbitraria)

Tsim = 300; % Orizzonte di simulazione

% Vincoli sugli stati
Hx = [eye(ordineStati); -eye(ordineStati)];
hx_max = [ub; posizioneAngolareCaricoMax; ub; ub];
hx_min = [lb; posizioneAngolareCaricoMin; lb; lb];
hx = [hx_max - x_eq; - hx_min + x_eq];

% Vincoli sugli ingressi
Hu = [eye(ordineIngressi); -eye(ordineIngressi)];
hu_max = [coppiaMotoreMax; tau_l; ub; ub];
hu_min = [coppiaMotoreMin; tau_l; lb; lb];
hu = [hu_max - u_eq; - hu_min + u_eq];

controllerN = 3; % Orizzonte/passi di predizione dell'MPC
if controllerN == 1
    N = 5;
elseif controllerN == 2
    N = 8;
elseif controllerN == 3
    N = 10;
elseif controllerN == 4
    N = 12;
elseif controllerN == 5
    N = 14;
else
    fprintf("Errore nella scelta del parametro N.\n")
    return
end

% Valori per controllerQR:
% _ 0 -> azione di controllo imparziale;
% _ 1 -> azione di controllo aggressiva;
% _ 2 -> azione di controllo cauta;
% _ 3 -> azione di controllo molto aggressiva;
% _ 4 -> azione di controllo molto cauta;
% _ 5 -> azione di controllo estremamente aggressiva;
% _ 6 -> azione di controllo estremamente cauta.
controllerQR = 0;
if controllerQR == 0
    Q = eye(ordineStati);
    R = eye(ordineIngressi);
elseif controllerQR == 1
    Q = 100 * eye(ordineStati);
    R = eye(ordineIngressi);
elseif controllerQR == 2
    Q = eye(ordineStati);
    R = 100 * eye(ordineIngressi);
elseif controllerQR == 3
    Q = 1000 * eye(ordineStati);
    R = eye(ordineIngressi);
elseif controllerQR == 4
    Q = eye(ordineStati);
    R = 1000 * eye(ordineIngressi);
elseif controllerQR == 5
    Q = 10000 * eye(ordineStati);
    R = eye(ordineIngressi);
elseif controllerQR == 6
    Q = eye(ordineStati);
    R = 10000 * eye(ordineIngressi);
else
    fprintf("Errore nella scelta dei parametri Q e R.\n")
    return
end

% Cambio di notazione opzionale: utilizzo x_ref e u_ref anzichè x_eq e u_eq
% per identificare lo stato di equilibrio e l'ingresso di equilibrio nelle
% coordinate del sistema originario.
x_ref = x_eq; % [425.0286; 1.7; 0; 0]
u_ref = u_eq; % [0.02002; 5.005; 0; 0]

% Stato di equilibrio e ingresso di equilibrio nelle coordinate del sistema
% linearizzato
x_ref_lin = x_ref - x_eq; % [0; 0; 0; 0]
u_ref_lin = u_ref - u_eq; % [0; 0; 0; 0]

%% 8. Design MPC.

mpc = designMPC(sys_d, Hx, hx, Hu, hu, x_ref_lin, u_ref_lin, Q, R, N);

%% Sezione di cancellazione.

% Obiettivo della sezione: pulizia della workspace dalle variabili inutili
% per il proseguo.

clear sys_d
clear x_eq u_eq
clear coppiaMotoreMin coppiaMotoreMax
clear tau_l
clear posizioneAngolareCaricoMin posizioneAngolareCaricoMax
clear lb ub
clear controllerQR controllerN Q R N
clear hx_min hx_max Hx hx
clear hu_min hu_max Hu hu

%% 9. Simulazione.

% Obiettivo della sezione: simulazione del sistema reale sfruttando l'MPC.

x_log = zeros(ordineStati, Tsim + 1); % Vettore in cui salvo, durante la simulazione, il valore dello stato
u_log = zeros(ordineIngressi, Tsim); % Vettore in cui salvo, durante la simulazione, il valore dell'ingresso

% Condizione iniziale nelle coordinate del sistema originario
x_log(:,1) = zeros(ordineStati,1);

for tt = 1:Tsim
    
    x_lin = x_log(:,tt) - x_ref; % Condizione iniziale per MPC
    x_linn = x_lin - x_ref_lin;
    
    f1 = x_linn' * mpc.f1_cost;
    f_qp = f1 + mpc.f2;
    
    b_x0 = mpc.b_cost * x_linn;
    b_qp = [mpc.htilde_u; mpc.htilde_x - b_x0];
    
    warning('off', 'all');
    options = optimoptions('quadprog', 'Display', 'off');
    [du, ~, exitflag] = quadprog(mpc.H_qp, f_qp, mpc.A_qp, b_qp, [], [], [], [], [], options);
    warning('on', 'all');
    
    if exitflag ~= 1 % exitflag = 1 indica che il minimo è stato trovato correttamente
        fprintf('Problemi all''iterazione %d/%d -> esito: %d.\n', tt, Tsim, exitflag);
        % https://it.mathworks.com/help/optim/ug/quadprog.html#mw_bd42ef06-6096-4303-afaa-7b3cb9c539b6
    end

    u_log(:,tt) = u_ref + du(1:4,1);
    
    dxdt = @(t, x) MyServomechanism(t, x, u_log(:,tt), K_el, D_el, J_m, J_l, n, tauf_segnato);
    [~, xx] = ode45(dxdt, [0 Ts], x_log(:,tt));
   
    x_log(:,tt+1) = xx(end,:)';

end

%% Sezione di cancellazione.

% Obiettivo della sezione: pulizia della workspace dalle variabili inutili
% per il proseguo.

clear J_m J_l K_el D_el tauf_segnato n
clear Ts Tsim
clear x_lin x_linn
clear x_ref x_ref_lin
clear u_ref u_ref_lin
clear mpc sys_d tt
clear f1 f_qp b_x0 b_qp
clear du exitflag options
clear dxdt xx

%% 10. Visualizzazione.

descrizioneStati = [ ...
    "Posizione angolare motore"; ...
    "Posizione angolare carico"; ...
    "Velocità angolare motore"; ...
    "Velocità angolare carico" ...
];

descrizioneIngressi = [ ...
    "Coppia motore"; ...
    "Coppia carico"; ...
    "Ingresso fittizio"; ...
    "Ingresso fittizio" ...
];

% Valori per visual:
% _ 1  -> Plot stati con Q ed R imparziali;
% _ 2  -> Plot ingressi con Q ed R imparziali;
% _ 3  -> Plot stati al variare di Q, con R costante;
% _ 4  -> Plot ingressi al variare di Q, con R costante;
% _ 5  -> Plot stati al variare di R, con Q costante;
% _ 6  -> Plot ingressi al variare di R, con Q costante;
% _ 7  -> Plot stati al variare di N mantenendo Q ed R imparziali;
% _ 8  -> Plot ingressi al variare di N mantenendo Q ed R imparziali.
visual = 1;

if visual == 1
    for i = 1:ordineStati
        subplot(2,2,i)
        plot(x_log(i,:))
        title(descrizioneStati(i));
        grid on
        axis tight
        hold on
        xlabel('t');
        ylabel(sprintf("x_%d(t)",i));
    end

elseif visual == 2
    for i = 1:ordineIngressi
        subplot(2,2,i)
        plot(u_log(i,:))
        title(descrizioneIngressi(i));
        grid on
        axis tight
        xlabel('t');
        ylabel(sprintf("u_%d(t)",i));
    end

elseif visual == 3
    x1 = load('QR_caso0.mat');
    x2 = load('QR_caso1.mat');
    x3 = load('QR_caso3.mat');
    x4 = load('QR_caso5.mat');
    x = [x1; x2; x3; x4];
    for i = 1:size(x,1)
        subplot(2,2,i)
        hold on
        grid on
        axis tight
        title(descrizioneStati(i));
        for j = 1:size(x,1)
            plot(x(j).x_log(i,:))
        end
        xlabel('t');
        ylabel(sprintf("x_%d(t)",i));
        legend({'Q=I','Q=100I','Q=1000I','Q=10000I'});
    end

elseif visual == 4
    x1 = load('QR_caso0.mat');
    x2 = load('QR_caso1.mat');
    x3 = load('QR_caso3.mat');
    x4 = load('QR_caso5.mat');
    x = [x1; x2; x3; x4];
    for i = 1:size(x,1)
        subplot(2,2,i)
        hold on
        grid on
        axis tight
        title(descrizioneIngressi(i));
        for j = 1:size(x,1)
            plot(x(j).u_log(i,:))
        end
        xlabel('t');
        ylabel(sprintf("u_%d(t)",i));
        legend({'Q=I','Q=100I','Q=1000I','Q=10000I'});
    end
    
elseif visual == 5
    x1 = load('QR_caso0.mat');
    x2 = load('QR_caso2.mat');
    x3 = load('QR_caso4.mat');
    x4 = load('QR_caso6.mat');
    x = [x1; x2; x3; x4];
    for i = 1:size(x,1)
        subplot(2,2,i)
        hold on
        grid on
        axis tight
        title(descrizioneStati(i));
        for j = 1:size(x,1)
            plot(x(j).x_log(i,:))
        end
        xlabel('t');
        ylabel(sprintf("x_%d(t)",i));
        legend({'R=I','R=100I','R=1000I','R=10000I'});
    end

elseif visual == 6
    x1 = load('QR_caso0.mat');
    x2 = load('QR_caso2.mat');
    x3 = load('QR_caso4.mat');
    x4 = load('QR_caso6.mat');
    x = [x1; x2; x3; x4];
    for i = 1:size(x,1)
        subplot(2,2,i)
        hold on
        grid on
        axis tight
        title(descrizioneIngressi(i));
        for j = 1:size(x,1)
            plot(x(j).u_log(i,:))
        end
        xlabel('t');
        ylabel(sprintf("u_%d(t)",i));
        legend({'R=I','R=100I','R=1000I','R=10000I'});
    end

elseif visual == 7
    x1 = load('N_caso1.mat');
    x2 = load('N_caso2.mat');
    x3 = load('N_caso3.mat');
    x4 = load('N_caso4.mat');
    x5 = load('N_caso5.mat');
    x = [x1; x2; x3; x4; x5];
    for i = 1:ordineStati
        subplot(2,2,i)
        title(descrizioneStati(i));
        hold on
        grid on
        axis tight
        for j = 1:size(x,1)
            plot(x(j).x_log(i,:))
        end
        xlabel('t');
        ylabel(sprintf("x_%d(t)",i));
        legend({'N=5','N=8','N=10','N=12','N=14'});
    end

elseif visual == 8
    x1 = load('N_caso1.mat');
    x2 = load('N_caso2.mat');
    x3 = load('N_caso3.mat');
    x4 = load('N_caso4.mat');
    x5 = load('N_caso5.mat');
    x = [x1; x2; x3; x4; x5];
    for i = 1:ordineIngressi
        subplot(2,2,i)
        title(descrizioneIngressi(i));
        hold on
        grid on
        axis tight
        for j = 1:size(x,1)
            plot(x(j).u_log(i,:))
        end
        xlabel('t');
        ylabel(sprintf("u_%d(t)",i));
        legend({'N=5','N=8','N=10','N=12','N=14'});
    end

else
    fprintf("Errore nella scelta del parametro di visualizzazione.\n");

end

clear;
