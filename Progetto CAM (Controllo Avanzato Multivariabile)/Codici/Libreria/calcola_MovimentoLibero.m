
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Autore: Fabio Filippo Mandalari.                      %
%                       Matricola: 1047426.                               %
%       Corso: CAM (Controllo Avanzato Multivariabile), UniBg.            %
%           Docenti: Prof. Antonio Ferramosca, Ing. Marco Polver.         %
%   Progetto relativo al controllo di un servomeccanismo mediante MPC.    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Funzione.

% Input:
% _ A -> Matrice relativa alla dinamica degli stati del sistema;
% _ N -> Orizzonte di predizione.

% Output:
% _ A_calligrafica -> Matrice contenente il movimento libero dello stato (pag. 12 exe. 4).

function A_calligrafica = calcola_MovimentoLibero(A, N)

    % 1^ iterazione (inizializzazione)
    A_calligrafica = eye(size(A));
    
    % N iterazioni
    for i = 1 : N
        A_calligrafica = [A_calligrafica; A^i];
    end

    % N+1 iterazioni totali, da 0 a N.

end
