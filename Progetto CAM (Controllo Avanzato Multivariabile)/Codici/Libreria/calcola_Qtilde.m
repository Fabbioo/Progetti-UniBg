
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Autore: Fabio Filippo Mandalari.                      %
%                       Matricola: 1047426.                               %
%       Corso: CAM (Controllo Avanzato Multivariabile), UniBg.            %
%           Docenti: Prof. Antonio Ferramosca, Ing. Marco Polver.         %
%   Progetto relativo al controllo di un servomeccanismo mediante MPC.    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Funzione.

% Input:
% _ Q -> Matrice per il peso degli stati;
% _ P -> Peso dello stato finale (soluzione dell'equazione di Riccati);
% _ N -> Orizzonte di predizione.

% Output:
% _ Qtilde -> Matrice che pesa gli stati (pag. 18 exe. 4).

function Qtilde = calcola_Qtilde(Q, P, N)
    
    % 1^ + 2^ iterazione (inizializzazione)
    zeri = zeros(size(Q));
    Qtilde = [Q zeri; zeri Q];
    
    % N-2 iterazioni
    for i = 3 : N
        zeri = [zeri zeros(size(Q))];
        Qtilde = [Qtilde zeri'; zeri Q];
    end
    
    % (N+1)-esima iterazione
    zeri = [zeri zeros(size(Q))];
    Qtilde = [Qtilde zeri'; zeri P];
    
    % N+1 iterazioni totali, da 0 a N-1 per inserire sulla diagonale di
    % Qtilde le matrici Q + ultima iterazione per l'inserimento
    % nell'ultimo elemento della diagonale di Qtilde della matrice P.

end
