
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Autore: Fabio Filippo Mandalari.                      %
%                       Matricola: 1047426.                               %
%       Corso: CAM (Controllo Avanzato Multivariabile), UniBg.            %
%           Docenti: Prof. Antonio Ferramosca, Ing. Marco Polver.         %
%   Progetto relativo al controllo di un servomeccanismo mediante MPC.    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Funzione.

% Input:
% _ R -> Matrice per il peso degli ingressi;
% _ N -> Orizzonte di predizione.

% Output:
% _ Rtilde -> Matrice che pesa gli ingressi (pag. 18 exe. 4).

function Rtilde = calcola_Rtilde(R, N)

    % 1^ iterazione (inizializzazione)
    zeri = [];
    Rtilde = R;
    
    % N-1 iterazioni
    for i = 1 : N-1
        zeri = [zeri zeros(size(R))];
        Rtilde = [Rtilde zeri'; zeri R];
    end

    % N iterazioni totali, da 0 a N.

end
