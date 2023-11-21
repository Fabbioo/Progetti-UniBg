
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
% _ B -> Matrice relativa alla dinamica degli ingressi del sistema;
% _ N -> Orizzonte di predizione.

% Output:
% _ B_calligrafica -> Matrice contenente il movimento forzato (pag. 12 exe. 4).

function B_calligrafica = calcola_MovimentoForzato(A, B, N)
    
    % 1^ iterazione (inizializzazione)
    zeri = zeros(size(B));
    B_calligrafica = [zeros(size(B)); B];
    
    % N-1 iterazioni
    for i = 1 : N-1
        
        riga = [];
        
        for k = i:-1:1
            riga = [riga A^k*B];
        end

        B_calligrafica = [B_calligrafica; riga];
        
        zeri = [zeros(size(B)); zeri];
        B_calligrafica = [B_calligrafica [zeri; B]];

    end

    % N iterazioni totali, da 0 a N-1.

end
