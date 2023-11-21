
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Autore: Fabio Filippo Mandalari.                      %
%                       Matricola: 1047426.                               %
%       Corso: CAM (Controllo Avanzato Multivariabile), UniBg.            %
%           Docenti: Prof. Antonio Ferramosca, Ing. Marco Polver.         %
%   Progetto relativo al controllo di un servomeccanismo mediante MPC.    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Funzione.

% Input:
% _ equilibri -> Vettore di equilibri da impilare;
% _ N         -> Orizzonte di predizione.

% Output:
% _ vettoreEquilibri -> Vettore contenente gli equilibri impilati N volte.

function vettoreEquilibri = calcola_xuRefLin(equilibri, N)
    
    vettoreEquilibri = [];
    
    % N iterazioni
    for i = 0:N
        vettoreEquilibri = [vettoreEquilibri; equilibri];
    end

end
