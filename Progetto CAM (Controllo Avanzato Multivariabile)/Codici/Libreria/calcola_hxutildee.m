
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Autore: Fabio Filippo Mandalari.                      %
%                       Matricola: 1047426.                               %
%       Corso: CAM (Controllo Avanzato Multivariabile), UniBg.            %
%           Docenti: Prof. Antonio Ferramosca, Ing. Marco Polver.         %
%   Progetto relativo al controllo di un servomeccanismo mediante MPC.    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Funzione.

% Input:
% _ h -> Vettore dei termini noti dei vincoli:
%   _ se hx -> Vettore dei termini noti dei vincoli relativi agli stati;
%   _ se hu -> Vettore dei termini noti dei vincoli relativi agli ingressi;
% _ N -> Orizzonte di predizione.

% Output:
% _ h_tilde -> Vettore "incrementato" dei termini noti dei vincoli:
%   _ se input = hx -> output = htilde_x -> Vettore "incrementato" dei
%     termini noti dei vincoli relativi agli stati (pag. 14 exe. 4);
%   _ se input = hu -> output = htilde_u -> Vettore "incrementato" dei
%     termini noti dei vincoli relativi agli ingressi (pag. 13 exe. 4).

function h_tilde = calcola_hxutildee(h, N)

    h_tilde = [];
    
    % N+1 iterazioni
    for i = 0:N
        h_tilde = [h_tilde; h];
    end

end
