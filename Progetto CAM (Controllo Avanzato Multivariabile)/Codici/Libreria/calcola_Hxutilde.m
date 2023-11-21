
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Autore: Fabio Filippo Mandalari.                      %
%                       Matricola: 1047426.                               %
%       Corso: CAM (Controllo Avanzato Multivariabile), UniBg.            %
%           Docenti: Prof. Antonio Ferramosca, Ing. Marco Polver.         %
%   Progetto relativo al controllo di un servomeccanismo mediante MPC.    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Funzione.

% Input:
% _ H -> Matrice dei vincoli:
%   _ se Hx -> Matrice dei vincoli relativi agli stati;
%   _ se Hu -> Matrice dei vincoli relativi agli ingressi;
% _ N -> Orizzonte di predizione.

% Output:
% _ H_tilde -> Matrice "incrementata" dei vincoli:
%   _ se input = Hx -> output = Htilde_x -> Matrice "incrementata" dei
%     vincoli relativi agli stati (pag. 14 exe. 4);
%   _ se input = Hu -> output = Htilde_u -> Matrice "incrementata" dei
%     vincoli relativi agli ingressi (pag. 13 exe. 4).

function H_tilde = calcola_Hxutilde(H, N)
    
    % 1^ e 2^ iterazione (inizializzazione)
    zeri = zeros(size(H));
    H_tilde = [H zeri; zeri H];
    
    % N-2 iterazioni
    for i = 3:N
        zeri = [zeri zeros(size(zeri)); zeros(size(zeri)) zeri];
        H_tilde = [H_tilde zeri; zeri H_tilde];
    end
    
    % Selezione delle righe e delle colonne utili (equivalente di N+1 iterazioni)
    H_tilde = H_tilde(1:(size(H,1)*(N+1)),1:(size(H,2)*(N+1)));

end
