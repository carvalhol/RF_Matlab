function sig_tap=cosine_taper(sig,tx)
% function sig_tap=cosine_taper(sig,tx)
%
%   Applique une pond�ration de type "cosine taper" � un 
%   signal temps (utilisation : � appliquer avant une FFT
%   afin d'�viter les effets de bords).
%
%   sig     :   signal entrant (signaux class�s en colonne)
%   tx      :   taux de la pond�ration (si 10%, entrer 0.1)
%   sig_tap :   signal trait�
%
% F. Hollender, le 04/02/2008.

[n,ns]=size(sig);
if n<ns
    sig = sig';
    N = ns;
    NS = n;
else
    N = n;
    NS = ns;
end
M=round((N-2)*tx);
w=0.5*(1-cos((0:M)*pi/(M+1)));
taper=[w ones(1,N-2*M-2) fliplr(w)]';
taper=taper*ones(1,NS);
%plot(taper)
%pause
sig_tap=sig.*taper;
if n<ns
    sig_tap = sig_tap';
end
%plot(sig_tap)
%pause

