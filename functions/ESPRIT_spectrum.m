function [A S D doa]= ESPRIT_spectrum(X1,X2,F)
% X1 = AB^t
% X2 = ADB^t
% F: number of sources

%% IF (ESPRIT)
[m1,L]=size(X1);
[m2,L]=size(X2);

X = [X1;X2];
[Us S V]=svd(X);
U = Us(:,1:F);

U1 = U(1:m1,:); U2=U(m1+1:end,:);

E = [(U1)';(U2)']*[U1,U2];
[VE EE]=eig((E+E')/2);
Ediag=diag(EE);
[Ediag sort_ord] = sort(Ediag,'descend');
VE=VE(:,sort_ord);

V12 = VE(1:F,F+1:end);
V22 = VE(F+1:end,F+1:end);

Psi_TLS = (-V12)*inv(V22);
[invT Phi]= eig(Psi_TLS);

A = U1*(invT);
S = pinv(A)*X1; 
D = Phi;%inv(invT)*Psi_TLS*invT;
w= angle(Phi);
doa=asin(diag(w/(pi)))*180/pi;
D = Phi;%inv(invT)*Psi_TLS*invT;
