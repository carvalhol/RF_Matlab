clc
close all

N=10000;

eta = linspace(0, 3, N);
zeta = linspace(1e-2, 10, N);
ell_c = 4;

%eta = r/ell_c;
%zeta = k/2/pi/ell_c;

r=eta*ell_c;
k=2*pi*ell_c*zeta;

%Auto Correlation Function (ACF)
figure (1)
hold on

%Exponential
ACF_exp = exp(-2*eta);
plot(eta, ACF_exp);
%Power-Law
ACF_pwrlaw = 1./(1+(pi^2*eta.^2)./4).^2;
plot(eta, ACF_pwrlaw);
%Gaussian
ACF_gauss = exp(-pi*eta.^2);
plot(eta, ACF_gauss);
%Triangular
coef = 2*pi*eta;
ACF_tri = 12*(2-2*cos(coef) - coef.*sin(coef))./(coef).^4;
plot(eta, ACF_tri);
%Low-pass white noise
coef = 3*pi*eta/2;
ACF_lowpass = 3*(sin(coef)-coef.*cos(coef))./(coef).^3;
plot(eta, ACF_lowpass);

hold off

%Power Spectral Density Function (PSDF)
figure (2)
hold on

%Exponential
PSDF_exp = 1/8/pi^2./(1+zeta.^2/4).^2;
plot(zeta, PSDF_exp);
%Power-Law
PSDF_pwrlaw = exp(-2.*zeta/pi)./pi^4;
plot(zeta, PSDF_pwrlaw);
%Gaussian
PSDF_gauss = 1/8/pi^3*exp(-zeta.^2./4/pi);
plot(zeta, PSDF_gauss);
%Triangular
PSDF_tri = 3/8/pi^4*(1-zeta/2/pi).*(zeta<2*pi);
plot(zeta, PSDF_tri);
%Low-pass white noise
PSDF_lowpass = 2/9/pi^4.*(zeta<3*pi/2);
plot(zeta, PSDF_lowpass);

ylim([10^(-4), 10^(-1)])
set(gca, 'YScale', 'log', 'XScale', 'log')
hold off


ell_c_out=2*trapz(r, ACF_gauss);