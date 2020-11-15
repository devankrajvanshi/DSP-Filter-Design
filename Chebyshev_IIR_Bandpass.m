%Chebyschev Discrete-Time Bandpass IIR Filter

%Un-normalised speifications
fs1 = 45.1;
fp1 = 49.1;
fp2 = 69.1;
fs2 = 73.1;
f_samp = 330;
delta = 0.15;

%Transformed specs using Bilinear Transformation
ws1 = tan(fs1/f_samp*pi);          
wp1 = tan(fp1/f_samp*pi);
wp2 = tan(fp2/f_samp*pi);
ws2 = tan(fs2/f_samp*pi);

%Parameters for Bandpass Transformation
W0 = sqrt(wp1*wp2);
B = wp2-wp1;
oml_s1 = ((ws1)^2 - (W0)^2)/(B*ws1);
oml_s2 = ((ws2)^2 - (W0)^2)/(B*ws2);
oml_s = min(abs(oml_s1),abs(oml_s2));

%Chebyshev LPF parameters
D1 = 1/((1-delta)^2)-1;
D2 = 1/(delta)^2 - 1;
epsilon = sqrt(D1);

N = acosh(sqrt(D2/D1))/acosh(oml_s); %Chebyschev approximation
N = ceil(N);

% Poles of the Chebyshev Polynomial of order N
% following sub-part is written for N=4, can be generalised later for any N
p1 = -sin(pi/(2*N))*sinh(asinh(1/epsilon)/N)+i*cos(pi/(2*N))*cosh(asinh(1/epsilon)/N);
p2 = -sin(pi/(2*N))*sinh(asinh(1/epsilon)/N)-i*cos(pi/(2*N))*cosh(asinh(1/epsilon)/N);
p3 = -sin(3*pi/(2*N))*sinh(asinh(1/epsilon)/N)+i*cos(3*pi/(2*N))*cosh(asinh(1/epsilon)/N);
p4 = -sin(3*pi/(2*N))*sinh(asinh(1/epsilon)/N)-i*cos(3*pi/(2*N))*cosh(asinh(1/epsilon)/N);        

%evaluating the Transfer function of Chebyshev Analog LPF
n1 = [1 -p1-p2 p1*p2];
n2 = [1 -p3-p4 p3*p4];
den = conv(n1,n2);          %multiply n1 and n2, which are the two quadratic factors in the denominator
num = [den(5)*sqrt(1/(1+epsilon*epsilon))];        % even order, DC Gain set as 1/(1+ epsilon^2)^0.5

%Evaluating Frequency Response of Final Filter
syms s z;
analog_lpf(s) = poly2sym(num,s)/poly2sym(den,s);    %analog lpf transfer function
analog_bpf(s) = analog_lpf((s*s +W0*W0)/(B*s));     %bandpass transformation
discrete_bpf(z) = analog_bpf((z-1)/(z+1));          %bilinear transformation

%coeffs of analog BPF
[ns, ds] = numden(analog_bpf(s));        %numerical simplification
ns = sym2poly(expand(ns));                          
ds = sym2poly(expand(ds));               %collect coeffs into matrix form
k = ds(1);    
ds = ds/k;
ns = ns/k;

%coeffs of discrete BPF
[nz, dz] = numden(discrete_bpf(z));      %numerical simplification
nz = sym2poly(expand(nz));                          
dz = sym2poly(expand(dz));               %collect coeffs into matrix form
kd = dz(1); %normalisation factor
kn = nz(1); %normalisation factor
dz = dz/kd;
nz = nz/kn;
k=kn/kd;
nz=nz*k;
fvtool(nz,dz,'Analysis','freq');         %frequency response in dB

[~,p,~]=tf2zp(nz,dz);
figure;
plot(real(p),imag(p),'rX');
title("Poles of the filter transfer function")
xlabel("Re(z)")
ylabel("Im(z)")
axis equal
grid on
t = linspace(0,2*pi,1000);
hold on
plot(cos(t),sin(t),'b-') 

%magnitude plot (not in log scale) 
[H,f] = freqz(nz,dz,1024*1024, f_samp);
figure;
plot(f,abs(H),'LineWidth',1);
hold on;
title("Magnitude Response")

%ylim([0,1.1]);
xlabel('Frequency (in kHz)');
ylabel('Magnitude Response');
xline(fp1,'--m');
xline(fs1,'--g');
yline(delta,'r');
xline(fp2,'--m');
xline(fs2,'--g');
yline(1+delta,'r');
yline(1-delta,'r');
grid
legend('Magnitude Response','Passband edge','Stopband edge','Tolerances','location','northeast')