%Butterworth Discrete-Time Band-Stop IIR Filter

%Un-normalised speifications
fp1 = 42.5;
fs1 = 46.5;
fs2 = 66.5;
fp2 = 70.5;
f_samp = 260;
delta = 0.15;

%Transformed specs using Bilinear Transformation       
wp1 = tan(fp1/f_samp*pi);
ws1 = tan(fs1/f_samp*pi); 
ws2 = tan(fs2/f_samp*pi);
wp2 = tan(fp2/f_samp*pi);

%Parameters for Bandstop Transformation
W0 = sqrt(wp1*wp2);
B = wp2-wp1;
oml_s1 = (B*ws1)/((W0)^2 - (ws1)^2);
oml_s2 = (B*ws2)/((W0)^2 - (ws2)^2);
oml_s = min(abs(oml_s1),abs(oml_s2));

%Butterworth Analog LPF parameters

D1 = 1/((1-delta)^2)-1;
D2 = 1/(delta)^2 - 1;
epsilon = sqrt(D1);

N = log(sqrt(D2/D1))/log(oml_s); %Butterworth approximation
N = ceil(N);         %order
Wc = 1.08;           %cut-off frequency

%poles of Butterworth polynomial of degree 8
% following sub-part is written for N=8, can be generalised later for any N
p1 = Wc*cos(pi/2 + pi/16) + i*Wc*sin(pi/2 + pi/16);
p2 = Wc*cos(pi/2 + pi/16) - i*Wc*sin(pi/2 + pi/16);
p3 = Wc*cos(pi/2 + pi/16+pi/8) + i*Wc*sin(pi/2 + pi/16+pi/8);
p4 = Wc*cos(pi/2 + pi/16+pi/8) - i*Wc*sin(pi/2 + pi/16+pi/8);
p5 = Wc*cos(pi/2 + pi/16+2*pi/8) + i*Wc*sin(pi/2 + pi/16+2*pi/8);
p6 = Wc*cos(pi/2 + pi/16+2*pi/8) - i*Wc*sin(pi/2 + pi/16+2*pi/8);
p7 = Wc*cos(pi/2 + pi/16+3*pi/8) + i*Wc*sin(pi/2 + pi/16+3*pi/8);
p8 = Wc*cos(pi/2 + pi/16+3*pi/8) - i*Wc*sin(pi/2 + pi/16+3*pi/8);


[num,den] = zp2tf([],[p1 p2 p3 p4 p5 p6 p7 p8],Wc^N);   %TF with poles p1-p8 and numerator Wc^N and no zeroes
                                                        %numerator chosen to make the DC Gain = 1

%Evaluating Frequency Response of Final Filter
syms s z;
analog_lpf(s) = poly2sym(num,s)/poly2sym(den,s);        %analog LPF Transfer Function
analog_bsf(s) = analog_lpf((B*s)/(s*s + W0*W0));        %bandstop transformation
discrete_bsf(z) = analog_bsf((z-1)/(z+1));              %bilinear transformation

%coeffs of analog bsf
[ns, ds] = numden(analog_bsf(s));                   %numerical simplification to collect coeffs
ns = sym2poly(expand(ns));                          
ds = sym2poly(expand(ds));                          %collect coeffs into matrix form
k = ds(1);    
ds = ds/k;
ns = ns/k;

%coeffs of discrete bsf
[nz, dz] = numden(discrete_bsf(z));                     %numerical simplification to collect coeffs                    
nz = sym2poly(expand(nz));
dz = sym2poly(expand(dz));                              %coeffs to matrix form
kd = dz(1);                                             %normalisation factor
k = dz(1);    
dz = dz/k;
nz = nz/k;
fvtool(nz,dz,'Analysis','freq');                        %frequency response
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%magnitude plot (not in log scale) 
[H,f] = freqz(nz,dz,1024*1024, f_samp);
figure;
plot(f,abs(H),'LineWidth',1);
hold on;
title("Magnitude Response")
xlabel('Frequency (in kHz)');
ylabel('Magnitude Response');
xline(fp1,'--m');
xline(fs1,'--g');
yline(delta,'r');
xline(fs2,'--g');
xline(fp2,'--m');
yline(1+delta,'r');
yline(1-delta,'r');
grid
legend('Magnitude Response','Passband edge','Stopband edge','Tolerances','location','northeast')
