% VX_robot.m [DAFXbook, 2nd ed., chapter 7]
%===== this program performs a robotization of a sound
clear; clf
%----- user data -----
n1 = 441; % analysis step [samples]
n2 = n1; % synthesis step [samples]
s_win = 1024; % analysis window length [samples]
[DAFx_in,FS] = audioread('426810__pax11__psalm10.wav');
%----- initialize windows, arrays, etc -----
w1 = hann(s_win, 'periodic'); % analysis window
w2 = w1; % synthesis window
L = length(DAFx_in);
DAFx_in = [zeros(s_win, 2); DAFx_in; zeros(s_win-mod(L,n1),2)] / max(abs(DAFx_in));
DAFx_out = zeros(length(DAFx_in),1);
tic
%UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
pin = 0;
pout = 0;
pend = length(DAFx_in)-s_win;
while pin<pend
grain = DAFx_in(pin+1:pin+s_win).* w1;
%===========================================
f = fft(grain);
r = abs(f);
grain = fftshift(real(ifft(r))).*w2;
% ===========================================
DAFx_out(pout+1:pout+s_win) = ...
DAFx_out(pout+1:pout+s_win) + grain;
pin = pin + n1;
pout = pout + n2;
end
%UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
toc
%----- listening and saving the output -----
% DAFx_in = DAFx_in(s_win+1:s_win+L);
DAFx_out = DAFx_out(s_win+1:s_win+L) / max(abs(DAFx_out));
%soundsc(DAFx_out, FS);
audiowrite('redwheel_robot.wav', DAFx_out, FS);