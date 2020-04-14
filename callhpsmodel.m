clear all
close all
n1 = 441; % analysis step [samples]
n2 = n1; % synthesis step [samples]
s_win = 1023; % analysis window length [samples]
[x,fs] = audioread('426810__pax11__psalm10.wav');
%----- initialize windows, arrays, etc -----
w1 = hann(s_win, 'periodic'); % analysis window
w2 = w1; % synthesis window
L = length(x);
DAFx_in = [zeros(s_win, 2); x; zeros(s_win-mod(L,n1),2)] / max(abs(x));
DAFx_out = zeros(length(x),1);
soundlength= length(x);
%-----time mapping-----%
timemapping = [ 0 1; % input time (sec)
0 1 ]; % output time (sec)
timemapping = timemapping*soundlength/fs;
tic
%UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
pin = 0;
pout = 0;
pend = length(x)-s_win;
w= w1;
N= 512;
t= -2000;
nH= 4;
minf0= 80;
maxf0= 260;
f0et= 5;
maxhd= 2;
stocf= 1;
[y, yh, ys]= hpsmodelparams(x, fs, w, N, t, nH, minf0, maxf0, f0et, maxhd, stocf, timemapping);
