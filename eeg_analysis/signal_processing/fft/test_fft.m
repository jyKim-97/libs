%% test signal
srate = 1024;
tmax = 20;
t = 0:1/srate:tmax;
y1 = 2*sin(2*pi*5*t+5);
y2 = 0.5*sin(2*pi*40*t+10);
y = y1 + y2;
y(t > tmax/2) = 0;

%% FFT 
[xffts, t_new, f_new] = getSTFFT(y, t, srate, "maxf", 100, "wbin", 1024, "mbin", 512, "nf", 500);

%% show
figure("Units", "normalized", "pos", [0.1, 0.1, 0.8, 0.6]);
subplot(211)
plot(t, y, "k", "LineWidth", 0.1);
xlabel("time (s)")
ylabel("V")

ax = subplot(212);
imagesc(t, f_new, imresize(xffts', 10));
colormap jet
ax.YDir = "normal";
xlim([0, 20])
ylim([0, 60])
colorbar("southoutside")