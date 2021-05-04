% Yuval and Bar doin' zero padding
% signals and systems final project
%      Bars's ID : 205684426
%      Yuval's ID : 312238207


%>>>>>>>>>>>constants<<<<<<<<<<<
N1 = 100;
N2 = 10000;


%>>>>>>>>>>>Q1<<<<<<<<<<<
n = 1:1:N1;

% prooving BIBO on the given system
x = zeros(N1,1);
y = zeros(N1,1);
for i = 1:N1
    x(i) = 1;
end

y(1) = 20*x(1);
y(2) = 90*x(1) + 10*x(2);
for i = 3:N1
    y(i) = 4*y(i-1)-4*y(i-2)+20*x(1)+10*x(i-1);
end

figure (101);
subplot(2,1,1);
plot (n,x(n));
grid on
title('inputs x[n] = 1');
subplot(2,1,2);
plot(n,y(n));
hold on 
title('output y(x[n]) of a bounded input');
grid on


%>>>>>>>>>>>Q2<<<<<<<<<<<
K = -N2:1:N2;
l = 1:81:1;
X = cos(0.4*pi*K)+cos(0.2*pi*K);
W_m = 0.4*pi;
W_ny= 2*W_m;
Fs = 0.5;
freq = linspace (0,2*pi, length(X));

% Q2_b
figure(201);
title('filters amplitude as a function of frequency');
plot_Q2amp(h2, 1); 
title('H_2');
plot_Q2amp(h3, 2); 
title('H_3');
plot_Q2amp(h4, 3); 
title('H_4');
plot_Q2amp(h6, 4); 
title('H_6');

figure(202);
title('filters phase as a function of frequency');
hold on
plot_Q2phase(h2, 1); 
title('H_2');
plot_Q2phase(h3, 2); 
title('H_3');
plot_Q2phase(h4, 3); 
title('H_4');
plot_Q2phase(h6, 4); 
title('H_6');

figure(203);
plot_X_jw(X);
grid on
title('X(e^jw)');
xlabel('frequency');
ylabel('Amplitude');
set (gca, 'XTick', 0:pi/4:2*pi);
set (gca, 'XTickLabel', {'0','0.25pi', '0.5pi', '0.75pi', 'pi', '1.25pi', '1.5pi', '1.75pi', '2pi'});

%Q2_d and e
figure(204);
plot_Q2d(X,h2,'H_2');

figure(205);
plot_Q2d(X,h3,'H_3');

figure(206);
plot_Q2d(X,h4,'H_4');

figure(207);
plot_Q2d(X,h6,'H_6');


%>>>>>>>>>>>Q3<<<<<<<<<<<
%that pice of code is just for visualization of our signals
k=0:1:200;
x1=sinc(k/6);
x2=(sinc(k/12)).^2;
x3=cos(k*pi/12);
x4=cos(k*pi/12)+sin(k*pi/6);

figure(301)
title('Q_3 signals');
plot_Q3_signals(k,x1,1,"x_1(t)");
plot_Q3_signals(k,x2,2,"x_2(t)");
plot_Q3_signals(k,x3,3,"x_3(t)");
plot_Q3_signals(k,x4,4,"x_4(t)");

%this glorius loop contains entire Q3
i=1;
T=[1,4,8];  %the vector of given T's

while i<4
    K=10000;
    ts=T(i);
    t=-K:ts:K;
    
    x1=sinc(t/6);
    x2=(sinc(t/12)).^2;
    x3=cos(t*pi/12);
    x4=cos(t*pi/12)+sin(t*pi/6);

    X1=fft(x1);
    X2=fft(x2);
    X3=fft(x3);
    X4=fft(x4);
    
    %so we can have some nice titles and good looking plots
    figure(301+i)
    str = 'T=%d';
    A=T(i);
    tit = sprintf(str,A);
    
    plot_Q3(X1,1,'X_1(jw),' ,tit);
    plot_Q3(X2,2,'X_2(jw), ',tit);
    plot_Q3(X3,3,'X_3(jw), ',tit);
    plot_Q3(X4,4,'X_4(jw), ',tit);

    i = i+1;
end



%>>>>>>>>>>>functions<<<<<<<<<<<

function plot_Y_jw(X, h)
%plots the DFT of the output of given x(t) signal
    Y= conv(h, X);                 %solving to find y(x)
    L = length(Y);                 %length of filter vector
    F = linspace(0,2*pi, L);       %the frequencies on the range
    FT = abs (fft(Y));
    plot (F, FT); 
    xlabel('frequency');
    ylabel('Amplitude');
    set (gca, 'XTick', 0:pi/4:2*pi);
    set (gca, 'XTickLabel', {'0','0.25pi', '0.5pi', '0.75pi', 'pi', '1.25pi', '1.5pi', '1.75pi', '2pi'});
end

function plot_X_jw(X)
%plots the DFT of the given x(t) signal
    L = length(X);                %length of filter vector
    F = linspace(0,2*pi, L);      %the frequencies on the range
    plot (F, abs (fft(X))); 
    xlabel('frequency');
    ylabel('Amplitude');
    set (gca, 'XTick', 0:pi/4:2*pi);
    set (gca, 'XTickLabel', {'0','0.25pi', '0.5pi', '0.75pi', 'pi', '1.25pi', '1.5pi', '1.75pi', '2pi'});
end

function Bode_plot_amplification(h)
%bode plot of the given filter amplification
    L = length(h);                   %length of filter vector
    F = linspace(0,pi, L/2);         %the frequencies on the range
    fourier = abs(fft(h));
    plot (F, fourier(1:L/2)); 
    xlabel('X');
    ylabel('Amplitude');
    set (gca, 'XTick', 0:pi/4:pi);
    set (gca, 'XTickLabel', {'0','0.25pi', '0.5pi', '0.75pi', 'pi'});
end

function Bode_plot_phase(h)
%bode plot of the given filter phase shift
    L = length(h);               %length of filter vector
    F = linspace(0,2*pi, L);     %the frequencies on the range
    fourier = abs(fft(h));
    plot (F, fourier); 
    xlabel('X');
    ylabel('Phase');
    set (gca, 'XTick', 0:pi/4:2*pi);
    set (gca, 'XTickLabel', {'0','0.25pi', '0.5pi', '0.75pi', 'pi', '1.25pi', '1.5pi', '1.75pi', '2pi'});
end

function plot_Q2amp(h, i) 
%plots Q2 bode plot of amplitude
    subplot(2,2,i);
    Bode_plot_amplification(h);
    grid on 
    title(i);
end

function plot_Q2phase(h, i) 
%plots Q2 bode plot of phase
    subplot(2,2,i);
    Bode_plot_phase(h);
    grid on 
    title(i);
end

function plot_Q2d(X,h,i)
%plots Q2_d figures
    title(i);
    subplot(2,1,1);
    plot_y_t(X,h);
    hold on
    xlabel('n');
    grid on
    plot_x_t(X);
    legend('y','x');
    grid on
    subplot(2,1,2);
    plot_Y_jw(X,h);
    title('Y(jw)');
    hold on
    grid on
end

function plot_Q3_signals(t,x,i, tit)
%plots Q3 signals
    subplot(2,2,i);
    plot(t,x),title(tit);
    xlabel('t');
    ylabel('x');
    grid on
end

function plot_Q3(X,i,tit1, tit2)
%plots Q3 figures
    subplot(2,2,i);
    plot_ft(X);
    s = strcat(tit1,tit2);
    title(s);
    xlabel('w');
    ylabel('amplitude');
    grid on
end

function plot_y_t(x,h)
%Plots y(n) - the output in time domain
    y = conv(x,h); % Result of convolution
    n=100:200;
    plot(n,y(100:200));    
end

function plot_x_t(x)
% Plots x(n) - the input in time domain
    n=100:200;
    plot(n,x(100:200));  
end

function plot_ft(x)
% Plots ft{x(n)} - the input in frequency domain
    w=linspace(0,2*pi,length(x));
    plot(w,abs(x));
    set(gca,'XTick',0:pi/2:2*pi); 
    set(gca,'XTickLabel',{'0','pi/2','pi','3*pi/2','2*pi'});
end