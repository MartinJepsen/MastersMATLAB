function fft_plot(output, input, t)
% FFT_PLOT(output, input, t)
% output: system output
% input: system input. chose input=1 to not normalise by the input


N = numel(t);                           % number of samples
bins = [0:N-1];                         % frequency bins
N2 = ceil(N/2);                         % half of N
fs = (t(2)-t(1))^-1;                    % sampling frequency
faxis=bins*fs/N;                        % frequency axis

f = figure;
f.Units = 'centimeters';
f.Position = [2, 2, 15, 8];
fs = 12;

f_output = fft(output);
f_input = fft(input);

f_plot = f_output./f_input;

plot(faxis(1:N2),abs(f_plot(1:N2))/(N2),'k')
% plot([0:1:N-1]*fs/N-fs/2,abs(f_plot)/fs,'k')
grid minor
set(gca,'XScale','linear')
set(gca,'YScale','log')
xlabel('Frequency [Hz]','Interpreter','latex','FontSize',fs)
ylabel('Amplitude','Interpreter','latex','FontSize',fs)
