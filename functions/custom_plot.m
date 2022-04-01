close all
t = 0:0.01:60;
u = cos(0.05*(2*pi)*t);
f = figure;
f.Units = 'centimeters';
f.Position = [2, 2, 15, 5];
fs = 12;

grid minor
[t_s, u_s] = shift(300, t, u(:));

hold on
set(gca, {'XMinorTick' ,'YMinorTick'},{'off', 'off'})
plot(t, u(:),'k')
stairs(t_s, u_s ,'r')
text(t_s(5), u_s(5)+0.1, '$\overbrace{~\quad\quad~}$','Interpreter','Latex','FontSize',fs)
text(t_s(5)+(t_s(6)-t_s(5))/2, u_s(5)+0.3, '$\Delta t$','Interpreter',...
    'Latex','FontSize',fs, 'HorizontalAlignment','center')



legend('$u(t)$','$u_k$','Interpreter','Latex','FontSize',fs)
xlabel('$t~[s]$','Interpreter','Latex','FontSize',fs,'Color','k')
xlim([0, 30])
xticks([0:5:40])
% xticklabels([])
yticklabels([])
% ylabel('$t,~k$','Interpreter','Latex','FontSize',fs)

% exportgraphics(f,'ZOH.pdf','Resolution',2)



%%
f = figure;
f.Units = 'centimeters';
f.Position = [2, 2, 15, 5];
fs = 12;
hold on
grid on

plot(t,d,'k');
plot(t,u,'r');

% numel(t(1:end-10))
% numel(s(2,1:end-9))

% plot(t,d-u,'color',[0,0.7,0])
% plot(t,s(1,:))
plot(t(500:end),x,'color',[0,0.7,0])

legend('Displacement response','Step input','Transient response','Interpreter','Latex','FontSize',10,'location','south east')
xlabel('$t~[s]$','Interpreter','Latex','FontSize',fs,'Color','k')
% exportgraphics(f,'stepinput.pdf','Resolution',100)

function [t_s, u_s] = shift(shiftcols, t, u)
c = shiftcols;
t_s = t(1:c:end-c);
u_s = u(ceil(c/2):c:end-ceil(c/2));
end