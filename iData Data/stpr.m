zy = spclow0;
zmy = modellow0;
zx = xaxislow0;
ox = xaxislowpt1;
omy = modellowpt1;
oy = spclowpt1;
fx = xaxislowpt5;
fmy = modellowpt5;
fy = spclowpt5;

plot(zx,zy,zx,zmy,'linewidth',1.5);hold on
plot(ox,oy,ox,omy,'linewidth',1.5)
plot(fx,fy,fx,fmy,'linewidth',1.5);hold off
legend({'H_{2}O','Fit','0.1 D-gly','Fit','0.5 D-gly','Fit'},'location','best')
xlabel('time (ns)')
ylabel('FID intensity')
yticks([])
xlim([-500 max(xaxis0)*1.1])
ax = gca;
ax.FontSize = 24;