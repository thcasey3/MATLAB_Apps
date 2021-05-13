clear,clc

[B,spc,params] = eprload('PREPR_5mM_TEMPO_pt5DGLY.DTA');
spcR = real(spc);
spcI = imag(spc);

m_phase = (180/pi)*acos(dot(spcR,spcI)/(sqrt(sum(spcR.^2))*sqrt(sum(spcI.^2))))

corr = pi/2-m_phase*pi/90;

sf = sin(corr);

cf = 1/cos(corr);

spcC = (spcR-spcI*sf)*cf;

phases = linspace(-pi/2, pi/2, 100);
spc2 = spc.*exp(-1i*phases);
for k = 1:100
    sum_sqrs(k) = sum(real(spc2(:,k)).^2) / sum(imag(spc2(:,k)).^2);
end
[~,b] = max(sum_sqrs);

plot(B,real(spc2(:,b)),B,spcC);legend('method1','method2')


