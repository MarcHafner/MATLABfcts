
DoseResp_fct = @(a,b,c,x) a + (1-a) ./ ( 1 + (x./b).^c);


get_newfigure(1,[50 50 450 400])
for i=1:4
    
    x = 10.^(-3:.5:1.5);
    x2 = 10.^(-4:.1:2);
    y = DoseResp_fct(.7-.3*i, 10^(.5*i-1.5), 1.3+.2*i, x) + (.15-.3*rand(size(x)));
        
    [xI50, Hill, Einf, xMax, Area, r2, EC50, fit_final, p, log, flag] = ICcurve_fit(x,y,'GR50');
    
    subplot(2,2,i)
    semilogx(x,y,'o')
    hold on
    plot(x2, fit_final(x2), '-k')
       
    plot(xI50*[1 1 1e-4], [-1 .5 .5],'r-')
    plot([1 1e3],Einf*[1 1], '-b')
    plot([1 1e3],xMax*[1 1], '-c')
    
    
    xlim([1e-4 1e2])
    ylim([-1 1.2])
end
