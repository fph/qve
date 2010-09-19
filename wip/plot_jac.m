1;
jrs=[];
ks=linspace(0.8,2,50);
for k=ks
    [a,M,b]=make_mbt('latouche',k);
    xn=qve_newton(a,M,b,1.e-14);
    [Jac dyad p1 dyad2 p2] = perron_jacobian(a,b,'eig',e-xn);
    jr=[svd(dyad)];
    jrs=[jrs jr];
end
plot(ks,jrs)