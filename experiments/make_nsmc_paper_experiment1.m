function times=make_nsmc_paper_experiment1()
experiment='rand';
size=100;
lambdas=[4000:10:5000];
times=[];
for lambda=lambdas
    [a,M,b]=make_mbt(experiment,size,lambda);
    f=@() qve_newton(a,M,b,1.e-13,1e8);
    fp=@() perron_iteration(a,b,'eigs',1.e-13,1e8);
    fpn=@() perron_newton(a,b,'eigs',1.e-13,1e8);    
    times=[times [timeit(f);timeit(fp);timeit(fpn)]];
end
save 'nsmc_paper_last_experiment1' experiment size times lambdas;
plot(lambdas,times(1,:),'-',lambdas,times(2,:),'--',lambdas,times(3,:),'-.');