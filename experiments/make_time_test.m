function [times iters]=make_time_test(experiment,lambdas)
iters=[];
times=[];
tries1=100;
tries2=100;
tries3=100;
%tries=1;

for lambda=lambdas
    [a,M,b]=make_mbt(experiment,lambda);
    tic;
    for k=1:tries1
        [x1 res1 iter1]=qve_newton_opt(a,b,1.e-13,1e8);
    end
    time1=toc;
    
    tic;
    for k=1:tries2
        [x2 res2 iter2]=perron_iteration_eig_norm1_opt(a,b,1.e-13,1e8);
    end
    time2=toc;

    for k=1:tries3
        [x3 res3 iter3]=perron_newton_eig_norm1_opt(a,b,1.e-13,1e8);
    end
    time3=toc;

    
    
    times=[times [time1/tries1;time2/tries2;time3/tries3]];
    iters=[iters [iter1;iter2;iter3]];
end
save 'last_time_experiments_notimeit' iters times lambdas;
plot(lambdas,times(1,:),'-',lambdas,times(2,:),'--',lambdas,times(3,:),'-.');