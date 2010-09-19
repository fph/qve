function make_time_test(experiment,lambdas)
iters=[];
tries=100;
%tries=1;

for lambda=lambdas
    [a,M,b]=make_mbt(experiment,lambda);
    tic;
    for k=1:tries
        x1=qve_newton(a,M,b,1.e-13);
    end
    time1=toc;
    
    tic;
    for k=1:tries
        x2=perron_iteration_eig(a,b,1.e-13,1000);
    end
    time2=toc;

    iters=[iters [time1;time2]/tries];
end
save 'last_time_experiments_notimeit' iters;
plot(lambdas,iters(1,:),'-',lambdas,iters(2,:),'--');