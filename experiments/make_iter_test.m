function make_iter_test(experiment,lambdas)
% makes a series of tests counting iteration numbers
%
iters=[];
for lambda=lambdas
    [a,M,b]=make_mbt(experiment,lambda);
    [x1 res1 iter1]=qve_newton(a,M,b,1.e-13);
    [x2 res2 iter2]=perron_iteration_newnorm_eig(a,b,1.e-13,1000000);
    iters=[iters [iter1;iter2;]];
end

plot(lambdas,iters(1,:),'-',lambdas,iters(2,:),'--');
