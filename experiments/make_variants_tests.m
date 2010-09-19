function [iters]=make_variants_tests(experiment,lambdas,variants)
k=length(variants);
iters=[];

for lambda=lambdas
    [a,M,b]=make_mbt(experiment,lambda);
    currentiters=[];
    for i=1:k
        [x res iter]=perron_iteration_eig_norm1_opt(a,juggleb(b,variants(i)),1.e-13,1e8);
        currentiters=[currentiters; iter];
    end
    iters=[iters currentiters];
end
save 'variants_experiments' iters lambdas variants;
