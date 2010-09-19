function G=overrelaxation(A,B,C,mchoice,iters)
n=size(B,1);
AA=(eye(n)-B)\A;
CC=(eye(n)-B)\C;
G=0;
for iter=1:iters
if(strcmp(mchoice,'GC'))
M=G*CC;
elseif(strcmp(mchoice,'halfGC'))
M=0.5*G*CC;
elseif(strcmp(mchoice,'heuristic'))
if(iter>1)
tot=min((eye(n)-CC*G)*ones(n,1));
if(tot<0)
'heuristic failed'
M=0;
%error('heuristic failed')
endif
tot
M=tot*ones(n)/n;
else
M=0;
endif
elseif(strcmp(mchoice,'0'))
M=0;
else
error('Invalid mchoice');
endif
Gnew=(eye(n)-CC*G-M)\(AA-M*G);
G=Gnew;
res=norm(G-A-B*G-C*G*G,'fro');
iter,res
endfor
