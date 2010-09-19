1;
n=64;
alpha=1e-9;
c=1-1e-9;
[delta d q]=make_lu(n,c,alpha);
nsteps=20;
[x reshistory]=lu_newton_slow(delta,d,q,0,nsteps);
[x reshistorymod]=lu_mod_newton_slow(delta,d,q,0,nsteps);
v=[reshistory; reshistorymod];
save "../1919.dat" v
