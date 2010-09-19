function [A B C X]=make_qbd(type,n,parameter)
% [A B C]=make_qbd(type,n,parameter)
% returns a QBD process X=A+BX+CX^2

% NARE from bot transformed with Ram's reduction
if(strcmp(type,'botnare')) 
    % [bot, example 4] weakly transient process
    AA=-[-0.003 0.0001;0.0001 -0.0030];
    BB=[0.0019 0.001;0.0019 0.001];
    CC=[0.0015 0.0015;0.0029 0.0001];
    DD=-[-0.003 0;0 -0.003];
    XX=[19/30 1/3;19/30 1/3];

    %ramaswami's red
    t=1/max(max(diag(AA)),max(diag(DD)));
    A=[eye(2)-t*DD zeros(2); t*BB zeros(2)];
    B=[zeros(2) t*CC; zeros(2) -t*AA];
    C=[zeros(2) zeros(2); zeros(2) eye(2)];
    X=[eye(2)-t*DD+t*CC*XX zeros(2); XX zeros(2)];
elseif(strcmp(type,'cosi'))
% Queueing networks and Markov chains: modeling and performance evaluation ...
% Di Gunter Bolch, Stefan Greiner, Hermann de Meer, Kishor Shridharbhai Trivedi
    n=4;
    A=0.15*eye(n);
    B=-0.75*eye(n);B(1,1)=-0.65;B(1,2)=0.5;
    B=B+eye(n);
    C=diag(0.6*ones(n-1,1),-1);
else
  error('no such qbd'); 
end