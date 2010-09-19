function [A,B,C,D,X] = make_nare (num,parin)

% Mriccatix generates examples for the equation
%               XCX-AX-XD+B=0
% whose coefficient form an M-matrix.
% Examples arise in fluid queues.
%
% Input: num is the number of the example and
%        parin are the relative parameters
%
% Example 1. Guo
%           parin: n (dimension)[10], alfa (real number)[0]
%
% Output: the coefficients A,B,C,D and the solution X (if provided)

% Author: Bruno Iannazzo

if num==1,
    % [bot, example 1]  null recurrent process
    A=-[-0.003 0.001;0.001 -0.003];
    B=[0.001 0.001;0.001 0.001];
    C=[0.001 0.001;0.001 0.001];
    D=-[-0.003 0.001;0.001 -0.003];
    X=[0.5 0.5;0.5 0.5];

elseif num==2,
    % [bot, example 2] null recurrent process
    A=-[-100.002 100;100 -100.002];
    B=[0.001 0.001;0.001 0.001];
    C=[0.001 0.001;0.001 0.001];
    D=-[-0.003 0.001;0.001 -0.003];
    X=[0.5 0.5;0.5 0.5];

elseif num==3,
    % [bot, example 3] strongly positive recurrent process
    m=2;n=18;
    A=0.018*eye(m);
    B=0.001*ones(m,n);
    C=0.001*ones(n,m);
    D=-10*ones(n)+180.002*eye(n);
    X=1/18*ones(2,18);

elseif num==4,
    % [bot, example 4] weakly transient process
    A=-[-0.003 0.0001;0.0001 -0.0030];
    B=[0.0019 0.001;0.0019 0.001];
    C=[0.0015 0.0015;0.0029 0.0001];
    D=-[-0.003 0;0 -0.003];
    X=[19/30 1/3;19/30 1/3];

elseif num==5,
    % [guo] random choose of a singular M-matrix with Me=e
    if nargin<2 n=10;alfa=0;
    else n=parin(1); alfa=parin(2);
    end;
    R=rand(2*n);e=ones([2*n 1]);
    % the matrix must have non zero elements
    for k=1:n
        for h=1:n
            if R(h,k)==0 error('Chosen a bad random matrix, please retry.'); end
        end
    end,
    W=diag(R*e)-R;I=eye(n);O=zeros(n);
    W=alfa*[I O;O I]+W;
    D=W(1:n,1:n);C=-W(1:n,n+1:2*n);B=-W(n+1:2*n,1:n);A=W(n+1:2*n,n+1:2*n);

elseif num==11
    % Bai, Guo, Xu. Esempio 1
    m=parin;n=m^2;
    ev=ones(n,1);
    S=1/50*ev*ev';
    v=zeros(m,1);v(2)=-1;T=toeplitz(v);
    for k=1:m
        T(k,k)=4+200/(k+1)^2;
    end
    v=zeros(n,1);v(1)=2;v(2)=1;C=1/50*toeplitz(v);
    v=zeros(n,1);v(m+1)=-1;A=toeplitz(v);
    for k=0:m-1
        A(k*m+1:k*m+m,k*m+1:k*m+m)=T;
    end
    D=A;B=S*C*S-A*S-S*D;

elseif num==13
    % Bai, Guo, Xu. Esempio 3
    n=parin(1);csi=parin(2);
    v=zeros(n,1);v(1)=3;v(2)=-1;
    A=gallery('circul',v);D=A;
    B=eye(n);C=csi*eye(n);

elseif num==14
    % Bai, Guo, Xu. Esempio 4
    n=parin(1);csi=parin(2);
    ev=ones(2*n,1);

    R=rand(2*n);W=diag(R*ev)-R;k=10;I=eye(n);
    D=W(1:n,1:n)+k*I;A=W(n+1:2*n,n+1:2*n)+k*I;
    B=-W(n+1:2*n,1:n);C=-csi*W(1:n,n+1:2*n);
elseif num==15
    % [jl]
    n=parin(1);c=parin(2);alpha=parin(3);
    if(mod(n,4)~=0)
        error('n must be a multiple of 4');
    end
    omega=zeros(n,1);
    ci=zeros(n,1);
    tmp_omega=zeros(n,1);
    tmp_ci=zeros(n,1);
    gauss_x=[-0.861136312 -0.339981044 0.339981044 0.861136312];
    gauss_w=[0.347854845 0.652145155 0.652145155 0.347854845];
    for i=0:4:n-1 %points on 0..1
        a=i/n;
        b=(i+4)/n;
        tmp_ci(i+1:i+4)=(b-a)/2*gauss_w;
        tmp_omega(i+1:i+4)=(b-a)/2*gauss_x+(a+b)/2;
    end

    omega(1:n)=tmp_omega(n:-1:1); %reverses
    ci(1:n)=tmp_ci(n:-1:1);

    delta=1./(c*omega*(1+alpha));
    d=1./(c*omega*(1-alpha));
    q=0.5*ci./omega;
    A=diag(delta)-ones(n,1)*q';
    B=ones(n,n);
    C=q*q';
    D=diag(d)-q*ones(1,n);
end

% References
% [bot] Nigel G. Bean, Malgorzata M. O'Reilly and P. G. Taylor,
% "Algorithms for Return Probabilities for Stochastic Fluid Flows", preprint
%
% [guo] C.-H. Guo,
% "Nonsymmetric algebraic Riccati equations and Wiener-Hopf  factorization
% for M-matrices" SIAM J. Matrix Anal. Appl. 23-1 (2001), pp. 225-242
%
% [jl]
%Juang, Jonq; Lin, Wen-Wei
% Nonsymmetric algebraic Riccati equations and Hamiltonian-like matrices.
% SIAM J. Matrix Anal. Appl.  20  (1999),  no. 1, 228--243 (electronic).
%