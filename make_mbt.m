function [a M BM]=make_mbt(kind,param);
if(strcmp(kind,'latouche')==1)
  delta=param;
  D0=zeros(9);
  v=[6 6 6 1 0 6 6 6];
  D0=D0+diag(v,1);
  D0(4,1)=6;
  D0(4,6)=1;
  D0(9,1)=1;
  D0(9,5)=1;
  D0(9,6)=6;
  D0=1.e-3*D0;
  v=[delta delta delta delta 5 4 4 4 4];
  D1=1.e-2*diag(v,0);
  d=zeros(9,1);
  d(5)=1;
  P1=zeros(9);
  P1(1:4,1:4)=eye(4);
  P1(5,5)=0.9;
  P1(5,1)=0.1;
  P1(6:9,5)=1;
  P2=zeros(9);
  P2(6:9,6:9)=eye(4);
  P2(1:4,5)=1;
  P2(5,1)=0.1;
  P2(5,5)=0.9;
  v=d+sum(D0,2)+sum(D1,2);
  v=-v;
  D0=D0+diag(v,0);
  DI=-inv(D0);
  a=DI*d;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % BM e' in forma di matrix nx(n*n)
  BM=zeros(9,9*9);
  for i=1:9
    for j=1:9
      for k=1:9
        BM(i,k+9*(j-1))=D1(i,i)*P1(i,j)*P2(i,k);
      end
    end
  end
  BM=DI*BM;
  % il primo indice di B e' quello in alto tra parentesi 
  % nelle nostre notazioni
%  B=reshape(BM,9,9,9);
  M=eye(9);
elseif(strcmp(kind,'latouche2')==1)
  p=param;
  D0=-10*eye(3);D0(3,2)=1;
  d=[1 1 9]';
  R=zeros(3,9);
  R(2,1)=9*p;
  R(1,5)=9*(1-p);
  R(1,7)=4.5*p;
  R(1,9)=4.5*p;
  R(2,9)=9*(1-p);
  a=-D0\d;
  BM=-D0\R;
  M=eye(3);
else
  error('Unknown kind');
end
