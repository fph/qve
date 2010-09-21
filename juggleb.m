function newb=juggleb(b,mode)
%performs various de-symmetrizations on b without
%affecting the quadratic form t->b(t,t)
%
% newb=juggleb(oldb,mode)
%
% mode may be: 'nothing', 'transpose','symmetrize','desymmetrize1','desymmetrize2'
%
% see [Meini, Poloni '10, Arxiv]
% (c) f.poloni@sns.it 2009-2010

n=size(b,1);
if(strcmp(mode,'nothing'))
  newb=b;
elseif(strcmp(mode,'transpose'))
  B=reshape(b,n,n,n);
  B=permute(B,[1 3 2]);
  newb=reshape(B,n,n*n);
elseif(strcmp(mode,'symmetrize'))
  newb=0.5*(b+juggleb(b,'transpose'));
elseif(strcmp(mode,'desymmetrize1'))
  B=reshape(b,n,n,n);
  for j=1:n
    for k=1:j-1
      B(:,j,k)=B(:,j,k)+B(:,k,j);
      B(:,k,j)=0;
    end
  end
  newb=reshape(B,n,n*n);
elseif(strcmp(mode,'desymmetrize2'))
  B=reshape(b,n,n,n);
  for j=1:n
    for k=1:j-1
      B(:,k,j)=B(:,j,k)+B(:,k,j);
      B(:,j,k)=0;
    end
  end
  newb=reshape(B,n,n*n); 
else
  error('Unknown mode');
end
