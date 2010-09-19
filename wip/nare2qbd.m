function [A0 A1 A2 Y]=nare2qbd(reduction,A,B,C,D,X)
%converts a NARE to a QBD process
if(strcmp(reduction,'ramaswami'))
  n=size(A,1);
  m=size(D,1);
  t=1/max(max(diag(A)),max(diag(D)));
  A0=[eye(m)-t*D zeros(m,n); t*B zeros(n,n)];
  A1=[zeros(m,m) t*C; zeros(n,m) -t*A];
  A2=[zeros(m,m) zeros(m,n); zeros(n,m) eye(n,n)];
  if(exist('X'))
    Y=[eye(m,m)-t*D+t*C*X zeros(m,n); X zeros(n,n)];
  end
else
  error('no such reduction');
end