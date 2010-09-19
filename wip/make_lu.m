function [delta d q]=make_lu(n,c,alpha)
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
