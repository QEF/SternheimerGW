
Np = 3;

par0=[0.25 1  5 0.25 1 10 0.25 1 15];

[par] = fmins ( 'multipole', par0);

for i = 1:Np
  S(i)   = par(3*i-2);
  eta(i) = par(3*i-1);
  wp(i)  = par(3*i  );
end;

% compare input data and fitting function
% along imaginary axis

tmp=load('fort.50');
beta=tmp(:,2);              
inveps=tmp(:,3)+j*tmp(:,4); 
inveps_tilde = 1;
for i = 1:Np
  inveps_tilde = inveps_tilde - abs(S(i))^2./(1 + ((beta+eta(i))/wp(i)).^2 );
end;
figure(1);
plot(beta,inveps,'or-');
hold on;
plot(beta,inveps_tilde,'ob-');

% compare ab-initio inveps on real axis
% and analytically continued function

figure(2);
hold off;
w=linspace(0,30);
inveps_re = 1;
for i = 1:Np
  inveps_re = inveps_re - abs(S(i))^2 * wp(i)/2 * (1./(w+wp(i)+j*eta(i))-1./(w-wp(i)+j*eta(i)));
end;
plot(w,imag(inveps_re),'ob-');
hold on;
plot(w,real(inveps_re),'ob-');
hold on;
tmp=load('inveps_q=0_G=0.dat');
w=tmp(:,2);          
inveps=tmp(:,3)+j*tmp(:,4); 
plot(w,imag(inveps),'or');
plot(w,real(inveps),'or');
