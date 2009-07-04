  %
  % Fit the imaginary axis inveps(w) using a multi-plasmon-pole expansion
  %
  function sqre = multipole ( par )

  Np = 3; % Number of poles

  tmp=load('fort.50');
  beta=tmp(:,2);              % imaginary frequency, eV
  inveps=tmp(:,3)+j*tmp(:,4); % inv(eps) on imag axis

  for i = 1:Np
    S(i)   = par(3*i-2); % pole strenght
    eta(i) = par(3*i-1); % pole width
    w(i)   = par(3*i  ); % pole frequency
  end;

  inveps_tilde = 1;
  for i = 1:Np
    inveps_tilde = inveps_tilde - abs(S(i))^2./(1 + ((beta+eta(i))/w(i)).^2 );
  end;

% square error for least squares fit

  sqre=sum(abs(inveps-inveps_tilde).^2);
