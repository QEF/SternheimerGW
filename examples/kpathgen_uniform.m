
  dk = 0.02;

  highsymlist = [ 0.500 0.500 0.500;  % L
                  0.000 0.000 0.000;  % G
                  0.000 0.000 1.000;  % X
                  1.000 1.000 1.000]; % G


  np = max(max(size(highsymlist))); 

  klist = [];

  for ip = 1:np-1

    k1 = highsymlist(ip  ,:);
    k2 = highsymlist(ip+1,:);

    nk = ceil(norm(k1-k2)/dk);

    c1 = linspace(k1(1),k2(1),nk);
    c2 = linspace(k1(2),k2(2),nk);
    c3 = linspace(k1(3),k2(3),nk);

    klist = [klist; c1' c2' c3'];

  end;

  nnk = max(max(size(klist)));

  w = ones(1,nnk)/nnk;

  klist = [klist';w]';

  fid = fopen('klist.dat','w');
  fprintf(fid,'%10.5f %10.5f %10.5f %10.5f\n',klist');
  fclose(fid);

  s = 0;
  ss = [ s ];

  for ik = 2:nnk;
    s = s + norm(klist(ik,1:3)-klist(ik-1,1:3));
    ss = [ss; s];
  end;

  fid = fopen('kcoord.dat','w');
  fprintf(fid,'%12.7f\n',ss');
  fclose(fid);
