  fid = fopen('appb.dat');
  u = fread(fid,'double');
  status = fclose(fid);

  N = size(u,1);
  NP = 512;   % number of plot point for X direction
  s =  -pi:2*pi/NP:pi;
  app = zeros(NP+1,6);
  for m = 1:N
     app(:,1) = app(:,1) + u(m)*sin(m*s');    % u
     app(:,2) = app(:,2) - u(m)*sin(m*s');    % -u
     app(:,3) = app(:,3) + m*u(m)*cos(m*s');    % u'
     app(:,4) = app(:,4) + u(m)/m*cos(m*s');    % -\int u
     app(:,5) = app(:,5) - m^2*u(m)*sin(m*s');    % u''
     app(:,6) = app(:,6) + u(m)/(m^2)*sin(m*s');    % -\int\int u
  end
  figure;
  set(gcf, 'DefaultLineLineWidth', 2);
  tiles = tiledlayout(3,2);
  title(tiles, 'Approximate solution');
  
  ax1 = nexttile;
  plot(s,app(:,1));
  title('u(x)');

  ax2 = nexttile;
  plot(s,app(:,2));
  title("$-u(x)=\psi''(x)$", 'Interpreter', 'latex');

  ax3 = nexttile;
  plot(s,app(:,3));
  title("$u'(x)$", 'Interpreter', 'latex');

  ax4 = nexttile;
  plot(s,app(:,4));
  title("$\psi'(x)$", 'Interpreter', 'latex');

  ax5 = nexttile;
  plot(s,app(:,5));
  title("$u''(x)$", 'Interpreter', 'latex');

  ax6 = nexttile;
  plot(s,app(:,6));
  title("$\psi(x)$", 'Interpreter', 'latex');

  linkaxes([ax1, ax2, ax3, ax4, ax5, ax6], 'xy');
  axis([-pi pi -Inf Inf]);
  allAxes = findobj(tiles, 'Type', 'axes');
  set(allAxes, 'XGrid', 'on', 'YGrid', 'on');
  set(allAxes, 'FontSize', 12);

