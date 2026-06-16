  fid = fopen('../../k100/R10000/appb.dat');
  u1 = fread(fid,'double');
  status = fclose(fid);

  fid = fopen('../../k101/R10000/appb.dat');
  u2 = fread(fid,'double');
  status = fclose(fid);


  % fid = fopen('../../k6/R10000/appb.dat');
  % u4 = fread(fid,'double');
  % status = fclose(fid);

  N = size(u1,1);
  NP = 2048;   % number of plot point for X direction
  s =  -pi:2*pi/NP:pi;
  app1 = zeros(NP+1,2);
  app2 = zeros(NP+1,2);
  % app4 = zeros(NP+1,2);
  for m = 1:N
     app1(:,1) = app1(:,1) + u1(m)*sin(m*s');    % u
     app1(:,2) = app1(:,2) + u1(m)/(m^2)*sin(m*s');    % -\int\int u
     app2(:,1) = app2(:,1) + u2(m)*sin(m*s');    % u
     app2(:,2) = app2(:,2) + u2(m)/(m^2)*sin(m*s');    % -\int\int u
  end
 
  fig = figure;
  width_cm  = 23.5;
  height_cm = 6.5;
  fig.Units = 'centimeters';
  fig.Position(3:4) = [width_cm, height_cm];
  set(gcf, 'DefaultLineLineWidth', 1.);
  tiles = tiledlayout(1,2);
  title(tiles, "$R=10^4$", 'Interpreter', 'latex');
  
  % colors = get(gca, 'ColorOrder');
  colors = [
            0.902, 0.624, 0.000;  % 橙
            0.337, 0.706, 0.914;  % 空色
            0.000, 0.620, 0.451;  % 緑
            0.941, 0.894, 0.259;  % 黄
            0.000, 0.447, 0.698;  % 青
            0.835, 0.369, 0.000;  % 朱
            0.800, 0.475, 0.655;  % 紫
           ];

  ax1 = nexttile;
  plot(ax1, s,app1(:,1), 'Color', colors(1,:));
  hold on;
  plot(ax1, s,app2(:,1), 'Color', colors(2,:));
  hold off;
  grid on;
  title(ax1, "$u(x)$", 'Interpreter', 'latex');
  legend(ax1, "$k=100$", "$k=101$", 'Interpreter', 'latex', 'Location','southeast');
  ax1.FontName = 'Times New Roman';
  ax1.FontSize = 10;

  ax2 = nexttile;
  plot(ax2, s,app1(:,2), 'Color', colors(1,:));
  hold on;
  plot(ax2, s,app2(:,2), 'Color', colors(2,:));
  hold off;
  grid on;
  title(ax2, "$\psi(x)$", 'Interpreter', 'latex');
  legend(ax2, "$k=100$", "$k=101$", 'Interpreter', 'latex', 'Location','southeast');
  ax2.FontName = 'Times New Roman';
  ax2.FontSize = 10;
  
  linkaxes([ax1, ax2], 'xy');
  axis([-pi pi -Inf Inf]);

  exportgraphics(fig, 'figure4.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');