function PlotDesignRatio(eta, b, L, x, y, tit)

% The function PlotDesignRatio plots the design factor as surface on 
% CLT-Concrete
% Input data
%   eta   = matrix of design ratio by coordinates x and y
%   b     = width, m
%   L     = span, m
%   x     = vector of coordinates along x-axis, only positive ones, m
%   y     = vector of coordinates along y-axis, m
%   tit   = title of the plot

% % Reshape x and y to ensure they are column vectors
%     x = reshape(x, [], 1);
%     y = reshape(y, 1, []);

% Find coordenates of maximum value
% [maxValue, maxIndex] = max(eta(:));
% [maxRow, maxCol] = ind2sub(size(eta), maxIndex);
% maxX = x(maxRow, maxCol);
% maxY = y(maxRow, maxCol);
% maxZ = eta(maxRow, maxCol);
% etamax = maxValue;
% widthWindow = 200*L;  
% heithWigndow = 800*b;   
[X, Y] = meshgrid(x, y);

figure(1)
%figure('Position', [100, 100, widthWindow, heigthWindow])
hold on
surf(X, Y, eta', 'EdgeColor', 'none', 'FaceColor', 'interp');
surf(-X, Y, eta', 'EdgeColor', 'none', 'FaceColor', 'interp');
view(45, 30);
xlabel('Length (m)');
ylabel('Width (m)');
zlabel('\eta');
title(tit);
xlim([-L/2 L/2]);
ylim([-b/2 b/2]);
zlim([0 max(eta(:))]);
grid on
box on
colormap(jet);
caxis([0.0, 1.0]);
colorbar; 
axis equal;
set(gcf, 'WindowState', 'maximized');

%if maxValue<=1
%    colortext = 'blue';
%else
%    colortext = 'red';
%end
%txt = text(-maxX, -maxY, maxZ, sprintf('Max Value: %.2f', maxValue), ...
%    'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom',...
%    'FontSize', 10, 'Color', colortext,'BackgroundColor', 'white');
%set(txt, 'EdgeColor', colortext, 'BackgroundColor', [1, 1, 1, 0.5]);
end