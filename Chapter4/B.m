I = imread("amsterdam.bmp");
% imshow(I)
% size(I,[1,2])
% center = (700/2,734/2), top = (700/2,0)
% from " Nieuwmarkt ”(near the center of the picture) to “ Station Amsterdam Centraal ”
% nimg = insertShape(I, 'Line', [350 367 350 0], 'Color', 'r');
height = size(I,1);
width = size(I,2);
%%%
imshow(I);
hold on;  % Enable hold to overlay the star on the image

line_vertices = [350, 367; 350, 0];
star_vertices = [350, 400; 370, 500;300, 440; 400, 440; 330, 500;350,400];

plot(line_vertices(:, 1), line_vertices(:, 2), 'r', 'LineWidth', 2);
plot(star_vertices(:, 1), star_vertices(:, 2), 'm', 'LineWidth', 2);

hold off; % Release the hold
%%%

%%%
for color = 1:3
    [~, max_indices] = max(I(:,:,color), [], 'all', 'linear');    
    [row, col] = ind2sub([height, width], max_indices);

    random_index = randi(length(row));
    max_row = row(random_index);
    max_col = col(random_index);
    
    max_color = I(max_row, max_col, :);
    
    circle_radius = 10;
%     disp(max_row-circle_radius)
%     hold on;
%     rectangle('Position', [350, 367, 2*circle_radius, 2*circle_radius],'Curvature', [1, 1], 'EdgeColor', max_color/255, 'LineWidth', 5);
%     hold off;
    hold on;
    center = [350; 367];
    radius = 50;
    x = center(1) - radius;
    y = center(2) - radius;
    w = radius * 2;
    h = radius * 2;
    pos = [x y w h];
    cur = [1 1];
%     imshow(I)
    rectangle('Position',pos,'Curvature',cur,'EdgeColor',max_color,'LineWidth',2)
    hold off;
end
%%%