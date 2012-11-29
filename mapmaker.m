close all;
clear all;
%% Load data
% source: http://pds-geosciences.wustl.edu/missions/mgs/megdr.html
f=fopen('megt90n000cb.img');
mapinfo = struct('west',0,'east',360,'north',90,'south',-90);
W = 1440;
H = 720;
jump = 5;
% source: http://toolserver.org/~geohack/geohack.php?pagename=Mars_Science_Laboratory&params=4.5895_S_137.4417_E_globe:Mars
curiosity = struct('north',-4.5895,'east',137.4417);
A=flipud(fread(f,[W,H],'*int16',0,'b')');
fclose(f);
A=A(1:jump:end,1:jump:end);
W = W/jump;
H = H/jump;
x0 = [
    floor((curiosity.east - mapinfo.west) / (mapinfo.east - mapinfo.west) * W)
    floor((curiosity.north - mapinfo.south) / (mapinfo.north - mapinfo.south) * H)];
figure(1);
imagesc(A);
figure(gcf);
set(gca,'YDir','normal');
hold on;
plot(x0(1),x0(2),'sr');
hold off;

%% Policy Controls:
%    812
%    703
%    654

%% Init grid (direct) and costs (euclidean)
cost = zeros(H,W);
control = zeros(H,W);
for i = 1:W % E/W
    for j = 1:H % N/S
        ci = mod(floor(atan2(x0(2)-j, x0(1)-i) / 2 / pi * 8 + 4.5),8)+1;
        costi = hypot(j-x0(2), i-x0(1));
        control(j,i) = ci;
        cost(j,i) = costi;
    end
end

%% Show stats
figure(2);
subplot(2,1,1);
imagesc(control);
figure(gcf);
set(gca,'YDir','normal');
subplot(2,1,2);
imagesc(cost);
figure(gcf);
set(gca,'YDir','normal');

%% Prep iteration
multiplier = [1,sqrt(2),1,sqrt(2),1,sqrt(2),1,sqrt(2)]';
nextx = zeros(8,H,W);
u_cost = zeros(8,H,W);
for i = 1:W
    for j = 1:H
        h = A(j,i);
        nx_xs = [
            i,            mod(j-2,H)+1 % down
            mod(i-2,W)+1, mod(j-2,H)+1 % down left
            mod(i-2,W)+1, j % left
            mod(i-2,W)+1, mod(j,H)+1 % up left
            i,            mod(j,H)+1 % up
            mod(i,W)+1,   mod(j,H)+1 % up right
            mod(i,W)+1,   j % right
            mod(i,W)+1,   mod(j-2,H)+1 % down right
            ];
        nx_x = nx_xs(:,2)+(nx_xs(:,1)-1)*H;
        nextx(:,j,i) = nx_x;
        next_h = A(nx_x);
        % small linear cost fcn gives nice cost map
        %u_cost(:,j,i) = (abs(double(next_h-h)/200) + 1).*multiplier;
        
        % medium square cost fcn gives nicer flow field
        u_cost(:,j,i) = ((double(next_h-h)/1000).^2 + 1).*multiplier;
    end
end

%% Value iteration of the actual costs
for iter = 1:50
    prevcost = cost;
    for i = 1:W
        for j = 1:H
            if x0(2) == j && x0(1) == i
                control(j,i) = 0;
                cost(j,i) = 0;
            else
                costs = u_cost(:,j,i) + prevcost(nextx(:,j,i));
                [costi, ci] = min(costs);
                control(j,i) = ci;
                cost(j,i) = costi;
            end
        end
    end
    if mod(iter,jump)==0 || true
        %% Show stats
        figure(2);
        subplot(2,1,1);
        imagesc(control);
        figure(gcf);
        set(gca,'YDir','normal');
        subplot(2,1,2);
        imagesc(cost);
        figure(gcf);
        set(gca,'YDir','normal');
    end
end
