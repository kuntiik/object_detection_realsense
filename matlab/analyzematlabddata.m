close all
HEIGHT = 480;
WIDTH = 848;
cx = 424.9;
cy = 244.22;
fx =428.8;
fy = 428.8;
n_pix = HEIGHT*WIDTH;
b = importdata('data.txt');
c = importdata('corrected.txt');
depth = [];
pic_norm = [0.601 0.22 0.77]';
p_norm_mm = pic_norm*1000;
d = 1.066;
d_mm = d*1000;
[XX,YY] = ndgrid(-2000:20:2000,-2000:20:2000);
z = (-p_norm_mm(1)*XX - p_norm_mm(2)*YY - d_mm)/p_norm_mm(3);
for i = 1:480
    %size(b{i})
    depth = vertcat(depth,str2num(b{i}));
end
%[X,Y] = meshgrid(1:2:848,1:2:480);
%Z = depth(Y,X);
hold on
%mesh(depth)
surf(XX,YY,z)
%surf(-depth)
%surface(depth)
rp = zeros(480*848,3);
for i = 1:HEIGHT
    for j = 1:WIDTH
        rp(i*480+j,1) = (i - cx)/fx*depth(i,j);
        rp(i*480+j,2) = (j - cy)/fy*depth(i,j);
        rp(i*480+j,3) = depth(i,j);
    end
end
%scatter3(rp(:,:,1),rp(:,:,2),rp(:,:,3),"o");
pcshow(-rp);
set(gca, 'XColor', [0.15 0.15 0.15], 'YColor', [0.15 0.15 0.15], 'ZColor', [0.15 0.15 0.15])
set(gca,'color','w');
hold off
