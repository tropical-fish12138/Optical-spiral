clear;
clc;
close all;

%% ================================================================= %%
%% 1. 基础物理常数与仿真参数设置
%% ================================================================= %%
% 单位换算
um = 1e-6;
nm = 1e-9;
mm = 1e-3;
deg = pi/180;

% 网格参数
N        = 800;             % 初始光场网格尺寸
N1       = 800;             % 输出光场网格尺寸 (为 padarray 预留)
N_slice  = 800;             % 光阑尺寸参考
Ts       = 8 * um;          % 空间采样间隔 (像素尺寸)

% 光学参数
wavelength = 0.6328 * um;   % 入射光波长
f          = 10 * mm;       % 衍射距离参数 (用于初始 classical_diffraction 验证)
k          = 2 * pi / wavelength; % 空间波数

%% ================================================================= %%
%% 2. 空间坐标网格初始化
%% ================================================================= %%
% 生成源平面坐标网格 (尺寸 N x N)
x_vec = (-N/2 + 0.5 : 1 : N/2 - 0.5) * Ts;
[x, y] = meshgrid(x_vec);
[theta, r] = cart2pol(x, y);

% 生成观察平面坐标网格 (尺寸 N1 x N1)
x1_vec = (-N1/2 + 0.5 : 1 : N1/2 - 0.5) * Ts;
[x1, y1] = meshgrid(x1_vec);

%% ================================================================= %%
%% 3. 孔径光阑与入射光场设计
%% ================================================================= %%
% 3.1 方形孔径光阑生成 (使用矢量化逻辑运算代替低效的 for 循环)
r1 = 0.5 * N_slice * Ts;  
aperture = double(abs(x) <= r1 & abs(y) <= r1); % 在范围内的透过率为1，否则为0
aperture = padarray(aperture, [(N1-N)/2, (N1-N)/2]);

% 3.2 环形高斯入射光束 (消除贝塞尔光束的无限旁瓣)
r_ring = 0.15 * mm;   % 环形光束的中心半径 (匹配 Z 轴观察窗口)
w_ring = 0.1 * mm;    % 环形光束的圆环宽度(束腰)

% 计算环形高斯光振幅分布
input_light = exp(-((r - r_ring).^2) ./ (w_ring^2));
input_light = padarray(input_light, [(N1-N)/2, (N1-N)/2]);

% 绘制入射光振幅图
figure('Name', '入射光场', 'Color', 'white');
mesh(x1.*1e3, y1.*1e3, input_light);
title("入射光束振幅 (环形高斯光)");
xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Amplitude');

%% ================================================================= %%
%% 4. 双涡旋贝塞尔相位设计
%% ================================================================= %%
% 光束参数设置
l1 = 1;                 % 光束1的拓扑电荷
l2 = -1;                % 光束2的拓扑电荷
alpha1 = 0.2 * deg;     % 等效轴锥镜1的折射角/圆锥角 (决定径向波矢)
alpha2 = 0.4 * deg;     % 等效轴锥镜2的折射角/圆锥角 (决定径向波矢)

kr1 = k * sin(alpha1);  % 光束1向内的径向波矢
kr2 = k * sin(alpha2);  % 光束2向内的径向波矢

% 理想纯相位公式：径向汇聚相位 + 角向涡旋相位
Phase_ideal1 = -kr1 * r + l1 * theta;
Phase_ideal2 = -kr2 * r + l2 * theta;

% 将相位包裹到 0~2pi 之间，并进行尺寸匹配补零
Phase_wrapped1 = padarray(mod(Phase_ideal1, 2*pi), [(N1-N)/2, (N1-N)/2]);
Phase_wrapped2 = padarray(mod(Phase_ideal2, 2*pi), [(N1-N)/2, (N1-N)/2]);

%% ================================================================= %%
%% 5. 初始单切面衍射计算 (验证用)
%% ================================================================= %%
% 叠加总输入光场 = 入射光振幅 * 孔径 * (两束光的复振幅之和)
% 注：避免使用 MATLAB 内置函数名 'input'，修改为 'input_field'
input_field = input_light .* aperture .* (exp(1i*Phase_wrapped1) + exp(1i*Phase_wrapped2));

% 调用外部函数计算单平面衍射
output_field = classical_diffraction(input_field, x1, y1, f, wavelength, 1/Ts, 1/Ts, "a");

% 计算光强并归一化
output_int = abs(output_field).^2;
output_int_norm = output_int ./ max(output_int(:));

% 绘制单切面仿真结果图
figure('Name', '单切面衍射验证', 'Color', 'white');
imshow(output_int_norm, []);
title("仿真结果图");

%% ================================================================= %%
%% 6. 不同 Z 轴传播距离下的横截面演化 (2D 切片图)
%% ================================================================= %%
z_observe_array = [10, 20, 30, 40, 50, 60] * mm;

figure('Name', '不同 Z 轴传播距离下的光强演化对比', 'Color', 'white', 'Position', [100 100 1200 700]);

for idx = 1:length(z_observe_array)
    current_z = z_observe_array(idx);
    
    % 计算当前距离 Z 的衍射光场
    output_z = classical_diffraction(input_field, x1, y1, current_z, wavelength, 1/Ts, 1/Ts, "a");
    
    % 归一化当前平面光强
    int_z = abs(output_z).^2;
    int_z_norm = int_z ./ max(int_z(:));
    
    % 子图绘制
    subplot(2, 3, idx);
    imagesc(x1(1,:)*1e3, y1(:,1)*1e3, int_z_norm);
    colormap(hot); % 激光热力伪彩
    axis square;
    
    % 限制视场范围，放大中心区域观测 (±0.8mm)
    xlim([-0.8 0.8]);
    ylim([-0.8 0.8]);
    
    title(sprintf('Z = %d mm', current_z * 1000), 'FontSize', 12);
    xlabel('X (mm)');
    ylabel('Y (mm)');
    
    % 添加参考十字虚线，辅助观察光瓣旋转
    hold on;
    plot([-1.5 1.5], [0 0], 'w:', 'LineWidth', 1);
    plot([0 0], [-1.5 1.5], 'w:', 'LineWidth', 1);
    hold off;
end
sgtitle(sprintf('不同传播距离下光学弹簧的横截面演化 (l_1=%d, l_2=%d)', l1, l2), 'FontSize', 16, 'FontWeight', 'bold');

%% ================================================================= %%
%% 7. 生成并绘制三维圆柱螺旋光学弹簧结构 (Volume Rendering)
%% ================================================================= %%
disp('正在计算三维光场，这可能需要一些时间，请稍候...');

% 7.1 三维采样参数设置
z_start    = 10 * mm;      
z_end      = 40 * mm;        
num_slices = 60;        % Z 轴采样切片数
Z_3D_array = linspace(z_start, z_end, num_slices);

% 7.2 确定感兴趣区域 (ROI) 以节省内存
roi_radius_mm = 0.6;    
roi_pixels = round((roi_radius_mm * mm) / Ts);
center_idx = N1 / 2;
ROI_range = (center_idx - roi_pixels) : (center_idx + roi_pixels);

% 7.3 初始化三维体积矩阵并逐层填充
Vol_Int = zeros(length(ROI_range), length(ROI_range), num_slices);

for idx = 1:num_slices
    current_z = Z_3D_array(idx);
    output_z = classical_diffraction(input_field, x1, y1, current_z, wavelength, 1/Ts, 1/Ts, "a");
    
    int_z = abs(output_z).^2;
    int_z_roi = int_z(ROI_range, ROI_range);
    
    % 局部平面归一化后存入体积矩阵
    Vol_Int(:, :, idx) = int_z_roi ./ max(int_z_roi(:));
end

% 7.4 构建三维物理坐标网格 (单位：mm)
[X_3D, Y_3D, Z_3D] = meshgrid(x1(1, ROI_range)*1e3, y1(ROI_range, 1)*1e3, Z_3D_array*1e3);

% 7.5 三维可视化渲染
figure('Name', '三维光学弹簧结构', 'Color', 'white', 'Position', [150 150 800 800]);
iso_threshold = 0.4; 

% 提取等值面数据并绘制
[faces, verts] = isosurface(X_3D, Y_3D, Z_3D, Vol_Int, iso_threshold);
p = patch('Faces', faces, 'Vertices', verts);

% 计算法线以确保光照正确反射
isonormals(X_3D, Y_3D, Z_3D, Vol_Int, p);

% 材质与色彩设定 (珊瑚红反光材质)
p.FaceColor = [1.0 0.35 0.25]; 
p.EdgeColor = 'none';          
p.AmbientStrength = 0.4;       
p.DiffuseStrength = 0.8;       
p.SpecularStrength = 0.9;      
p.SpecularExponent = 25;       

% 光照与视角设定
view(-40, 30);                 
lighting phong;                
camlight('headlight');         
camlight('left');              

% 坐标轴格式化与比例压缩
axis([-roi_radius_mm, roi_radius_mm, -roi_radius_mm, roi_radius_mm, z_start*1e3, z_end*1e3]);
grid on; box on;
daspect([1 1 12]); % 调整显示长宽比，防止 Z 轴过度拉伸

xlabel('X 轴 (mm)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Y 轴 (mm)', 'FontSize', 12, 'FontWeight', 'bold');
zlabel('Z 轴传播距离 (mm)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('涡旋光干涉的三维螺旋结构 (等值面 = %.1f)', iso_threshold), 'FontSize', 16, 'FontWeight', 'bold');

rotate3d on; 
disp('三维光场计算与渲染完成！');