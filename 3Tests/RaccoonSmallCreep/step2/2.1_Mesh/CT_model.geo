// CT模型的gmsh几何文件
// 尺寸参数（单位：m）
lc = 0.0005;  // 网格尺寸，从1mm换算为0.001m

// 外边界尺寸
W = 0.025;    // 总宽度，从25mm换算为0.025m
H = 0.024;    // 总高度，从24mm换算为0.024m
w = 0.020;    // 内部宽度，从20mm换算为0.020m
h = 0.0071;   // 内部高度到圆孔中心的距离，从7.1mm换算为0.0071m
hole_y_top = H/2+h;
hole_y_bottom = H/2-h;
// 圆孔参数
d = 0.006;     // 圆孔直径，从6mm换算为0.006m
r = d/2;       // 圆孔半径

// 圆孔水平位置（距离左边界）
hole_x = (W - w);  // 位于左侧

// 缺口参数
tip_x = 0.0085+hole_x;  // 缺口深度，从8.5mm换算为0.0085m
notch_width = 0.001;  // 缺口宽度，从1mm换算为0.001m
//预制裂纹长度
a0 = 0.010+hole_x;  // 从10mm换算为0.010m
precrack_length = a0-tip_x;

// 修改缺口几何为更真实的预制裂纹形状
// 计算裂纹尖端的y坐标，使其与上圆孔中心的垂直距离是7.1mm
crack_tip_y = H/2;  // 距离上圆孔中心7.1mm

// 重新定义外边界点，包含缺口
Point(1) = {0, 0, 0, lc};
Point(2) = {W, 0, 0, lc};
Point(3) = {W, H, 0, lc};
Point(4) = {0, H, 0, lc};

// 缺口点（从左边界开始）
Point(15) = {0, crack_tip_y + notch_width/2, 0, lc/4};  // 缺口上边界点
Point(16) = {0, crack_tip_y - notch_width/2, 0, lc/4};  // 缺口下边界点

// 缺口中间点（略微扩展）
Point(17) = {tip_x*0.95, crack_tip_y + notch_width/2, 0, lc/4};
Point(18) = {tip_x*0.95, crack_tip_y - notch_width/2, 0, lc/4};

// 缺口尖端（裂纹前沿）- 与上圆孔中心垂直距离7.1mm
Point(19) = {tip_x, crack_tip_y + 1e-6, 0, lc/6};
Point(20) = {tip_x, crack_tip_y - 1e-6, 0, lc/6};
// 预制裂纹定义
// 预制裂纹长度：从缺口尖端到圆孔中心的x方向距离
precrack_end_x = tip_x + precrack_length;  // 预制裂纹终点x坐标

// 预制裂纹终点的上下表面
Point(23) = {precrack_end_x, crack_tip_y + 1e-6, 0, lc/6};  // 预制裂纹终点上表面，从1e-6调整为1e-5
Point(24) = {precrack_end_x, crack_tip_y - 1e-6, 0, lc/6};  // 预制裂纹终点下表面，从1e-6调整为1e-5

// 重新定义外边界线，包含缺口
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 15};   // 从左上角到缺口上边界
Line(13) = {15, 17};  // 缺口上边界
Line(14) = {17, 19};  // 缺口上边界到裂纹尖端左
Line(15) = {19, 23};  // 裂纹尖端左到裂纹尖端右
Line(16) = {23, 24};  // 裂纹尖端右到裂纹尖端右
Line(17) = {24, 20};   // 从裂纹尖端右到裂纹尖端左

// 预制裂纹线（上表面和下表面）
Line(18) = {20, 18};  // 从裂纹尖端左到缺口下边界
Line(19) = {18, 16};  // 从缺口下边界到外圈交点
Line(20) = {16, 1};  // 从缺口下边界到外圈交点
// 修改后的外边界循环，在裂纹尖端闭合
Line Loop(1) = {1, 2, 3, 4, 13, 14, 15, 16, 17,18,19,20};

// 上圆孔中心 (距离左边界 hole_x, 距离底边 hole_y_top)
Point(5) = {hole_x, hole_y_top, 0, lc/2};

// 上圆孔的边界点（简化为4个基本方向）
Point(6) = {hole_x + r, hole_y_top, 0, lc/2};  // 右
//Point(7) = {hole_x, hole_y_top + r, 0, lc/2};  // 上
Point(8) = {hole_x - r, hole_y_top, 0, lc/2};  // 左
Point(9) = {hole_x, hole_y_top - r, 0, lc/2};  // 下

// 上圆孔的圆弧（简化为4段）
//Circle(5) = {6, 5, 7};   // 右到上
//Circle(6) = {7, 5, 8};   // 上到左
Circle(7) = {8, 5, 9};   // 左到下
Circle(8) = {9, 5, 6};   // 下到右
Line(50) = {6,8};
// 上圆孔循环
Line Loop(2) = {7,8,50};

// 下圆孔中心 (距离左边界 hole_x, 距离底边 hole_y_bottom)
Point(10) = {hole_x, hole_y_bottom, 0, lc/2};

// 下圆孔的边界点（简化为4个基本方向）
Point(11) = {hole_x + r, hole_y_bottom, 0, lc/2};  // 右
Point(12) = {hole_x, hole_y_bottom + r, 0, lc/2};  // 上
Point(13) = {hole_x - r, hole_y_bottom, 0, lc/2};  // 左
//Point(14) = {hole_x, hole_y_bottom - r, 0, lc/2};  // 下

// 下圆孔的圆弧（简化为4段）
Circle(9) = {11, 10, 12};  // 右到上
Circle(10) = {12, 10, 13}; // 上到左
//Circle(11) = {13, 10, 14}; // 左到下
//Circle(12) = {14, 10, 11}; // 下到右
Line(51) = {13,11};
// 下圆孔循环
Line Loop(3) = {51, 9, 10};

// 定义平面表面（从修改后的外边界减去两个圆孔和预制裂纹）
Plane Surface(1) = {1, 2, 3};

// 定义边界ID
Physical Line("bottom") = {1};
Physical Line("right") = {2};
Physical Line("top") = {3};
Physical Line("left") = {4, 17};  // 左边界包含缺口两侧
Physical Line("upper_hole") = {50};  // 上圆孔完整边界
Physical Line("lower_hole") = {51};  // 下圆孔完整边界
Physical Line("notch") = {13, 14, 15, 16};  // 缺口边界
Physical Line("precrack_upper") = {18};  // 预制裂纹上表面
Physical Line("precrack_lower") = {20};  // 预制裂纹下表面
Physical Line("precrack_tip") = {19};    // 预制裂纹尖端连接线

// 定义受力点
Physical Point("upper_load_point") = {7};  // 上圆孔顶部受力点
Physical Point("lower_load_point") = {14};  // 下圆孔底部受力点

// 定义区域ID
Physical Surface("body") = {1};

// 网格参数
Mesh.Algorithm = 6;  // 回到稳定的Frontal-Delaunay算法
Mesh.ElementOrder = 1;
Mesh.CharacteristicLengthMax = lc;
Mesh.CharacteristicLengthMin = lc/5;  // 放宽最小网格尺寸

// 简化的四边形网格设置
Mesh.RecombineAll = 1;  // 将三角形重新组合成四边形
Mesh.RecombinationAlgorithm = 1;  // 使用标准重组算法

// 定义裂纹尖端网格加密
// 裂纹尖端周围的网格尺寸
lc_tip = lc/5;  // 裂纹尖端最细网格尺寸
lc_transition = lc;  // 过渡区域网格尺寸

// 重新设置裂纹尖端相关点的网格尺寸
Characteristic Length {19} = lc_tip;  // 裂纹尖端
Characteristic Length {23, 24} = lc_tip;  // 预制裂纹终点

// 简化的长方体加密区域
// 内部细网格区域
Field[1] = Box;
Field[1].XMin = tip_x - 0.002;           // 扩大边界，减少网格梯度
Field[1].XMax = precrack_end_x + 0.002;  // 扩大边界，减少网格梯度
Field[1].YMin = crack_tip_y - 0.002;     // 扩大边界，减少网格梯度
Field[1].YMax = crack_tip_y + 0.002;     // 扩大边界，减少网格梯度
Field[1].ZMin = 0;                  // 长方体前边界
Field[1].ZMax = 0;                   // 长方体后边界
Field[1].VIn = lc_tip;               // 长方体内部网格尺寸
Field[1].VOut = lc;                  // 长方体外部网格尺寸

// 设置背景网格场
Background Field = 1;

// 添加表面重组为四边形
Recombine Surface {1}; 