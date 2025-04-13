#pragma once

/*这个头文件专门放控制代码运行方式的宏定义指令，所有宏定义集中在一个文件里方便不同文件共用，同时修改起来也方便*/

//#define print_debug//打印调试
//#define fixed_debug//固定随机矩阵调试
#define covar_eye//协方差矩阵固定为单位阵
//#define new_Hx//查表新方式计算Hx

#define performance_cal//计算性能
//#define real_detect//实数检测，后面应该用不到
#define accel//加速

//#define generate_data
//#define test
#define run_detection
//#define file_op//文件操作