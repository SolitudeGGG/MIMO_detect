MHGD_accel_hw.cpp包含顶层函数和部分函数的硬件实现版本
main_hw.cpp包含testbench
matrix_mul.sv 提供参数化的定点矩阵/标量乘法模块：
- 支持方阵或矩阵与行/列向量相乘，`scalar_mode` 可做矩阵与标量的逐元素乘法；
- `INT_W`/`FRAC_W` 参数化定点位宽；
- 时间复用乘累加避免宽位宽脉动阵列的 LUT 暴涨，`SERIAL_CHUNK`/`CHUNK_W` 将乘法器切块到 DSP 友好的宽度以降低 DSP/LUT 消耗。
