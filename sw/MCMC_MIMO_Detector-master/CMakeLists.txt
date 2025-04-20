cmake_minimum_required(VERSION 3.10)  # 设定最低CMake版本
project(MCMC_MIMO_Detector C)  # 设定项目名称，使用C语言

set(CMAKE_C_STANDARD 99)  # 设定C标准为C99

# 查找当前目录下所有C文件
file(GLOB SOURCES "*.c")  # 获取所有 .c 文件
file(GLOB HEADERS "*.h")  # 获取所有 .h 文件

# 添加可执行目标
add_executable(MCMC_MIMO_Detector ${SOURCES} ${HEADERS})

# 设置输出目录
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})