
set(EXEC_NAME height_detection.cpp)
set(OUTPUT_NAME find_objects)
cmake_minimum_required(VERSION 3.10.2)
project( ${OUTPUT_NAME})

add_subdirectory(functions)
include_directories(functions)

find_package( OpenCV REQUIRED )
add_executable( ${OUTPUT_NAME}  ${EXEC_NAME})
target_link_libraries(${OUTPUT_NAME} ${OpenCV_LIBS} functions make_groups)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_BUILD_TYPE RelWithDebInfo)
