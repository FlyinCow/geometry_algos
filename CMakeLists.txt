cmake_minimum_required(VERSION 3.30)
project(cmake_tempalte CXX)

add_executable(${PROJECT_NAME})
file(GLOB_RECURSE sources src/*.cpp)
target_sources(${PROJECT_NAME} PRIVATE ${sources})
target_include_directories(${PROJECT_NAME} PRIVATE
    "${CMAKE_SOURCE_DIR}/src"
)