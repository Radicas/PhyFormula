project(redcgl_benchmark)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(
        ${HEADERS_PATH}
        ${SOURCES_PATH}/include
        ${SOURCES_PATH}/cjson
        "/usr/local/opt/google-benchmark/include"
)

file(GLOB_RECURSE SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/*.c
        ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
        ${SOURCES_PATH}/*.c
)

# 查找GTest
find_package(Benchmark REQUIRED)

message(STATUS "project name " ${PROJECT_NAME})

add_executable(${PROJECT_NAME}
        ${SOURCES}
)

# 添加库搜索路径
target_link_directories(${PROJECT_NAME} PRIVATE
        ${LIBS_PATH}
        $ENV{REDEDA_PLUS_DEV_DIST}
        "/usr/local/opt/google-benchmark/lib"
)

if(UNIX)
        target_link_libraries(${PROJECT_NAME} PRIVATE
                m
                benchmark
        )
else()
        target_link_libraries(${PROJECT_NAME} PRIVATE
                benchmark
        )
endif()

# 安装目标文件
install(TARGETS ${PROJECT_NAME}
        RUNTIME DESTINATION ${INSTALL_PATH}/bin
        LIBRARY DESTINATION ${INSTALL_PATH}/lib
        ARCHIVE DESTINATION ${INSTALL_PATH}/lib
)
