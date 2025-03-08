project(slip_system_benchmark)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(
        ${HEADERS_PATH}
        ${SOURCES_PATH}/include
)

file(GLOB_RECURSE SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
        ${SOURCES_PATH}/*.cpp
)

# 查找GTest
find_package(Benchmark REQUIRED)

add_executable(${PROJECT_NAME}
        ${SOURCES}
)

target_link_libraries(${PROJECT_NAME} PRIVATE
        benchmark
)

if(UNIX)
        target_link_libraries(${PROJECT_NAME} PRIVATE
                m
                benchmark
        )
endif()

# 安装目标文件
install(TARGETS ${PROJECT_NAME}
        RUNTIME DESTINATION ${INSTALL_PATH}/bin
        LIBRARY DESTINATION ${INSTALL_PATH}/lib
        ARCHIVE DESTINATION ${INSTALL_PATH}/lib
)
