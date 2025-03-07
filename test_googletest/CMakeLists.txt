project(formula_test)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(
        ${HEADERS_PATH}
        ${SOURCES_PATH}/include
        ${SOURCES_PATH}/cjson
)

file(GLOB_RECURSE SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
        ${SOURCES_PATH}/*.cpp
)

# 查找GTest
find_package(GTest)
include(GoogleTest)

add_executable(${PROJECT_NAME}
        ${SOURCES}
)

# 添加库搜索路径
target_link_directories(${PROJECT_NAME} PRIVATE
        ${LIBS_PATH}
)

target_link_libraries(${PROJECT_NAME} PRIVATE
        GTest::gtest
)

if(UNIX)
        target_link_libraries(${PROJECT_NAME} PRIVATE
                m
        )
endif()

# 安装目标文件
install(TARGETS ${PROJECT_NAME}
        RUNTIME DESTINATION ${INSTALL_PATH}/bin
        LIBRARY DESTINATION ${INSTALL_PATH}/lib
        ARCHIVE DESTINATION ${INSTALL_PATH}/lib
)