# 收集源文件
file(GLOB_RECURSE SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
)

# 添加头文件目录
include_directories(
        ${HEADERS_PATH}
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# 编译成库
add_library(${PROJECT_NAME} SHARED
        ${SOURCES}
)

if (UNIX)
        target_link_libraries(${PROJECT_NAME} PRIVATE
                m
        )
endif ()

# 安装目标文件
install(TARGETS ${PROJECT_NAME}
        RUNTIME DESTINATION ${INSTALL_PATH}/bin
        LIBRARY DESTINATION ${INSTALL_PATH}/bin
        ARCHIVE DESTINATION ${INSTALL_PATH}/bin
)
