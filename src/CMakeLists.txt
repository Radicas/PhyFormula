# 收集源文件
file(GLOB_RECURSE SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
)

# 添加头文件目录
include_directories(
        ${HEADERS_PATH}
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/cjson
)

# 编译成库
add_library(${PROJECT_NAME} SHARED
        ${SOURCES}
)

# 添加库搜索路径
target_link_directories(${PROJECT_NAME} PRIVATE
        ${LIBS_PATH}
)

if (UNIX)
        target_link_libraries(${PROJECT_NAME} PRIVATE
                m
        )
endif ()

# 获取所有依赖的动态库
file(GLOB DYLIBS ${LIBS_PATH}/*)

# 将动态库安装到bin目录
install(FILES ${DYLIBS} DESTINATION
        ${INSTALL_PATH}/bin
)

# 安装目标文件
install(TARGETS ${PROJECT_NAME}
        RUNTIME DESTINATION ${INSTALL_PATH}/bin
        LIBRARY DESTINATION ${INSTALL_PATH}/bin
        ARCHIVE DESTINATION ${INSTALL_PATH}/bin
)

# ####### 自定义安装 ########
# # 获取当前动态库
# file(GLOB_RECURSE DLL ${CMAKE_HOME_PATH}/install/bin/*.dll)

# message(STATUS "REDEDA_DEV_DIST path:" ${REDEDA_DEV_DIST})

# if(REDEDA_DEV_DIST)
# # 将当前动态库安装到DEV_DIST
# install(FILES ${DLL} DESTINATION
# ${REDEDA_DEV_DIST}
# )

# # 获取头文件
# set(INC ${CMAKE_HOME_PATH}/inc/redcgl)

# # 将头文件安装到DEV_DIST
# install(DIRECTORY ${INC} DESTINATION
# ${REDEDA_DEV_DIST}/rededa/include/rededa
# )
# endif ()
# #########################