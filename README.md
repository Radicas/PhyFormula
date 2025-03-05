RedEDA Computational Geometry Library

## 1. 依赖

- boolean 库
- cJSON 库

### 1.1. boolean

将 libboolean.dll 放在 libs 目录下，仓库会默认提供一个。如果过期或者需要跨平台编译，找 Radica 获取支持。

### 1.2. cJSON

#### 1.2.1. windows

```bash
pacman -S mingw-w64-x86_64-cjson
```

#### 1.2.2. Unix

待补充

#### 1.2.3. MacOS

```bash
brew install cjson
```

## 2. 编译

终端进入项目根目录后：

```bash
mkdir build && cd build
cmake ..
make
```

## 3. 使用

### 3.1. 外部使用

将 inc 目录下的 redcgl 文件夹拷贝至目标项目中，将构建后的 libredcgl.dll 和 libboolean.dll 链接至目标项目即可。

### 3.2. 创建对象

#### 3.2.1. 创建边

```C
#include "redcgl/edge.h"

/* 创建一条边 */
int edge_size = 1;
redcgl::Edge* edges = redcgl::create_empty_edge_arrary(edge_size);

redcgl::Edge* e0 = redcgl::edge_at(edges, 0);
redcgl::edge_set_st_x(e0, 0.0);
redcgl::edge_set_st_y(e0, -5.0);
redcgl::edge_set_et_x(e0, -5.0);
redcgl::edge_set_et_y(e0, 0.0);
redcgl::edge_set_ct_x(e0, 0.0);
redcgl::edge_set_ct_y(e0, 0.0);
redcgl::edge_set_radius(e0, -5.0);
redcgl::edge_set_sa(e0, M_PI_2 + M_PI_4);
redcgl::edge_set_ea(e0, M_PI);

/* 创建多条边 */
edge_size = 4;
redcgl::Edge* edges2 = redcgl::create_empty_edge_arrary(edge_size);
for (int i = 0; i < edge_size; i++)
{
    // specific logic
}

/* 释放内存 */
free_edge(edges)
free_edge(edges2);
```

#### 3.2.2. 创建多边形

```C
#include "redcgl/edge.h"
#include "redcgl/polygon.h"

edge_size = 4;
redcgl::Edge* edges = redcgl::create_empty_edge_arrary(edge_size);
for (int i = 0; i < edge_size; i++)
{
    // 具体处理边的逻辑
}

/* 创建不带孔的多边形 */
redcgl::Polygon* poly_without_hole = redcgl::create_polygon(edges, edge_size, NULL);

/* 使用 */
redcgl::Polygon* next = redcgl::polygon_get_next(poly_without_hole);
redcgl::print_polygon(next);

/* 一旦边传递给多边形，边的释放交由多边形处理 */
redcgl::free_polygon(poly_without_hole);
```

## 4. 开发规范

1. 函数命名采用小写字母加下划线的形式
2. 函数声明要加上完备的 doxygen 格式的注释
3. 算法在源文件的前面需要加上详细的原理说明
4. 函数名根据功能作为前缀
5. 每添加一个函数，必须要有对应的单元测试
6. 头文件函数、源文件函数、单元测试 case 需要保持相对顺序一致
7. 实现不同版本但相同功能的函数，在函数名后面加\_n，比如实现第二种，就\_2，以此类推
8. 如果形参全为基础类型，则加上\_v

头文件要包上以下内容：

```C
#ifdef __cplusplus·
// clang-format on
extern "C"
{
    namespace redcgl
    {
#endif

// 头文件内容

#ifdef __cplusplus
    }  // namespace redcgl
}
// clang-format on
#endif
```

源文件包上以下内容：

```C
#ifdef __cplusplus
namespace redcgl
{
#endif

#ifdef __cplusplus
}  // namespace redcgl
#endif
```
