#!/bin/bash

# 检查 clang-format 是否存在
if ! command -v clang-format &>/dev/null; then
    echo "clang-format is not installed. Please install clang-format first."
    exit 1
fi

# 要格式化的目录
DIRECTORIES=("../src" "../test" "../inc")

# 支持的文件扩展名
EXTENSIONS=("*.cpp" "*.c" "*.h")

# 遍历每个目录并格式化文件
for DIR in "${DIRECTORIES[@]}"; do
    for EXT in "${EXTENSIONS[@]}"; do
        find "$DIR" -type f -name "$EXT" | while read -r FILE; do
            echo "Formatting: $FILE"
            clang-format -i "$FILE"
        done
    done
done

echo "Formatting complete."
