param (
    [string]$srcDir = "../src",
    [string]$testDir = "../test"
)

# 检查 clang-format 是否存在
if (-not (Get-Command clang-format -ErrorAction SilentlyContinue)) {
    Write-Host "clang-format is not installed. Please install clang-format first." -ForegroundColor Red
    exit 1
}

# 要格式化的目录
$directories = @($srcDir, $testDir)

# 支持的文件扩展名
$extensions = @("*.cpp", "*.c", "*.h")

# 遍历每个目录并格式化文件
foreach ($dir in $directories) {
    foreach ($ext in $extensions) {
        Get-ChildItem -Path $dir -Recurse -Filter $ext | ForEach-Object {
            Write-Host "Formatting: $_.FullName"
            clang-format -i $_.FullName
        }
    }
}

Write-Host "Formatting complete."
