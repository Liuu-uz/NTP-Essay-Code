#!/bin/bash
# GPU进程清理和环境检查脚本

echo "🔍 GPU状态检查和清理工具"
echo "=================================="

# 1. 检查当前GPU使用情况
echo "📊 当前GPU使用情况:"
if command -v nvidia-smi &> /dev/null; then
    nvidia-smi --query-gpu=index,name,memory.used,memory.total,utilization.gpu --format=csv,noheader,nounits
    echo ""
    
    echo "🔍 GPU进程详情:"
    nvidia-smi pmon -c 1 2>/dev/null || echo "   没有检测到GPU进程"
    echo ""
else
    echo "   nvidia-smi 不可用，无法检查GPU状态"
fi

# 2. 查找可能占用GPU的Python进程
echo "🐍 查找Python进程（可能占用GPU）:"
python_procs=$(ps aux | grep python | grep -v grep | grep -v $$ | head -10)
if [ -n "$python_procs" ]; then
    echo "$python_procs"
    echo ""
    
    # 查找包含protenix、torch、cuda的进程
    echo "🎯 查找相关进程 (protenix/torch/cuda):"
    ps aux | grep -E "(protenix|torch|cuda)" | grep -v grep | head -5
else
    echo "   没有发现Python进程"
fi
echo ""

# 3. 检查当前环境变量
echo "🌍 当前CUDA相关环境变量:"
echo "   CUDA_VISIBLE_DEVICES = '${CUDA_VISIBLE_DEVICES}'"
echo "   CUDA_DEVICE_ORDER = '${CUDA_DEVICE_ORDER}'"
echo "   PYTORCH_CUDA_ALLOC_CONF = '${PYTORCH_CUDA_ALLOC_CONF}'"
echo "   TORCH_CUDA_ENABLE = '${TORCH_CUDA_ENABLE}'"
echo "   FORCE_CPU = '${FORCE_CPU}'"
echo ""

# 4. 提供清理选项
echo "🧹 清理选项:"
echo "1. 清理环境变量并强制CPU模式"
echo "2. 杀死可疑的Python进程"
echo "3. 重启Python环境"
echo "4. 查看详细的GPU内存使用"
echo ""

read -p "选择操作 (1-4, 回车跳过): " choice

case $choice in
    1)
        echo "🛡️ 设置强制CPU环境..."
        export CUDA_VISIBLE_DEVICES=""
        export TORCH_CUDA_ENABLE=0
        export FORCE_CPU=1
        export DS_SKIP_CUDA_CHECK=1
        
        echo "✅ 环境变量已设置:"
        echo "   CUDA_VISIBLE_DEVICES = '${CUDA_VISIBLE_DEVICES}'"
        echo "   TORCH_CUDA_ENABLE = '${TORCH_CUDA_ENABLE}'"
        echo "   FORCE_CPU = '${FORCE_CPU}'"
        echo ""
        
        # 测试PyTorch
        echo "🧪 测试PyTorch CUDA状态:"
        python3 -c "
import torch
print(f'CUDA available: {torch.cuda.is_available()}')
print(f'CUDA device count: {torch.cuda.device_count()}')
if torch.cuda.is_available():
    print('⚠️  PyTorch仍能访问CUDA')
else:
    print('✅ PyTorch将使用CPU')
" 2>/dev/null || echo "❌ PyTorch测试失败"
        ;;
        
    2)
        echo "🔍 查找可疑进程..."
        suspicious_procs=$(ps aux | grep -E "(protenix|python.*torch)" | grep -v grep | awk '{print $2}')
        
        if [ -n "$suspicious_procs" ]; then
            echo "发现可疑进程: $suspicious_procs"
            echo "是否要杀死这些进程? (y/N)"
            read -p "> " kill_confirm
            
            if [[ $kill_confirm =~ ^[Yy]$ ]]; then
                echo "$suspicious_procs" | xargs kill -9
                echo "✅ 进程已清理"
            fi
        else
            echo "✅ 没有发现可疑进程"
        fi
        ;;
        
    3)
        echo "🔄 重启Python环境建议:"
        echo "   1. 退出当前shell: exit"
        echo "   2. 重新登录或新开终端"
        echo "   3. 重新激活conda/virtualenv环境"
        echo "   4. 设置CPU环境变量再运行"
        ;;
        
    4)
        echo "🔍 详细GPU内存使用:"
        if command -v nvidia-smi &> /dev/null; then
            nvidia-smi
        else
            echo "nvidia-smi 不可用"
        fi
        ;;
        
    *)
        echo "跳过清理操作"
        ;;
esac

echo ""
echo "💡 建议的下一步操作:"
echo "1. 在新的终端中运行:"
echo "   export CUDA_VISIBLE_DEVICES=''"
echo "   export TORCH_CUDA_ENABLE=0"
echo "   export FORCE_CPU=1"
echo ""
echo "2. 然后运行强制CPU版本的脚本"
echo ""
echo "3. 如果还有问题，考虑重启系统清理所有GPU进程"