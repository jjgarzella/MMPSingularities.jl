julia_exec=$(which julia)

$julia_exec startup.jl
$julia_exec --project -e 'using Pkg; Pkg.add(url="https://github.com/alexp616/CudaNTTs.jl"); Pkg.update()'
$julia_exec --project --threads 2 --gcthreads=1 experiment/run.jl &

wait