julia_exec=$(which julia)

$julia_exec startup.jl
$julia_exec --project -e 'using Pkg; Pkg.update()'
$julia_exec --project --threads 2 --gcthreads=1 -e "include(\"experiment/run.jl\")" &

wait
