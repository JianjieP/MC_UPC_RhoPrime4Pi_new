cmake -S . -B build
cmake --build build --target generate_grid_opt_job merge_grid_opt_jobs -j4

for j in $(seq 0 19); do
  ./build/generate_grid_opt_job PbPb 5360 NoTag build/grid_PbPb5360_Notag_smoke.root 40 100 1 4 1 20 400 1 "$j" 20 > "build/split20_job${j}.log" 2>&1 &
done
wait

./build/merge_grid_opt_jobs build/grid_PbPb5360_Notag_smoke.root 20