cmake -S ../../ -B ../../build -G Ninja -DCOMPILER_OPTION=clang -DFLOAT_PREC=2 -DPREALLOC_SIZE="1<<15" -DARRAY_DFT_METHOD=ARRAY_DFT_ASYNC -DNUM_DFT_THREADS=8 -DFFT_SPEED=ON -DFFT_ACCURACY=ON
cmake --build ../../build
read -p "Done!"
