cd build/Release/generators
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE=conan_toolchain.cmake ../../..
cmake --build . --config Release 
./main