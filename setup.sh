conan profile detect --force
if [ -d build ]; then
    rm -r build
fi
conan install . --output-folder=. --build=missing --settings=build_type=Release

