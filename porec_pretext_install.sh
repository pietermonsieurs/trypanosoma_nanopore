module purge
module load calcua/2024a
module load Meson
module load Clang

mkdir cd data/software/PretextMap/
cd cd data/software/PretextMap/
env CXX=clang meson setup --buildtype=release --unity on --prefix=/user/antwerpen/205/vsc20587/data/software/PretextMap  builddir


cd builddir/
meson compile
meson test
meson install

