mkdir build
cp *.h build/
cat groupie.f endfio.f scratchb.f timer.f seconds.f > build/x.f
cd build
gfortran -o groupie x.f
