mkdir build
cp *.f build
cd build
gfortran -o sixlin   sixlin.f
gfortran -o gamlin   gamlin.f
gfortran -o doace    doace.f
gfortran -o dotsl    dotsl.f
gfortran -o dodos    dodos.f
gfortran -o dophn    dophn.f
gfortran -o acemaker acemaker.f
