all: BeQ2L

RMSE: computeRMSE.cpp matmethods.hpp h5objcpp.hpp
	g++ -O0 -std=c++11 -g computeRMSE.cpp   -o RMSE  -I/usr/local/boost_1_47_0 -I/usr/include -I/opt/intel/composerxe-2011.5.220/mkl/include/ -L/opt/intel/composerxe-2011.5.220/mkl/lib/intel64/ -L/usr/local/lib -I/home/nwk2/include -L/home/nwk2/lib -L/home/nwk2/usr/lib64 -I/home/nwk2/usr/include -larmadillo -lhdf5_cpp -lhdf5_hl_cpp -lssl -lcrypto -lhdf5 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm

BeQ2L: main.cpp matmethods.hpp parameter_read.hpp importdata.hpp h5objcpp.hpp
	icpc -std=c++11 -ipo -prof-gen -prof-dir=/home/nwk2/ -O3 -g main.cpp -o BeQ2L -I/usr/local/include -I/usr/local/boost_1_47_0 -I/usr/local/boost_1_47_0/boost -L/usr/local/lib  -larmadillo -lhdf5_cpp -lhdf5_hl_cpp -lssl -lcrypto -lhdf5 -lboost_timer

cubetest: cubetest.cpp 
	g++ -std=c++11 -O3 -g cubetest.cpp -o testcube -I/usr/local/include -I/usr/local/boost_1_47_0 -L/usr/local/lib -larmadillo

findtest: findtest.cpp
	g++ -std=c++11 -O3 -g findtest.cpp -o findtest -I/usr/local/include -L/usr/local/lib -larmadillo

test_import: import_data.cpp h5objcpp.hpp
	g++ -std=c++11 -O3 -g import_data.cpp -o test_import -I/usr/local/include -I/usr/local/boost_1_47_0 -L/usr/local/lib -I/home/nwk2/include -L/home/nwk2/lib -L/home/nwk2/usr/lib64 -I/home/nwk2/usr/include -larmadillo -lhdf5_cpp -lhdf5_hl_cpp -lssl -lcrypto -lhdf5

demomat: demomat.cpp h5objcpp.hpp
	g++ -std=c++11 -O3 -g demomat.cpp -o demomat -I/usr/local/include -I/usr/local/boost_1_47_0 -L/usr/local/lib -I/home/nwk2/include -L/home/nwk2/lib -L/home/nwk2/usr/lib64 -I/home/nwk2/usr/include -larmadillo -lhdf5_cpp -lhdf5_hl_cpp -lssl -lcrypto -lhdf5

divtest: divtest.cpp
	g++ -std=c++11 -O3 -g divtest.cpp -o divtest -I/usr/local/include -I/home/nwk2/include -L/home/nwk2/lib -L/home/nwk2/usr/lib64 -I/home/nwk2/usr/include -larmadillo

chunktest: chunktest.cpp h5objcpp.hpp
	g++ -std=c++11 -O0 -g chunktest.cpp -o chunktest -I/usr/local/include -I/usr/local/boost_1_47_0 -L/usr/local/lib -I/home/nwk2/include -L/home/nwk2/lib -L/home/nwk2/usr/lib64 -I/home/nwk2/usr/include -larmadillo -lhdf5_cpp -lhdf5_hl_cpp -lssl -lcrypto -lhdf5

RMSEtest: RMSEtest.cpp matmethods.hpp
	g++ -std=c++11 -O3 -g RMSEtest.cpp -o RMSEtest -I/usr/local/include -I/usr/local/boost_1_47_0 -L/usr/local/lib -I/home/nwk2/include -L/home/nwk2/lib -L/home/nwk2/usr/lib64 -I/home/nwk2/usr/include -larmadillo -lhdf5_cpp -lhdf5_hl_cpp -lssl -lcrypto -lhdf5

Gen_Submat: generate_submat.cpp matmethods.hpp h5objcpp.hpp
	/opt/gcc-4.9.1/bin/g++ -std=c++11 -O3 -g generate_submat.cpp -o Gen_Submat -I/usr/local/include -I/usr/local/boost_1_47_0 -L/usr/local/lib -I/home/nwk2/include -L/home/nwk2/lib -L/home/nwk2/usr/lib64 -I/home/nwk2/usr/include -larmadillo -lhdf5_cpp -lhdf5_hl_cpp -lssl -lcrypto -lhdf5

mkstreamtest: mkstreamtest.cpp
	icc  -O3 -g mkstreamtest.cpp -o mstreamtest  -L/usr/local/lib  -L/home/nwk2/usr/lib64 -L/opt/intel/composerxe/mkl/lib/intel64/ -I/opt/intel/composerexe-2011.5.220/mkl/include  -L/opt/intel/composerexe-2011.5.220/mkl/lib/intel64/ -I/home/nwk2/usr/include -larmadillo -mkl


clean:
	rm -f BeQ2L
