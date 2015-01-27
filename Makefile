all: BeQ2L

RMSE: computeRMSE.cpp matmethods.hpp h5objcpp.hpp
	g++ -std=c++11 -O3 -g computeRMSE.cpp -o RMSE -I/usr/local/include -I/usr/local/boost_1_47_0 -L/usr/local/lib -larmadillo -lhdf5_cpp -lhdf5_hl_cpp -lssl -lcrypto -lhdf5	

BeQ2L: main.cpp matmethods.hpp parameter_read.hpp importdata.hpp h5objcpp.hpp
	g++ -std=c++11 -O3 -g main.cpp -o BeQ2L -I/usr/local/include -I/usr/local/boost_1_47_0 -I/usr/local/boost_1_47_0/boost -L/usr/local/lib  -larmadillo -lhdf5_cpp -lhdf5_hl_cpp -lssl -lcrypto -lhdf5 -lboost_timer

cubetest: cubetest.cpp 
	g++ -std=c++11 -O3 -g cubetest.cpp -o testcube -I/usr/local/include -I/usr/local/boost_1_47_0 -L/usr/local/lib -larmadillo

findtest: findtest.cpp
	g++ -std=c++11 -O3 -g findtest.cpp -o findtest -I/usr/local/include -L/usr/local/lib -larmadillo

test_import: import_data.cpp h5objcpp.hpp
	g++ -std=g++11 -O3 -g import_data.cpp -o test_import -I/usr/local/include -I/usr/local/boost_1_47_0 -L/usr/local/lib -larmadillo -lhdf5_cpp -lhdf5_hl_cpp -lhdf5 -lssl -lcrypto



clean:
	rm -f BeQ2L
