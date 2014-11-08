BeQ2L
=====

Free and Open source implementation of Bootstrap eQTL

By making use of the libre licensed OpenBLAS library and boost libraries in place of the proprietary Intel MKL.
The entire codebase is rewritten to improve clarity of the code as well as performance.

The software consists of the following components

1) import_data
   This program takes as arguments two files corresponding to the matrices A&B as well as a destination output file, which is a single HDF5 file which will contain the data for each matrix.
   Each matrix must have the same number of columns, and each row must be labeled with a row name.


Classes
------

Param
Fields
	const char* paramfile
	size_t rownum
	vector<size_t> Acols
	vector<size_t> Bcols
	bool streaming
	size_t bootnum
	int chunkstart
	int chunknum
	const char* Amatfield
	const char* Bmatfield
	vector<string> Arows
	vector<string> Brows

