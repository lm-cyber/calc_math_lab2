



do:
	g++ -O3 main.cpp util.cpp secant_method.cpp fixed_point_iteration_metod_sys.cpp  fixed_point_iteration_metod.cpp -o prog -lblas -llapack -DHAVE_CBLAS=1
