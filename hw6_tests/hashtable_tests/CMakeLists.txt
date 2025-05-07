include_directories(. ../hashtable_tests)

add_header_problem(
	NAME hashtable 
	TEST_SOURCE 
		hashtable_tests.cpp)


#		RUNTIME_TEST_SOURCE
#		hashtable_stress_tests.cpp)
