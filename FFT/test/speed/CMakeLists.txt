add_executable("${PROJECT_NAME}_speed_test")

target_sources("${PROJECT_NAME}_speed_test" PRIVATE 
	"${CMAKE_CURRENT_LIST_DIR}/test_speed.cpp" 
	"${CMAKE_SOURCE_DIR}/FFT/src/matrix_vector_mapper.cpp")

target_include_directories("${PROJECT_NAME}_speed_test" PRIVATE "${CMAKE_SOURCE_DIR}/FFT/include")