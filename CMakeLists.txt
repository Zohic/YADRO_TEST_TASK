cmake_minimum_required(VERSION 3.18)
include(CMakePrintHelpers)

set(CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS}")

IF(DEFINED COMPILER_OPTION)
	IF(${COMPILER_OPTION} STREQUAL "clang")
		set(CMAKE_CXX_COMPILER "clang++.exe")
		set(CMAKE_C_COMPILER "clang.exe")
		message("compiler is set to clang, make sure c/с++ compilers binaries are specified in PATH")
	ELSEIF(${COMPILER_OPTION} STREQUAL "gcc")
		set(CMAKE_CXX_COMPILER "g++.exe")
		set(CMAKE_C_COMPILER "gcc.exe")
		message("compiler is set to gcc, make sure c/с++ compilers binaries are specified in PATH")
	ELSE()
		message("WARNING: SPECIFIED UNKNOW COMPILER")
		message("Compiler will be chosen automatically. On Windows it might be MSVC which can only work from powershell terminal")
	ENDIF()
ELSE()
	message("WARNING: Compiler will be chosen automatically (you may set \"COMPILER_OPTION\" to \"clang\" or \"gcc\"). On Windows it might be MSVC which can only work from powershell terminal")
ENDIF()
	

IF(DEFINED FLOAT_PREC)
	add_compile_definitions("prec=${FLOAT_PREC}")
ELSE()
	message(FATAL_ERROR "ERROR: \"FLOAT_PREC\" IS NOT DEFINED. Use -D FLOAT_PREC=1|2")
ENDIF()

IF(DEFINED PREALLOC_SIZE)
	add_compile_definitions("PREALLOC_SIZE_DEF=${PREALLOC_SIZE}")
ELSE()
	set(PREALLOC_SIZE "1<<15")
	add_compile_definitions("PREALLOC_SIZE_DEF=${PREALLOC_SIZE}")
	message("WARNING: amount of preallocated memory isn't specified (use \"PREALLOC_SIZE\"). Set to default 2^15")
ENDIF()

IF(DEFINED ARRAY_DFT_METHOD)
	IF(${ARRAY_DFT_METHOD} STREQUAL "ARRAY_DFT_ASYNC" OR ${ARRAY_DFT_METHOD} STREQUAL "ARRAY_DFT_THREADS")
		add_compile_definitions("ARRAY_DFT_METHOD=${ARRAY_DFT_METHOD}")
	ELSE()
		message(FATAL_ERROR "UNKNOW parameter for \"ARRAY_DFT_METHOD\" use ARRAY_DFT_ASYNC or ARRAY_DFT_THREADS (for ARRAY_DFT_THREADS you must also specify NUM_DFT_THREADS)")
	ENDIF()
ELSE()
	message(FATAL_ERROR "ERROR: \"ARRAY_DFT_METHOD\" IS NOT DEFINED. use ARRAY_DFT_ASYNC or ARRAY_DFT_THREADS (for ARRAY_DFT_THREADS you must also specify NUM_DFT_THREADS)")
ENDIF()

IF(DEFINED NUM_DFT_THREADS)
	add_compile_definitions("NUM_DFT_THREADS_DEF=${NUM_DFT_THREADS}")
ELSE()
	set(NUM_DFT_THREADS 4)
	add_compile_definitions("NUM_DFT_THREADS_DEF=${NUM_DFT_THREADS}")
	message("WARNING: number of threads isn't specified (use \"NUM_DFT_THREADS\"). Set to default 4")
ENDIF()

set(CMAKE_BUILD_TYPE Release)

project(YADRO_FFT
  VERSION 1.0
  DESCRIPTION "test task for YADRO"     
  HOMEPAGE_URL "https://github.com/Zohic/YADRO_TEST_TASK" 
  LANGUAGES C CXX                      
)

option(FFT_SPEED "build speed test" OFF)
option(FFT_ACCURACY "build accuracy test" OFF)

IF(NOT FFT_ACCURACY AND NOT FFT_SPEED)
	message(FATAL_ERROR "at least one of the options \"FFT_SPEED\" or \"FFT_ACCURACY\" must be ON in order to build")
ENDIF()

if(FFT_ACCURACY)
    add_subdirectory("FFT/test/accuracy")
endif()

if(FFT_SPEED)
    add_subdirectory("FFT/test/speed")
endif()

unset(COMPILER_OPTION CACHE)

unset(FFT_SPEED CACHE)
unset(FFT_ACCURACY CACHE)

unset(ARRAY_DFT_METHOD CACHE)
unset(NUM_DFT_THREADS CACHE)

unset(PREALLOC_SIZE CACHE)

unset(FLOAT_PREC CACHE)
