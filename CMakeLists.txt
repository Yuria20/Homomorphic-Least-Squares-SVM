cmake_minimum_required(VERSION 3.8.0)
project(homsvp VERSION 0.1.0 LANGUAGES C CXX)

include_directories(/usr/local/include)
link_directories(/usr/local/lib)


include(CTest)
enable_testing()


#add_executable(homsvp main.cpp)
#add_executable(homsvp main.cpp poly.h poly.cpp CustomEvaluator.cpp CustomEvaluator.h lssvp.h lssvp.cpp)
add_executable(homsvp test2.cpp poly.h poly.cpp custom_eval.h custom_eval.cpp HEmatrix.h lssvm.h hom_lssvm.h)

target_link_libraries(homsvp PRIVATE /usr/local/lib/libHEaaN.so /usr/local/lib/libntl.a /usr/local/lib/libgmp.a /usr/local/lib/libgmp.so)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
set(CMAKE_CXX_STANDARD 11)
include(CPack)
