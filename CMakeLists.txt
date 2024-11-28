cmake_minimum_required(VERSION 3.29)
project(CGALStudy)

set(CMAKE_CXX_STANDARD 17)
set(CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE TRUE)

find_package(CGAL CONFIG REQUIRED COMPONENTS Qt6)
#find_package(CGAL CONFIG REQUIRED)
add_executable(CGALStudy
    main.cpp
    Kriging.cpp
    Kriging.h
    UsingDeclaration.h
    MyBbox.h
)

if(CGAL_Qt6_FOUND)
    target_link_libraries(CGALStudy PUBLIC CGAL::CGAL_Basic_viewer)
else()
    message(STATUS "NOTICE: The example 'draw_triangulation_3' requires Qt6, and will not be compiled.")
endif()
target_link_libraries(CGALStudy PRIVATE CGAL::CGAL)

include_directories("D:/geoLIbs/eigen-3.4.0")