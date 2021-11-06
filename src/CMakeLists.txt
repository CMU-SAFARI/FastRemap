if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_RELEASE "-O3")

cmake_minimum_required (VERSION 3.4)
project (seqan3_tutorial CXX)
set(CMAKE_CXX_STANDARD 17)  

#set(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} "/usr/local/opt/boost\@1.76") 
#SET(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} "/usr/local/opt/boost\@1.76/lib") 

#set(BOOST_ROOT "/usr/local/opt/boost\@1.76") 
#set(_boost_INCLUDE_SEARCH_DIRS ${_boost_INCLUDE_SEARCH_DIRS} "/usr/local/opt/boost\@1.76") 
#set(Boost_DEBUG 1) 
#set(Boost_USE_STATIC_LIBS OFF) 
#set(Boost_USE_MULTITHREADED ON)  
#set(Boost_USE_STATIC_RUNTIME OFF) 
#find_package(Boost COMPONENTS) 


# add seqan3 to search path
list(APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../seqan2")

# require seqan3 with a version between >=3.0.0 and <4.0.0
find_package (seqan2 2.0 REQUIRED)
 
# build app with seqan3
#include_directories(${Boost_INCLUDE_DIR}) 
add_executable (FastRemap main.cpp utils.cpp mapbam.cpp) 
target_link_libraries (FastRemap seqan2::seqan2)  # ${Boost_LIBRARIES}) 

