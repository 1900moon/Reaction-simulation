cmake_minimum_required(VERSION 3.20)
project(lise_practice CXX)

find_package(ROOT 6.22 CONFIG REQUIRED)
if(ROOT_FOUND)
  message(STATUS "ROOT_INCLUDE_DIRS: ${ROOT_INCLUDE_DIRS}")
endif()
include("${ROOT_USE_FILE}")

add_subdirectory(./lib)
add_executable(a.out main.cpp)

target_link_libraries(a.out ROOT::Core ROOT::Gpad ROOT::Graf3d ROOT::Graf
  ROOT::Hist ROOT::Imt ROOT::MathCore ROOT::Matrix ROOT::MultiProc ROOT::Net
  ROOT::Physics ROOT::Postscript ROOT::RIO ROOT::ROOTDataFrame
  ROOT::ROOTVecOps ROOT::Rint ROOT::Thread ROOT::TreePlayer ROOT::Tree)


target_link_libraries(a.out refer)
