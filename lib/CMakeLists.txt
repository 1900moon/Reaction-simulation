cmake_minimum_required(VERSION 3.20)
project(lise_lib
    VERSION 1.0.0
    DESCRIPTION "okawa practice"
    # (reference "https://kamino.hatenablog.com/entry/cmake_tutorial2")
    # HOMEPAGE_URL "https://example.com"
    LANGUAGES CXX)

add_library(refer STATIC 
    ./source/function.cpp
    ./source/beam.cpp)

target_link_libraries(refer ROOT::Core ROOT::Gpad ROOT::Graf3d ROOT::Graf
  ROOT::Hist ROOT::Imt ROOT::MathCore ROOT::Matrix ROOT::MultiProc ROOT::Net
  ROOT::Physics ROOT::Postscript ROOT::RIO ROOT::ROOTDataFrame
  ROOT::ROOTVecOps ROOT::Rint ROOT::Thread ROOT::TreePlayer ROOT::Tree)
target_compile_features(refer PUBLIC cxx_std_11)
target_include_directories(refer INTERFACE ./include)
set_target_properties(refer
    PROPERTIES
    VERSION ${PROJECT_VERSION})





