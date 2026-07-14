cmake_minimum_required(VERSION 3.16)
project(RhoPrimeGenerator LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(ROOT REQUIRED COMPONENTS Core Hist MathCore RIO Tree Physics)
find_package(OpenMP)

add_library(rhoprime_core
  src/EventGeneratorBose.cxx
  src/UpcRhoPrimeModel.cxx
)

target_include_directories(rhoprime_core PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${ROOT_INCLUDE_DIRS}
)

target_link_libraries(rhoprime_core PUBLIC
  ROOT::Core
  ROOT::Hist
  ROOT::MathCore
  ROOT::RIO
  ROOT::Tree
  ROOT::Physics
)

if(OpenMP_CXX_FOUND)
  set_property(SOURCE src/EventGeneratorBose.cxx APPEND PROPERTY
    COMPILE_OPTIONS ${OpenMP_CXX_FLAGS}
  )
  target_link_libraries(rhoprime_core PUBLIC $<LINK_ONLY:OpenMP::OpenMP_CXX>)
endif()

add_executable(generate_bose
  src/main_generate_bose.cxx
)

target_link_libraries(generate_bose PRIVATE rhoprime_core)

add_executable(generate_grid
  src/main_rhoprime_dsigma_dy.cxx
)

target_link_libraries(generate_grid PRIVATE rhoprime_core)

add_executable(generate_grid_opt
  src/main_generate_grid_opt.cxx
  src/UpcRhoPrimeModelOpt.cxx
)

target_link_libraries(generate_grid_opt PRIVATE rhoprime_core)

add_executable(generate_grid_opt_job
  src/main_generate_grid_opt_job.cxx
  src/UpcRhoPrimeModelOpt.cxx
)

target_link_libraries(generate_grid_opt_job PRIVATE rhoprime_core)

add_executable(merge_grid_opt_jobs
  src/main_merge_grid_opt_jobs.cxx
)

target_link_libraries(merge_grid_opt_jobs PRIVATE ROOT::Core ROOT::Hist ROOT::MathCore ROOT::RIO)

add_executable(analyze_bose
  analyze_bose.cxx
)

target_link_libraries(analyze_bose PRIVATE ROOT::Core ROOT::Hist ROOT::MathCore ROOT::RIO ROOT::Tree ROOT::Physics)

add_executable(analyze_nonbose
  analyze_nonbose.cxx
)

target_link_libraries(analyze_nonbose PRIVATE ROOT::Core ROOT::Hist ROOT::MathCore ROOT::RIO ROOT::Tree ROOT::Physics)

add_executable(compare_bose_modulations
  compare_bose_modulations.cxx
)

target_link_libraries(compare_bose_modulations PRIVATE ROOT::Core ROOT::MathCore ROOT::RIO ROOT::Tree ROOT::Physics)
