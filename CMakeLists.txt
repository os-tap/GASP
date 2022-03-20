cmake_minimum_required(VERSION 3.0.0)
project(GASP VERSION 0.1.0)

set (CMAKE_CXX_STANDARD 17)


# include(CTest)
# enable_testing()


add_subdirectory(library)

add_executable(GASP src/main.cpp
    src/GAS.cpp
    src/Params.cpp
    src/Field.cpp
    src/Particle.cpp
    src/Segments.cpp
    # src/Frontline.cpp
    src/Screen.cpp
)

# set(CPACK_PROJECT_NAME ${PROJECT_NAME})
# set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
# include(CPack)


target_link_libraries(GASP
  PRIVATE  
    fmt::fmt 
    nlohmann_json::nlohmann_json
    SDL2main SDL2-static m
    tbb
    gsl
)