project(GASP VERSION 0.1.0)
set (CMAKE_CXX_STANDARD 20)

add_executable(GASP main.cpp
    GAS.cpp
    Params.cpp
    Field.cpp
    Particle.cpp
    Segments.cpp
    Screen.cpp
    # Frontline.cpp
)

target_link_libraries(GASP
  # PRIVATE  
    fmt::fmt
    nlohmann_json::nlohmann_json
    SDL3::SDL3
    # gsl
    # tbb
    TBB::tbb
    # Threads::Threads
)