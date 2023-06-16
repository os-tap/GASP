include(FetchContent)

# --------------------------------------------------------------------

# Offline mode (uncomment next line to enable)
# set(FETCHCONTENT_FULLY_DISCONNECTED ON)

# set(FETCHCONTENT_QUIET OFF)

# --------------------------------------------------------------------

# Libraries

# --------------------------------------------------------------------

message(STATUS "FetchContent: SDL")

FetchContent_Declare(
        SDL
        GIT_REPOSITORY https://github.com/libsdl-org/SDL.git
        GIT_TAG 739f783
)
FetchContent_MakeAvailable(SDL)

# FetchContent_GetProperties(SDL2)
# if (NOT sdl2_POPULATED)
#     FetchContent_Populate(SDL2)
#     add_subdirectory(${sdl2_SOURCE_DIR} ${sdl2_BINARY_DIR})
# endif ()

# --------------------------------------------------------------------

message(STATUS "FetchContent: fmtlib")

FetchContent_Declare(
        fmtlib
        GIT_REPOSITORY https://github.com/fmtlib/fmt.git
        GIT_TAG 8.0.1
)
FetchContent_MakeAvailable(fmtlib)

# --------------------------------------------------------------------

message(STATUS "FetchContent: json")

FetchContent_Declare(json
  GIT_REPOSITORY https://github.com/ArthurSonzogni/nlohmann_json_cmake_fetchcontent
  GIT_TAG v3.10.3)

FetchContent_GetProperties(json)
if(NOT json_POPULATED)
  FetchContent_Populate(json)
  add_subdirectory(${json_SOURCE_DIR} ${json_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()

# --------------------------------------------------------------------

# message(STATUS "FetchContent: splinter")

# FetchContent_Declare(
#         splinter
#         GIT_REPOSITORY https://github.com/bgrimstad/splinter.git
#         GIT_TAG v3.0
# )
# FetchContent_MakeAvailable(splinter)

# --------------------------------------------------------------------

# message(STATUS "FetchContent: gsl")

# FetchContent_Declare(
#         gsl
#         GIT_REPOSITORY https://github.com/ampl/gsl.git
#         GIT_TAG v2.7.0
# )
# FetchContent_MakeAvailable(gsl)

# --------------------------------------------------------------------

message(STATUS "FetchContent: one tbb")

FetchContent_Declare(
        oneTBB
        GIT_REPOSITORY https://github.com/oneapi-src/oneTBB.git
        GIT_TAG v2021.9.0
)
FetchContent_MakeAvailable(oneTBB)

# --------------------------------------------------------------------
