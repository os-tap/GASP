#pragma once
#include <cstdio>
#include "Params.h"
#include "Particle.h"
#include "Segments.h"
#include "Frontline.h"
#include "Screen.h"

class GAS {
public:
    GAS();

private:
    ps::Params params;

    ps::Segments main_swarm{ params };
    ps::Frontline front_line{ params };
    ps::Screen screen{ params };

    unsigned ii = 0;
    unsigned swarm_step = 0;
    void MainLoop();
    void Iteration();
    void ReadSDLEvents();
    void IterateSwarm();
    void BuildFrontline();
    void DrawScreen();
    void PrintFiles();
    void CalcFPS();


    int print_step_counter = 0;
    void CountFiles();
    void ClearFiles();

    Uint32 startTime = 0;
    Uint32 endTime = 0;
    Uint32 delta = 0;
    short fps, max_fps = 60, fps_sum = 0;
    short timePerFrame = 1000 / max_fps; // milliseconds
    char buffer[8] = "FPS: 00";


    int mouse_x, mouse_y;
    double burn_x, burn_y;

    std::string burn_radius_str{ "t, burn_radius" };
    void PrintBurnRadius();


    bool key_pressed = false;

    struct Input {
        bool sdl_quit = 0,
            set_burn = 0,
            lights_out = 0,
            clear = 0,
            pause = 0,
            step = 0,
            toggle_fill = 0,
            print_step = 0,
            print_denisty = 0,
            clear_csv = 0,
            refill = 0,
            update_curve = 0;
    } input;

    struct State {
        bool move = true,
            burn = true,
            pause = false,
            display_swarm = true,
            display_line = false,
            bold_points = false,
            blur = false,
            update_line = true,
            test = false,
            svm = false,
            printer = false;
    } state;
};