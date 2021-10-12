#include "GAS.h"

GAS::GAS()
{
    CountFiles();
    startTime = SDL_GetTicks();
    std::cout << params.burn_radius / params.burn_radius_cross;
    params.Print();
    MainLoop();
}

void GAS::MainLoop()
{
    for (; !input.sdl_quit;)
    {
        Iteration();
        ++ii;
        CalcFPS();
    }
}

void GAS::Iteration()
{
    ReadSDLEvents();

        
    IterateSwarm();


    if (state.display_line && state.update_line && (!state.pause || input.step || key_pressed))
    {
        BuildFrontline();
    }


    if (!state.pause || input.step || key_pressed)
    {
        DrawScreen();
    }
        

    if (input.print_step || (state.printer && !state.pause))
    {
        PrintFiles();
    }


    if (input.print_denisty)
    {
        main_swarm.Density_Grid();
        main_swarm.Density_Radius();
        std::cout << "\nPrint Denisty";
    }

    if (input.clear_csv)
    {
        ClearFiles();
    }
        
}

void GAS::ReadSDLEvents()
{
    SDL_Event e;
    input = {};
    key_pressed = false;
    while (SDL_PollEvent(&e))
    {
        switch (e.type)
        {
        case SDL_QUIT: input.sdl_quit = true;
            break;

        case SDL_KEYDOWN:
            switch (e.key.keysym.sym)
            {
            case SDLK_1: params.SetStream(1); break;
            case SDLK_2: params.SetStream(2); break;
            case SDLK_3: params.SetStream(3); break;
            case SDLK_4: params.SetStream(4); break;
            case SDLK_SPACE: input.lights_out = true; key_pressed = 1; break;
            case SDLK_DELETE: input.clear = true; key_pressed = 1; break;
            case SDLK_p: state.pause = !state.pause;
                std::cout << (state.pause ? "\n- Pause\n-" : "\n- Play\n-");
                break;
            case SDLK_s: input.step = true; key_pressed = 1; break;
            case SDLK_q: state.bold_points = !state.bold_points; key_pressed = 1; break;
            case SDLK_b: state.blur = !state.blur; key_pressed = 1; break;
            case SDLK_v: state.svm = !state.svm; key_pressed = 1; break;
            case SDLK_y: std::cout << ((state.printer = !state.printer) ? "\nPrinter On" : "\nPrinter Off"); key_pressed = 1; break;

            case SDLK_t: state.test = !state.test; break;

            case SDLK_h: std::cout << ((state.burn = !state.burn) ? "\n- Burn On" : "\n- Burn Off"); break;
            case SDLK_m: state.move = !state.move;
                std::cout << (state.move ? "\n- Move" : "\n- Freeze");
                break;
            case SDLK_w: state.display_swarm = !state.display_swarm; key_pressed = 1; break;
            case SDLK_e: state.display_line = !state.display_line; key_pressed = 1; params.clear_curve_spline();  break;
            case SDLK_l: state.update_line = !state.update_line; key_pressed = 1; break;
            case SDLK_f: main_swarm.Toggle_Fill();
                std::cout << "\n- Toggle Fill";
                break;
            case SDLK_c: input.print_step = true;
                break;
            case SDLK_d: input.print_denisty = true;
                break;
            case SDLK_x: input.clear_csv = true;
                break;
            case SDLK_r: params.refill = !params.refill;
                break;
            case SDLK_k: params.scale_burn = !params.scale_burn;
                params.curve_burn_coef = 0;
                break;
            case SDLK_i: input.update_curve = true;
                break;
            case SDLK_u:
                key_pressed = 1;
                params.Read();
                params.Print();
                main_swarm.Update();
                front_line.Init();
                screen.calc_refract_points();
                std::cout << main_swarm.all_list.size() << " - size\n";
                std::cout << front_line.spline_points.size() << " - error\n";
                break;
            }
            break;

            break;
        }
        // break SDL_PollEvent
        if (input.sdl_quit) return;
    }

    if (SDL_BUTTON(SDL_BUTTON_LEFT) && SDL_GetMouseState(&mouse_x, &mouse_y)) {
        key_pressed = true;
        input.set_burn = true;
        burn_x = params.screen_to_area_x(mouse_x);
        burn_y = params.screen_to_area_y(mouse_y);
    }
}

void GAS::IterateSwarm()
{
    if (!state.pause || input.step) {
        if (state.move) {
            //if (State.pause) std::cout << "\nMove";
            if (!params.refill) main_swarm.Emit();
            main_swarm.MoveParticles();
        }
    }

    if (!state.pause || input.step) {
        main_swarm.ClearParticles();
        main_swarm.UpdateSegments();
        if (params.refill) main_swarm.Refill();
        main_swarm.CalcBurnRadius(params.grid_curve_calc);
        main_swarm.PlaceBurned();
    }



    if (input.set_burn) main_swarm.BurnSegmentByPoint(burn_x, burn_y);
    if (input.lights_out) main_swarm.LightsOut();
    if (input.clear) main_swarm.EraseParticles();


    if (!state.pause || input.step) {
        //if (State.pause) std::cout << "\nCross";
        if (state.burn)
        {
            main_swarm.CrossParticles();
            main_swarm.StepParticles();

        }        //main_swarm.RefractParticles();
        //main_swarm.FinalLoop();

    }
    

}

void GAS::BuildFrontline()
{
    front_line.Calc(main_swarm.all_will_burn);
    

    if (params.calc_cross)
    {
        //std::cout << "cross";
        main_swarm.UpdateSegments();
        main_swarm.CalcFrontlineRadius(front_line.spline_points);

        //front_line.SetCrosses(main_swarm.front_crosses);
        front_line.BuildCrossesSpline(main_swarm.front_crosses);
    }
    else {
        front_line.BuildCurvatureSpline();
    }

    if (params.scale_burn)
    {
        //if (input.update_curve)
        front_line.BuildCurvatureSpline();
        params.set_curve_spline(front_line.curve_spline, front_line.curve_start, front_line.curve_end);
    }
    //front_line.Finalize();
}

void GAS::DrawScreen()
{

    screen.clear();


    if (state.display_swarm) screen.load_swarm(main_swarm.all_list, state.bold_points);
    if (state.blur) screen.box_blur();
    screen.UpdateTexture();

    //screen.draw_grid(main_swarm.grid_count_x, main_swarm.grid_count_z);


    if (state.display_line) screen.draw_frontline(front_line.spline_points);
    //if (params.display_kinks) screen.draw_hlines(front_line.kinks);


    screen.Render();
}

void GAS::PrintBurnRadius() {

    burn_radius_str += fmt::format("\n{},{}", print_step_counter, params.get_burn_radius(0));

    std::ofstream csv(params.csv_folder + "burn_radius.csv." + std::to_string(print_step_counter));
    csv << burn_radius_str;
    csv.close();
}


void GAS::PrintFiles()
{

    main_swarm.PrintCurvature(print_step_counter);
    main_swarm.PrintCount(print_step_counter, params.print_count);
    front_line.Print(print_step_counter);
    PrintBurnRadius();
    std::cout << "\nprint - " << print_step_counter;
    ++print_step_counter;
}

void GAS::CalcFPS()
{
    endTime = SDL_GetTicks();
    delta = endTime - startTime;

    startTime = endTime;

    if (delta < timePerFrame) {
        SDL_Delay(timePerFrame - delta);
        fps = max_fps;
    }
    else fps = 1000 / delta;

    fps_sum += fps;

    if (ii % 10 == 0)
    {
        sprintf_s(buffer, "FPS: %d", fps);
        screen.SetTitle(buffer);
        fps_sum = 0;
    }
}

void GAS::ClearFiles()
{
    bool gas_out = 0, line_out = 0;

    for (int i = 0;; ++i)
    {
        if (!gas_out && std::remove((params.csv_folder + "gas.csv." + std::to_string(i)).c_str()) != 0) {
            gas_out = 1;
        }
        if (!line_out && std::remove((params.csv_folder + "line.csv." + std::to_string(i)).c_str()) != 0) {
            line_out = 1;
        }

        if (gas_out && line_out) break;

    }

    print_step_counter = 0;
    burn_radius_str = "t, burn_radius";
    burn_radius_str += fmt::format("\n{},{}", -1, params.get_burn_radius(0));
    std::cout << "\nclear files";
}

void GAS::CountFiles()
{
    std::ifstream check_file;
    for (;;) {
        check_file.open((params.csv_folder + "line.csv." + std::to_string(print_step_counter)).c_str());
        if (check_file.fail()) break;
        check_file.close();
        ++print_step_counter;
    }
}