#include "Params.h"



ps::Params::Params() {
    Read();
}

void ps::Params::Read()
{
    Load(GetFromFile());
    frontline_kinks.clear();
}


json ps::Params::GetFromFile() {
    std::ifstream file("params.json");
    assert(file);

    nlohmann::json j;
    file >> j;
    file.close();

    return j;
}



void ps::Params::Load(json j)
{

    csv_folder = j["csv_folder"];

    area_beg = j["area_beg"];
    area_end = j["area_end"];
    area_size = area_end - area_beg;
    area_height = area_size / screen_width * screen_height;

    curve_start = area_end; curve_end = area_beg;

    screen_proportion = screen_width / area_size;
    screen_x_proportion = area_size / screen_width;
    screen_y_proportion = area_height / screen_height;

    stream_beg = j["stream_beg"];
    stream_end = j["stream_end"];
    stream_width = stream_end - stream_beg;
    stream_radius = stream_width / 2;

    area_center = stream_beg + stream_radius;
    L = stream_width;

    scale = j["scale"];
    DSR = L / scale;



    burn_radius = (double)j["burn_radius"] * DSR;
    burn_radius_2 = burn_radius * burn_radius;

    base_particles = (int)j["base_particles"];
    particles_dist = burn_radius / (int)j["particles_dist"];

    //burn_radius_cross = burn_radius + burn_radius * (double)j["burn_fix"];
    //burn_radius_cross = burn_radius * (1 + (pow(base_particles, (double)j["burn_fix"])));
    burn_fix_a = (double)j["burn_fix_a"];
    burn_fix_b = (double)j["burn_fix_b"];
    burn_radius_cross = make_radius_cross_fix(burn_radius);
    burn_radius_2_cross = burn_radius_cross * burn_radius_cross;

    emitter_begin = (double)j["emitter_begin"] * burn_radius;


    base_speed = (double)j["base_speed"] * burn_radius;

    burn_speed = (double)j["burn_speed"] / burn_radius;

    const_speed = (double)j["const_speed"] * burn_radius;

    iterations = (int)j["iterations"];

    iterate_speed = base_speed / iterations;
    iterate_const = const_speed / iterations;
    iterate_particles = (int)floor(base_particles * particle_speed(area_center) / burn_radius_2_cross / M_PI * L);

    burn_time = (int)j["burn_time"] * iterations;
    sage_time = (int)j["sage_time"] * iterations;
    wave_time = (int)j["wave_time"] * iterations;

    frontline_spline_steps = (int)j["frontline_spline_steps"];
    frontline_spline_alpha = (double)j["frontline_spline_alpha"];
    move_normal = (bool)j["move_normal"];
    move_speed = (bool)j["move_speed"];
    curve_by_radius = (bool)j["curve_by_radius"];


    frontline_window_steps = (int)j["frontline_window_steps"];
    frontline_window_size = (double)j["frontline_window_size"];

    frontline_stencil_h = (int)j["frontline_stencil_h"];
    frontline_radius_h = (int)j["frontline_radius_h"];

    refract_coef = (double)j["refract_coef"];
    refract_offset = (double)j["refract_offset"];

    svm_count = (int)j["svm_count"];
    print_count = (int)j["print_count"];
    display_count = (int)j["display_count"];

    calc_cross = (bool)j["calc_cross"];
    frontline_cross_chunk = (bool)j["frontline_cross_chunk"];
    frontline_cross_multipler = (double)j["frontline_cross_multipler"];
    frontline_cross_area = (double)j["frontline_cross_area"];
    frontline_cross_radius = (stream_width * frontline_cross_area);
    frontline_cross_radius_2 = frontline_cross_radius * frontline_cross_radius;

    frontline_cross_border = (double)j["frontline_cross_border"];
    frontline_cross_spline_alpha = (double)j["frontline_cross_spline_alpha"];

    
    who_cross = (bool)j["who_cross"];
    scale_burn = (bool)j["scale_burn"];
    scale_burn_multipler = (double)j["scale_burn_multipler"];
    scale_burn_size = stream_width * (double)j["scale_burn_size"];
    scale_burn_pow = (double)j["scale_burn_pow"];
    scale_burn_condition = (double)j["scale_burn_condition"];

    curve_burn_coef = (double)j["curve_burn_coef"];

}


void ps::Params::Print()
{
    std::cout << "\n-------------\n";
    std::cout << stream_function(area_center) * base_speed << " - base speed\n";
    std::cout << particle_speed(area_center) << " - max delta\n";
    std::cout << burn_radius << " - burn radius\n";
    std::cout << burn_radius_2 << " - burn radius 2\n";
    std::cout << burn_radius_2_cross << " - burn radius 2 fix\n";
    std::cout << burn_speed << " - burn speed\n";
    std::cout << area_height << " - area_height\n";
}
/*

std::string ps::Params::frontline_params() const
{
	return fmt::format("#steps {}\n#window 1/{}\n#h {}\n", front_line_steps, front_line_windows, front_line_h);
}
*/


