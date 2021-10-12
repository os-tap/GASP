#include "Segments.h"




namespace ps {

    Segments::Segments(const Params& aParams) : P(&aParams)
    {
        Update();
    }
    void Segments::Update()
    {
        SetSegmentsGrid(P->burn_radius_cross);
        ResetFillGrid();
        refill = P->refill;
    }


    void Segments::Toggle_Fill() {
        _fill_one = !_fill_one;
    }

    void Segments::FillParticles() {
        _fill_one ? Fill_Grid() : Fill_Sampling();
    }



    void Segments::CreateParticle(double x_cord, double z_cord) {
        all_list.emplace_back(x_cord, z_cord);
    }

    void Segments::SetSegmentsGrid(double seg_size)
    {
        //double seg_size = P->burn_radius_cross;
        grid_count_x = (int)floor(P->area_size / seg_size);
        grid_count_z = (int)floor(P->area_height / seg_size);
        grid_count = grid_count_x * grid_count_z;

        grid_x_size = P->area_size / grid_count_x;
        grid_z_size = P->area_height / grid_count_z;
        grid_min_size = grid_x_size * (grid_x_size <= grid_z_size) + grid_z_size * (grid_z_size < grid_x_size);
        grid_max_size = grid_x_size * (grid_x_size >= grid_z_size) + grid_z_size * (grid_z_size > grid_x_size);

        grid_count_x_percent = grid_count_x / P->area_size;
        grid_count_z_percent = grid_count_z / P->area_height;


        grid.resize(grid_count);

        for (int i = 0; i < grid_count; ++i) {
            grid[i].x = i % grid_count_x;
            grid[i].z = i / grid_count_x;            
            
            grid[i].seg_start_x = grid[i].x * grid_x_size + P->area_beg;
            grid[i].seg_start_z = grid[i].z * grid_z_size;
        }

        grid_particles_count = (double)P->base_particles / P->burn_radius_2 / M_PI * grid_x_size * grid_z_size;
        grid_particles_min = grid_particles_count / (1 + P->grid_count_gap);
        grid_particles_max = grid_particles_count * (1 + P->grid_count_gap);


        delete[]gxa;
        delete[]gya;

        gxa = new double[grid_count_x - P->grid_curve_calc*2];
        gya = new double[grid_count_z - P->grid_curve_calc*2];

        for (size_t i = P->grid_curve_calc; i < grid_count_x - P->grid_curve_calc; i++) {
            double xai = (i + 0.5) * grid_x_size + P->area_beg;
            gxa[i - P->grid_curve_calc] = xai;
        }
        for (size_t i = P->grid_curve_calc; i < grid_count_z - P->grid_curve_calc; i++) {
            gya[i - P->grid_curve_calc] = (i + 0.5) * grid_z_size;
        }

        //spline2d.init(gxa, gya, grid_count_x - P->grid_curve_calc * 2, grid_count_z - P->grid_curve_calc * 2);

        CircleSquareSampling(P->grid_curve_area);

    }



    void Segments::ResetFillGrid()
    {
        fill_grid_count = (int)floor(P->stream_width / P->particles_dist);
        last_particles.resize(fill_grid_count);
        for (int i = 0; i < fill_grid_count; i++) {
            last_particles[i] = 0;
        }
    }



    void Segments::PrintSwarm(int num)
    {
        std::string output;
        output += "x,z,burn";
        for (auto& particle : all_list) {
            output += fmt::format("\n{},{},{}", particle._x(), particle._z(), particle.burn_counter);
        }
        std::ofstream csv(P->csv_folder + "gas.csv." + std::to_string(num));
        csv << output;
        csv.close();
    }


    void Segments::PrintCount(int num, int count)
    {

        std::string output;

        //output += "#base_particles " + std::to_string(P->base_particles);

        output += "x,z,burn";

        int step = all_list.size() / count;
        for (size_t i = 0; i < all_list.size(); i += step) {
            output += fmt::format("\n{},{},{}", all_list[i].x, all_list[i].z, all_list[i].burn_counter);
        }

        std::ofstream csv(P->csv_folder + "gas.csv." + std::to_string(num));
        csv << output;
        csv.close();

    }

    void Segments::PrintCurvature(int num)
    {

        std::string output;

        //output += "#base_particles " + std::to_string(P->base_particles);

        output += "x,z,c,c_ok,c_bu";
        int g = P->grid_curve_calc+1;

        for (const auto& seg : grid) {

            //if (seg.x >= g && seg.x < grid_count_x - g && seg.z >= g && seg.z < grid_count_z - g)
            {

                double x = seg.x / (double)grid_count_x * P->area_size + P->area_beg + grid_x_size / 2;
                double z = seg.z / (double)grid_count_z * P->area_height + grid_z_size / 2;
                output += fmt::format("\n{},{},{},{},{}", x, z, seg.curvature,seg.ok_list.size(),seg.b_list.size());
                //output += fmt::format("\n{},{},{},{},{},{}", x, z, seg.curvature,spline2d.eval(x,z),seg.ok_list.size(),seg.b_list.size());

            }
        }

        std::ofstream csv(P->csv_folder + "gas_curv.csv." + std::to_string(num));
        csv << output;
        csv.close();



    }

    void Segments::PrintSVM(int num, int count)
    {
        //all_will_burn.clear();

        std::string output;

        //output += "#base_particles " + std::to_string(P->base_particles);

        //output += "x,z,burn";
        int step = all_list.size() / count;
        for (size_t i = 0; i < all_list.size(); i+= step)
        {
            output += fmt::format("{} 1:{} 2:{}\n", 1 + 2*!all_list[i].isOk(), all_list[i].x, 1 - all_list[i].z);
        }

        std::ofstream csv(P->csv_folder + "svm.csv." + std::to_string(num));
        csv << output;
        csv.close();



    }

    void Segments::Emit()
    {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> dist;

        double p_x_cord, p_z_cord;
        double emitter_height = P->particle_speed_z(P->area_center, 0)*1.1;
        int particles_per_step = round((double)P->base_particles / P->burn_radius_2 / M_PI * P->stream_width * emitter_height);

        for (int i = 0; i < particles_per_step; i++)
        {
            p_x_cord = dist(gen) * P->stream_width + P->stream_beg;
            p_z_cord = -dist(gen) * emitter_height;
            CreateParticle(p_x_cord, p_z_cord);
        }
    }


    void Segments::Fill_Grid() {

        for (int i = 0; i < fill_grid_count; i++)
        {
            double x_cord = P->stream_beg + P->particles_dist / 2 + P->particles_dist * i;
            double p_speed = P->particle_speed(x_cord);
            for (double z_cord = last_particles[i] - P->particles_dist; z_cord >= 0; z_cord -= P->particles_dist)
            {
                CreateParticle(x_cord, z_cord);
                last_particles[i] = z_cord;
            }
            last_particles[i] += p_speed;
        }


    }


    void Segments::Fill_Sampling() {

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> dist;
        //std::uniform_real_distribution<double> dist_x(P->stream_beg, P->stream_end), dist_z(0, P->particle_speed(P->area_center));

        double p_x_cord, p_z_cord, p_speed, p_burn_radius, fabs_x, max_z = P->particle_speed(P->area_center);

        for (int pi = P->iterate_particles; pi; --pi)
        {
            p_x_cord = dist(gen) * P->stream_width + P->stream_beg;
            p_z_cord = dist(gen) * max_z;
            p_speed = P->particle_speed(p_x_cord);

            if (p_z_cord < p_speed) {
                CreateParticle(p_x_cord, p_z_cord);
            }
        }

    }

    void Segments::DoSegment(Segment& seg) {
        //if (!(seg.burn_list.empty()))
        for (auto& ok_i : seg.ok_list) {
            for (auto& burn_i : seg.burn_list) {
                if (ok_i.Cross(burn_i)) {
                    will_burn_index.push_back(ok_i.index);
                    break;
                }
            }
        }
    }

    void Segments::CrossParticles() {

        std::for_each(std::execution::par,
            grid.begin(), grid.end(), [this](Segment& seg) {
                if (!(seg.burn_list.empty())) {
                    DoSegment(seg);
                }
            });

        all_will_burn.resize(will_burn_index.size());
        tbb::parallel_for(size_t(0),all_will_burn.size(), [=](size_t i) {
            BurnParticle(all_list[will_burn_index[i]]);
            all_will_burn[i] = &all_list[will_burn_index[i]];
        });

    }




    void Segments::CircleSquareSampling(double r)
    {
        //circle_squares.clear();

        unsigned g = r / grid_x_size;
        circle_squares.resize(g+1);
        circle_squares[0] = g;

        double r2 = r * r;
        double grid_z_size_2 = grid_z_size * grid_z_size;
        for (size_t i = 1; i <= g; i++)
        {
            circle_squares[i] = sqrt(r2 - i * i * grid_z_size_2) / grid_x_size;
        }


    }


    void Segments::StepParticles()
    {

        std::for_each(std::execution::par, all_list.begin(), all_list.end(), [this](Particle& p) {

            CounterType::reference thread_sage_counter = SageCounter.local();

            if (p.state == Particle::State::WARM && ++p.warm_counter >= P->iterations)
            {
                p.state = Particle::State::BURN;
                //++burn_counter;
            }
            else if (p.state == Particle::State::BURN && ++p.burn_counter > P->burn_time && P->burn_time)
            {
                p.state = P->sage_time && ++thread_sage_counter == 5 ? Particle::State::SAGE : Particle::State::DIED;
                thread_sage_counter *= thread_sage_counter < 5;
            }
            else if (p.state == Particle::State::SAGE && ++p.burn_counter >= P->sage_time)
            {
                p.state = Particle::State::DIED;
            }

            });

    }

    void Segments::MoveParticles()
    {
        std::for_each(std::execution::par, all_list.begin(), all_list.end(), [this](Particle& p) {
            MoveParticle(p);
        });

    }



    inline int Segments::GetSegmentX(double x_cord) const
    {
        int x = ceil((x_cord - P->area_beg) * grid_count_x_percent) - 1;
        return x * (x > 0);
    }
    inline int Segments::GetSegmentZ(double z_cord) const
    {
        int z = ceil(z_cord * grid_count_z_percent) - 1;
        return z * (z > 0);
    }
    Segments::Segment& Segments::SegmentByPoint(double x, double z)
    {
        return grids(GetSegmentX(x), GetSegmentZ(z));
    }
    void Segments::BurnSegmentByPoint(double x, double z)
    {
        BurnSegment(grids(GetSegmentX(x), GetSegmentZ(z)));
    }


    void Segments::ClearSegments()
    {
        will_burn_index.clear();
        std::for_each(std::execution::par,
            grid.begin(), grid.end(), [this](Segment& seg) {
                seg.burn_list.clear();
                seg.b_list.clear();
                seg.ok_list.clear();
                seg.burn_indexes.clear();
                seg.front_points.clear();
            });
    }


    void Segments::UpdateSegments()
    {
        ClearSegments();


        tbb::parallel_for(size_t(0),all_list.size(), [=] (size_t i) {
                ParticleToSegment(all_list[i], i);
            }
        );

       std::for_each(std::execution::par,
            grid.begin(), grid.end(), [this](Segment& seg) {
               seg.ok_size = seg.ok_list.size();
               seg.b_size = seg.b_list.size();
            });


    }


    void Segments::ParticleToSegment(Particle& p, size_t index) {
        int seg_x = GetSegmentX(p.x);
        int seg_z = GetSegmentZ(p.z);
        auto &seg = grids(seg_x, seg_z);

        if (refill && seg.ok_list.size() + seg.b_list.size() > grid_particles_max) {
            p.state = Particle::State::DIED;
            return;
        }

        if (p.state == Particle::State::OK) {
            seg.ok_list.emplace_back(p, index);
        }
        else if (p.state == Particle::State::BURN) {
            seg.b_list.emplace_back(p, index);
        }
    }








    void Segments::LightsOut() {
        all_will_burn.clear();
        for (auto& p : all_list) {
            if (p.state != Particle::State::OK)
            {
                p.state = Particle::State::SAGE;
            }
        }
        ClearSegments();
        UpdateSegments();
        return;

        /*for (auto& seg : burn_segments)
        {
            seg->burn_list.clear();
        }*/
    }
    void Segments::EraseParticles() {

        all_will_burn.clear();
        for (auto& p : all_list) {
            p.state = Particle::State::SAGE;
        }
        ClearSegments();
        //        UpdateSegments();
    }
    void Segments::Density_Grid()
    {
        std::string output = "count";
        for (const auto& seg : grid)
        {
            int size = seg.ok_list.size()+ seg.b_list.size();
            if (size) {
                output += fmt::format("\n{}", size);
            }
        }

        std::ofstream csv(P->csv_folder + "d_grid.csv");
        csv << output;
        csv.close();
    }

   


    void Segments::MoveParticle(Particle& p) {
        double half_x = p.x + P->particle_speed_x(p.x, p.z)/2;
        double half_z = p.z + P->particle_speed_z(p.x, p.z)/2;

        p.x += P->particle_speed_x(half_x, half_z);
        p.z += P->particle_speed_z(half_x, half_z);

        if (p.z < 0 || p.z >= P->area_height || p.x > P->area_end || p.x < P->area_beg) p.state = Particle::State::DIED;
    }

    void Segments::StepParticle(Particle& p) {
        //p.Step();

        /*if (p.state == Particle::State::WAVE &&
            ++p.wave_counter >= P->wave_time)
        {
            p.state = Particle::State::DIED;
        }
        else */if (p.state == Particle::State::WARM &&
            ++p.warm_counter >= P->iterations)
        {
            p.state = Particle::State::BURN;
            //++burn_counter;
        }
        else if (p.state == Particle::State::BURN &&
            ++p.burn_counter > P->burn_time)
        {
            p.state = P->sage_time ? Particle::State::SAGE : Particle::State::DIED;
        }

        else if (p.state == Particle::State::SAGE &&
            ++p.burn_counter >= P->sage_time)
        {
            p.state = Particle::State::DIED;
        }
    }
    void Segments::BurnParticle(Particle& particle) {
        particle.setBurn();
    }


    void Segments::BurnSegment(Segment& segment) {
        
        std::for_each(std::execution::par,
            segment.ok_list.begin(), segment.ok_list.end(), [this](SegPointOk &p) {
                BurnParticle(all_list[p.index]);
            });
        segment.ok_list.clear();
    }


    void Segments::ClearParticles() {

        auto erase_it = std::remove_if(std::execution::par, all_list.begin(), all_list.end(), [](Particle& p) {
            return p.state == Particle::State::DIED;
        });
        all_list.erase(erase_it, all_list.end());


    }

    void Segments::CalcBurnRadius(int g)
    {
        std::for_each(std::execution::par, grid.begin(), grid.end(), [&,g](Segment& seg) {
            seg.curvature = 0;
            seg.c_ok = seg.c_b = 0;
            if (seg.b_list.size() == 0) return;
            //if (seg.x >= g && seg.x < grid_count_x - g && seg.z >= g && seg.z < grid_count_z - g) 
            {

                int nbr = 0;
                for (int j = -circle_squares[0], zi = seg.z - circle_squares[0]; zi <= seg.z + circle_squares[0]; zi++,j++)
                {
                    int gi = circle_squares[abs(j)];
                    for (int xi = seg.x - gi; xi <= seg.x + gi; xi++)
                    {
                        if (zi < 0) {
                            seg.c_ok += grid_particles_count;
                            nbr++;
                        }
                        else if (zi >= grid_count_z || xi < 0 || xi >= grid_count_x)
                        {

                            //seg.c_b += grid_particles_count;
                        }
                        else {
                            nbr++;
                            auto& seg_i = grids(xi, zi);
                            seg.c_ok += seg_i.ok_size;
                            seg.c_b += seg_i.b_size;
                        }
                    }
                }



                if (seg.c_ok && seg.c_b)  
                {
                    //seg.curvature = (0.5 - (double)seg.c_ok / (seg.c_b + seg.c_ok)) * M_PI * M_PI / ( P->grid_curve_calc * grid_x_size) * P->curve_burn_coef;
                    //seg.curvature = (0.5 - (double)seg.c_ok / nbr / grid_particles_count);
                    seg.curvature = (0.5 - (double)seg.c_ok / nbr / grid_particles_count) * M_PI * M_PI / P->grid_curve_area * P->curve_burn_coef;
                    //seg.curvature *= seg.curvature *(seg.curvature > 0);
                    seg.curvature *= seg.curvature > 0;
                    //seg.curvature *= std::fabs(seg.curvature);
                }
            }
        });
        SplineGrid();
    }

    void Segments::SplineGrid() {
        /*for (size_t j = 0 ; j < grid_count_z - P->grid_curve_calc*2; j++)
        {
            for (size_t i = 0; i < grid_count_x - P->grid_curve_calc*2; i++)
            {
                spline2d.set(i, j, grids(i + P->grid_curve_calc, j + P->grid_curve_calc).curvature);
            }
        }
        spline2d.fit();*/
    }



    void Segments::PlaceBurned()
    {
        std::for_each(std::execution::par, grid.begin(), grid.end(), [&](Segment& seg) {
            if (seg.c_ok && seg.c_b)
            {
                for (auto& bp : seg.b_list) {


                    //double br = P->make_radius_cross_fix(P->burn_radius * (1 + seg.curvature));

                    //double br = P->burn_radius_cross;
                    double br = P->burn_radius_cross * (1 + seg.curvature);
                    //double br = P->burn_radius_cross * (1 + spline2d.eval(bp.x,bp.z));
                    double br2 = br * br;
                    //bp.r2 = br2;
                    int grids_calc = ceil(br / grid_min_size);

                    int seg_x_start = (seg.x - grids_calc) * (seg.x >= grids_calc);
                    int seg_z_start = (seg.z - grids_calc) * (seg.z >= grids_calc);

                    int seg_x_end = seg.x + grids_calc + 1;
                    if (seg_x_end > grid_count_x) seg_x_end = grid_count_x;
                    int seg_z_end = seg.z + grids_calc + 1;
                    if (seg_z_end > grid_count_z) seg_z_end = grid_count_z;

                    for (int xi = seg_x_start; xi < seg_x_end; xi++) {
                        for (int zi = seg_z_start; zi < seg_z_end; zi++) {
                            if (!grids(xi, zi).ok_list.empty()) {
                                grids(xi, zi).burn_list.emplace_back(bp, br2);
                            }
                        }
                    }
                }
            }
        });
    }



    void Segments::Refill()
    {

        std::for_each(std::execution::par, grid.begin(), grid.end(), [this](Segment& seg) {
            GenType::reference thread_gen = ParticleGenerator.local();

            if (seg.b_list.empty()) {
                int sz = seg.ok_list.size();

                if (sz < grid_particles_min) {
                    for (size_t i = 0; i < grid_particles_count - sz; i++)
                    {
                        double p_x_cord = thread_gen.first(thread_gen.second) * grid_x_size + seg.seg_start_x;
                        double p_z_cord = thread_gen.first(thread_gen.second) * grid_z_size + seg.seg_start_z;
                        Particle p(p_x_cord, p_z_cord);
                        auto it = refilled.push_back(p);
                        seg.ok_list.emplace_back(p, all_list.size() + (it - refilled.begin()));
                    }
                }
            }
            else if (seg.ok_list.empty()) {
                int sz = seg.b_list.size();
                if (sz < grid_particles_min) {
                    for (size_t i = 0; i < grid_particles_count - sz; i++)
                    {
                        double p_x_cord = thread_gen.first(thread_gen.second) * grid_x_size + seg.seg_start_x;
                        double p_z_cord = thread_gen.first(thread_gen.second) * grid_z_size + seg.seg_start_z;
                        Particle p(p_x_cord, p_z_cord);
                        p.state = Particle::State::BURN;
                        p.burn_counter = 1;
                        auto it = refilled.push_back(p);
                        seg.b_list.emplace_back(p, 0);
                    }
                }

            }
            else {

            }

            });


        all_list.insert(all_list.end(), refilled.begin(), refilled.end());
        refilled.clear();

    }

}