#include "Segments.h"




namespace ps {

    Segments::Segments(const Params& aParams) : P(&aParams)
    {
        Update();
    }
    void Segments::Update()
    {
        //SetSegmentsGrid(P->burn_radius_cross * (1 + pow(1 - fabs(0) / 0.25, 4)));
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



    void Segments::CreateParticle(double x_cord, double z_cord, double p_speed) {
        //double p_burn_radius = P->get_burn_radius(x_cord);
        double p_burn_radius = P->burn_radius_cross;
        //p_burn_radius *= 1 + (fabs(x_cord) < 0.2);
        //p_burn_radius *= 1 + (fabs(x_cord) < 0.1) * pow(1.15 - fabs(x_cord) / 0.1 / 1.15, 4);
        //p_burn_radius *= 1 + (fabs(x_cord) < 0.2) * pow(1 - fabs(x_cord) / 0.2 * 1, 2) * 1.8;
        //p_burn_radius *= 1 + (fabs(x_cord) < 0.1);
        all_list.emplace_back(x_cord, z_cord, p_speed, p_burn_radius);
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

        //output += "#base_particles " + std::to_string(P->base_particles);

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
        //std::uniform_real_distribution<double> dist_x(P->stream_beg, P->stream_end), dist_z(0, P->particle_speed(P->area_center));

        double p_x_cord, p_z_cord, p_speed, p_burn_radius, fabs_x, max_z = P->particle_speed_z(P->area_center, 0);

        int particles_per_step = round((double)P->base_particles / P->burn_radius_2 / M_PI * P->stream_width * P->stream_width / 10);

        for (int pi = particles_per_step; pi; --pi)
        {
            p_x_cord = dist(gen) * P->stream_width + P->stream_beg;
            p_z_cord = -dist(gen) * P->stream_width / 10;
            p_speed = P->particle_speed(p_x_cord);

            //if (p_z_cord < p_speed) 
            {
                CreateParticle(p_x_cord, p_z_cord, p_speed);
            }
        }
    }


    void Segments::Fill_Grid() {

        for (int i = 0; i < fill_grid_count; i++)
        {
            double x_cord = P->stream_beg + P->particles_dist / 2 + P->particles_dist * i;
            double p_speed = P->particle_speed(x_cord);
            for (double z_cord = last_particles[i] - P->particles_dist; z_cord >= 0; z_cord -= P->particles_dist)
            {
                CreateParticle(x_cord, z_cord, p_speed);
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
                CreateParticle(p_x_cord, p_z_cord, p_speed);
            }
        }

    }



    void Segments::Refill()
    {

        std::for_each(pstl::execution::par, grid.begin(), grid.end(), [this](Segment& seg) {
            GenType::reference thread_gen = ParticleGenerator.local();

            if (seg.b_list.empty()) {
                int sz = seg.ok_list.size();

                if (sz < grid_particles_min) {
                    for (size_t i = 0; i < grid_particles_count - sz; i++)
                    {
                        double p_x_cord = thread_gen.first(thread_gen.second) * grid_x_size + seg.seg_start_x;
                        double p_z_cord = thread_gen.first(thread_gen.second) * grid_z_size + seg.seg_start_z;
                        Particle p(p_x_cord, p_z_cord, 0, 0);
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
                        Particle p(p_x_cord, p_z_cord, 0, 0);
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



    //void Segments::Fill_Sampling_2() {

    //	std::random_device rd;
    //	std::uniform_real_distribution<double> dist_x(P->stream_beg, P->stream_end), dist_z(0, P->particle_speed(P->area_center));

    //	double p_x_cord, p_z_cord, p_speed, p_burn_radius, fabs_x;

    //	for (int pi = P->iterate_particles; pi; --pi)
    //	{
    //		p_x_cord = dist_x(rd);
    //		p_z_cord = dist_z(rd);
    //		p_speed = P->particle_speed(p_x_cord);

    //		if (p_z_cord < p_speed) {
    //			CreateParticle(p_x_cord, p_z_cord, p_speed);
    //		}
    //	}

    //}

    int Segments::Line_Count()
    {
        int sum = 0;
        for (int xi = 0; xi < grid_count_x; xi++)
        {
            sum += (int)grids(xi, 0).ok_list.size();
        }
        return sum;
    }






    void Segments::DoSegment(Segment& seg) {
        //if (!(seg.burn_list.empty()))
        {
            for (auto& ok_i : seg.ok_list) {
                for (auto& burn_i : seg.burn_list) {
                    if (ok_i.Cross(burn_i)) {
                    //if ((P->who_cross ? burn_i.Cross(ok_i) : ok_i.Cross(burn_i)) ) {
                        //seg.burn_indexes.push_back(ok_i.index);
                        //BurnParticle(all_list[ok_i.index]);
                        will_burn_index.push_back(ok_i.index);
                        break;
                    }
                }
            }
            /*for(auto &i : seg.burn_indexes)
                will_burn_index.push_back(i);*/

        }
    }






    void Segments::CrossParticles() {


        //all_will_burn_concurrent.clear();


        std::for_each(pstl::execution::par,
            grid.begin(), grid.end(), [this](Segment& seg) {
                if (!(seg.burn_list.empty()))
                    DoSegment(seg);
            });

        //all_will_burn.clear();
        //all_will_burn.reserve(will_burn_index.size());
        all_will_burn.resize(will_burn_index.size());
        tbb::parallel_for(size_t(0),all_will_burn.size(), [=](size_t i) {
            BurnParticle(all_list[will_burn_index[i]]);
            all_will_burn[i] = &all_list[will_burn_index[i]];
        });


        //std::for_each(pstl::execution::par,
        //    will_burn_index.begin(), will_burn_index.end(), [this](size_t index) {
        //        BurnParticle(all_list[index]);
        //        //all_will_burn_concurrent.push_back(all_list[index]);
        //    });

        

        //all_will_burn.clear();
        //all_will_burn.insert(all_will_burn.end(), all_will_burn_concurrent.begin(), all_will_burn_concurrent.end());



    }




    bool Segments::CheckSegmentNeighborsBurn(int seg_x, int seg_z)
    {

        for (int xi = seg_x ? seg_x - 1 : 0; xi < (seg_x < grid_count_x - 1 ? seg_x + 2 : grid_count_x); xi++)
        {
            for (int zi = seg_z ? seg_z - 1 : 0; zi < (seg_z < grid_count_z - 1 ? seg_z + 2 : grid_count_z); zi++)
            {
                if (!(grids(xi, zi).burn_list.empty()))
                {
                    return true;
                }
            }
        }
        return false;

    }

    void Segments::CircleSquareSampling(double r)
    {
        //circle_squares.clear();

        int g = r / grid_x_size;
        circle_squares.resize(g+1);
        circle_squares[0] = g;

        double r2 = r * r;
        double grid_z_size_2 = grid_z_size * grid_z_size;
        for (size_t i = 1; i <= g; i++)
        {
            circle_squares[i] = sqrt(r2 - i * i * grid_z_size_2) / grid_x_size;
        }


    }

    //bool Segments::ParticleInBurnSegments(Particle* particle, int seg_x, int seg_z)
    //{
    //    for (int xi = seg_x ? seg_x - 1 : 0; xi < (seg_x < grid_count_x - 1 ? seg_x + 2 : grid_count_x); xi++)
    //    {
    //        for (int zi = seg_z ? seg_z - 1 : 0; zi < (seg_z < grid_count_z - 1 ? seg_z + 2 : grid_count_z); zi++)
    //        {
    //            //auto seg = grids(xi,zi);
    //            //if (!seg.burn_list.empty()) 
    //            for (auto& burn_i : grids(xi, zi).burn_list)
    //            {
    //                if (particle.Cross(burn_i))
    //                {
    //                    return true;
    //                }
    //            }
    //        }
    //    }
    //    return false;
    //}

    void Segments::StepParticles()
    {

        std::for_each(pstl::execution::par, all_list.begin(), all_list.end(), [this](Particle& p) {

            CounterType::reference thread_sage_counter = SageCounter.local();

            p.Step();

            /*if (p.state == Particle::State::WAVE &&
                ++p.wave_counter >= P->wave_time)
            {
                p.state = Particle::State::DIED;
            }
            else */if (p.state == Particle::State::WARM && ++p.warm_counter >= P->iterations)
            {
                p.state = Particle::State::BURN;
                //++burn_counter;
            }
            //else if (p.state == Particle::State::BURN && ++p.burn_counter > P->burn_time)
            else if (p.state == Particle::State::BURN && ++p.burn_counter > P->burn_time && P->burn_time)
            {
                p.state = P->sage_time && ++thread_sage_counter == 5 ? Particle::State::SAGE : Particle::State::DIED;
                thread_sage_counter *= thread_sage_counter < 5;
            }

            else if (p.state == Particle::State::SAGE && ++p.burn_counter >= P->sage_time)
            {
                p.state = Particle::State::DIED;
            }

            //StepParticle(p);
            });

    }

    void Segments::MoveParticles()
    {
        std::for_each(pstl::execution::par, all_list.begin(), all_list.end(), [this](Particle& p) {
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
        std::for_each(pstl::execution::par,
            grid.begin(), grid.end(), [this](Segment& seg) {
                seg.burn_list.clear();
                seg.b_list.clear();
                seg.ok_list.clear();
                seg.burn_indexes.clear();
                seg.front_points.clear();
            });
        //burn_segments.clear();
    }


    void Segments::UpdateSegments()
    {
        ClearSegments();

        //        auto put_particle = [](Particle &p) { ParticleToSegment(&p); };

        /*std::for_each(pstl::execution::par,
            all_list.begin(), all_list.end(), [this](Particle& p) {
                ParticleToSegment(p);
            }
        );*/

        tbb::parallel_for(size_t(0),all_list.size(), [=] (size_t i) {
                ParticleToSegment(all_list[i], i);
            }
        );

       std::for_each(pstl::execution::par,
            grid.begin(), grid.end(), [this](Segment& seg) {
               seg.ok_size = seg.ok_list.size();
               seg.b_size = seg.b_list.size();
            });


    }

    //void Segments::ParticleToSegment(Particle& p) {
    //    int seg_x = GetSegmentX(p.x);
    //    int seg_z = GetSegmentZ(p.z);
    //    //Segment segment = grids(seg_x,seg_z);

    //    if (p.state == Particle::State::OK) {
    //        grids(seg_x, seg_z).ok_list.push_back(&p);
    //    }
    //    else if (p.state == Particle::State::BURN) {

    //        int grids_calc = ceil(p.burn_radius / grid_min_size);

    //        int seg_x_start = (seg_x - grids_calc) * (seg_x >= grids_calc);
    //        int seg_z_start = (seg_z - grids_calc) * (seg_z >= grids_calc);

    //        int seg_x_end = seg_x + grids_calc + 1 <= grid_count_x ? seg_x + grids_calc + 1 : grid_count_x;
    //        int seg_z_end = seg_z + grids_calc + 1 <= grid_count_z ? seg_z + grids_calc + 1 : grid_count_z;

    //        for (int xi = seg_x_start; xi < seg_x_end; xi++) {
    //            for (int zi = seg_z_start; zi < seg_z_end; zi++) {
    //                grids(xi, zi).burn_list.push_back(&p);
    //            }
    //        }

    //        //segment.burn_list.push_back(p);
    //    }
    //}


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
            /*int grids_calc = ceil(p.burn_radius / grid_min_size);

            int seg_x_start = (seg_x - grids_calc) * (seg_x >= grids_calc);
            int seg_z_start = (seg_z - grids_calc) * (seg_z >= grids_calc);

            int seg_x_end = seg_x + grids_calc + 1;
            if (seg_x_end > grid_count_x) seg_x_end = grid_count_x;
            int seg_z_end = seg_z + grids_calc + 1;
            if (seg_z_end > grid_count_z) seg_z_end = grid_count_z;

            for (int xi = seg_x_start; xi < seg_x_end; xi++) {
                for (int zi = seg_z_start; zi < seg_z_end; zi++) {
                    grids(xi, zi).burn_list.emplace_back(p, index);
                }
            }*/

            //segment.burn_list.push_back(p);
        }
    }

    void Segments::ParticleToSegment(size_t index) {
        auto& p = all_list[index];
        int seg_x = GetSegmentX(p.x);
        int seg_z = GetSegmentZ(p.z);

        if (p.state == Particle::State::OK) {
            grids(seg_x, seg_z).ok_list.emplace_back(p, index);
        }
        else if (p.state == Particle::State::BURN) {
            /*int grids_calc = ceil(p.burn_radius / grid_min_size);

            int seg_x_start = (seg_x - grids_calc) * (seg_x >= grids_calc);
            int seg_z_start = (seg_z - grids_calc) * (seg_z >= grids_calc);

            int seg_x_end = seg_x + grids_calc + 1;
            if (seg_x_end > grid_count_x) seg_x_end = grid_count_x;
            int seg_z_end = seg_z + grids_calc + 1;
            if (seg_z_end > grid_count_z) seg_z_end = grid_count_z;

            for (int xi = seg_x_start; xi < seg_x_end; xi++) {
                for (int zi = seg_z_start; zi < seg_z_end; zi++) {
                    grids(xi, zi).burn_list.emplace_back(p, index);
                }
            }*/
            grids(seg_x, seg_z).b_list.emplace_back(p, index);

            //segment.burn_list.push_back(p);
        }
    }

    void Segments::CalcFrontlineRadius(std::vector <Point>& points) {
        //frontline_points = points;


        double cross_radius = P->frontline_cross_radius;
        double cross_radius_2 = P->frontline_cross_radius_2;
        //double cross_radius = P->burn_radius_cross;
        //double cross_radius_2 = P->burn_radius_2_cross;
        int grids_calc = ceil(cross_radius / grid_min_size);

        //tbb::parallel_for(size_t(0), frontline_points.size(), [=](size_t i)
        
        for (size_t i = 0; i < points.size(); i++)
        {
                //ParticleToSegment(all_list[i], i);
                //auto p = frontline_points[i];
                if (points[i].x > P->stream_beg + cross_radius && points[i].x < P->stream_end - cross_radius)
                    FrontPointToSegment(points[i], i, grids_calc, cross_radius_2);
        }
        //);



        std::for_each(pstl::execution::par_unseq, grid.begin(), grid.end(), [this](Segment& seg) {
        //for(auto& seg : grid)
            for (auto& sp : seg.front_points) {
                for (const auto& op : seg.ok_list) {
                    if (sp.CrossOk(op)) {
                        ++sp.count;
                    }
                }
            }
        });


        front_crosses.resize(points.size());
        std::fill(front_crosses.begin(), front_crosses.end(), 0);

        for (const auto& seg : grid) {
            for (auto& sp : seg.front_points) {
                front_crosses[sp.index] += sp.count;
            }
        }


        /*for (int i = 0; i < front_crosses.size(); i++) {
            if (front_crosses[i] > 0) {
                for (int j = i-1; j >= 0; j--) {
                    front_crosses[j] = front_crosses[i];
                }
                break;
            }
        }
        for (int i = front_crosses.size() - 1; i >= 0; i--) {
            if (front_crosses[i] > 0) {
                for (int j = i+1; j < front_crosses.size(); j++) {
                    front_crosses[j] = front_crosses[i];
                }
                break;
            }
        }*/

        for (auto& c : front_crosses)
        {
            c = (P->frontline_cross_border - c / P->base_particles / P->frontline_cross_radius_2 * P->burn_radius_2_cross) * M_PI * M_PI / P->frontline_cross_radius;
            //c *= c < 0.5F;
            //if (c) c = 0.5F - c;

        }



    }

    void Segments::FrontPointToSegment(Point& p, size_t index, int grids_calc, double cross_radius_2) {
        int seg_x = GetSegmentX(p.x);
        int seg_z = GetSegmentZ(p.z);

        //int seg_x_start = (seg_x - grids_calc) * (seg_x >= grids_calc);
        //int seg_z_start = (seg_z - grids_calc) * (seg_z >= grids_calc);

        int seg_x_start = seg_x - grids_calc;
        if (seg_x_start < 0) return;
        int seg_z_start = seg_z - grids_calc;
        if (seg_z_start < 0) return;

        int seg_x_end = seg_x + grids_calc + 1;
        if (seg_x_end > grid_count_x) return;// seg_x_end = grid_count_x;
        int seg_z_end = seg_z + grids_calc + 1;
        if (seg_z_end > grid_count_z) return;// seg_z_end = grid_count_z;

        for (int xi = seg_x_start; xi < seg_x_end; xi++) {
            for (int zi = seg_z_start; zi < seg_z_end; zi++) {
                grids(xi, zi).front_points.emplace_back(p, cross_radius_2, index);
            }
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
    void Segments::Density_Radius()
    {

        std::map <int, int> denisty_radius;
        for (int zi = 1; zi < grid_count_z - 1; zi++)
        {
            for (int xi = 1; xi < grid_count_x - 1; xi++)
            {
                for (auto& particle_1 : grids(xi, zi).ok_list)
                {
                    int crossed = 0;
                    for (int i = xi - 1; i < xi + 2; i++)
                    {
                        for (int j = zi - 1; j < zi + 2; j++)
                        {
                            for (auto& particle_2 : grids(i, j).ok_list)
                            {
                                if (particle_1.Cross(particle_2))
                                {
                                    ++crossed;
                                }
                            }
                        }
                    }
                    denisty_radius[crossed] += 1;
                }

            }
        }

        std::string output = "radius_count, particles_count";
        for (auto const& [key, val] : denisty_radius) {
            output += fmt::format("\n{},{}", val, key);
        }

        std::ofstream csv(P->csv_folder + "d_radius.csv");
        csv << output;
        csv.close();
    }
    void Segments::Max_Radius()
    {
        std::string output = "count";
        double max = 0, distance;
        for (int zi = 1; zi < grid_count_z - 1; zi++)
        {
            for (int xi = 1; xi < grid_count_x - 1; xi++)
            {
                for (auto& particle_1 : grids(xi, zi).ok_list)
                {
                    max = 0;
                    for (int i = xi - 1; i < xi + 2; i++)
                    {
                        for (int j = zi - 1; j < zi + 2; j++)
                        {
                            for (auto& particle_2 : grids(i, j).ok_list)
                            {
                                if (particle_1.Cross(particle_2))
                                {
                                    /*distance = particle_1.Distance(particle_2);
                                    if (distance < P->burn_radius_2 && distance > max)
                                    {
                                        max = distance;
                                    }*/
                                }
                            }
                        }
                    }
                    //					output+= fmt::format("\n{}", sqrt(max));
                }

            }
        }

        std::ofstream csv(P->csv_folder + "d_distance.csv");
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
        p.Step();

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
        particle.SetBurnRadius(P->get_burn_radius(particle.x));
        //all_will_burn.push_back(particle);
    }
    void Segments::BurnParticle(Particle* particle) {
        particle->setBurn();
        particle->SetBurnRadius(P->get_burn_radius(particle->x));
        //all_will_burn.push_back(particle);
    }


    void Segments::BurnSegment(Segment& segment) {
        /*for (auto& particle : segment.ok_list) {
            will_burn_index.push_back(particle.index);
        }*/
        
        std::for_each(pstl::execution::par,
            segment.ok_list.begin(), segment.ok_list.end(), [this](SegPoint &p) {
                BurnParticle(all_list[p.index]);
                //all_will_burn_concurrent.push_back(all_list[index]);
            });
        segment.ok_list.clear();
    }


    void Segments::ClearParticles() {

        auto erase_it = std::remove_if(pstl::execution::par, all_list.begin(), all_list.end(), [](Particle& p) {
            return p.state == Particle::State::DIED;
        });
        all_list.erase(erase_it, all_list.end());


    }

    void Segments::CalcBurnRadius(int g)
    {
        std::for_each(pstl::execution::par, grid.begin(), grid.end(), [&,g](Segment& seg) {
            seg.curvature = 0;
            seg.c_ok = seg.c_b = 0;
            if (seg.b_list.size() == 0) return;
            //if (seg.x >= g && seg.x < grid_count_x - g && seg.z >= g && seg.z < grid_count_z - g) 
            {
                /*for (size_t zi = seg.z-g; zi <= seg.z+g ; zi++)
                {
                    for (size_t xi = seg.x-g; xi <= seg.x+g; xi++)
                    {
                        seg.c_ok += grids(xi, zi).ok_list.size();
                        seg.c_b += grids(xi, zi).b_list.size();
                    }
                }*/

                //int xs[] = { 1,3,4,5,6,6,7,7,7,6,6,5,4,3,1 };
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
        std::for_each(pstl::execution::par, grid.begin(), grid.end(), [&](Segment& seg) {
            if (seg.c_ok && seg.c_b)
            {
                for (auto& bp : seg.b_list) {


                    //double br = P->make_radius_cross_fix(P->burn_radius * (1 + seg.curvature));

                    //double br = P->burn_radius_cross;
                    double br = P->burn_radius_cross * (1 + seg.curvature);
                    //double br = P->burn_radius_cross * (1 + spline2d.eval(bp.x,bp.z));
                    double br2 = br * br;
                    bp.r2 = br2;
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
                                grids(xi, zi).burn_list.push_back(bp);
                            }
                        }
                    }
                }
            }
        });
    }


    void Segments::RefractParticles()
    {
        for (auto& p : all_list) {

            if (p.state == Particle::State::OK && p.z >= P->refract_func(p.x)) {
                p.state = Particle::State::WAVE;
            }
            StepParticle(p);

        }
    }

    void Segments::CalcFrontlineRadius2(std::vector <Point>& points) {
        double cross_radius = P->frontline_cross_radius;
        int grids_calc = ceil(cross_radius / grid_min_size);
        /*double cross_radius = P->burn_radius * P->frontline_cross_multipler;
        int grids_calc = ceil(P->frontline_cross_multipler);*/
        //std::cout << cross_radius;

        int forntline_start = 0, frontline_end = P->frontline_window_steps - 1;

        if (P->frontline_cross_chunk)
        {
            while ( points[forntline_start].z < cross_radius && forntline_start < P->frontline_window_steps - 1 )
                ++forntline_start;

            while ( points[frontline_end].z < cross_radius && frontline_end > 1 )
                --frontline_end;
        }





        for (int i = forntline_start; i <= frontline_end; ++i) {

            int counter = 0;
            int seg_x = GetSegmentX(points[i].x);
            int seg_z = GetSegmentZ(points[i].z);

            Particle particle(points[i].x, points[i].z, 0, cross_radius);

            for (const auto& p : all_list) {
                //counter+= p.state == Particle::State::OK && particle.Cross(p);

                int particle_front_seg = ceil((p.x - P->stream_beg) / P->stream_width * P->frontline_window_steps) - 1;
                particle_front_seg*= particle_front_seg > 0;

                if (p.state == Particle::State::OK && particle.Cross(p))
                {
                    ++counter;
                    /*int particle_front_seg = ceil((p.x - P->stream_beg) / P->stream_width * P->frontline_window_steps) - 1;
                    particle_front_seg = particle_front_seg * (particle_front_seg > 0);
                    counter += (p.z <= points[particle_front_seg].z);*/
                }
            }
            //for (int xi = (seg_x - grids_calc >= 0) * (seg_x - grids_calc);
            //    xi <= (seg_x + grids_calc) * (seg_x + grids_calc < grid_count_x) +
            //    (grid_count_x - 1) * (seg_x + grids_calc >= grid_count_x);
            //    xi++)
            //{
            //    for (int zi = (seg_z - grids_calc >= 0) * (seg_z - grids_calc);
            //        zi <= (seg_z + grids_calc) * (seg_z + grids_calc < grid_count_z) +
            //        (grid_count_z - 1) * (seg_x + grids_calc >= grid_count_z);
            //        zi++)
            //    {
            //        for (const auto& ok_particle : grids(xi, zi).ok_list)
            //        {
            //            counter += particle.Cross(ok_particle) && ok_particle->state == Particle::State::OK;

            //            /*if (particle.Cross(ok_particle) && ok_particle->state == Particle::State::OK) {
            //                ++counter;
            //            }*/
            //            /*int particle_front_seg = ceil((ok_particle->x - P->stream_beg) / P->stream_width * P->frontline_window_steps) - 1;
            //            particle_front_seg = particle_front_seg * (particle_front_seg > 0);
            //            counter += particle.Cross(ok_particle) *(ok_particle->z <= points[particle_front_seg].z);*/
            //        }
            //    }
            //}
            //points[i].cross = counter / P->base_particles / P->frontline_cross_multipler / P->frontline_cross_multipler;
            //points[i].cross = counter / P->base_particles / cross_radius / cross_radius * P->burn_radius_2_cross;

            //points[i].cross = counter / P->base_particles / cross_radius / cross_radius * P->burn_radius_2_cross;

            //points[i].cross = counter;
            /*points[i].cross -= 0.5;
            points[i].cross *= 100;*/
        }
    }


}