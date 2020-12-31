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
    }


    void Segments::Toggle_Fill() {
        _fill_one = !_fill_one;
    }

    void Segments::FillParticles() {
        _fill_one ? Fill_Grid() : Fill_Sampling();
    }



    void Segments::CreateParticle(double x_cord, double z_cord, double p_speed) {
        double p_burn_radius = P->burn_radius_cross;
        //p_burn_radius *= 1 + (fabs(x_cord) < 0.25) * pow(1.2 - fabs(x_cord) / 0.25 * 1.2, 3);
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
        }




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

    void Segments::PrintSVM(int num, int count)
    {
        //all_will_burn.clear();

        std::string output;

        //output += "#base_particles " + std::to_string(P->base_particles);

        //output += "x,z,burn";
        int step = all_list.size() / count;
        for (size_t i = 0; i < all_list.size(); i+= step)
        {
            output += fmt::format("{} 1:{} 2:{}\n", !all_list[i].isOk() + 1, all_list[i].x, 1 - all_list[i].z);
        }

        std::ofstream csv(P->csv_folder + "svm.csv." + std::to_string(num));
        csv << output;
        csv.close();



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




    /*void Segments::DoSegment(Segment &seg) {
        if (CheckSegmentBurn(seg.x, seg.z))
            for (auto& particle_i : seg.ok_list) {
                if (ParticleInBurnSegment(particle_i, seg.x, seg.z)) {
                    BurnParticle(particle_i);
//                    seg.will_burn.push_back(particle_i);
                    all_will_burn_concurrent.push_back(particle_i);
                }
            }

    }*/

    void Segments::DoSegment(Segment& seg) {
        //if (!(seg.burn_list.empty()))
        {
            for (auto& ok_i : seg.ok_list) {
                for (auto& burn_i : seg.burn_list) {
                    if (ok_i.Cross(burn_i)) {
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


        std::for_each(std::execution::par,
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


        //std::for_each(std::execution::par,
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
        std::for_each(std::execution::par, all_list.begin(), all_list.end(), [this](Particle& p) {
            StepParticle(p);
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
                seg.ok_list.clear();
                seg.burn_indexes.clear();
            });
        //burn_segments.clear();
    }


    void Segments::UpdateSegments()
    {
        ClearSegments();

        //        auto put_particle = [](Particle &p) { ParticleToSegment(&p); };

        /*std::for_each(std::execution::par,
            all_list.begin(), all_list.end(), [this](Particle& p) {
                ParticleToSegment(p);
            }
        );*/

        tbb::parallel_for(size_t(0),all_list.size(), [=] (size_t i) {
                ParticleToSegment(all_list[i], i);
            }
        );

 /*       std::for_each(std::execution::par,
            grid.begin(), grid.end(), [this](Segment& seg) {
                if (!(seg.burn_list.empty())) burn_segments.push_back(&seg);
            });*/


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

        if (p.state == Particle::State::OK) {
            grids(seg_x, seg_z).ok_list.emplace_back(p, index);
        }
        else if (p.state == Particle::State::BURN) {

            int grids_calc = ceil(p.burn_radius / grid_min_size);

            int seg_x_start = (seg_x - grids_calc) * (seg_x >= grids_calc);
            int seg_z_start = (seg_z - grids_calc) * (seg_z >= grids_calc);

            int seg_x_end = seg_x + grids_calc + 1 <= grid_count_x ? seg_x + grids_calc + 1 : grid_count_x;
            int seg_z_end = seg_z + grids_calc + 1 <= grid_count_z ? seg_z + grids_calc + 1 : grid_count_z;

            for (int xi = seg_x_start; xi < seg_x_end; xi++) {
                for (int zi = seg_z_start; zi < seg_z_end; zi++) {
                    grids(xi, zi).burn_list.emplace_back(p, index);
                }
            }

            //segment.burn_list.push_back(p);
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
        int size = 0;
        for (int zi = 0; zi < grid_count_z; zi++)
        {
            for (int xi = 0; xi < grid_count_x; xi++)
            {
                size = (int)grids(xi, zi).ok_list.size();
                if (size) {
                    //					output += fmt::format("\n{}", size);
                }

            }
        }

        std::ofstream csv(P->csv_folder + "d_grid.csv");
        csv << output;
        csv.close();
    }
    void Segments::Density_Radius()
    {
        std::string output = "radius_count, particles_count";
        int crossed = 0;
        std::unordered_map <int, int> denisty_radius;
        for (int zi = 1; zi < grid_count_z - 1; zi++)
        {
            for (int xi = 1; xi < grid_count_x - 1; xi++)
            {
                for (auto& particle_1 : grids(xi, zi).ok_list)
                {
                    crossed = 0;
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
                    //output += fmt::format("\n{}", crossed);
                }

            }
        }
        //        for (auto const& [key, val] : denisty_radius) {
        //			output += fmt::format("\n{},{}", val, key);
        //        }

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
        p.Move();
        if (p.z >= P->area_height) p.state = Particle::State::DIED;
    }

    void Segments::StepParticle(Particle& p) {

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
            ++p.sage_counter >= P->sage_time)
        {
            p.state = Particle::State::DIED;
        }
    }
    void Segments::BurnParticle(Particle& particle) {
        particle.setBurn();
        //all_will_burn.push_back(particle);
    }
    void Segments::BurnParticle(Particle* particle) {
        particle->setBurn();
        //all_will_burn.push_back(particle);
    }


    void Segments::BurnSegment(Segment& segment) {
        /*for (auto& particle : segment.ok_list) {
            will_burn_index.push_back(particle.index);
        }*/
        
        std::for_each(std::execution::par,
            segment.ok_list.begin(), segment.ok_list.end(), [this](SegPoint &p) {
                BurnParticle(all_list[p.index]);
                //all_will_burn_concurrent.push_back(all_list[index]);
            });
        segment.ok_list.clear();
    }


    void Segments::ClearParticles() {

        auto erase_it = std::remove_if(std::execution::par, all_list.begin(), all_list.end(), [](Particle& p) {
            return p.state == Particle::State::DIED;
            });
        all_list.erase(erase_it, all_list.end());


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

    void Segments::CalcFrontlineRadius(std::vector <Point>& points) {
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