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

    void Segments::Fill() {
        _fill_one ? Fill_Grid() : Fill_Sampling();
    }



    void Segments::CreateParticle(double x_cord, double z_cord, double p_speed) {
        double p_burn_radius = P->burn_radius_cross;
        //if (fabs(x_cord) < 0.25) p_burn_radius *= 1 + pow(1 - fabs(x_cord) / 0.25, 4);
        all_list.emplace_back(x_cord, z_cord, p_speed, p_burn_radius);
//        ++particles_count;
    }

    void Segments::SetSegmentsGrid(double seg_size)
    {
        //double seg_size = P->burn_radius_cross;
        grid_count_x = (int)floor(P->area_size / seg_size);
        grid_count_z = (int)floor(P->area_height / seg_size);
        grid_count = grid_count_x * grid_count_z;

        grid_count_x_percent = grid_count_x / P->area_size;
        grid_count_z_percent = grid_count_z / P->area_height;


        grids = std::vector<Segment>(grid_count);
//        grids.resize(grid_count);

        for (int i = 0; i < grid_count; ++i) {
            grids[i].x = i % grid_count_x;
            grids[i].z = i / grid_count_x;
        }

//        grids.resize(grid_count);
//        delete [] mutex_mem;
//        mutex_mem = new segment_mutex[grid_count];

//        grids.clear();
//        for (int i = 0; i < grid_count; ++i) {
//            grids[i].mtx = &mutex_mem[i];
//            grids.push_back({});
//            grids.emplace_back();
//        }

        delete[]grid;
        delete[]grid_mem;

        grid = new Segment* [grid_count_x];
        grid_mem = new Segment [grid_count];

        for (size_t i = 0; i < grid_count_x; i++) {
            grid[i] = grid_mem + i * grid_count_z;
        }

        /*for (int i = 0; i < grid_count; ++i) {
            grid_mem[i].mtx = &mutex_mem[i];
        }*/


    }

    Segments::Segment* Segments::GetSeg(int x, int z) {
        return &(grids[z * grid_count_x + x]);
    }

    void Segments::ResetFillGrid() {

        fill_grid_count = (int)floor(P->stream_width / P->particles_dist);
        delete[]last_particles;
        last_particles = new double[fill_grid_count];
        for (int i = 0; i < fill_grid_count; i++) {
            last_particles[i] = 0;
        }
    }



    void Segments::Print(int num)
    {
        //all_will_burn.clear();

        std::string output;

        //output += "#base_particles " + std::to_string(P->base_particles);

        output += "x,z,burn";
        for (auto& particle : all_list) {
			output += fmt::format("\n{},{},{}", particle._x(), particle._z(), particle.burn_counter);

            /*if (particle.burn_counter == 1 || particle.getState() == Particle::State::WARM)
            {
                all_will_burn.push_back(&particle);
            }*/
        }

        std::ofstream csv(P->csv_folder + "gas.csv." + std::to_string(num));
        csv << output;
        csv.close();



    }



    void Segments::Fill_Grid() {

        for (int i = 0; i < fill_grid_count; i++)
        {
            double x_cord = P->stream_beg + P->particles_dist / 2 + P->particles_dist*i;
            double p_speed = P->particle_speed(x_cord);
            for (double z_cord = last_particles[i] - P->particles_dist; z_cord >= 0; z_cord -= P->particles_dist)
            {
                CreateParticle(x_cord, z_cord, p_speed);
                last_particles[i] = z_cord;
            }
            last_particles[i]+= p_speed;
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
            sum+= (int)grid[xi][0].ok_list.size();
        }
        return sum;
    }


    void Segments::DoSegment(int i) {
        int x = i % grid_count_x, z = i / grid_count_x;
        if (CheckSegmentBurn(x, z))
        {
            for (auto& particle_i : GetSeg(x,z)->ok_list)
//            for (auto& particle_i : grid[x][z].ok_list)
            {
                ParticleInBurnSegment(particle_i, x, z);
            }
        }
    }

    void Segments::DoSegment(Segment *seg) {
        for (auto &particle_i : seg->ok_list) {
            for (auto &burn_i : seg->burn_list) {
                if (particle_i->CrossBurn(*burn_i)) {
                    BurnParticle(particle_i);
                    //seg.will_burn.push_back(particle_i);
                    all_will_burn_concurrent.push_back(particle_i);
                    break;
                }
            }
        }
    }
    void Segments::DoSegment(Segment &seg) {
//        if (!(seg.burn_list.empty()))
        {
            for (auto &particle_i : seg.ok_list) {
                for (auto &burn_i : seg.burn_list) {
                    if (particle_i->CrossBurn(*burn_i)) {
                        BurnParticle(particle_i);
                        //seg.will_burn.push_back(particle_i);
                        all_will_burn_concurrent.push_back(particle_i);
                        break;
                    }
                }
            }
        }
    }
/*    void Segments::DoSegment(Segment &seg) {
        if (CheckSegmentBurn(seg.x, seg.z))
            for (auto& particle_i : seg.ok_list) {
                if (ParticleInBurnSegment(particle_i, seg.x, seg.z)) {
                    BurnParticle(particle_i);
//                    seg.will_burn.push_back(particle_i);
                    all_will_burn_concurrent.push_back(particle_i);
                }
            }

    }*/

    void Segments::DoSegments(int beg, int end) {
        for (int i = beg; i < end; ++i) {
            DoSegment(i);
        }
    }




    void Segments::CrossParticles() {


        all_will_burn_concurrent.clear();


        /*std::for_each(std::execution::par,
            grids.begin(), grids.end(),
            [this](Segment &seg){
           DoSegment(seg);
        });*/
        
        std::for_each(std::execution::par,
            burn_segments.begin(), burn_segments.end(),
            [this](Segment *seg){
           DoSegment(seg);
        });

        all_will_burn.clear();
        all_will_burn.insert(all_will_burn.end(), all_will_burn_concurrent.begin(), all_will_burn_concurrent.end());


//        for(auto &pi : all_will_burn_concurrent) all_will_burn.push_back(pi);
//        all_will_burn.clear();
//        for(auto& seg : grids) all_will_burn.insert(all_will_burn.end(), seg.will_burn.begin(), seg.will_burn.end());


    }




    bool Segments::CheckSegmentBurn(int seg_x, int seg_z)
    {

        for (int xi = seg_x ? seg_x - 1 : 0; xi < (seg_x < grid_count_x - 1 ? seg_x + 2 : grid_count_x); xi++)
        {
            for (int zi = seg_z ? seg_z - 1 : 0; zi < (seg_z < grid_count_z - 1 ? seg_z + 2 : grid_count_z); zi++)
            {
//                if (!(grid[xi][zi].burn_list.empty()))
                if (!(GetSeg(xi,zi)->burn_list.empty()))
                {
                    return true;
                }
            }
        }
        return false;

    }

    bool Segments::ParticleInBurnSegment(Particle* particle, int seg_x, int seg_z)
    {
        //int xi = seg_x;
        for (int xi = seg_x ? seg_x - 1 : 0; xi < (seg_x < grid_count_x - 1 ? seg_x + 2 : grid_count_x); xi++)
        {
            for (int zi = seg_z ? seg_z - 1 : 0; zi < (seg_z < grid_count_z - 1 ? seg_z + 2 : grid_count_z); zi++)
            {
//                auto seg = &grid[xi][zi];
                auto seg = GetSeg(xi,zi);
                if (!seg->burn_list.empty()) for(auto& burn_i : seg->burn_list)
                    {
                        if (particle->CrossBurn(*burn_i))
                        {
                            return true;
                        }
                    }
            }
        }
        return false;
    }

    void Segments::StepParticles()
    {

//        all_will_burn_concurrent.clear();

        std::for_each(std::execution::par, all_list.begin(), all_list.end(), [this](Particle &p) {
            StepParticle(&p);
//            if (p.burn_counter == 1) all_will_burn_concurrent.push_back(&p);
        });

//        all_will_burn.clear();
//        all_will_burn.insert(all_will_burn.end(), all_will_burn_concurrent.begin(), all_will_burn_concurrent.end());


        /*all_will_burn.clear();
        for (auto& particle_i : all_list) {
            StepParticle(&particle_i);
            if (particle_i.burn_counter == 1)
            {
                all_will_burn.push_back(&particle_i);
            }
        }*/
    }
    void Segments::MoveParticles()
    {

        std::for_each(std::execution::par, all_list.begin(), all_list.end(), [this](Particle &p) { MoveParticle(&p); });

        /*for (auto& particle_i : all_list) {
            MoveParticle(&particle_i);
        }*/
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
    Segments::Segment* Segments::SegmentByPoint(double x, double z)
    {
        return GetSeg(GetSegmentX(x),GetSegmentZ(z));
    }
    /*Segments::Segment* Segments::GetSegment(double x, double z) const
    {
        return &grid[GetSegmentX(x)][GetSegmentZ(z)];
    }*/


    void Segments::ClearSegments()
    {
        std::for_each(std::execution::par,
        grids.begin(), grids.end(), [this](Segment &seg){
            seg.has_burn = false;
            seg.burn_list.clear();
            seg.ok_list.clear();
            seg.will_burn.clear();
        });
        burn_segments.clear();
/*        for (int i = 0; i < grid_count; i++)
        {
            grid_mem[i].has_burn = false;
            grid_mem[i].burn_list.clear();
            grid_mem[i].ok_list.clear();
        }*/
    }


    void Segments::UpdateSegments()
    {
        ClearSegments();

//        auto put_particle = [](Particle &p) { ParticleToSegment(&p); };

        std::for_each(
                std::execution::par,
                all_list.begin(), all_list.end(), [this](Particle &p) { ParticleToSegment(&p); });

        std::for_each(std::execution::par,
        grids.begin(), grids.end(),
        [this](Segment &seg){
          if (!(seg.burn_list.empty())) burn_segments.push_back(&seg);
        });




        /*for (auto& p : all_list) {
            ParticleToSegment(&p);
        }*/
    }/*
    void Segments::UpdateSegments()
    {
        burn_segments.clear();
        ClearSegments();

        for (auto& p : all_list) {

            unsigned seg_x = GetSegmentX(p.x);
            unsigned seg_z = GetSegmentZ(p.z);
            Segment *segment = &grid[seg_x][seg_z];
            //Segment *segment = GetSegment(p.x, p.z);

            if (p.state == Particle::State::BURN) {
                burn_segments.emplace_back(seg_x, seg_z);
                segment->has_burn = 1;
                segment->burn_list.push_back(&p);
            }
            else if (p.state == Particle::State::OK) {
                segment->ok_list.push_back(&p);
            }
        }
    }*/

    void Segments::ParticleToSegment(Particle *p) {
        int seg_x = GetSegmentX(p->x);
        int seg_z = GetSegmentZ(p->z);
//        Segment *segment = &grid[seg_x][seg_z];
        Segment *segment = GetSeg(seg_x,seg_z);

        if (p->state == Particle::State::OK) {
            segment->ok_list.push_back(p);
        }
        else if (p->state == Particle::State::BURN) {

            for (int xi = seg_x ? seg_x - 1 : 0; xi < (seg_x < grid_count_x - 1 ? seg_x + 2 : grid_count_x); xi++) {
                for (int zi = seg_z ? seg_z - 1 : 0; zi < (seg_z < grid_count_z - 1 ? seg_z + 2 : grid_count_z); zi++) {
                    GetSeg(xi, zi)->burn_list.push_back(p);
                }
            }

//            segment->burn_list.push_back(p);
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

        for (auto& seg : burn_segments)
        {
            seg->burn_list.clear();
        }
    }
    void Segments::Clear() {

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
//                size = (int)grid[xi][zi].ok_list.size();
                size = (int)GetSeg(xi,zi)->ok_list.size();
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
        for (int zi = 1; zi < grid_count_z-1; zi++)
        {
            for (int xi = 1; xi < grid_count_x-1; xi++)
            {
                for (auto& particle_1 : grid[xi][zi].ok_list)
                {
                    crossed = 0;
                    for (int i = xi - 1; i < xi + 2; i++)
                    {
                        for (int j = zi - 1; j < zi + 2; j++)
                        {
                            for (auto& particle_2 : grid[i][j].ok_list)
                            {
                                if (particle_1->Cross(*particle_2))
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
        for (int zi = 1; zi < grid_count_z-1; zi++)
        {
            for (int xi = 1; xi < grid_count_x-1; xi++)
            {
                for (auto& particle_1 : grid[xi][zi].ok_list)
                {
                    max = 0;
                    for (int i = xi - 1; i < xi + 2; i++)
                    {
                        for (int j = zi - 1; j < zi + 2; j++)
                        {
                            for (auto& particle_2 : grid[i][j].ok_list)
                            {
                                if (particle_1->Cross(*particle_2))
                                {
                                    distance = particle_1->Distance(*particle_2);
                                    if (distance < P->burn_radius_2 && distance > max)
                                    {
                                        max = distance;
                                    }
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


    void Segments::MoveParticle(Particle* p) {
        p->Move();
        if (p->z >= P->area_height) p->state = Particle::State::DIED;
    }

    void Segments::StepParticle(Particle* p) {

        if (p->state == Particle::State::WAVE &&
            ++p->wave_counter >= P->wave_time)
        {
            p->state = Particle::State::DIED;
        }
        else if (p->state == Particle::State::WARM &&
                 ++p->warm_counter >= P->iterations)
        {
            p->state = Particle::State::BURN;
            //++burn_counter;
        }
        else if (p->state == Particle::State::BURN &&
                 ++p->burn_counter > P->burn_time)
        {
            p->state = P->sage_time ? Particle::State::SAGE : Particle::State::DIED;
        }

        else if (p->state == Particle::State::SAGE &&
                 ++p->sage_counter >= P->sage_time)
        {
            p->state = Particle::State::DIED;
        }
    }
    void Segments::BurnParticle(Particle* particle) {
        particle->setBurn();
        //all_will_burn.push_back(particle);
    }


    void Segments::BurnSegment(Segment* segment) {
        for (auto& particle : segment->ok_list) {
            BurnParticle(particle);
        }
    }


    void Segments::ClearParticles() {

        auto erase_it = std::remove_if(std::execution::par, all_list.begin(), all_list.end(), [](Particle p) {
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
            StepParticle(&p);

        }
    }

    void Segments::CalcFrontlineRadius(std::vector <Frontline::front_line_point> &points) {
        int grids_calc = ceil(P->frontline_cross_multipler);
        double cross_radius = P->burn_radius_cross*P->frontline_cross_multipler;

        for (int i = 0; i < P->front_line_steps; ++i) {

            int counter = 0;
            int seg_x = GetSegmentX(points[i].x);
            int seg_z = GetSegmentZ(points[i].z);

            Particle particle(points[i].x, points[i].z, 0, cross_radius);

            /*for (const auto& p : all_list) {
                if (p.state == Particle::State::OK && particle.Cross(p))
                {
                    ++counter;
                }
            }*/
            for (int xi = (seg_x - grids_calc >= 0)*(seg_x - grids_calc);
                xi <= (seg_x + grids_calc)*(seg_x + grids_calc < grid_count_x) +
                    (grid_count_x - 1)*(seg_x + grids_calc >= grid_count_x);
                xi++)
            {
                for (int zi = (seg_z - grids_calc >= 0)*(seg_z - grids_calc);
                    zi <= (seg_z + grids_calc)*(seg_z + grids_calc < grid_count_z) +
                        (grid_count_z - 1)*(seg_x + grids_calc >= grid_count_z);
                    zi++)
                {
//                    for (const auto& ok_particle : grid[xi][zi].ok_list)
                    for (const auto& ok_particle : GetSeg(xi,zi)->ok_list)
                    {
                        int particle_front_seg = ceil((ok_particle->_x() - P->stream_beg)/P->stream_width * P->front_line_steps) - 1;
                        particle_front_seg = particle_front_seg * (particle_front_seg > 0);
                        counter+= particle.Cross(*ok_particle) *
                                (ok_particle->_z() <= points[particle_front_seg].z);
                    }
                }
            }
            points[i].cross = counter / P->base_particles / P->frontline_cross_multipler / P->frontline_cross_multipler;
        }
    }


}