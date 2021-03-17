#include "Classifier.h"
#include <random>

namespace ps {



	void Classifier::Init()
	{
		SetParams(P->svm_params.c_str());
		//SetParams("-t 2 -c 100000");

		grid_width = P->screen_width;
		grid_x.resize(grid_width);

		grid_height = P->screen_height;
		grid_z.resize(grid_height);

		grid_size = grid_width * grid_height;
		grid_predicted.resize(grid_size);

		
		for (int i = 0; i < grid_width; i++)
		{
			//grid_x[i] = i / P->screen_x_proportion;
			grid_x[i] = (double)i / grid_width * P->area_size + P->area_beg;
		}
		for (int i = 0; i < grid_height; i++)
		{
			//grid_z[i] = P->area_height - i / P->screen_y_proportion;
			grid_z[i] = P->area_height - (double)i / grid_height * P->area_height;
		}

		int dbg = 0;
		std::cout << "dbd";
		
	}


	void Classifier::SetParams(const char* input_params)
	{

		// default values
		param.svm_type = C_SVC;
		param.kernel_type = RBF;
		param.degree = 3;
		param.gamma = 0;
		param.coef0 = 0;
		param.nu = 0.5;
		param.cache_size = 100;
		param.C = 1;
		param.eps = 1e-3;
		param.p = 0.1;
		param.shrinking = 1;
		param.probability = 0;
		param.nr_weight = 0;
		param.weight_label = NULL;
		param.weight = NULL;

		// parse options
		const char* p = input_params;

		while (1) {
			while (*p && *p != '-')
				p++;

			if (*p == '\0')
				break;

			p++;
			switch (*p++) {
			case 's':
				param.svm_type = atoi(p);
				break;
			case 't':
				param.kernel_type = atoi(p);
				break;
			case 'd':
				param.degree = atoi(p);
				break;
			case 'g':
				param.gamma = atof(p);
				break;
			case 'r':
				param.coef0 = atof(p);
				break;
			case 'n':
				param.nu = atof(p);
				break;
			case 'm':
				param.cache_size = atof(p);
				break;
			case 'c':
				param.C = atof(p);
				break;
			case 'e':
				param.eps = atof(p);
				break;
			case 'p':
				param.p = atof(p);
				break;
			case 'h':
				param.shrinking = atoi(p);
				break;
			case 'b':
				param.probability = atoi(p);
				break;
			case 'w':
				++param.nr_weight;
				param.weight_label = (int*)realloc(param.weight_label, sizeof(int) * param.nr_weight);
				param.weight = (double*)realloc(param.weight, sizeof(double) * param.nr_weight);
				param.weight_label[param.nr_weight - 1] = atoi(p);
				while (*p && !isspace(*p)) ++p;
				param.weight[param.nr_weight - 1] = atof(p);
				break;
			}
		}
		free(param.weight_label);
		free(param.weight);
	}


	void Classifier::Train(std::vector<Particle>& particles_list)
	{
		// guard
		if (particles_list.empty()) return;

		//svm_free_and_destroy_model(&model);

		

		// build problem
		

		prob.l = P->svm_count;
		prob.y = new double[prob.l];

		if (param.kernel_type == PRECOMPUTED)
		{
		}
		else if (param.svm_type == EPSILON_SVR || param.svm_type == NU_SVR)
		{
			//if (param.gamma == 0) param.gamma = 1;
			//svm_node* x_space = new svm_node[2 * prob.l];
			//prob.x = new svm_node * [prob.l];

			//for (int i = 0; i < particles_list.size(); i++)
			//{
			//	x_space[2 * i].index = 1;
			//	x_space[2 * i].value = particles_list[i].x;
			//	x_space[2 * i + 1].index = -1;
			//	prob.x[i] = &x_space[2 * i];
			//	prob.y[i] = particles_list[i].z;
			//}

			//// build model & classify
			//svm_model* model = svm_train(&prob, &param);
			//svm_node x[2];
			//x[0].index = 1;
			//x[1].index = -1;
			//int* j = new int[XLEN];

			//for (i = 0; i < XLEN; i++)
			//{
			//	x[0].value = (double)i / XLEN;
			//	j[i] = (int)(YLEN * svm_predict(model, x));
			//}

			//buffer_painter.setPen(colors[0]);
			//buffer_painter.drawLine(0, 0, 0, YLEN - 1);

			//int p = (int)(param.p * YLEN);
			//for (i = 1; i < XLEN; i++)
			//{
			//	buffer_painter.setPen(colors[0]);
			//	buffer_painter.drawLine(i, 0, i, YLEN - 1);

			//	buffer_painter.setPen(colors[5]);
			//	buffer_painter.drawLine(i - 1, j[i - 1], i, j[i]);

			//	if (param.svm_type == EPSILON_SVR)
			//	{
			//		buffer_painter.setPen(colors[2]);
			//		buffer_painter.drawLine(i - 1, j[i - 1] + p, i, j[i] + p);

			//		buffer_painter.setPen(colors[2]);
			//		buffer_painter.drawLine(i - 1, j[i - 1] - p, i, j[i] - p);
			//	}
			//}

			//svm_free_and_destroy_model(&model);
			//delete[] j;
			//delete[] x_space;
			//delete[] prob.x;
			//delete[] prob.y;
		}
		else
		{
			if (param.gamma == 0) param.gamma = 0.5;
			svm_node* x_space = new svm_node[3 * prob.l];
			prob.x = new svm_node * [prob.l];

			int step = particles_list.size() / P->svm_count;
			for (int i = 0; i < P->svm_count; i++)
			{
				/*if (particles_list[i * step].state == Particle::State::DIED) {
					--i;
					continue;
				}*/
				//if (particles_list[i * step].state == Particle::State::OK || particles_list[i * step].state == Particle::State::BURN)
				{
					x_space[3 * i].index = 1;
					x_space[3 * i].value = particles_list[i * step].x;
					x_space[3 * i + 1].index = 2;
					x_space[3 * i + 1].value = particles_list[i * step].z;
					x_space[3 * i + 2].index = -1;
					prob.x[i] = &x_space[3 * i];
					prob.y[i] = static_cast<int>(particles_list[i * step].state);
				}
			}

			// build model & classify
			model = svm_train(&prob, &param);

			//

			svm_node x[3];
			x[0].index = 1;
			x[1].index = 2;
			x[2].index = -1;



			for (int i = 0; i < grid_width; i++) {
				for (int j = 0; j < grid_height; j++) {
					x[0].value = grid_x[i];
					x[1].value = grid_z[j];
					double d = svm_predict(model, x);
					if (param.svm_type == ONE_CLASS && d < 0) d = 2;
					grid_predicted[j * grid_width + i] = (int)d;
				}
			}


			particles_list.clear();
			std::random_device rd;
			std::mt19937 gen(rd());
			std::uniform_real_distribution<double> dist;

			////for (int pi = 1000; pi; --pi)
			for (int pi = P->full_particles_count; pi; --pi)
			{
				x[0].value = dist(gen) * P->area_size + P->area_beg;
				x[1].value = dist(gen) * P->area_height;
				Particle::State state = static_cast<Particle::State>((int)(svm_predict(model, x)));


				//if (state == Particle::State::OK || state == Particle::State::BURN) 
				{
					particles_list.emplace_back(x[0].value, x[1].value, state);
				}
			}
			delete[] x_space;
		}
		
	}

	int Classifier::predict(double x, double z) {
		svm_node node[3];

		node[0].index = 1;
		node[1].index = 2;
		node[2].index = -1;

		node[0].value = x;
		node[1].value = z;
		double d = 0;
		try
		{
			d = svm_predict(model, node);
		}
		catch (const std::exception&)
		{
			std::cout << "WTF";
		}
		
		if (param.svm_type == ONE_CLASS && d < 0) d = 2;
		return (int)d;
	}

}