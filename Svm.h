#pragma once
#include <svm.h>
#include <cctype>
#include <vector>

#include "Params.h"
#include "Particle.h"

namespace ps {
	class Svm
	{
	public:
		Svm(const Params& _P) : P(&_P) {
			Init();
		};
		void FreeModel() {
			delete[] prob.x;
			delete[] prob.y;
			svm_free_and_destroy_model(&model);
		};
		~Svm() {
			//delete[] prob.x;
			//delete[] prob.y;
			//svm_free_and_destroy_model(&model);
		};

		void Init();

		void SetParams(const char* input_params);
		void Train(std::vector<Particle>& particles_list);
		void PredictGrid();
		int predict(double x, double z);

		std::vector<double>grid_x;
		std::vector<double>grid_z;
		std::vector<int>grid_predicted;


	private:
		const Params* P;
		svm_parameter param;
		svm_problem prob;
		svm_model* model{ nullptr };
		unsigned grid_width = 0;
		unsigned grid_height = 0;
		unsigned grid_size = 0;
	};

}