#pragma once
#include <svm.h>
#include <cctype>
#include <vector>

#include "Particle.h"

namespace ps {
	class SVMClassifier
	{
	public:
		SVMClassifier();
		~SVMClassifier();

		void SetParams(const char* input_params);
		void Train(std::vector<Particle> particles_list);


	private:

	};

}