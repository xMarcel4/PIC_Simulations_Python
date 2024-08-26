#pragma once
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "Fields_.h"
#include "World.h"
#include "Species.h"

namespace Output {
	void fields(World& world, std::vector<Species>& species);
	void screenOutput(World& world, std::vector<Species>& species);
	void diagOutput(World& world, std::vector<Species>& species);
}