#include "RandomGenerator.h"

RandomGenerator* RandomGenerator::s_instance = new RandomGenerator();

RandomGenerator::RandomGenerator()
{
    std::random_device device;
    m_generator.seed(device());
}

RandomGenerator::~RandomGenerator() {}