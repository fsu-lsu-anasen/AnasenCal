#ifndef RANDOM_GENERATOR_H
#define RANDOM_GENERATOR_H

#include <random>

//Singleton method, always have a RandomGenerator
class RandomGenerator
{
public:
    ~RandomGenerator();

    std::mt19937& GetGenerator() { return m_generator; }
    static RandomGenerator& GetInstance() { return *s_instance; }

private:
    RandomGenerator();
    static RandomGenerator* s_instance;
    std::mt19937 m_generator; //mersenne-twister
};

#endif