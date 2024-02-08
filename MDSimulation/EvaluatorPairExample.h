#include <iostream>
#include <vector>
#include <cmath>

struct Vector3
{
    float x, y, z;
};

class EvaluatorPairExample
{
public:
    struct param_type
    {
        Scalar k;     //!< Spring constant
        Scalar sigma; //!< Minima of the spring
    };

    DEVICE EvaluatorPairExample(Scalar _rsq, Scalar _rcutsq, const param_type& _params)
        : rsq(_rsq), rcutsq(_rcutsq), k(_params.k), sigma(_params.sigma)
    {
    }

    DEVICE bool evalForceAndEnergy(Scalar& force_divr, Scalar& pair_eng, bool energy_shift)
    {
        if (rsq < rcutsq)
        {
            Scalar r = fast::sqrt(rsq);
            Scalar rinv = 1 / r;
            Scalar overlap = sigma - r;

            force_divr = k * overlap * rinv;
            pair_eng = Scalar(0.5) * k * overlap * overlap;

            if (energy_shift)
            {
                Scalar rcut = fast::sqrt(rcutsq);
                Scalar cut_overlap = sigma - rcut;
                pair_eng -= Scalar(0.5) * k * cut_overlap * cut_overlap;
            }
            return true;
        }
        else
        {
            return false;
        }
    }

protected:
    Scalar rsq;
    Scalar rcutsq;
    Scalar k;
    Scalar sigma;
};

// Potential gradient function
DEVICE void potential_gradient(Scalar rsq, Scalar rcutsq, const EvaluatorPairExample::param_type& params,
                               Scalar& force_divr)
{
    Scalar r = fast::sqrt(rsq);
    Scalar rinv = 1 / r;
    Scalar overlap = params.sigma - r;

    force_divr = params.k * overlap * rinv;

    if (rsq < rcutsq)
    {
        Scalar rcut = fast::sqrt(rcutsq);
        Scalar cut_overlap = params.sigma - rcut;
        force_divr -= params.k * cut_overlap / rcut;
    }
}

// Energy function
DEVICE Scalar energy(Scalar rsq, Scalar rcutsq, const EvaluatorPairExample::param_type& params)
{
    if (rsq < rcutsq)
    {
        Scalar r = fast::sqrt(rsq);
        Scalar overlap = params.sigma - r;

        Scalar pair_eng = Scalar(0.5) * params.k * overlap * overlap;

        Scalar rcut = fast::sqrt(rcutsq);
        Scalar cut_overlap = params.sigma - rcut;
        pair_eng -= Scalar(0.5) * params.k * cut_overlap * cut_overlap;

        return pair_eng;
    }
    else
    {
        return Scalar(0.0);
    }
}

// Exponential force and energy function
std::pair<Vector3, float> exponentialForceParticleEnergy(const std::vector<Vector3>& positions,
                                                        const std::vector<std::vector<int>>& bonds,
                                                        float l0, float b, int particleIndex)
{
    int numParticles = positions.size();
    Vector3 force = {0.0f, 0.0f, 0.0f};
    float energy = 0.0f;

    // Calculate repulsion force and energy for the specified particle
    for (int i = 0; i < numParticles; ++i)
    {
        if (i == particleIndex)
            continue;

        float distance = std::abs(positions[i].x - positions[particleIndex].x) +
                         std::abs(positions[i].y - positions[particleIndex].y) +
                         std::abs(positions[i].z - positions[particleIndex].z);

        force.x -= 4.0f / l0 * std::signbit(positions[i].x - positions[particleIndex].x) *
                    std::exp(-4.0f * distance / l0) / l0;
        force.y -= 4.0f / l0 * std::signbit(positions[i].y - positions[particleIndex].y) *
                    std::exp(-4.0f * distance / l0) / l0;
        force.z -= 4.0f / l0 * std::signbit(positions[i].z - positions[particleIndex].z) *
                    std::exp(-4.0f * distance / l0) / l0;

        energy += std::exp(-4.0f * distance / l0);
    }

    // Calculate bonding forces and energy
    for (const auto& bond : bonds)
    {
        if (particleIndex == bond[0] || particleIndex == bond[1])
        {
            int partner = (particleIndex == bond[0]) ? bond[1] : bond[0];
            float bondForce = 4.0f * b / l0 * std::signbit(positions[particleIndex].x - positions[partner].x) *
                                   std::exp(-4.0f * b * std::abs
