#include "parameters.h"

Parameters::Parameters()
{
}

Parameters::Parameters( double maxEnergy ) : maxEnergy(maxEnergy)
{
}

void Parameters::setPotential( std::function<double(double)> potential )
{
    this->potential = potential;
}

void Parameters::findTurningPoints()
{
    // мерсенновский генератор псевдопростых чисел
    std::random_device rd;
    std::mt19937 gen(rd());
    // создаем равномерные генераторы в областях [LOWERBOUND, 0.0] и [0.0, UPPERBOUND]
    std::uniform_real_distribution<double> lowerDistribution(lowerBound, 0.0);
    std::uniform_real_distribution<double> upperDistribution(0.0, upperBound);

    // найдем точки слева и справа от корня в правой половине потенциала (UPPERBOUND)
    double positivePoint = upperDistribution(gen);
    double negativePoint = upperDistribution(gen);

    // генерируем точки до тех пор пока не найдем точку, в которой потенциал больше 0
    // она у нас будет в positivePoint.
    // точка с отрицательным потенциалом будет в negativePoint
    std::cout << "(looking for positivePoint)..." << std::endl;
    for (int i = 0; potential(positivePoint) - maxEnergy < 0.0; i++ )
    {
        positivePoint = upperDistribution(gen);
        if ( i > 1e6 )
            throw std::invalid_argument("Too small bounds!");
    }

    std::cout << "(found positivepoint)" << std::endl;
    std::cout << "(looking for negativePoint)..." << std::endl;
    while (potential(negativePoint) - maxEnergy > 0.0)
        negativePoint = upperDistribution(gen);
    std::cout << "(found negativePoint)" << std::endl;

    // стандартный алгоритм бисекции
    for ( ;; )
    {
        double midPoint = (positivePoint + negativePoint) / 2.0;
        double potentialMidPoint = potential(midPoint);

        if ( fabs(potentialMidPoint - maxEnergy) < epsilon )
        {
            //std::cout << std::fixed << std::setprecision(8);
            //std::cout << "abs(...) = " << fabs(potentialMidPoint - maxEnergy) << std::endl;
            rightTurningPoint = midPoint;
            break;
        }

        if ( (potentialMidPoint - maxEnergy) >= 0.0 )
            positivePoint = midPoint;
        else
            negativePoint = midPoint;
    }

    std::cout << "rightTurningPoint: " << rightTurningPoint << std::endl;

    positivePoint = lowerDistribution(gen);
    negativePoint = lowerDistribution(gen);

    while ( potential(positivePoint) - maxEnergy < 0.0 )
        positivePoint = lowerDistribution(gen);
    while( potential(negativePoint) - maxEnergy > 0.0 )
        negativePoint = lowerDistribution(gen);

    for ( ;; )
    {
        double midPoint = 0.5 * (negativePoint + positivePoint);
        double potentialMidPoint = potential(midPoint);

        if ( fabs(potentialMidPoint - maxEnergy) < epsilon )
        {
            leftTurningPoint = midPoint;
            break;
        }

        if ( (potentialMidPoint - maxEnergy) >= 0.0 )
            positivePoint = midPoint;
        else
            negativePoint = midPoint;
    }

    std::cout << "leftTurningPoint = " << leftTurningPoint << std::endl;
}

void Parameters::setGridParameters()
{
    d = 1.0 / sqrt(2 * mass * maxEnergy);
    N = std::round(2 * (rightTurningPoint / d + 4 * M_PI));

    std::cout << "d = " << d << std::endl;
    std::cout << "N = " << N << std::endl;
}

void Parameters::show()
{
    std::cout << "########################" << std::endl;
    std::cout << "parameters.maxEnergy = " << maxEnergy << std::endl;
    std::cout << "parameters.epsilon = " << epsilon << std::endl;
    std::cout << "parameters.lowerBound = " << lowerBound << std::endl;
    std::cout << "parameters.upperBound = " << upperBound << std::endl;
    std::cout << "########################" << std::endl;
}

