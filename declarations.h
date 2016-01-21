#ifndef DECLARATIONS_H
#define DECLARATIONS_H
#include <QDir>
#include <iostream>
#include <typeinfo>
#include <valarray>
#include <vector>
#include <QString>
#include <QVector>
#include <quadmath.h>
#include <QDebug>

//typedef __float128 real_type;
typedef long double real_type;
enum MethodOfIntegration {euler, rkutta4, merson5, butcher6, fehlberg8, gragg_bulirsch_stoer};
enum ChaoticIndicator {noIndicator, FLI, OFLI, MEGNO, OMEGNO, SOMEGNO, SALI, APLE, allIndicators, R};
enum TimeIntervalPartition {noTimePartition, uniformPartition, logarithmicPartition, poincarePartition};

const real_type thresholdValues[] = {0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.01, 0.005};
const int nThresholds = 8;

struct Threshold
{
    Threshold(){}
    Threshold(real_type i_value)
    {
        value = i_value;
        timeOfReaching = 0;
    }
    real_type value;
    real_type timeOfReaching;
};
//const real_type mu = 0.5;

//const std::valarray<real_type>rb0({0.0, 0.0, -mu, 0.0});
//const std::valarray<real_type>rb1({0.0, 0.0, 1-mu, 0.0});

void concatValarrays(std::valarray<real_type>& val1, const std::valarray<real_type>& val2);
real_type cosBetweenVectors(const std::valarray<real_type>& val1, const std::valarray<real_type>& val2);
QString getIndicatorName(ChaoticIndicator indicator);

class MyString128 : public QString
{
public:
       static QString number(__float128 x, char f, int prec) {
           char buf[128];
           int succ = quadmath_snprintf(buf,128,"%+-#46.*QE", prec, x);
           QString result(buf);

           return result.trimmed();
       }
};

real_type mySqrt(real_type x);
real_type myLog(real_type x);
real_type myPow(real_type x, real_type p);
real_type myFabs(real_type x);
#endif // DECLARATIONS_H
