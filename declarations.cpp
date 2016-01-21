#include "declarations.h"
void concatValarrays(std::valarray<real_type>& val1, const std::valarray<real_type>& val2)
{
    std::valarray<real_type> temp(val1.size()+val2.size());
    for (int i=0; i<  temp.size(); ++i) {
       temp[i] = i < val1.size() ?  val1[i] : val2[i-val1.size()];
    }
    val1 = temp;
}

real_type cosBetweenVectors(const std::valarray<real_type>& val1, const std::valarray<real_type>& val2)
{
    if (val1.size()!=val2.size()) return NAN;
    return (val1*val2).sum()/mySqrt( (val1*val1).sum() * (val2*val2).sum());
}

QString getIndicatorName(ChaoticIndicator indicator) {
    switch (indicator) {
    case FLI:
        return "FLI";
    case OFLI:
        return "OFLI";
    case MEGNO:
        return "MEGNO";
    case OMEGNO:
        return "OMEGNO";
    case SALI:
        return "SALI";
    case APLE:
        return "APLE";
    case R:
        return "R";
    case allIndicators:
        return "ALL";
    }
    return "";
}

real_type mySqrt(real_type x) {
    if (typeid(x).name() == typeid(real_type).name()) {
        return sqrtq(x);
    }
    if (typeid(x).hash_code() == typeid (long double ).hash_code() || typeid(x).hash_code() == typeid (double ).hash_code()) {
        return sqrt(x);
    }
}

real_type myLog(real_type x) {
    if (typeid(x).hash_code() == typeid (real_type).hash_code()) {
        return logq(x);
    }
    if (typeid(x).hash_code() == typeid (long double ).hash_code() || typeid(x).hash_code() == typeid (double ).hash_code()) {
        return log(x);
    }
}

real_type myPow(real_type x, real_type p) {
    if (typeid(x).hash_code() == typeid (__float128).hash_code()) {
        return powq(x, p);
    }
    if (typeid(x).hash_code() == typeid (long double ).hash_code() || typeid(x).hash_code() == typeid (double ).hash_code()) {
        return pow(x,p );
    }
}
real_type myFabs(real_type x) {
    if (typeid(x).name() == typeid(real_type).name()) {
        return fabsq(x);
    }
    if (typeid(x).hash_code() == typeid (long double ).hash_code() || typeid(x).hash_code() == typeid (double ).hash_code()) {
        return fabs(x);
    }
}
