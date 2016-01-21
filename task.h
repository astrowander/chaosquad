#ifndef TASK_H
#define TASK_H
#include "declarations.h"
#include <random>
#include <chrono>



struct Task
{

        std::valarray<real_type> initialVector;
        MethodOfIntegration method;
        ChaoticIndicator selectIndicator;
        real_type tmin,tmax, h, eps;
        real_type mu;
        std::vector<real_type> timeCounts;
        int nCounts=20;
        TimeIntervalPartition partition;
        Task(MethodOfIntegration i_method, real_type i_tmin, real_type i_tmax, real_type i_h, real_type i_eps, real_type i_mu, TimeIntervalPartition partition): method(i_method), tmin(i_tmin), tmax(i_tmax), h(i_h), eps(i_eps), mu(i_mu)
        {



        }

        static inline real_type norm3(const std::valarray<real_type> &x) {
            real_type result = x[2]*x[2]+x[3]*x[3];
            return result*mySqrt(result);

        }

        void getInitialVector(real_type x0, real_type C0, bool randomDeviationVector)
        {
            //std::vector<real_type> tempVector;

            initialVector = {0.0, mySqrt(x0*x0 + 2.0*(1.0-mu)/(myFabs(x0+mu)) +2.0*mu/(myFabs(-1.0+mu+x0)) - C0) , x0, 0.0};
            std::valarray<real_type> rb0({0.0,0.0,-mu,0.0});
            std::valarray<real_type> rb1({0.0,0.0,1-mu,0.0});

            //std::valarray<real_type> delta = {0.0, -initialVector[1], initialVector[2]- (initialVector[2]+mu)*(1.0-mu)/norm3(initialVector-rb0) - (initialVector[2]-1+mu)*mu/norm3(initialVector-rb1), 0.0 };
            //real_type squareBracketsExpression = x0 - (1-mu)*(x0+mu)/norm3(initialVector-rb0)-mu*(x0-1+mu)/norm3(initialVector-rb1);
            std::valarray<real_type> delta = {0.0, 0.0,1.0,0.0};
            std::cout << "delta vector: (";
            for (int i=0; i< delta.size(); ++i) {
                std::cout << MyString128::number(delta[i],'g', 10).toStdString() << "; ";
            }
            std::cout<< ")" <<std::endl;
            /*real_type delta_norm = mySqrt(delta[1]*delta[1]+delta[2]*delta[2]);
            std::cout<< "delta norm = "  << MyString128::number(delta_norm,'g', 10).toStdString() << std::endl;
            delta = delta/delta_norm;*/

            if (selectIndicator != ChaoticIndicator::noIndicator && selectIndicator != ChaoticIndicator::R) {
                if (randomDeviationVector) {
                    concatValarrays(initialVector, getRandomDeviationVector());
                }
                else {
                    concatValarrays(initialVector, delta);
                }
            }

            if (selectIndicator == ChaoticIndicator::MEGNO or selectIndicator == ChaoticIndicator::OMEGNO  or selectIndicator == ChaoticIndicator::SOMEGNO) {
                concatValarrays(initialVector, {0.0, 0.0});
            }

            if (selectIndicator == ChaoticIndicator::SALI /*or selectIndicator==SelectIndicator::allIndicators*/) {
               /* if (randomDeviationVector) {
                    concatValarrays(initialVector, getRandomDeviationVector());
                }*/
                //else {
                if (randomDeviationVector) {
                    real_type a3 = initialVector[6], a4 = initialVector[7], b4 = 1/mySqrt(1+(a4/a3)*(a4/a3)), b3 = -a4*b4/a3;
                    concatValarrays(initialVector, {0.0, 0.0, b3, b4});
                }
                else
                {
                    concatValarrays(initialVector, {0.0, 0.0, 1.0, 0.0});
                }
               // }
            }

            if (selectIndicator == ChaoticIndicator::allIndicators) {
                concatValarrays(initialVector, {0.0, 0.0});
            }
          //  initialVector = std::valarray<real_type>(tempVector.data(), tempVector.size());
            std::cout << "Initial Vector: (";
            for (int i=0; i< initialVector.size(); ++i) {
                std::cout << MyString128::number(initialVector[i],'g', 10).toStdString() << "; ";
            }
            std::cout<< ")" <<std::endl;
        }
        void getInitialVectorHH(real_type h, real_type y, real_type y_dot, bool randomDeviationVector)
        {
            initialVector = {mySqrt(2*h - y*y - 2*y*y*y/3 - y_dot*y_dot),y_dot,0.0,y};
            if (selectIndicator != ChaoticIndicator::noIndicator) {
                if (randomDeviationVector) {
                    concatValarrays(initialVector, getRandomDeviationVector());
                }
                else {
                    concatValarrays(initialVector, {1.0, 0.0, 0.0, 0.0});
                }
            }
            if (selectIndicator == ChaoticIndicator::MEGNO or selectIndicator == ChaoticIndicator::OMEGNO) {
                concatValarrays(initialVector, {0.0, 0.0});
            }
            if (selectIndicator == ChaoticIndicator::SALI) {
                if (randomDeviationVector) {
                    concatValarrays(initialVector, getRandomDeviationVector());
                }
                else {
                    concatValarrays(initialVector, {0.0, 0.0, 0.0, 1.0});
                }
            }
        }

        void buildTimePartition()
        {
            timeCounts.clear();
            switch(partition) {
            case noTimePartition:
                timeCounts.push_back(tmax);
                break;
            case uniformPartition:
                for (int i=0; i<nCounts; ++i) {
                    timeCounts.push_back(tmin + (i+1)*(tmax-tmin)/(nCounts+0.0));
                }
                break;
            case logarithmicPartition:
            {
                real_type base = 10*h, time=0;
                int mantiss = 0;
                do {
                  mantiss++;
                  if (mantiss == 100) {
                        mantiss = 10;
                        base *=10;
                  }
                  time = base*mantiss;
                  timeCounts.push_back(time);
                } while (time < tmax - h);
                break;
            }
             default:
                break;
            }
        }

private:
        const std::valarray<real_type> getRandomDeviationVector()
        {
            unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();

            std::default_random_engine generator(seed1);
            std::normal_distribution<double> distribution(0.0, 1.0);
            real_type x[4], r=0;
            for (int i=0; i<4; ++i) {
                x[i]=distribution(generator);
                r+=x[i]*x[i];
            }
            r=mySqrt(r);
            for (int i=0; i<4; ++i) {
                x[i] = x[i] / r;
            }
            return {x[0], x[1], x[2], x[3]};
        }
};



#endif // TASK_H
