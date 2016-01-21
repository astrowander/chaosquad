#include "multithreadsolver.h"

void MultiThreadSolver::solveOneStep()
{

//    timeValues.push_back(timeValues.back());
//    real_type& t= timeValues.back();
//    std::valarray<real_type>& x = solution.back();
//    solution.push_back(((*this).*(pointerToMethod))(t,x,h));
    x_prev=x_current;
    t_prev = t_current;
    prevIndicatorValue = currentIndicatorValue;
    x_current=((*this).*(pointerToMethod))(t_current,x_current,h);

    if (selectIndicator!=ChaoticIndicator::noIndicator) ((*this).*(pointerToIndicator))();
  //  qDebug() << "t = " << MyString128::number(t_current, 'g', 10) << "; ";
    ++size;
}

void MultiThreadSolver::solve() {
    int i=0;
    x_current=initialVector;
    x_prev=initialVector;
    t_current=tmin;

    int poincarePoints = 0;
    int nextTimeCount = 0;
    while (t_current < tmax) {
        if ((x_prev[3]*x_current[3]<=0 && (x_current[1] >= 0 || x_prev[1]>=0) && (i > 1 ))) {
            ++poincarePoints;

            if (nReturns!=0) {
                h = - x_prev[3]/x_prev[1];
                forbidEnlargingStep = true;
                while (myFabs(x_current[3])>1e-17) {
                    x_current=x_prev;
                    t_current=t_prev;
                    qDebug() << "h1 = " << MyString128::number(h, 'g', 10);
                    solveOneStep();
                    qDebug() << "h2 = " << MyString128::number(h, 'g', 10);
                    qDebug() << "t = " << MyString128::number(t_current, 'g', 10);
                    qDebug() << "y = " << MyString128::number(x_current[3], 'g', 10);
                    qDebug() << "yprev = " << MyString128::number(x_prev[3], 'g', 10);
                    qDebug() << "";
                    h = h - x_current[3]/x_current[1];
                }
                timeValues.push_back(t_current);
                solution.push_back(x_current);
                if (selectIndicator!=ChaoticIndicator::noIndicator) ((*this).*(pointerToIndicator))();
                indicatorValues.push_back(currentIndicatorValue);
                if (poincarePoints == nReturns) break;
                forbidEnlargingStep = false;
                solveOneStep();
                solveOneStep();
                continue;
            }
        }
        if (!poincare && t_current+h>=timeValues[nextTimeCount]) {
             h = timeValues[nextTimeCount] - t_current;
             std::cout << "t = " << MyString128::number(timeValues[nextTimeCount++],'g', 10).toStdString() << std::endl;
             solveOneStep();             
             continue;
        }

        solveOneStep();
        ++i;

       /* if (x_prev[3]*x_current[3] < 0 ) {
            qDebug() << "Breakpoint";
        }*/
    }
    real_type nominator = 0, denominator = 0;
    for (int j=0; j<timeValues.size(); ++j) {
        if (timeValues[j]>=1.0) {
            nominator += timeValues[j]*timeValues[j];
            denominator += timeValues[j]/indicatorValues[j][0];
        }
    }
    std::cout << "Относительная ошибка интеграла Якоби: " << MyString128::number(myFabs(jacobiIntegral(x_current)-C0)/C0, 'g', 15).toStdString() << std::endl;
    std::cout << "x0 = " << MyString128::number(initialVector[2], 'g', 10).toStdString() << std::endl << "Интегрирование завершено за " << i << " итераций" << std::endl;
    std::cout << "Коэффициент аппроксимации A =  " << MyString128::number(nominator/denominator, 'g', 15).toStdString() << std::endl;
    std::cout << std::endl << "Данные о порогах:" << std::endl;
    for (int j=0; j<thresholds.size(); ++j) {
        std::cout << "Порог = " << MyString128::number(thresholds[j].value,'g',10).toStdString() << "; Время достижения = " << MyString128::number(thresholds[j].timeOfReaching,'g',10).toStdString() << std::endl;
    }
    std::cout << std::endl;
}

void MultiThreadSolver::saveSolution(const QDir &workingDirectory, QTextStream &thresholdsTextStream, int prec)
{
    ChildSolver::saveSolution(workingDirectory, prec);
    thresholdsTextStream << MyString128::number(initialVector[2],'g',10) <<";"<<MyString128::number(cosa,'g',10)<<";"<<MyString128::number(cosb,'g',10);
    for (int i=0; i<thresholds.size(); ++i) {
        thresholdsTextStream << ";" <<QString::number(thresholds[i].timeOfReaching,'g',10);
    }
    thresholdsTextStream << endl;
}

/*std::valarray<real_type> MultiThreadSolver::forceNoIndicator(real_type t, const std::valarray<real_type> &r)
{
    std::valarray<real_type> x(4);
   // std::valarray<real_type>r0(r-rb0);
    //std::valarray<real_type>r1(r-rb1);
    x[0]=   -r[2]-2*r[2]*r[3];
    x[1]= - r[3] -r[2]*r[2] +r[3]*r[3];
    x[2] = r[0];
    x[3] = r[1];
    return x;
}

std::valarray<real_type> MultiThreadSolver::forceFLI(real_type t, const std::valarray<real_type> &r)
{
    std::valarray<real_type> x(dim);
   // std::valarray<real_type>r0(r-rb0);
    //std::valarray<real_type>r1(r-rb1);
    x[0]=   -r[2]-2*r[2]*r[3];
    x[1]= - r[3] -r[2]*r[2] +r[3]*r[3];
    x[2] = r[0];
    x[3] = r[1];

    x[4] = -(1+2*r[3])*r[6] -2*r[2]*r[7];
    x[5] = -2*r[2]*r[6]-(1-2*r[3])*r[7];
    x[6] = r[4];
    x[7] = r[5];

    return x;
}*/
