#ifndef CHILDSOLVER_H
#define CHILDSOLVER_H
#include "odesolver.h"

#include<QFile>
#include<QTextStream>
#include "qcpcolormapplotdata.h"

QString getIndicatorName(ChaoticIndicator indicator);

class ChildSolver : public OdeSolver
{

public:
    ChildSolver(const std::valarray<real_type> &i_solution, MethodOfIntegration i_method, real_type i_tmin,
                                real_type i_tmax, real_type i_h, real_type i_eps, real_type i_mu, std::vector<real_type> i_timeValues, ChaoticIndicator i_select)
                                                    : OdeSolver(i_solution, i_method, i_tmin, i_tmax, i_h, i_eps, i_timeValues)
    {
       // initialVector = i_solution;
        selectIndicator=i_select;
        h_initial=h;
        globalMegnoParameter=0;
        globalOmegnoParameter=0;
        C0=4;
        mu=i_mu;

        if (i_select==ChaoticIndicator::allIndicators) {
            currentIndicatorValue = QVector<real_type>(6,0);
        }
        else {
            currentIndicatorValue = QVector<real_type>(1,0);
        }

        smallStepCrashFlag = false;

        //if (i_select!=noIndicator) indicatorValues.push_back(0);
        switch (i_select) {
            case ChaoticIndicator::noIndicator:
                pointerToIndicator = nullptr;
                pointerToForce = (&ChildSolver::forceNoIndicator);
                break;
            case ChaoticIndicator::FLI:
                pointerToIndicator = (&ChildSolver::indicatorFLI);
                pointerToForce = (&ChildSolver::forceFLI);
                break;
            case ChaoticIndicator::OFLI:
                pointerToIndicator = (&ChildSolver::indicatorOFLI);
                pointerToForce = (&ChildSolver::forceFLI);
                break;
            case ChaoticIndicator::MEGNO:
                pointerToIndicator = (&ChildSolver::indicatorMEGNO);
                pointerToForce = (&ChildSolver::forceMEGNO);
                break;
            case ChaoticIndicator::OMEGNO:
                pointerToIndicator = (&ChildSolver::indicatorOMEGNO);
                pointerToForce = (&ChildSolver::forceMEGNO);
                break;
        case ChaoticIndicator::SOMEGNO:
            pointerToIndicator = (&ChildSolver::indicatorSOMEGNO);
            pointerToForce = (&ChildSolver::forceMEGNO);
            break;
        case ChaoticIndicator::SALI:
            pointerToIndicator = (&ChildSolver::indicatorSALI);
            pointerToForce = (&ChildSolver::forceSALI);
                break;
        case ChaoticIndicator::APLE:
            pointerToIndicator = (&ChildSolver::indicatorAPLE);
            pointerToForce = (&ChildSolver::forceFLI);
                break;
        case ChaoticIndicator::R:
            pointerToForce = (&ChildSolver::forceNoIndicator);
            pointerToIndicator = (&ChildSolver::indicatorR);
            break;
        case ChaoticIndicator::allIndicators:
            pointerToIndicator = (&ChildSolver::indicatorALL);
            pointerToForce = (&ChildSolver::forceALL);
                break;
        }
       rb0={0.0,0.0,-mu,0.0};
       rb1={0.0,0.0,1-mu,0.0};
    }
    QVector <real_type> getIndicatorValue(int n) {
        return indicatorValues[n];
    }

    QVector <real_type> getCurrentIndicatorValue(int n) {
        return currentIndicatorValue;
    }

    std::vector<QVector <real_type> > getIndicatorValues() {
        return indicatorValues;
    }

    void reverseSolve(); //overloaded
    bool isSmallStepCrash() {
        return smallStepCrashFlag;
    }
    ChaoticIndicator getIndicator()
    {
        return selectIndicator;
    }
    void saveSolution(const QDir& dir, int prec = 16);

protected:
    real_type jacobiIntegral(const std::valarray<real_type> &r);
    real_type C0;
    bool smallStepCrashFlag;
    ChaoticIndicator selectIndicator;
    std::vector< QVector<real_type> > indicatorValues;
    QVector<real_type> currentIndicatorValue;
    void (ChildSolver::*pointerToIndicator)();
    real_type cosa, cosb;

    void solveOneStep();
    void saveSolvingStep();
    void solve();

private:

    //std::valarray<real_type>rb0, rb1;

    std::valarray<real_type> (ChildSolver::*pointerToForce)(real_type, const std::valarray<real_type> &);

    real_type h_initial, mu;
    real_type globalMegnoParameter, globalOmegnoParameter, somegnoParam;


    //real_type mu;
    std::valarray<real_type> rb0;
    std::valarray<real_type> rb1;



    std::valarray<real_type> f(real_type t, const std::valarray<real_type> &x);
    virtual std::valarray<real_type> forceNoIndicator(real_type t, const std::valarray<real_type> &r);
    virtual std::valarray<real_type> forceFLI(real_type t, const std::valarray<real_type> &r);
    std::valarray<real_type> forceMEGNO(real_type t, const std::valarray<real_type> &r);
    std::valarray<real_type> forceSALI(real_type t, const std::valarray<real_type> &r);
    std::valarray<real_type> forceALL(real_type t, const std::valarray<real_type> &r);

    void indicatorR();
    void indicatorFLI();
    void indicatorOFLI();
    void indicatorAPLE();
    void indicatorMEGNO();
    void indicatorOMEGNO();
    void indicatorSOMEGNO();
    void indicatorSALI();
    void indicatorALL();


};

#endif // CHILDSOLVER_H
