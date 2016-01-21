#ifndef ODESOLVER_H
#define ODESOLVER_H

#include <QMessageBox>
#include "declarations.h"

class OdeSolver
{
private:





    MethodOfIntegration method;
    std::valarray<real_type> euler(real_type &t, const std::valarray<real_type>& x, real_type &dt);
    std::valarray<real_type> rkutta4(real_type &t, const std::valarray<real_type>& x, real_type &dt);
    std::valarray<real_type> merson5(real_type &t, const std::valarray<real_type> &x, real_type &dt);
    std::valarray<real_type> butcher6(real_type &t, const std::valarray<real_type> &x, real_type &h);
    std::valarray<real_type> fehlberg8(real_type &t, const std::valarray<real_type> &x, real_type &h);
    std::valarray<real_type> gragg_bulirsch_stoer(real_type &x, const std::valarray<real_type> &y0, real_type &h);

    std::valarray<real_type> graggs_method( const std::valarray<real_type> &y0, real_type x0,real_type x, int number_of_steps );
    int Polynomial_Extrapolation_to_Zero(  std::valarray<real_type> &fzero, std::vector<std::valarray<real_type> > &tableau, real_type x[], const std::valarray<real_type> &f, int n );

protected:
    std::vector< std::valarray<real_type> > solution;
    std::valarray<real_type> initialVector, x_current, x_mid;
    std::valarray<real_type> (OdeSolver::*pointerToMethod)(real_type &, const std::valarray<real_type>&, real_type&);
    std::vector<real_type> timeValues;
    real_type h;
    real_type tmin, tmax;
    real_type t_current;
    std::valarray<real_type> x_prev;
    real_type t_prev;
    real_type savingStep;
    real_type ode_eps;

    bool forbidEnlargingStep;
    int dim;
    int nIterations;
    int size;
    int max_order;
    virtual void solveOneStep();
    virtual void saveSolvingStep();
    inline virtual std::valarray <real_type> f(real_type t, const std::valarray<real_type> &r);    

public:
    OdeSolver(const std::valarray<real_type> &i_solution, MethodOfIntegration i_method, real_type i_tmin, real_type i_tmax, real_type i_h, real_type i_eps, std::vector<real_type> i_timeValues):  ode_eps(i_eps), method(i_method),  h(i_h), tmin(i_tmin), tmax(i_tmax)
    {
        //solution.push_back(i_solution);
       // max_order=0;
        x_current=i_solution;
        t_current=0;
        initialVector = i_solution;
        timeValues = i_timeValues;
        //timeValues.push_back(tmin);
        dim=i_solution.size();
        size=0;
        nIterations=0;
        forbidEnlargingStep = false;
        switch (method) {
            case MethodOfIntegration::euler:
                pointerToMethod = (&OdeSolver::euler);
                break;
            case MethodOfIntegration::rkutta4:
                pointerToMethod = (&OdeSolver::rkutta4);
                break;
            case MethodOfIntegration::merson5:
                pointerToMethod = (&OdeSolver::merson5);
                break;
            case MethodOfIntegration::butcher6:
                pointerToMethod = (&OdeSolver::butcher6);
                break;
            case MethodOfIntegration::fehlberg8:
                pointerToMethod = (&OdeSolver::fehlberg8);
                break;
            case MethodOfIntegration::gragg_bulirsch_stoer:
                pointerToMethod = (&OdeSolver::gragg_bulirsch_stoer);
                break;
        }
    }

    virtual ~OdeSolver() {
       // std::cout << "Максимальный порядок интегрирования = " << 2*max_order+2 << std::endl;
    }
    virtual void solve();
    virtual void reverseSolve();
    real_type getTimeValue(int n) { return timeValues[n];}
    std::valarray<real_type> getSolutionValue(int n) { return solution[n];}
    std::valarray<real_type> getFinalSolution() { return solution.back();}
    std::valarray<real_type> getInitialVector() {return initialVector;}
    std::vector<real_type> getTimeValues() {
        return timeValues;
    }


    int getSolutionSize() { return solution.size();}
    real_type getH() {return h;}

    void setTmax(real_type i_tmax) {tmax=i_tmax;}
    void setH(real_type i_h) {h=i_h;}
    MyString128 printFinalSolution() {
        MyString128 result;
        result += "Solution: ";
        for (int i=0; i<dim; ++i) {
            result += MyString128::number(solution.back()[i],'g',10) + ", ";
        }
        result.chop(2);
        result+=")";
        return result;
    }

};

#endif // ODESOLVER_H
