#ifndef MULTITHREADSOLVER_H
#define MULTITHREADSOLVER_H
#include <QThread>
#include "childsolver.h"
#include "task.h"
enum State {occupied, vacant};
class MultiThreadSolver : public ChildSolver, public QThread
{
private:
    bool over;
    int id, nReturns;
    bool poincare;
    
    QVector<real_type> prevIndicatorValue;
    
    QVector<Threshold> thresholds;
protected:
    void solveOneStep();
public:
    MultiThreadSolver(const Task & task, int i_id = 0) : ChildSolver(task.initialVector,task.method,task.tmin,task.tmax, task.h,task.eps, task.mu, task.timeCounts, task.selectIndicator)
    {
        id = i_id;
        over=false;
        nReturns=0;
        for (int i=0; i<nThresholds; ++i) {
            thresholds.push_back(Threshold(thresholdValues[i]));
        }
        if (task.partition == poincarePartition) poincare = true;
        else poincare = false;
    }
    void run() {
        over=false;
        solve();
        over=true;
    }
    bool isOver() {
        return over;
    }

    int getID() {
        return id;
    }

    void setNReturns(int i_nreturns) {
        nReturns = i_nreturns;
    }

    Threshold getThresholdData(int i) {
        return thresholds[i];
    }

    int getNOfThresholds() {
        return thresholds.size();
    }

    void solve();
    void saveSolution(const QDir& workingDirectory, QTextStream& thresholdsTextStream, int prec = 16);
   // std::valarray<real_type> forceNoIndicator(real_type t, const std::valarray<real_type> &r);
  //  std::valarray<real_type> forceFLI(real_type t, const std::valarray<real_type> &r);
};

#endif // MULTITHREADSOLVER_H
