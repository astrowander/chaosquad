#include "childsolver.h"
static real_type n=1.0;

static inline real_type norm(const std::valarray<real_type> &x) {
    return mySqrt(x[2]*x[2]+x[3]*x[3]);
}

static inline real_type norm3(const std::valarray<real_type> &x) {
    real_type result = x[2]*x[2]+x[3]*x[3];
    return result*mySqrt(result);

}

static inline real_type norm2(const std::valarray<real_type> &x) {
    return x[2]*x[2]+x[3]*x[3];
}

static inline real_type norm5(const std::valarray<real_type> &x) {
    real_type result = x[2]*x[2]+x[3]*x[3];
    return result*result*mySqrt(result);
}

static inline std::valarray <real_type> project(const std::valarray<real_type> &r, const std::valarray<real_type> f)
{
    real_type f_module =f[0]*f[0]+f[1]*f[1]+f[2]*f[2]+f[3]*f[3];
    return r - f * ((r*f).sum()/f_module);
}

std::valarray<real_type> ChildSolver::f(real_type t, const std::valarray<real_type> &x)
{
   return ((*this).*(pointerToForce))(t,x);
}

std::valarray<real_type> ChildSolver::forceNoIndicator(real_type t, const std::valarray<real_type> &r)
{
    std::valarray<real_type> x(4);
    std::valarray<real_type>r0(r-rb0);
    std::valarray<real_type>r1(r-rb1);
    x[0]=   2 * r[1] *n + r[2] - (1 - mu)*r0[2]/norm3(r0) - mu*r1[2]/norm3(r1);
    x[1]= - 2 * r[0]*n + r[3] - (1-mu)*r0[3]/norm3(r0) - mu*r1[3]/norm3(r1);
    x[2] = r[0];
    x[3] = r[1];
    return x;
}

std::valarray<real_type> ChildSolver::forceFLI(real_type t, const std::valarray<real_type> &r)
{
    std::valarray<real_type> f(dim);
    std::valarray<real_type>r0(r-rb0);
    std::valarray<real_type>r1(r-rb1);

   /* std::cout<< "r: " << r[0] << " "<< r[1] << " "<< r[2] << " "<< r[3] << " " <<std::endl;
    std::cout<< "rb0: " << rb0[0] << " "<< rb0[1] << " "<< rb0[2] << " "<< rb0[3] << " " <<std::endl;
    std::cout<< "rb1: " << rb1[0] << " "<< rb1[1] << " "<< rb1[2] << " "<< rb1[3] << " " <<std::endl;
    std::cout<< "r0: " << r0[0] << " "<< r0[1] << " "<< r0[2] << " "<< r0[3] << " " <<std::endl;
    std::cout<< "r1: " << r1[0] << " "<< r1[1] << " "<< r1[2] << " "<< r1[3] << " " <<std::endl;*/

    f[0]=   real_type(2.0) * r[1] *n + r[2] - (real_type(1.0) - mu)*r0[2]/norm3(r0) - mu*r1[2]/norm3(r1);
    f[1]= - real_type(2.0) * r[0]*n + r[3] - (real_type(1.0)-mu)*r0[3]/norm3(r0) - mu*r1[3]/norm3(r1);
    f[2] = r[0];
    f[3] = r[1];

    real_type D11= real_type(1.0) - (real_type(1.0)-mu)*(norm2(r0) - real_type(3.0)*(r[2]+mu)*(r[2]+mu))/norm5(r0) - mu*( norm2(r1) - real_type(3.0) * (r[2]-real_type(1.0)+mu)*(r[2]-real_type(1.0)+mu)) /norm5(r1);
    real_type D12 = (real_type(1.0)-mu)*real_type(3.0)*r[3]*(r[2]+mu)/norm5(r0) + mu*real_type(3.0)*r[3]*(r[2]-real_type(1.0)+mu)/norm5(r1);
    real_type D21= (real_type(1.0)-mu)*real_type(3.0)*r[3]*(r[2]+mu)/norm5(r0) + mu*real_type(3.0)*r[3]*(r[2]-real_type(1.0)+mu)/norm5(r1);
    real_type D22 = real_type(1.0) - (real_type(1.0)-mu)*(norm2(r0) - real_type(3.0)*r[3]*r[3])/norm5(r0) - mu*(norm2(r1) - real_type(3.0)*r[3]*r[3])/norm5(r1);

    real_type D14 = real_type(2.0)*n;
    real_type D23 = -real_type(2.0)*n;

    f[4] = D11*r[6]+D12*r[7]+D14*r[5];
    f[5] = D21*r[6]+D22*r[7]+D23*r[4];
    f[6] = r[4];
    f[7] = r[5];

    return f;
}

std::valarray<real_type> ChildSolver::forceMEGNO(real_type t, const std::valarray<real_type> &r)
{
    std::valarray<real_type> f(dim);
    f = forceFLI(t,r);
    f[8]=  globalMegnoParameter; //teta
    f[9]= (t==0) ? real_type(0.0) : r[8]/t;
    return f;
}

std::valarray<real_type> ChildSolver::forceSALI(real_type t, const std::valarray<real_type> &r)
{
    std::valarray<real_type> f(dim);
    std::valarray<real_type>r0(r-rb0);
    std::valarray<real_type>r1(r-rb1);
    real_type D11= real_type(1.0) - (real_type(1.0)-mu)*(norm2(r0) - real_type(3.0)*(r[2]+mu)*(r[2]+mu))/norm5(r0) - mu*( norm2(r1) - real_type(3.0) * (r[2]-real_type(1.0)+mu)*(r[2]-real_type(1.0)+mu)) /norm5(r1);
    real_type D12 = (real_type(1.0)-mu)*real_type(3.0)*r[3]*(r[2]+mu)/norm5(r0) + mu*real_type(3.0)*r[3]*(r[2]-real_type(1.0)+mu)/norm5(r1);
    real_type D21= (real_type(1.0)-mu)*real_type(3.0)*r[3]*(r[2]+mu)/norm5(r0) + mu*real_type(3.0)*r[3]*(r[2]-real_type(1.0)+mu)/norm5(r1);
    real_type D22 = real_type(1.0) - (real_type(1.0)-mu)*(norm2(r0) - real_type(3.0)*r[3]*r[3])/norm5(r0) - mu*(norm2(r1) - real_type(3.0)*r[3]*r[3])/norm5(r1);

    real_type D14 = real_type(2.0)*n;
    real_type D23 = -real_type(2.0)*n;

    f[0]=   real_type(2.0) * r[1] *n + r[2] - (real_type(1.0) - mu)*r0[2]/norm3(r0) - mu*r1[2]/norm3(r1);
    f[1]= - real_type(2.0) * r[0]*n + r[3] - (real_type(1.0)-mu)*r0[3]/norm3(r0) - mu*r1[3]/norm3(r1);
    f[2] = r[0];
    f[3] = r[1];
    f[4] = D11*r[6]+D12*r[7]+D14*r[5];
    f[5] = D21*r[6]+D22*r[7]+D23*r[4];
    f[6] = r[4];
    f[7] = r[5];

    f[8] = D11*r[10]+D12*r[11]+D14*r[9];
    f[9] = D21*r[10]+D22*r[11]+D23*r[8];
    f[10] = r[8];
    f[11] = r[9];




    return f;
}

std::valarray<real_type> ChildSolver::forceALL(real_type t, const std::valarray<real_type> &r)
{
    std::valarray<real_type> f(dim);
    std::valarray<real_type>r0(r-rb0);
    std::valarray<real_type>r1(r-rb1);
    real_type D11= real_type(1.0) - (real_type(1.0)-mu)*(norm2(r0) - real_type(3.0)*(r[2]+mu)*(r[2]+mu))/norm5(r0) - mu*( norm2(r1) - real_type(3.0) * (r[2]-real_type(1.0)+mu)*(r[2]-real_type(1.0)+mu)) /norm5(r1);
    real_type D12 = (real_type(1.0)-mu)*real_type(3.0)*r[3]*(r[2]+mu)/norm5(r0) + mu*real_type(3.0)*r[3]*(r[2]-real_type(1.0)+mu)/norm5(r1);
    real_type D21= (real_type(1.0)-mu)*real_type(3.0)*r[3]*(r[2]+mu)/norm5(r0) + mu*real_type(3.0)*r[3]*(r[2]-real_type(1.0)+mu)/norm5(r1);
    real_type D22 = real_type(1.0) - (real_type(1.0)-mu)*(norm2(r0) - real_type(3.0)*r[3]*r[3])/norm5(r0) - mu*(norm2(r1) - real_type(3.0)*r[3]*r[3])/norm5(r1);

    real_type D14 = real_type(2.0)*n;
    real_type D23 = -real_type(2.0)*n;

    f[0]=   real_type(2.0) * r[1] *n + r[2] - (real_type(1.0) - mu)*r0[2]/norm3(r0) - mu*r1[2]/norm3(r1);
    f[1]= - real_type(2.0) * r[0]*n + r[3] - (real_type(1.0)-mu)*r0[3]/norm3(r0) - mu*r1[3]/norm3(r1);
    f[2] = r[0];
    f[3] = r[1];
    f[4] = D11*r[6]+D12*r[7]+D14*r[5];
    f[5] = D21*r[6]+D22*r[7]+D23*r[4];
    f[6] = r[4];
    f[7] = r[5];

   /* f[8] = D11*r[10]+D12*r[11]+D14*r[9];
    f[9] = D21*r[10]+D22*r[11]+D23*r[8];
    f[10] = r[8];
    f[11] = r[9];

    f[12] = globalMegnoParameter;
    f[13]= (t==0) ? real_type(0.0) : r[12]/t;*/

    f[8] = globalOmegnoParameter;
    f[9]= (t==0) ? real_type(0.0) : r[8]/t;

    return f;
}

void ChildSolver::indicatorALL()
{
    std::valarray<real_type>& x = x_current;
    real_type& t = t_current;
    std::valarray<real_type> delta1( {x[4], x[5], x[6], x[7]} );
   // std::valarray<real_type> delta2( {x[8], x[9], x[10], x[11]} );

    //real_type dev1 = mySqrt((delta1*delta1).sum());
    //real_type dev2 = mySqrt((delta2*delta2).sum());
    std::valarray<real_type> f(forceNoIndicator(0.0,x));
    real_type odev = mySqrt((project(delta1,f)*project(delta1,f)).sum());
    //FLI
    /*if (myLog(dev1) > currentIndicatorValue[0]) {
       currentIndicatorValue[0] = myLog(dev1);
    }*/
    //OFLI
    currentIndicatorValue[0] = myLog(odev);
    if (myLog(odev) > currentIndicatorValue[1]) {
       currentIndicatorValue[1] = myLog(odev);
    }
    //APLE
    //currentIndicatorValue[2] = myLog(dev1)/myLog(t);
    //MEGNO
    /*globalMegnoParameter=myLog(dev1);
    currentIndicatorValue[3] =  real_type(2.0)*(x[12]-x[13])/t;*/
    //OMEGNO
    globalOmegnoParameter=myLog(odev);
    currentIndicatorValue[4] = real_type(2.0)*(x[14]-x[15])/t;
    //SALI
    /*real_type d1 = mySqrt((x[8]/dev2-x[4]/dev1)*(x[8]/dev2-x[4]/dev1) + (x[9]/dev2-x[5]/dev1)*(x[9]/dev2-x[5]/dev1) + (x[10]/dev2-x[6]/dev1)*(x[10]/dev2-x[6]/dev1) + (x[11]/dev2-x[7]/dev1)*(x[11]/dev2-x[7]/dev1)   );
    real_type d2 = mySqrt((x[8]/dev2+x[4]/dev1)*(x[8]/dev2+x[4]/dev1) + (x[9]/dev2+x[5]/dev1)*(x[9]/dev2+x[5]/dev1) + (x[10]/dev2+x[6]/dev1)*(x[10]/dev2+x[6]/dev1) + (x[11]/dev2+x[7]/dev1)*(x[11]/dev2+x[7]/dev1)   );
    currentIndicatorValue[5] =  (d1<=d2) ? d1 : d2 ;*/

}

void ChildSolver::indicatorR()
{
    if (solution.size()==0) return;
    real_type x1, vx1, x2, vx2;
    x1 = initialVector[2];
    vx1 = initialVector[0];

    x2 = solution.back()[2];
    vx2 = solution.back()[0];

    //currentIndicatorValue[0] = mySqrt((x2-x1)*(x2-x1) + (vx2-vx1)*(vx2-vx1));
    currentIndicatorValue[0] = (x2-x1)*(x2-x1) + (vx2-vx1)*(vx2-vx1);
   /* qDebug() << "R = " << MyString128::number(currentIndicatorValue[0], 'g', 10);
    qDebug() << "x1 = " << MyString128::number(x1, 'g', 10);
    qDebug() << "vx1 = " << MyString128::number(vx1, 'g', 10);
    qDebug() << "x2 = " << MyString128::number(x2, 'g', 10);
    qDebug() << "vx2 = " << MyString128::number(vx2, 'g', 10);*/

}
 void ChildSolver::indicatorFLI()
 {
     std::valarray<real_type>& x = x_current;
     real_type dev = mySqrt(x[4]*x[4]+x[5]*x[5]+x[6]*x[6]+x[7]*x[7]);
     if (myLog(dev) > currentIndicatorValue[0]) {
        currentIndicatorValue[0] = myLog(dev);
     }
 }

 void ChildSolver::indicatorOFLI()
 {
     std::valarray<real_type>& x = x_current;
     std::valarray<real_type> delta( {x[4], x[5], x[6], x[7]} );
     std::valarray<real_type> f(forceNoIndicator(real_type(0.0),x));
     real_type dev = mySqrt((project(delta,f)*project(delta,f)).sum());

     if (myLog(dev) > currentIndicatorValue[0]) {
         currentIndicatorValue[0] = myLog(dev);
     }
 }

 void ChildSolver::indicatorAPLE()
 {
     std::valarray<real_type>& x = x_current;
     real_type& t = t_current;
     real_type dev = mySqrt(x[4]*x[4]+x[5]*x[5]+x[6]*x[6]+x[7]*x[7]);
     currentIndicatorValue[0] = myLog(dev)/myLog(t);
 }

 void ChildSolver::indicatorMEGNO()
 {
     std::valarray<real_type>& x = x_current;
     real_type& t = t_current;
     real_type dev = mySqrt(x[4]*x[4]+x[5]*x[5]+x[6]*x[6]+x[7]*x[7]);

     globalMegnoParameter=myLog(dev);
     currentIndicatorValue[0] =  real_type(2.0)*(x[8]-x[9])/t;
 }

 void ChildSolver::indicatorSOMEGNO()
 {
    std::valarray<real_type>& x = x_current;
    real_type& t = t_current;
    std::valarray<real_type> delta( {x[4], x[5], x[6], x[7]} );

    std::valarray<real_type> f(forceNoIndicator(t,x));
    globalMegnoParameter=myLog(mySqrt((project(delta,f)*project(delta,f)).sum()));

    real_type result = x[8]-x[9];
    if (result > somegnoParam) somegnoParam = result;

    currentIndicatorValue[0] = 2*somegnoParam/t;
 }

 void ChildSolver::indicatorOMEGNO()
 {
    std::valarray<real_type>& x = x_current;
    real_type& t = t_current;
    std::valarray<real_type> delta( {x[4], x[5], x[6], x[7]} );

    std::valarray<real_type> f(forceNoIndicator(t,x));
    globalMegnoParameter=myLog(mySqrt((project(delta,f)*project(delta,f)).sum()));
    currentIndicatorValue[0] = real_type(2.0)*(x[8]-x[9]);
 }

 void ChildSolver::indicatorSALI()
 {
     std::valarray<real_type>& r = x_current;
     real_type dev1 = mySqrt(r[4]*r[4]+r[5]*r[5]+r[6]*r[6]+r[7]*r[7]);
     real_type dev2 = mySqrt(r[8]*r[8]+r[9]*r[9]+r[10]*r[10]+r[11]*r[11]);

     real_type d1 = mySqrt((r[8]/dev2-r[4]/dev1)*(r[8]/dev2-r[4]/dev1) + (r[9]/dev2-r[5]/dev1)*(r[9]/dev2-r[5]/dev1) + (r[10]/dev2-r[6]/dev1)*(r[10]/dev2-r[6]/dev1) + (r[11]/dev2-r[7]/dev1)*(r[11]/dev2-r[7]/dev1)   );
     real_type d2 = mySqrt((r[8]/dev2+r[4]/dev1)*(r[8]/dev2+r[4]/dev1) + (r[9]/dev2+r[5]/dev1)*(r[9]/dev2+r[5]/dev1) + (r[10]/dev2+r[6]/dev1)*(r[10]/dev2+r[6]/dev1) + (r[11]/dev2+r[7]/dev1)*(r[11]/dev2+r[7]/dev1)   );
     currentIndicatorValue[0] =  (d1<=d2) ? d1 : d2 ;
 }

 void ChildSolver::solveOneStep()
 {
    OdeSolver::solveOneStep();
    if (myFabs(h)<myFabs(1.1e-17*t_current)) {
        ode_eps*=real_type(10.0);
    }

   if (selectIndicator!=ChaoticIndicator::noIndicator) ((*this).*(pointerToIndicator))();
 /*  qDebug() << "t = " << MyString128::number(t_current, 'g', 10) << "; ";
   for (int i=0; i<dim; ++i) {
       qDebug() << "r[" << i << "] = " << MyString128::number(x_current[i], 'g', 10) << "; ";
   }
   qDebug() << getIndicatorName(selectIndicator) << " = " << MyString128::number(currentIndicatorValue[0], 'g', 10) << endl;*/
 }

 void ChildSolver::reverseSolve()
 {
    //n=-1;
    h=-h_initial;
    OdeSolver::reverseSolve();
 }
void ChildSolver::saveSolvingStep() {
    OdeSolver::saveSolvingStep();    
    indicatorValues.push_back(currentIndicatorValue);
}

 void ChildSolver::saveSolution(const QDir& dir, int prec)
 {
     std::valarray<real_type> delta({initialVector[4], initialVector[5], initialVector[6], initialVector[7]});
     std::valarray<real_type> x({initialVector[0], initialVector[1], initialVector[2], initialVector[3]});
     std::valarray<real_type> n({-x[0], -x[1], x[2]- (x[2]+mu)*(1-mu)/norm3(x-rb0) - (x[2]-1+mu)*mu/norm3(x-rb1), x[3]-x[3]*(1-mu)/norm3(x-rb0) - x[3]*mu/norm3(x-rb1)});
     std::valarray<real_type> f( forceNoIndicator(0.0,x));
     std::cout << "Вектор нормали к инт. поверхности" << double(n[0]) <<"; " << double(n[1]) << "; " << double(n[2]) << "; " << double(n[3]) << std::endl;
     std::cout << "Вектор фазового потока" << double(f[0]) <<"; " << double(f[1]) << "; " << double(f[2]) << "; " << double(f[3]) << std::endl;
     cosa = double(cosBetweenVectors(delta,f));
     cosb = double(cosBetweenVectors(delta,n));



     QString fullFileName = dir.path() + "/evol_" + MyString128::number(double(initialVector[2]),'g',prec) + "_" + MyString128::number(double(initialVector[4]),'g',prec) + ".csv";
     QFile fout(fullFileName);
     if (!fout.open(QIODevice::WriteOnly)) {
         QMessageBox msgBox;
         msgBox.setText("Ошибка! Невозможно записать файл.");
         msgBox.exec();
         return;
     }
     QTextStream outStream(&fout);

     int nColumns = initialVector.size()+3;

     outStream<<"Name;"<< MyString128::number(initialVector[2],'g',prec) << endl;
     outStream << 0 << ";";
     for(int j=0; j<dim; j++) {
         outStream<< MyString128::number(initialVector[j],'g',prec) << ";";
     }

     if (selectIndicator!=ChaoticIndicator::noIndicator) {
         for (real_type item : currentIndicatorValue) {
            ++nColumns;
            outStream << 0 << ";";
         }
     }
     outStream << MyString128::number(cosa,'g',16) << ";" << MyString128::number(cosa,'g',16) << endl;
     for (int i=0; i<solution.size(); ++i) {
         outStream << MyString128::number(timeValues[i],'g',prec) << ";";
         for(int j=0; j<dim; j++) {
             outStream<< MyString128::number(solution[i][j],'g',prec) << ";";
         }
         if (selectIndicator!=ChaoticIndicator::noIndicator) {
             for (int j=0; j<indicatorValues[i].size(); j++) {
                outStream<<MyString128::number(indicatorValues[i][j],'g',prec) << ";";
             }
         }
         outStream << MyString128::number(cosa,'g',16) << ";" << MyString128::number(cosa,'g',16) << endl;
     }
     fout.close();
 }

 real_type ChildSolver::jacobiIntegral(const std::valarray<real_type> &r)
 {
     return r[2]*r[2]+r[3]*r[3] + real_type(2.0)*(real_type(1.0)-mu)/norm(r-rb0) + real_type(2.0)*mu/norm(r-rb1) - r[0]*r[0] -r[1]*r[1];
 }

 void ChildSolver::solve()
 {
     std::cout<< "Интегрирование начато, x0 = " << MyString128::number(initialVector[2], 'g', 10).toStdString() << std::endl;
     OdeSolver::solve();
     std::cout << "Орбита х = " << MyString128::number(initialVector[2],'g', 15).toStdString() << "Относительная ошибка интеграла Якоби: " << MyString128::number(myFabs(jacobiIntegral(x_current)-C0)/C0, 'g', 10).toStdString() << std::endl;
    // reverseSolve();
     std::cout << "Орбита х = " << MyString128::number(initialVector[2],'g', 15).toStdString() << ". Интегрирование завершено за " << nIterations << " итераций" << std::endl;
    // std::cout << "Сдвиг вектора положения после интергрирования вперед-назад: " << double(mySqrt(myPow(initialVector[2]-x_current[2],real_type(2.0)) + myPow(initialVector[3] - x_current[3],real_type(2.0))) )<< std::endl;
    // std::cout << "Сдвиг вектора скорости после интергрирования вперед-назад: " << double(mySqrt(myPow(initialVector[1]-x_current[1],real_type(2.0)) + myPow(initialVector[0]-x_current[0],real_type(2.0)))) << std::endl;
 }
