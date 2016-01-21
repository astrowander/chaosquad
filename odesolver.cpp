#include "odesolver.h"
static int number_of_steps[] = { 2,4,6,8,12,16,24,32,48,64,96,128, 192, 256, 384, 512, 768, 1024};

#define ATTEMPTS sizeof(number_of_steps)/sizeof(number_of_steps[0])

inline std::valarray <real_type> OdeSolver::f(real_type t, const std::valarray<real_type> &x)
{
   //(*temp)[0]=sin(t);
   // return *temp;
    std::valarray<real_type> temp(x.size());
    temp[0]=1.5 + 1.5*cosq(4.83*t);
    return temp;
}


std::valarray<real_type> OdeSolver::euler(real_type &t, const std::valarray<real_type> &x, real_type &h)
{
        t+=h;
        return x + f(t,x)*h;
}

std::valarray<real_type> OdeSolver::rkutta4(real_type &t, const std::valarray<real_type> &x, real_type &h)
{
    int dim=x.size();
    std::valarray<real_type> k1(dim), k2(dim), k3(dim), k4(dim);
    k1=f(t,x)*h;
    t+=h/real_type(2.0);
    k2=f(t,x+k1*real_type(0.5))*h;
    k3=f(t,x+k2*real_type(0.5))*h;
    t+=h/real_type(2.0);
    k4=f(t,x+k3)*h;
    return x + (k1 + k2*real_type(2.0) + k3*real_type(2.0) + k4)/real_type(6.0);
}

std::valarray<real_type> OdeSolver::merson5(real_type &t, const std::valarray<real_type> &x, real_type &h)
{
    int dim=x.size();
    std::valarray<real_type> k1(dim), k2(dim), k3(dim), k4(dim), k5(dim), delta;
    h*=real_type(2.0);
    do
    {
        h/=real_type(2.0);
        k1=f(t,x)*h/real_type(3.0);
        k2=f(t+h/real_type(3.0),x+k1)* h/real_type(3.0);
        k3=f(t+h/real_type(3.0),x+k1/real_type(2.0)+k2/real_type(2.0))* h/real_type(3.0);
        k4=f(t+h/real_type(2.0),x+k1*real_type(0.375)+k3*real_type(1.125))* h/real_type(3.0);
        k5=f(t+h/real_type(2.0), x+k1*real_type(1.5) - k2*real_type(4.5) +k4*real_type(6.0)) *h/real_type(3.0);
        delta = k1 - k2*real_type(4.5) + k4*real_type(4.0) -k5*real_type(0.5);
    }
    while ((delta*delta).sum()>5*ode_eps);
    t+=h;
    if ((delta*delta).sum() < 5*ode_eps/32) h*=2;
    return x+(k1+k4*real_type(4.0)+k5)/real_type(2.0);
}

std::valarray<real_type> OdeSolver::butcher6(real_type &t, const std::valarray<real_type> &x, real_type &h)
{
    int dim=x.size();
    std::valarray<real_type> k1(dim), k2(dim), k3(dim), k4(dim), k5(dim), k6(dim), k7(dim);
    k1=f(t, x)*h;
    k2=f(t +h/real_type(2.0), 	x +k1*(real_type(1.0)/real_type(2.0)))*h;
    k3=f(t +real_type(2.0)*h/real_type(3.0), x +k1*(real_type(2.0)/real_type(9.0)) 	 +k2*(real_type(4.0)/real_type(9.0)))*h;
    k4=f(t +h/real_type(3.0), 	x +k1*(real_type(7.0)/real_type(36.0))	 +k2*(real_type(2.0)/real_type(9.0))   -k3*(real_type(1.0)/real_type(12.0)))*h;
    k5=f(t +real_type(5.0)*h/real_type(6.0), x -k1*(real_type(35.0)/real_type(144.0)) -k2*(real_type(55.0)/real_type(36.0)) +k3*(real_type(35.0)/real_type(48.0))  +k4*(real_type(15.0)/real_type(8.0)))*h;
    k6=f(t +h/real_type(6.0), 	x -k1*(real_type(1.0)/real_type(360.0))  -k2*(real_type(11.0)/real_type(36.0)) -k3*(real_type(1.0)/real_type(8.0))    +k4*(real_type(1.0)/real_type(2.0))    +k5*(real_type(1.0)/real_type(10.0)))*h;
    k7=f(t +h, 	x -k1*(real_type(41.0)/real_type(260.0)) +k2*(real_type(22.0)/real_type(13.0)) +k3*(real_type(43.0)/real_type(156.0)) -k4*(real_type(118.0)/real_type(39.0)) +k5*(real_type(32.0)/real_type(195.0)) +k6*(real_type(80.0)/real_type(39.0)))*h;

    t+=h;
    return x + k1*(real_type(13.0)/real_type(200.0))+k3*(real_type(11.0)/real_type(40.0))+k4*(real_type(11.0)/real_type(40.0))+k5*(real_type(4.0)/real_type(25.0))+k6*(real_type(4.0)/real_type(25.0))+k7*(real_type(13.0)/real_type(200.0));
}

std::valarray<real_type> OdeSolver::fehlberg8(real_type &t, const std::valarray<real_type> &x, real_type &h)
{
    int dim=x.size();
    std::valarray<real_type> k1(dim), k2(dim), k3(dim), k4(dim), k5(dim), k6(dim),
        k7(dim), k8(dim), k9(dim), k10(dim), k11(dim), k12(dim), k13(dim), err(dim);
     const real_type c_1_11 = real_type(41.0) / real_type(840.0);
    const real_type c6 = real_type(34.0) / real_type(105.0);
    const real_type c_7_8= real_type(9.0) / real_type(35.0);
    const real_type c_9_10 = real_type(9.0) / real_type(280.0);

    const real_type a2 = real_type(2.0) / real_type(27.0);
    const real_type a3 = real_type(1.0) / real_type(9.0);
    const real_type a4 = real_type(1.0) / real_type(6.0);
    const real_type a5 = real_type(5.0) / real_type(12.0);
    const real_type a6 = real_type(1.0) / real_type(2.0);
    const real_type a7 = real_type(5.0) / real_type(6.0);
    const real_type a8 = real_type(1.0) / real_type(6.0);
    const real_type a9 = real_type(2.0) / real_type(3.0);
    const real_type a10 = real_type(1.0) / real_type(3.0);

    const real_type b31 = real_type(1.0) / real_type(36.0);
    const real_type b32 = real_type(3.0) / real_type(36.0);
    const real_type b41 = real_type(1.0) / real_type(24.0);
    const real_type b43 = real_type(3.0) / real_type(24.0);
    const real_type b51 = real_type(20.0) / real_type(48.0);
    const real_type b53 = -real_type(75.0) / real_type(48.0);
    const real_type b54 = real_type(75.0) / real_type(48.0);
    const real_type b61 = real_type(1.0) / real_type(20.0);
    const real_type b64 = real_type(5.0) / real_type(20.0);
    const real_type b65 = real_type(4.0) / real_type(20.0);
    const real_type b71 = -real_type(25.0) / real_type(108.0);
    const real_type b74 =  real_type(125.0) / real_type(108.0);
    const real_type b75 = -real_type(260.0) / real_type(108.0);
    const real_type b76 =  real_type(250.0) / real_type(108.0);
    const real_type b81 = real_type(31.0)/real_type(300.0);
    const real_type b85 = real_type(61.0)/real_type(225.0);
    const real_type b86 = -real_type(2.0)/real_type(9.0);
    const real_type b87 = real_type(13.0)/real_type(900.0);
    const real_type b91 = real_type(2.0);
    const real_type b94 = -real_type(53.0)/real_type(6.0);
    const real_type b95 = real_type(704.0) / real_type(45.0);
    const real_type b96 = -real_type(107.0) / real_type(9.0);
    const real_type b97 = real_type(67.0) / real_type(90.0);
    const real_type b98 = real_type(3.0);
    const real_type b10_1 = -real_type(91.0) / real_type(108.0);
    const real_type b10_4 = real_type(23.0) / real_type(108.0);
    const real_type b10_5 = -real_type(976.0) / real_type(135.0);
    const real_type b10_6 = real_type(311.0) / real_type(54.0);
    const real_type b10_7 = -real_type(19.0) / real_type(60.0);
    const real_type b10_8 = real_type(17.0) / real_type(6.0);
    const real_type b10_9 = -real_type(1.0) / real_type(12.0);
    const real_type b11_1 = real_type(2383.0) / real_type(4100.0);
    const real_type b11_4 = -real_type(341.0) / real_type(164.0);
    const real_type b11_5 = real_type(4496.0) / real_type(1025.0);
    const real_type b11_6 = -real_type(301.0) / real_type(82.0);
    const real_type b11_7 = real_type(2133.0) / real_type(4100.0);
    const real_type b11_8 = real_type(45.0) / real_type(82.0);
    const real_type b11_9 = real_type(45.0) / real_type(164.0);
    const real_type b11_10 = real_type(18.0) / real_type(41.0);
    const real_type b12_1 = real_type(3.0) / real_type(205.0);
    const real_type b12_6 = - real_type(6.0) / real_type(41.0);
    const real_type b12_7 = - real_type(3.0) / real_type(205.0);
    const real_type b12_8 = - real_type(3.0) / real_type(41.0);
    const real_type b12_9 = real_type(3.0) / real_type(41.0);
    const real_type b12_10 = real_type(6.0) / real_type(41.0);
    const real_type b13_1 = -real_type(1777.0) / real_type(4100.0);
    const real_type b13_4 = -real_type(341.0) / real_type(164.0);
    const real_type b13_5 = real_type(4496.0) / real_type(1025.0);
    const real_type b13_6 = -real_type(289.0) / real_type(82.0);
    const real_type b13_7 = real_type(2193.0) / real_type(4100.0);
    const real_type b13_8 = real_type(51.0) / real_type(82.0);
    const real_type b13_9 = real_type(33.0) / real_type(164.0);
    const real_type b13_10 = real_type(12.0) / real_type(41.0);

    const real_type err_factor  = -real_type(41.0) / real_type(840.0);

    real_type module_err;
    do {
    real_type h2_7 = a2 * h;

   k1 = f(t, x);
   k2 = f(t+h2_7, x + h2_7 * k1);
   k3 = f(t+a3*h, x + h * ( b31*k1 + b32*k2) );
   k4 = f(t+a4*h, x + h * ( b41*k1 + b43*k3) );
   k5 = f(t+a5*h, x + h * ( b51*k1 + b53*k3 + b54*k4) );
   k6 = f(t+a6*h, x + h * ( b61*k1 + b64*k4 + b65*k5) );
   k7 = f(t+a7*h, x + h * ( b71*k1 + b74*k4 + b75*k5 + b76*k6) );
   k8 = f(t+a8*h, x + h * ( b81*k1 + b85*k5 + b86*k6 + b87*k7) );
   k9 = f(t+a9*h, x + h * ( b91*k1 + b94*k4 + b95*k5 + b96*k6
                                                          + b97*k7 + b98*k8) );
   k10 = f(t+a10*h, x + h * ( b10_1*k1 + b10_4*k4 + b10_5*k5 + b10_6*k6
                                          + b10_7*k7 + b10_8*k8 + b10_9*k9 ) );
   k11 = f(t+h, x + h * ( b11_1*k1 + b11_4*k4 + b11_5*k5 + b11_6*k6
                           + b11_7*k7 + b11_8*k8 + b11_9*k9 + b11_10 * k10 ) );
   k12 = f(t, x + h * ( b12_1*k1 + b12_6*k6 + b12_7*k7 + b12_8*k8
                                                 + b12_9*k9 + b12_10 * k10 ) );
   k13 = f(t+h, x + h * ( b13_1*k1 + b13_4*k4 + b13_5*k5 + b13_6*k6
                + b13_7*k7 + b13_8*k8 + b13_9*k9 + b13_10*k10 + k12 ) );

   err =  err_factor * (k1 + k11 - k12 - k13);

   module_err = mySqrt((err*err).sum())/mySqrt((x*x).sum());
   if (module_err>real_type(5.0)*ode_eps) h *=real_type(0.5)*myPow(ode_eps/module_err,real_type(0.125));
   } while (module_err>real_type(5.0)*ode_eps);
    t+=h;
    std::valarray<real_type> add(dim);
    add = h * (c_1_11 * (k1 + k11)  + c6 * k6 + c_7_8 * (k7 + k8)
                                          + c_9_10 * (k9 + k10) );
    if (module_err < real_type(5.0)*ode_eps/real_type(32.0) && !forbidEnlargingStep) h*=real_type(2.0);
    return x + add;
}


void OdeSolver::solveOneStep()
{

//    timeValues.push_back(timeValues.back());
//    real_type& t= timeValues.back();
//    std::valarray<real_type>& x = solution.back();
//    solution.push_back(((*this).*(pointerToMethod))(t,x,h));
    x_current=((*this).*(pointerToMethod))(t_current,x_current,h);
    ++size;
}
void OdeSolver::saveSolvingStep() {
   // timeValues.push_back(t_current);
    solution.push_back(x_current);
}

void OdeSolver::solve() {
   x_current=initialVector;
   t_current=tmin;
   /* while (timeValues.back() < tmax) {
       if (timeValues.back()+h>tmax) {
           h=tmax-timeValues.back();
       }
       solveOneStep();
       i++;
   }*/
   //real_type savingTimeStep = h;
   //int savingTimeMultiplier = 0, interm = 1;
   int nextTimeCount = 0;
   while (t_current <tmax ) {
       t_prev = t_current;
       x_prev = x_current;

       if (t_current>=timeValues[nextTimeCount]) {
            h = timeValues[nextTimeCount++] - t_current;
            std::cout << "t = " << MyString128::number(timeValues[nextTimeCount],'g', 10).toStdString() << std::endl;
            solveOneStep();
            saveSolvingStep();
            continue;
       }
       solveOneStep();
       ++nIterations;
   }

}
void OdeSolver::reverseSolve() {
   while (t_current > tmin) {
       if (t_current+h<tmin) {
           h=tmin-t_current;
       }
       solveOneStep();
   }
}

std::valarray<real_type> OdeSolver::gragg_bulirsch_stoer(real_type &x, const std::valarray<real_type> &y0, real_type &h)
{
      real_type step_size2[ATTEMPTS];
      std::vector< std::valarray<real_type> > tableau(ATTEMPTS+1);
      for (unsigned int i=0; i<ATTEMPTS+1; i++) {
          tableau[i].resize(dim);
      }
      real_type dum;
      real_type h_new=h;
      std::valarray<real_type> est;
      std::valarray<real_type> old_est;
      std::valarray<real_type> yscale(dim);
      yscale=y0;
      for (int i=0; i<dim; i++) {
          if (myFabs(yscale[i])<real_type(1e-33)) yscale[i]=1;
      }
      std::valarray<real_type> y1(dim);
      //int (*Extrapolate)(real_type*,real_type*,real_type*,real_type,int);

      int err;

      /* Perform the first estimate of y(x+h), store the step size which */
      /* was used squared and the estimate for subsequent called to the  */
      /* rational function approximation for the value of y(x+h) as the  */
      /* square of the step size used tends to zero.                     */
      real_type sum;
      for (int i=0; i<dim; ++i) sum += myFabs(yscale[i]);
      if (sum == real_type(0.0)) return std::valarray<real_type>(NAN,dim);

      //if (rational_extrapolate) Extrapolate = Rational_Extrapolation_to_Zero;
      //else
     // Extrapolate = Polynomial_Extrapolation_to_Zero;

      est = graggs_method(y0, x, x+h, number_of_steps[0] );
      step_size2[0] = (dum = h / (real_type) number_of_steps[0], dum * dum);
      y1 = est;

      //if (
              Polynomial_Extrapolation_to_Zero(y1, tableau, step_size2, est, 0) ;//< 0
              //) return err-1;

       /* Continue using Gragg's method with smaller step sizes, followed    */
       /* by an estimate of y(x+h) using the rational function approximation */
       /* as the number of steps becomes large or step size used in Gragg's  */
       /* method tends to zero, and then halt if the absolute difference     */
       /* between the current estimate and the previous estimate is less     */
       /* than the user-specified tolerance.  If an attempt to divide by     */
       /* zero was made in the rational function approximation, then halt    */
       /* with a return code of -2. If the procedure fails to converge       */
       /* within ATTEMPTS - 1 iterations, then halt with a return code of -1.*/
       /* If an attempt is made with a step size of 0, then halt with a      */
       /* return code of -3.                                                 */

      for (unsigned i = 1; i < ATTEMPTS; i++) {
         old_est = y1;
         est = graggs_method( y0, x, x+h, number_of_steps[i] );
         step_size2[i] = (dum = h / (real_type) number_of_steps[i], dum * dum);

         if ((err = Polynomial_Extrapolation_to_Zero(y1, tableau, step_size2, est, i) < 0)) {
             return y0;
         }
         std::valarray<real_type> temp = y1 / yscale - old_est / yscale;
        real_type module_err = mySqrt((temp*temp).sum());
        //for (int j=0; j<dim; ++j) module_err += fabsq(temp[i]);
         if ( module_err < ode_eps) {
            if (i > 1 && !forbidEnlargingStep) h_new = real_type(8.0) * h / (real_type) number_of_steps[i-1];
            else h_new = h;
            x+=h;
            h=h_new;
            return y1;
         }
         if (i>max_order) max_order=i;
      }
      //return std::valarray<real_type>(NAN,dim);
      x+=h;
      if (!forbidEnlargingStep) {
          h=h_new;
      }
      return y1;
}

std::valarray<real_type> OdeSolver::graggs_method( const std::valarray<real_type> &y_old, real_type x0,real_type x, int number_of_steps )
{
    real_type h = (x - x0) / (real_type) number_of_steps;
    real_type h2 =  h + h;
    std::valarray<real_type> y0(y_old);
    std::valarray<real_type> y1 = y0 + h * f(x0,y0);
    std::valarray<real_type> y2;

    while ( --number_of_steps ) {
          x0 += h;
          y2 = y0 + h2 * f(x0,y1);
          y0 = y1;
          y1 = y2;
    }
    return  real_type(0.5) * ( y0 + y1 + h * f(x,y1) );
}

int OdeSolver::Polynomial_Extrapolation_to_Zero( std::valarray<real_type> &fzero, std::vector<std::valarray<real_type> > &tableau,
                                               real_type x[], const std::valarray<real_type> &f, int n ) {

   std::valarray<real_type> back_two_columns(dim);    //  T[row,col-2];
   std::valarray<real_type> old_aux(dim);             //  T[row-1,col];
              //  T[row,col];
   std::valarray<real_type>  vertical_diff(dim);       //  T[row,col]-T[row-1,col]
   std::valarray<real_type>  backslant_diff(dim);      //  T[row,col]-T[row,col-1]
   std::valarray<real_type>  forwardslant_diff(dim);   //  T[row,col]-T[row-1,col-1];
   std::valarray<real_type> denominator(dim);
   int i;

   if (n == 0) { tableau[0] = f; return 0; }
   if ( x[n] == real_type(0.0) ) { fzero = f; return -2; }

   //back_two_columns = real_type(0.0);
   old_aux = tableau[0];
   tableau[0] = f;
   for (i = 0; i < n; i++) {
      vertical_diff = tableau[i] - old_aux;
      backslant_diff = tableau[i] - back_two_columns;
      forwardslant_diff = backslant_diff - vertical_diff;
      denominator = (x[n-i-1]/x[n]) * forwardslant_diff - backslant_diff;
      //if (abs(denominator).sum() == real_type(0.0)) return -1;
      back_two_columns = old_aux;
      old_aux = tableau[i+1];
    /*  for (int j=0; j<dim; j++) {
        if (fabsq(denominator[j]) >1e-14) {
            tableau[i+1][j] = tableau[i][j] + vertical_diff[j] * backslant_diff[j] / denominator[j];
        }
      }*/
      tableau[i+1] = tableau[i];//+ vertical_diff* backslant_diff / denominator;
      for (int j=0; j<dim; j++) {
              if (myFabs(denominator[j])!=0) {
                  tableau[i+1][j] += vertical_diff[j] * backslant_diff[j] / denominator[j];
              }
            }
  }


   fzero = tableau[n];
   return 0;
}
