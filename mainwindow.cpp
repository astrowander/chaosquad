#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <fstream>
#include <QTime>
MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    //solver=nullptr;
    //plot = nullptr;
    ui->progressBar->setValue(0);
    fileName.clear();
    workingDirectory.setPath("/home/astrowander/Qt/Projects/chaosDetector2/outputs");
}

MainWindow::~MainWindow()
{
    //plots.clear();
    delete ui;
}

void MainWindow::on_pushButton_clicked()
{
    real_type xmin, xmax, dx, tmin, tmax, h, eps, C0=4.0l;
    QTime timer;
    timer.start();
    MethodOfIntegration method;

    initialConditions.clear();
    ui->progressBar->setValue(0);
    //TODO: переименовать текстовые строки
    xmin = static_cast<real_type> (ui->lineEdit_xmin->text().toDouble());
    xmax = static_cast<real_type> (ui->lineEdit_xmax->text().toDouble());
    dx = static_cast<real_type> (ui->lineEdit_dx->text().toDouble());
    tmin = static_cast<real_type> (ui->lineEdit_tmin->text().toDouble());
    tmax = static_cast<real_type> (ui->lineEdit_tmax->text().toDouble());
    h = static_cast<real_type> (ui->lineEdit_h->text().toDouble());
    eps=static_cast<real_type> (ui->lineEdit_eps->text().toDouble());

    //QDateTime now = QDateTime::currentDateTime();
    //QString nowString = QDateTime::currentDateTime().toString("yyyy-MM-dd hh:mm:ss");
    std::cout << double(tmax) << " " << double(h)  << std::endl << std::endl;



    if (ui->radioButton->isChecked()) {
        method=MethodOfIntegration::euler;
    }
    else if (ui->radioButton_2->isChecked()) {
        method=MethodOfIntegration::rkutta4;
    }
    else if (ui->radioButton_3->isChecked()) {
        method=MethodOfIntegration::merson5;
    }
    else if (ui->radioButton_4->isChecked()) {
        method=MethodOfIntegration::butcher6;
    }
    else if (ui->radioButton_5->isChecked()) {
        method=MethodOfIntegration::fehlberg8;
    }
    else if (ui->radioButton_13->isChecked()) {
        method=MethodOfIntegration::gragg_bulirsch_stoer;
    }
    else {
        return;
    }

    int n = ui->listWidget->count(), nRepeats = ui->spinBox->value();
    for (int i=0; i< n; ++i){
        initialConditions.append(ui->listWidget->item(i)->text().toDouble());
    }
    real_type x0=xmin;
    if (dx!=0)
    while (x0 < xmax + dx / 2) {
        initialConditions.append(x0);
        x0+=dx;
    }
    if (initialConditions.isEmpty()) return;
    ChaoticIndicator chaoticIndicator;
    if (ui->radioButton_6->isChecked()) {
       chaoticIndicator=ChaoticIndicator::noIndicator;

    }
    else if (ui->radioButton_7->isChecked()) {
       chaoticIndicator=ChaoticIndicator::FLI;

    }
    else if (ui->radioButton_8->isChecked()) {
       chaoticIndicator=ChaoticIndicator::OFLI;

    }
    else if (ui->radioButton_9->isChecked()) {
        chaoticIndicator=ChaoticIndicator::APLE;

    }
    else if (ui->radioButton_10->isChecked()) {
       chaoticIndicator=ChaoticIndicator::MEGNO;

    }
    else if (ui->radioButton_11->isChecked()) {
       chaoticIndicator=ChaoticIndicator::OMEGNO;

    }
    else if (ui->radioButton_12->isChecked()) {
       chaoticIndicator=ChaoticIndicator::SALI;

    }
    else if (ui->radioButton_16->isChecked()) {
        chaoticIndicator=ChaoticIndicator::allIndicators;
    }
    else if (ui->radioButton_16->isChecked()) {
        chaoticIndicator=ChaoticIndicator::SOMEGNO;
    }
    else if (ui->radioButton_indicator_R->isChecked()) {
        chaoticIndicator=ChaoticIndicator::R;
    }

    int N = initialConditions.size()*nRepeats;
    for (int i=0; i<N; ++i) {
        Task task(method,tmin,tmax,h,eps, ui->lineEdit_mu->text().toDouble(), partition);
        x0=initialConditions[i/nRepeats];
        task.selectIndicator=chaoticIndicator;

        if (ui->radioButton_14->isChecked()) {
            if (ui->lineEdit_ncounts->text()=="") return;
            task.partition = uniformPartition;
            task.nCounts = ui->lineEdit_ncounts->text().toInt();
        }
        else if (ui->radioButton_15->isChecked()) {
            task.partition = logarithmicPartition;
        }
        else if (ui->radioButton_17->isChecked()) {
            task.partition = noTimePartition;
        }
        else if (ui->radioButton_19->isChecked()) {
            task.partition = poincarePartition;
        }
        task.buildTimePartition();
        task.getInitialVector(x0, C0, ui->randomDeviationCheckBox->isChecked());
        //task.getInitialVectorHH(1.0/6.0,x0, 0.0,ui->randomDeviationCheckBox->isChecked());
        tasksQueue.enqueue(task);
    }

    workingDirectory.setPath("/home/astrowander/Qt/Projects/chaosDetector2/outputs");

    QStringList nowStrings = QDateTime::currentDateTime().toString("yyyy-MM-dd hh:mm:ss").split(" ");
    if ( !workingDirectory.cd(nowStrings[0]) ) {
        workingDirectory.mkdir(nowStrings[0]);
        workingDirectory.cd(nowStrings[0]);
    }
    QString subdir = getIndicatorName(chaoticIndicator) + "_" + QString::number(N) + "_" + nowStrings[1];
    workingDirectory.mkdir(subdir);
    workingDirectory.cd(subdir);
    QString mainLogName = workingDirectory.path() + "/main_log.csv";
    QFile mainLogFile(mainLogName);

    if (!mainLogFile.open(QIODevice::WriteOnly)) {
        QMessageBox msgBox;
        msgBox.setText("Ошибка! Невозможно записать файл.");
        msgBox.exec();
        return;
    }

    QTextStream mainLogStream(&mainLogFile);

    mainLogStream<<"Title;" << getIndicatorName(tasksQueue.back().selectIndicator) << " map" << endl;
    mainLogStream<< "Number of data sets;" << tasksQueue.size() << endl;
    mainLogStream << "Names of columns" << endl << "t;";
    for(int j=0; j<tasksQueue.back().initialVector.size(); j++) {
        mainLogStream<< "r[" << j << "];";
    }

    if (tasksQueue.back().selectIndicator == ChaoticIndicator::allIndicators) {
        mainLogStream << getIndicatorName(ChaoticIndicator::FLI)    << ";"  << getIndicatorName(ChaoticIndicator::OFLI) << ";"
                      << getIndicatorName(ChaoticIndicator::APLE)   << ";"  << getIndicatorName(ChaoticIndicator::MEGNO) << ";"
                      << getIndicatorName(ChaoticIndicator::OMEGNO) << ";"  << getIndicatorName(ChaoticIndicator::SALI)
                      <<  "; cos(α); cos(β) ";
    }
    else if (tasksQueue.back().selectIndicator != ChaoticIndicator::noIndicator) {
        mainLogStream << getIndicatorName(tasksQueue.back().selectIndicator) << "; cos(α); cos(β) ";
    }

    mainLogStream << endl << ";" << endl;

    QFile thresholdsLogFile(workingDirectory.path() + "/thresholds_log.csv");
    if (!thresholdsLogFile.open(QIODevice::WriteOnly)) {
        QMessageBox msgBox;
        msgBox.setText("Ошибка! Невозможно записать файл thresholds_log.csv.");
        msgBox.exec();
        return;
    }
    QTextStream thresholdsLogStream(&thresholdsLogFile);
    thresholdsLogStream << "Орбита\Порог;" << "cos(a);" << "cos(b)";
    for (int j=0; j<nThresholds; ++j) {
        thresholdsLogStream << ";" << QString::number(thresholdValues[j],'g',10);
    }
    thresholdsLogStream << endl;
    
    ui->progressBar->setMaximum(tasksQueue.size());
    std::cout<<"Queue of tasks: " << std::endl;
    do {
        if (solvers.size()<maxThreads && (!tasksQueue.isEmpty())) {
            solvers.append(new MultiThreadSolver(tasksQueue.dequeue()));
            solvers.last()->setNReturns(ui->lineEdit_nreturns->text().toInt());
            solvers.last()->start();            
        }
        for (int i=0; i<solvers.size(); i++) {
            if (solvers[i]->isOver()) {                
                    ui->progressBar->setValue(ui->progressBar->value()+1);
                    solvers[i]->saveSolution(workingDirectory, thresholdsLogStream);
                    solvers[i]->quit();
                    solvers[i]->wait();
                    delete solvers[i];
                    solvers.erase(solvers.begin()+i);                
            }
        }        
    } while (!solvers.isEmpty());
   // mainLogFile.close();

    QStringList filesList = workingDirectory.entryList(QDir::Files);
    
    for (QString entryName : filesList) {
        if (entryName=="main_log.csv"|| entryName=="thresholds_log.csv") continue;
        QFile tempFile(workingDirectory.path()+"/"+entryName);
        if(!tempFile.open(QIODevice::ReadOnly)) continue;
        QTextStream tempStream(&tempFile);
        QString line;
        while(!tempStream.atEnd()) {
            line = tempStream.readLine();
            mainLogStream << line <<endl;
        }
        tempFile.remove();
    }    
    mainLogFile.close();
    thresholdsLogFile.close();
    std::cout<< "Расчет завершен  Прошло " << timer.elapsed() / 1000 << "c" << " " << timer.elapsed() % 1000 << "мс" << std::endl;
    //workingDirectory.cd("../..");
}

void MainWindow::on_pushButton_2_clicked()
{




   real_type tmin = static_cast<real_type> (ui->lineEdit_tmin->text().toDouble());
   real_type tmax = static_cast<real_type> (ui->lineEdit_tmax->text().toDouble());
   real_type h = static_cast<real_type> (ui->lineEdit_h->text().toDouble());
   real_type eps=static_cast<real_type> (ui->lineEdit_eps->text().toDouble());
   real_type C0 = 4;

   //QDateTime now = QDateTime::currentDateTime();
   //QString nowString = QDateTime::currentDateTime().toString("yyyy-MM-dd hh:mm:ss");

   ChaoticIndicator chaoticIndicator;
   MethodOfIntegration method;
   if (ui->radioButton_6->isChecked()) {
      chaoticIndicator=ChaoticIndicator::noIndicator;

   }
   else if (ui->radioButton_7->isChecked()) {
      chaoticIndicator=ChaoticIndicator::FLI;

   }
   else if (ui->radioButton_8->isChecked()) {
      chaoticIndicator=ChaoticIndicator::OFLI;

   }
   else if (ui->radioButton_9->isChecked()) {
       chaoticIndicator=ChaoticIndicator::APLE;

   }
   else if (ui->radioButton_10->isChecked()) {
      chaoticIndicator=ChaoticIndicator::MEGNO;

   }
   else if (ui->radioButton_11->isChecked()) {
      chaoticIndicator=ChaoticIndicator::OMEGNO;

   }
   else if (ui->radioButton_12->isChecked()) {
      chaoticIndicator=ChaoticIndicator::SALI;

   }
   else if (ui->radioButton_16->isChecked()) {
       chaoticIndicator=ChaoticIndicator::allIndicators;
   }
   else if (ui->radioButton_indicator_R->isChecked()) {
       chaoticIndicator=ChaoticIndicator::R;
   }


   if (ui->radioButton->isChecked()) {
       method=MethodOfIntegration::euler;
   }
   else if (ui->radioButton_2->isChecked()) {
       method=MethodOfIntegration::rkutta4;
   }
   else if (ui->radioButton_3->isChecked()) {
       method=MethodOfIntegration::merson5;
   }
   else if (ui->radioButton_4->isChecked()) {
       method=MethodOfIntegration::butcher6;
   }
   else if (ui->radioButton_5->isChecked()) {
       method=MethodOfIntegration::fehlberg8;
   }
   else if (ui->radioButton_13->isChecked()) {
       method=MethodOfIntegration::gragg_bulirsch_stoer;
   }
   else {
       return;


   }
/*
   Task tasks[3] { Task (method,tmin,tmax,h,eps, ui->lineEdit_mu->text().toDouble(), partition) ,
                   Task (method,tmin,tmax,h,eps, ui->lineEdit_mu->text().toDouble(), partition) ,
                   Task (method,tmin,tmax,h,eps, ui->lineEdit_mu->text().toDouble(), partition)  };

   real_type x[4] = {ui->lineEdit_xmin->text().toDouble(), 0, ui->lineEdit_xmax->text().toDouble(), 0};
   x[1] = (x[0]+x[2])/2;

   real_type indValues[4];
   real_type tol = 5e-9, x_pred=1e6;

   for (int i=0; i<3; ++i) {
       tasks[i].selectIndicator=chaoticIndicator;
       tasks[i].partition = noTimePartition;
       tasks[i].buildTimePartition();
       tasks[i].getInitialVector(x[i], C0, false);
       solvers.append(new MultiThreadSolver(tasks[i], i));
       solvers.last()->start();
   }

   do {
       for (int i=0; i<solvers.size(); ++i) {
        if (solvers[i]->isOver()) {
           indValues[solvers[i]->getID()] = solvers[i]->getIndicatorValue(0).last();
           std::cout << "x = " << QString::number(double(solvers[i]-> getInitialVector()[2]),'g', 10).toStdString() << "; " << getIndicatorName(solvers[i]->getIndicator()).toStdString() << " = " << QString::number(double(solvers[i]->getIndicatorValue(0).last()),'g', 10).toStdString() << std::endl << std::endl;
           solvers[i]->quit();
           solvers[i]->wait();
           delete solvers[i];
           solvers.erase(solvers.begin()+i);
        }
       }
   } while (!solvers.isEmpty());

   int nIterations=0;
   while (true) {
       ++nIterations;
       real_type a = indValues[0]*(x[2]-x[1]) - indValues[1]*(x[2]-x[0]) + indValues[2]*(x[1]-x[0]);
       real_type b = -indValues[0]*(x[2]*x[2]-x[1]*x[1]) + indValues[1]*(x[2]*x[2]-x[0]*x[0]) - indValues[2]*(x[1]*x[1]-x[0]*x[0]);
       x[3] = -b/(2.0l*a);
       if (fabsq(x[3]-x[1])<tol) break;


       tasks[0].getInitialVector(x[3],C0, false);
       solvers.append(new MultiThreadSolver(tasks[0])
               );
       solvers[0]->start();

       do {
           if (solvers[0]->isOver()) {
               indValues[3] = solvers[0]->getIndicatorValue(0).last();
               std::cout << "x = " << QString::number(double(solvers[0]-> getInitialVector()[2]),'g', 10).toStdString() << "; " << getIndicatorName(solvers[0]->getIndicator()).toStdString() << " = " << QString::number(double(solvers[0]->getIndicatorValue(0).last()),'g', 10).toStdString() << std::endl << std::endl;
               solvers[0]->quit();
               solvers[0]->wait();
               delete solvers[0];
               solvers.clear();
           }
       } while (!solvers.isEmpty());

       if (x[3] > x[1]) {
           x[0] = x[1];
           indValues[0] = indValues[1];
       }
       else {
           x[2] = x[1];
           indValues[2] = indValues[1];
       }
       x[1] = x[3];
       indValues[1] = indValues[3];
/*
       if (max_i!=0) {
        std::swap(x[max_i],x[0]);
        std::swap(indValues[max_i], indValues[0]);
        if (min_i==0) min_i=max_i;
       }

       if (min_i!=1) {
            std::swap(x[min_i],x[1]);
            std::swap(indValues[min_i], indValues[1]);
       }
*/
  // }

/*
   Task tasks[2] { Task (method,tmin,tmax,h,eps, ui->lineEdit_mu->text().toDouble(), partition) ,
                   Task (method,tmin,tmax,h,eps, ui->lineEdit_mu->text().toDouble(), partition)  };
   real_type x[4] = {ui->lineEdit_xmin->text().toDouble(), 0, 0, ui->lineEdit_xmax->text().toDouble()};
   real_type phi = (real_type(1.0) + mySqrt(real_type(5.0)))/real_type(2.0);
   x[1] = x[3] - (x[3]-x[0])/phi;
   x[2] = x[0] + (x[3]-x[0])/phi;
   real_type indValues[4];
   real_type tol=ui->lineEdit_dx->text().toDouble();

   for (int i=0; i<2; i++) {
        tasks[i].selectIndicator=chaoticIndicator;
        tasks[i].partition = noTimePartition;
        tasks[i].buildTimePartition();
        tasks[i].getInitialVector(x[i+1], C0, false);
        solvers.append(new MultiThreadSolver(tasks[i],i));
        solvers.last()->start();
   }

   do {
       for (int i=0; i<solvers.size(); ++i) {
        if (solvers[i]->isOver()) {
           indValues[solvers[i]->getID()+1] = solvers[i]->getIndicatorValue(solvers[i]->getSolutionSize()-1).at(0);

           std::cout << "solution back element number = " << solvers[i]->getSolutionSize()-1 << "; " << "solution back element" <<  MyString128::number(indValues[solvers[i]->getID()+1],'g', 10).toStdString() << std::endl << "x = " << MyString128::number(double(solvers[i]-> getInitialVector()[2]),'g', 10).toStdString() << "; " << getIndicatorName(solvers[i]->getIndicator()).toStdString() << " = " << MyString128::number(double(solvers[i]->getIndicatorValue(solvers[i]->getSolutionSize()-1).at(0)),'g', 10).toStdString() << std::endl << std::endl;
           solvers[i]->quit();
           solvers[i]->wait();
           delete solvers[i];
           solvers.erase(solvers.begin()+i);
        }
       }
   } while (!solvers.isEmpty());

   int nIterations=0;
   bool x1IsComputing = false;

   if (indValues[1]<indValues[2]) {
       x[3] = x[2];
       indValues[3] = indValues[2];

       x[2] = x[1];
       indValues[2] = indValues[1];

       x[1] = x[1] = x[3] - (x[3]-x[0])/phi;
       x1IsComputing = true;
   }
   else {
       x[0] = x[1];
       indValues[0] = indValues[1];
       x[1] = x[2];
       indValues[1] = indValues[2];

       x[2] = x[0] + (x[3]-x[0])/phi;
   }

   while (fabsq(x[1]-x[2]) > tol) {
       ++nIterations;

       //x[2] = x[1] + (x[1]-x[0])*(indValues[1]-indValues[0]) / (2* (2* indValues[1] - indValues[0]));

       if (x1IsComputing) {
        tasks[0].getInitialVector(x[1],C0, false);
       }
       else {
        tasks[0].getInitialVector(x[2],C0, false);
       }
       solvers.append(new MultiThreadSolver(tasks[0]));
       solvers[0]->start();

       do {
           if (solvers[0]->isOver()) {
               if (x1IsComputing) {
                indValues[1] = solvers[0]->getIndicatorValue(solvers[0]->getSolutionSize()-1).at(0);
               }
               else {
                indValues[2] = solvers[0]->getIndicatorValue(solvers[0]->getSolutionSize()-1).at(0);
               }

               std::cout << "x = " << QString::number(double(solvers[0]-> getInitialVector()[2]),'g', 10).toStdString() << "; " << getIndicatorName(solvers[0]->getIndicator()).toStdString() << " = " << QString::number(double(solvers[0]->getIndicatorValue(solvers[0]->getSolutionSize()-1).at(0)),'g', 10).toStdString() << std::endl << std::endl;
               solvers[0]->quit();
               solvers[0]->wait();
               delete solvers[0];
               solvers.clear();
           }          
       } while (!solvers.isEmpty());

       if (indValues[1]<indValues[2]) {
           x[3] = x[2];
           indValues[3] = indValues[2];

           x[2] = x[1];
           indValues[2] = indValues[1];

           x[1] = x[1] = x[3] - (x[3]-x[0])/phi;
           x1IsComputing = true;
       }
       else {
           x[0] = x[1];
           indValues[0] = indValues[1];
           x[1] = x[2];
           indValues[1] = indValues[2];

           x[2] = x[0] + (x[3]-x[0])/phi;
       }

   }
*/
/*
   Task tasks[8] { Task (method,tmin,tmax,h,eps, ui->lineEdit_mu->text().toDouble(), partition) ,
                   Task (method,tmin,tmax,h,eps, ui->lineEdit_mu->text().toDouble(), partition) ,
                   Task (method,tmin,tmax,h,eps, ui->lineEdit_mu->text().toDouble(), partition) ,
                   Task (method,tmin,tmax,h,eps, ui->lineEdit_mu->text().toDouble(), partition) ,
                   Task (method,tmin,tmax,h,eps, ui->lineEdit_mu->text().toDouble(), partition) ,
                   Task (method,tmin,tmax,h,eps, ui->lineEdit_mu->text().toDouble(), partition) ,
                   Task (method,tmin,tmax,h,eps, ui->lineEdit_mu->text().toDouble(), partition) ,
                   Task (method,tmin,tmax,h,eps, ui->lineEdit_mu->text().toDouble(), partition)
                 };
   real_type x[10] = {ui->lineEdit_xmin->text().toDouble(), 0, 0, 0, 0, 0, 0, 0, 0, ui->lineEdit_xmax->text().toDouble()};
   real_type indValues[10];
   real_type tol=ui->lineEdit_dx->text().toDouble();
   int minPos, nIterations=0;

   while (fabsq(x[9]-x[0]) > tol) {
       ++nIterations;

       for (int i=1; i<9; ++i) {
           x[i] = x[0] + i*(x[9]-x[0])/real_type(9.0);
       }

       for (int i=1; i<9; i++) {
            tasks[i-1].selectIndicator=chaoticIndicator;
            tasks[i-1].partition = noTimePartition;
            tasks[i-1].buildTimePartition();
            tasks[i-1].getInitialVector(x[i+1], C0, false);
            solvers.append(new MultiThreadSolver(tasks[i-1],i));
            solvers.last()->start();
       }

       do {
           for (int i=0; i<solvers.size(); ++i) {
            if (solvers[i]->isOver()) {
               int xpos = solvers[i]->getID();
               indValues[xpos] = solvers[i]->getIndicatorValue(solvers[i]->getSolutionSize()-1).at(0);

               //std::cout << "solution back element number = " << solvers[i]->getSolutionSize()-1 << "; " << "solution back element" <<  MyString128::number(indValues[solvers[i]->getID()+1],'g', 10).toStdString() << std::endl << "x = " << MyString128::number(double(solvers[i]-> getInitialVector()[2]),'g', 10).toStdString() << "; " << getIndicatorName(solvers[i]->getIndicator()).toStdString() << " = " << MyString128::number(double(solvers[i]->getIndicatorValue(solvers[i]->getSolutionSize()-1).at(0)),'g', 10).toStdString() << std::endl << std::endl;
               solvers[i]->quit();
               solvers[i]->wait();
               delete solvers[i];
               solvers.erase(solvers.begin()+i);
            }
           }
       } while (!solvers.isEmpty());


       real_type minValue = 1e100;
       for (int i=1; i<9; ++i) {
           if (indValues[i]<minValue) {
               minValue = indValues[i];
               minPos = i;
           }
       }

       x[0] = x[minPos-1];
       x[9] = x[minPos+1];
       qDebug() << "x[0] = " << MyString128::number(x[0], 'g', 15);
       qDebug() << "x[9] = " << MyString128::number(x[9], 'g', 15);
   }

  // real_type left, right, middleLeft, middleRight;
  */
/*

   Task tasks[2] { Task (method,tmin,tmax,h,eps, ui->lineEdit_mu->text().toDouble(), partition) ,
                   Task (method,tmin,tmax,h,eps, ui->lineEdit_mu->text().toDouble(), partition)  };

   const real_type phi = (1.0l + sqrt(5.0l))/2.0l;
   real_type x[4] = {ui->lineEdit_xmin->text().toDouble(), 0, 0, ui->lineEdit_xmax->text().toDouble()};
   x[1] = x[3] - (x[3]-x[0])/phi;
   x[2] = x[0] +(x[3]-x[0])/phi;

   real_type indValues[4];
   real_type tol = 5e-9;

   for (int i=0; i<2; i++) {
        tasks[i].selectIndicator=chaoticIndicator;
        tasks[i].partition = noTimePartition;
        tasks[i].buildTimePartition();
        tasks[i].getInitialVector(x[i+1], C0, false);
        solvers.append(new MultiThreadSolver(tasks[i]));
        solvers.last()->start();
   }

   do {
       for (int i=0; i<solvers.size(); ++i) {
        if (solvers[i]->isOver()) {
           indValues[i+1] = solvers[i]->getIndicatorValue(0).last();
           std::cout << "x = " << QString::number(double(solvers[i]-> getInitialVector()[2]),'g', 10).toStdString() << "; " << getIndicatorName(solvers[i]->getIndicator()).toStdString() << " = " << QString::number(double(solvers[i]->getIndicatorValue(0).last()),'g', 10).toStdString() << std::endl << std::endl;
           solvers[i]->quit();
           solvers[i]->wait();
           delete solvers[i];
           solvers.erase(solvers.begin()+i);
        }
       }
   } while (!solvers.isEmpty());
    int nIterations=0;
   while (fabsq(x[3]-x[0]) > tol) {
       //Task task1;
       ++nIterations;
       bool x2IsComputing;
       if (indValues[1] >=indValues[2]) {

           x[0] = x[1];
           indValues[0] = indValues[1];
           x[1] = x[2];
           indValues[1] = indValues[2];

           x[2] = x[0] +(x[3]-x[0])/phi;
           x2IsComputing = true;
           tasks[0].getInitialVector(x[2], C0, false);
       }
       else {
           x[3] =x[2];
           indValues[3] = indValues[2];

           x[2] = x[1];
           indValues[2] = indValues[1];

           x[1] = x[3] - (x[3]-x[0])/phi;
           x2IsComputing = false;
           tasks[0].getInitialVector(x[1], C0, false);
       }

       solvers.append(new MultiThreadSolver(tasks[0]));
       solvers[0]->start();

       do {
           if (solvers[0]->isOver()) {
               if (!x2IsComputing) {
                    indValues[1] = solvers[0]->getIndicatorValue(0).last();
               }
               else {
                   indValues[2] = solvers[0]->getIndicatorValue(0).last();
               }
               std::cout << "x = " << QString::number(double(solvers[0]-> getInitialVector()[2]),'g', 10).toStdString() << "; " << getIndicatorName(solvers[0]->getIndicator()).toStdString() << " = " << QString::number(double(solvers[0]->getIndicatorValue(0).last()),'g', 10).toStdString() << std::endl << std::endl;
               solvers[0]->quit();
               solvers[0]->wait();
               delete solvers[0];
               solvers.clear();
           }
       } while (!solvers.isEmpty());
 }
 */
/*
    real_type x, x_pred=0, dx=ui->lineEdit_dx->text().toDouble(), tol = 5*dx;
    QVector<real_type> firstDerivativeValues, dxValues;
    real_type vals[5];
    x = ui->lineEdit_xmin->text().toDouble();
    int nIterations=0;

    while (true) {
        nIterations++;
        real_type args[5] {x-2*dx, x - dx, x, x + dx,x + 2*dx};

        Task* tasks[5];

        for (int i=0; i<5; ++i) {
            tasks[i] = new Task(method,tmin,tmax,h,eps, ui->lineEdit_mu->text().toDouble(), partition);
            tasks[i]->selectIndicator=chaoticIndicator;
            tasks[i]->partition = noTimePartition;
            tasks[i]->buildTimePartition();
            tasks[i]->getInitialVector(args[i], C0, false);
            solvers.append(new MultiThreadSolver(*tasks[i], i));
            solvers.last()->start();
        }

        do {
            for (int i=0; i<solvers.size(); ++i) {
             if (solvers[i]->isOver()) {
                vals[solvers[i]->getID()] = solvers[i]->getIndicatorValue(0).last();
                std::cout << "x = " << QString::number(double(solvers[i]-> getInitialVector()[2]),'g', 15).toStdString() << "; " << getIndicatorName(solvers[i]->getIndicator()).toStdString() << "R =" << " = " << MyString128::number(vals[solvers[i]->getID()],'g', 15).toStdString() << std::endl << std::endl;
                solvers[i]->quit();
                solvers[i]->wait();
                delete solvers[i];
                solvers.erase(solvers.begin()+i);
             }
            }
        } while (!solvers.isEmpty());

        real_type firstDerivative = (-vals[4] + 8*vals[3] - 8*vals[1] + vals[0])/(12*dx);
        firstDerivativeValues.append(firstDerivative);
        dxValues.append(dx);
        //real_type secondDerivative = (-vals[4] + 16*vals[3] - 30 * vals[2] + 16*vals[1] -vals[0])/(12*dx*dx);
        qDebug()<<"dx = " << MyString128::number(dx,'g',2);
        qDebug() << "f'(x) = " << MyString128::number(firstDerivative,'g',10);
        int fdvSize = firstDerivativeValues.size();
        if (fdvSize>=2) qDebug()<< "|D'(n) - D'(n-1)| = " <<  MyString128::number(fabsq(firstDerivativeValues[fdvSize-2]-firstDerivativeValues[fdvSize-1]),'g',2);
        if (fdvSize>=3) qDebug()<< "|D'(n) - D'(n-1)| = " <<  MyString128::number(fabsq(firstDerivativeValues[fdvSize-3]-firstDerivativeValues[fdvSize-2]),'g',2);

        if(fdvSize < 3 || fabsq(firstDerivativeValues[fdvSize-1]-firstDerivativeValues[fdvSize-2]) < fabsq(firstDerivativeValues[fdvSize-2]-firstDerivativeValues[fdvSize-3])) {
            dx/=10;
            continue;
        }
        break;
        }
       qDebug() << "Оптимальный шаг = " << MyString128::number(dx, 'g', 2);
*/
   real_type dx=1e-9, tol =real_type(ui->lineEdit_dx->text().toDouble());
   QVector<real_type> xValues;
   real_type vals[5];
   xValues.append( real_type(ui->lineEdit_xmin->text().toDouble()));
   int nIterations=0;

   while (true) {
       nIterations++;
       real_type x = xValues.last();
       real_type args[5] {x-2*dx, x - dx, x, x + dx,x + 2*dx};

       Task* tasks[5];

       for (int i=0; i<5; ++i) {
           tasks[i] = new Task(method,tmin,tmax,h,eps, real_type(ui->lineEdit_mu->text().toDouble()), partition);
           tasks[i]->selectIndicator=chaoticIndicator;
           tasks[i]->partition = logarithmicPartition;
           tasks[i]->buildTimePartition();
           tasks[i]->getInitialVector(args[i], C0, false);
           solvers.append(new MultiThreadSolver(*tasks[i], i));
           solvers.last()->setNReturns(ui->lineEdit_nreturns->text().toInt());
           solvers.last()->start();
       }

       do {
           for (int i=0; i<solvers.size(); ++i) {
            if (solvers[i]->isOver()) {
               vals[solvers[i]->getID()] = solvers[i]->getIndicatorValue(0).last();
               std::cout << "x = " << QString::number(double(solvers[i]-> getInitialVector()[2]),'g', 15).toStdString() << "; " << getIndicatorName(solvers[i]->getIndicator()).toStdString() << "R =" << " = " << MyString128::number(vals[solvers[i]->getID()],'g', 15).toStdString() << std::endl << std::endl;
               solvers[i]->quit();
               solvers[i]->wait();
               delete solvers[i];
               solvers.erase(solvers.begin()+i);
            }
           }
       } while (!solvers.isEmpty());

       real_type firstDerivative = (-vals[4] + 8*vals[3] - 8*vals[1] + vals[0])/(12*dx);
       //firstDerivativeValues.append(firstDerivative);
       //dxValues.append(dx);
       xValues.append(xValues.last()-vals[2]/firstDerivative);
       //real_type secondDerivative = (-vals[4] + 16*vals[3] - 30 * vals[2] + 16*vals[1] -vals[0])/(12*dx*dx);
       //qDebug()<<"dx = " << MyString128::number(dx,'g',2);
       qDebug() << "f'(x) = " << MyString128::number(firstDerivative,'g',10);
       qDebug() << "x[n] = " << MyString128::number(xValues.last(),'g',10);
       //int xSize = firstDerivativeValues.size();
       if (nIterations>=1) qDebug()<< "|x(n) - x(n-1) = " <<  MyString128::number(myFabs(xValues[nIterations-1]-xValues[nIterations]),'g',10);
       if (nIterations>=2) qDebug()<< "|x(n-1) - x(n-2)| = " <<  MyString128::number(myFabs(xValues[nIterations-2]-xValues[nIterations-1]),'g',10);

       if(myFabs(xValues[nIterations]-xValues[nIterations-1]) < tol || (nIterations >= 3 && myFabs(xValues[nIterations]-xValues[nIterations-1]) > myFabs(xValues[nIterations-1]-xValues[nIterations-2]))) {
           break;
       }

       }
      qDebug() << "Периодическая орбита " << MyString128::number(xValues.last(), 'g', 15);
      qDebug() <<"Погрешность локализации = " << MyString128::number(myFabs(xValues[nIterations]-xValues[nIterations-1]), 'g', 15);
   /*
   real_type x, x_pred=0, tol = ui->lineEdit_dx->text().toDouble();
   std::valarray<real_type> solution;
   real_type indValue;
   x = ui->lineEdit_xmin->text().toDouble();
   int nIterations=0;

   while (fabsq(x-x_pred)>tol) {
       nIterations++;
       Task* task;
       task = new Task(method,tmin,tmax,h,eps, ui->lineEdit_mu->text().toDouble(), partition);
       task->selectIndicator=chaoticIndicator;
       task->getInitialVector(x, C0, false);
       solvers.append(new MultiThreadSolver(*task));
       solvers.last()->start();

       do {

            if (solvers[0]->isOver()) {
               indValue = solvers[0]->getCurrentIndicatorValue(0)[0];
               solution = solvers[0]->getFinalSolution();
               std::cout << solvers[0]->printFinalSolution().toStdString() << std::endl;
               std::cout << "x = " << QString::number(double(solvers[0]-> getInitialVector()[2]),'g', 15).toStdString() << "; " << getIndicatorName(solvers[0]->getIndicator()).toStdString() << "R =" << " = " << MyString128::number(indValue,'g', 15).toStdString() << std::endl << std::endl;
               solvers[0]->quit();
               solvers[0]->wait();
               delete solvers[0];
               solvers.erase(solvers.begin()+0);
            }

       } while (!solvers.isEmpty());

       real_type firstDerivative = 2*(x - solution[2])*(1-solution[6]) + 2 * solution[0]*solution[4];
       std::cout<<"D(x) = " << MyString128::number(indValue,'g',15).toStdString() << std::endl;
       std::cout<<"D'(x) = " << MyString128::number(firstDerivative,'g',15).toStdString() << std::endl;
       //real_type secondDerivative = (-vals[4] + 16*vals[3] - 30 * vals[2] + 16*vals[1] -vals[0])/(12*dx*dx);
       //std::cout << "f'(x) = " << QString::number(double(firstDerivative),'g',10).toStdString() << std::endl << QString::number(double(secondDerivative),'g',10).toStdString() << std::endl;
       x_pred = x;
       x = x - indValue/firstDerivative;

   }
   */
/*
    if (x1IsComputing) {
        std::cout << "Периодическая орбита x = " << MyString128::number(x[2], 'g', 15).toStdString() << std::endl << "R(x) = " << MyString128::number(indValues[2],'g', 15).toStdString() << std::endl;
    }
    else    {
        std::cout << "Периодическая орбита x = " << MyString128::number(x[1], 'g', 15).toStdString() << std::endl << "R(x) = " << MyString128::number(indValues[1],'g', 15).toStdString() << std::endl;
    }
*/
   //std::cout << "Периодическая орбита x = " << MyString128::number(x[minPos], 'g', 15).toStdString() << std::endl << "R(x) = " << MyString128::number(indValues[minPos],'g', 15).toStdString() << std::endl;
   //std::cout << "Периодическая орбита x = " << MyString128::number(x, 'g', 15).toStdString() << std::endl << "R(x) = " << MyString128::number(indValue,'g', 15).toStdString() << std::endl;
   //std::cout << "Число итераций = " << nIterations << std::endl;
}

void MainWindow::on_pushButton_3_clicked()
{
    QString ss = QString::number(QInputDialog::getDouble(this,"Chaos Detector","Enter initial condition:",-0.2506555,-0.5,0.5,16),'g',16);
    ui->listWidget->addItem(ss);
}

void MainWindow::on_pushButton_4_clicked()
{
    qDeleteAll(ui->listWidget->selectedItems());
}

void MainWindow::on_radioButton_14_clicked(bool checked)
{
    if (checked) {
        ui->lineEdit_ncounts->setEnabled(true);
        return;
    }
    ui->lineEdit_ncounts->clear();
    ui->lineEdit_ncounts->setEnabled(false);
}

void MainWindow::on_pushButton_5_clicked()
{

}

void MainWindow::on_checkBox_autoscale_stateChanged(int checkState)
{

}



void MainWindow::on_pushButton_6_clicked()
{

}
