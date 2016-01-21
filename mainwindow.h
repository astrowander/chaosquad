#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

#include "multithreadsolver.h"
#include "plotwindow.h"
#include <QFileDialog>

#include <QQueue>

#include <QInputDialog>
#include <QStandardPaths>

const int maxThreads = 8;
namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private slots:
    void on_pushButton_clicked();

    void on_pushButton_2_clicked();

   // void on_checkBox_3_stateChanged(int arg1);

    void on_pushButton_3_clicked();

    void on_pushButton_4_clicked();

    void on_radioButton_14_clicked(bool checked);

    void on_pushButton_5_clicked();

   // void on_checkBox_autoscale_clicked();

    void on_checkBox_autoscale_stateChanged(int arg1); 

    void on_pushButton_6_clicked();
 
private:
    int dim;
    Ui::MainWindow *ui;    
    QVector<real_type> initialConditions;
    QQueue<Task> tasksQueue;
    QList<MultiThreadSolver*> solvers;
    QString fileName;
    TimeIntervalPartition partition;
    QDir workingDirectory;
};

#endif // MAINWINDOW_H
