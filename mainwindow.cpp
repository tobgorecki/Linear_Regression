#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "fit_model/fitting_linear_model.h"
#include <stdio.h>
#include <QDebug>
#include <QString>
#include <QtMath>
#include <vector>
#include <iostream>
#include <plot/qcustomplot.h>
#include <QFile>
MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);


    QString windowTitle("Linear Regression v1.0");
    this->setWindowTitle(windowTitle);
    this->setStyleSheet("background-color: qlineargradient(spread:pad, x1:0, y1:0, x2:1, y2:1, stop:0" " #159957, stop:1 #155799)");
    //this->setStyleSheet("background-color: qlineargradient(spread:pad, x1:0, y1:0, x2:1, y2:1, stop:0" " #44a08d, stop:1 #093637)");



    ui->LSM->setStyleSheet("background-color: rgb(213, 226, 188);");
    ui->LSM_w->setStyleSheet("background-color: rgb(213, 226, 188);");
    ui->GD->setStyleSheet("background-color: rgb(213, 226, 188);");
    ui->MCMC->setStyleSheet("background-color: rgb(213, 226, 188);");
    ui->Progress->setStyleSheet("background-color: rgb(213, 226, 188);");

    ui->upload_run_plot->setStyleSheet("background-color: rgb(220,220,220);");
    ui->exit_aplication->setStyleSheet("background-color: rgb(220,220,220) ;");

    this->setWindowTitle(windowTitle);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::update_progress(QString progress){

    ui->progress_value->setText(progress);

}

void MainWindow::plot_result(std::vector<double> x_v,std::vector<double> y_v,std::vector<double> y_err,
                             double para_A_classic, double para_B_classic,
                             double para_A_weighted, double para_B_weighted,
                             double para_A_gradient_descent, double para_B_gradient_descent,
                             double para_A_mcmc, double para_B_mcmc,
                             int num_lines)
{


    QVector<double> x(num_lines), y(num_lines),y_error(num_lines),
                     model_classic(num_lines),
                     model_weighted(num_lines),
                     model_gradient_descent(num_lines),
                     model_mcmc(num_lines);


    for (int iter=0; iter<num_lines; ++iter)
    {
      x[iter] = x_v[iter];
      y[iter] = y_v[iter];
      y_error[iter] = y_err[iter];
      model_classic[iter] = para_A_classic * x_v[iter] + para_B_classic;
      model_weighted[iter] =  para_A_weighted * x_v[iter] + para_B_weighted;
      model_gradient_descent[iter] = para_A_gradient_descent * x_v[iter] + para_B_gradient_descent;
      model_mcmc[iter] = para_A_mcmc * x_v[iter] + para_B_mcmc;

    }
    // remove all items in legend
     ui->plot->legend->clearItems();
     ui->plot->clearGraphs();

    // add legend
    ui->plot->legend->setVisible(true);
    ui->plot->legend->setFont(QFont("Helvetica", 13));
    ui->plot->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignTop |Qt::AlignLeft); //legend position
    QStringList lineNames;
    lineNames <<"Data" <<"Least square classic" << "Least square weighted" << "Gradient descent" << "MCMC" ;



    // ploting original data (from file )
    ui->plot->addGraph();
    ui->plot->graph()->setData(x, y);
    ui->plot->graph()->setLineStyle(QCPGraph::lsNone);
    ui->plot->graph()->setName(lineNames.at(0));
    ui->plot->graph()->setPen(QPen(Qt::black));
    ui->plot->graph()->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 5));


    QCPErrorBars *errorBars = new QCPErrorBars(ui->plot->xAxis,ui->plot->yAxis);
    errorBars->removeFromLegend();
    errorBars->setDataPlottable(ui->plot->graph());
    errorBars->setData((y_error));


    // change axis size

    QFont font_axis("Helvetica [Cronyx]", 14);
    font_axis.setPointSize(14);
    ui->plot->xAxis->setTickLabelFont(font_axis);
    ui->plot->yAxis->setTickLabelFont(font_axis);
    ui->plot->xAxis->setLabelFont(font_axis);
    ui->plot->yAxis->setLabelFont(font_axis);
    ui->plot->xAxis->setLabel("x");
    ui->plot->yAxis->setLabel("y");

    // set axes ranges
    auto  minmax_x = std::minmax_element(x_v.begin(), x_v.end());
    auto  minmax_y = std::minmax_element(y_v.begin(), y_v.end());
    int delta =3 ;
    ui->plot->xAxis->setRange(*minmax_x.first-delta , *minmax_x.second+delta );
    ui->plot->yAxis->setRange(*minmax_y.first-5*delta , *minmax_y.second+5*delta);
    ui->plot->replot();


    // plotting model classic
    ui->plot->addGraph();
    ui->plot->graph()->setData(x, model_classic);
    ui->plot->graph()->setName(lineNames.at(1));
    ui->plot->graph()->setPen(QPen(Qt::red));
    ui->plot->replot();

    // plotting model weighted
    ui->plot->addGraph();
    ui->plot->graph()->setData(x,model_weighted);
    ui->plot->graph()->setName(lineNames.at(2));
    ui->plot->graph()->setPen(QPen(Qt::blue));
    ui->plot->replot();

    // plotting model gradient descent
    ui->plot->addGraph();
    ui->plot->graph()->setData(x, model_gradient_descent);
    ui->plot->graph()->setName(lineNames.at(3));
    ui->plot->graph()->setPen(QPen(Qt::green));
    ui->plot->replot();

    // plotting  model mcmc
    ui->plot->addGraph();
    ui->plot->graph()->setData(x, model_mcmc);
    ui->plot->graph()->setName(lineNames.at(4));
    ui->plot->graph()->setPen(QPen(QColorConstants::Svg::purple ));


    // update legend
   ui->plot->legend->updateLayout();
    ui->plot->replot();
    ui->plot->update();

   //delete errorBars;

}

void MainWindow::on_upload_run_plot_clicked()
{
    // read data
    QMessageBox msgBox;
    msgBox.setText("Format must be x-axis, y-axis, error_y");
    msgBox.setIcon(QMessageBox::Icon::Warning);
    msgBox.exec();


    QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"), "Linear_Regression");

    if (fileName.isEmpty())
    {
        QMessageBox::information(this, tr("ERROR"), tr("Data not upload"));
        return;
    }

    int number_lines=0;
     number_lines =count_number_lines_data(fileName );


     double *x  = new double[number_lines];
     double *y  = new double[number_lines];
     double *err_y  = new double[number_lines];


    read_file_with_data( fileName,x ,y,err_y,number_lines);

     std::vector<double> x_value;
     std::vector<double> y_value;
     std::vector<double> y_error;

     for (auto i =0; i < number_lines; ++i){
            x_value.push_back(x[i]);
            y_value.push_back(y[i]);
            y_error.push_back(err_y[i]);
     }
      //####################### computing model parameters ##############################

     double para_A_classic=0.0,  para_B_classic=0.0, para_A_classic_error=0.0,para_B_classic_error=0.0;
     double para_A_weighted=0.0, para_B_weighted=0.0,para_A_weighted_error=0.0, para_B_weighted_error=0.0;
     double para_A_gradient_descent=0.0,  para_B_gradient_descent=0.0;
     double para_A_mcmc=0.0,  para_B_mcmc = 0.0;
     double correlation_coefficient=0.0;
     Fitting_linear_model  model(0,0);

     QString   progress = "Start...\n";

     // least square classic method
       model.Least_sqaure_classic( number_lines, x_value,y_value,
                                    &para_A_classic , &para_B_classic,
                                    &para_A_classic_error , &para_B_classic_error,&correlation_coefficient);


       ui->lsm_A_result->setText( QString::number(para_A_classic));
       ui->lsm_B_result->setText(QString::number(para_B_classic));
       ui->lsm_A_error_result->setText(QString::number(para_A_classic_error));
       ui->lsm_B_error_result->setText(QString::number(para_B_classic_error));
       ui->lsm_r_result->setText(QString::number(correlation_coefficient));
         progress += "Least sqaure-- > OK \n";

   // least square weighted method
     model.Least_sqaure_weighted( number_lines, x_value,y_value,y_error ,
                                  &para_A_weighted , &para_B_weighted,
                                  &para_A_weighted_error , &para_B_weighted_error,&correlation_coefficient);



     ui->lsm_w_A_result->setText( QString::number(para_A_weighted));
     ui->lsm_w_B_result->setText(QString::number(para_B_weighted));
     ui->lsm_w_A_error_result->setText(QString::number(para_A_weighted_error));
     ui->lsm_w_B_error_result->setText(QString::number(para_B_weighted_error));
     ui->lsm_w_r_result->setText(QString::number(correlation_coefficient));
     progress += "Least sqaure weighted-- > OK \n";


    // gradient descent
     model.Gradinet_descent(number_lines, x_value,y_value,y_error ,
                            &para_A_gradient_descent, &para_B_gradient_descent, para_A_weighted, para_B_weighted);

     ui->gd_A_result->setText( QString::number(para_A_gradient_descent));
     ui->gd_B_result->setText(QString::number(para_B_gradient_descent));
     progress += "Gradinet descent -- > OK\n";
     update_progress(progress);


     // mcmc
     model.Monte_Carto_Markov_Chains(number_lines, x_value,y_value,y_error ,
                                     &para_A_mcmc , &para_B_mcmc,para_A_weighted, para_B_weighted);

     ui->mcmc_A_result->setText( QString::number(para_A_mcmc));
     ui->mcmc_B_result->setText(QString::number(para_B_mcmc));
     progress += "MCMC -- > OK\n";
     update_progress(progress);

    //######################### plotting part #############################################


     plot_result(x_value,y_value,y_error ,
                   para_A_classic, para_B_classic,
                   para_A_weighted, para_B_weighted,
                   para_A_gradient_descent, para_B_gradient_descent,
                   para_A_mcmc, para_B_mcmc,
                   number_lines);

    // deallocated memory
    delete[] x;
    delete[] y;
    delete [] err_y;

}



void  MainWindow::read_file_with_data(const QString filename, double *x , double *y , double  *err_y, int num_lines)
{
    QFile file(filename);
    if(!file.open(QFile::ReadOnly |
                  QFile::Text))
    {
        qDebug() << " Could not open the file for reading";
        return;
    }

    QTextStream in(&file);
    int iter=0;

       while ( iter < num_lines) {
        QString n1,n2,n3;
         in >> n1 >> n2>>n3;
         x[iter] = n1.toDouble();
         y[iter] = n2.toDouble();
         err_y[iter] = n3.toDouble();
        // qDebug() <<QString::number(x[iter])<<QString::number(y[iter])<<QString::number(err_y[iter]);
        iter +=1;
      }

    file.close();

}



int MainWindow::count_number_lines_data(const QString filename)//,double *X_value, double *Y_value
{
    QFile file(filename);
    if(!file.open(QFile::ReadOnly |
                  QFile::Text))
    {
        qDebug() << " Could not open the file for counting lines";
        return 1;
    }
    QTextStream in(&file);
    int iter=0;
    while (!in.atEnd()) {

          QString line = in.readLine();
            if (line == "") //read until white line occurring
                   break;
          iter +=1;
      }


    file.close();
    return iter;
}



void MainWindow::on_exit_aplication_clicked()
{
    close();
}




