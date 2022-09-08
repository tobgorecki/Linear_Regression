#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();


private slots:
    void on_upload_run_plot_clicked();
    void on_exit_aplication_clicked();
    void read_file_with_data(QString filename, double *x, double *y , double  *err_y, int num_lines);
    int count_number_lines_data(QString filename);
    void plot_result(std::vector<double> x_v, std::vector<double> y_v, std::vector<double> y_err,
                                 double para_A_classic, double para_B_classic,
                                 double para_A_weighted, double para_B_weighted,
                                 double para_A_gradient_descent, double para_B_gradient_descent,
                                 double para_A_mcmc, double para_B_mcmc, int num_lines);
    void update_progress(QString progress);



private:
    Ui::MainWindow *ui;


};
#endif // MAINWINDOW_H
