#ifndef class_APF
#define class_APF
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <opencv2/opencv.hpp>
#include <ctime>

using namespace std;
using namespace cv;
void random_obs(double(*obs)[2], int n);//������������
double compute_attfield(double* Current_point, double* goal_point, double Attraction_K, double d);
double compute_repfield(double* Current_point, int n, double* goal_point, double Repulsion_K, double Obstacles_dis);
class APF
{
private:
    double Attraction_K;//�����߶�����
    double Repulsion_K;//�����߶�����
    double Obstacles_dis;//�ϰ����÷���
    double a;//������������ϵ��
    double d;//������������
    double step;//����
    double* start_point;//��ʼ��
    double* goal_point;//Ŀ���
    double(*obstacles)[2];//�ϰ�������
public:
    APF(double Attraction_K, double Repulsion_K, double Obstacles_dis, double a, double d, double step, double* start_point, double* goal_point, double(*obstacles)[2]) //����������
    {
        this->Attraction_K = Attraction_K;
        this->Repulsion_K = Repulsion_K;
        this->Obstacles_dis = Obstacles_dis;
        this->a = a;
        this->d = d;
        this->step = step;
        this->start_point = start_point;
        this->goal_point = goal_point;
        this->obstacles = obstacles;
    }//���캯��
    double* compute_angle(double* Current_point, int n);
    double* compute_attraction(double* Current_point, double* angle, double d);
    double* compute_repulsion(double* Current_point, double* angle, int n);
    double sum(double* p1, int n);
    double output_step();
    double output_goal_pointx();
    double output_goal_pointy();
};
/*************************************************
�������ƣ�APF::compute_angle
���ͣ����㺯��
���ã����㵱ǰ����Ŀ����Լ������ϰ���֮��ļн�
������Current_point(��ǰ��)��n(�ϰ������)
*************************************************/
double* APF::compute_angle(double* Current_point, int n)
{
    double* Y = new double[++n];//Y���ڴ洢�н�
    double deltax, deltay, r;
    for (int i = 0; i < n; i++)
    {
        if (i != 0)
        {
            deltax = this->obstacles[i - 1][0] - Current_point[0];
            deltay = this->obstacles[i - 1][1] - Current_point[1];
        }
        else
        {
            deltax = this->goal_point[0] - Current_point[0];
            deltay = this->goal_point[1] - Current_point[1];
        }
        r = sqrt(deltax * deltax + deltay * deltay);
        if (deltay > 0)
            Y[i] = acos(deltax / r);
        else
            Y[i] = -acos(deltax / r);
    }//��нǵĹ���
    return Y;
};
/*************************************************
�������ƣ�APF::compute_attraction
���ͣ����㺯��
���ã����㵱ǰ���ܵ�������
������Current_point(��ǰ��)��angle(�н�����)��d(������������)
*************************************************/
double* APF::compute_attraction(double* Current_point, double* angle, double d)
{
    double R = pow((this->goal_point[0] - Current_point[0]), 2) + pow((this->goal_point[1] - Current_point[1]), 2);
    double r = sqrt(R);
    static double Yatt[2];
    if (r < d)
    {
        Yatt[0] = this->Attraction_K * r * cos(angle[0]);
        Yatt[1] = this->Attraction_K * r * sin(angle[0]);
    }
    else
    {
        Yatt[0] = this->Attraction_K * d * cos(angle[0]);
        Yatt[1] = this->Attraction_K * d * sin(angle[0]);
    }//��������
    return Yatt;
}
/*************************************************
�������ƣ�APF::compute_repulsion
���ͣ����㺯��
���ã����㵱ǰ���ܵ��ĳ���
������Current_point(��ǰ��)��angle(�н�����)��n(�ϰ������)
*************************************************/
double* APF::compute_repulsion(double* Current_point, double* angle, int n)
{
    double* YY = new double[4];//�洢���ճ���
    double Rat = pow((Current_point[0] - this->goal_point[0]), 2) + pow((Current_point[1] - this->goal_point[1]), 2);
    double rat = sqrt(Rat);
    double Rre, rre, Yrer, Yata;
    double* Yrerx = new double[n], * Yrery = new double[n], * Yatax = new double[n], * Yatay = new double[n];
    for (int i = 0; i < n; i++)
    {
        Rre = pow((Current_point[0] - this->obstacles[i][0]), 2) + pow((Current_point[1] - this->obstacles[i][1]), 2);
        rre = sqrt(Rre);
        if (rre > this->Obstacles_dis)
        {
            Yrerx[i] = 0;
            Yrery[i] = 0;
            Yatax[i] = 0;
            Yatay[i] = 0;
        }
        else
        {
            Yrer = this->Repulsion_K * (1 / rre - 1 / this->Obstacles_dis) * (1 / Rre) * pow(rat, a);//���Fre1(����)����
            Yata = a * this->Repulsion_K * pow((1 / rre - 1 / this->Obstacles_dis), 2) * pow(rat, a - 1) / 2;//���Fre2(����)����
            Yrerx[i] = Yrer * cos(angle[i + 1] + 3.1415);
            Yrery[i] = Yrer * sin(angle[i + 1] + 3.1415);
            Yatax[i] = Yata * cos(angle[0]);
            Yatay[i] = Yata * sin(angle[0]);
        }
    }
    YY[0] = sum(Yrerx, n);
    YY[1] = sum(Yrery, n);
    YY[2] = sum(Yatax, n);
    YY[3] = sum(Yatay, n);
    return YY;
}
/*************************************************
�������ƣ�APF::sum
���ͣ���������
���ã���double�����Ա���
������p1(������)��n(�����Ա����)
*************************************************/
double APF::sum(double* p1, int n)
{
    double sum = 0;
    for (int i = 0; i < n; i++)
    {
        sum += p1[i];
    }
    return sum;
}
/*************************************************
�������ƣ�APF::output_step
���ͣ���������
���ã����˽�����ݳ�Աstep
*************************************************/
double APF::output_step()
{
    return step;
}
/*************************************************
�������ƣ�APF::goal_pointx
���ͣ���������
���ã����˽�����ݳ�Աgoal_point[0]
*************************************************/
double APF::output_goal_pointx()
{
    return goal_point[0];
}
/*************************************************
�������ƣ�APF::output_pointy
���ͣ���������
���ã����˽�����ݳ�Աgoal_point[1]
*************************************************/
double APF::output_goal_pointy()
{
    return goal_point[1];
}
#endif 

