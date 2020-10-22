#ifndef class_APF
#define class_APF
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <opencv2/opencv.hpp>
#include <ctime>
#include <windows.h>
using namespace std;
using namespace cv;
void random_obs(double(*obs)[2], int n);//辅助函数声明
double compute_attfield(double* Current_point, double* goal_point, double Attraction_K, double d);
double compute_repfield(double* Current_point, int n, double* goal_point, double Repulsion_K, double Obstacles_dis);
class APF
{
private:
    double Attraction_K;//引力尺度因子
    double Repulsion_K;//斥力尺度因子
    double Obstacles_dis;//障碍作用阀域
    double a;//修正斥力函数系数
    double d;//引力修正阀域
    double step;//步长
    double* start_point;//起始点
    double* goal_point;//目标点
    double(*obstacles)[2];//障碍物数组
public:
    APF(double Attraction_K, double Repulsion_K, double Obstacles_dis, double a, double d, double step, double* start_point, double* goal_point, double(*obstacles)[2]) //将各个参数
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
    }//构造函数
    double* compute_angle(double* Current_point, int n);
    double* compute_attraction(double* Current_point, double* angle, double d);
    double* compute_repulsion(double* Current_point, double* angle, int n);
    double sum(double* p1, int n);
    double output_step();
    double output_goal_pointx();
    double output_goal_pointy();
};
/*************************************************
函数名称：APF::compute_angle
类型：计算函数
作用：计算当前点与目标点以及各个障碍物之间的夹角
参数：Current_point(当前点)，n(障碍物个数)
*************************************************/
double* APF::compute_angle(double* Current_point, int n)
{
    double* Y = new double[++n];//Y用于存储夹角
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
    }//求夹角的过程
    return Y;
};
/*************************************************
函数名称：APF::compute_attraction
类型：计算函数
作用：计算当前点受到的引力
参数：Current_point(当前点)，angle(夹角数组)，d(引力修正阀域)
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
    }//引力修正
    return Yatt;
}
/*************************************************
函数名称：APF::compute_repulsion
类型：计算函数
作用：计算当前点受到的斥力
参数：Current_point(当前点)，angle(夹角数组)，n(障碍物个数)
*************************************************/
double* APF::compute_repulsion(double* Current_point, double* angle, int n)
{
    double* YY = new double[4];//存储最终斥力
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
            Yrer = this->Repulsion_K * (1 / rre - 1 / this->Obstacles_dis) * (1 / Rre) * pow(rat, a);//求的Fre1(斥力)向量
            Yata = a * this->Repulsion_K * pow((1 / rre - 1 / this->Obstacles_dis), 2) * pow(rat, a - 1) / 2;//求的Fre2(引力)向量
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
函数名称：APF::sum
类型：辅助函数
作用：对double数组成员求和
参数：p1(数组名)，n(数组成员个数)
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
函数名称：APF::output_step
类型：辅助函数
作用：输出私有数据成员step
*************************************************/
double APF::output_step()
{
    return step;
}
/*************************************************
函数名称：APF::goal_pointx
类型：辅助函数
作用：输出私有数据成员goal_point[0]
*************************************************/
double APF::output_goal_pointx()
{
    return goal_point[0];
}
/*************************************************
函数名称：APF::output_pointy
类型：辅助函数
作用：输出私有数据成员goal_point[1]
*************************************************/
double APF::output_goal_pointy()
{
    return goal_point[1];
}
#endif 

