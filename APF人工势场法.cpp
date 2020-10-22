#include "APF_class.h"

int main()
{
    double obs[7][2] = { {2.5,1.5},{2,2},{1.5,2.5},{1,2},{2,1},{5.5,5.5},{8,8} };
    /*double obs[100][2];
    for (int i = 0; i < 5; i++)
    {
        obs[i][0] = double(2 * i);
        obs[i][1] = 6.5 + sin(obs[i][0]);
    }
    for (int i = 5; i < 10; i++)
    {
        obs[i][0] = double(2 * i) - 10;
        obs[i][1] = 3.5 + sin(obs[i][0]);
    }*/
    int n = sizeof(obs) / sizeof(double) / 2;//障碍物个数
    double Attraction_K = 200;//引力尺度因子
    double Repulsion_K = 4;//斥力尺度因子
    double Obstacles_dis = 4;//障碍物作用阀域
    double a = 2;////修正斥力函数系数
    double d = 5;//引力修正阀域
    double step = 0.1;//步长
    double start_point[2] = { 0,0 }, goal_point[2] = { 10,10 };//起始点和目标点
    double* angle_re, * Yatt, * Y;//分别用于存储角度，引力，斥力
    double velocity = 1.0;//预设速度
    double acceleration = 0;//加速度
    double m = 1000;//预设质量
    APF APF1(Attraction_K, Repulsion_K, Obstacles_dis, a, d, step, start_point, goal_point, obs);
    const int iterator = 2000;//迭代次数
    double path[iterator][2];
    int k1 = 0;
    double Current_point[2] = { 0,0 };//当前点
    //**************************参数的初始化**************************
    //random_obs(obs, n);
    Mat img(700, 700, CV_8UC3, Scalar(255, 255, 255));
    Point o(0, 1000);
    for (int i = 0; i < n; i++)
    {
        o.x = 100 + int(obs[i][0] * 50);
        o.y = 600 - int(obs[i][1] * 50);
        circle(img, o, 5, Scalar(0, 0, 0), -1);//画出障碍物，半径为5，颜色为黑色
    }
    o.x = int(goal_point[0]) * 50 + 100;
    o.y = 600 - int(goal_point[1]) * 50;
    circle(img, o, 5, Scalar(255, 0, 0), -1);//画出目标点，半径为5，颜色为蓝色
    Point p(0, 1000);//opencv坐标系的原点在左上角
    for (int j = 0; j < iterator - 1; j++)
    {
        path[j][0] = Current_point[0];
        path[j][1] = Current_point[1];
        angle_re = APF1.compute_angle(Current_point, n);
        Yatt = APF1.compute_attraction(Current_point, angle_re, d);
        Y = APF1.compute_repulsion(Current_point, angle_re, n);
        double Fsumyj = Yatt[1] + Y[1] + Y[3];
        double Fsumxj = Yatt[0] + Y[0] + Y[2];
        double S[2] = { Current_point[0] ,Current_point[1] };
        if ((Fsumyj == 0 && Fsumxj == 0) || sqrt(pow(path[j][0] - path[j - 2][0], 2) + pow(path[j][1] - path[j - 2][0], 2)) < 0.3 * step)//判断是否陷入局部最优解
        {
            double T = 10;
            while (0.85 * T)
            {
                srand(time(NULL));
                double x[2] = { Current_point[0] ,Current_point[1] };
                double dx = double(rand() % int(step * 1000)) / 1000;
                double theta = double(rand() % 728) / 100;
                double x1[2];
                x1[0] = x[0] + dx * cos(theta);
                x1[1] = x[1] + dx * sin(theta);
                double Ux = compute_attfield(x, goal_point, Attraction_K, d) - compute_repfield(x, n, goal_point, Repulsion_K, Obstacles_dis);
                double Ux1 = compute_attfield(x1, goal_point, Attraction_K, d) - compute_repfield(x1, n, goal_point, Repulsion_K, Obstacles_dis);
                double dUx = Ux - Ux1;
                if (dUx <= 0)
                {
                    Current_point[0] = x1[0];
                    Current_point[1] = x1[1];
                    Point s(0, 1000);
                    s.x = 100 + int(Current_point[0] * 50);
                    s.y = 600 - int(Current_point[1] * 50);
                    circle(img, s, 1, Scalar(0, 255, 0), -1);
                    if (sqrt(pow(Current_point[0] - S[0], 2) + pow(Current_point[1] - S[1], 2)) > 8 * step)
                        break;
                }
                else
                {
                    double P = pow(2.7, -dUx / T);
                    double a = double(rand() % 100) / 100;
                    if (P > a)
                    {
                        Current_point[0] = x1[0];
                        Current_point[1] = x1[1];
                        Point s(0, 1000);
                        s.x = 100 + int(Current_point[0] * 50);
                        s.y = 600 - int(Current_point[1] * 50);
                        circle(img, s, 1, Scalar(0, 255, 0), -1);
                        if (sqrt(pow(Current_point[0] - S[0], 2) + pow(Current_point[1] - S[1], 2)) > 8 * step)
                            break;
                    }
                }
                waitKey(20);
                imshow("APF", img);
            }
        }//模拟退火算法
        double Position_angle;
        if (Fsumxj > 0)//确定力的方向
            Position_angle = atan(Fsumyj / Fsumxj);//atan范围-pi/2到pi/2;
        else
            Position_angle = atan(Fsumyj / Fsumxj) + 3.1415;
        Current_point[0] += APF1.output_step() * cos(Position_angle);
        Current_point[1] += APF1.output_step() * sin(Position_angle);
        a = sqrt(Fsumyj * Fsumyj + Fsumxj * Fsumxj) / m;

        p.x = 100 + int(path[j][0] * 50);
        p.y = 600 - int(path[j][1] * 50);//opencv坐标系与惯性坐标系的转换
        circle(img, p, 1, Scalar(0, 0, 255), -1);//画出路径点，半径为1，颜色为红色
        if (fabs(Current_point[0] - APF1.output_goal_pointx()) < 0.1 && fabs(Current_point[1] - APF1.output_goal_pointy()) < 0.1)//在接近目标点时停止迭代
        {
            path[j + 1][0] = Current_point[0];
            path[j + 1][1] = Current_point[1];
            k1 = j;
            Point s(0, 1000);
            s.x = 100 + int(path[j + 1][0] * 50);
            s.y = 600 - int(path[j + 1][1] * 50);
            circle(img, s, 1, Scalar(0, 0, 255), -1);
            imshow("APF", img);
            cout << "迭代次数" << k1 << endl;
            break;
        }

        double t = 1000 * step / velocity;//单位是毫秒
        velocity += t / 1000 * a;
        imshow("APF", img);
        waitKey(int(t));
    }
    waitKey(0);
    return 0;
}
/*************************************************
函数名称：random_obs
类型：辅助函数
作用：随机生成障碍物坐标
参数：obs(障碍物坐标数组)，n(障碍物个数)
*************************************************/
void random_obs(double(*obs)[2], int n)
{
    srand(time(NULL));
    for (int i = 0; i < n; i++)
    {
        obs[i][0] = double(rand() % 1000) / 100;
        obs[i][1] = double(rand() % 1000) / 100;
    }
}
/*************************************************
函数名称：compute_attfield
类型：计算函数
作用：计算引力势
参数：current_point(当前点)，goal_point(目标点),Attraction_K(引力尺度因子),d(引力修正阀域)
*************************************************/
double compute_attfield(double* Current_point, double* goal_point, double Attraction_K, double d)
{
    double r = sqrt(pow((goal_point[0] - Current_point[0]), 2) + pow((goal_point[1] - Current_point[1]), 2));
    if (r < d)
        return 0.5 * Attraction_K * r * r;
    else
        return 0.5 * Attraction_K * r * d;
}
/*************************************************
函数名称：compute_repfield
类型：计算函数
作用：计算斥力势
参数：current_point(当前点)，n(障碍物个数),goal_point(目标点),Repulsion_K(斥力尺度因子),Obstacles_dis(障碍物作用阀域)
*************************************************/
double compute_repfield(double* Current_point, int n, double* goal_point, double Repulsion_K, double Obstacles_dis)
{
    double sum = 0;
    for (int i = 0; i < n; i++)
    {

        double r = sqrt(pow((goal_point[0] - Current_point[0]), 2) + pow((goal_point[1] - Current_point[1]), 2));
        if (r < Obstacles_dis)
            sum += 0.5 * Repulsion_K * pow((1 / r - 1 / Obstacles_dis), 2);
        else
            ;
    }
    return sum;
}