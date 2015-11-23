
/*************************************************
模型：子网格，文献SIW_filter，波导过渡，14个圆孔(d=2mm)，孔径固定2*dx_coarse，缝隙宽s=3mm,3*dx_coarse
      计算电场Ey的总场值，并利用所得场值计算S11和S21
	  重叠型DDM：重叠区域为1个均匀网格
	  子网格法，二次加细，3:2的空隙圆柱比。
参考程序： 樊瑞SIW_filter(矩阵剖分）
尺寸：a=17mm，d0=18mm，d1=9mm，w=2mm
频率：7.49GHz;  波长：27mm;  介电常数：2.2
激励：沿z向打入TE波主模TE10：sin(pi*x/a)*exp(-j*kz*dx)
*************************************************/
#include<iostream>
#include <math.h>
#include <time.h>
#include <fstream>
#include <iomanip>
#include "BasicFunc.h"
#include "Complex.h"
#include "Vector.h"
using namespace std;

//Fundamental constants
const double pi=3.1415926;
//   freq:   6.86;  6.96; 7.04;  7.06;  7.12;  7.14; 7.17;  7.20;  7.21;  7.22;  7.23;  7.31 ;  7.40
//      N:  118.0; 116.3; 114.9; 114.6; 113.7; 113.4; 112.9; 112.4; 112.3; 112.1; 112.0; 110.75; 109.4
const double N_dx=112.9; //可修改!!!
const double dx=1.0/N_dx;     //步长,均匀网格   dx_coarse=1mm, dx_sub1=0.5mm, dx_sub2=0.25mm
const double a=(17.0*4)/N_dx;       //输入波导x向归一化宽度,y向不计,得出主模波长介于(a,2a)之间 a=17/... !!!!!!!
const double k=2*pi;        //传播常数k
const double kz=sqrt ( k*k - (pi/a)*(pi/a) );   //传播常数kz
const int N_d0=18;      //x向外部空间长度  d0/4dx=18   d0=18mm, 为2/3*lambda
const int N_d1=9;      //z向波导过渡长度  d1/4dx=9   d1=9mm，为1/3*lambd
const int N_w=2;        //周期金属圆孔直径  w/4dx=2   w=2mm
const int N_s=3;        //周期缝隙z向长度  s/4dx=3
const int N_slot=1;     //ddm缝隙较长部分  4(可调)/5
const int N_slot2=0;    //ddm缝隙较短部分  1(可调)/2
//const int N_slot3=2;    //ddm缝隙较短部分  1(可调)/2,new!!!
const int N_a=16;       //波导端口x向长度  a/4dx-1=16  a=17mm
const int Nx_side=52;      //波导段整体x向总长  2N_d0+N_a=53
const int Nx_1=26;      //..中间段整体x1向总长  2N_d0+N_a+2N_w+2=63
const int Nx_2=74;      //..中间段整体x1向总长  2N_d0+N_a+2N_w+2=63
const int Nx_3=30;      //..中间段整体x2向总长  2N_d0+N_a+2N_w+2=63
const bool Debug = true;

//基本圆孔区域参数*******************************************************
//circle layer_nums
const int Num1=74, Num2=26, Num3=46, Num4=30, Num5=80, Num6=16, 
          Num7=28, Num8=12, Num9=72;
//1circle last node num圆孔部分每行第一个圆最下一个节点编号
const int Node_1=30, Node_5=33, Node_9=29;
//1circle last node num圆孔部分每行第二个圆最下一个节点编号
const int Node_11=57, Node_55=63, Node_99=55;
//1circle last node num圆孔部分靠近园表面的几个特殊点位置
const int Node11=20,  Node12=26,  Node13=47,  Node14=53, 
          Node41=4,   Node42=10,  Node43=19,  Node44=25, 
          Node51=24,  Node52=54, 
		  Node61=3,   Node62=11, 
		  Node71=6,   Node72=20, 
          Node91=22,  Node92=48;
//total num
//两边过渡波导节点数 Nx_side*(N_d1+1)+N_slot*(Nx_2+Nx_3)+Nx_2+Nx_1=724
const int M_end=Nx_side*(N_d1+1)+N_slot*(Nx_2+Nx_3)+Nx_2+Nx_1;  
//中间含圆孔周期节点数 2*(Num1+Num2+Num3+Num4+Num5+Num6+Num7+Num8)+Num9+2*N_slot*(Nx_2+Nx_3)=904
const int M_middle=2*(Num1+Num2+Num3+Num4+Num5+Num6+Num7+Num8)+Num9+2*N_slot*(Nx_2+Nx_3);


//含post1圆孔区域参数*******************************************************
//circle layer_nums
const int Num1_p1=74,  Num2_p1=26,  Num3_p1=46,  Num4_p1=30,  Num5_p1=80, 
          Num6_p1=16,  Num7_p1=28,  Num8_p1=12,  Num9_p1=72,  Num10_p1=12,
          Num11_p1=35, Num12_p1=16, Num13_p1=85, Num14_p1=30, Num15_p1=55,
          Num16_p1=31, Num17_p1=83, Num18_p1=4,  Num19_p1=40;
//1circle last node num圆孔部分每行第一个圆最下一个节点编号
const int Node_1_p1=30,  Node_5_p1=33,  Node_9_p1=29,  Node_13_p1=33,    Node_17_p1=30;
//1circle last node num圆孔部分每行第二个圆最下一个节点编号
const int Node_11_p1=57, Node_55_p1=63, Node_99_p1=55, Node_1313_p1=68,  Node_1717_p1=66;
//1circle last node num圆孔部分靠近园表面的几个特殊点位置
const int Node11_p1=20,  Node12_p1=26,  Node13_p1=47,  Node14_p1=53, 
          Node41_p1=4,   Node42_p1=10,  Node43_p1=19,  Node44_p1=25, 
          Node51_p1=24,  Node52_p1=54, 
          Node61_p1=3,   Node62_p1=11,
          Node71_p1=6,   Node72_p1=20, 
          Node91_p1=22,  Node92_p1=34,  Node93_p1=37,  Node94_p1=48,
          Node111_p1=6,  Node112_p1=13, Node113_p1=21, Node114_p1=27,
          Node121_p1=3,  Node122_p1=11, 
          Node131_p1=24, Node132_p1=37, Node133_p1=47, Node134_p1=59, 
          Node141_p1=4,  Node142_p1=10, Node143_p1=19, Node144_p1=25, 
          Node151_p1=22, Node152_p1=26, Node153_p1=32,  
          Node161_p1=12, Node162_p1=18,  
          Node171_p1=20, Node172_p1=26, Node173_p1=34, Node174_p1=48, Node175_p1=56, Node176_p1=62,
          Node191_p1=14, Node192_p1=19, Node193_p1=25;
//total num
//含小post1的圆孔周期节点数!!:      // M_post1=1718
const int M_post1=2*(Num1_p1+Num2_p1+Num3_p1+Num4_p1+Num5_p1+Num6_p1+Num7_p1+Num8_p1+Num9_p1+Num10_p1
					 +Num11_p1+Num12_p1+Num13_p1+Num14_p1+Num15_p1+Num16_p1+Num17_p1+Num18_p1)
					 +Num19_p1+2*N_slot*(Nx_2+Nx_3);

const int M_post1_q1=Nx_2+Nx_3+Num1_p1+Num2_p1+Num3_p1+Num4_p1+Num5_p1+Num6_p1+Num7_p1+Num8_p1+Num9_p1;

const int M_post1_q2=M_post1_q1+Num10_p1+Num11_p1+Num12_p1+Num13_p1+Num14_p1
                    +Num15_p1+Num16_p1+Num17_p1+Num18_p1+Num19_p1;
const int M_post1_q3=M_post1_q2+Num9_p1+ Num10_p1+Num11_p1+Num12_p1+Num13_p1
                    +Num14_p1+Num15_p1+Num16_p1+Num17_p1+Num18_p1;


//含post2圆孔区域参数*******************************************************
//circle layer_nums
const int Num1_p2=74,  Num2_p2=26,  Num3_p2=46,  Num4_p2=30,  Num5_p2=80, 
          Num6_p2=16,  Num7_p2=39,  Num8_p2=12,  Num9_p2=79,  Num10_p2=12,
          Num11_p2=41, Num12_p2=29, Num13_p2=95, Num14_p2=45, Num15_p2=62,
          Num16_p2=34, Num17_p2=80, Num18_p2=6,  Num19_p2=42;
//1circle last node num圆孔部分每行第一个圆最下一个节点编号
const int Node_1_p2=30,  Node_5_p2=33,  Node_9_p2=29,  Node_13_p2=33,    Node_17_p2=30;
//1circle last node num圆孔部分每行第二个圆最下一个节点编号
const int Node_11_p2=57, Node_55_p2=63, Node_99_p2=62, Node_1313_p2=78,  Node_1717_p2=63;
//1circle last node num圆孔部分靠近园表面的几个特殊点位置
const int Node11_p2=20,  Node12_p2=26,  Node13_p2=47,  Node14_p2=53, 
          Node41_p2=4,   Node42_p2=10,  Node43_p2=19,  Node44_p2=25, 
          Node51_p2=24,  Node52_p2=37,  Node53_p2=42,  Node54_p2=54, 
          Node61_p2=3,   Node62_p2=11,
          Node71_p2=6,   Node72_p2=13,  Node73_p2=25,  Node74_p2=31, 
          Node91_p2=22,  Node92_p2=33,  Node93_p2=45,  Node94_p2=55,
          Node111_p2=6,  Node112_p2=13, Node113_p2=27, Node114_p2=33,
          Node121_p2=3,  Node122_p2=7,  Node123_p2=21, Node124_p2=24, 
          Node131_p2=24, Node132_p2=37, Node133_p2=57, Node134_p2=69, 
          Node141_p2=4,  Node142_p2=10, Node143_p2=19, Node144_p2=25, Node145_p2=34, Node146_p2=40, 
          Node151_p2=22, Node152_p2=30, Node153_p2=39,  
          Node161_p2=12, Node162_p2=16, Node163_p2=21,  
          Node171_p2=20, Node172_p2=26, Node173_p2=34, Node174_p2=45, Node175_p2=53, Node176_p2=59,
          Node191_p2=14, Node192_p2=20, Node193_p2=27;
//total num
//含大post2的圆孔周期节点数!!:      // M_post2=1862 !     506/952/1435
const int M_post2=2*(Num1_p2+Num2_p2+Num3_p2+Num4_p2+Num5_p2+Num6_p2+Num7_p2+Num8_p2+Num9_p2+Num10_p2
					+Num11_p2+Num12_p2+Num13_p2+Num14_p2+Num15_p2+Num16_p2+Num17_p2+Num18_p2)
					+Num19_p2+2*N_slot*(Nx_2+Nx_3);

const int M_post2_q1=Nx_2+Nx_3+Num1_p2+Num2_p2+Num3_p2+Num4_p2+Num5_p2+Num6_p2+Num7_p2+Num8_p2+Num9_p2;

const int M_post2_q2=M_post2_q1+Num10_p2+Num11_p2+Num12_p2+Num13_p2+Num14_p2
                               +Num15_p2+Num16_p2+Num17_p2+Num18_p2+Num19_p2;
const int M_post2_q3=M_post2_q2+Num9_p2+ Num10_p2+Num11_p2+Num12_p2+Num13_p2
                               +Num14_p2+Num15_p2+Num16_p2+Num17_p2+Num18_p2;


//region total=10;  这里计算10周期（缝隙加孔）,重叠区域为两个sub1子网格
const int N_circle=12;
const int N_end=2;
const int N_middle=8;
const int N_post1=2;
const int N_post2=1;
const int N_c1=5,N_c2=7,N_c3=9;
//node total=

const int N_matrix=8;                //均匀网格，5点格式
complex Einc0[N_a],Einc1[N_a];         //前两层点的入射波，source中初始化

void Source()	//4dx!!!
{
	int j=0;
	for( j=0; j<N_a; j++ )
	{
		Einc0[j].real=sin( pi*(j+1)*4*dx/a );//边界
		Einc0[j].image=0;
		Einc1[j].real=sin( pi*(j+1)*4*dx/a ) * cos( kz*4*dx );//第一层
		Einc1[j].image=-sin( pi*(j+1)*4*dx/a ) * sin( kz*4*dx );
	}
	
};

//边界节点参数的扫描//

//一类过度节点参数
const double A1_1=3*(10+k*k*4*dx*dx),    A1_2=3*(10+k*k*dx*dx);
const double A2_1=6-k*k*4*dx*dx,         A2_2=6-k*k*dx*dx;
const double A3_1=96-32*k*k*4*dx*dx,     A3_2=96-32*k*k*dx*dx;
const double A4_1=48*(2-k*k*4*dx*dx),    A4_2=48*(2-k*k*dx*dx);
const double A5_1=-360+240*k*k*4*dx*dx-30*k*k*4*dx*dx*k*k*4*dx*dx, A5_2=-360+240*k*k*dx*dx-30*k*k*dx*dx*k*k*dx*dx;
//二类过渡节点参数
//*************新六点格式
const double B1_1=20*k*k*4*dx*dx-35, B1_2=20*k*k*dx*dx-35;
const double B2=16;
const double B3=5;
const double B4=10;
const double B5=-1;
//三类过渡节点参数
const double C1_1=-16+12*k*k*4*dx*dx-2*k*k*4*dx*dx*k*k*4*dx*dx, C1_2=-16+12*k*k*dx*dx-2*k*k*dx*dx*k*k*dx*dx;
const double C2_1=4-2*k*k*4*dx*dx,    C2_2=4-2*k*k*dx*dx;
const double C3_1=5-2*k*k*4*dx*dx,    C3_2=5-2*k*k*dx*dx;
const double C4=1;

//several repeated type of node
//special 1st kind node
void kind1_node( int j, int i, int i_0, int i_1, int i_2, int i_3, int i_4, int i_5, int i_6, int i_7, 
double w0, double w1, double w2, double w3, double w4, double w5, double w6, double w7, complex** Coef, int** Coef_location )
{
	Coef[j][0].real=w0; Coef[j][1].real=w1; Coef[j][2].real=w2; Coef[j][3].real=w3;
	Coef[j][4].real=w4; Coef[j][5].real=w5; Coef[j][6].real=w6; Coef[j][7].real=w7;
	Coef_location[j][0]=i+i_0; Coef_location[j][1]=i+i_1; Coef_location[j][2]=i+i_2; Coef_location[j][3]=i+i_3;
	Coef_location[j][4]=i+i_4; Coef_location[j][5]=i+i_5; Coef_location[j][6]=i+i_6; Coef_location[j][7]=i+i_7;
};
//special 2nd kind node  二类点，新六点格式
void kind2_node( int j, int i, int i_0, int i_1, int i_2, int i_3, int i_4, int i_5,
	double w0, double w1, double w2, double w3, double w4, double w5, complex** Coef, int** Coef_location )
{
	Coef[j][0].real=w0; Coef[j][1].real=w1; Coef[j][2].real=w2;
	Coef[j][3].real=w3; Coef[j][4].real=w4; Coef[j][5].real=w5; 
	Coef_location[j][0]=i+i_0; Coef_location[j][1]=i+i_1; Coef_location[j][2]=i+i_2;
	Coef_location[j][3]=i+i_3; Coef_location[j][4]=i+i_4; Coef_location[j][5]=i+i_5;
};
//special 3rd kind node
void kind3_node( int j, int i, int i_0, int i_1, int i_2, int i_3, int i_4, int i_5,
				double w0, double w1, double w2, double w3, double w4, double w5, complex** Coef, int** Coef_location )
{
	Coef[j][0].real=w0; Coef[j][1].real=w1; Coef[j][2].real=w2;
	Coef[j][3].real=w3; Coef[j][4].real=w4; Coef[j][5].real=w5; 
	Coef_location[j][0]=i+i_0; Coef_location[j][1]=i+i_1; Coef_location[j][2]=i+i_2;
	Coef_location[j][3]=i+i_3; Coef_location[j][4]=i+i_4; Coef_location[j][5]=i+i_5;
};
//special 5-node
void other_5_node( int j, int i, int i_0, int i_1, int i_2, int i_3, int i_4,
				double w0, double w1, double w2, double w3, double w4, complex** Coef, int** Coef_location )
{
	Coef[j][0].real=w0; Coef[j][1].real=w1; Coef[j][2].real=w2;
	Coef[j][3].real=w3; Coef[j][4].real=w4; 
	Coef_location[j][0]=i+i_0; Coef_location[j][1]=i+i_1; Coef_location[j][2]=i+i_2;
	Coef_location[j][3]=i+i_3; Coef_location[j][4]=i+i_4;
};

// inner node, 5_node_format//j,i--node_num; spacing--node_type; i_0,i_4--diff
void inner_5node( int j, int i, int spacing, int i_0, int i_4, complex** Coef, int** Coef_location )
{
	Coef[j][1].real=Coef[j][3].real
		=Coef[j][0].real=Coef[j][4].real=1;
	Coef[j][2].real=k*k*spacing*spacing*dx*dx-4;// spacing=1
	Coef_location[j][0]=i+i_0;
	Coef_location[j][1]=i-1;
	Coef_location[j][2]=i;
	Coef_location[j][3]=i+1;
	Coef_location[j][4]=i+i_4;
};
void inner_4node( int j, int i, int i_0, int i_1, int i_2, int i_3, 
				 double w0, double w1, double w2, double w3, complex** Coef, int** Coef_location )
{
	Coef[j][0].real=w0; Coef[j][1].real=w1; Coef[j][2].real=w2; Coef[j][3].real=w3;
	Coef_location[j][0]=i+i_0;	Coef_location[j][1]=i+i_1;
	Coef_location[j][2]=i+i_2;	Coef_location[j][3]=i+i_3;
};
void inner_3node( int j, int i, int i_0, int i_1, int i_2, 
				 double w0, double w1, double w2, complex** Coef, int** Coef_location )
{
	Coef[j][0].real=w0; Coef[j][1].real=w1; Coef[j][2].real=w2;
	Coef_location[j][0]=i+i_0;	Coef_location[j][1]=i+i_1;
	Coef_location[j][2]=i+i_2;
};
// boundary node(up and down)//j,i--node_num; type--up or down; i_0,i_4--diff
void boundary_node( int j, int i, int type, int i_0, int i_4, complex** Coef, int** Coef_location )
{
	Coef[j][2].real=(2-2*k*k*16*dx*dx)*cos(k*4*dx)-2;
	Coef[j][2].image=(2-2*k*k*16*dx*dx)*sin(k*4*dx);
	Coef[j][0].real=Coef[j][4].real=1-cos( k*4*dx );
	Coef[j][0].image=Coef[j][4].image=-sin( k*4*dx );
	Coef_location[j][2]=i;
	Coef_location[j][0]=i+i_0;
	Coef_location[j][4]=i+i_4;
	if ( type == 1 )//up
	{
		Coef[j][3].real=2*k*k*16*dx*dx;
		Coef_location[j][3]=i+1;
	} 
	else// 2--down
	{
		Coef[j][1].real=2*k*k*16*dx*dx;
		Coef_location[j][1]=i-1;
	}
};

//*********************************************************************************//
//输入端口区域1参数扫描：n代表重叠区域号，Ey[n][M]代表每次迭代所求的解
//*********************************************************************************//
void scan_coef_Matrix_end1(complex** Coef, int num_total, int ** Coef_location, int sizeM, int sizeN,int offset)
{
	int i = 0, j = 0, m = 0;
	for (j = 0; j<sizeM; j++)         //系数矩阵初始置零
	{
		for (i = 0; i<N_matrix; i++)
		{
			Coef[j][i].real = 0;
			Coef[j][i].image = 0;
			Coef_location[j][i] = -1;//初始化为负值
		}
	}
	//参数扫描
	j = 0;
	for (i = 0; i<Nx_side; i++)//left截断左边界 1阶mur吸收边界条件 均匀网格：4dx
	{
		switch (i)
		{
		case 0:          //角点处理：场平均吸收边条:kz/k!?!??!!?
			Coef[j][0].real = -3;
			Coef[j][2].image = Coef[j][1].image = -sin(k * 4 * dx);
			Coef[j][2].real = Coef[j][1].real = cos(k * 4 * dx);
			Coef[j][3].image = -sin(sqrt(2.0)*k * 4 * dx);
			Coef[j][3].real = cos(sqrt(2.0)*k * 4 * dx);
			Coef_location[j][0] = i;
			Coef_location[j][1] = i + 1;
			Coef_location[j][2] = i + Nx_side;
			Coef_location[j][3] = i + Nx_side + 1;
			break;
		case (Nx_side - 1) :
			Coef[j][1].real = -3;
			Coef[j][3].image = Coef[j][0].image = -sin(k * 4 * dx);
			Coef[j][3].real = Coef[j][0].real = cos(k * 4 * dx);
			Coef[j][2].image = -sin(sqrt(2.0)*k * 4 * dx);
			Coef[j][2].real = cos(sqrt(2.0)*k * 4 * dx);
			Coef_location[j][0] = i - 1;
			Coef_location[j][1] = i;
			Coef_location[j][2] = i + Nx_side - 1;
			Coef_location[j][3] = i + Nx_side;
			break;
		case (N_d0 - 1) :
			Coef[j][2].real = (2 - 2 * k*k * 16 * dx*dx)*cos(k * 4 * dx) - 2;
			Coef[j][2].image = (2 - 2 * k*k * 16 * dx*dx)*sin(k * 4 * dx);
			Coef[j][1].real = 1 - cos(k * 4 * dx);
			Coef[j][1].image = -sin(k * 4 * dx);
			Coef[j][4].real = 2 * k*k * 16 * dx*dx;
			Coef_location[j][2] = i;
			Coef_location[j][1] = i - 1;
			Coef_location[j][4] = i + Nx_side;
			break;
		case (N_d0 + N_a) :
			Coef[j][2].real = (2 - 2 * k*k * 16 * dx*dx)*cos(k * 4 * dx) - 2;
			Coef[j][2].image = (2 - 2 * k*k * 16 * dx*dx)*sin(k * 4 * dx);
			Coef[j][3].real = 1 - cos(k * 4 * dx);
			Coef[j][3].image = -sin(k * 4 * dx);
			Coef[j][4].real = 2 * k*k * 16 * dx*dx;
			Coef_location[j][2] = i;
			Coef_location[j][3] = i + 1;
			Coef_location[j][4] = i + Nx_side;
			break;
		default:
			if ((i >= N_d0) && (i < (N_d0 + N_a)))// add source
			{
				Coef[j][2].real = -1;
				Coef[j][4].image = -sin(kz * 4 * dx);
				Coef[j][4].real = cos(kz * 4 * dx);
				Coef_location[j][2] = i;
				Coef_location[j][4] = i + Nx_side;
			}
			else
			{
				Coef[j][2].real = (2 - 2 * k*k * 16 * dx*dx)*cos(k * 4 * dx) - 2;
				Coef[j][2].image = (2 - 2 * k*k * 16 * dx*dx)*sin(k * 4 * dx);
				Coef[j][1].real = Coef[j][3].real = 1 - cos(k * 4 * dx);
				Coef[j][1].image = Coef[j][3].image = -sin(k * 4 * dx);
				Coef[j][4].real = 2 * k*k * 16 * dx*dx;
				Coef_location[j][2] = i;
				Coef_location[j][1] = i - 1;
				Coef_location[j][3] = i + 1;
				Coef_location[j][4] = i + Nx_side;
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_side); i<(2 * Nx_side); i++)//source node
	{
		switch (i - Nx_side)
		{
		case 0:
			boundary_node(j, j - offset, 1, -Nx_side, Nx_side, Coef, Coef_location);
			break;
		case (Nx_side - 1) :
			boundary_node(j, j - offset, 2, -Nx_side, Nx_side, Coef, Coef_location);
			break;
		case (N_d0 - 1) :
			inner_4node(j, j - offset, -Nx_side, -1, 0, Nx_side, 1, 1, k*k * 16 * dx*dx - 4, 1, Coef, Coef_location);
			break;
		case N_d0:
			Coef[j][3].real = Coef[j][0].real = Coef[j][4].real = 1;
			Coef[j][2].real = k*k * 16 * dx*dx - 4;// spacing=2
			Coef_location[j][0] = i - Nx_side;
			Coef_location[j][2] = i;
			Coef_location[j][3] = i + 1;
			Coef_location[j][4] = i + Nx_side;
			break;
		case (N_d0 + N_a - 1) :
			Coef[j][1].real = Coef[j][0].real = Coef[j][4].real = 1;
			Coef[j][2].real = k*k * 16 * dx*dx - 4;// spacing=2
			Coef_location[j][0] = i - Nx_side;
			Coef_location[j][2] = i;
			Coef_location[j][1] = i - 1;
			Coef_location[j][4] = i + Nx_side;
			break;
		case (N_d0 + N_a) :
			inner_4node(j, j - offset, -Nx_side, 0, 1, Nx_side, 1, k*k * 16 * dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		default:
			if (((i - Nx_side) >= N_d0) && ((i - Nx_side) < (N_d0 + N_a)))// add source
			{
				Coef[j][1].real = Coef[j][3].real
					= Coef[j][0].real = Coef[j][4].real = 1;
				Coef[j][2].real = k*k * 16 * dx*dx - 4;// spacing=1
				Coef_location[j][0] = i - Nx_side;
				Coef_location[j][1] = i - 1;
				Coef_location[j][2] = i;
				Coef_location[j][3] = i + 1;
				Coef_location[j][4] = i + Nx_side;
			}
			else
			{
				inner_5node(j, j - offset, 4, -Nx_side, Nx_side, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (2 * Nx_side); i<(Nx_side*N_d1); i++)//middle node  4dx
	{
		switch (i%Nx_side)
		{
		case 0:
			boundary_node(j, j - offset, 1, -Nx_side, Nx_side, Coef, Coef_location);
			break;
		case (Nx_side - 1) :
			boundary_node(j, j - offset, 2, -Nx_side, Nx_side, Coef, Coef_location);
			break;
		case (N_d0 - 1) :
			inner_4node(j, j - offset, -Nx_side, -1, 0, Nx_side, 1, 1, k*k * 16 * dx*dx - 4, 1, Coef, Coef_location);
			break;
		case N_d0:
			inner_4node(j, j - offset, -Nx_side, 0, 1, Nx_side, 1, k*k * 16 * dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (N_d0 + N_a - 1) :
			inner_4node(j, j - offset, -Nx_side, -1, 0, Nx_side, 1, 1, k*k * 16 * dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (N_d0 + N_a) :
			inner_4node(j, j - offset, -Nx_side, 0, 1, Nx_side, 1, k*k * 16 * dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		default:
			inner_5node(j, j - offset, 4, -Nx_side, Nx_side, Coef, Coef_location);
			break;
		}
		j = j + 1;
	}
	//****************************************************************************************!!!
	for (i = (Nx_side*N_d1); i<(Nx_side*N_d1 + Nx_side); i++)//middle right ,one line
	{
		switch (i - Nx_side*N_d1)
		{
		case 0:
			boundary_node(j, j - offset, 1, -Nx_side, Nx_side + Nx_1, Coef, Coef_location);
			break;
		case (Nx_side - 1) :
			boundary_node(j, j - offset, 2, -Nx_side, Nx_2 + Nx_1, Coef, Coef_location);
			break;
		case (N_d0 - 2) :// 2 lei
			kind2_node(j, j - offset, -Nx_side * 2, -Nx_side, -1, 0, 1, Nx_side - N_d0 + 2,
			B5, B4, B3, B1_1, B3, B2, Coef, Coef_location);
			break;
		case (N_d0 - 1) :// 2 lei缺一点
			other_5_node(j, j - offset, -2 * Nx_side, -Nx_side, -1, 0, (Nx_side - N_d0 + 1 + 2),
			B5, B4, B3, B1_1, B2, Coef, Coef_location);
			break;
		case N_d0:// 2 lei缺一点
			other_5_node(j, j - offset, -2 * Nx_side, -Nx_side, 0, 1, (Nx_side - N_d0 + Nx_1 / 2 - 3),
				B5, B4, B1_1, B3, B2, Coef, Coef_location);
			break;
		case (N_d0 + 1) :// 2 lei
			kind2_node(j, j - offset, -Nx_side * 2, -Nx_side, -1, 0, 1, Nx_side - N_d0 - 1 + Nx_1 / 2 - 1,
			B5, B4, B3, B1_1, B3, B2, Coef, Coef_location);
			break;
		case (N_d0 + N_a - 2) :// 2 lei
			kind2_node(j, j - offset, -Nx_side * 2, -Nx_side, -1, 0, 1, (Nx_side - N_d0 - N_a + 2 + Nx_1 / 2),
			B5, B4, B3, B1_1, B3, B2, Coef, Coef_location);
			break;
		case (N_d0 + N_a - 1) :// 2 lei缺一点
			other_5_node(j, j - offset, -2 * Nx_side, -Nx_side, -1, 0, (Nx_side - N_d0 - N_a + 1 + Nx_1 / 2 + 2),
			B5, B4, B3, B1_1, B2, Coef, Coef_location);
			break;
		case (N_d0 + N_a) :// 2 lei缺一点
			other_5_node(j, j - offset, -2 * Nx_side, -Nx_side, 0, 1, (Nx_side - N_d0 - N_a + Nx_1 - 3),
			B5, B4, B1_1, B3, B2, Coef, Coef_location);
			break;
		case (N_d0 + N_a + 1) :// 2 lei
			kind2_node(j, j - offset, -Nx_side * 2, -Nx_side, -1, 0, 1, (Nx_side - N_d0 - N_a - 1 + Nx_1 - 1),
			B5, B4, B3, B1_1, B3, B2, Coef, Coef_location);
			break;
		default:// 4dx
			if ((i - Nx_side*N_d1) < (N_d0 - 2))//up
			{
				inner_5node(j, j - offset, 4, -Nx_side, Nx_side + Nx_1, Coef, Coef_location);
			}
			else if ((i - Nx_side*N_d1) > (N_d0 + N_a + 1))//down
			{
				inner_5node(j, j - offset, 4, -Nx_side, Nx_2 + Nx_1, Coef, Coef_location);
			}
			else//middle
			{
				inner_5node(j, j - offset, 4, -Nx_side, Nx_side - N_d0 - 1 + Nx_1 + Node_1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//**********************************************************************************************
	for (i = (Nx_side*N_d1 + Nx_side); i<(Nx_side*N_d1 + Nx_side + Nx_1 / 2); i++)//right1 Nx_1: 26 node part1
	{
		switch (i - Nx_side*N_d1 - Nx_side)
		{
		case 0:// 3 lei
			kind3_node(j, j - offset, -(Nx_side - N_d0 + 2) - 1, -(Nx_side - N_d0 + 2), 0, 1, Nx_1 + N_d0 - 3, Nx_1 + N_d0 - 1,
				C4, C3_1, C1_1, C2_1, C4, C3_1, Coef, Coef_location);
			break;
		case (Nx_1 / 2 - 1) ://node2: 3 lei
			kind3_node(j, j - offset, -(Nx_side - N_d0 - 1 + Nx_1 / 2 - 1), -(Nx_side - N_d0 - 1 + Nx_1 / 2 - 1) + 1, -1, 0, Nx_1 + N_d0 - 1, Nx_1 + N_d0 - 1 + 2,
			C3_1, C4, C2_1, C1_1, C3_1, C4, Coef, Coef_location);
			break;
		case 1:// 1 lei
			kind1_node(j, j - offset, -(Nx_side * 2 - N_d0 + 2 + 1), -(Nx_side * 2 - N_d0 + 2 + 1) + 1, -(Nx_side - N_d0 + 2 + 1), -(Nx_side - N_d0 + 2), -1, 0, 1, Nx_1 + N_d0 - 1,
				A2_1, A2_1, A1_1, A1_1, A3_1, A5_1, A3_1, A4_1, Coef, Coef_location);
			break;
		case 2:// 5 node 2dx
			inner_5node(j, j - offset, 2, -(Nx_side - N_d0 + 1 + 2), Nx_1 + N_d0 - 1, Coef, Coef_location);
			break;
		case 3:// 1 lei que 2 dian
			kind2_node(j, j - offset, -(Nx_side - N_d0 + 1 + 3) - Nx_side, -(Nx_side - N_d0 + 1 + 3), -1, 0, 1,
				Nx_1 + N_d0 - 1, A2_1, A1_1, A3_1, A5_1, A3_1, A4_1, Coef, Coef_location);
			break;
		case (Nx_1 / 2 - 4) :// 1 lei que 2 dian
			kind2_node(j, j - offset, -(Nx_side - N_d0 + Nx_1 / 2 - 4) - Nx_side, -(Nx_side - N_d0 + Nx_1 / 2 - 4), -1, 0, 1,
			Nx_1 + N_d0 - 1, A2_1, A1_1, A3_1, A5_1, A3_1, A4_1, Coef, Coef_location);
			break;
		case (Nx_1 / 2 - 3) :// 5 node 2dx
			inner_5node(j, j - offset, 2, -(Nx_side - N_d0 + Nx_1 / 2 - 3), Nx_1 + N_d0 - 1, Coef, Coef_location);
			break;
		case (Nx_1 / 2 - 2) :// 1 lei
			kind1_node(j, j - offset, -(Nx_side * 2 - N_d0 + Nx_1 / 2 - 2), -(Nx_side * 2 - N_d0 + Nx_1 / 2 - 2) + 1, -(Nx_side - N_d0 + Nx_1 / 2 - 2), -(Nx_side - N_d0 + Nx_1 / 2 - 2) + 1,
			-1, 0, 1, Nx_1 + N_d0 - 1,
			A2_1, A2_1, A1_1, A1_1, A3_1, A5_1, A3_1, A4_1, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_side*N_d1 - Nx_side) % 2 == 0)
			{
				inner_4node(j, j - offset, -1, 0, 1, Nx_1 + N_d0 - 1, 1, k*k * 4 * dx*dx - 4, 1, 1, Coef, Coef_location);
			}
			else// 1 lei que 4 dian
			{
				inner_4node(j, j - offset, -1, 0, 1, Nx_1 + N_d0 - 1, A3_1, A5_1, A3_1, A4_1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_side*N_d1 + Nx_side + Nx_1 / 2); i<(Nx_side*N_d1 + Nx_side + Nx_1); i++)// Nx_1 :  26 node! part2
	{
		switch (i - Nx_side*N_d1 - Nx_side - Nx_1 / 2)
		{
		case 0:// 3 lei
			kind3_node(j, j - offset, -(Nx_side - N_d0 - N_a + 2 + Nx_1 / 2) - 1, -(Nx_side - N_d0 - N_a + 2 + Nx_1 / 2), 0, 1, Node_11 - 2, Node_11,
				C4, C3_1, C1_1, C2_1, C4, C3_1, Coef, Coef_location);
			break;
		case (Nx_1 / 2 - 1) :// 3 lei
			kind3_node(j, j - offset, -(Nx_side - N_d0 - N_a + Nx_1 - 2), -(Nx_side - N_d0 - N_a + Nx_1 - 2) + 1, -1, 0, Node_11, Node_11 + 2,
			C3_1, C4, C2_1, C1_1, C3_1, C4, Coef, Coef_location);
			break;
		case 1:// 1 lei
			kind1_node(j, j - offset, -(Nx_side * 2 - N_d0 - N_a + 2 + Nx_1 / 2) - 1, -(Nx_side * 2 - N_d0 - N_a + 2 + Nx_1 / 2), -(Nx_side - N_d0 - N_a + 2 + Nx_1 / 2) - 1, -(Nx_side - N_d0 - N_a + 2 + Nx_1 / 2),
				-1, 0, 1, Node_11,
				A2_1, A2_1, A1_1, A1_1, A3_1, A5_1, A3_1, A4_1, Coef, Coef_location);
			break;
		case 2:// 5 node 2dx
			inner_5node(j, j - offset, 2, -(Nx_side - N_d0 - N_a + 1 + Nx_1 / 2 + 2), Node_11, Coef, Coef_location);
			break;
		case 3:// 1 lei que 2 dian
			kind2_node(j, j - offset, -(Nx_side - N_d0 - N_a + 1 + Nx_1 / 2 + 3) - Nx_side, -(Nx_side - N_d0 - N_a + 1 + Nx_1 / 2 + 3), -1, 0, 1,
				Node_11, A2_1, A1_1, A3_1, A5_1, A3_1, A4_1, Coef, Coef_location);
			break;
		case (Nx_1 / 2 - 4) :// 1 lei que 2 dian
			kind2_node(j, j - offset, -(Nx_side - N_d0 - N_a + Nx_1 - 4) - Nx_side, -(Nx_side - N_d0 - N_a + Nx_1 - 4), -1, 0, 1, Node_11,
			A2_1, A1_1, A3_1, A5_1, A3_1, A4_1, Coef, Coef_location);
			break;
		case (Nx_1 / 2 - 3) :// 5 node 2dx
			inner_5node(j, j - offset, 2, -(Nx_side - N_d0 - N_a + Nx_1 - 3), Node_11, Coef, Coef_location);
			break;
		case (Nx_1 / 2 - 2) :// 1 lei
			kind1_node(j, j - offset, -(Nx_side * 2 - N_d0 - N_a + Nx_1 - 2), -(Nx_side * 2 - N_d0 - N_a + Nx_1 - 2) + 1, -(Nx_side - N_d0 - N_a + Nx_1 - 2), -(Nx_side - N_d0 - N_a + Nx_1 - 2) + 1,
			-1, 0, 1, Node_11, A2_1, A2_1, A1_1, A1_1, A3_1, A5_1, A3_1, A4_1, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_side*N_d1 - Nx_side - Nx_1 / 2) % 2 == 0)
			{
				inner_4node(j, j - offset, -1, 0, 1, Node_11, 1, k*k * 4 * dx*dx - 4, 1, 1, Coef, Coef_location);
			}
			else//1 lei
			{
				inner_4node(j, j - offset, -1, 0, 1, Node_11, A3_1, A5_1, A3_1, A4_1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//***************************************************************************************
	for (i = (Nx_side*N_d1 + Nx_side + Nx_1); i<(Nx_side*N_d1 + Nx_side + Nx_1 + Nx_2); i++)//right 2 Nx_2 :  74 node! 
	{
		switch (i - Nx_side*N_d1 - Nx_1 - Nx_side)
		{
		case 0:
			boundary_node(j, j - offset, 1, -Nx_side - Nx_1, Nx_2 + Nx_3, Coef, Coef_location);
			break;
		case (Nx_2 - 1) :
			boundary_node(j, j - offset, 2, -Nx_1 - Nx_2, Nx_2 + Nx_3, Coef, Coef_location);
			break;
		case (N_d0 - 3) ://2 lei shang
			kind2_node(j, j - offset, -Nx_side - Nx_1, -2, -1, 0, 1, Nx_2 + Nx_3,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (N_d0 - 2) ://3 lei shang
			kind3_node(j, j - offset, -Nx_side - Nx_1 - 1, -Nx_side - Nx_1, -1, 0, 1, Nx_2 - N_d0 + 2,
			C4, C4, C3_1, C1_1, C3_1, C2_1, Coef, Coef_location);
			break;
		case Node_1://3 lei xia
			kind3_node(j, j - offset, -(Nx_side - N_d0 - 1 + Nx_1 + Node_1), -(Nx_side - N_d0 - 1 + Nx_1 + Node_1) + 1, -1, 0, 1, Nx_2 - N_d0 + 2,
				C4, C4, C3_1, C1_1, C3_1, C2_1, Coef, Coef_location);
			break;
		case (Node_1 + 1) ://2 lei xia
			kind2_node(j, j - offset, -(Nx_side - N_d0 - 1 + Nx_1 + Node_1), -1, 0, 1, 2, Nx_2 + Nx_3,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_1 + N_a - 4) ://2 lei shang
			kind2_node(j, j - offset, -(Nx_side - N_d0 - 1 + Nx_1 + Node_1), -2, -1, 0, 1, Nx_2 + Nx_3,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_1 + N_a - 3) ://3 lei shang
			kind3_node(j, j - offset, -(Nx_side - N_d0 - 1 + Nx_1 + Node_1) - 1, -(Nx_side - N_d0 - 1 + Nx_1 + Node_1), -1, 0, 1, Nx_2 - Node_11 + Nx_3 - 1,
			C4, C4, C3_1, C1_1, C3_1, C2_1, Coef, Coef_location);
			break;
		case Node_11://3 lei xia
			kind3_node(j, j - offset, -Nx_1 - Nx_2, -Nx_1 - Nx_2 + 1, -1, 0, 1, Nx_2 - Node_11 + Nx_3 - 1,
				C4, C4, C3_1, C1_1, C3_1, C2_1, Coef, Coef_location);
			break;
		case (Node_11 + 1) ://2 lei xia
			kind2_node(j, j - offset, -Nx_1 - Nx_2, -1, 0, 1, 2, Nx_2 + Nx_3,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_side*N_d1 - Nx_1 - Nx_side) < (N_d0 - 3))// 4dx
			{
				inner_5node(j, j - offset, 4, -Nx_side - Nx_1, Nx_2 + Nx_3, Coef, Coef_location);
			}
			else if (((i - Nx_side*N_d1 - Nx_1 - Nx_side) >= (N_d0 - 1)) && ((i - Nx_side*N_d1 - Nx_1 - Nx_side) < Node_1))// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Nx_1 + N_d0 - 1), Nx_2 - N_d0 + 2, Coef, Coef_location);
			}
			else if (((i - Nx_side*N_d1 - Nx_1 - Nx_side) > (Node_1 + 1)) && ((i - Nx_side*N_d1 - Nx_1 - Nx_side) < (Node_1 + N_a - 4)))// 4dx!
			{
				inner_5node(j, j - offset, 4, -(Nx_side - N_d0 - 1 + Nx_1 + Node_1), Nx_2 + Nx_3, Coef, Coef_location);
			}
			else if (((i - Nx_side*N_d1 - Nx_1 - Nx_side) >= (Node_1 + N_a - 2)) && ((i - Nx_side*N_d1 - Nx_1 - Nx_side) < Node_11))// 2dx!
			{
				inner_5node(j, j - offset, 2, -Node_11, Nx_2 - Node_11 + Nx_3 - 1, Coef, Coef_location);
			}
			else// 4dx!
			{
				inner_5node(j, j - offset, 4, -Nx_1 - Nx_2, Nx_2 + Nx_3, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//***************************************************************************************
	//cout << "right 3:" << j << endl;
	for (i = (Nx_side*N_d1 + Nx_side + Nx_1 + Nx_2); i<(Nx_side*N_d1 + Nx_side + Nx_1 + Nx_2 + Nx_3); i++)// right 3 :  30 node! 
	{
		switch (i - Nx_side*N_d1 - Nx_side - Nx_1 - Nx_2)
		{
		case 0:// 1 lei shang
			kind1_node(j, j - offset, -(Nx_2 - N_d0 + 2) - 2, -(Nx_2 - N_d0 + 2) - 1, -(Nx_2 - N_d0 + 2), 0, 1, Nx_3 + N_d0 - 2 - 2, Nx_3 + N_d0 - 2 - 1, Nx_3 + N_d0 - 2,
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case (Nx_3 / 2 - 1) :// 1 lei xia
			kind1_node(j, j - offset, -(Nx_2 - N_d0 + 2), -(Nx_2 - N_d0 + 2) + 1, -(Nx_2 - N_d0 + 2) + 2, -1, 0,
			Nx_3 + N_d0 - 2, Nx_3 + N_d0 - 2 + 1, Nx_3 + N_d0 - 2 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		case (Nx_3 / 2) :// 1 lei shang
			kind1_node(j, j - offset, -(Nx_2 - Node_11 - 1 + Nx_3) - 2, -(Nx_2 - Node_11 - 1 + Nx_3) - 1, -(Nx_2 - Node_11 - 1 + Nx_3), 0, 1,
			1 + Node_11 - 2, 1 + Node_11 - 1, 1 + Node_11,
			A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case (Nx_3 - 1) :// 1 lei xia
			kind1_node(j, j - offset, -(Nx_2 - Node_11 - 1 + Nx_3), -(Nx_2 - Node_11 - 1 + Nx_3) + 1, -(Nx_2 - Node_11 - 1 + Nx_3) + 2, -1, 0,
			1 + Node_11, 1 + Node_11 + 1, 1 + Node_11 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_side*N_d1 - Nx_side - Nx_1 - Nx_2) < (Nx_3 / 2 - 1))
			{
				inner_5node(j, j - offset, 2, -(Nx_2 - N_d0 + 2), Nx_3 + N_d0 - 2, Coef, Coef_location);
			}
			else
			{
				inner_5node(j, j - offset, 2, -(Nx_2 - Node_11 - 1 + Nx_3), 1 + Node_11, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//***************************************************************************************
	//最后一行 本节点向下的区域没有影响
	//cout << "circle num1:" << j << endl;
	for (i = (Nx_2 + Nx_3); i<(Nx_2 + Nx_3 + Node_1 + 2); i++)// circle num1 Nx_2:  74 node!  part1 !
	{
		switch (i - Nx_2 - Nx_3)
		{
		case 0:
			boundary_node(j, j - offset, 1, -Nx_2 - Nx_3, Num1 + Num2 + Num3 + Num4, Coef, Coef_location);
			break;
		case (N_d0 - 3) ://2 lei shang 1j
			kind2_node(j, j - offset, -Nx_2 - Nx_3, -2, -1, 0, 1, Num1 + Num2 + Num3 + Num4,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_1 + 1) ://2 lei xia 1j
			kind2_node(j, j - offset, -(Nx_2 + Nx_3), -1, 0, 1, 2, Num1 - Node_1 + Num2 + Num3 + Num4 + Node_5,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3) < (N_d0 - 3))// 4dx!
			{
				inner_5node(j, j - offset, 4, -Nx_3 - Nx_2, Num1 + Num2 + Num3 + Num4, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3) > (N_d0 - 3)) && ((i - Nx_2 - Nx_3) < Node11))// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Nx_3 + N_d0 - 2), Num1 - N_d0 + 2 + Num2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3) >= Node11) && ((i - Nx_2 - Nx_3) <= Node12))// 2 lei left 2j !
			{
				kind2_node(j, j - offset, -Nx_3 - Nx_2, -(Nx_3 + N_d0 - 2), -1, 0, 1, Num1 - (i - Nx_2 - Nx_3) + ((i - Nx_2 - Nx_3) - Node11) * 2,
					B5, B4, B3, B1_2, B3, B2, Coef, Coef_location);
			}
			else// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Nx_3 + N_d0 - 2), Num1 - Node_1 + Num2 + Num3 / 2 - 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Node_1 + 2); i<(Nx_2 + Nx_3 + Num1); i++)// circle num1 Nx_2:  74 node!  part2 !
	{
		switch (i - Nx_2 - Nx_3)
		{
		case (Num1 - 1) :
			boundary_node(j, j - offset, 2, -Nx_2 - Nx_3, Num5 + Num2 + Num3 + Num4, Coef, Coef_location);
			break;
		case (Node_1 + N_a - 4) ://2 lei shang 1j
			kind2_node(j, j - offset, -Nx_2 - Nx_3, -2, -1, 0, 1, Num1 - Node_1 + Num2 + Num3 + Num4 + Node_5,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_11 + 1) ://2 lei xia 1j
			kind2_node(j, j - offset, -(Nx_2 + Nx_3), -1, 0, 1, 2, Num5 + Num2 + Num3 + Num4,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3) < (Node_1 + N_a - 4))// 4dx!
			{
				inner_5node(j, j - offset, 4, -Nx_3 - Nx_2, Num1 - Node_1 + Num2 + Num3 + Num4 + Node_5, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3) > (Node_1 + N_a - 4)) && ((i - Nx_2 - Nx_3) < Node13))// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Node_11 + 1), Num1 - Node_1 - N_a + 3 + Num2 + Num3 / 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3) >= Node13) && ((i - Nx_2 - Nx_3) <= Node14))// 2 lei left 2j !
			{
				kind2_node(j, j - offset, -Nx_3 - Nx_2, -(Node_11 + 1), -1, 0, 1, Num1 - (i - Nx_2 - Nx_3) + Num2 / 2 + ((i - Nx_2 - Nx_3) - Node13) * 2,
					B5, B4, B3, B1_2, B3, B2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3) >= Node14) && ((i - Nx_2 - Nx_3) <= Node_11))// 2dx
			{
				inner_5node(j, j - offset, 2, -(Node_11 + 1), Num1 - Node_11 + Num2 + Num3 - 1, Coef, Coef_location);
			}
			else
			{
				inner_5node(j, j - offset, 4, -(Nx_3 + Nx_2), Num5 + Num2 + Num3 + Num4, Coef, Coef_location);// 4dx!
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!!!!!!!!!
	//cout << "circle num 2:" << j << endl;
	for (i = (Nx_2 + Nx_3 + Num1); i<(Nx_2 + Nx_3 + Num1 + Num2 / 2); i++)// num 2 :  26 node! part1
	{
		switch (i - Nx_2 - Nx_3 - Num1)
		{
		case 0:// 3 lei 2j shang
			kind3_node(j, j - offset, -(Num1 - Node11) - 1, -(Num1 - Node11), 0, 1, Num2 + 5 - 2, Num2 + 5,
				C4, C3_2, C1_2, C2_2, C4, C3_2, Coef, Coef_location);
			break;
		case (Num2 / 2 - 1) ://node2: 3 lei 2j xia
			kind3_node(j, j - offset, -(Num1 - Node12 + Num2 / 2 - 1), -(Num1 - Node12 + Num2 / 2 - 1) + 1, -1, 0, Num2 + 5, Num2 + 5 + 2,
			C3_2, C4, C2_2, C1_2, C3_2, C4, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1) % 2 == 0)// 1dx
			{
				inner_5node(j, j - offset, 1, -((i - Nx_2 - Nx_3 - Num1) + Num1 - Node11 - (i - Nx_2 - Nx_3 - Num1) / 2), Num2 + 5, Coef, Coef_location);
			}
			else//1 lei  2j left
			{
				kind1_node(j, j - offset, -((i - Nx_2 - Nx_3 - Num1) + Num1 - Node11 - (i - Nx_2 - Nx_3 - Num1) / 2) - (Nx_3 + N_d0 - 2),
					-((i - Nx_2 - Nx_3 - Num1) + Num1 - Node11 - (i - Nx_2 - Nx_3 - Num1) / 2) - (Nx_3 + N_d0 - 2) + 1,
					-((i - Nx_2 - Nx_3 - Num1) + Num1 - Node11 - (i - Nx_2 - Nx_3 - Num1) / 2),
					-((i - Nx_2 - Nx_3 - Num1) + Num1 - Node11 - (i - Nx_2 - Nx_3 - Num1) / 2) + 1,
					-1, 0, 1, 5 + Num2, A2_2, A2_2, A1_2, A1_2, A3_2, A5_2, A3_2, A4_2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 / 2); i<(Nx_2 + Nx_3 + Num1 + Num2); i++)// num 2 :  26 node! part2
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 / 2)
		{
		case 0:// 3 lei 2j shang
			kind3_node(j, j - offset, -(Num1 - Node13 + Num2 / 2) - 1, -(Num1 - Node13 + Num2 / 2), 0, 1, Num3 - 5 - 2, Num3 - 5,
				C4, C3_2, C1_2, C2_2, C4, C3_2, Coef, Coef_location);
			break;
		case (Num2 / 2 - 1) :// 3 lei 2j xia
			kind3_node(j, j - offset, -(Num1 - Node14 + Num2 - 1), -(Num1 - Node14 + Num2 - 1) + 1, -1, 0, Num3 - 5, Num3 - 5 + 2,
			C3_2, C4, C2_2, C1_2, C3_2, C4, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 / 2) % 2 == 0)// 1dx
			{
				inner_5node(j, j - offset, 1, -((i - Nx_2 - Nx_3 - Num1 - Num2 / 2) + Num1 - Node13 - (i - Nx_2 - Nx_3 - Num1 - Num2 / 2) / 2) - Num2 / 2,
					Num3 - 5, Coef, Coef_location);
			}
			else//1 lei 2j left
			{
				kind1_node(j, j - offset, -((i - Nx_2 - Nx_3 - Num1 - Num2 / 2) + Num1 - Node13 - (i - Nx_2 - Nx_3 - Num1 - Num2 / 2) / 2) - Num2 / 2 - (Node_11 + 1),
					-((i - Nx_2 - Nx_3 - Num1 - Num2 / 2) + Num1 - Node13 - (i - Nx_2 - Nx_3 - Num1 - Num2 / 2) / 2) - Num2 / 2 - (Node_11 + 1) + 1,
					-((i - Nx_2 - Nx_3 - Num1 - Num2 / 2) + Num1 - Node13 - (i - Nx_2 - Nx_3 - Num1 - Num2 / 2) / 2) - Num2 / 2,
					-((i - Nx_2 - Nx_3 - Num1 - Num2 / 2) + Num1 - Node13 - (i - Nx_2 - Nx_3 - Num1 - Num2 / 2) / 2) - Num2 / 2 + 1,
					-1, 0, 1, Num3 - 5, A2_2, A2_2, A1_2, A1_2, A3_2, A5_2, A3_2, A4_2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!!!
	//cout << "circle num 3:" << j << endl;
	for (i = (Nx_2 + Nx_3 + Num1 + Num2); i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 / 2); i++)// num 3 :  46 node! part1
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2)
		{
		case 0:// 1 lei 1j shang
			kind1_node(j, j - offset, -(Num1 - N_d0 + 2 + Num2) - 2, -(Num1 - N_d0 + 2 + Num2) - 1, -(Num1 - N_d0 + 2 + Num2), 0, 1, Num3 + Num4 + N_d0 - 2 - 2,
				Num3 + Num4 + N_d0 - 2 - 1, Num3 + Num4 + N_d0 - 2, A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case 3:// 2 lei 2j shang
			kind2_node(j, j - offset, -(Num1 - N_d0 + 2 + Num2), -2, -1, 0, 1, Num3 + Num4 + N_d0 - 2,
				B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case 4:// 3 lei 2j shang
			kind3_node(j, j - offset, -(Num1 - N_d0 + 2 + Num2) - 1, -(Num1 - N_d0 + 2 + Num2), -1, 0, 1, Num3 - 4,
				C4, C4, C3_2, C1_2, C3_2, C2_2, Coef, Coef_location);
			break;
		case (Num3 / 2 - 5) :// 3 lei 2j xia
			kind3_node(j, j - offset, -(Num1 - Node_1 + Num2 + Num3 / 2 - 1), -(Num1 - Node_1 + Num2 + Num3 / 2 - 1) + 1, -1, 0, 1, Num3 - 4,
			C4, C4, C3_2, C1_2, C3_2, C2_2, Coef, Coef_location);
			break;
		case (Num3 / 2 - 4) :// 2 lei 2j xia
			kind2_node(j, j - offset, -(Num1 - Node_1 + Num2 + Num3 / 2 - 1), -1, 0, 1, 2, Num3 / 2 + 1 + Num4 + Node_5,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num3 / 2 - 1) :// 1 lei 1j xia
			kind1_node(j, j - offset, -(Num1 - Node_1 + Num2 + Num3 / 2 - 1), -(Num1 - Node_1 + Num2 + Num3 / 2 - 1) + 1, -(Num1 - Node_1 + Num2 + Num3 / 2 - 1) + 2, -1, 0,
			Num3 / 2 + 1 + Num4 + Node_5, Num3 / 2 + 1 + Num4 + Node_5 + 1, Num3 / 2 + 1 + Num4 + Node_5 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2) < 3)// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num1 - N_d0 + 2 + Num2), Num3 + Num4 + N_d0 - 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2) > 4) && ((i - Nx_2 - Nx_3 - Num1 - Num2) < (Num3 / 2 - 5)))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num2 + 5), Num3 - 4, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num1 - Node_1 + Num2 + Num3 / 2 - 1), Num3 / 2 + 1 + Num4 + Node_5, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 / 2); i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3); i++)// num 3 :  46 node! part2
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 / 2)
		{
		case 0:// 1 lei 1j shang
			kind1_node(j, j - offset, -(Num1 - Node_1 - N_a + 3 + Num2 + Num3 / 2) - 2, -(Num1 - Node_1 - N_a + 3 + Num2 + Num3 / 2) - 1, -(Num1 - Node_1 - N_a + 3 + Num2 + Num3 / 2), 0, 1,
				Num3 / 2 + Num4 + Node_5 + N_a - 3 - 2, Num3 / 2 + Num4 + Node_5 + N_a - 3 - 1, Num3 / 2 + Num4 + Node_5 + N_a - 3,
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case 3:// 2 lei 2j shang
			kind2_node(j, j - offset, -(Num1 - Node_1 - N_a + 3 + Num2 + Num3 / 2), -2, -1, 0, 1, Num3 / 2 + Num4 + Node_5 + N_a - 3,
				B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case 4:// 3 lei 2j shang
			kind3_node(j, j - offset, -(Num1 - Node_1 - N_a + 3 + Num2 + Num3 / 2) - 1, -(Num1 - Node_1 - N_a + 3 + Num2 + Num3 / 2), -1, 0, 1, 5 + Num4 - 1,
				C4, C4, C3_2, C1_2, C3_2, C2_2, Coef, Coef_location);
			break;
		case (Num3 / 2 - 5) :// 3 lei 2j xia
			kind3_node(j, j - offset, -(Num1 - Node_11 + Num2 + Num3 - 1), -(Num1 - Node_11 + Num2 + Num3 - 1) + 1, -1, 0, 1, 5 + Num4 - 1,
			C4, C4, C3_2, C1_2, C3_2, C2_2, Coef, Coef_location);
			break;
		case (Num3 / 2 - 4) :// 2 lei 2j xia
			kind2_node(j, j - offset, -(Num1 - Node_11 + Num2 + Num3 - 1), -1, 0, 1, 2, 1 + Num4 + Node_55,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num3 / 2 - 1) :// 1 lei 1j xia
			kind1_node(j, j - offset, -(Num1 - Node_11 + Num2 + Num3 - 1), -(Num1 - Node_11 + Num2 + Num3 - 1) + 1, -(Num1 - Node_11 + Num2 + Num3 - 1) + 2, -1, 0,
			1 + Num4 + Node_55, 1 + Num4 + Node_55 + 1, 1 + Num4 + Node_55 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 / 2) < 3)// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num1 - Node_1 - N_a + 3 + Num2 + Num3 / 2), Num3 / 2 + Num4 + Node_5 + N_a - 3, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 / 2) > 4) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 / 2) < (Num3 / 2 - 5)))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num3 - 5), 5 + Num4 - 1, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num1 - Node_11 + Num2 + Num3 - 1), 1 + Num4 + Node_55, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	//cout << "circle num 4:" << j << endl;
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3); i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 / 2); i++)// num 4 :  30 node! part1
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3)
		{
		case 0:// 1 lei 2j shang
			kind1_node(j, j - offset, -(Num3 - 4) - 2, -(Num3 - 4) - 1, -(Num3 - 4), 0, 1, Num4 + N_d0 + 2 - 2, Num4 + N_d0 + 2 - 1, Num4 + N_d0 + 2,
				A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case (Num4 / 2 - 1) :// 1 lei 2j xia
			kind1_node(j, j - offset, -(Num3 - 4), -(Num3 - 4) + 1, -(Num3 - 4) + 2, -1, 0,
			Num4 + N_d0 + 2 - N_c1, Num4 + N_d0 + 2 - N_c1 + 1, Num4 + N_d0 + 2 - N_c1 + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3) <= Node41)// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num3 - 4), Num4 + N_d0 + 2, Coef, Coef_location);
			}
			else if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3) >= Node42) // 1dx
			{
				inner_5node(j, j - offset, 1, -(Num3 - 4), Num4 + N_d0 + 2 - N_c1, Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_4node(j, j - offset, -(Num3 - 4), -1, 0, 1, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 / 2); i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4); i++)// num 4 :  30 node! part2
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3)
		{
		case (Num4 / 2) :// 1 lei 2j shang
			kind1_node(j, j - offset, -(5 + Num4 - 1) - 2, -(5 + Num4 - 1) - 1, -(5 + Num4 - 1), 0, 1,
			(Node_55 - 4 + 1 + N_c1) - 2, (Node_55 - 4 + 1 + N_c1) - 1, (Node_55 - 4 + 1 + N_c1),
			A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case (Num4 - 1) :// 1 lei 2j xia
			kind1_node(j, j - offset, -(5 + Num4 - 1), -(5 + Num4 - 1) + 1, -(5 + Num4 - 1) + 2, -1, 0,
			(Node_55 - 4 + 1), (Node_55 - 4 + 1) + 1, (Node_55 - 4 + 1) + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3) <= Node43)// 1dx
			{
				inner_5node(j, j - offset, 1, -(5 + Num4 - 1), (Node_55 - 4 + 1 + N_c1), Coef, Coef_location);
			}
			else if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3) >= Node44) // 1dx
			{
				inner_5node(j, j - offset, 1, -(5 + Num4 - 1), (Node_55 - 4 + 1), Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_4node(j, j - offset, -(5 + Num4 - 1), -1, 0, 1, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	//cout << "circle num5:" << j << endl;
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4); i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Node_5 + 2); i++)// circle num5 :  80 node!  part1 !
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4)
		{
		case 0:
			boundary_node(j, j - offset, 1, -Num1 - Num2 - Num3 - Num4, Num5 + Num6 + Num7 + Num8, Coef, Coef_location);
			break;
		case (N_d0 - 3) ://2 lei shang 1j
			kind2_node(j, j - offset, -Num1 - Num2 - Num3 - Num4, -2, -1, 0, 1, Num5 + Num6 + Num7 + Num8,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (N_d0 + 1) :// 2 lei 2j shang
			kind2_node(j, j - offset, -(Num3 + Num4 + N_d0 - 2), -2, -1, 0, 1, Num5 - N_d0 + 2 + Num6,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node51:// inner 3 node
			inner_3node(j, j - offset, -(Num4 + N_d0 + 2), -1, 0, 1, 1, k*k*dx*dx - 4, Coef, Coef_location);
			break;
		case (Node51 + 1) :// inner 3 node
			inner_3node(j, j - offset, -(Num4 + N_d0 + 2) + N_c1, 0, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node_5 - 3) ://2 lei xia 2j
			kind2_node(j, j - offset, -(Num3 / 2 + 1 + Num4 + Node_5), -1, 0, 1, 2, Num5 - Node_5 + Num6 + Num7 / 2 - 1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_5 + 1) ://2 lei xia 1j
			kind2_node(j, j - offset, -(Num1 - Node_1 + Num2 + Num3 + Num4 + Node_5), -1, 0, 1, 2, Num5 - Node_5 + Num6 + Num7 + Num8 + Node_9,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < (N_d0 - 3))// 4dx!
			{
				inner_5node(j, j - offset, 4, -Num1 - Num2 - Num3 - Num4, Num5 + Num6 + Num7 + Num8, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) > (N_d0 - 3)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < (N_d0 + 1)))// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Num3 + Num4 + N_d0 - 2), Num5 - N_d0 + 2 + Num6, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) > (N_d0 + 1)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < Node51))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num4 + N_d0 + 2), Num5 - N_d0 - 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) > (Node51 + 1)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < (Node_5 - 3)))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num4 + N_d0 + 2) + N_c1, Num5 - N_d0 - 2 - 2, Coef, Coef_location);
			}
			else// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Num3 / 2 + 1 + Num4 + Node_5), Num5 - Node_5 + Num6 + Num7 / 2 - 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Node_5 + 2); i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5); i++)// circle num5 :  80 node!  part2 !
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4)
		{
		case (Num5 - 1) :
			boundary_node(j, j - offset, 2, -Num2 - Num3 - Num4 - Num5, Num6 + Num7 + Num8 + Num9, Coef, Coef_location);
			break;
		case (Node_5 + N_a - 4) ://2 lei shang 1j
			kind2_node(j, j - offset, -(Num1 - Node_1 + Num2 + Num3 + Num4 + Node_5), -2, -1, 0, 1, Num5 - Node_5 + Num6 + Num7 + Num8 + Node_9,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_5 + N_a) :// 2 lei 2j shang
			kind2_node(j, j - offset, -(Num3 / 2 + Num4 + Node_5 + N_a - 3), -2, -1, 0, 1, Num5 - Node_5 - N_a + 3 + Num6 + Num7 / 2,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node52:// inner 3 node
			inner_3node(j, j - offset, -(Node_55 - 4 + 1) - N_c1, -1, 0, 1, 1, k*k*dx*dx - 4, Coef, Coef_location);
			break;
		case (Node52 + 1) :// inner 3 node
			inner_3node(j, j - offset, -(Node_55 - 4 + 1), 0, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node_55 - 3) ://2 lei xia 2j
			kind2_node(j, j - offset, -(Num4 + Node_55 + 1), -1, 0, 1, 2, Num5 - Node_55 + Num6 + Num7 - 1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_55 + 1) ://2 lei xia 1j
			kind2_node(j, j - offset, -Num2 - Num3 - Num4 - Num5, -1, 0, 1, 2, Num6 + Num7 + Num8 + Num9,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < (Node_5 + N_a - 4))// 4dx!
			{
				inner_5node(j, j - offset, 4, -(Num1 - Node_1 + Num2 + Num3 + Num4 + Node_5), Num5 - Node_5 + Num6 + Num7 + Num8 + Node_9, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) > (Node_5 + N_a - 4)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < (Node_5 + N_a)))// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Num3 / 2 + Num4 + Node_5 + N_a - 3), Num5 - Node_5 - N_a + 3 + Num6 + Num7 / 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) > (Node_5 + N_a)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < Node52))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Node_55 - 4 + 1) - N_c1, (Num5 - Node_55 + 4 + Num6 - 1) + 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) > (Node52 + 1)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < (Node_55 - 3)))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Node_55 - 4 + 1), (Num5 - Node_55 + 4 + Num6 - 1), Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) > (Node_55 - 3)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < (Node_55 + 1)))// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num4 + Node_55 + 1), Num5 - Node_55 + Num6 + Num7 - 1, Coef, Coef_location);
			}
			else// 4dx!
			{
				inner_5node(j, j - offset, 4, -Num2 - Num3 - Num4 - Num5, Num6 + Num7 + Num8 + Num9, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	cout << "end1区域节点数：" << j << endl;
	if (j == sizeM)
		cout << "end1 passed..." << endl;
	else
		cout << "end1 failed..." << endl;
	//修正矩阵，将不在本区域的节点刨除
	for (int i = 0; i != sizeM; i++)
		for (int j = 0; j != N_matrix; j++){
		if (Coef_location[i][j] < 0 || Coef_location[i][j] >= sizeN)
			Coef_location[i][j] = -1;
		}
	if (Debug)
	{
		ofstream outFile_Coef("end1.txt");
		if (!outFile_Coef)
		{
			cout << "Can't write the coefficient matrix into file !" << endl;
			exit(1);//Terminate with fault
		}
		else
		{
			int i, j;
			for (i = 0; i < sizeM; i++)
			{
				for (j = 0; j < N_matrix; j++)
				{
					if ((Coef[i][j].real != 0 || Coef[i][j].image != 0) && Coef_location[i][j] >= 0)
						outFile_Coef << Coef[i][j].real << "+j*" << Coef[i][j].image << " (" << i << "," << Coef_location[i][j] << ") ";
				}
				outFile_Coef << endl << endl;
			}
		}
		outFile_Coef.close();
	}

};
//*********************************************************************************//
//输出端口区域2参数扫描：n代表重叠区域号，Ey[n][M]代表每次迭代所求的解
//*********************************************************************************//
void scan_coef_Matrix_end2(complex** Coef, int num_total, int ** Coef_location, int sizeM, int sizeN, int offset)
{
	int i = 0, j = 0, m = 0;
	for (j = 0; j<sizeM; j++)         //系数矩阵初始置零
	{
		for (i = 0; i<N_matrix; i++)
		{
			Coef[j][i].real = 0;
			Coef[j][i].image = 0;
			Coef_location[j][i] = -1;//初始化为负值
		}
	}
	//参数扫描
	j = 0;
	//******************************************************************************************
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9);
		i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Node_5 + 2); i++)// circle num13(5) :  80 node!  part1 !
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9)
		{
		case 0:
			boundary_node(j, j - offset, 1, -Num9 - Num8 - Num7 - Num6, Num5 + Num4 + Num3 + Num2, Coef, Coef_location);
			break;
		case (N_d0 - 3) ://2 lei shang 1j
			kind2_node(j, j - offset, -Num9 - Num8 - Num7 - Num6, -2, -1, 0, 1, Num5 + Num4 + Num3 + Num2,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (N_d0 + 1) :// 2 lei 2j shang
			kind2_node(j, j - offset, -(Num7 + Num6 + N_d0 - 2), -2, -1, 0, 1, Num5 - N_d0 + 2 + Num4,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node51:// inner 3 node
			inner_3node(j, j - offset, -1, 0, (Num5 - N_d0 - 2), 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node51 + 1) :// inner 3 node
			inner_3node(j, j - offset, 0, 1, (Num5 - N_d0 - 2) + N_c1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Node_5 - 3) ://2 lei xia 2j
			kind2_node(j, j - offset, -(Num7 / 2 + 1 + Num6 + Node_5), -1, 0, 1, 2, Num5 - Node_5 + Num4 + Num3 / 2 - 1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_5 + 1) ://2 lei xia 1j
			kind2_node(j, j - offset, -(Num9 - Node_9 + Num8 + Num7 + Num6 + Node_5), -1, 0, 1, 2, Num5 - Node_5 + Num4 + Num3 + Num2 + Node_1,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < (N_d0 - 3))// 4dx!
			{
				inner_5node(j, j - offset, 4, -Num9 - Num8 - Num7 - Num6, Num5 + Num4 + Num3 + Num2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) > (N_d0 - 3)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < (N_d0 + 1)))// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Num7 + Num6 + N_d0 - 2), Num5 - N_d0 + 2 + Num4, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) > (N_d0 + 1)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < Node51))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num6 + N_d0 + 2), (Num5 - N_d0 - 2), Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) > (Node51 + 1)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < (Node_5 - 3)))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num6 + N_d0 + 2) - 2, (Num5 - N_d0 - 2) + N_c1, Coef, Coef_location);
			}
			else// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Num7 / 2 + 1 + Num6 + Node_5), Num5 - Node_5 + Num4 + Num3 / 2 - 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Node_5 + 2);
		i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9); i++)// circle num13(5) :  80 node!  part2 !
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9)
		{
		case (Num5 - 1) :
			boundary_node(j, j - offset, 2, -Num8 - Num7 - Num6 - Num5, Num4 + Num3 + Num2 + Num1, Coef, Coef_location);
			break;
		case (Node_5 + N_a - 4) ://2 lei shang 1j
			kind2_node(j, j - offset, -(Num9 - Node_9 + Num8 + Num7 + Num6 + Node_5), -2, -1, 0, 1, Num5 - Node_5 + Num4 + Num3 + Num2 + Node_1,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_5 + N_a) :// 2 lei 2j shang
			kind2_node(j, j - offset, -(Num7 / 2 + Num6 + Node_5 + N_a - 3), -2, -1, 0, 1, Num5 - Node_5 - N_a + 3 + Num4 + Num3 / 2,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node52:// inner 3 node
			inner_3node(j, j - offset, -1, 0, (Num5 - Node_55 + 4 + Num4 - 1) - N_c1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node52 + 1) :// inner 3 node
			inner_3node(j, j - offset, 0, 1, (Num5 - Node_55 + 4 + Num4 - 1), k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Node_55 - 3) ://2 lei xia 2j
			kind2_node(j, j - offset, -(Num6 + Node_55 + 1), -1, 0, 1, 2, Num5 - Node_55 + Num4 + Num3 - 1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_55 + 1) ://2 lei xia 1j
			kind2_node(j, j - offset, -Num8 - Num7 - Num6 - Num5, -1, 0, 1, 2, Num4 + Num3 + Num2 + Num1,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < (Node_5 + N_a - 4))// 4dx!
			{
				inner_5node(j, j - offset, 4, -(Num9 - Node_9 + Num8 + Num7 + Num6 + Node_5), Num5 - Node_5 + Num4 + Num3 + Num2 + Node_1, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) > (Node_5 + N_a - 4)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < (Node_5 + N_a)))// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Num7 / 2 + Num6 + Node_5 + N_a - 3), Num5 - Node_5 - N_a + 3 + Num4 + Num3 / 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) > (Node_5 + N_a)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < Node52))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Node_55 - 4 + 1) + 2, (Num5 - Node_55 + 4 + Num4 - 1) - N_c1, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) > (Node52 + 1)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < (Node_55 - 3)))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Node_55 - 4 + 1), (Num5 - Node_55 + 4 + Num4 - 1), Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) > (Node_55 - 3)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < (Node_55 + 1)))// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num6 + Node_55 + 1), Num5 - Node_55 + Num4 + Num3 - 1, Coef, Coef_location);
			}
			else// 4dx!
			{
				inner_5node(j, j - offset, 4, -Num8 - Num7 - Num6 - Num5, Num4 + Num3 + Num2 + Num1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9);
		i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Num4 / 2); i++)// num 14(4) :  30 node! part1
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9)
		{
		case 0:// 1 lei 2j shang
			kind1_node(j, j - offset, -(Num5 - N_d0 - 2) - 2, -(Num5 - N_d0 - 2) - 1, -(Num5 - N_d0 - 2), 0, 1, Num4 + 4 - 2, Num4 + 4 - 1, Num4 + 4,
				A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case (Num4 / 2 - 1) :// 1 lei 2j xia
			kind1_node(j, j - offset, -(Num5 - N_d0 - 2) - N_c1, -(Num5 - N_d0 - 2) - N_c1 + 1, -(Num5 - N_d0 - 2) - N_c1 + 2, -1, 0,
			Num4 + 4, Num4 + 4 + 1, Num4 + 4 + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) <= Node41)// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num5 - N_d0 - 2), Num4 + 4, Coef, Coef_location);
			}
			else if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) >= Node42) // 1dx
			{
				inner_5node(j, j - offset, 1, -(Num5 - N_d0 - 2) - N_c1, Num4 + 4, Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_4node(j, j - offset, -1, 0, 1, Num4 + 4, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Num4 / 2);
		i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9); i++)// num14(4) :  30 node! part2
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9)
		{
		case (Num4 / 2) :// 1 lei 2j shang
			kind1_node(j, j - offset, -(Num5 - Node_55 + 4 + Num4 - 1) + N_c1 - 2, -(Num5 - Node_55 + 4 + Num4 - 1) + N_c1 - 1, -(Num5 - Node_55 + 4 + Num4 - 1) + N_c1, 0, 1,
			(Num3 - 5 + 1) - 2, (Num3 - 5 + 1) - 1, (Num3 - 5 + 1),
			A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case (Num4 - 1) :// 1 lei 2j xia
			kind1_node(j, j - offset, -(Num5 - Node_55 + 4 + Num4 - 1), -(Num5 - Node_55 + 4 + Num4 - 1) + 1, -(Num5 - Node_55 + 4 + Num4 - 1) + 2, -1, 0,
			(Num3 - 5 + 1), (Num3 - 5 + 1) + 1, (Num3 - 5 + 1) + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) <= Node43)// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num5 - Node_55 + 4 + Num4 - 1) + N_c1, (Num3 - 5 + 1), Coef, Coef_location);
			}
			else if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) >= Node44) // 1dx
			{
				inner_5node(j, j - offset, 1, -(Num5 - Node_55 + 4 + Num4 - 1), (Num3 - 5 + 1), Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_4node(j, j - offset, -1, 0, 1, (Num3 - 5 + 1), 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9);
		i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Num3 / 2); i++)// num 15(3) :  46 node! part1
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9)
		{
		case 0:// 1 lei 1j shang
			kind1_node(j, j - offset, -(Num5 - N_d0 + 2 + Num4) - 2, -(Num5 - N_d0 + 2 + Num4) - 1, -(Num5 - N_d0 + 2 + Num4), 0, 1,
				Num3 + Num2 + N_d0 - 2 - 2, Num3 + Num2 + N_d0 - 2 - 1, Num3 + Num2 + N_d0 - 2,
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case 3:// 2 lei 2j shang
			kind2_node(j, j - offset, -(Num5 - N_d0 + 2 + Num4), -2, -1, 0, 1, Num3 + Num2 + N_d0 - 2,
				B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case 4:// 3 lei 2j shang
			kind3_node(j, j - offset, -(Num4 + 4), -1, 0, 1, (Num3 + Num2 + N_d0 - 2) - 1, (Num3 + Num2 + N_d0 - 2),
				C2_2, C3_2, C1_2, C3_2, C4, C4, Coef, Coef_location);
			break;
		case (Num3 / 2 - 5) :// 3 lei 2j xia
			kind3_node(j, j - offset, -(Num4 + 4), -1, 0, 1, Num3 / 2 + 1 + Num2 + Node_1, Num3 / 2 + 1 + Num2 + Node_1 + 1,
			C2_2, C3_2, C1_2, C3_2, C4, C4, Coef, Coef_location);
			break;
		case (Num3 / 2 - 4) :// 2 lei 2j xia
			kind2_node(j, j - offset, -(Num5 - Node_5 + Num4 + Num3 / 2 - 1), -1, 0, 1, 2, Num3 / 2 + 1 + Num2 + Node_1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num3 / 2 - 1) :// 1 lei 1j xia
			kind1_node(j, j - offset, -(Num5 - Node_5 + Num4 + Num3 / 2 - 1), -(Num5 - Node_5 + Num4 + Num3 / 2 - 1) + 1, -(Num5 - Node_5 + Num4 + Num3 / 2 - 1) + 2, -1, 0,
			Num3 / 2 + 1 + Num2 + Node_1, Num3 / 2 + 1 + Num2 + Node_1 + 1, Num3 / 2 + 1 + Num2 + Node_1 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < 3)// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num5 - N_d0 + 2 + Num4), Num3 + Num2 + N_d0 - 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) > 4) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < (Num3 / 2 - 5)))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num4 + 4), Num3 - 5, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num5 - Node_5 + Num4 + Num3 / 2 - 1), Num3 / 2 + 1 + Num2 + Node_1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Num3 / 2);
		i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 * 2 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9); i++)// num 15(3) :  46 node! part2
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num3 / 2)
		{
		case 0:// 1 lei 1j shang
			kind1_node(j, j - offset, -(Num5 - Node_5 - N_a + 3 + Num4 + Num3 / 2) - 2, -(Num5 - Node_5 - N_a + 3 + Num4 + Num3 / 2) - 1, -(Num5 - Node_5 - N_a + 3 + Num4 + Num3 / 2), 0, 1,
				Num3 / 2 + Num2 + Node_1 + N_a - 3 - 2, Num3 / 2 + Num2 + Node_1 + N_a - 3 - 1, Num3 / 2 + Num2 + Node_1 + N_a - 3,
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case 3:// 2 lei 2j shang
			kind2_node(j, j - offset, -(Num5 - Node_5 - N_a + 3 + Num4 + Num3 / 2), -2, -1, 0, 1, Num3 / 2 + Num2 + Node_1 + N_a - 3,
				B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case 4:// 3 lei 2j shang
			kind3_node(j, j - offset, -(Num3 - 5 + 1), -1, 0, 1, (Num3 / 2 + Num2 + Node_1 + N_a - 3) - 1, (Num3 / 2 + Num2 + Node_1 + N_a - 3),
				C2_2, C3_2, C1_2, C3_2, C4, C4, Coef, Coef_location);
			break;
		case (Num3 / 2 - 5) :// 3 lei 2j xia
			kind3_node(j, j - offset, -(Num3 - 5 + 1), -1, 0, 1, 1 + Num2 + Node_11, 1 + Num2 + Node_11 + 1,
			C2_2, C3_2, C1_2, C3_2, C4, C4, Coef, Coef_location);
			break;
		case (Num3 / 2 - 4) :// 2 lei 2j xia
			kind2_node(j, j - offset, -(Num5 - Node_55 + Num4 + Num3 - 1), -1, 0, 1, 2, 1 + Num2 + Node_11,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num3 / 2 - 1) :// 1 lei 1j xia
			kind1_node(j, j - offset, -(Num5 - Node_55 + Num4 + Num3 - 1), -(Num5 - Node_55 + Num4 + Num3 - 1) + 1, -(Num5 - Node_55 + Num4 + Num3 - 1) + 2, -1, 0,
			1 + Num2 + Node_11, 1 + Num2 + Node_11 + 1, 1 + Num2 + Node_11 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num3 / 2) < 3)// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num5 - Node_5 - N_a + 3 + Num4 + Num3 / 2), Num3 / 2 + Num2 + Node_1 + N_a - 3, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num3 / 2) > 4) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num3 / 2) < (Num3 / 2 - 5)))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num3 - 5 + 1), 6 + Num2 - 1, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num5 - Node_55 + Num4 + Num3 - 1), 1 + Num2 + Node_11, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 * 2 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9);
		i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 * 2 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Num2 / 2); i++)// num 16(2) :  26 node! part1
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9)
		{
		case 0:// 3 lei 2j shang
			kind3_node(j, j - offset, -(Num3 - 3), -(Num3 - 3) + 2, 0, 1, Num2 + Node11 - 1, Num2 + Node11,
				C4, C3_2, C1_2, C2_2, C4, C3_2, Coef, Coef_location);
			break;
		case (Num2 / 2 - 1) ://node2: 3 lei 2j xia
			kind3_node(j, j - offset, -(Num3 - 5), -(Num3 - 5) + 2, -1, 0, Num2 / 2 + Node_1 - 4 + 1, Num2 / 2 + Node_1 - 4 + 1 + 1,
			C3_2, C4, C2_2, C1_2, C3_2, C4, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) % 2 == 0)// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num3 - 5),
					(Num2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) + Node11 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) / 2), Coef, Coef_location);
			}
			else//1 lei  2j left
			{
				kind1_node(j, j - offset, -(Num3 - 5), -1, 0, 1,
					(Num2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) + Node11 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) / 2),
					(Num2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) + Node11 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) / 2) + 1,
					(Num2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) + Node11 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) / 2) + (Num1 - N_d0 + 2),
					(Num2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) + Node11 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) / 2) + 1 + (Num1 - N_d0 + 2),
					A4_2, A3_2, A5_2, A3_2, A1_2, A1_2, A2_2, A2_2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 * 2 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Num2 / 2);
		i<(Nx_2 + Nx_3 + Num1 + Num2 * 2 + Num3 * 2 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9); i++)// num 16(2) :  26 node! part2
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2)
		{
		case 0:// 3 lei 2j shang
			kind3_node(j, j - offset, -(6 + Num2 - 1) - 2, -(6 + Num2 - 1), 0, 1, Num2 / 2 + Node13 - 1, Num2 / 2 + Node13,
				C4, C3_2, C1_2, C2_2, C4, C3_2, Coef, Coef_location);
			break;
		case (Num2 / 2 - 1) :// 3 lei 2j xia
			kind3_node(j, j - offset, -(6 + Num2 - 1), -(6 + Num2 - 1) + 2, -1, 0, Node_11 - 4 + 1, Node_11 - 4 + 1 + 1,
			C3_2, C4, C2_2, C1_2, C3_2, C4, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) % 2 == 0)// 1dx
			{
				inner_5node(j, j - offset, 1, -(6 + Num2 - 1), Num2 / 2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) + Node13 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) / 2,
					Coef, Coef_location);
			}
			else//1 lei 2j left
			{
				kind1_node(j, j - offset, -(6 + Num2 - 1), -1, 0, 1,
					(Num2 / 2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) + Node13 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) / 2),
					(Num2 / 2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) + Node13 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) / 2) + 1,
					(Num2 / 2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) + Node13 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) / 2) + (Num1 - Node_11 + Nx_3 - 1),
					(Num2 / 2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) + Node13 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) / 2) + 1 + (Num1 - Node_11 + Nx_3 - 1),
					A4_2, A3_2, A5_2, A3_2, A1_2, A1_2, A2_2, A2_2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (num_total - Nx_2 - Nx_3 - Num1); i<(num_total - Nx_2 - Nx_3 - Num1 + Node_1 + 2); i++)// circle num17(1) Nx_2:  74 node!  part1 !
	{
		switch (i - num_total + Nx_2 + Nx_3 + Num1)
		{
		case 0:
			boundary_node(j, j - offset, 1, -Num5 - Num4 - Num3 - Num2, Nx_2 + Nx_3, Coef, Coef_location);
			break;
		case (N_d0 - 3) ://2 lei shang 1j
			kind2_node(j, j - offset, -Num5 - Num4 - Num3 - Num2, -2, -1, 0, 1, Nx_2 + Nx_3,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_1 + 1) ://2 lei xia 1j
			kind2_node(j, j - offset, -(Num5 - Node_5 + Num4 + Num3 + Num2 + Node_1), -1, 0, 1, 2, Nx_2 + Nx_3,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - num_total + Nx_2 + Nx_3 + Num1) < (N_d0 - 3))// 4dx!
			{
				inner_5node(j, j - offset, 4, -Num5 - Num4 - Num3 - Num2, Nx_2 + Nx_3, Coef, Coef_location);
			}
			else if (((i - num_total + Nx_2 + Nx_3 + Num1) > (N_d0 - 3)) && ((i - num_total + Nx_2 + Nx_3 + Num1) < Node11))// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Num3 + Num2 + N_d0 - 2), Num1 - N_d0 + 2, Coef, Coef_location);
			}
			else if (((i - num_total + Nx_2 + Nx_3 + Num1) >= Node11) && ((i - num_total + Nx_2 + Nx_3 + Num1) <= Node12))// 2 lei right 2j !
			{
				kind2_node(j, j - offset, -((i - num_total + Nx_2 + Nx_3 + Num1) + Num2 - (i - num_total + Nx_2 + Nx_3 + Num1 - Node11) * 2),
					-1, 0, 1, Num1 - N_d0 + 2, Nx_2 + Nx_3, B2, B3, B1_2, B3, B4, B5, Coef, Coef_location);
			}
			else// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Num3 / 2 + 1 + Num2 + Node_1), Num1 - N_d0 + 2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (num_total - Nx_2 - Nx_3 - Num1 + Node_1 + 2); i<(num_total - Nx_2 - Nx_3); i++)// circle num17(1) Nx_2:  74 node!  part2 !
	{
		switch (i - num_total + Nx_2 + Nx_3 + Num1)
		{
		case (Num1 - 1) :
			boundary_node(j, j - offset, 2, -Num1 - Num2 - Num3 - Num4, Nx_2 + Nx_3, Coef, Coef_location);
			break;
		case (Node_1 + N_a - 4) ://2 lei shang 1j
			kind2_node(j, j - offset, -(Num5 - Node_5 + Num4 + Num3 + Num2 + Node_1), -2, -1, 0, 1, Nx_2 + Nx_3,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_11 + 1) ://2 lei xia 1j
			kind2_node(j, j - offset, -Num1 - Num2 - Num3 - Num4, -1, 0, 1, 2, Nx_2 + Nx_3,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - num_total + Nx_2 + Nx_3 + Num1) < (Node_1 + N_a - 4))// 4dx!
			{
				inner_5node(j, j - offset, 4, -(Num5 - Node_5 + Num4 + Num3 + Num2 + Node_1), Nx_2 + Nx_3, Coef, Coef_location);
			}
			else if (((i - num_total + Nx_2 + Nx_3 + Num1) > (Node_1 + N_a - 4)) && ((i - num_total + Nx_2 + Nx_3 + Num1) < Node13))// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Num3 / 2 + Num2 + Node_1 + N_a - 3), Num1 - Node_11 + Nx_3 - 1, Coef, Coef_location);
			}
			else if (((i - num_total + Nx_2 + Nx_3 + Num1) >= Node13) && ((i - num_total + Nx_2 + Nx_3 + Num1) <= Node14))// 2 lei left 2j !
			{
				kind2_node(j, j - offset, -((i - num_total + Nx_2 + Nx_3 + Num1) + Num2 / 2 - (i - num_total + Nx_2 + Nx_3 + Num1 - Node13) * 2),
					-1, 0, 1, Num1 - Node_11 + Nx_3 - 1, Nx_2 + Nx_3, B2, B3, B1_2, B3, B4, B5, Coef, Coef_location);
			}
			else if (((i - num_total + Nx_2 + Nx_3 + Num1) >= Node14) && ((i - num_total + Nx_2 + Nx_3 + Num1) <= Node_11))// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num2 + Node_11 + 1), Num1 - Node_11 + Nx_3 - 1, Coef, Coef_location);
			}
			else
			{
				inner_5node(j, j - offset, 4, -Num1 - Num2 - Num3 - Num4, Nx_2 + Nx_3, Coef, Coef_location);// 4dx!
			}
			break;
		}
		j = j + 1;
	}
	//*****************************************************************************************
	for (i = (Nx_2); i<(Nx_2 + Nx_3); i++)//left 3 Nx_3:  30 node! 
	{
		switch (i - Nx_2)
		{
		case 0:// 1 lei shang
			kind1_node(j, j - offset, -(Nx_2 - N_d0 + 2) - 2, -(Nx_2 - N_d0 + 2) - 1, -(Nx_2 - N_d0 + 2), 0, 1, Nx_3 + N_d0 - 2 - 2, Nx_3 + N_d0 - 2 - 1, Nx_3 + N_d0 - 2,
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case (Nx_3 / 2 - 1) :// 1 lei xia
			kind1_node(j, j - offset, -(Nx_2 - N_d0 + 2), -(Nx_2 - N_d0 + 2) + 1, -(Nx_2 - N_d0 + 2) + 2, -1, 0,
			Nx_3 + N_d0 - 2, Nx_3 + N_d0 - 2 + 1, Nx_3 + N_d0 - 2 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		case (Nx_3 / 2) :// 1 lei shang
			kind1_node(j, j - offset, -(Nx_2 - Node_11 - 1 + Nx_3) - 2, -(Nx_2 - Node_11 - 1 + Nx_3) - 1, -(Nx_2 - Node_11 - 1 + Nx_3), 0, 1,
			1 + Node_11 - 2, 1 + Node_11 - 1, 1 + Node_11,
			A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case (Nx_3 - 1) :// 1 lei xia
			kind1_node(j, j - offset, -(Nx_2 - Node_11 - 1 + Nx_3), -(Nx_2 - Node_11 - 1 + Nx_3) + 1, -(Nx_2 - Node_11 - 1 + Nx_3) + 2, -1, 0,
			1 + Node_11, 1 + Node_11 + 1, 1 + Node_11 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2) < (Nx_3 / 2 - 1))
			{
				inner_5node(j, j - offset, 2, -(Nx_2 - N_d0 + 2), Nx_3 + N_d0 - 2, Coef, Coef_location);
			}
			else
			{
				inner_5node(j, j - offset, 2, -(Nx_2 - Node_11 - 1 + Nx_3), 1 + Node_11, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//***************************************************************************************
	for (i = (Nx_2 + Nx_3); i<(Nx_2 * 2 + Nx_3); i++)//left 2 Nx_2:  77 node! 
	{
		switch (i - Nx_2 - Nx_3)
		{
		case 0:
			boundary_node(j, j - offset, 1, -Nx_2 - Nx_3, Nx_2 + Nx_1, Coef, Coef_location);
			break;
		case (Nx_2 - 1) :
			boundary_node(j, j - offset, 2, -Nx_3 - Nx_2, Nx_1 + Nx_side, Coef, Coef_location);
			break;
		case (N_d0 - 3) ://2 lei shang
			kind2_node(j, j - offset, -Nx_2 - Nx_3, -2, -1, 0, 1, Nx_2 + Nx_1,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (N_d0 - 2) ://3 lei shang
			kind3_node(j, j - offset, -(Nx_3 + N_d0 - 2), -1, 0, 1, Nx_2 + Nx_1 - 1, Nx_2 + Nx_1,
			C2_1, C3_1, C1_1, C3_1, C4, C4, Coef, Coef_location);
			break;
		case Node_1://3 lei xia
			kind3_node(j, j - offset, -(Nx_3 + N_d0 - 2), -1, 0, 1, Nx_2 - Node_1 + Nx_1 + N_d0 + 1, Nx_2 - Node_1 + Nx_1 + N_d0 + 1 + 1,
				C2_1, C3_1, C1_1, C3_1, C4, C4, Coef, Coef_location);
			break;
		case (Node_1 + 1) ://2 lei xia
			kind2_node(j, j - offset, -(Nx_2 + Nx_3), -1, 0, 1, 2, Nx_2 - Node_1 + Nx_1 + N_d0 + 1,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_1 + N_a - 4) ://2 lei shang
			kind2_node(j, j - offset, -(Nx_2 + Nx_3), -2, -1, 0, 1, Nx_2 - Node_1 + Nx_1 + N_d0 + 1,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_1 + N_a - 3) ://3 lei shang
			kind3_node(j, j - offset, -(Node_11 + 1), -1, 0, 1, Nx_2 - Node_1 + Nx_1 + N_d0 + 1 - 1, Nx_2 - Node_1 + Nx_1 + N_d0 + 1,
			C2_1, C3_1, C1_1, C3_1, C4, C4, Coef, Coef_location);
			break;
		case Node_11://3 lei xia
			kind3_node(j, j - offset, -(Node_11 + 1), -1, 0, 1, Nx_side + Nx_1, Nx_side + Nx_1 + 1,
				C2_1, C3_1, C1_1, C3_1, C4, C4, Coef, Coef_location);
			break;
		case (Node_11 + 1) ://2 lei xia
			kind2_node(j, j - offset, -Nx_3 - Nx_2, -1, 0, 1, 2, Nx_1 + Nx_side,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3) < (N_d0 - 3))// 2dx!
			{
				inner_5node(j, j - offset, 4, -Nx_3 - Nx_2, Nx_2 + Nx_1, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3) >= (N_d0 - 1)) && ((i - Nx_2 - Nx_3) < Node_1))
			{
				inner_5node(j, j - offset, 2, -(Nx_3 + N_d0 - 2), Nx_2 - N_d0 + 1, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3) > (Node_1 + 1)) && ((i - Nx_2 - Nx_3) < (Node_1 + N_a - 4)))// 2dx!
			{
				inner_5node(j, j - offset, 4, -(Nx_2 + Nx_3), Nx_2 - Node_1 + Nx_1 + N_d0 + 1, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3) >= (Node_1 + N_a - 2)) && ((i - Nx_2 - Nx_3) < Node_11))
			{
				inner_5node(j, j - offset, 2, -(Node_11 + 1), Nx_1 + Nx_2 - Node_11, Coef, Coef_location);
			}
			else// 2dx!
			{
				inner_5node(j, j - offset, 4, -Nx_3 - Nx_2, Nx_side + Nx_1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (Nx_2 * 2 + Nx_3); i<(Nx_2 * 2 + Nx_3 + Nx_1 / 2); i++)// left 1 Nx_1: 26 node! part1
	{
		switch (i - Nx_2 * 2 - Nx_3)
		{
		case 0:// 3 lei
			kind3_node(j, j - offset, -(Nx_2 - N_d0 + 1) - 2, -(Nx_2 - N_d0 + 1), 0, 1, Nx_1 + N_d0 - 3, Nx_1 + N_d0 - 2,
				C4, C3_1, C1_1, C2_1, C4, C3_1, Coef, Coef_location);
			break;
		case (Nx_1 / 2 - 1) ://node2: 3 lei
			kind3_node(j, j - offset, -(Nx_2 - N_d0 + 1), -(Nx_2 - N_d0 + 1) + 2, -1, 0, Nx_1 / 2 + 1 + N_d0 + 1, Nx_1 / 2 + 1 + N_d0 + 2,
			C3_1, C4, C2_1, C1_1, C3_1, C4, Coef, Coef_location);
			break;
		case 1:// 1 lei
			kind1_node(j, j - offset, -(Nx_2 - N_d0 + 1), -1, 0, 1, (Nx_1 - 1 + N_d0 - 2), (Nx_1 - 1 + N_d0 - 2) + 1, (Nx_1 - 1 + N_d0 - 2) + Nx_side, (Nx_1 - 1 + N_d0 - 2) + 1 + Nx_side,
				A4_1, A3_1, A5_1, A3_1, A1_1, A1_1, A2_1, A2_1, Coef, Coef_location);
			break;
		case 2:// 5 node 2dx
			inner_5node(j, j - offset, 2, -(Nx_2 - N_d0 + 1), Nx_1 - 2 + N_d0 - 1, Coef, Coef_location);
			break;
		case 3:// 1 lei que 2 dian
			kind2_node(j, j - offset, -(Nx_2 - N_d0 + 1), -1, 0, 1, Nx_1 - 3 + N_d0 - 1, Nx_1 - 3 + N_d0 - 1 + Nx_side,
				A4_1, A3_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		case (Nx_1 / 2 - 4) :// 1 lei que 2 dian
			kind2_node(j, j - offset, -(Nx_2 - N_d0 + 1), -1, 0, 1, N_d0 + Nx_1 / 2 + 4, N_d0 + Nx_1 / 2 + 4 + Nx_side,
			A4_1, A3_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		case (Nx_1 / 2 - 3) :// 5 node 2dx
			inner_5node(j, j - offset, 2, -(Nx_2 - N_d0 + 1), Nx_1 / 2 + 3 + N_d0, Coef, Coef_location);
			break;
		case (Nx_1 / 2 - 2) :// 1 lei
			kind1_node(j, j - offset, -(Nx_2 - N_d0 + 1), -1, 0, 1, (Nx_1 / 2 + 2 + N_d0), (Nx_1 / 2 + 2 + N_d0) + 1, (Nx_1 / 2 + 2 + N_d0) + Nx_side, (Nx_1 / 2 + 2 + N_d0) + 1 + Nx_side,
			A4_1, A3_1, A5_1, A3_1, A1_1, A1_1, A2_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 * 2 - Nx_3) % 2 == 0)
			{
				inner_4node(j, j - offset, -(Nx_2 - N_d0 + 1), -1, 0, 1, 1, 1, k*k * 4 * dx*dx - 4, 1, Coef, Coef_location);
			}
			else//1 lei
			{
				inner_4node(j, j - offset, -(Nx_2 - N_d0 + 1), -1, 0, 1, A4_1, A3_1, A5_1, A3_1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 * 2 + Nx_3 + Nx_1 / 2); i<(Nx_2 * 2 + Nx_3 + Nx_1); i++)//left 1 Nx_1: 26 node! part2
	{
		switch (i - Nx_2 * 2 - Nx_3 - Nx_1 / 2)
		{
		case 0:// 3 lei
			kind3_node(j, j - offset, -(Nx_2 - Node_11 + Nx_1) - 2, -(Nx_2 - Node_11 + Nx_1), 0, 1, Nx_1 / 2 + N_d0 + N_a - 3, Nx_1 / 2 + N_d0 + N_a - 2,
				C4, C3_1, C1_1, C2_1, C4, C3_1, Coef, Coef_location);
			break;
		case (Nx_1 / 2 - 1) ://node2: 3 lei
			kind3_node(j, j - offset, -(Nx_2 - Node_11 + Nx_1), -(Nx_2 - Node_11 + Nx_1) + 2, -1, 0, N_d0 + N_a + 1 + 1, N_d0 + N_a + 2 + 1,
			C3_1, C4, C2_1, C1_1, C3_1, C4, Coef, Coef_location);
			break;
		case 1:// 1 lei
			kind1_node(j, j - offset, -(Nx_2 - Node_11 + Nx_1), -1, 0, 1, (Nx_1 / 2 - 1 + N_d0 + N_a - 2), (Nx_1 / 2 - 1 + N_d0 + N_a - 2) + 1,
				(Nx_1 / 2 - 1 + N_d0 + N_a - 2) + Nx_side, (Nx_1 / 2 - 1 + N_d0 + N_a - 2) + 1 + Nx_side,
				A4_1, A3_1, A5_1, A3_1, A1_1, A1_1, A2_1, A2_1, Coef, Coef_location);
			break;
		case 2:// 5 node 2dx
			inner_5node(j, j - offset, 2, -(Nx_2 - Node_11 + Nx_1), Nx_1 / 2 - 2 + N_d0 + N_a - 1, Coef, Coef_location);
			break;
		case 3:// 1 lei que 2 dian
			kind2_node(j, j - offset, -(Nx_2 - Node_11 + Nx_1), -1, 0, 1, Nx_1 / 2 - 3 + N_d0 + N_a - 1, Nx_1 / 2 - 3 + N_d0 + N_a - 1 + Nx_side,
				A4_1, A3_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		case (Nx_1 / 2 - 4) :// 1 lei que 2 dian
			kind2_node(j, j - offset, -(Nx_2 - Node_11 + Nx_1), -1, 0, 1, N_d0 + N_a + 4, N_d0 + N_a + 4 + Nx_side,
			A4_1, A3_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		case (Nx_1 / 2 - 3) :// 5 node 2dx
			inner_5node(j, j - offset, 2, -(Nx_2 - Node_11 + Nx_1), 3 + N_d0 + N_a, Coef, Coef_location);
			break;
		case (Nx_1 / 2 - 2) :// 1 lei
			kind1_node(j, j - offset, -(Nx_2 - Node_11 + Nx_1), -1, 0, 1, (2 + N_d0 + N_a), (2 + N_d0 + N_a) + 1, (2 + N_d0 + N_a) + Nx_side, (2 + N_d0 + N_a) + 1 + Nx_side,
			A4_1, A3_1, A5_1, A3_1, A1_1, A1_1, A2_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 * 2 - Nx_3 - Nx_1 / 2) % 2 == 0)// 2dx
			{
				inner_4node(j, j - offset, -(Nx_2 - Node_11 + Nx_1), -1, 0, 1, 1, 1, k*k * 4 * dx*dx - 4, 1, Coef, Coef_location);
			}
			else//1 lei
			{
				inner_4node(j, j - offset, -(Nx_2 - Node_11 + Nx_1), -1, 0, 1, A4_1, A3_1, A5_1, A3_1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!!!!
	for (i = (Nx_2 * 2 + Nx_3 + Nx_1); i<(num_total - Nx_side*N_d1); i++)//middle left ,one line
	{
		switch (i - Nx_2 * 2 - Nx_3 - Nx_1)
		{
		case 0:
			boundary_node(j, j - offset, 1, -Nx_2 - Nx_1, Nx_side, Coef, Coef_location);
			break;
		case (Nx_side - 1) :
			boundary_node(j, j - offset, 2, -Nx_1 - Nx_side, Nx_side, Coef, Coef_location);
			break;
		case (N_d0 - 2) :// 2 lei
			kind2_node(j, j - offset, -(Nx_1 + N_d0 - 2), -1, 0, 1, Nx_side, Nx_side * 2, B2, B3, B1_1, B3, B4, B5, Coef, Coef_location);
			break;
		case (N_d0 - 1) :// 2 lei缺一点
			other_5_node(j, j - offset, -(Nx_1 - 2 + N_d0 - 1), -1, 0, Nx_side, 2 * Nx_side, B2, B3, B1_1, B4, B5, Coef, Coef_location);
			break;
		case N_d0:// 2 lei缺一点
			other_5_node(j, j - offset, -(Nx_1 / 2 + 3 + N_d0), 0, 1, Nx_side, 2 * Nx_side, B2, B1_1, B3, B4, B5, Coef, Coef_location);
			break;
		case (N_d0 + 1) :// 2 lei
			kind2_node(j, j - offset, -(Nx_1 / 2 + 1 + N_d0 + 1), -1, 0, 1, Nx_side, Nx_side * 2, B2, B3, B1_1, B3, B4, B5, Coef, Coef_location);
			break;
		case (N_d0 + N_a - 2) :// 2 lei
			kind2_node(j, j - offset, -(Nx_1 / 2 + N_a + N_d0 - 2), -1, 0, 1, Nx_side, Nx_side * 2, B2, B3, B1_1, B3, B4, B5, Coef, Coef_location);
			break;
		case (N_d0 + N_a - 1) :// 2 lei缺一点
			other_5_node(j, j - offset, -(Nx_1 / 2 - 2 + N_d0 + N_a - 1), -1, 0, Nx_side, 2 * Nx_side, B2, B3, B1_1, B4, B5, Coef, Coef_location);
			break;
		case (N_d0 + N_a) :// 2 lei缺一点
			other_5_node(j, j - offset, -(N_d0 + N_a + 3), 0, 1, Nx_side, 2 * Nx_side, B2, B1_1, B3, B4, B5, Coef, Coef_location);
			break;
		case (N_d0 + N_a + 1) :// 2 lei
			kind2_node(j, j - offset, -(1 + N_a + N_d0 + 1), -1, 0, 1, Nx_side, Nx_side * 2, B2, B3, B1_1, B3, B4, B5, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 * 2 - Nx_3 - Nx_1) < (N_d0 - 2))//up
			{
				inner_5node(j, j - offset, 4, -Nx_2 - Nx_1, Nx_side, Coef, Coef_location);
			}
			else if ((i - Nx_2 * 2 - Nx_3 - Nx_1) > (N_d0 + N_a + 1))//down
			{
				inner_5node(j, j - offset, 4, -Nx_1 - Nx_side, Nx_side, Coef, Coef_location);
			}
			else//middle
			{
				inner_5node(j, j - offset, 4, -(Nx_2 - Node_1 + Nx_1 + N_d0 + 1), Nx_side, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (num_total - Nx_side*N_d1); i<(num_total - Nx_side); i++)//middle
	{
		switch ((i - num_total + Nx_side*N_d1) % Nx_side)
		{
		case 0:
			boundary_node(j, j - offset, 1, -Nx_side, Nx_side, Coef, Coef_location);
			break;
		case (Nx_side - 1) :
			boundary_node(j, j - offset, 2, -Nx_side, Nx_side, Coef, Coef_location);
			break;
		case (N_d0 - 1) :
			inner_4node(j, j - offset, -Nx_side, -1, 0, Nx_side, 1, 1, k*k * 16 * dx*dx - 4, 1, Coef, Coef_location);
			break;
		case N_d0:
			inner_4node(j, j - offset, -Nx_side, 0, 1, Nx_side, 1, k*k * 16 * dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (N_d0 + N_a - 1) :
			inner_4node(j, j - offset, -Nx_side, -1, 0, Nx_side, 1, 1, k*k * 16 * dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (N_d0 + N_a) :
			inner_4node(j, j - offset, -Nx_side, 0, 1, Nx_side, 1, k*k * 16 * dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		default:
			inner_5node(j, j - offset, 4, -Nx_side, Nx_side, Coef, Coef_location);
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (num_total - Nx_side); i<num_total; i++)//截断右边界 1阶mur吸收边界条件
	{
		switch (i - num_total + Nx_side)
		{
		case 0:          //角点处理：场平均吸收边条
			Coef[j][2].real = -3;
			Coef[j][0].image = Coef[j][3].image = -sin(k * 4 * dx);
			Coef[j][0].real = Coef[j][3].real = cos(k * 4 * dx);
			Coef[j][1].image = -sin(sqrt(2.0)*k * 4 * dx);
			Coef[j][1].real = cos(sqrt(2.0)*k * 4 * dx);
			Coef_location[j][0] = j - offset - Nx_side;
			Coef_location[j][1] = j - offset - Nx_side + 1;
			Coef_location[j][2] = j - offset;
			Coef_location[j][3] = j - offset + 1;
			break;
		case (Nx_side - 1) :
			Coef[j][3].real = -3;
			Coef[j][1].image = Coef[j][2].image = -sin(k * 4 * dx);
			Coef[j][1].real = Coef[j][2].real = cos(k * 4 * dx);
			Coef[j][0].image = -sin(sqrt(2.0)*k * 4 * dx);
			Coef[j][0].real = cos(sqrt(2.0)*k * 4 * dx);
			Coef_location[j][0] = j - offset - Nx_side - 1;
			Coef_location[j][1] = j - offset - Nx_side;
			Coef_location[j][2] = j - offset - 1;
			Coef_location[j][3] = j - offset;
			break;
		case (N_d0 - 1) :
			Coef[j][2].real = (2 - 2 * k*k * 16 * dx*dx)*cos(k * 4 * dx) - 2;
			Coef[j][2].image = (2 - 2 * k*k * 16 * dx*dx)*sin(k * 4 * dx);
			Coef[j][1].real = 1 - cos(k * 4 * dx);
			Coef[j][1].image = -sin(k * 4 * dx);
			Coef[j][0].real = 2 * k*k * 16 * dx*dx;
			Coef_location[j][2] = j - offset;
			Coef_location[j][1] = j - offset - 1;
			Coef_location[j][0] = j - offset - Nx_side;
			break;
		case (N_d0 + N_a) :
			Coef[j][2].real = (2 - 2 * k*k * 16 * dx*dx)*cos(k * 4 * dx) - 2;
			Coef[j][2].image = (2 - 2 * k*k * 16 * dx*dx)*sin(k * 4 * dx);
			Coef[j][3].real = 1 - cos(k * 4 * dx);
			Coef[j][3].image = -sin(k * 4 * dx);
			Coef[j][0].real = 2 * k*k * 16 * dx*dx;
			Coef_location[j][2] = j - offset;
			Coef_location[j][3] = j - offset + 1;
			Coef_location[j][0] = j - offset - Nx_side;
			break;
		default:
			if (((i - num_total + Nx_side) >= N_d0) && ((i - num_total + Nx_side) < (N_d0 + N_a)))
			{
				Coef[j][2].real = -1;
				Coef[j][0].image = -sin(kz * 4 * dx);
				Coef[j][0].real = cos(kz * 4 * dx);
				Coef_location[j][2] = j - offset;
				Coef_location[j][0] = j - offset - Nx_side;
			}
			else
			{
				Coef[j][2].real = (2 - 2 * k*k * 16 * dx*dx)*cos(k * 4 * dx) - 2;
				Coef[j][2].image = (2 - 2 * k*k * 16 * dx*dx)*sin(k * 4 * dx);
				Coef[j][1].real = Coef[j][3].real = 1 - cos(k * 4 * dx);
				Coef[j][1].image = Coef[j][3].image = -sin(k * 4 * dx);
				Coef[j][0].real = 2 * k*k * 16 * dx*dx;
				Coef_location[j][2] = j - offset;
				Coef_location[j][1] = j - offset - 1;
				Coef_location[j][3] = j - offset + 1;
				Coef_location[j][0] = j - offset - Nx_side;
			}
			break;
		}
		j = j + 1;
	}
	cout << "end2区域节点数：" << j << endl;
	if (j == sizeM)
		cout << "end2 passed..." << endl;
	else
		cout << "end2 failed..." << endl;
	//修正矩阵，将不在本区域的节点刨除
	for (int i = 0; i != sizeM; i++)
		for (int j = 0; j != N_matrix; j++){
		if (Coef_location[i][j] < 0 || Coef_location[i][j] >= sizeN)
			Coef_location[i][j] = -1;
		}
	if (Debug)
	{
		ofstream outFile_Coef("end2.txt");
		if (!outFile_Coef)
		{
			cout << "Can't write the coefficient matrix into file !" << endl;
			exit(1);//Terminate with fault
		}
		else
		{
			int i, j;
			for (i = 0; i < sizeM; i++)
			{
				for (j = 0; j < N_matrix; j++)
				{
					if ((Coef[i][j].real != 0 || Coef[i][j].image != 0) && Coef_location[i][j] >= 0)
						outFile_Coef << Coef[i][j].real << "+j*" << Coef[i][j].image << " (" << i << "," << Coef_location[i][j] << ") ";
				}
				outFile_Coef << endl << endl;
			}
		}
		outFile_Coef.close();
	}

}
//*********************************************************************************//
//中间圆孔区域扫参（共20个圆孔）：n代表重叠区域号，Ey[n][M]代表每次迭代所求的解
//*********************************************************************************//
//具有相同区域填充的封装函数(方便矩阵填充 本函数填充的是整体矩阵）
void scan_coef_Matrix_middle_wrapper(int &j, complex** Coef, int** Coef_location,int num_total, int num, int offset){
	int i;
	//*****************************************************************************************
	for (i = (Nx_2); i<(Nx_2 + Nx_3); i++)//left 2 Nx_3:  30 node! 
	{
		switch (i - Nx_2)
		{
		case 0:// 1 lei shang
			kind1_node(j, i, -(Nx_2 - N_d0 + 2) - 2, -(Nx_2 - N_d0 + 2) - 1, -(Nx_2 - N_d0 + 2), 0, 1, Nx_3 + N_d0 - 2 - 2, Nx_3 + N_d0 - 2 - 1, Nx_3 + N_d0 - 2,
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case (Nx_3 / 2 - 1) :// 1 lei xia
			kind1_node(j, i, -(Nx_2 - N_d0 + 2), -(Nx_2 - N_d0 + 2) + 1, -(Nx_2 - N_d0 + 2) + 2, -1, 0,
			Nx_3 + N_d0 - 2, Nx_3 + N_d0 - 2 + 1, Nx_3 + N_d0 - 2 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		case (Nx_3 / 2) :// 1 lei shang
			kind1_node(j, i, -(Nx_2 - Node_11 - 1 + Nx_3) - 2, -(Nx_2 - Node_11 - 1 + Nx_3) - 1, -(Nx_2 - Node_11 - 1 + Nx_3), 0, 1,
			1 + Node_11 - 2, 1 + Node_11 - 1, 1 + Node_11,
			A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case (Nx_3 - 1) :// 1 lei xia
			kind1_node(j, i, -(Nx_2 - Node_11 - 1 + Nx_3), -(Nx_2 - Node_11 - 1 + Nx_3) + 1, -(Nx_2 - Node_11 - 1 + Nx_3) + 2, -1, 0,
			1 + Node_11, 1 + Node_11 + 1, 1 + Node_11 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2) < (Nx_3 / 2 - 1))
			{
				inner_5node(j, i, 2, -(Nx_2 - N_d0 + 2), Nx_3 + N_d0 - 2, Coef, Coef_location);
			}
			else
			{
				inner_5node(j, i, 2, -(Nx_2 - Node_11 - 1 + Nx_3), 1 + Node_11, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (Nx_2 + Nx_3); i<(Nx_2 + Nx_3 + Node_1 + 2); i++)// circle num1 Nx_2:  74 node!  part1 !
	{
		switch (i - Nx_2 - Nx_3)
		{
		case 0:
			boundary_node(j, i, 1, -Nx_2 - Nx_3, Num1 + Num2 + Num3 + Num4, Coef, Coef_location);
			break;
		case (N_d0 - 3) ://2 lei shang 1j
			kind2_node(j, i, -Nx_2 - Nx_3, -2, -1, 0, 1, Num1 + Num2 + Num3 + Num4,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_1 + 1) ://2 lei xia 1j
			kind2_node(j, i, -(Nx_2 + Nx_3), -1, 0, 1, 2, Num1 - Node_1 + Num2 + Num3 + Num4 + Node_5,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3) < (N_d0 - 3))// 4dx!
			{
				inner_5node(j, i, 4, -Nx_3 - Nx_2, Num1 + Num2 + Num3 + Num4, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3) > (N_d0 - 3)) && ((i - Nx_2 - Nx_3) < Node11))// 2dx!
			{
				inner_5node(j, i, 2, -(Nx_3 + N_d0 - 2), Num1 - N_d0 + 2 + Num2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3) >= Node11) && ((i - Nx_2 - Nx_3) <= Node12))// 2 lei left 2j !
			{
				kind2_node(j, i, -Nx_3 - Nx_2, -(Nx_3 + N_d0 - 2), -1, 0, 1, Num1 - (i - Nx_2 - Nx_3) + ((i - Nx_2 - Nx_3) - Node11) * 2,
					B5, B4, B3, B1_2, B3, B2, Coef, Coef_location);
			}
			else// 2dx!
			{
				inner_5node(j, i, 2, -(Nx_3 + N_d0 - 2), Num1 - Node_1 + Num2 + Num3 / 2 - 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Node_1 + 2); i<(Nx_2 + Nx_3 + Num1); i++)// circle num1 Nx_2:  74 node!  part2 !
	{
		switch (i - Nx_2 - Nx_3)
		{
		case (Num1 - 1) :
			boundary_node(j, i, 2, -Nx_2 - Nx_3, Num5 + Num2 + Num3 + Num4, Coef, Coef_location);
			break;
		case (Node_1 + N_a - 4) ://2 lei shang 1j
			kind2_node(j, i, -Nx_2 - Nx_3, -2, -1, 0, 1, Num1 - Node_1 + Num2 + Num3 + Num4 + Node_5,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_11 + 1) ://2 lei xia 1j
			kind2_node(j, i, -(Nx_2 + Nx_3), -1, 0, 1, 2, Num5 + Num2 + Num3 + Num4,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3) < (Node_1 + N_a - 4))// 4dx!
			{
				inner_5node(j, i, 4, -Nx_3 - Nx_2, Num1 - Node_1 + Num2 + Num3 + Num4 + Node_5, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3) > (Node_1 + N_a - 4)) && ((i - Nx_2 - Nx_3) < Node13))// 2dx!
			{
				inner_5node(j, i, 2, -(Node_11 + 1), Num1 - Node_1 - N_a + 3 + Num2 + Num3 / 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3) >= Node13) && ((i - Nx_2 - Nx_3) <= Node14))// 2 lei left 2j !
			{
				kind2_node(j, i, -Nx_3 - Nx_2, -(Node_11 + 1), -1, 0, 1, Num1 - (i - Nx_2 - Nx_3) + Num2 / 2 + ((i - Nx_2 - Nx_3) - Node13) * 2,
					B5, B4, B3, B1_2, B3, B2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3) >= Node14) && ((i - Nx_2 - Nx_3) <= Node_11))// 2dx
			{
				inner_5node(j, i, 2, -(Node_11 + 1), Num1 - Node_11 + Num2 + Num3 - 1, Coef, Coef_location);
			}
			else
			{
				inner_5node(j, i, 4, -(Nx_3 + Nx_2), Num5 + Num2 + Num3 + Num4, Coef, Coef_location);// 4dx!
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!!!!!!!!!
	for (i = (Nx_2 + Nx_3 + Num1); i<(Nx_2 + Nx_3 + Num1 + Num2 / 2); i++)// num 2 :  26 node! part1
	{
		switch (i - Nx_2 - Nx_3 - Num1)
		{
		case 0:// 3 lei 2j shang
			kind3_node(j, i, -(Num1 - Node11) - 1, -(Num1 - Node11), 0, 1, Num2 + 5 - 2, Num2 + 5,
				C4, C3_2, C1_2, C2_2, C4, C3_2, Coef, Coef_location);
			break;
		case (Num2 / 2 - 1) ://node2: 3 lei 2j xia
			kind3_node(j, i, -(Num1 - Node12 + Num2 / 2 - 1), -(Num1 - Node12 + Num2 / 2 - 1) + 1, -1, 0, Num2 + 5, Num2 + 5 + 2,
			C3_2, C4, C2_2, C1_2, C3_2, C4, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1) % 2 == 0)// 1dx
			{
				inner_5node(j, i, 1, -((i - Nx_2 - Nx_3 - Num1) + Num1 - Node11 - (i - Nx_2 - Nx_3 - Num1) / 2), Num2 + 5, Coef, Coef_location);
			}
			else//1 lei  2j left
			{
				kind1_node(j, i, -((i - Nx_2 - Nx_3 - Num1) + Num1 - Node11 - (i - Nx_2 - Nx_3 - Num1) / 2) - (Nx_3 + N_d0 - 2),
					-((i - Nx_2 - Nx_3 - Num1) + Num1 - Node11 - (i - Nx_2 - Nx_3 - Num1) / 2) - (Nx_3 + N_d0 - 2) + 1,
					-((i - Nx_2 - Nx_3 - Num1) + Num1 - Node11 - (i - Nx_2 - Nx_3 - Num1) / 2),
					-((i - Nx_2 - Nx_3 - Num1) + Num1 - Node11 - (i - Nx_2 - Nx_3 - Num1) / 2) + 1,
					-1, 0, 1, 5 + Num2, A2_2, A2_2, A1_2, A1_2, A3_2, A5_2, A3_2, A4_2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 / 2); i<(Nx_2 + Nx_3 + Num1 + Num2); i++)// num 2 :  26 node! part2
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 / 2)
		{
		case 0:// 3 lei 2j shang
			kind3_node(j, i, -(Num1 - Node13 + Num2 / 2) - 1, -(Num1 - Node13 + Num2 / 2), 0, 1, Num3 - 5 - 2, Num3 - 5,
				C4, C3_2, C1_2, C2_2, C4, C3_2, Coef, Coef_location);
			break;
		case (Num2 / 2 - 1) :// 3 lei 2j xia
			kind3_node(j, i, -(Num1 - Node14 + Num2 - 1), -(Num1 - Node14 + Num2 - 1) + 1, -1, 0, Num3 - 5, Num3 - 5 + 2,
			C3_2, C4, C2_2, C1_2, C3_2, C4, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 / 2) % 2 == 0)// 1dx
			{
				inner_5node(j, i, 1, -((i - Nx_2 - Nx_3 - Num1 - Num2 / 2) + Num1 - Node13 - (i - Nx_2 - Nx_3 - Num1 - Num2 / 2) / 2) - Num2 / 2,
					Num3 - 5, Coef, Coef_location);
			}
			else//1 lei 2j left
			{
				kind1_node(j, i, -((i - Nx_2 - Nx_3 - Num1 - Num2 / 2) + Num1 - Node13 - (i - Nx_2 - Nx_3 - Num1 - Num2 / 2) / 2) - Num2 / 2 - (Node_11 + 1),
					-((i - Nx_2 - Nx_3 - Num1 - Num2 / 2) + Num1 - Node13 - (i - Nx_2 - Nx_3 - Num1 - Num2 / 2) / 2) - Num2 / 2 - (Node_11 + 1) + 1,
					-((i - Nx_2 - Nx_3 - Num1 - Num2 / 2) + Num1 - Node13 - (i - Nx_2 - Nx_3 - Num1 - Num2 / 2) / 2) - Num2 / 2,
					-((i - Nx_2 - Nx_3 - Num1 - Num2 / 2) + Num1 - Node13 - (i - Nx_2 - Nx_3 - Num1 - Num2 / 2) / 2) - Num2 / 2 + 1,
					-1, 0, 1, Num3 - 5, A2_2, A2_2, A1_2, A1_2, A3_2, A5_2, A3_2, A4_2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!!!
	for (i = (Nx_2 + Nx_3 + Num1 + Num2); i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 / 2); i++)// num 3 :  46 node! part1
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2)
		{
		case 0:// 1 lei 1j shang
			kind1_node(j, i, -(Num1 - N_d0 + 2 + Num2) - 2, -(Num1 - N_d0 + 2 + Num2) - 1, -(Num1 - N_d0 + 2 + Num2), 0, 1, Num3 + Num4 + N_d0 - 2 - 2,
				Num3 + Num4 + N_d0 - 2 - 1, Num3 + Num4 + N_d0 - 2, A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case 3:// 2 lei 2j shang
			kind2_node(j, i, -(Num1 - N_d0 + 2 + Num2), -2, -1, 0, 1, Num3 + Num4 + N_d0 - 2,
				B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case 4:// 3 lei 2j shang
			kind3_node(j, i, -(Num1 - N_d0 + 2 + Num2) - 1, -(Num1 - N_d0 + 2 + Num2), -1, 0, 1, Num3 - 4,
				C4, C4, C3_2, C1_2, C3_2, C2_2, Coef, Coef_location);
			break;
		case (Num3 / 2 - 5) :// 3 lei 2j xia
			kind3_node(j, i, -(Num1 - Node_1 + Num2 + Num3 / 2 - 1), -(Num1 - Node_1 + Num2 + Num3 / 2 - 1) + 1, -1, 0, 1, Num3 - 4,
			C4, C4, C3_2, C1_2, C3_2, C2_2, Coef, Coef_location);
			break;
		case (Num3 / 2 - 4) :// 2 lei 2j xia
			kind2_node(j, i, -(Num1 - Node_1 + Num2 + Num3 / 2 - 1), -1, 0, 1, 2, Num3 / 2 + 1 + Num4 + Node_5,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num3 / 2 - 1) :// 1 lei 1j xia
			kind1_node(j, i, -(Num1 - Node_1 + Num2 + Num3 / 2 - 1), -(Num1 - Node_1 + Num2 + Num3 / 2 - 1) + 1, -(Num1 - Node_1 + Num2 + Num3 / 2 - 1) + 2, -1, 0,
			Num3 / 2 + 1 + Num4 + Node_5, Num3 / 2 + 1 + Num4 + Node_5 + 1, Num3 / 2 + 1 + Num4 + Node_5 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2) < 3)// 2dx
			{
				inner_5node(j, i, 2, -(Num1 - N_d0 + 2 + Num2), Num3 + Num4 + N_d0 - 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2) > 4) && ((i - Nx_2 - Nx_3 - Num1 - Num2) < (Num3 / 2 - 5)))// 1dx
			{
				inner_5node(j, i, 1, -(Num2 + 5), Num3 - 4, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, i, 2, -(Num1 - Node_1 + Num2 + Num3 / 2 - 1), Num3 / 2 + 1 + Num4 + Node_5, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 / 2); i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3); i++)// num 3 :  46 node! part2
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 / 2)
		{
		case 0:// 1 lei 1j shang
			kind1_node(j, i, -(Num1 - Node_1 - N_a + 3 + Num2 + Num3 / 2) - 2, -(Num1 - Node_1 - N_a + 3 + Num2 + Num3 / 2) - 1, -(Num1 - Node_1 - N_a + 3 + Num2 + Num3 / 2), 0, 1,
				Num3 / 2 + Num4 + Node_5 + N_a - 3 - 2, Num3 / 2 + Num4 + Node_5 + N_a - 3 - 1, Num3 / 2 + Num4 + Node_5 + N_a - 3,
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case 3:// 2 lei 2j shang
			kind2_node(j, i, -(Num1 - Node_1 - N_a + 3 + Num2 + Num3 / 2), -2, -1, 0, 1, Num3 / 2 + Num4 + Node_5 + N_a - 3,
				B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case 4:// 3 lei 2j shang
			kind3_node(j, i, -(Num1 - Node_1 - N_a + 3 + Num2 + Num3 / 2) - 1, -(Num1 - Node_1 - N_a + 3 + Num2 + Num3 / 2), -1, 0, 1, 5 + Num4 - 1,
				C4, C4, C3_2, C1_2, C3_2, C2_2, Coef, Coef_location);
			break;
		case (Num3 / 2 - 5) :// 3 lei 2j xia
			kind3_node(j, i, -(Num1 - Node_11 + Num2 + Num3 - 1), -(Num1 - Node_11 + Num2 + Num3 - 1) + 1, -1, 0, 1, 5 + Num4 - 1,
			C4, C4, C3_2, C1_2, C3_2, C2_2, Coef, Coef_location);
			break;
		case (Num3 / 2 - 4) :// 2 lei 2j xia
			kind2_node(j, i, -(Num1 - Node_11 + Num2 + Num3 - 1), -1, 0, 1, 2, 1 + Num4 + Node_55,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num3 / 2 - 1) :// 1 lei 1j xia
			kind1_node(j, i, -(Num1 - Node_11 + Num2 + Num3 - 1), -(Num1 - Node_11 + Num2 + Num3 - 1) + 1, -(Num1 - Node_11 + Num2 + Num3 - 1) + 2, -1, 0,
			1 + Num4 + Node_55, 1 + Num4 + Node_55 + 1, 1 + Num4 + Node_55 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 / 2) < 3)// 2dx
			{
				inner_5node(j, i, 2, -(Num1 - Node_1 - N_a + 3 + Num2 + Num3 / 2), Num3 / 2 + Num4 + Node_5 + N_a - 3, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 / 2) > 4) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 / 2) < (Num3 / 2 - 5)))// 1dx
			{
				inner_5node(j, i, 1, -(Num3 - 5), 5 + Num4 - 1, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, i, 2, -(Num1 - Node_11 + Num2 + Num3 - 1), 1 + Num4 + Node_55, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3); i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 / 2); i++)// num 4 :  30 node! part1
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3)
		{
		case 0:// 1 lei 2j shang
			kind1_node(j, i, -(Num3 - 4) - 2, -(Num3 - 4) - 1, -(Num3 - 4), 0, 1, Num4 + N_d0 + 2 - 2, Num4 + N_d0 + 2 - 1, Num4 + N_d0 + 2,
				A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case (Num4 / 2 - 1) :// 1 lei 2j xia
			kind1_node(j, i, -(Num3 - 4), -(Num3 - 4) + 1, -(Num3 - 4) + 2, -1, 0,
			Num4 + N_d0 + 2 - N_c1, Num4 + N_d0 + 2 - N_c1 + 1, Num4 + N_d0 + 2 - N_c1 + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3) <= Node41)// 1dx
			{
				inner_5node(j, i, 1, -(Num3 - 4), Num4 + N_d0 + 2, Coef, Coef_location);
			}
			else if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3) >= Node42) // 1dx
			{
				inner_5node(j, i, 1, -(Num3 - 4), Num4 + N_d0 + 2 - N_c1, Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_4node(j, i, -(Num3 - 4), -1, 0, 1, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 / 2); i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4); i++)// num 4 :  30 node! part2
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3)
		{
		case (Num4 / 2) :// 1 lei 2j shang
			kind1_node(j, i, -(5 + Num4 - 1) - 2, -(5 + Num4 - 1) - 1, -(5 + Num4 - 1), 0, 1,
			(Node_55 - 4 + 1 + N_c1) - 2, (Node_55 - 4 + 1 + N_c1) - 1, (Node_55 - 4 + 1 + N_c1),
			A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case (Num4 - 1) :// 1 lei 2j xia
			kind1_node(j, i, -(5 + Num4 - 1), -(5 + Num4 - 1) + 1, -(5 + Num4 - 1) + 2, -1, 0,
			(Node_55 - 4 + 1), (Node_55 - 4 + 1) + 1, (Node_55 - 4 + 1) + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3) <= Node43)// 1dx
			{
				inner_5node(j, i, 1, -(5 + Num4 - 1), (Node_55 - 4 + 1 + N_c1), Coef, Coef_location);
			}
			else if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3) >= Node44) // 1dx
			{
				inner_5node(j, i, 1, -(5 + Num4 - 1), (Node_55 - 4 + 1), Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_4node(j, i, -(5 + Num4 - 1), -1, 0, 1, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4); i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Node_5 + 2); i++)// circle num5 :  80 node!  part1 !
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4)
		{
		case 0:
			boundary_node(j, i, 1, -Num1 - Num2 - Num3 - Num4, Num5 + Num6 + Num7 + Num8, Coef, Coef_location);
			break;
		case (N_d0 - 3) ://2 lei shang 1j
			kind2_node(j, i, -Num1 - Num2 - Num3 - Num4, -2, -1, 0, 1, Num5 + Num6 + Num7 + Num8,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (N_d0 + 1) :// 2 lei 2j shang
			kind2_node(j, i, -(Num3 + Num4 + N_d0 - 2), -2, -1, 0, 1, Num5 - N_d0 + 2 + Num6,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node51:// inner 3 node
			inner_3node(j, i, -(Num4 + N_d0 + 2), -1, 0, 1, 1, k*k*dx*dx - 4, Coef, Coef_location);
			break;
		case (Node51 + 1) :// inner 3 node
			inner_3node(j, i, -(Num4 + N_d0 + 2) + N_c1, 0, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node_5 - 3) ://2 lei xia 2j
			kind2_node(j, i, -(Num3 / 2 + 1 + Num4 + Node_5), -1, 0, 1, 2, Num5 - Node_5 + Num6 + Num7 / 2 - 1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_5 + 1) ://2 lei xia 1j
			kind2_node(j, i, -(Num1 - Node_1 + Num2 + Num3 + Num4 + Node_5), -1, 0, 1, 2, Num5 - Node_5 + Num6 + Num7 + Num8 + Node_9,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < (N_d0 - 3))// 4dx!
			{
				inner_5node(j, i, 4, -Num1 - Num2 - Num3 - Num4, Num5 + Num6 + Num7 + Num8, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) > (N_d0 - 3)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < (N_d0 + 1)))// 2dx!
			{
				inner_5node(j, i, 2, -(Num3 + Num4 + N_d0 - 2), Num5 - N_d0 + 2 + Num6, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) > (N_d0 + 1)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < Node51))// 1dx
			{
				inner_5node(j, i, 1, -(Num4 + N_d0 + 2), Num5 - N_d0 - 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) > (Node51 + 1)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < (Node_5 - 3)))// 1dx
			{
				inner_5node(j, i, 1, -(Num4 + N_d0 + 2) + N_c1, Num5 - N_d0 - 2 - 2, Coef, Coef_location);
			}
			else// 2dx!
			{
				inner_5node(j, i, 2, -(Num3 / 2 + 1 + Num4 + Node_5), Num5 - Node_5 + Num6 + Num7 / 2 - 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Node_5 + 2); i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5); i++)// circle num5 :  80 node!  part2 !
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4)
		{
		case (Num5 - 1) :
			boundary_node(j, i, 2, -Num2 - Num3 - Num4 - Num5, Num6 + Num7 + Num8 + Num9, Coef, Coef_location);
			break;
		case (Node_5 + N_a - 4) ://2 lei shang 1j
			kind2_node(j, i, -(Num1 - Node_1 + Num2 + Num3 + Num4 + Node_5), -2, -1, 0, 1, Num5 - Node_5 + Num6 + Num7 + Num8 + Node_9,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_5 + N_a) :// 2 lei 2j shang
			kind2_node(j, i, -(Num3 / 2 + Num4 + Node_5 + N_a - 3), -2, -1, 0, 1, Num5 - Node_5 - N_a + 3 + Num6 + Num7 / 2,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node52:// inner 3 node
			inner_3node(j, i, -(Node_55 - 4 + 1) - N_c1, -1, 0, 1, 1, k*k*dx*dx - 4, Coef, Coef_location);
			break;
		case (Node52 + 1) :// inner 3 node
			inner_3node(j, i, -(Node_55 - 4 + 1), 0, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node_55 - 3) ://2 lei xia 2j
			kind2_node(j, i, -(Num4 + Node_55 + 1), -1, 0, 1, 2, Num5 - Node_55 + Num6 + Num7 - 1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_55 + 1) ://2 lei xia 1j
			kind2_node(j, i, -Num2 - Num3 - Num4 - Num5, -1, 0, 1, 2, Num6 + Num7 + Num8 + Num9,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < (Node_5 + N_a - 4))// 4dx!
			{
				inner_5node(j, i, 4, -(Num1 - Node_1 + Num2 + Num3 + Num4 + Node_5), Num5 - Node_5 + Num6 + Num7 + Num8 + Node_9, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) > (Node_5 + N_a - 4)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < (Node_5 + N_a)))// 2dx!
			{
				inner_5node(j, i, 2, -(Num3 / 2 + Num4 + Node_5 + N_a - 3), Num5 - Node_5 - N_a + 3 + Num6 + Num7 / 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) > (Node_5 + N_a)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < Node52))// 1dx
			{
				inner_5node(j, i, 1, -(Node_55 - 4 + 1) - N_c1, (Num5 - Node_55 + 4 + Num6 - 1) + 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) > (Node52 + 1)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < (Node_55 - 3)))// 1dx
			{
				inner_5node(j, i, 1, -(Node_55 - 4 + 1), (Num5 - Node_55 + 4 + Num6 - 1), Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) > (Node_55 - 3)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < (Node_55 + 1)))// 2dx
			{
				inner_5node(j, i, 2, -(Num4 + Node_55 + 1), Num5 - Node_55 + Num6 + Num7 - 1, Coef, Coef_location);
			}
			else// 4dx!
			{
				inner_5node(j, i, 4, -Num2 - Num3 - Num4 - Num5, Num6 + Num7 + Num8 + Num9, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5); i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6); i++)// num 6 :  16 node!
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5)
		{
		case 0:// 1 lei 2j shang
			kind1_node(j, i, -(Num5 - N_d0 - 2) - 2, -(Num5 - N_d0 - 2) - 1, -(Num5 - N_d0 - 2), 0, 1, Num6 + 4 - 2, Num6 + 4 - 1, Num6 + 4,
				A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case Node61:// inner 3 node
			inner_3node(j, i, -(Num5 - N_d0 - 2), -1, 0, 1, 1, k*k*dx*dx - 4, Coef, Coef_location);
			break;
		case (Node61 + 1) :// inner 3 node
			inner_3node(j, i, -(Num5 - N_d0 - 2 - 2), 0, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Num6 / 2 - 1) :// 1 lei 2j xia
			kind1_node(j, i, -(Num5 - N_d0 - 2 - 2), -(Num5 - N_d0 - 2 - 2) + 1, -(Num5 - N_d0 - 2 - 2) + 2, -1, 0,
			Num6 + 4 - 2, Num6 + 4 - 2 + 1, Num6 + 4 - 2 + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		case (Num6 / 2) :// 1 lei 2j shang
			kind1_node(j, i, -(Num5 - Node_55 + 4 + Num6 - 1) - 2 - 2, -(Num5 - Node_55 + 4 + Num6 - 1) - 2 - 1, -(Num5 - Node_55 + 4 + Num6 - 1) - 2, 0, 1,
			Num7 - 4 + 2 - 2, Num7 - 4 + 2 - 1, Num7 - 4 + 2, A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case Node62:// inner 3 node
			inner_3node(j, i, -(Num5 - Node_55 + 4 + Num6 - 1) - 2, -1, 0, 1, 1, k*k*dx*dx - 4, Coef, Coef_location);
			break;
		case (Node62 + 1) :// inner 3 node
			inner_3node(j, i, -(Num5 - Node_55 + 4 + Num6 - 1), 0, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Num6 - 1) :// 1 lei 2j xia
			kind1_node(j, i, -(Num5 - Node_55 + 4 + Num6 - 1), -(Num5 - Node_55 + 4 + Num6 - 1) + 1, -(Num5 - Node_55 + 4 + Num6 - 1) + 2, -1, 0,
			Num7 - 4, Num7 - 4 + 1, Num7 - 4 + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5) < Node61)// 1dx
			{
				inner_5node(j, i, 1, -(Num5 - N_d0 - 2), Num6 + 4, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5) > (Node61 + 1)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5) < (Num6 / 2 - 1)))// 2dx
			{
				inner_5node(j, i, 1, -(Num5 - N_d0 - 2 - 2), Num6 + 4 - 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5) > (Num6 / 2)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5) < Node62))// 2dx
			{
				inner_5node(j, i, 1, -(Num5 - Node_55 + 4 + Num6 - 1) - 2, Num7 - 4 + 2, Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_5node(j, i, 1, -(Num5 - Node_55 + 4 + Num6 - 1), Num7 - 4, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!!
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6); i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 + Num7 / 2); i++)// num 7 :  28 node! part1
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6)
		{
		case 0:// 1 lei 1j shang
			kind1_node(j, i, -(Num5 - N_d0 + 2 + Num6) - 2, -(Num5 - N_d0 + 2 + Num6) - 1, -(Num5 - N_d0 + 2 + Num6), 0, 1,
				Num7 + Num8 + N_d0 - 2 - 2, Num7 + Num8 + N_d0 - 2 - 1, Num7 + Num8 + N_d0 - 2,
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case 3:// 2 lei 2j shang
			kind2_node(j, i, -(Num5 - N_d0 + 2 + Num6), -2, -1, 0, 1, Num7 + Num8 + N_d0 - 2,
				B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node71://inner 4 node
			inner_4node(j, i, -(Num6 + 4), -1, 0, Num7 - 4, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node71 + 1) ://inner 4 node
			inner_4node(j, i, -(Num6 + 4) + 2, 0, 1, Num7 - 4, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Num7 / 2 - 4) :// 2 lei 2j xia
			kind2_node(j, i, -(Num5 - Node_5 + Num6 + Num7 / 2 - 1), -1, 0, 1, 2, Num7 / 2 + 1 + Num8 + Node_9,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num7 / 2 - 1) :// 1 lei 1j xia
			kind1_node(j, i, -(Num5 - Node_5 + Num6 + Num7 / 2 - 1), -(Num5 - Node_5 + Num6 + Num7 / 2 - 1) + 1, -(Num5 - Node_5 + Num6 + Num7 / 2 - 1) + 2, -1, 0,
			Num7 / 2 + 1 + Num8 + Node_9, Num7 / 2 + 1 + Num8 + Node_9 + 1, Num7 / 2 + 1 + Num8 + Node_9 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6) < 3)// 2dx
			{
				inner_5node(j, i, 2, -(Num5 - N_d0 + 2 + Num6), Num7 + Num8 + N_d0 - 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6) > 3) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6) < (Node71)))// 1dx
			{
				inner_5node(j, i, 1, -(Num6 + 4), Num7 - 4, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6) > (Node71 + 1)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6) < (Num7 / 2 - 4)))// 1dx
			{
				inner_5node(j, i, 1, -(Num6 + 4) + 2, Num7 - 4, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, i, 2, -(Num5 - Node_5 + Num6 + Num7 / 2 - 1), Num7 / 2 + 1 + Num8 + Node_9, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 + Num7 / 2);
		i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 + Num7); i++)// num 7 :  28 node! part2
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6)
		{
		case (Num7 / 2) :// 1 lei 1j shang
			kind1_node(j, i, -(Num5 - Node_5 - N_a + 3 + Num6 + Num7 / 2) - 2, -(Num5 - Node_5 - N_a + 3 + Num6 + Num7 / 2) - 1, -(Num5 - Node_5 - N_a + 3 + Num6 + Num7 / 2), 0, 1,
			Num7 / 2 + Num8 + Node_9 + N_a - 3 - 2, Num7 / 2 + Num8 + Node_9 + N_a - 3 - 1, Num7 / 2 + Num8 + Node_9 + N_a - 3,
			A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case (Num7 / 2 + 3) :// 2 lei 2j shang
			kind2_node(j, i, -(Num5 - Node_5 - N_a + 3 + Num6 + Num7 / 2), -2, -1, 0, 1, Num7 / 2 + Num8 + Node_9 + N_a - 3,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node72://inner 4 node
			inner_4node(j, i, -(Num7 - 5 + 1) - 2, -1, 0, 5 + Num8 - 1, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node72 + 1) ://inner 4 node
			inner_4node(j, i, -(Num7 - 5 + 1), 0, 1, 5 + Num8 - 1, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Num7 - 4) :// 2 lei 2j xia
			kind2_node(j, i, -(Num5 - Node_55 + Num6 + Num7 - 1), -1, 0, 1, 2, 1 + Num8 + Node_99,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num7 - 1) :// 1 lei 1j xia
			kind1_node(j, i, -(Num5 - Node_55 + Num6 + Num7 - 1), -(Num5 - Node_55 + Num6 + Num7 - 1) + 1, -(Num5 - Node_55 + Num6 + Num7 - 1) + 2, -1, 0,
			1 + Num8 + Node_99, 1 + Num8 + Node_99 + 1, 1 + Num8 + Node_99 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6) < (Num7 / 2 + 3))// 2dx
			{
				inner_5node(j, i, 2, -(Num5 - Node_5 - N_a + 3 + Num6 + Num7 / 2), Num7 / 2 + Num8 + Node_9 + N_a - 3, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6) > (Num7 / 2 + 3)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6) < (Node72)))// 1dx
			{
				inner_5node(j, i, 1, -(Num7 - 5 + 1) - 2, 5 + Num8 - 1, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6) > (Node72 + 1)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6) < (Num7 - 4)))// 1dx
			{
				inner_5node(j, i, 1, -(Num7 - 5 + 1), 5 + Num8 - 1, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, i, 2, -(Num5 - Node_55 + Num6 + Num7 - 1), 1 + Num8 + Node_99, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 + Num7);
		i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 + Num7 + Num8); i++)// num 8 :  12 node! 
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 - Num7)
		{
		case 0:// 1 lei shang 2j
			kind1_node(j, i, -(Num7 - 4) - 2, -(Num7 - 4) - 1, -(Num7 - 4), 0, 1, Num8 + N_d0 + 2 - 2, Num8 + N_d0 + 2 - 1, Num8 + N_d0 + 2,
				A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case 1:
			inner_5node(j, i, 1, -(Num7 - 4), Num8 + N_d0 + 2, Coef, Coef_location);
			break;
		case 2://inner 4 node
			inner_4node(j, i, -(Num7 - 4), -1, 0, Num8 + N_d0 + 2, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case 3://inner 4 node
			inner_4node(j, i, -(Num7 - 4), 0, 1, Num8 + N_d0 + 2, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case 4:
			inner_5node(j, i, 1, -(Num7 - 4), Num8 + N_d0 + 2, Coef, Coef_location);
			break;
		case (Num8 / 2 - 1) :// 1 lei xia 2j
			kind1_node(j, i, -(Num7 - 4), -(Num7 - 4) + 1, -(Num7 - 4) + 2, -1, 0, Num8 + N_d0 + 2, Num8 + N_d0 + 2 + 1, Num8 + N_d0 + 2 + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		case (Num8 / 2) :// 1 lei shang 2j
			kind1_node(j, i, -(5 + Num8 - 1) - 2, -(5 + Num8 - 1) - 1, -(5 + Num8 - 1), 0, 1,
			Node_99 - 4 + 1 - 2, Node_99 - 4 + 1 - 1, Node_99 - 4 + 1,
			A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case 7:
			inner_5node(j, i, 1, -(5 + Num8 - 1), Node_99 - 4 + 1, Coef, Coef_location);
			break;
		case 8://inner 4 node
			inner_4node(j, i, -(5 + Num8 - 1), -1, 0, Node_99 - 4 + 1, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case 9://inner 4 node
			inner_4node(j, i, -(5 + Num8 - 1), 0, 1, Node_99 - 4 + 1, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case 10:
			inner_5node(j, i, 1, -(5 + Num8 - 1), Node_99 - 4 + 1, Coef, Coef_location);
			break;
		case (Num8 - 1) :// 1 lei xia 2j
			kind1_node(j, i, -(5 + Num8 - 1), -(5 + Num8 - 1) + 1, -(5 + Num8 - 1) + 2, -1, 0,
			Node_99 - 4 + 1, Node_99 - 4 + 1 + 1, Node_99 - 4 + 1 + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 + Num7 + Num8);
		i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 + Num7 + Num8 + Node_9 + 2); i++)// circle num9 :  72 node!  part1 !
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 - Num7 - Num8)
		{
		case 0:
			boundary_node(j, i, 1, -Num5 - Num6 - Num7 - Num8, Num9 + Num8 + Num7 + Num6, Coef, Coef_location);
			break;
		case (N_d0 - 3) ://2 lei shang 1j
			kind2_node(j, i, -Num5 - Num6 - Num7 - Num8, -2, -1, 0, 1, Num9 + Num8 + Num7 + Num6,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (N_d0 + 1) :// 2 lei 2j shang
			kind2_node(j, i, -(Num7 + Num8 + N_d0 - 2), -2, -1, 0, 1, Num9 - N_d0 + 2 + Num8,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node91:// inner 4 node
			inner_4node(j, i, -(Num8 + N_d0 + 2), -1, 0, Num9 - N_d0 - 2, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node91 + 1) :// inner 4 node
			inner_4node(j, i, -(Num8 + N_d0 + 2), 0, 1, Num9 - N_d0 - 2, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Node_9 - 3) ://2 lei xia 2j
			kind2_node(j, i, -(Num7 / 2 + 1 + Num8 + Node_9), -1, 0, 1, 2, Num9 - Node_9 + Num8 + Num7 / 2 - 1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_9 + 1) ://2 lei xia 1j
			kind2_node(j, i, -(Num5 - Node_5 + Num6 + Num7 + Num8 + Node_9), -1, 0, 1, 2, Num9 - Node_9 + Num8 + Num7 + Num6 + Node_5,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 - Num7 - Num8) < (N_d0 - 3))// 4dx!
			{
				inner_5node(j, i, 4, -Num5 - Num6 - Num7 - Num8, Num9 + Num8 + Num7 + Num6, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 - Num7 - Num8) > (N_d0 - 3)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 - Num7 - Num8) < (N_d0 + 1)))// 2dx!
			{
				inner_5node(j, i, 2, -(Num7 + Num8 + N_d0 - 2), Num9 - N_d0 + 2 + Num8, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 - Num7 - Num8) > (N_d0 + 1)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 - Num7 - Num8) < Node91))// 1dx
			{
				inner_5node(j, i, 1, -(Num8 + N_d0 + 2), Num9 - N_d0 - 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 - Num7 - Num8) > (Node91 + 1)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 - Num7 - Num8) < (Node_9 - 3)))// 1dx
			{
				inner_5node(j, i, 1, -(Num8 + N_d0 + 2), Num9 - N_d0 - 2, Coef, Coef_location);
			}
			else// 2dx!
			{
				inner_5node(j, i, 2, -(Num7 / 2 + 1 + Num8 + Node_9), Num9 - Node_9 + Num8 + Num7 / 2 - 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 + Num7 + Num8 + Node_9 + 2);
		i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 + Num7 + Num8 + Num9); i++)// circle num9 :  72 node!  part2 !
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 - Num7 - Num8)
		{
		case (Num9 - 1) :
			boundary_node(j, i, 2, -Num6 - Num7 - Num8 - Num9, Num8 + Num7 + Num6 + Num5, Coef, Coef_location);
			break;
		case (Node_9 + N_a - 4) ://2 lei shang 1j
			kind2_node(j, i, -(Num5 - Node_5 + Num6 + Num7 + Num8 + Node_9), -2, -1, 0, 1, Num9 - Node_9 + Num8 + Num7 + Num6 + Node_5,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_9 + N_a) :// 2 lei 2j shang
			kind2_node(j, i, -(Num7 / 2 + Num8 + Node_9 + N_a - 3), -2, -1, 0, 1, Num9 - Node_9 - N_a + 3 + Num8 + Num7 / 2,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node92:// inner 4 node
			inner_4node(j, i, -(Node_99 - 4 + 1), -1, 0, Num9 - Node_99 + 4 + Num8 - 1, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node92 + 1) :// inner 4 node
			inner_4node(j, i, -(Node_99 - 4 + 1), 0, 1, Num9 - Node_99 + 4 + Num8 - 1, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Node_99 - 3) ://2 lei xia 2j
			kind2_node(j, i, -(Num8 + Node_99 + 1), -1, 0, 1, 2, Num9 - Node_99 + Num8 + Num7 - 1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_99 + 1) ://2 lei xia 1j
			kind2_node(j, i, -Num6 - Num7 - Num8 - Num9, -1, 0, 1, 2, Num8 + Num7 + Num6 + Num5,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 - Num7 - Num8) < (Node_9 + N_a - 4))// 4dx!
			{
				inner_5node(j, i, 4, -(Num5 - Node_5 + Num6 + Num7 + Num8 + Node_9), Num9 - Node_9 + Num8 + Num7 + Num6 + Node_5, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 - Num7 - Num8) > (Node_9 + N_a - 4)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 - Num7 - Num8) < (Node_9 + N_a)))// 2dx!
			{
				inner_5node(j, i, 2, -(Num7 / 2 + Num8 + Node_9 + N_a - 3), Num9 - Node_9 - N_a + 3 + Num8 + Num7 / 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 - Num7 - Num8) > (Node_9 + N_a)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 - Num7 - Num8) < Node92))// 1dx
			{
				inner_5node(j, i, 1, -(Node_99 - 4 + 1), Num9 - Node_99 + 4 + Num8 - 1, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 - Num7 - Num8) > (Node92 + 1)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 - Num7 - Num8) < (Node_99 - 3)))// 1dx
			{
				inner_5node(j, i, 1, -(Node_99 - 4 + 1), Num9 - Node_99 + 4 + Num8 - 1, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 - Num7 - Num8) > (Node_99 - 3)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 - Num7 - Num8) < (Node_99 + 1)))// 2dx
			{
				inner_5node(j, i, 2, -(Num8 + Node_99 + 1), Num9 - Node_99 + Num8 + Num7 - 1, Coef, Coef_location);
			}
			else// 4dx!
			{
				inner_5node(j, i, 4, -Num6 - Num7 - Num8 - Num9, Num8 + Num7 + Num6 + Num5, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 + Num7 + Num8 + Num9);
		i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 + Num7 + Num8 * 2 + Num9); i++)// num 10(8) :  12 node! 
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 - Num7 - Num8 - Num9)
		{
		case 0:// 1 lei shang 2j
			kind1_node(j, i, -(Num9 - N_d0 - 2) - 2, -(Num9 - N_d0 - 2) - 1, -(Num9 - N_d0 - 2), 0, 1, Num8 + 4 - 2, Num8 + 4 - 1, Num8 + 4,
				A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case 1:
			inner_5node(j, i, 1, -(Num9 - N_d0 - 2), Num8 + 4, Coef, Coef_location);
			break;
		case 2://inner 4 node
			inner_4node(j, i, -(Num9 - N_d0 - 2), -1, 0, Num8 + 4, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case 3://inner 4 node
			inner_4node(j, i, -(Num9 - N_d0 - 2), 0, 1, Num8 + 4, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case 4:
			inner_5node(j, i, 1, -(Num9 - N_d0 - 2), Num8 + 4, Coef, Coef_location);
			break;
		case (Num8 / 2 - 1) :// 1 lei xia 2j
			kind1_node(j, i, -(Num9 - N_d0 - 2), -(Num9 - N_d0 - 2) + 1, -(Num9 - N_d0 - 2) + 2, -1, 0, Num8 + 4, Num8 + 4 + 1, Num8 + 4 + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		case (Num8 / 2) :// 1 lei shang 2j
			kind1_node(j, i, -(Num9 - Node_99 + 4 + Num8 - 1) - 2, -(Num9 - Node_99 + 4 + Num8 - 1) - 1, -(Num9 - Node_99 + 4 + Num8 - 1), 0, 1,
			Num7 - 4 - 2, Num7 - 4 - 1, Num7 - 4, A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case 7:
			inner_5node(j, i, 1, -(Num9 - Node_99 + 4 + Num8 - 1), Num7 - 4, Coef, Coef_location);
			break;
		case 8://inner 4 node
			inner_4node(j, i, -(Num9 - Node_99 + 4 + Num8 - 1), -1, 0, Num7 - 4, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case 9://inner 4 node
			inner_4node(j, i, -(Num9 - Node_99 + 4 + Num8 - 1), 0, 1, Num7 - 4, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case 10:
			inner_5node(j, i, 1, -(Num9 - Node_99 + 4 + Num8 - 1), Num7 - 4, Coef, Coef_location);
			break;
		case (Num8 - 1) :// 1 lei xia 2j
			kind1_node(j, i, -(Num9 - Node_99 + 4 + Num8 - 1), -(Num9 - Node_99 + 4 + Num8 - 1) + 1, -(Num9 - Node_99 + 4 + Num8 - 1) + 2, -1, 0,
			Num7 - 4, Num7 - 4 + 1, Num7 - 4 + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 + Num7 + Num8 * 2 + Num9);
		i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 + Num7 + Num8 * 2 + Num9 + Num7 / 2); i++)// num 11(7) :  28 node! part1
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 - Num7 - Num8 * 2 - Num9)
		{
		case 0:// 1 lei 1j shang
			kind1_node(j, i, -(Num9 - N_d0 + 2 + Num8) - 2, -(Num9 - N_d0 + 2 + Num8) - 1, -(Num9 - N_d0 + 2 + Num8), 0, 1,
				Num7 + Num6 + N_d0 - 2 - 2, Num7 + Num6 + N_d0 - 2 - 1, Num7 + Num6 + N_d0 - 2,
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case 3:// 2 lei 2j shang
			kind2_node(j, i, -(Num9 - N_d0 + 2 + Num8), -2, -1, 0, 1, Num7 + Num6 + N_d0 - 2,
				B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node71://inner 4 node
			inner_4node(j, i, -(Num8 + 4), -1, 0, Num7 - 4, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node71 + 1) ://inner 4 node
			inner_4node(j, i, -(Num8 + 4), 0, 1, Num7 - 4 + 2, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Num7 / 2 - 4) :// 2 lei 2j xia
			kind2_node(j, i, -(Num9 - Node_9 + Num8 + Num7 / 2 - 1), -1, 0, 1, 2, Num7 / 2 + 1 + Num6 + Node_5,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num7 / 2 - 1) :// 1 lei 1j xia
			kind1_node(j, i, -(Num9 - Node_9 + Num8 + Num7 / 2 - 1), -(Num9 - Node_9 + Num8 + Num7 / 2 - 1) + 1, -(Num9 - Node_9 + Num8 + Num7 / 2 - 1) + 2, -1, 0,
			Num7 / 2 + 1 + Num6 + Node_5, Num7 / 2 + 1 + Num6 + Node_5 + 1, Num7 / 2 + 1 + Num6 + Node_5 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 - Num7 - Num8 * 2 - Num9) < 3)// 2dx
			{
				inner_5node(j, i, 2, -(Num9 - N_d0 + 2 + Num8), Num7 + Num6 + N_d0 - 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 - Num7 - Num8 * 2 - Num9) > 3) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 - Num7 - Num8 * 2 - Num9) < (Node71)))// 1dx
			{
				inner_5node(j, i, 1, -(Num8 + 4), Num7 - 4, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 - Num7 - Num8 * 2 - Num9) > (Node71 + 1)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 - Num7 - Num8 * 2 - Num9) < (Num7 / 2 - 4)))// 1dx
			{
				inner_5node(j, i, 1, -(Num8 + 4), Num7 - 4 + 2, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, i, 2, -(Num9 - Node_9 + Num8 + Num7 / 2 - 1), Num7 / 2 + 1 + Num6 + Node_5, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 + Num7 + Num8 * 2 + Num9 + Num7 / 2);
		i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 + Num7 * 2 + Num8 * 2 + Num9); i++)// num 11(7) :  28 node! part2
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 - Num7 - Num8 * 2 - Num9)
		{
		case (Num7 / 2) :// 1 lei 1j shang
			kind1_node(j, i, -(Num9 - Node_9 - N_a + 3 + Num8 + Num7 / 2) - 2, -(Num9 - Node_9 - N_a + 3 + Num8 + Num7 / 2) - 1, -(Num9 - Node_9 - N_a + 3 + Num8 + Num7 / 2), 0, 1,
			Num7 / 2 + Num6 + Node_5 + N_a - 3 - 2, Num7 / 2 + Num6 + Node_5 + N_a - 3 - 1, Num7 / 2 + Num6 + Node_5 + N_a - 3,
			A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case (Num7 / 2 + 3) :// 2 lei 2j shang
			kind2_node(j, i, -(Num9 - Node_9 - N_a + 3 + Num8 + Num7 / 2), -2, -1, 0, 1, Num7 / 2 + Num6 + Node_5 + N_a - 3,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node72://inner 4 node
			inner_4node(j, i, -(Num7 - 5 + 1), -1, 0, 5 + Num6 - 1 - 2, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node72 + 1) ://inner 4 node
			inner_4node(j, i, -(Num7 - 5 + 1), 0, 1, 5 + Num6 - 1, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Num7 - 4) :// 2 lei 2j xia
			kind2_node(j, i, -(Num9 - Node_99 + Num8 + Num7 - 1), -1, 0, 1, 2, 1 + Num6 + Node_55,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num7 - 1) :// 1 lei 1j xia
			kind1_node(j, i, -(Num9 - Node_99 + Num8 + Num7 - 1), -(Num9 - Node_99 + Num8 + Num7 - 1) + 1, -(Num9 - Node_99 + Num8 + Num7 - 1) + 2, -1, 0,
			1 + Num6 + Node_55, 1 + Num6 + Node_55 + 1, 1 + Num6 + Node_55 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 - Num7 - Num8 * 2 - Num9) < (Num7 / 2 + 3))// 2dx
			{
				inner_5node(j, i, 2, -(Num9 - Node_9 - N_a + 3 + Num8 + Num7 / 2), Num7 / 2 + Num6 + Node_5 + N_a - 3, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 - Num7 - Num8 * 2 - Num9) > (Num7 / 2 + 3)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 - Num7 - Num8 * 2 - Num9) < (Node72)))// 1dx
			{
				inner_5node(j, i, 1, -(Num7 - 5 + 1), 5 + Num6 - 1 - 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 - Num7 - Num8 * 2 - Num9) > (Node72 + 1)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 - Num7 - Num8 * 2 - Num9) < (Num7 - 4)))// 1dx
			{
				inner_5node(j, i, 1, -(Num7 - 5 + 1), 5 + Num6 - 1, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, i, 2, -(Num9 - Node_99 + Num8 + Num7 - 1), 1 + Num6 + Node_55, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 + Num7 * 2 + Num8 * 2 + Num9);
		i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9); i++)// num 12(6) :  16 node!
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 - Num7 * 2 - Num8 * 2 - Num9)
		{
		case 0:// 1 lei 2j shang
			kind1_node(j, i, -(Num7 - 4) - 2, -(Num7 - 4) - 1, -(Num7 - 4), 0, 1, Num6 + N_d0 + 2 - 2, Num6 + N_d0 + 2 - 1, Num6 + N_d0 + 2,
				A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case Node61:// inner 3 node
			inner_3node(j, i, -1, 0, (Num6 + N_d0 + 2), 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node61 + 1) :// inner 3 node
			inner_3node(j, i, 0, 1, (Num6 + N_d0 + 2) + 2, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Num6 / 2 - 1) :// 1 lei 2j xia
			kind1_node(j, i, (-(Num7 - 4) - 2), (-(Num7 - 4) - 2) + 1, (-(Num7 - 4) - 2) + 2, -1, 0,
			(Num6 + N_d0 + 2) + 2, (Num6 + N_d0 + 2) + 2 + 1, (Num6 + N_d0 + 2) + 2 + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		case (Num6 / 2) :// 1 lei 2j shang
			kind1_node(j, i, -(5 + Num6 - 1) + 2 - 2, -(5 + Num6 - 1) + 2 - 1, -(5 + Num6 - 1) + 2, 0, 1,
			(Node_55 - 4 + 1) - 2 - 2, (Node_55 - 4 + 1) - 2 - 1, (Node_55 - 4 + 1) - 2,
			A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case Node62:// inner 3 node
			inner_3node(j, i, -1, 0, (Node_55 - 4 + 1) - 2, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node62 + 1) :// inner 3 node
			inner_3node(j, i, 0, 1, (Node_55 - 4 + 1), k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Num6 - 1) :// 1 lei 2j xia
			kind1_node(j, i, -(5 + Num6 - 1), -(5 + Num6 - 1) + 1, -(5 + Num6 - 1) + 2, -1, 0,
			(Node_55 - 4 + 1), (Node_55 - 4 + 1) + 1, (Node_55 - 4 + 1) + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 - Num7 * 2 - Num8 * 2 - Num9) < Node61)// 1dx
			{
				inner_5node(j, i, 1, -(Num7 - 4), Num6 + N_d0 + 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 - Num7 * 2 - Num8 * 2 - Num9) > (Node61 + 1)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 - Num7 * 2 - Num8 * 2 - Num9) < (Num6 / 2 - 1)))// 2dx
			{
				inner_5node(j, i, 1, -(Num7 - 4) - 2, (Num6 + N_d0 + 2) + 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 - Num7 * 2 - Num8 * 2 - Num9) > (Num6 / 2)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 - Num7 * 2 - Num8 * 2 - Num9) < Node62))// 2dx
			{
				inner_5node(j, i, 1, -(5 + Num6 - 1) + 2, (Node_55 - 4 + 1) - 2, Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_5node(j, i, 1, -(5 + Num6 - 1), (Node_55 - 4 + 1), Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9);
		i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Node_5 + 2); i++)// circle num13(5) :  80 node!  part1 !
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9)
		{
		case 0:
			boundary_node(j, i, 1, -Num9 - Num8 - Num7 - Num6, Num5 + Num4 + Num3 + Num2, Coef, Coef_location);
			break;
		case (N_d0 - 3) ://2 lei shang 1j
			kind2_node(j, i, -Num9 - Num8 - Num7 - Num6, -2, -1, 0, 1, Num5 + Num4 + Num3 + Num2,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (N_d0 + 1) :// 2 lei 2j shang
			kind2_node(j, i, -(Num7 + Num6 + N_d0 - 2), -2, -1, 0, 1, Num5 - N_d0 + 2 + Num4,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node51:// inner 3 node
			inner_3node(j, i, -1, 0, (Num5 - N_d0 - 2), 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node51 + 1) :// inner 3 node
			inner_3node(j, i, 0, 1, (Num5 - N_d0 - 2) + N_c1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Node_5 - 3) ://2 lei xia 2j
			kind2_node(j, i, -(Num7 / 2 + 1 + Num6 + Node_5), -1, 0, 1, 2, Num5 - Node_5 + Num4 + Num3 / 2 - 1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_5 + 1) ://2 lei xia 1j
			kind2_node(j, i, -(Num9 - Node_9 + Num8 + Num7 + Num6 + Node_5), -1, 0, 1, 2, Num5 - Node_5 + Num4 + Num3 + Num2 + Node_1,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < (N_d0 - 3))// 4dx!
			{
				inner_5node(j, i, 4, -Num9 - Num8 - Num7 - Num6, Num5 + Num4 + Num3 + Num2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) > (N_d0 - 3)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < (N_d0 + 1)))// 2dx!
			{
				inner_5node(j, i, 2, -(Num7 + Num6 + N_d0 - 2), Num5 - N_d0 + 2 + Num4, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) > (N_d0 + 1)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < Node51))// 1dx
			{
				inner_5node(j, i, 1, -(Num6 + N_d0 + 2), (Num5 - N_d0 - 2), Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) > (Node51 + 1)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < (Node_5 - 3)))// 1dx
			{
				inner_5node(j, i, 1, -(Num6 + N_d0 + 2) - 2, (Num5 - N_d0 - 2) + N_c1, Coef, Coef_location);
			}
			else// 2dx!
			{
				inner_5node(j, i, 2, -(Num7 / 2 + 1 + Num6 + Node_5), Num5 - Node_5 + Num4 + Num3 / 2 - 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Node_5 + 2);
		i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9); i++)// circle num13(5) :  80 node!  part2 !
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9)
		{
		case (Num5 - 1) :
			boundary_node(j, i, 2, -Num8 - Num7 - Num6 - Num5, Num4 + Num3 + Num2 + Num1, Coef, Coef_location);
			break;
		case (Node_5 + N_a - 4) ://2 lei shang 1j
			kind2_node(j, i, -(Num9 - Node_9 + Num8 + Num7 + Num6 + Node_5), -2, -1, 0, 1, Num5 - Node_5 + Num4 + Num3 + Num2 + Node_1,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_5 + N_a) :// 2 lei 2j shang
			kind2_node(j, i, -(Num7 / 2 + Num6 + Node_5 + N_a - 3), -2, -1, 0, 1, Num5 - Node_5 - N_a + 3 + Num4 + Num3 / 2,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node52:// inner 3 node
			inner_3node(j, i, -1, 0, (Num5 - Node_55 + 4 + Num4 - 1) - N_c1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node52 + 1) :// inner 3 node
			inner_3node(j, i, 0, 1, (Num5 - Node_55 + 4 + Num4 - 1), k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Node_55 - 3) ://2 lei xia 2j
			kind2_node(j, i, -(Num6 + Node_55 + 1), -1, 0, 1, 2, Num5 - Node_55 + Num4 + Num3 - 1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_55 + 1) ://2 lei xia 1j
			kind2_node(j, i, -Num8 - Num7 - Num6 - Num5, -1, 0, 1, 2, Num4 + Num3 + Num2 + Num1,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < (Node_5 + N_a - 4))// 4dx!
			{
				inner_5node(j, i, 4, -(Num9 - Node_9 + Num8 + Num7 + Num6 + Node_5), Num5 - Node_5 + Num4 + Num3 + Num2 + Node_1, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) > (Node_5 + N_a - 4)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < (Node_5 + N_a)))// 2dx!
			{
				inner_5node(j, i, 2, -(Num7 / 2 + Num6 + Node_5 + N_a - 3), Num5 - Node_5 - N_a + 3 + Num4 + Num3 / 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) > (Node_5 + N_a)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < Node52))// 1dx
			{
				inner_5node(j, i, 1, -(Node_55 - 4 + 1) + 2, (Num5 - Node_55 + 4 + Num4 - 1) - N_c1, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) > (Node52 + 1)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < (Node_55 - 3)))// 1dx
			{
				inner_5node(j, i, 1, -(Node_55 - 4 + 1), (Num5 - Node_55 + 4 + Num4 - 1), Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) > (Node_55 - 3)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < (Node_55 + 1)))// 2dx
			{
				inner_5node(j, i, 2, -(Num6 + Node_55 + 1), Num5 - Node_55 + Num4 + Num3 - 1, Coef, Coef_location);
			}
			else// 4dx!
			{
				inner_5node(j, i, 4, -Num8 - Num7 - Num6 - Num5, Num4 + Num3 + Num2 + Num1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9);
		i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Num4 / 2); i++)// num 14(4) :  30 node! part1
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9)
		{
		case 0:// 1 lei 2j shang
			kind1_node(j, i, -(Num5 - N_d0 - 2) - 2, -(Num5 - N_d0 - 2) - 1, -(Num5 - N_d0 - 2), 0, 1, Num4 + 4 - 2, Num4 + 4 - 1, Num4 + 4,
				A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case (Num4 / 2 - 1) :// 1 lei 2j xia
			kind1_node(j, i, -(Num5 - N_d0 - 2) - N_c1, -(Num5 - N_d0 - 2) - N_c1 + 1, -(Num5 - N_d0 - 2) - N_c1 + 2, -1, 0,
			Num4 + 4, Num4 + 4 + 1, Num4 + 4 + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) <= Node41)// 1dx
			{
				inner_5node(j, i, 1, -(Num5 - N_d0 - 2), Num4 + 4, Coef, Coef_location);
			}
			else if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) >= Node42) // 1dx
			{
				inner_5node(j, i, 1, -(Num5 - N_d0 - 2) - N_c1, Num4 + 4, Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_4node(j, i, -1, 0, 1, Num4 + 4, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Num4 / 2);
		i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9); i++)// num14(4) :  30 node! part2
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9)
		{
		case (Num4 / 2) :// 1 lei 2j shang
			kind1_node(j, i, -(Num5 - Node_55 + 4 + Num4 - 1) + N_c1 - 2, -(Num5 - Node_55 + 4 + Num4 - 1) + N_c1 - 1, -(Num5 - Node_55 + 4 + Num4 - 1) + N_c1, 0, 1,
			(Num3 - 5 + 1) - 2, (Num3 - 5 + 1) - 1, (Num3 - 5 + 1),
			A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case (Num4 - 1) :// 1 lei 2j xia
			kind1_node(j, i, -(Num5 - Node_55 + 4 + Num4 - 1), -(Num5 - Node_55 + 4 + Num4 - 1) + 1, -(Num5 - Node_55 + 4 + Num4 - 1) + 2, -1, 0,
			(Num3 - 5 + 1), (Num3 - 5 + 1) + 1, (Num3 - 5 + 1) + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) <= Node43)// 1dx
			{
				inner_5node(j, i, 1, -(Num5 - Node_55 + 4 + Num4 - 1) + N_c1, (Num3 - 5 + 1), Coef, Coef_location);
			}
			else if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) >= Node44) // 1dx
			{
				inner_5node(j, i, 1, -(Num5 - Node_55 + 4 + Num4 - 1), (Num3 - 5 + 1), Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_4node(j, i, -1, 0, 1, (Num3 - 5 + 1), 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9);
		i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Num3 / 2); i++)// num 15(3) :  46 node! part1
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9)
		{
		case 0:// 1 lei 1j shang
			kind1_node(j, i, -(Num5 - N_d0 + 2 + Num4) - 2, -(Num5 - N_d0 + 2 + Num4) - 1, -(Num5 - N_d0 + 2 + Num4), 0, 1,
				Num3 + Num2 + N_d0 - 2 - 2, Num3 + Num2 + N_d0 - 2 - 1, Num3 + Num2 + N_d0 - 2,
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case 3:// 2 lei 2j shang
			kind2_node(j, i, -(Num5 - N_d0 + 2 + Num4), -2, -1, 0, 1, Num3 + Num2 + N_d0 - 2,
				B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case 4:// 3 lei 2j shang
			kind3_node(j, i, -(Num4 + 4), -1, 0, 1, (Num3 + Num2 + N_d0 - 2) - 1, (Num3 + Num2 + N_d0 - 2),
				C2_2, C3_2, C1_2, C3_2, C4, C4, Coef, Coef_location);
			break;
		case (Num3 / 2 - 5) :// 3 lei 2j xia
			kind3_node(j, i, -(Num4 + 4), -1, 0, 1, Num3 / 2 + 1 + Num2 + Node_1, Num3 / 2 + 1 + Num2 + Node_1 + 1,
			C2_2, C3_2, C1_2, C3_2, C4, C4, Coef, Coef_location);
			break;
		case (Num3 / 2 - 4) :// 2 lei 2j xia
			kind2_node(j, i, -(Num5 - Node_5 + Num4 + Num3 / 2 - 1), -1, 0, 1, 2, Num3 / 2 + 1 + Num2 + Node_1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num3 / 2 - 1) :// 1 lei 1j xia
			kind1_node(j, i, -(Num5 - Node_5 + Num4 + Num3 / 2 - 1), -(Num5 - Node_5 + Num4 + Num3 / 2 - 1) + 1, -(Num5 - Node_5 + Num4 + Num3 / 2 - 1) + 2, -1, 0,
			Num3 / 2 + 1 + Num2 + Node_1, Num3 / 2 + 1 + Num2 + Node_1 + 1, Num3 / 2 + 1 + Num2 + Node_1 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < 3)// 2dx
			{
				inner_5node(j, i, 2, -(Num5 - N_d0 + 2 + Num4), Num3 + Num2 + N_d0 - 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) > 4) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < (Num3 / 2 - 5)))// 1dx
			{
				inner_5node(j, i, 1, -(Num4 + 4), Num3 - 5, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, i, 2, -(Num5 - Node_5 + Num4 + Num3 / 2 - 1), Num3 / 2 + 1 + Num2 + Node_1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Num3 / 2);
		i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 * 2 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9); i++)// num 15(3) :  46 node! part2
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num3 / 2)
		{
		case 0:// 1 lei 1j shang
			kind1_node(j, i, -(Num5 - Node_5 - N_a + 3 + Num4 + Num3 / 2) - 2, -(Num5 - Node_5 - N_a + 3 + Num4 + Num3 / 2) - 1, -(Num5 - Node_5 - N_a + 3 + Num4 + Num3 / 2), 0, 1,
				Num3 / 2 + Num2 + Node_1 + N_a - 3 - 2, Num3 / 2 + Num2 + Node_1 + N_a - 3 - 1, Num3 / 2 + Num2 + Node_1 + N_a - 3,
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case 3:// 2 lei 2j shang
			kind2_node(j, i, -(Num5 - Node_5 - N_a + 3 + Num4 + Num3 / 2), -2, -1, 0, 1, Num3 / 2 + Num2 + Node_1 + N_a - 3,
				B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case 4:// 3 lei 2j shang
			kind3_node(j, i, -(Num3 - 5 + 1), -1, 0, 1, (Num3 / 2 + Num2 + Node_1 + N_a - 3) - 1, (Num3 / 2 + Num2 + Node_1 + N_a - 3),
				C2_2, C3_2, C1_2, C3_2, C4, C4, Coef, Coef_location);
			break;
		case (Num3 / 2 - 5) :// 3 lei 2j xia
			kind3_node(j, i, -(Num3 - 5 + 1), -1, 0, 1, 1 + Num2 + Node_11, 1 + Num2 + Node_11 + 1,
			C2_2, C3_2, C1_2, C3_2, C4, C4, Coef, Coef_location);
			break;
		case (Num3 / 2 - 4) :// 2 lei 2j xia
			kind2_node(j, i, -(Num5 - Node_55 + Num4 + Num3 - 1), -1, 0, 1, 2, 1 + Num2 + Node_11,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num3 / 2 - 1) :// 1 lei 1j xia
			kind1_node(j, i, -(Num5 - Node_55 + Num4 + Num3 - 1), -(Num5 - Node_55 + Num4 + Num3 - 1) + 1, -(Num5 - Node_55 + Num4 + Num3 - 1) + 2, -1, 0,
			1 + Num2 + Node_11, 1 + Num2 + Node_11 + 1, 1 + Num2 + Node_11 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num3 / 2) < 3)// 2dx
			{
				inner_5node(j, i, 2, -(Num5 - Node_5 - N_a + 3 + Num4 + Num3 / 2), Num3 / 2 + Num2 + Node_1 + N_a - 3, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num3 / 2) > 4) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num3 / 2) < (Num3 / 2 - 5)))// 1dx
			{
				inner_5node(j, i, 1, -(Num3 - 5 + 1), 6 + Num2 - 1, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, i, 2, -(Num5 - Node_55 + Num4 + Num3 - 1), 1 + Num2 + Node_11, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 * 2 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9);
		i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 * 2 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Num2 / 2); i++)// num 16(2) :  26 node! part1
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9)
		{
		case 0:// 3 lei 2j shang
			kind3_node(j, i, -(Num3 - 3), -(Num3 - 3) + 2, 0, 1, Num2 + Node11 - 1, Num2 + Node11,
				C4, C3_2, C1_2, C2_2, C4, C3_2, Coef, Coef_location);
			break;
		case (Num2 / 2 - 1) ://node2: 3 lei 2j xia
			kind3_node(j, i, -(Num3 - 5), -(Num3 - 5) + 2, -1, 0, Num2 / 2 + Node_1 - 4 + 1, Num2 / 2 + Node_1 - 4 + 1 + 1,
			C3_2, C4, C2_2, C1_2, C3_2, C4, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) % 2 == 0)// 1dx
			{
				inner_5node(j, i, 1, -(Num3 - 5),
					(Num2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) + Node11 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) / 2), Coef, Coef_location);
			}
			else//1 lei  2j left
			{
				kind1_node(j, i, -(Num3 - 5), -1, 0, 1,
					(Num2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) + Node11 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) / 2),
					(Num2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) + Node11 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) / 2) + 1,
					(Num2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) + Node11 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) / 2) + (Num1 - N_d0 + 2),
					(Num2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) + Node11 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) / 2) + 1 + (Num1 - N_d0 + 2),
					A4_2, A3_2, A5_2, A3_2, A1_2, A1_2, A2_2, A2_2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 * 2 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Num2 / 2);
		i<(Nx_2 + Nx_3 + Num1 + Num2 * 2 + Num3 * 2 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9); i++)// num 16(2) :  26 node! part2
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2)
		{
		case 0:// 3 lei 2j shang
			kind3_node(j, i, -(6 + Num2 - 1) - 2, -(6 + Num2 - 1), 0, 1, Num2 / 2 + Node13 - 1, Num2 / 2 + Node13,
				C4, C3_2, C1_2, C2_2, C4, C3_2, Coef, Coef_location);
			break;
		case (Num2 / 2 - 1) :// 3 lei 2j xia
			kind3_node(j, i, -(6 + Num2 - 1), -(6 + Num2 - 1) + 2, -1, 0, Node_11 - 4 + 1, Node_11 - 4 + 1 + 1,
			C3_2, C4, C2_2, C1_2, C3_2, C4, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) % 2 == 0)// 1dx
			{
				inner_5node(j, i, 1, -(6 + Num2 - 1), Num2 / 2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) + Node13 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) / 2,
					Coef, Coef_location);
			}
			else//1 lei 2j left
			{
				kind1_node(j, i, -(6 + Num2 - 1), -1, 0, 1,
					(Num2 / 2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) + Node13 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) / 2),
					(Num2 / 2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) + Node13 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) / 2) + 1,
					(Num2 / 2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) + Node13 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) / 2) + (Num1 - Node_11 + Nx_3 - 1),
					(Num2 / 2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) + Node13 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) / 2) + 1 + (Num1 - Node_11 + Nx_3 - 1),
					A4_2, A3_2, A5_2, A3_2, A1_2, A1_2, A2_2, A2_2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (num_total - Nx_2 - Nx_3 - Num1); i<(num_total - Nx_2 - Nx_3 - Num1 + Node_1 + 2); i++)// circle num17(1) Nx_2:  74 node!  part1 !
	{
		switch (i - num_total + Nx_2 + Nx_3 + Num1)
		{
		case 0:
			boundary_node(j, i, 1, -Num5 - Num4 - Num3 - Num2, Nx_2 + Nx_3, Coef, Coef_location);
			break;
		case (N_d0 - 3) ://2 lei shang 1j
			kind2_node(j, i, -Num5 - Num4 - Num3 - Num2, -2, -1, 0, 1, Nx_2 + Nx_3,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_1 + 1) ://2 lei xia 1j
			kind2_node(j, i, -(Num5 - Node_5 + Num4 + Num3 + Num2 + Node_1), -1, 0, 1, 2, Nx_2 + Nx_3,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - num_total + Nx_2 + Nx_3 + Num1) < (N_d0 - 3))// 4dx!
			{
				inner_5node(j, i, 4, -Num5 - Num4 - Num3 - Num2, Nx_2 + Nx_3, Coef, Coef_location);
			}
			else if (((i - num_total + Nx_2 + Nx_3 + Num1) > (N_d0 - 3)) && ((i - num_total + Nx_2 + Nx_3 + Num1) < Node11))// 2dx!
			{
				inner_5node(j, i, 2, -(Num3 + Num2 + N_d0 - 2), Num1 - N_d0 + 2, Coef, Coef_location);
			}
			else if (((i - num_total + Nx_2 + Nx_3 + Num1) >= Node11) && ((i - num_total + Nx_2 + Nx_3 + Num1) <= Node12))// 2 lei right 2j !
			{
				kind2_node(j, i, -((i - num_total + Nx_2 + Nx_3 + Num1) + Num2 - (i - num_total + Nx_2 + Nx_3 + Num1 - Node11) * 2),
					-1, 0, 1, Num1 - N_d0 + 2, Nx_2 + Nx_3, B2, B3, B1_2, B3, B4, B5, Coef, Coef_location);
			}
			else// 2dx!
			{
				inner_5node(j, i, 2, -(Num3 / 2 + 1 + Num2 + Node_1), Num1 - N_d0 + 2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (num_total - Nx_2 - Nx_3 - Num1 + Node_1 + 2); i<(num_total - Nx_2 - Nx_3); i++)// circle num17(1) Nx_2:  74 node!  part2 !
	{
		switch (i - num_total + Nx_2 + Nx_3 + Num1)
		{
		case (Num1 - 1) :
			boundary_node(j, i, 2, -Num1 - Num2 - Num3 - Num4, Nx_2 + Nx_3, Coef, Coef_location);
			break;
		case (Node_1 + N_a - 4) ://2 lei shang 1j
			kind2_node(j, i, -(Num5 - Node_5 + Num4 + Num3 + Num2 + Node_1), -2, -1, 0, 1, Nx_2 + Nx_3,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_11 + 1) ://2 lei xia 1j
			kind2_node(j, i, -Num1 - Num2 - Num3 - Num4, -1, 0, 1, 2, Nx_2 + Nx_3,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - num_total + Nx_2 + Nx_3 + Num1) < (Node_1 + N_a - 4))// 4dx!
			{
				inner_5node(j, i, 4, -(Num5 - Node_5 + Num4 + Num3 + Num2 + Node_1), Nx_2 + Nx_3, Coef, Coef_location);
			}
			else if (((i - num_total + Nx_2 + Nx_3 + Num1) > (Node_1 + N_a - 4)) && ((i - num_total + Nx_2 + Nx_3 + Num1) < Node13))// 2dx!
			{
				inner_5node(j, i, 2, -(Num3 / 2 + Num2 + Node_1 + N_a - 3), Num1 - Node_11 + Nx_3 - 1, Coef, Coef_location);
			}
			else if (((i - num_total + Nx_2 + Nx_3 + Num1) >= Node13) && ((i - num_total + Nx_2 + Nx_3 + Num1) <= Node14))// 2 lei left 2j !
			{
				kind2_node(j, i, -((i - num_total + Nx_2 + Nx_3 + Num1) + Num2 / 2 - (i - num_total + Nx_2 + Nx_3 + Num1 - Node13) * 2),
					-1, 0, 1, Num1 - Node_11 + Nx_3 - 1, Nx_2 + Nx_3, B2, B3, B1_2, B3, B4, B5, Coef, Coef_location);
			}
			else if (((i - num_total + Nx_2 + Nx_3 + Num1) >= Node14) && ((i - num_total + Nx_2 + Nx_3 + Num1) <= Node_11))// 2dx
			{
				inner_5node(j, i, 2, -(Num2 + Node_11 + 1), Num1 - Node_11 + Nx_3 - 1, Coef, Coef_location);
			}
			else
			{
				inner_5node(j, i, 4, -Num1 - Num2 - Num3 - Num4, Nx_2 + Nx_3, Coef, Coef_location);// 4dx!
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (num_total - Nx_2 - Nx_3); i<(num_total - Nx_2); i++)//right1 Nx_3:  30 node! 
	{
		switch (i - num_total + Nx_2 + Nx_3)
		{
		case 0:// 1 lei shang
			kind1_node(j, i, -(Nx_2 - N_d0 + 2) - 2, -(Nx_2 - N_d0 + 2) - 1, -(Nx_2 - N_d0 + 2), 0, 1, Nx_3 + N_d0 - 2 - 2, Nx_3 + N_d0 - 2 - 1, Nx_3 + N_d0 - 2,
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case (Nx_3 / 2 - 1) :// 1 lei xia
			kind1_node(j, i, -(Nx_2 - N_d0 + 2), -(Nx_2 - N_d0 + 2) + 1, -(Nx_2 - N_d0 + 2) + 2, -1, 0,
			Nx_3 + N_d0 - 2, Nx_3 + N_d0 - 2 + 1, Nx_3 + N_d0 - 2 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		case (Nx_3 / 2) :// 1 lei shang
			kind1_node(j, i, -(Nx_2 - Node_11 - 1 + Nx_3) - 2, -(Nx_2 - Node_11 - 1 + Nx_3) - 1, -(Nx_2 - Node_11 - 1 + Nx_3), 0, 1,
			1 + Node_11 - 2, 1 + Node_11 - 1, 1 + Node_11,
			A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case (Nx_3 - 1) :// 1 lei xia
			kind1_node(j, i, -(Nx_2 - Node_11 - 1 + Nx_3), -(Nx_2 - Node_11 - 1 + Nx_3) + 1, -(Nx_2 - Node_11 - 1 + Nx_3) + 2, -1, 0,
			1 + Node_11, 1 + Node_11 + 1, 1 + Node_11 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - num_total + Nx_2 + Nx_3) < (Nx_3 / 2 - 1))
			{
				inner_5node(j, i, 2, -(Nx_2 - N_d0 + 2), Nx_3 + N_d0 - 2, Coef, Coef_location);
			}
			else
			{
				inner_5node(j, i, 2, -(Nx_2 - Node_11 - 1 + Nx_3), 1 + Node_11, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//cout << "offset1 start:" << j << endl;
	

}
void scan_coef_Matrix_middle(complex** Coef, int num_total, int ** Coef_location, int sizeM, int sizeN, int offset)
{
	int i=0,j=0,m=0;
	for( j=0; j<num_total; j++ )         //系数矩阵初始置零
	{
		for( i=0; i<N_matrix; i++ )
		{
			Coef[j][i].real=0;
			Coef[j][i].image=0;
			Coef_location[j][i]=-1;//初始化为负值
		}
	}
	j=0;
    //系数矩阵扫描
	//******************************************************************************************
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9);
		i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Node_5 + 2); i++)// circle num13(5) :  80 node!  part1 !
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9)
		{
		case 0:
			boundary_node(j, j - offset, 1, -Num9 - Num8 - Num7 - Num6, Num5 + Num4 + Num3 + Num2, Coef, Coef_location);
			break;
		case (N_d0 - 3) ://2 lei shang 1j
			kind2_node(j, j - offset, -Num9 - Num8 - Num7 - Num6, -2, -1, 0, 1, Num5 + Num4 + Num3 + Num2,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (N_d0 + 1) :// 2 lei 2j shang
			kind2_node(j, j - offset, -(Num7 + Num6 + N_d0 - 2), -2, -1, 0, 1, Num5 - N_d0 + 2 + Num4,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node51:// inner 3 node
			inner_3node(j, j - offset, -1, 0, (Num5 - N_d0 - 2), 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node51 + 1) :// inner 3 node
			inner_3node(j, j - offset, 0, 1, (Num5 - N_d0 - 2) + N_c1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Node_5 - 3) ://2 lei xia 2j
			kind2_node(j, j - offset, -(Num7 / 2 + 1 + Num6 + Node_5), -1, 0, 1, 2, Num5 - Node_5 + Num4 + Num3 / 2 - 1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_5 + 1) ://2 lei xia 1j
			kind2_node(j, j - offset, -(Num9 - Node_9 + Num8 + Num7 + Num6 + Node_5), -1, 0, 1, 2, Num5 - Node_5 + Num4 + Num3 + Num2 + Node_1,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < (N_d0 - 3))// 4dx!
			{
				inner_5node(j, j - offset, 4, -Num9 - Num8 - Num7 - Num6, Num5 + Num4 + Num3 + Num2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) > (N_d0 - 3)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < (N_d0 + 1)))// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Num7 + Num6 + N_d0 - 2), Num5 - N_d0 + 2 + Num4, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) > (N_d0 + 1)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < Node51))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num6 + N_d0 + 2), (Num5 - N_d0 - 2), Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) > (Node51 + 1)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < (Node_5 - 3)))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num6 + N_d0 + 2) - 2, (Num5 - N_d0 - 2) + N_c1, Coef, Coef_location);
			}
			else// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Num7 / 2 + 1 + Num6 + Node_5), Num5 - Node_5 + Num4 + Num3 / 2 - 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Node_5 + 2);
		i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9); i++)// circle num13(5) :  80 node!  part2 !
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9)
		{
		case (Num5 - 1) :
			boundary_node(j, j - offset, 2, -Num8 - Num7 - Num6 - Num5, Num4 + Num3 + Num2 + Num1, Coef, Coef_location);
			break;
		case (Node_5 + N_a - 4) ://2 lei shang 1j
			kind2_node(j, j - offset, -(Num9 - Node_9 + Num8 + Num7 + Num6 + Node_5), -2, -1, 0, 1, Num5 - Node_5 + Num4 + Num3 + Num2 + Node_1,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_5 + N_a) :// 2 lei 2j shang
			kind2_node(j, j - offset, -(Num7 / 2 + Num6 + Node_5 + N_a - 3), -2, -1, 0, 1, Num5 - Node_5 - N_a + 3 + Num4 + Num3 / 2,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node52:// inner 3 node
			inner_3node(j, j - offset, -1, 0, (Num5 - Node_55 + 4 + Num4 - 1) - N_c1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node52 + 1) :// inner 3 node
			inner_3node(j, j - offset, 0, 1, (Num5 - Node_55 + 4 + Num4 - 1), k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Node_55 - 3) ://2 lei xia 2j
			kind2_node(j, j - offset, -(Num6 + Node_55 + 1), -1, 0, 1, 2, Num5 - Node_55 + Num4 + Num3 - 1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_55 + 1) ://2 lei xia 1j
			kind2_node(j, j - offset, -Num8 - Num7 - Num6 - Num5, -1, 0, 1, 2, Num4 + Num3 + Num2 + Num1,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < (Node_5 + N_a - 4))// 4dx!
			{
				inner_5node(j, j - offset, 4, -(Num9 - Node_9 + Num8 + Num7 + Num6 + Node_5), Num5 - Node_5 + Num4 + Num3 + Num2 + Node_1, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) > (Node_5 + N_a - 4)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < (Node_5 + N_a)))// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Num7 / 2 + Num6 + Node_5 + N_a - 3), Num5 - Node_5 - N_a + 3 + Num4 + Num3 / 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) > (Node_5 + N_a)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < Node52))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Node_55 - 4 + 1) + 2, (Num5 - Node_55 + 4 + Num4 - 1) - N_c1, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) > (Node52 + 1)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < (Node_55 - 3)))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Node_55 - 4 + 1), (Num5 - Node_55 + 4 + Num4 - 1), Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) > (Node_55 - 3)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < (Node_55 + 1)))// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num6 + Node_55 + 1), Num5 - Node_55 + Num4 + Num3 - 1, Coef, Coef_location);
			}
			else// 4dx!
			{
				inner_5node(j, j - offset, 4, -Num8 - Num7 - Num6 - Num5, Num4 + Num3 + Num2 + Num1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9);
		i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Num4 / 2); i++)// num 14(4) :  30 node! part1
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9)
		{
		case 0:// 1 lei 2j shang
			kind1_node(j, j - offset, -(Num5 - N_d0 - 2) - 2, -(Num5 - N_d0 - 2) - 1, -(Num5 - N_d0 - 2), 0, 1, Num4 + 4 - 2, Num4 + 4 - 1, Num4 + 4,
				A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case (Num4 / 2 - 1) :// 1 lei 2j xia
			kind1_node(j, j - offset, -(Num5 - N_d0 - 2) - N_c1, -(Num5 - N_d0 - 2) - N_c1 + 1, -(Num5 - N_d0 - 2) - N_c1 + 2, -1, 0,
			Num4 + 4, Num4 + 4 + 1, Num4 + 4 + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) <= Node41)// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num5 - N_d0 - 2), Num4 + 4, Coef, Coef_location);
			}
			else if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) >= Node42) // 1dx
			{
				inner_5node(j, j - offset, 1, -(Num5 - N_d0 - 2) - N_c1, Num4 + 4, Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_4node(j, j - offset, -1, 0, 1, Num4 + 4, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Num4 / 2);
		i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9); i++)// num14(4) :  30 node! part2
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9)
		{
		case (Num4 / 2) :// 1 lei 2j shang
			kind1_node(j, j - offset, -(Num5 - Node_55 + 4 + Num4 - 1) + N_c1 - 2, -(Num5 - Node_55 + 4 + Num4 - 1) + N_c1 - 1, -(Num5 - Node_55 + 4 + Num4 - 1) + N_c1, 0, 1,
			(Num3 - 5 + 1) - 2, (Num3 - 5 + 1) - 1, (Num3 - 5 + 1),
			A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case (Num4 - 1) :// 1 lei 2j xia
			kind1_node(j, j - offset, -(Num5 - Node_55 + 4 + Num4 - 1), -(Num5 - Node_55 + 4 + Num4 - 1) + 1, -(Num5 - Node_55 + 4 + Num4 - 1) + 2, -1, 0,
			(Num3 - 5 + 1), (Num3 - 5 + 1) + 1, (Num3 - 5 + 1) + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) <= Node43)// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num5 - Node_55 + 4 + Num4 - 1) + N_c1, (Num3 - 5 + 1), Coef, Coef_location);
			}
			else if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) >= Node44) // 1dx
			{
				inner_5node(j, j - offset, 1, -(Num5 - Node_55 + 4 + Num4 - 1), (Num3 - 5 + 1), Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_4node(j, j - offset, -1, 0, 1, (Num3 - 5 + 1), 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9);
		i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Num3 / 2); i++)// num 15(3) :  46 node! part1
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9)
		{
		case 0:// 1 lei 1j shang
			kind1_node(j, j - offset, -(Num5 - N_d0 + 2 + Num4) - 2, -(Num5 - N_d0 + 2 + Num4) - 1, -(Num5 - N_d0 + 2 + Num4), 0, 1,
				Num3 + Num2 + N_d0 - 2 - 2, Num3 + Num2 + N_d0 - 2 - 1, Num3 + Num2 + N_d0 - 2,
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case 3:// 2 lei 2j shang
			kind2_node(j, j - offset, -(Num5 - N_d0 + 2 + Num4), -2, -1, 0, 1, Num3 + Num2 + N_d0 - 2,
				B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case 4:// 3 lei 2j shang
			kind3_node(j, j - offset, -(Num4 + 4), -1, 0, 1, (Num3 + Num2 + N_d0 - 2) - 1, (Num3 + Num2 + N_d0 - 2),
				C2_2, C3_2, C1_2, C3_2, C4, C4, Coef, Coef_location);
			break;
		case (Num3 / 2 - 5) :// 3 lei 2j xia
			kind3_node(j, j - offset, -(Num4 + 4), -1, 0, 1, Num3 / 2 + 1 + Num2 + Node_1, Num3 / 2 + 1 + Num2 + Node_1 + 1,
			C2_2, C3_2, C1_2, C3_2, C4, C4, Coef, Coef_location);
			break;
		case (Num3 / 2 - 4) :// 2 lei 2j xia
			kind2_node(j, j - offset, -(Num5 - Node_5 + Num4 + Num3 / 2 - 1), -1, 0, 1, 2, Num3 / 2 + 1 + Num2 + Node_1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num3 / 2 - 1) :// 1 lei 1j xia
			kind1_node(j, j - offset, -(Num5 - Node_5 + Num4 + Num3 / 2 - 1), -(Num5 - Node_5 + Num4 + Num3 / 2 - 1) + 1, -(Num5 - Node_5 + Num4 + Num3 / 2 - 1) + 2, -1, 0,
			Num3 / 2 + 1 + Num2 + Node_1, Num3 / 2 + 1 + Num2 + Node_1 + 1, Num3 / 2 + 1 + Num2 + Node_1 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < 3)// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num5 - N_d0 + 2 + Num4), Num3 + Num2 + N_d0 - 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) > 4) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < (Num3 / 2 - 5)))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num4 + 4), Num3 - 5, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num5 - Node_5 + Num4 + Num3 / 2 - 1), Num3 / 2 + 1 + Num2 + Node_1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Num3 / 2);
		i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 * 2 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9); i++)// num 15(3) :  46 node! part2
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num3 / 2)
		{
		case 0:// 1 lei 1j shang
			kind1_node(j, j - offset, -(Num5 - Node_5 - N_a + 3 + Num4 + Num3 / 2) - 2, -(Num5 - Node_5 - N_a + 3 + Num4 + Num3 / 2) - 1, -(Num5 - Node_5 - N_a + 3 + Num4 + Num3 / 2), 0, 1,
				Num3 / 2 + Num2 + Node_1 + N_a - 3 - 2, Num3 / 2 + Num2 + Node_1 + N_a - 3 - 1, Num3 / 2 + Num2 + Node_1 + N_a - 3,
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case 3:// 2 lei 2j shang
			kind2_node(j, j - offset, -(Num5 - Node_5 - N_a + 3 + Num4 + Num3 / 2), -2, -1, 0, 1, Num3 / 2 + Num2 + Node_1 + N_a - 3,
				B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case 4:// 3 lei 2j shang
			kind3_node(j, j - offset, -(Num3 - 5 + 1), -1, 0, 1, (Num3 / 2 + Num2 + Node_1 + N_a - 3) - 1, (Num3 / 2 + Num2 + Node_1 + N_a - 3),
				C2_2, C3_2, C1_2, C3_2, C4, C4, Coef, Coef_location);
			break;
		case (Num3 / 2 - 5) :// 3 lei 2j xia
			kind3_node(j, j - offset, -(Num3 - 5 + 1), -1, 0, 1, 1 + Num2 + Node_11, 1 + Num2 + Node_11 + 1,
			C2_2, C3_2, C1_2, C3_2, C4, C4, Coef, Coef_location);
			break;
		case (Num3 / 2 - 4) :// 2 lei 2j xia
			kind2_node(j, j - offset, -(Num5 - Node_55 + Num4 + Num3 - 1), -1, 0, 1, 2, 1 + Num2 + Node_11,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num3 / 2 - 1) :// 1 lei 1j xia
			kind1_node(j, j - offset, -(Num5 - Node_55 + Num4 + Num3 - 1), -(Num5 - Node_55 + Num4 + Num3 - 1) + 1, -(Num5 - Node_55 + Num4 + Num3 - 1) + 2, -1, 0,
			1 + Num2 + Node_11, 1 + Num2 + Node_11 + 1, 1 + Num2 + Node_11 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num3 / 2) < 3)// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num5 - Node_5 - N_a + 3 + Num4 + Num3 / 2), Num3 / 2 + Num2 + Node_1 + N_a - 3, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num3 / 2) > 4) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num3 / 2) < (Num3 / 2 - 5)))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num3 - 5 + 1), 6 + Num2 - 1, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num5 - Node_55 + Num4 + Num3 - 1), 1 + Num2 + Node_11, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 * 2 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9);
		i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 * 2 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Num2 / 2); i++)// num 16(2) :  26 node! part1
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9)
		{
		case 0:// 3 lei 2j shang
			kind3_node(j, j - offset, -(Num3 - 3), -(Num3 - 3) + 2, 0, 1, Num2 + Node11 - 1, Num2 + Node11,
				C4, C3_2, C1_2, C2_2, C4, C3_2, Coef, Coef_location);
			break;
		case (Num2 / 2 - 1) ://node2: 3 lei 2j xia
			kind3_node(j, j - offset, -(Num3 - 5), -(Num3 - 5) + 2, -1, 0, Num2 / 2 + Node_1 - 4 + 1, Num2 / 2 + Node_1 - 4 + 1 + 1,
			C3_2, C4, C2_2, C1_2, C3_2, C4, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) % 2 == 0)// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num3 - 5),
					(Num2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) + Node11 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) / 2), Coef, Coef_location);
			}
			else//1 lei  2j left
			{
				kind1_node(j, j - offset, -(Num3 - 5), -1, 0, 1,
					(Num2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) + Node11 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) / 2),
					(Num2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) + Node11 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) / 2) + 1,
					(Num2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) + Node11 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) / 2) + (Num1 - N_d0 + 2),
					(Num2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) + Node11 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) / 2) + 1 + (Num1 - N_d0 + 2),
					A4_2, A3_2, A5_2, A3_2, A1_2, A1_2, A2_2, A2_2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 * 2 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Num2 / 2);
		i<(Nx_2 + Nx_3 + Num1 + Num2 * 2 + Num3 * 2 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9); i++)// num 16(2) :  26 node! part2
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2)
		{
		case 0:// 3 lei 2j shang
			kind3_node(j, j - offset, -(6 + Num2 - 1) - 2, -(6 + Num2 - 1), 0, 1, Num2 / 2 + Node13 - 1, Num2 / 2 + Node13,
				C4, C3_2, C1_2, C2_2, C4, C3_2, Coef, Coef_location);
			break;
		case (Num2 / 2 - 1) :// 3 lei 2j xia
			kind3_node(j, j - offset, -(6 + Num2 - 1), -(6 + Num2 - 1) + 2, -1, 0, Node_11 - 4 + 1, Node_11 - 4 + 1 + 1,
			C3_2, C4, C2_2, C1_2, C3_2, C4, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) % 2 == 0)// 1dx
			{
				inner_5node(j, j - offset, 1, -(6 + Num2 - 1), Num2 / 2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) + Node13 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) / 2,
					Coef, Coef_location);
			}
			else//1 lei 2j left
			{
				kind1_node(j, j - offset, -(6 + Num2 - 1), -1, 0, 1,
					(Num2 / 2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) + Node13 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) / 2),
					(Num2 / 2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) + Node13 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) / 2) + 1,
					(Num2 / 2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) + Node13 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) / 2) + (Num1 - Node_11 + Nx_3 - 1),
					(Num2 / 2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) + Node13 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) / 2) + 1 + (Num1 - Node_11 + Nx_3 - 1),
					A4_2, A3_2, A5_2, A3_2, A1_2, A1_2, A2_2, A2_2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (num_total - Nx_2 - Nx_3 - Num1); i<(num_total - Nx_2 - Nx_3 - Num1 + Node_1 + 2); i++)// circle num17(1) Nx_2:  74 node!  part1 !
	{
		switch (i - num_total + Nx_2 + Nx_3 + Num1)
		{
		case 0:
			boundary_node(j, j - offset, 1, -Num5 - Num4 - Num3 - Num2, Nx_2 + Nx_3, Coef, Coef_location);
			break;
		case (N_d0 - 3) ://2 lei shang 1j
			kind2_node(j, j - offset, -Num5 - Num4 - Num3 - Num2, -2, -1, 0, 1, Nx_2 + Nx_3,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_1 + 1) ://2 lei xia 1j
			kind2_node(j, j - offset, -(Num5 - Node_5 + Num4 + Num3 + Num2 + Node_1), -1, 0, 1, 2, Nx_2 + Nx_3,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - num_total + Nx_2 + Nx_3 + Num1) < (N_d0 - 3))// 4dx!
			{
				inner_5node(j, j - offset, 4, -Num5 - Num4 - Num3 - Num2, Nx_2 + Nx_3, Coef, Coef_location);
			}
			else if (((i - num_total + Nx_2 + Nx_3 + Num1) > (N_d0 - 3)) && ((i - num_total + Nx_2 + Nx_3 + Num1) < Node11))// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Num3 + Num2 + N_d0 - 2), Num1 - N_d0 + 2, Coef, Coef_location);
			}
			else if (((i - num_total + Nx_2 + Nx_3 + Num1) >= Node11) && ((i - num_total + Nx_2 + Nx_3 + Num1) <= Node12))// 2 lei right 2j !
			{
				kind2_node(j, j - offset, -((i - num_total + Nx_2 + Nx_3 + Num1) + Num2 - (i - num_total + Nx_2 + Nx_3 + Num1 - Node11) * 2),
					-1, 0, 1, Num1 - N_d0 + 2, Nx_2 + Nx_3, B2, B3, B1_2, B3, B4, B5, Coef, Coef_location);
			}
			else// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Num3 / 2 + 1 + Num2 + Node_1), Num1 - N_d0 + 2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (num_total - Nx_2 - Nx_3 - Num1 + Node_1 + 2); i<(num_total - Nx_2 - Nx_3); i++)// circle num17(1) Nx_2:  74 node!  part2 !
	{
		switch (i - num_total + Nx_2 + Nx_3 + Num1)
		{
		case (Num1 - 1) :
			boundary_node(j, j - offset, 2, -Num1 - Num2 - Num3 - Num4, Nx_2 + Nx_3, Coef, Coef_location);
			break;
		case (Node_1 + N_a - 4) ://2 lei shang 1j
			kind2_node(j, j - offset, -(Num5 - Node_5 + Num4 + Num3 + Num2 + Node_1), -2, -1, 0, 1, Nx_2 + Nx_3,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_11 + 1) ://2 lei xia 1j
			kind2_node(j, j - offset, -Num1 - Num2 - Num3 - Num4, -1, 0, 1, 2, Nx_2 + Nx_3,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - num_total + Nx_2 + Nx_3 + Num1) < (Node_1 + N_a - 4))// 4dx!
			{
				inner_5node(j, j - offset, 4, -(Num5 - Node_5 + Num4 + Num3 + Num2 + Node_1), Nx_2 + Nx_3, Coef, Coef_location);
			}
			else if (((i - num_total + Nx_2 + Nx_3 + Num1) > (Node_1 + N_a - 4)) && ((i - num_total + Nx_2 + Nx_3 + Num1) < Node13))// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Num3 / 2 + Num2 + Node_1 + N_a - 3), Num1 - Node_11 + Nx_3 - 1, Coef, Coef_location);
			}
			else if (((i - num_total + Nx_2 + Nx_3 + Num1) >= Node13) && ((i - num_total + Nx_2 + Nx_3 + Num1) <= Node14))// 2 lei left 2j !
			{
				kind2_node(j, j - offset, -((i - num_total + Nx_2 + Nx_3 + Num1) + Num2 / 2 - (i - num_total + Nx_2 + Nx_3 + Num1 - Node13) * 2),
					-1, 0, 1, Num1 - Node_11 + Nx_3 - 1, Nx_2 + Nx_3, B2, B3, B1_2, B3, B4, B5, Coef, Coef_location);
			}
			else if (((i - num_total + Nx_2 + Nx_3 + Num1) >= Node14) && ((i - num_total + Nx_2 + Nx_3 + Num1) <= Node_11))// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num2 + Node_11 + 1), Num1 - Node_11 + Nx_3 - 1, Coef, Coef_location);
			}
			else
			{
				inner_5node(j, j - offset, 4, -Num1 - Num2 - Num3 - Num4, Nx_2 + Nx_3, Coef, Coef_location);// 4dx!
			}
			break;
		}
		j = j + 1;
	}
	//cout << "offset:" << j << endl;
	//====================================================================================================
	//调用封装函数
	//left 1 --- right 1
	scan_coef_Matrix_middle_wrapper(j, Coef, Coef_location, num_total, sizeN - 2 * Num1, offset);
	//=======================================================================================================================
	//***************************************************************************************
	//cout << "circle num1:" << j << endl;
	for (i = (Nx_2 + Nx_3); i<(Nx_2 + Nx_3 + Node_1 + 2); i++)// circle num1 Nx_2:  74 node!  part1 !
	{
		switch (i - Nx_2 - Nx_3)
		{
		case 0:
			boundary_node(j, j - offset, 1, -Nx_2 - Nx_3, Num1 + Num2 + Num3 + Num4, Coef, Coef_location);
			break;
		case (N_d0 - 3) ://2 lei shang 1j
			kind2_node(j, j - offset, -Nx_2 - Nx_3, -2, -1, 0, 1, Num1 + Num2 + Num3 + Num4,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_1 + 1) ://2 lei xia 1j
			kind2_node(j, j - offset, -(Nx_2 + Nx_3), -1, 0, 1, 2, Num1 - Node_1 + Num2 + Num3 + Num4 + Node_5,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3) < (N_d0 - 3))// 4dx!
			{
				inner_5node(j, j - offset, 4, -Nx_3 - Nx_2, Num1 + Num2 + Num3 + Num4, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3) > (N_d0 - 3)) && ((i - Nx_2 - Nx_3) < Node11))// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Nx_3 + N_d0 - 2), Num1 - N_d0 + 2 + Num2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3) >= Node11) && ((i - Nx_2 - Nx_3) <= Node12))// 2 lei left 2j !
			{
				kind2_node(j, j - offset, -Nx_3 - Nx_2, -(Nx_3 + N_d0 - 2), -1, 0, 1, Num1 - (i - Nx_2 - Nx_3) + ((i - Nx_2 - Nx_3) - Node11) * 2,
					B5, B4, B3, B1_2, B3, B2, Coef, Coef_location);
			}
			else// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Nx_3 + N_d0 - 2), Num1 - Node_1 + Num2 + Num3 / 2 - 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Node_1 + 2); i<(Nx_2 + Nx_3 + Num1); i++)// circle num1 Nx_2:  74 node!  part2 !
	{
		switch (i - Nx_2 - Nx_3)
		{
		case (Num1 - 1) :
			boundary_node(j, j - offset, 2, -Nx_2 - Nx_3, Num5 + Num2 + Num3 + Num4, Coef, Coef_location);
			break;
		case (Node_1 + N_a - 4) ://2 lei shang 1j
			kind2_node(j, j - offset, -Nx_2 - Nx_3, -2, -1, 0, 1, Num1 - Node_1 + Num2 + Num3 + Num4 + Node_5,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_11 + 1) ://2 lei xia 1j
			kind2_node(j, j - offset, -(Nx_2 + Nx_3), -1, 0, 1, 2, Num5 + Num2 + Num3 + Num4,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3) < (Node_1 + N_a - 4))// 4dx!
			{
				inner_5node(j, j - offset, 4, -Nx_3 - Nx_2, Num1 - Node_1 + Num2 + Num3 + Num4 + Node_5, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3) > (Node_1 + N_a - 4)) && ((i - Nx_2 - Nx_3) < Node13))// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Node_11 + 1), Num1 - Node_1 - N_a + 3 + Num2 + Num3 / 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3) >= Node13) && ((i - Nx_2 - Nx_3) <= Node14))// 2 lei left 2j !
			{
				kind2_node(j, j - offset, -Nx_3 - Nx_2, -(Node_11 + 1), -1, 0, 1, Num1 - (i - Nx_2 - Nx_3) + Num2 / 2 + ((i - Nx_2 - Nx_3) - Node13) * 2,
					B5, B4, B3, B1_2, B3, B2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3) >= Node14) && ((i - Nx_2 - Nx_3) <= Node_11))// 2dx
			{
				inner_5node(j, j - offset, 2, -(Node_11 + 1), Num1 - Node_11 + Num2 + Num3 - 1, Coef, Coef_location);
			}
			else
			{
				inner_5node(j, j - offset, 4, -(Nx_3 + Nx_2), Num5 + Num2 + Num3 + Num4, Coef, Coef_location);// 4dx!
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!!!!!!!!!
	//cout << "circle num 2:" << j << endl;
	for (i = (Nx_2 + Nx_3 + Num1); i<(Nx_2 + Nx_3 + Num1 + Num2 / 2); i++)// num 2 :  26 node! part1
	{
		switch (i - Nx_2 - Nx_3 - Num1)
		{
		case 0:// 3 lei 2j shang
			kind3_node(j, j - offset, -(Num1 - Node11) - 1, -(Num1 - Node11), 0, 1, Num2 + 5 - 2, Num2 + 5,
				C4, C3_2, C1_2, C2_2, C4, C3_2, Coef, Coef_location);
			break;
		case (Num2 / 2 - 1) ://node2: 3 lei 2j xia
			kind3_node(j, j - offset, -(Num1 - Node12 + Num2 / 2 - 1), -(Num1 - Node12 + Num2 / 2 - 1) + 1, -1, 0, Num2 + 5, Num2 + 5 + 2,
			C3_2, C4, C2_2, C1_2, C3_2, C4, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1) % 2 == 0)// 1dx
			{
				inner_5node(j, j - offset, 1, -((i - Nx_2 - Nx_3 - Num1) + Num1 - Node11 - (i - Nx_2 - Nx_3 - Num1) / 2), Num2 + 5, Coef, Coef_location);
			}
			else//1 lei  2j left
			{
				kind1_node(j, j - offset, -((i - Nx_2 - Nx_3 - Num1) + Num1 - Node11 - (i - Nx_2 - Nx_3 - Num1) / 2) - (Nx_3 + N_d0 - 2),
					-((i - Nx_2 - Nx_3 - Num1) + Num1 - Node11 - (i - Nx_2 - Nx_3 - Num1) / 2) - (Nx_3 + N_d0 - 2) + 1,
					-((i - Nx_2 - Nx_3 - Num1) + Num1 - Node11 - (i - Nx_2 - Nx_3 - Num1) / 2),
					-((i - Nx_2 - Nx_3 - Num1) + Num1 - Node11 - (i - Nx_2 - Nx_3 - Num1) / 2) + 1,
					-1, 0, 1, 5 + Num2, A2_2, A2_2, A1_2, A1_2, A3_2, A5_2, A3_2, A4_2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 / 2); i<(Nx_2 + Nx_3 + Num1 + Num2); i++)// num 2 :  26 node! part2
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 / 2)
		{
		case 0:// 3 lei 2j shang
			kind3_node(j, j - offset, -(Num1 - Node13 + Num2 / 2) - 1, -(Num1 - Node13 + Num2 / 2), 0, 1, Num3 - 5 - 2, Num3 - 5,
				C4, C3_2, C1_2, C2_2, C4, C3_2, Coef, Coef_location);
			break;
		case (Num2 / 2 - 1) :// 3 lei 2j xia
			kind3_node(j, j - offset, -(Num1 - Node14 + Num2 - 1), -(Num1 - Node14 + Num2 - 1) + 1, -1, 0, Num3 - 5, Num3 - 5 + 2,
			C3_2, C4, C2_2, C1_2, C3_2, C4, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 / 2) % 2 == 0)// 1dx
			{
				inner_5node(j, j - offset, 1, -((i - Nx_2 - Nx_3 - Num1 - Num2 / 2) + Num1 - Node13 - (i - Nx_2 - Nx_3 - Num1 - Num2 / 2) / 2) - Num2 / 2,
					Num3 - 5, Coef, Coef_location);
			}
			else//1 lei 2j left
			{
				kind1_node(j, j - offset, -((i - Nx_2 - Nx_3 - Num1 - Num2 / 2) + Num1 - Node13 - (i - Nx_2 - Nx_3 - Num1 - Num2 / 2) / 2) - Num2 / 2 - (Node_11 + 1),
					-((i - Nx_2 - Nx_3 - Num1 - Num2 / 2) + Num1 - Node13 - (i - Nx_2 - Nx_3 - Num1 - Num2 / 2) / 2) - Num2 / 2 - (Node_11 + 1) + 1,
					-((i - Nx_2 - Nx_3 - Num1 - Num2 / 2) + Num1 - Node13 - (i - Nx_2 - Nx_3 - Num1 - Num2 / 2) / 2) - Num2 / 2,
					-((i - Nx_2 - Nx_3 - Num1 - Num2 / 2) + Num1 - Node13 - (i - Nx_2 - Nx_3 - Num1 - Num2 / 2) / 2) - Num2 / 2 + 1,
					-1, 0, 1, Num3 - 5, A2_2, A2_2, A1_2, A1_2, A3_2, A5_2, A3_2, A4_2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!!!
	//cout << "circle num 3:" << j << endl;
	for (i = (Nx_2 + Nx_3 + Num1 + Num2); i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 / 2); i++)// num 3 :  46 node! part1
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2)
		{
		case 0:// 1 lei 1j shang
			kind1_node(j, j - offset, -(Num1 - N_d0 + 2 + Num2) - 2, -(Num1 - N_d0 + 2 + Num2) - 1, -(Num1 - N_d0 + 2 + Num2), 0, 1, Num3 + Num4 + N_d0 - 2 - 2,
				Num3 + Num4 + N_d0 - 2 - 1, Num3 + Num4 + N_d0 - 2, A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case 3:// 2 lei 2j shang
			kind2_node(j, j - offset, -(Num1 - N_d0 + 2 + Num2), -2, -1, 0, 1, Num3 + Num4 + N_d0 - 2,
				B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case 4:// 3 lei 2j shang
			kind3_node(j, j - offset, -(Num1 - N_d0 + 2 + Num2) - 1, -(Num1 - N_d0 + 2 + Num2), -1, 0, 1, Num3 - 4,
				C4, C4, C3_2, C1_2, C3_2, C2_2, Coef, Coef_location);
			break;
		case (Num3 / 2 - 5) :// 3 lei 2j xia
			kind3_node(j, j - offset, -(Num1 - Node_1 + Num2 + Num3 / 2 - 1), -(Num1 - Node_1 + Num2 + Num3 / 2 - 1) + 1, -1, 0, 1, Num3 - 4,
			C4, C4, C3_2, C1_2, C3_2, C2_2, Coef, Coef_location);
			break;
		case (Num3 / 2 - 4) :// 2 lei 2j xia
			kind2_node(j, j - offset, -(Num1 - Node_1 + Num2 + Num3 / 2 - 1), -1, 0, 1, 2, Num3 / 2 + 1 + Num4 + Node_5,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num3 / 2 - 1) :// 1 lei 1j xia
			kind1_node(j, j - offset, -(Num1 - Node_1 + Num2 + Num3 / 2 - 1), -(Num1 - Node_1 + Num2 + Num3 / 2 - 1) + 1, -(Num1 - Node_1 + Num2 + Num3 / 2 - 1) + 2, -1, 0,
			Num3 / 2 + 1 + Num4 + Node_5, Num3 / 2 + 1 + Num4 + Node_5 + 1, Num3 / 2 + 1 + Num4 + Node_5 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2) < 3)// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num1 - N_d0 + 2 + Num2), Num3 + Num4 + N_d0 - 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2) > 4) && ((i - Nx_2 - Nx_3 - Num1 - Num2) < (Num3 / 2 - 5)))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num2 + 5), Num3 - 4, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num1 - Node_1 + Num2 + Num3 / 2 - 1), Num3 / 2 + 1 + Num4 + Node_5, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 / 2); i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3); i++)// num 3 :  46 node! part2
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 / 2)
		{
		case 0:// 1 lei 1j shang
			kind1_node(j, j - offset, -(Num1 - Node_1 - N_a + 3 + Num2 + Num3 / 2) - 2, -(Num1 - Node_1 - N_a + 3 + Num2 + Num3 / 2) - 1, -(Num1 - Node_1 - N_a + 3 + Num2 + Num3 / 2), 0, 1,
				Num3 / 2 + Num4 + Node_5 + N_a - 3 - 2, Num3 / 2 + Num4 + Node_5 + N_a - 3 - 1, Num3 / 2 + Num4 + Node_5 + N_a - 3,
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case 3:// 2 lei 2j shang
			kind2_node(j, j - offset, -(Num1 - Node_1 - N_a + 3 + Num2 + Num3 / 2), -2, -1, 0, 1, Num3 / 2 + Num4 + Node_5 + N_a - 3,
				B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case 4:// 3 lei 2j shang
			kind3_node(j, j - offset, -(Num1 - Node_1 - N_a + 3 + Num2 + Num3 / 2) - 1, -(Num1 - Node_1 - N_a + 3 + Num2 + Num3 / 2), -1, 0, 1, 5 + Num4 - 1,
				C4, C4, C3_2, C1_2, C3_2, C2_2, Coef, Coef_location);
			break;
		case (Num3 / 2 - 5) :// 3 lei 2j xia
			kind3_node(j, j - offset, -(Num1 - Node_11 + Num2 + Num3 - 1), -(Num1 - Node_11 + Num2 + Num3 - 1) + 1, -1, 0, 1, 5 + Num4 - 1,
			C4, C4, C3_2, C1_2, C3_2, C2_2, Coef, Coef_location);
			break;
		case (Num3 / 2 - 4) :// 2 lei 2j xia
			kind2_node(j, j - offset, -(Num1 - Node_11 + Num2 + Num3 - 1), -1, 0, 1, 2, 1 + Num4 + Node_55,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num3 / 2 - 1) :// 1 lei 1j xia
			kind1_node(j, j - offset, -(Num1 - Node_11 + Num2 + Num3 - 1), -(Num1 - Node_11 + Num2 + Num3 - 1) + 1, -(Num1 - Node_11 + Num2 + Num3 - 1) + 2, -1, 0,
			1 + Num4 + Node_55, 1 + Num4 + Node_55 + 1, 1 + Num4 + Node_55 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 / 2) < 3)// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num1 - Node_1 - N_a + 3 + Num2 + Num3 / 2), Num3 / 2 + Num4 + Node_5 + N_a - 3, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 / 2) > 4) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 / 2) < (Num3 / 2 - 5)))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num3 - 5), 5 + Num4 - 1, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num1 - Node_11 + Num2 + Num3 - 1), 1 + Num4 + Node_55, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	//cout << "circle num 4:" << j << endl;
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3); i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 / 2); i++)// num 4 :  30 node! part1
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3)
		{
		case 0:// 1 lei 2j shang
			kind1_node(j, j - offset, -(Num3 - 4) - 2, -(Num3 - 4) - 1, -(Num3 - 4), 0, 1, Num4 + N_d0 + 2 - 2, Num4 + N_d0 + 2 - 1, Num4 + N_d0 + 2,
				A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case (Num4 / 2 - 1) :// 1 lei 2j xia
			kind1_node(j, j - offset, -(Num3 - 4), -(Num3 - 4) + 1, -(Num3 - 4) + 2, -1, 0,
			Num4 + N_d0 + 2 - N_c1, Num4 + N_d0 + 2 - N_c1 + 1, Num4 + N_d0 + 2 - N_c1 + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3) <= Node41)// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num3 - 4), Num4 + N_d0 + 2, Coef, Coef_location);
			}
			else if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3) >= Node42) // 1dx
			{
				inner_5node(j, j - offset, 1, -(Num3 - 4), Num4 + N_d0 + 2 - N_c1, Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_4node(j, j - offset, -(Num3 - 4), -1, 0, 1, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 / 2); i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4); i++)// num 4 :  30 node! part2
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3)
		{
		case (Num4 / 2) :// 1 lei 2j shang
			kind1_node(j, j - offset, -(5 + Num4 - 1) - 2, -(5 + Num4 - 1) - 1, -(5 + Num4 - 1), 0, 1,
			(Node_55 - 4 + 1 + N_c1) - 2, (Node_55 - 4 + 1 + N_c1) - 1, (Node_55 - 4 + 1 + N_c1),
			A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case (Num4 - 1) :// 1 lei 2j xia
			kind1_node(j, j - offset, -(5 + Num4 - 1), -(5 + Num4 - 1) + 1, -(5 + Num4 - 1) + 2, -1, 0,
			(Node_55 - 4 + 1), (Node_55 - 4 + 1) + 1, (Node_55 - 4 + 1) + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3) <= Node43)// 1dx
			{
				inner_5node(j, j - offset, 1, -(5 + Num4 - 1), (Node_55 - 4 + 1 + N_c1), Coef, Coef_location);
			}
			else if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3) >= Node44) // 1dx
			{
				inner_5node(j, j - offset, 1, -(5 + Num4 - 1), (Node_55 - 4 + 1), Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_4node(j, j - offset, -(5 + Num4 - 1), -1, 0, 1, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	//cout << "circle num5:" << j << endl;
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4); i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Node_5 + 2); i++)// circle num5 :  80 node!  part1 !
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4)
		{
		case 0:
			boundary_node(j, j - offset, 1, -Num1 - Num2 - Num3 - Num4, Num5 + Num6 + Num7 + Num8, Coef, Coef_location);
			break;
		case (N_d0 - 3) ://2 lei shang 1j
			kind2_node(j, j - offset, -Num1 - Num2 - Num3 - Num4, -2, -1, 0, 1, Num5 + Num6 + Num7 + Num8,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (N_d0 + 1) :// 2 lei 2j shang
			kind2_node(j, j - offset, -(Num3 + Num4 + N_d0 - 2), -2, -1, 0, 1, Num5 - N_d0 + 2 + Num6,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node51:// inner 3 node
			inner_3node(j, j - offset, -(Num4 + N_d0 + 2), -1, 0, 1, 1, k*k*dx*dx - 4, Coef, Coef_location);
			break;
		case (Node51 + 1) :// inner 3 node
			inner_3node(j, j - offset, -(Num4 + N_d0 + 2) + N_c1, 0, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node_5 - 3) ://2 lei xia 2j
			kind2_node(j, j - offset, -(Num3 / 2 + 1 + Num4 + Node_5), -1, 0, 1, 2, Num5 - Node_5 + Num6 + Num7 / 2 - 1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_5 + 1) ://2 lei xia 1j
			kind2_node(j, j - offset, -(Num1 - Node_1 + Num2 + Num3 + Num4 + Node_5), -1, 0, 1, 2, Num5 - Node_5 + Num6 + Num7 + Num8 + Node_9,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < (N_d0 - 3))// 4dx!
			{
				inner_5node(j, j - offset, 4, -Num1 - Num2 - Num3 - Num4, Num5 + Num6 + Num7 + Num8, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) > (N_d0 - 3)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < (N_d0 + 1)))// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Num3 + Num4 + N_d0 - 2), Num5 - N_d0 + 2 + Num6, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) > (N_d0 + 1)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < Node51))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num4 + N_d0 + 2), Num5 - N_d0 - 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) > (Node51 + 1)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < (Node_5 - 3)))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num4 + N_d0 + 2) + N_c1, Num5 - N_d0 - 2 - 2, Coef, Coef_location);
			}
			else// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Num3 / 2 + 1 + Num4 + Node_5), Num5 - Node_5 + Num6 + Num7 / 2 - 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Node_5 + 2); i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5); i++)// circle num5 :  80 node!  part2 !
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4)
		{
		case (Num5 - 1) :
			boundary_node(j, j - offset, 2, -Num2 - Num3 - Num4 - Num5, Num6 + Num7 + Num8 + Num9, Coef, Coef_location);
			break;
		case (Node_5 + N_a - 4) ://2 lei shang 1j
			kind2_node(j, j - offset, -(Num1 - Node_1 + Num2 + Num3 + Num4 + Node_5), -2, -1, 0, 1, Num5 - Node_5 + Num6 + Num7 + Num8 + Node_9,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_5 + N_a) :// 2 lei 2j shang
			kind2_node(j, j - offset, -(Num3 / 2 + Num4 + Node_5 + N_a - 3), -2, -1, 0, 1, Num5 - Node_5 - N_a + 3 + Num6 + Num7 / 2,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node52:// inner 3 node
			inner_3node(j, j - offset, -(Node_55 - 4 + 1) - N_c1, -1, 0, 1, 1, k*k*dx*dx - 4, Coef, Coef_location);
			break;
		case (Node52 + 1) :// inner 3 node
			inner_3node(j, j - offset, -(Node_55 - 4 + 1), 0, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node_55 - 3) ://2 lei xia 2j
			kind2_node(j, j - offset, -(Num4 + Node_55 + 1), -1, 0, 1, 2, Num5 - Node_55 + Num6 + Num7 - 1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_55 + 1) ://2 lei xia 1j
			kind2_node(j, j - offset, -Num2 - Num3 - Num4 - Num5, -1, 0, 1, 2, Num6 + Num7 + Num8 + Num9,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < (Node_5 + N_a - 4))// 4dx!
			{
				inner_5node(j, j - offset, 4, -(Num1 - Node_1 + Num2 + Num3 + Num4 + Node_5), Num5 - Node_5 + Num6 + Num7 + Num8 + Node_9, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) > (Node_5 + N_a - 4)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < (Node_5 + N_a)))// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Num3 / 2 + Num4 + Node_5 + N_a - 3), Num5 - Node_5 - N_a + 3 + Num6 + Num7 / 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) > (Node_5 + N_a)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < Node52))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Node_55 - 4 + 1) - N_c1, (Num5 - Node_55 + 4 + Num6 - 1) + 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) > (Node52 + 1)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < (Node_55 - 3)))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Node_55 - 4 + 1), (Num5 - Node_55 + 4 + Num6 - 1), Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) > (Node_55 - 3)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < (Node_55 + 1)))// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num4 + Node_55 + 1), Num5 - Node_55 + Num6 + Num7 - 1, Coef, Coef_location);
			}
			else// 4dx!
			{
				inner_5node(j, j - offset, 4, -Num2 - Num3 - Num4 - Num5, Num6 + Num7 + Num8 + Num9, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//**************************************************************
	cout << "middle区域节点数:" << j << endl;
	if (j == sizeM)
		cout << "middle passed..." << endl;
	else
		cout << "middle failed..." << endl;
	//修正矩阵，将不在本区域的节点刨除
	for (int i = 0; i != sizeM; i++)
		for (int j = 0; j != N_matrix; j++){
		if (Coef_location[i][j] < 0 || Coef_location[i][j] >= sizeN)
			Coef_location[i][j] = -1;
		}
	if (Debug)
	{
		ofstream outFile_Coef("middle.txt");
		if (!outFile_Coef)
		{
			cout << "Can't write the coefficient matrix into file !" << endl;
			exit(1);//Terminate with fault
		}
		else
		{
			int i, j;
			for (i = 0; i < sizeM; i++)
			{
				for (j = 0; j < N_matrix; j++)
				{
					if ((Coef[i][j].real != 0 || Coef[i][j].image != 0) && Coef_location[i][j] >= 0)
						outFile_Coef << Coef[i][j].real << "+j*" << Coef[i][j].image << " (" << i << "," << Coef_location[i][j] << ") ";
				}
				outFile_Coef << endl << endl;
			}
		}
		outFile_Coef.close();
	}

};
void scan_coef_Matrix_middle_1(complex** Coef, int num_total, int ** Coef_location, int sizeM, int sizeN, int offset)
{
	int i = 0, j = 0, m = 0;
	for (j = 0; j<num_total; j++)         //系数矩阵初始置零
	{
		for (i = 0; i<N_matrix; i++)
		{
			Coef[j][i].real = 0;
			Coef[j][i].image = 0;
			Coef_location[j][i] = -1;//初始化为负值
		}
	}
	j = 0;
	//系数矩阵扫描
	//****************************************************************************************!!!
	for (i = (Nx_side*N_d1); i<(Nx_side*N_d1 + Nx_side); i++)//middle right ,one line
	{
		switch (i - Nx_side*N_d1)
		{
		case 0:
			boundary_node(j, j - offset, 1, -Nx_side, Nx_side + Nx_1, Coef, Coef_location);
			break;
		case (Nx_side - 1) :
			boundary_node(j, j - offset, 2, -Nx_side, Nx_2 + Nx_1, Coef, Coef_location);
			break;
		case (N_d0 - 2) :// 2 lei
			kind2_node(j, j - offset, -Nx_side * 2, -Nx_side, -1, 0, 1, Nx_side - N_d0 + 2,
			B5, B4, B3, B1_1, B3, B2, Coef, Coef_location);
			break;
		case (N_d0 - 1) :// 2 lei缺一点
			other_5_node(j, j - offset, -2 * Nx_side, -Nx_side, -1, 0, (Nx_side - N_d0 + 1 + 2),
			B5, B4, B3, B1_1, B2, Coef, Coef_location);
			break;
		case N_d0:// 2 lei缺一点
			other_5_node(j, j - offset, -2 * Nx_side, -Nx_side, 0, 1, (Nx_side - N_d0 + Nx_1 / 2 - 3),
				B5, B4, B1_1, B3, B2, Coef, Coef_location);
			break;
		case (N_d0 + 1) :// 2 lei
			kind2_node(j, j - offset, -Nx_side * 2, -Nx_side, -1, 0, 1, Nx_side - N_d0 - 1 + Nx_1 / 2 - 1,
			B5, B4, B3, B1_1, B3, B2, Coef, Coef_location);
			break;
		case (N_d0 + N_a - 2) :// 2 lei
			kind2_node(j, j - offset, -Nx_side * 2, -Nx_side, -1, 0, 1, (Nx_side - N_d0 - N_a + 2 + Nx_1 / 2),
			B5, B4, B3, B1_1, B3, B2, Coef, Coef_location);
			break;
		case (N_d0 + N_a - 1) :// 2 lei缺一点
			other_5_node(j, j - offset, -2 * Nx_side, -Nx_side, -1, 0, (Nx_side - N_d0 - N_a + 1 + Nx_1 / 2 + 2),
			B5, B4, B3, B1_1, B2, Coef, Coef_location);
			break;
		case (N_d0 + N_a) :// 2 lei缺一点
			other_5_node(j, j - offset, -2 * Nx_side, -Nx_side, 0, 1, (Nx_side - N_d0 - N_a + Nx_1 - 3),
			B5, B4, B1_1, B3, B2, Coef, Coef_location);
			break;
		case (N_d0 + N_a + 1) :// 2 lei
			kind2_node(j, j - offset, -Nx_side * 2, -Nx_side, -1, 0, 1, (Nx_side - N_d0 - N_a - 1 + Nx_1 - 1),
			B5, B4, B3, B1_1, B3, B2, Coef, Coef_location);
			break;
		default:// 4dx
			if ((i - Nx_side*N_d1) < (N_d0 - 2))//up
			{
				inner_5node(j, j - offset, 4, -Nx_side, Nx_side + Nx_1, Coef, Coef_location);
			}
			else if ((i - Nx_side*N_d1) > (N_d0 + N_a + 1))//down
			{
				inner_5node(j, j - offset, 4, -Nx_side, Nx_2 + Nx_1, Coef, Coef_location);
			}
			else//middle
			{
				inner_5node(j, j - offset, 4, -Nx_side, Nx_side - N_d0 - 1 + Nx_1 + Node_1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//**********************************************************************************************
	for (i = (Nx_side*N_d1 + Nx_side); i<(Nx_side*N_d1 + Nx_side + Nx_1 / 2); i++)//right1 Nx_1: 26 node part1
	{
		switch (i - Nx_side*N_d1 - Nx_side)
		{
		case 0:// 3 lei
			kind3_node(j, j - offset, -(Nx_side - N_d0 + 2) - 1, -(Nx_side - N_d0 + 2), 0, 1, Nx_1 + N_d0 - 3, Nx_1 + N_d0 - 1,
				C4, C3_1, C1_1, C2_1, C4, C3_1, Coef, Coef_location);
			break;
		case (Nx_1 / 2 - 1) ://node2: 3 lei
			kind3_node(j, j - offset, -(Nx_side - N_d0 - 1 + Nx_1 / 2 - 1), -(Nx_side - N_d0 - 1 + Nx_1 / 2 - 1) + 1, -1, 0, Nx_1 + N_d0 - 1, Nx_1 + N_d0 - 1 + 2,
			C3_1, C4, C2_1, C1_1, C3_1, C4, Coef, Coef_location);
			break;
		case 1:// 1 lei
			kind1_node(j, j - offset, -(Nx_side * 2 - N_d0 + 2 + 1), -(Nx_side * 2 - N_d0 + 2 + 1) + 1, -(Nx_side - N_d0 + 2 + 1), -(Nx_side - N_d0 + 2), -1, 0, 1, Nx_1 + N_d0 - 1,
				A2_1, A2_1, A1_1, A1_1, A3_1, A5_1, A3_1, A4_1, Coef, Coef_location);
			break;
		case 2:// 5 node 2dx
			inner_5node(j, j - offset, 2, -(Nx_side - N_d0 + 1 + 2), Nx_1 + N_d0 - 1, Coef, Coef_location);
			break;
		case 3:// 1 lei que 2 dian
			kind2_node(j, j - offset, -(Nx_side - N_d0 + 1 + 3) - Nx_side, -(Nx_side - N_d0 + 1 + 3), -1, 0, 1,
				Nx_1 + N_d0 - 1, A2_1, A1_1, A3_1, A5_1, A3_1, A4_1, Coef, Coef_location);
			break;
		case (Nx_1 / 2 - 4) :// 1 lei que 2 dian
			kind2_node(j, j - offset, -(Nx_side - N_d0 + Nx_1 / 2 - 4) - Nx_side, -(Nx_side - N_d0 + Nx_1 / 2 - 4), -1, 0, 1,
			Nx_1 + N_d0 - 1, A2_1, A1_1, A3_1, A5_1, A3_1, A4_1, Coef, Coef_location);
			break;
		case (Nx_1 / 2 - 3) :// 5 node 2dx
			inner_5node(j, j - offset, 2, -(Nx_side - N_d0 + Nx_1 / 2 - 3), Nx_1 + N_d0 - 1, Coef, Coef_location);
			break;
		case (Nx_1 / 2 - 2) :// 1 lei
			kind1_node(j, j - offset, -(Nx_side * 2 - N_d0 + Nx_1 / 2 - 2), -(Nx_side * 2 - N_d0 + Nx_1 / 2 - 2) + 1, -(Nx_side - N_d0 + Nx_1 / 2 - 2), -(Nx_side - N_d0 + Nx_1 / 2 - 2) + 1,
			-1, 0, 1, Nx_1 + N_d0 - 1,
			A2_1, A2_1, A1_1, A1_1, A3_1, A5_1, A3_1, A4_1, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_side*N_d1 - Nx_side) % 2 == 0)
			{
				inner_4node(j, j - offset, -1, 0, 1, Nx_1 + N_d0 - 1, 1, k*k * 4 * dx*dx - 4, 1, 1, Coef, Coef_location);
			}
			else// 1 lei que 4 dian
			{
				inner_4node(j, j - offset, -1, 0, 1, Nx_1 + N_d0 - 1, A3_1, A5_1, A3_1, A4_1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_side*N_d1 + Nx_side + Nx_1 / 2); i<(Nx_side*N_d1 + Nx_side + Nx_1); i++)// Nx_1 :  26 node! part2
	{
		switch (i - Nx_side*N_d1 - Nx_side - Nx_1 / 2)
		{
		case 0:// 3 lei
			kind3_node(j, j - offset, -(Nx_side - N_d0 - N_a + 2 + Nx_1 / 2) - 1, -(Nx_side - N_d0 - N_a + 2 + Nx_1 / 2), 0, 1, Node_11 - 2, Node_11,
				C4, C3_1, C1_1, C2_1, C4, C3_1, Coef, Coef_location);
			break;
		case (Nx_1 / 2 - 1) :// 3 lei
			kind3_node(j, j - offset, -(Nx_side - N_d0 - N_a + Nx_1 - 2), -(Nx_side - N_d0 - N_a + Nx_1 - 2) + 1, -1, 0, Node_11, Node_11 + 2,
			C3_1, C4, C2_1, C1_1, C3_1, C4, Coef, Coef_location);
			break;
		case 1:// 1 lei
			kind1_node(j, j - offset, -(Nx_side * 2 - N_d0 - N_a + 2 + Nx_1 / 2) - 1, -(Nx_side * 2 - N_d0 - N_a + 2 + Nx_1 / 2), -(Nx_side - N_d0 - N_a + 2 + Nx_1 / 2) - 1, -(Nx_side - N_d0 - N_a + 2 + Nx_1 / 2),
				-1, 0, 1, Node_11,
				A2_1, A2_1, A1_1, A1_1, A3_1, A5_1, A3_1, A4_1, Coef, Coef_location);
			break;
		case 2:// 5 node 2dx
			inner_5node(j, j - offset, 2, -(Nx_side - N_d0 - N_a + 1 + Nx_1 / 2 + 2), Node_11, Coef, Coef_location);
			break;
		case 3:// 1 lei que 2 dian
			kind2_node(j, j - offset, -(Nx_side - N_d0 - N_a + 1 + Nx_1 / 2 + 3) - Nx_side, -(Nx_side - N_d0 - N_a + 1 + Nx_1 / 2 + 3), -1, 0, 1,
				Node_11, A2_1, A1_1, A3_1, A5_1, A3_1, A4_1, Coef, Coef_location);
			break;
		case (Nx_1 / 2 - 4) :// 1 lei que 2 dian
			kind2_node(j, j - offset, -(Nx_side - N_d0 - N_a + Nx_1 - 4) - Nx_side, -(Nx_side - N_d0 - N_a + Nx_1 - 4), -1, 0, 1, Node_11,
			A2_1, A1_1, A3_1, A5_1, A3_1, A4_1, Coef, Coef_location);
			break;
		case (Nx_1 / 2 - 3) :// 5 node 2dx
			inner_5node(j, j - offset, 2, -(Nx_side - N_d0 - N_a + Nx_1 - 3), Node_11, Coef, Coef_location);
			break;
		case (Nx_1 / 2 - 2) :// 1 lei
			kind1_node(j, j - offset, -(Nx_side * 2 - N_d0 - N_a + Nx_1 - 2), -(Nx_side * 2 - N_d0 - N_a + Nx_1 - 2) + 1, -(Nx_side - N_d0 - N_a + Nx_1 - 2), -(Nx_side - N_d0 - N_a + Nx_1 - 2) + 1,
			-1, 0, 1, Node_11, A2_1, A2_1, A1_1, A1_1, A3_1, A5_1, A3_1, A4_1, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_side*N_d1 - Nx_side - Nx_1 / 2) % 2 == 0)
			{
				inner_4node(j, j - offset, -1, 0, 1, Node_11, 1, k*k * 4 * dx*dx - 4, 1, 1, Coef, Coef_location);
			}
			else//1 lei
			{
				inner_4node(j, j - offset, -1, 0, 1, Node_11, A3_1, A5_1, A3_1, A4_1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//***************************************************************************************
	for (i = (Nx_side*N_d1 + Nx_side + Nx_1); i<(Nx_side*N_d1 + Nx_side + Nx_1 + Nx_2); i++)//right 2 Nx_2 :  74 node! 
	{
		switch (i - Nx_side*N_d1 - Nx_1 - Nx_side)
		{
		case 0:
			boundary_node(j, j - offset, 1, -Nx_side - Nx_1, Nx_2 + Nx_3, Coef, Coef_location);
			break;
		case (Nx_2 - 1) :
			boundary_node(j, j - offset, 2, -Nx_1 - Nx_2, Nx_2 + Nx_3, Coef, Coef_location);
			break;
		case (N_d0 - 3) ://2 lei shang
			kind2_node(j, j - offset, -Nx_side - Nx_1, -2, -1, 0, 1, Nx_2 + Nx_3,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (N_d0 - 2) ://3 lei shang
			kind3_node(j, j - offset, -Nx_side - Nx_1 - 1, -Nx_side - Nx_1, -1, 0, 1, Nx_2 - N_d0 + 2,
			C4, C4, C3_1, C1_1, C3_1, C2_1, Coef, Coef_location);
			break;
		case Node_1://3 lei xia
			kind3_node(j, j - offset, -(Nx_side - N_d0 - 1 + Nx_1 + Node_1), -(Nx_side - N_d0 - 1 + Nx_1 + Node_1) + 1, -1, 0, 1, Nx_2 - N_d0 + 2,
				C4, C4, C3_1, C1_1, C3_1, C2_1, Coef, Coef_location);
			break;
		case (Node_1 + 1) ://2 lei xia
			kind2_node(j, j - offset, -(Nx_side - N_d0 - 1 + Nx_1 + Node_1), -1, 0, 1, 2, Nx_2 + Nx_3,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_1 + N_a - 4) ://2 lei shang
			kind2_node(j, j - offset, -(Nx_side - N_d0 - 1 + Nx_1 + Node_1), -2, -1, 0, 1, Nx_2 + Nx_3,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_1 + N_a - 3) ://3 lei shang
			kind3_node(j, j - offset, -(Nx_side - N_d0 - 1 + Nx_1 + Node_1) - 1, -(Nx_side - N_d0 - 1 + Nx_1 + Node_1), -1, 0, 1, Nx_2 - Node_11 + Nx_3 - 1,
			C4, C4, C3_1, C1_1, C3_1, C2_1, Coef, Coef_location);
			break;
		case Node_11://3 lei xia
			kind3_node(j, j - offset, -Nx_1 - Nx_2, -Nx_1 - Nx_2 + 1, -1, 0, 1, Nx_2 - Node_11 + Nx_3 - 1,
				C4, C4, C3_1, C1_1, C3_1, C2_1, Coef, Coef_location);
			break;
		case (Node_11 + 1) ://2 lei xia
			kind2_node(j, j - offset, -Nx_1 - Nx_2, -1, 0, 1, 2, Nx_2 + Nx_3,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_side*N_d1 - Nx_1 - Nx_side) < (N_d0 - 3))// 4dx
			{
				inner_5node(j, j - offset, 4, -Nx_side - Nx_1, Nx_2 + Nx_3, Coef, Coef_location);
			}
			else if (((i - Nx_side*N_d1 - Nx_1 - Nx_side) >= (N_d0 - 1)) && ((i - Nx_side*N_d1 - Nx_1 - Nx_side) < Node_1))// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Nx_1 + N_d0 - 1), Nx_2 - N_d0 + 2, Coef, Coef_location);
			}
			else if (((i - Nx_side*N_d1 - Nx_1 - Nx_side) > (Node_1 + 1)) && ((i - Nx_side*N_d1 - Nx_1 - Nx_side) < (Node_1 + N_a - 4)))// 4dx!
			{
				inner_5node(j, j - offset, 4, -(Nx_side - N_d0 - 1 + Nx_1 + Node_1), Nx_2 + Nx_3, Coef, Coef_location);
			}
			else if (((i - Nx_side*N_d1 - Nx_1 - Nx_side) >= (Node_1 + N_a - 2)) && ((i - Nx_side*N_d1 - Nx_1 - Nx_side) < Node_11))// 2dx!
			{
				inner_5node(j, j - offset, 2, -Node_11, Nx_2 - Node_11 + Nx_3 - 1, Coef, Coef_location);
			}
			else// 4dx!
			{
				inner_5node(j, j - offset, 4, -Nx_1 - Nx_2, Nx_2 + Nx_3, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//====================================================================================================
	//调用封装函数
	//left 1 --- right 1
	scan_coef_Matrix_middle_wrapper(j, Coef, Coef_location, num_total, sizeN - 2 * Num1, offset);
	//=======================================================================================================================
	//***************************************************************************************
	//cout << "circle num1:" << j << endl;
	for (i = (Nx_2 + Nx_3); i<(Nx_2 + Nx_3 + Node_1 + 2); i++)// circle num1 Nx_2:  74 node!  part1 !
	{
		switch (i - Nx_2 - Nx_3)
		{
		case 0:
			boundary_node(j, j - offset, 1, -Nx_2 - Nx_3, Num1 + Num2 + Num3 + Num4, Coef, Coef_location);
			break;
		case (N_d0 - 3) ://2 lei shang 1j
			kind2_node(j, j - offset, -Nx_2 - Nx_3, -2, -1, 0, 1, Num1 + Num2 + Num3 + Num4,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_1 + 1) ://2 lei xia 1j
			kind2_node(j, j - offset, -(Nx_2 + Nx_3), -1, 0, 1, 2, Num1 - Node_1 + Num2 + Num3 + Num4 + Node_5,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3) < (N_d0 - 3))// 4dx!
			{
				inner_5node(j, j - offset, 4, -Nx_3 - Nx_2, Num1 + Num2 + Num3 + Num4, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3) > (N_d0 - 3)) && ((i - Nx_2 - Nx_3) < Node11))// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Nx_3 + N_d0 - 2), Num1 - N_d0 + 2 + Num2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3) >= Node11) && ((i - Nx_2 - Nx_3) <= Node12))// 2 lei left 2j !
			{
				kind2_node(j, j - offset, -Nx_3 - Nx_2, -(Nx_3 + N_d0 - 2), -1, 0, 1, Num1 - (i - Nx_2 - Nx_3) + ((i - Nx_2 - Nx_3) - Node11) * 2,
					B5, B4, B3, B1_2, B3, B2, Coef, Coef_location);
			}
			else// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Nx_3 + N_d0 - 2), Num1 - Node_1 + Num2 + Num3 / 2 - 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Node_1 + 2); i<(Nx_2 + Nx_3 + Num1); i++)// circle num1 Nx_2:  74 node!  part2 !
	{
		switch (i - Nx_2 - Nx_3)
		{
		case (Num1 - 1) :
			boundary_node(j, j - offset, 2, -Nx_2 - Nx_3, Num5 + Num2 + Num3 + Num4, Coef, Coef_location);
			break;
		case (Node_1 + N_a - 4) ://2 lei shang 1j
			kind2_node(j, j - offset, -Nx_2 - Nx_3, -2, -1, 0, 1, Num1 - Node_1 + Num2 + Num3 + Num4 + Node_5,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_11 + 1) ://2 lei xia 1j
			kind2_node(j, j - offset, -(Nx_2 + Nx_3), -1, 0, 1, 2, Num5 + Num2 + Num3 + Num4,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3) < (Node_1 + N_a - 4))// 4dx!
			{
				inner_5node(j, j - offset, 4, -Nx_3 - Nx_2, Num1 - Node_1 + Num2 + Num3 + Num4 + Node_5, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3) > (Node_1 + N_a - 4)) && ((i - Nx_2 - Nx_3) < Node13))// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Node_11 + 1), Num1 - Node_1 - N_a + 3 + Num2 + Num3 / 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3) >= Node13) && ((i - Nx_2 - Nx_3) <= Node14))// 2 lei left 2j !
			{
				kind2_node(j, j - offset, -Nx_3 - Nx_2, -(Node_11 + 1), -1, 0, 1, Num1 - (i - Nx_2 - Nx_3) + Num2 / 2 + ((i - Nx_2 - Nx_3) - Node13) * 2,
					B5, B4, B3, B1_2, B3, B2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3) >= Node14) && ((i - Nx_2 - Nx_3) <= Node_11))// 2dx
			{
				inner_5node(j, j - offset, 2, -(Node_11 + 1), Num1 - Node_11 + Num2 + Num3 - 1, Coef, Coef_location);
			}
			else
			{
				inner_5node(j, j - offset, 4, -(Nx_3 + Nx_2), Num5 + Num2 + Num3 + Num4, Coef, Coef_location);// 4dx!
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!!!!!!!!!
	//cout << "circle num 2:" << j << endl;
	for (i = (Nx_2 + Nx_3 + Num1); i<(Nx_2 + Nx_3 + Num1 + Num2 / 2); i++)// num 2 :  26 node! part1
	{
		switch (i - Nx_2 - Nx_3 - Num1)
		{
		case 0:// 3 lei 2j shang
			kind3_node(j, j - offset, -(Num1 - Node11) - 1, -(Num1 - Node11), 0, 1, Num2 + 5 - 2, Num2 + 5,
				C4, C3_2, C1_2, C2_2, C4, C3_2, Coef, Coef_location);
			break;
		case (Num2 / 2 - 1) ://node2: 3 lei 2j xia
			kind3_node(j, j - offset, -(Num1 - Node12 + Num2 / 2 - 1), -(Num1 - Node12 + Num2 / 2 - 1) + 1, -1, 0, Num2 + 5, Num2 + 5 + 2,
			C3_2, C4, C2_2, C1_2, C3_2, C4, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1) % 2 == 0)// 1dx
			{
				inner_5node(j, j - offset, 1, -((i - Nx_2 - Nx_3 - Num1) + Num1 - Node11 - (i - Nx_2 - Nx_3 - Num1) / 2), Num2 + 5, Coef, Coef_location);
			}
			else//1 lei  2j left
			{
				kind1_node(j, j - offset, -((i - Nx_2 - Nx_3 - Num1) + Num1 - Node11 - (i - Nx_2 - Nx_3 - Num1) / 2) - (Nx_3 + N_d0 - 2),
					-((i - Nx_2 - Nx_3 - Num1) + Num1 - Node11 - (i - Nx_2 - Nx_3 - Num1) / 2) - (Nx_3 + N_d0 - 2) + 1,
					-((i - Nx_2 - Nx_3 - Num1) + Num1 - Node11 - (i - Nx_2 - Nx_3 - Num1) / 2),
					-((i - Nx_2 - Nx_3 - Num1) + Num1 - Node11 - (i - Nx_2 - Nx_3 - Num1) / 2) + 1,
					-1, 0, 1, 5 + Num2, A2_2, A2_2, A1_2, A1_2, A3_2, A5_2, A3_2, A4_2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 / 2); i<(Nx_2 + Nx_3 + Num1 + Num2); i++)// num 2 :  26 node! part2
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 / 2)
		{
		case 0:// 3 lei 2j shang
			kind3_node(j, j - offset, -(Num1 - Node13 + Num2 / 2) - 1, -(Num1 - Node13 + Num2 / 2), 0, 1, Num3 - 5 - 2, Num3 - 5,
				C4, C3_2, C1_2, C2_2, C4, C3_2, Coef, Coef_location);
			break;
		case (Num2 / 2 - 1) :// 3 lei 2j xia
			kind3_node(j, j - offset, -(Num1 - Node14 + Num2 - 1), -(Num1 - Node14 + Num2 - 1) + 1, -1, 0, Num3 - 5, Num3 - 5 + 2,
			C3_2, C4, C2_2, C1_2, C3_2, C4, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 / 2) % 2 == 0)// 1dx
			{
				inner_5node(j, j - offset, 1, -((i - Nx_2 - Nx_3 - Num1 - Num2 / 2) + Num1 - Node13 - (i - Nx_2 - Nx_3 - Num1 - Num2 / 2) / 2) - Num2 / 2,
					Num3 - 5, Coef, Coef_location);
			}
			else//1 lei 2j left
			{
				kind1_node(j, j - offset, -((i - Nx_2 - Nx_3 - Num1 - Num2 / 2) + Num1 - Node13 - (i - Nx_2 - Nx_3 - Num1 - Num2 / 2) / 2) - Num2 / 2 - (Node_11 + 1),
					-((i - Nx_2 - Nx_3 - Num1 - Num2 / 2) + Num1 - Node13 - (i - Nx_2 - Nx_3 - Num1 - Num2 / 2) / 2) - Num2 / 2 - (Node_11 + 1) + 1,
					-((i - Nx_2 - Nx_3 - Num1 - Num2 / 2) + Num1 - Node13 - (i - Nx_2 - Nx_3 - Num1 - Num2 / 2) / 2) - Num2 / 2,
					-((i - Nx_2 - Nx_3 - Num1 - Num2 / 2) + Num1 - Node13 - (i - Nx_2 - Nx_3 - Num1 - Num2 / 2) / 2) - Num2 / 2 + 1,
					-1, 0, 1, Num3 - 5, A2_2, A2_2, A1_2, A1_2, A3_2, A5_2, A3_2, A4_2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!!!
	//cout << "circle num 3:" << j << endl;
	for (i = (Nx_2 + Nx_3 + Num1 + Num2); i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 / 2); i++)// num 3 :  46 node! part1
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2)
		{
		case 0:// 1 lei 1j shang
			kind1_node(j, j - offset, -(Num1 - N_d0 + 2 + Num2) - 2, -(Num1 - N_d0 + 2 + Num2) - 1, -(Num1 - N_d0 + 2 + Num2), 0, 1, Num3 + Num4 + N_d0 - 2 - 2,
				Num3 + Num4 + N_d0 - 2 - 1, Num3 + Num4 + N_d0 - 2, A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case 3:// 2 lei 2j shang
			kind2_node(j, j - offset, -(Num1 - N_d0 + 2 + Num2), -2, -1, 0, 1, Num3 + Num4 + N_d0 - 2,
				B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case 4:// 3 lei 2j shang
			kind3_node(j, j - offset, -(Num1 - N_d0 + 2 + Num2) - 1, -(Num1 - N_d0 + 2 + Num2), -1, 0, 1, Num3 - 4,
				C4, C4, C3_2, C1_2, C3_2, C2_2, Coef, Coef_location);
			break;
		case (Num3 / 2 - 5) :// 3 lei 2j xia
			kind3_node(j, j - offset, -(Num1 - Node_1 + Num2 + Num3 / 2 - 1), -(Num1 - Node_1 + Num2 + Num3 / 2 - 1) + 1, -1, 0, 1, Num3 - 4,
			C4, C4, C3_2, C1_2, C3_2, C2_2, Coef, Coef_location);
			break;
		case (Num3 / 2 - 4) :// 2 lei 2j xia
			kind2_node(j, j - offset, -(Num1 - Node_1 + Num2 + Num3 / 2 - 1), -1, 0, 1, 2, Num3 / 2 + 1 + Num4 + Node_5,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num3 / 2 - 1) :// 1 lei 1j xia
			kind1_node(j, j - offset, -(Num1 - Node_1 + Num2 + Num3 / 2 - 1), -(Num1 - Node_1 + Num2 + Num3 / 2 - 1) + 1, -(Num1 - Node_1 + Num2 + Num3 / 2 - 1) + 2, -1, 0,
			Num3 / 2 + 1 + Num4 + Node_5, Num3 / 2 + 1 + Num4 + Node_5 + 1, Num3 / 2 + 1 + Num4 + Node_5 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2) < 3)// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num1 - N_d0 + 2 + Num2), Num3 + Num4 + N_d0 - 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2) > 4) && ((i - Nx_2 - Nx_3 - Num1 - Num2) < (Num3 / 2 - 5)))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num2 + 5), Num3 - 4, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num1 - Node_1 + Num2 + Num3 / 2 - 1), Num3 / 2 + 1 + Num4 + Node_5, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 / 2); i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3); i++)// num 3 :  46 node! part2
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 / 2)
		{
		case 0:// 1 lei 1j shang
			kind1_node(j, j - offset, -(Num1 - Node_1 - N_a + 3 + Num2 + Num3 / 2) - 2, -(Num1 - Node_1 - N_a + 3 + Num2 + Num3 / 2) - 1, -(Num1 - Node_1 - N_a + 3 + Num2 + Num3 / 2), 0, 1,
				Num3 / 2 + Num4 + Node_5 + N_a - 3 - 2, Num3 / 2 + Num4 + Node_5 + N_a - 3 - 1, Num3 / 2 + Num4 + Node_5 + N_a - 3,
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case 3:// 2 lei 2j shang
			kind2_node(j, j - offset, -(Num1 - Node_1 - N_a + 3 + Num2 + Num3 / 2), -2, -1, 0, 1, Num3 / 2 + Num4 + Node_5 + N_a - 3,
				B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case 4:// 3 lei 2j shang
			kind3_node(j, j - offset, -(Num1 - Node_1 - N_a + 3 + Num2 + Num3 / 2) - 1, -(Num1 - Node_1 - N_a + 3 + Num2 + Num3 / 2), -1, 0, 1, 5 + Num4 - 1,
				C4, C4, C3_2, C1_2, C3_2, C2_2, Coef, Coef_location);
			break;
		case (Num3 / 2 - 5) :// 3 lei 2j xia
			kind3_node(j, j - offset, -(Num1 - Node_11 + Num2 + Num3 - 1), -(Num1 - Node_11 + Num2 + Num3 - 1) + 1, -1, 0, 1, 5 + Num4 - 1,
			C4, C4, C3_2, C1_2, C3_2, C2_2, Coef, Coef_location);
			break;
		case (Num3 / 2 - 4) :// 2 lei 2j xia
			kind2_node(j, j - offset, -(Num1 - Node_11 + Num2 + Num3 - 1), -1, 0, 1, 2, 1 + Num4 + Node_55,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num3 / 2 - 1) :// 1 lei 1j xia
			kind1_node(j, j - offset, -(Num1 - Node_11 + Num2 + Num3 - 1), -(Num1 - Node_11 + Num2 + Num3 - 1) + 1, -(Num1 - Node_11 + Num2 + Num3 - 1) + 2, -1, 0,
			1 + Num4 + Node_55, 1 + Num4 + Node_55 + 1, 1 + Num4 + Node_55 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 / 2) < 3)// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num1 - Node_1 - N_a + 3 + Num2 + Num3 / 2), Num3 / 2 + Num4 + Node_5 + N_a - 3, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 / 2) > 4) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 / 2) < (Num3 / 2 - 5)))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num3 - 5), 5 + Num4 - 1, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num1 - Node_11 + Num2 + Num3 - 1), 1 + Num4 + Node_55, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	//cout << "circle num 4:" << j << endl;
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3); i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 / 2); i++)// num 4 :  30 node! part1
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3)
		{
		case 0:// 1 lei 2j shang
			kind1_node(j, j - offset, -(Num3 - 4) - 2, -(Num3 - 4) - 1, -(Num3 - 4), 0, 1, Num4 + N_d0 + 2 - 2, Num4 + N_d0 + 2 - 1, Num4 + N_d0 + 2,
				A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case (Num4 / 2 - 1) :// 1 lei 2j xia
			kind1_node(j, j - offset, -(Num3 - 4), -(Num3 - 4) + 1, -(Num3 - 4) + 2, -1, 0,
			Num4 + N_d0 + 2 - N_c1, Num4 + N_d0 + 2 - N_c1 + 1, Num4 + N_d0 + 2 - N_c1 + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3) <= Node41)// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num3 - 4), Num4 + N_d0 + 2, Coef, Coef_location);
			}
			else if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3) >= Node42) // 1dx
			{
				inner_5node(j, j - offset, 1, -(Num3 - 4), Num4 + N_d0 + 2 - N_c1, Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_4node(j, j - offset, -(Num3 - 4), -1, 0, 1, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 / 2); i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4); i++)// num 4 :  30 node! part2
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3)
		{
		case (Num4 / 2) :// 1 lei 2j shang
			kind1_node(j, j - offset, -(5 + Num4 - 1) - 2, -(5 + Num4 - 1) - 1, -(5 + Num4 - 1), 0, 1,
			(Node_55 - 4 + 1 + N_c1) - 2, (Node_55 - 4 + 1 + N_c1) - 1, (Node_55 - 4 + 1 + N_c1),
			A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case (Num4 - 1) :// 1 lei 2j xia
			kind1_node(j, j - offset, -(5 + Num4 - 1), -(5 + Num4 - 1) + 1, -(5 + Num4 - 1) + 2, -1, 0,
			(Node_55 - 4 + 1), (Node_55 - 4 + 1) + 1, (Node_55 - 4 + 1) + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3) <= Node43)// 1dx
			{
				inner_5node(j, j - offset, 1, -(5 + Num4 - 1), (Node_55 - 4 + 1 + N_c1), Coef, Coef_location);
			}
			else if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3) >= Node44) // 1dx
			{
				inner_5node(j, j - offset, 1, -(5 + Num4 - 1), (Node_55 - 4 + 1), Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_4node(j, j - offset, -(5 + Num4 - 1), -1, 0, 1, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	//cout << "circle num5:" << j << endl;
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4); i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Node_5 + 2); i++)// circle num5 :  80 node!  part1 !
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4)
		{
		case 0:
			boundary_node(j, j - offset, 1, -Num1 - Num2 - Num3 - Num4, Num5 + Num6 + Num7 + Num8, Coef, Coef_location);
			break;
		case (N_d0 - 3) ://2 lei shang 1j
			kind2_node(j, j - offset, -Num1 - Num2 - Num3 - Num4, -2, -1, 0, 1, Num5 + Num6 + Num7 + Num8,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (N_d0 + 1) :// 2 lei 2j shang
			kind2_node(j, j - offset, -(Num3 + Num4 + N_d0 - 2), -2, -1, 0, 1, Num5 - N_d0 + 2 + Num6,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node51:// inner 3 node
			inner_3node(j, j - offset, -(Num4 + N_d0 + 2), -1, 0, 1, 1, k*k*dx*dx - 4, Coef, Coef_location);
			break;
		case (Node51 + 1) :// inner 3 node
			inner_3node(j, j - offset, -(Num4 + N_d0 + 2) + N_c1, 0, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node_5 - 3) ://2 lei xia 2j
			kind2_node(j, j - offset, -(Num3 / 2 + 1 + Num4 + Node_5), -1, 0, 1, 2, Num5 - Node_5 + Num6 + Num7 / 2 - 1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_5 + 1) ://2 lei xia 1j
			kind2_node(j, j - offset, -(Num1 - Node_1 + Num2 + Num3 + Num4 + Node_5), -1, 0, 1, 2, Num5 - Node_5 + Num6 + Num7 + Num8 + Node_9,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < (N_d0 - 3))// 4dx!
			{
				inner_5node(j, j - offset, 4, -Num1 - Num2 - Num3 - Num4, Num5 + Num6 + Num7 + Num8, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) > (N_d0 - 3)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < (N_d0 + 1)))// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Num3 + Num4 + N_d0 - 2), Num5 - N_d0 + 2 + Num6, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) > (N_d0 + 1)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < Node51))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num4 + N_d0 + 2), Num5 - N_d0 - 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) > (Node51 + 1)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < (Node_5 - 3)))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num4 + N_d0 + 2) + N_c1, Num5 - N_d0 - 2 - 2, Coef, Coef_location);
			}
			else// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Num3 / 2 + 1 + Num4 + Node_5), Num5 - Node_5 + Num6 + Num7 / 2 - 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Node_5 + 2); i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5); i++)// circle num5 :  80 node!  part2 !
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4)
		{
		case (Num5 - 1) :
			boundary_node(j, j - offset, 2, -Num2 - Num3 - Num4 - Num5, Num6 + Num7 + Num8 + Num9, Coef, Coef_location);
			break;
		case (Node_5 + N_a - 4) ://2 lei shang 1j
			kind2_node(j, j - offset, -(Num1 - Node_1 + Num2 + Num3 + Num4 + Node_5), -2, -1, 0, 1, Num5 - Node_5 + Num6 + Num7 + Num8 + Node_9,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_5 + N_a) :// 2 lei 2j shang
			kind2_node(j, j - offset, -(Num3 / 2 + Num4 + Node_5 + N_a - 3), -2, -1, 0, 1, Num5 - Node_5 - N_a + 3 + Num6 + Num7 / 2,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node52:// inner 3 node
			inner_3node(j, j - offset, -(Node_55 - 4 + 1) - N_c1, -1, 0, 1, 1, k*k*dx*dx - 4, Coef, Coef_location);
			break;
		case (Node52 + 1) :// inner 3 node
			inner_3node(j, j - offset, -(Node_55 - 4 + 1), 0, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node_55 - 3) ://2 lei xia 2j
			kind2_node(j, j - offset, -(Num4 + Node_55 + 1), -1, 0, 1, 2, Num5 - Node_55 + Num6 + Num7 - 1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_55 + 1) ://2 lei xia 1j
			kind2_node(j, j - offset, -Num2 - Num3 - Num4 - Num5, -1, 0, 1, 2, Num6 + Num7 + Num8 + Num9,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < (Node_5 + N_a - 4))// 4dx!
			{
				inner_5node(j, j - offset, 4, -(Num1 - Node_1 + Num2 + Num3 + Num4 + Node_5), Num5 - Node_5 + Num6 + Num7 + Num8 + Node_9, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) > (Node_5 + N_a - 4)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < (Node_5 + N_a)))// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Num3 / 2 + Num4 + Node_5 + N_a - 3), Num5 - Node_5 - N_a + 3 + Num6 + Num7 / 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) > (Node_5 + N_a)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < Node52))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Node_55 - 4 + 1) - N_c1, (Num5 - Node_55 + 4 + Num6 - 1) + 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) > (Node52 + 1)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < (Node_55 - 3)))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Node_55 - 4 + 1), (Num5 - Node_55 + 4 + Num6 - 1), Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) > (Node_55 - 3)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < (Node_55 + 1)))// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num4 + Node_55 + 1), Num5 - Node_55 + Num6 + Num7 - 1, Coef, Coef_location);
			}
			else// 4dx!
			{
				inner_5node(j, j - offset, 4, -Num2 - Num3 - Num4 - Num5, Num6 + Num7 + Num8 + Num9, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//**************************************************************
	cout << "middle_1区域节点数:" << j << endl;
	if (j == sizeM)
		cout << "middle_1 passed..." << endl;
	else
		cout << "middle_1 failed..." << endl;
	//修正矩阵，将不在本区域的节点刨除
	for (int i = 0; i != sizeM; i++)
		for (int j = 0; j != N_matrix; j++){
		if (Coef_location[i][j] < 0 || Coef_location[i][j] >= sizeN)
			Coef_location[i][j] = -1;
		}
	if (Debug)
	{
		ofstream outFile_Coef("middle_1.txt");
		if (!outFile_Coef)
		{
			cout << "Can't write the coefficient matrix into file !" << endl;
			exit(1);//Terminate with fault
		}
		else
		{
			int i, j;
			for (i = 0; i < sizeM; i++)
			{
				for (j = 0; j < N_matrix; j++)
				{
					if ((Coef[i][j].real != 0 || Coef[i][j].image != 0) && Coef_location[i][j] >= 0)
						outFile_Coef << Coef[i][j].real << "+j*" << Coef[i][j].image << " (" << i << "," << Coef_location[i][j] << ") ";
				}
				outFile_Coef << endl << endl;
			}
		}
		outFile_Coef.close();
	}

};
void scan_coef_Matrix_middle_11(complex** Coef, int num_total, int ** Coef_location, int sizeM, int sizeN, int offset)
{
	int i = 0, j = 0, m = 0;
	for (j = 0; j<num_total; j++)         //系数矩阵初始置零
	{
		for (i = 0; i<N_matrix; i++)
		{
			Coef[j][i].real = 0;
			Coef[j][i].image = 0;
			Coef_location[j][i] = -1;//初始化为负值
		}
	}
	j = 0;
	//系数矩阵扫描
	//******************************************************************************************
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9);
		i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Node_5 + 2); i++)// circle num13(5) :  80 node!  part1 !
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9)
		{
		case 0:
			boundary_node(j, j - offset, 1, -Num9 - Num8 - Num7 - Num6, Num5 + Num4 + Num3 + Num2, Coef, Coef_location);
			break;
		case (N_d0 - 3) ://2 lei shang 1j
			kind2_node(j, j - offset, -Num9 - Num8 - Num7 - Num6, -2, -1, 0, 1, Num5 + Num4 + Num3 + Num2,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (N_d0 + 1) :// 2 lei 2j shang
			kind2_node(j, j - offset, -(Num7 + Num6 + N_d0 - 2), -2, -1, 0, 1, Num5 - N_d0 + 2 + Num4,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node51:// inner 3 node
			inner_3node(j, j - offset, -1, 0, (Num5 - N_d0 - 2), 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node51 + 1) :// inner 3 node
			inner_3node(j, j - offset, 0, 1, (Num5 - N_d0 - 2) + N_c1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Node_5 - 3) ://2 lei xia 2j
			kind2_node(j, j - offset, -(Num7 / 2 + 1 + Num6 + Node_5), -1, 0, 1, 2, Num5 - Node_5 + Num4 + Num3 / 2 - 1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_5 + 1) ://2 lei xia 1j
			kind2_node(j, j - offset, -(Num9 - Node_9 + Num8 + Num7 + Num6 + Node_5), -1, 0, 1, 2, Num5 - Node_5 + Num4 + Num3 + Num2 + Node_1,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < (N_d0 - 3))// 4dx!
			{
				inner_5node(j, j - offset, 4, -Num9 - Num8 - Num7 - Num6, Num5 + Num4 + Num3 + Num2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) > (N_d0 - 3)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < (N_d0 + 1)))// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Num7 + Num6 + N_d0 - 2), Num5 - N_d0 + 2 + Num4, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) > (N_d0 + 1)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < Node51))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num6 + N_d0 + 2), (Num5 - N_d0 - 2), Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) > (Node51 + 1)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < (Node_5 - 3)))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num6 + N_d0 + 2) - 2, (Num5 - N_d0 - 2) + N_c1, Coef, Coef_location);
			}
			else// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Num7 / 2 + 1 + Num6 + Node_5), Num5 - Node_5 + Num4 + Num3 / 2 - 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Node_5 + 2);
		i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9); i++)// circle num13(5) :  80 node!  part2 !
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9)
		{
		case (Num5 - 1) :
			boundary_node(j, j - offset, 2, -Num8 - Num7 - Num6 - Num5, Num4 + Num3 + Num2 + Num1, Coef, Coef_location);
			break;
		case (Node_5 + N_a - 4) ://2 lei shang 1j
			kind2_node(j, j - offset, -(Num9 - Node_9 + Num8 + Num7 + Num6 + Node_5), -2, -1, 0, 1, Num5 - Node_5 + Num4 + Num3 + Num2 + Node_1,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_5 + N_a) :// 2 lei 2j shang
			kind2_node(j, j - offset, -(Num7 / 2 + Num6 + Node_5 + N_a - 3), -2, -1, 0, 1, Num5 - Node_5 - N_a + 3 + Num4 + Num3 / 2,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node52:// inner 3 node
			inner_3node(j, j - offset, -1, 0, (Num5 - Node_55 + 4 + Num4 - 1) - N_c1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node52 + 1) :// inner 3 node
			inner_3node(j, j - offset, 0, 1, (Num5 - Node_55 + 4 + Num4 - 1), k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Node_55 - 3) ://2 lei xia 2j
			kind2_node(j, j - offset, -(Num6 + Node_55 + 1), -1, 0, 1, 2, Num5 - Node_55 + Num4 + Num3 - 1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_55 + 1) ://2 lei xia 1j
			kind2_node(j, j - offset, -Num8 - Num7 - Num6 - Num5, -1, 0, 1, 2, Num4 + Num3 + Num2 + Num1,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < (Node_5 + N_a - 4))// 4dx!
			{
				inner_5node(j, j - offset, 4, -(Num9 - Node_9 + Num8 + Num7 + Num6 + Node_5), Num5 - Node_5 + Num4 + Num3 + Num2 + Node_1, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) > (Node_5 + N_a - 4)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < (Node_5 + N_a)))// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Num7 / 2 + Num6 + Node_5 + N_a - 3), Num5 - Node_5 - N_a + 3 + Num4 + Num3 / 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) > (Node_5 + N_a)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < Node52))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Node_55 - 4 + 1) + 2, (Num5 - Node_55 + 4 + Num4 - 1) - N_c1, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) > (Node52 + 1)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < (Node_55 - 3)))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Node_55 - 4 + 1), (Num5 - Node_55 + 4 + Num4 - 1), Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) > (Node_55 - 3)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < (Node_55 + 1)))// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num6 + Node_55 + 1), Num5 - Node_55 + Num4 + Num3 - 1, Coef, Coef_location);
			}
			else// 4dx!
			{
				inner_5node(j, j - offset, 4, -Num8 - Num7 - Num6 - Num5, Num4 + Num3 + Num2 + Num1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9);
		i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Num4 / 2); i++)// num 14(4) :  30 node! part1
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9)
		{
		case 0:// 1 lei 2j shang
			kind1_node(j, j - offset, -(Num5 - N_d0 - 2) - 2, -(Num5 - N_d0 - 2) - 1, -(Num5 - N_d0 - 2), 0, 1, Num4 + 4 - 2, Num4 + 4 - 1, Num4 + 4,
				A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case (Num4 / 2 - 1) :// 1 lei 2j xia
			kind1_node(j, j - offset, -(Num5 - N_d0 - 2) - N_c1, -(Num5 - N_d0 - 2) - N_c1 + 1, -(Num5 - N_d0 - 2) - N_c1 + 2, -1, 0,
			Num4 + 4, Num4 + 4 + 1, Num4 + 4 + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) <= Node41)// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num5 - N_d0 - 2), Num4 + 4, Coef, Coef_location);
			}
			else if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) >= Node42) // 1dx
			{
				inner_5node(j, j - offset, 1, -(Num5 - N_d0 - 2) - N_c1, Num4 + 4, Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_4node(j, j - offset, -1, 0, 1, Num4 + 4, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Num4 / 2);
		i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9); i++)// num14(4) :  30 node! part2
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9)
		{
		case (Num4 / 2) :// 1 lei 2j shang
			kind1_node(j, j - offset, -(Num5 - Node_55 + 4 + Num4 - 1) + N_c1 - 2, -(Num5 - Node_55 + 4 + Num4 - 1) + N_c1 - 1, -(Num5 - Node_55 + 4 + Num4 - 1) + N_c1, 0, 1,
			(Num3 - 5 + 1) - 2, (Num3 - 5 + 1) - 1, (Num3 - 5 + 1),
			A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case (Num4 - 1) :// 1 lei 2j xia
			kind1_node(j, j - offset, -(Num5 - Node_55 + 4 + Num4 - 1), -(Num5 - Node_55 + 4 + Num4 - 1) + 1, -(Num5 - Node_55 + 4 + Num4 - 1) + 2, -1, 0,
			(Num3 - 5 + 1), (Num3 - 5 + 1) + 1, (Num3 - 5 + 1) + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) <= Node43)// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num5 - Node_55 + 4 + Num4 - 1) + N_c1, (Num3 - 5 + 1), Coef, Coef_location);
			}
			else if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) >= Node44) // 1dx
			{
				inner_5node(j, j - offset, 1, -(Num5 - Node_55 + 4 + Num4 - 1), (Num3 - 5 + 1), Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_4node(j, j - offset, -1, 0, 1, (Num3 - 5 + 1), 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9);
		i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Num3 / 2); i++)// num 15(3) :  46 node! part1
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9)
		{
		case 0:// 1 lei 1j shang
			kind1_node(j, j - offset, -(Num5 - N_d0 + 2 + Num4) - 2, -(Num5 - N_d0 + 2 + Num4) - 1, -(Num5 - N_d0 + 2 + Num4), 0, 1,
				Num3 + Num2 + N_d0 - 2 - 2, Num3 + Num2 + N_d0 - 2 - 1, Num3 + Num2 + N_d0 - 2,
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case 3:// 2 lei 2j shang
			kind2_node(j, j - offset, -(Num5 - N_d0 + 2 + Num4), -2, -1, 0, 1, Num3 + Num2 + N_d0 - 2,
				B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case 4:// 3 lei 2j shang
			kind3_node(j, j - offset, -(Num4 + 4), -1, 0, 1, (Num3 + Num2 + N_d0 - 2) - 1, (Num3 + Num2 + N_d0 - 2),
				C2_2, C3_2, C1_2, C3_2, C4, C4, Coef, Coef_location);
			break;
		case (Num3 / 2 - 5) :// 3 lei 2j xia
			kind3_node(j, j - offset, -(Num4 + 4), -1, 0, 1, Num3 / 2 + 1 + Num2 + Node_1, Num3 / 2 + 1 + Num2 + Node_1 + 1,
			C2_2, C3_2, C1_2, C3_2, C4, C4, Coef, Coef_location);
			break;
		case (Num3 / 2 - 4) :// 2 lei 2j xia
			kind2_node(j, j - offset, -(Num5 - Node_5 + Num4 + Num3 / 2 - 1), -1, 0, 1, 2, Num3 / 2 + 1 + Num2 + Node_1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num3 / 2 - 1) :// 1 lei 1j xia
			kind1_node(j, j - offset, -(Num5 - Node_5 + Num4 + Num3 / 2 - 1), -(Num5 - Node_5 + Num4 + Num3 / 2 - 1) + 1, -(Num5 - Node_5 + Num4 + Num3 / 2 - 1) + 2, -1, 0,
			Num3 / 2 + 1 + Num2 + Node_1, Num3 / 2 + 1 + Num2 + Node_1 + 1, Num3 / 2 + 1 + Num2 + Node_1 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < 3)// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num5 - N_d0 + 2 + Num4), Num3 + Num2 + N_d0 - 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) > 4) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < (Num3 / 2 - 5)))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num4 + 4), Num3 - 5, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num5 - Node_5 + Num4 + Num3 / 2 - 1), Num3 / 2 + 1 + Num2 + Node_1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Num3 / 2);
		i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 * 2 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9); i++)// num 15(3) :  46 node! part2
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num3 / 2)
		{
		case 0:// 1 lei 1j shang
			kind1_node(j, j - offset, -(Num5 - Node_5 - N_a + 3 + Num4 + Num3 / 2) - 2, -(Num5 - Node_5 - N_a + 3 + Num4 + Num3 / 2) - 1, -(Num5 - Node_5 - N_a + 3 + Num4 + Num3 / 2), 0, 1,
				Num3 / 2 + Num2 + Node_1 + N_a - 3 - 2, Num3 / 2 + Num2 + Node_1 + N_a - 3 - 1, Num3 / 2 + Num2 + Node_1 + N_a - 3,
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case 3:// 2 lei 2j shang
			kind2_node(j, j - offset, -(Num5 - Node_5 - N_a + 3 + Num4 + Num3 / 2), -2, -1, 0, 1, Num3 / 2 + Num2 + Node_1 + N_a - 3,
				B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case 4:// 3 lei 2j shang
			kind3_node(j, j - offset, -(Num3 - 5 + 1), -1, 0, 1, (Num3 / 2 + Num2 + Node_1 + N_a - 3) - 1, (Num3 / 2 + Num2 + Node_1 + N_a - 3),
				C2_2, C3_2, C1_2, C3_2, C4, C4, Coef, Coef_location);
			break;
		case (Num3 / 2 - 5) :// 3 lei 2j xia
			kind3_node(j, j - offset, -(Num3 - 5 + 1), -1, 0, 1, 1 + Num2 + Node_11, 1 + Num2 + Node_11 + 1,
			C2_2, C3_2, C1_2, C3_2, C4, C4, Coef, Coef_location);
			break;
		case (Num3 / 2 - 4) :// 2 lei 2j xia
			kind2_node(j, j - offset, -(Num5 - Node_55 + Num4 + Num3 - 1), -1, 0, 1, 2, 1 + Num2 + Node_11,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num3 / 2 - 1) :// 1 lei 1j xia
			kind1_node(j, j - offset, -(Num5 - Node_55 + Num4 + Num3 - 1), -(Num5 - Node_55 + Num4 + Num3 - 1) + 1, -(Num5 - Node_55 + Num4 + Num3 - 1) + 2, -1, 0,
			1 + Num2 + Node_11, 1 + Num2 + Node_11 + 1, 1 + Num2 + Node_11 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num3 / 2) < 3)// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num5 - Node_5 - N_a + 3 + Num4 + Num3 / 2), Num3 / 2 + Num2 + Node_1 + N_a - 3, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num3 / 2) > 4) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num3 / 2) < (Num3 / 2 - 5)))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num3 - 5 + 1), 6 + Num2 - 1, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num5 - Node_55 + Num4 + Num3 - 1), 1 + Num2 + Node_11, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 * 2 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9);
		i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 * 2 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Num2 / 2); i++)// num 16(2) :  26 node! part1
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9)
		{
		case 0:// 3 lei 2j shang
			kind3_node(j, j - offset, -(Num3 - 3), -(Num3 - 3) + 2, 0, 1, Num2 + Node11 - 1, Num2 + Node11,
				C4, C3_2, C1_2, C2_2, C4, C3_2, Coef, Coef_location);
			break;
		case (Num2 / 2 - 1) ://node2: 3 lei 2j xia
			kind3_node(j, j - offset, -(Num3 - 5), -(Num3 - 5) + 2, -1, 0, Num2 / 2 + Node_1 - 4 + 1, Num2 / 2 + Node_1 - 4 + 1 + 1,
			C3_2, C4, C2_2, C1_2, C3_2, C4, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) % 2 == 0)// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num3 - 5),
					(Num2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) + Node11 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) / 2), Coef, Coef_location);
			}
			else//1 lei  2j left
			{
				kind1_node(j, j - offset, -(Num3 - 5), -1, 0, 1,
					(Num2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) + Node11 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) / 2),
					(Num2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) + Node11 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) / 2) + 1,
					(Num2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) + Node11 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) / 2) + (Num1 - N_d0 + 2),
					(Num2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) + Node11 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) / 2) + 1 + (Num1 - N_d0 + 2),
					A4_2, A3_2, A5_2, A3_2, A1_2, A1_2, A2_2, A2_2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 * 2 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Num2 / 2);
		i<(Nx_2 + Nx_3 + Num1 + Num2 * 2 + Num3 * 2 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9); i++)// num 16(2) :  26 node! part2
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2)
		{
		case 0:// 3 lei 2j shang
			kind3_node(j, j - offset, -(6 + Num2 - 1) - 2, -(6 + Num2 - 1), 0, 1, Num2 / 2 + Node13 - 1, Num2 / 2 + Node13,
				C4, C3_2, C1_2, C2_2, C4, C3_2, Coef, Coef_location);
			break;
		case (Num2 / 2 - 1) :// 3 lei 2j xia
			kind3_node(j, j - offset, -(6 + Num2 - 1), -(6 + Num2 - 1) + 2, -1, 0, Node_11 - 4 + 1, Node_11 - 4 + 1 + 1,
			C3_2, C4, C2_2, C1_2, C3_2, C4, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) % 2 == 0)// 1dx
			{
				inner_5node(j, j - offset, 1, -(6 + Num2 - 1), Num2 / 2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) + Node13 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) / 2,
					Coef, Coef_location);
			}
			else//1 lei 2j left
			{
				kind1_node(j, j - offset, -(6 + Num2 - 1), -1, 0, 1,
					(Num2 / 2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) + Node13 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) / 2),
					(Num2 / 2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) + Node13 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) / 2) + 1,
					(Num2 / 2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) + Node13 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) / 2) + (Num1 - Node_11 + Nx_3 - 1),
					(Num2 / 2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) + Node13 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) / 2) + 1 + (Num1 - Node_11 + Nx_3 - 1),
					A4_2, A3_2, A5_2, A3_2, A1_2, A1_2, A2_2, A2_2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (num_total - Nx_2 - Nx_3 - Num1); i<(num_total - Nx_2 - Nx_3 - Num1 + Node_1 + 2); i++)// circle num17(1) Nx_2:  74 node!  part1 !
	{
		switch (i - num_total + Nx_2 + Nx_3 + Num1)
		{
		case 0:
			boundary_node(j, j - offset, 1, -Num5 - Num4 - Num3 - Num2, Nx_2 + Nx_3, Coef, Coef_location);
			break;
		case (N_d0 - 3) ://2 lei shang 1j
			kind2_node(j, j - offset, -Num5 - Num4 - Num3 - Num2, -2, -1, 0, 1, Nx_2 + Nx_3,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_1 + 1) ://2 lei xia 1j
			kind2_node(j, j - offset, -(Num5 - Node_5 + Num4 + Num3 + Num2 + Node_1), -1, 0, 1, 2, Nx_2 + Nx_3,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - num_total + Nx_2 + Nx_3 + Num1) < (N_d0 - 3))// 4dx!
			{
				inner_5node(j, j - offset, 4, -Num5 - Num4 - Num3 - Num2, Nx_2 + Nx_3, Coef, Coef_location);
			}
			else if (((i - num_total + Nx_2 + Nx_3 + Num1) > (N_d0 - 3)) && ((i - num_total + Nx_2 + Nx_3 + Num1) < Node11))// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Num3 + Num2 + N_d0 - 2), Num1 - N_d0 + 2, Coef, Coef_location);
			}
			else if (((i - num_total + Nx_2 + Nx_3 + Num1) >= Node11) && ((i - num_total + Nx_2 + Nx_3 + Num1) <= Node12))// 2 lei right 2j !
			{
				kind2_node(j, j - offset, -((i - num_total + Nx_2 + Nx_3 + Num1) + Num2 - (i - num_total + Nx_2 + Nx_3 + Num1 - Node11) * 2),
					-1, 0, 1, Num1 - N_d0 + 2, Nx_2 + Nx_3, B2, B3, B1_2, B3, B4, B5, Coef, Coef_location);
			}
			else// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Num3 / 2 + 1 + Num2 + Node_1), Num1 - N_d0 + 2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (num_total - Nx_2 - Nx_3 - Num1 + Node_1 + 2); i<(num_total - Nx_2 - Nx_3); i++)// circle num17(1) Nx_2:  74 node!  part2 !
	{
		switch (i - num_total + Nx_2 + Nx_3 + Num1)
		{
		case (Num1 - 1) :
			boundary_node(j, j - offset, 2, -Num1 - Num2 - Num3 - Num4, Nx_2 + Nx_3, Coef, Coef_location);
			break;
		case (Node_1 + N_a - 4) ://2 lei shang 1j
			kind2_node(j, j - offset, -(Num5 - Node_5 + Num4 + Num3 + Num2 + Node_1), -2, -1, 0, 1, Nx_2 + Nx_3,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_11 + 1) ://2 lei xia 1j
			kind2_node(j, j - offset, -Num1 - Num2 - Num3 - Num4, -1, 0, 1, 2, Nx_2 + Nx_3,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - num_total + Nx_2 + Nx_3 + Num1) < (Node_1 + N_a - 4))// 4dx!
			{
				inner_5node(j, j - offset, 4, -(Num5 - Node_5 + Num4 + Num3 + Num2 + Node_1), Nx_2 + Nx_3, Coef, Coef_location);
			}
			else if (((i - num_total + Nx_2 + Nx_3 + Num1) > (Node_1 + N_a - 4)) && ((i - num_total + Nx_2 + Nx_3 + Num1) < Node13))// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Num3 / 2 + Num2 + Node_1 + N_a - 3), Num1 - Node_11 + Nx_3 - 1, Coef, Coef_location);
			}
			else if (((i - num_total + Nx_2 + Nx_3 + Num1) >= Node13) && ((i - num_total + Nx_2 + Nx_3 + Num1) <= Node14))// 2 lei left 2j !
			{
				kind2_node(j, j - offset, -((i - num_total + Nx_2 + Nx_3 + Num1) + Num2 / 2 - (i - num_total + Nx_2 + Nx_3 + Num1 - Node13) * 2),
					-1, 0, 1, Num1 - Node_11 + Nx_3 - 1, Nx_2 + Nx_3, B2, B3, B1_2, B3, B4, B5, Coef, Coef_location);
			}
			else if (((i - num_total + Nx_2 + Nx_3 + Num1) >= Node14) && ((i - num_total + Nx_2 + Nx_3 + Num1) <= Node_11))// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num2 + Node_11 + 1), Num1 - Node_11 + Nx_3 - 1, Coef, Coef_location);
			}
			else
			{
				inner_5node(j, j - offset, 4, -Num1 - Num2 - Num3 - Num4, Nx_2 + Nx_3, Coef, Coef_location);// 4dx!
			}
			break;
		}
		j = j + 1;
	}
	//cout << "offset:" << j << endl;
	//====================================================================================================
	//调用封装函数
	//left 1 --- right 1
	scan_coef_Matrix_middle_wrapper(j, Coef, Coef_location, num_total, sizeN - 2 * Num1, offset);
	//=======================================================================================================================
	//***************************************************************************************
	for (i = (Nx_2 + Nx_3); i<(Nx_2 * 2 + Nx_3); i++)//left 2 Nx_2:  74 node! 
	{
		switch (i - Nx_2 - Nx_3)
		{
		case 0:
			boundary_node(j, j - offset, 1, -Nx_2 - Nx_3, Nx_2 + Nx_1, Coef, Coef_location);
			break;
		case (Nx_2 - 1) :
			boundary_node(j, j - offset, 2, -Nx_3 - Nx_2, Nx_1 + Nx_side, Coef, Coef_location);
			break;
		case (N_d0 - 3) ://2 lei shang
			kind2_node(j, j - offset, -Nx_2 - Nx_3, -2, -1, 0, 1, Nx_2 + Nx_1,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (N_d0 - 2) ://3 lei shang
			kind3_node(j, j - offset, -(Nx_3 + N_d0 - 2), -1, 0, 1, Nx_2 + Nx_1 - 1, Nx_2 + Nx_1,
			C2_1, C3_1, C1_1, C3_1, C4, C4, Coef, Coef_location);
			break;
		case Node_1://3 lei xia
			kind3_node(j, j - offset, -(Nx_3 + N_d0 - 2), -1, 0, 1, Nx_2 - Node_1 + Nx_1 + N_d0 + 1, Nx_2 - Node_1 + Nx_1 + N_d0 + 1 + 1,
				C2_1, C3_1, C1_1, C3_1, C4, C4, Coef, Coef_location);
			break;
		case (Node_1 + 1) ://2 lei xia
			kind2_node(j, j - offset, -(Nx_2 + Nx_3), -1, 0, 1, 2, Nx_2 - Node_1 + Nx_1 + N_d0 + 1,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_1 + N_a - 4) ://2 lei shang
			kind2_node(j, j - offset, -(Nx_2 + Nx_3), -2, -1, 0, 1, Nx_2 - Node_1 + Nx_1 + N_d0 + 1,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_1 + N_a - 3) ://3 lei shang
			kind3_node(j, j - offset, -(Node_11 + 1), -1, 0, 1, Nx_2 - Node_1 + Nx_1 + N_d0 + 1 - 1, Nx_2 - Node_1 + Nx_1 + N_d0 + 1,
			C2_1, C3_1, C1_1, C3_1, C4, C4, Coef, Coef_location);
			break;
		case Node_11://3 lei xia
			kind3_node(j, j - offset, -(Node_11 + 1), -1, 0, 1, Nx_side + Nx_1, Nx_side + Nx_1 + 1,
				C2_1, C3_1, C1_1, C3_1, C4, C4, Coef, Coef_location);
			break;
		case (Node_11 + 1) ://2 lei xia
			kind2_node(j, j - offset, -Nx_3 - Nx_2, -1, 0, 1, 2, Nx_1 + Nx_side,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3) < (N_d0 - 3))// 2dx!
			{
				inner_5node(j, j - offset, 4, -Nx_3 - Nx_2, Nx_2 + Nx_1, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3) >= (N_d0 - 1)) && ((i - Nx_2 - Nx_3) < Node_1))
			{
				inner_5node(j, j - offset, 2, -(Nx_3 + N_d0 - 2), Nx_2 - N_d0 + 1, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3) > (Node_1 + 1)) && ((i - Nx_2 - Nx_3) < (Node_1 + N_a - 4)))// 2dx!
			{
				inner_5node(j, j - offset, 4, -(Nx_2 + Nx_3), Nx_2 - Node_1 + Nx_1 + N_d0 + 1, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3) >= (Node_1 + N_a - 2)) && ((i - Nx_2 - Nx_3) < Node_11))
			{
				inner_5node(j, j - offset, 2, -(Node_11 + 1), Nx_1 + Nx_2 - Node_11, Coef, Coef_location);
			}
			else// 2dx!
			{
				inner_5node(j, j - offset, 4, -Nx_3 - Nx_2, Nx_side + Nx_1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (Nx_2 * 2 + Nx_3); i<(Nx_2 * 2 + Nx_3 + Nx_1 / 2); i++)// left 1 Nx_1: 26 node! part1
	{
		switch (i - Nx_2 * 2 - Nx_3)
		{
		case 0:// 3 lei
			kind3_node(j, j - offset, -(Nx_2 - N_d0 + 1) - 2, -(Nx_2 - N_d0 + 1), 0, 1, Nx_1 + N_d0 - 3, Nx_1 + N_d0 - 2,
				C4, C3_1, C1_1, C2_1, C4, C3_1, Coef, Coef_location);
			break;
		case (Nx_1 / 2 - 1) ://node2: 3 lei
			kind3_node(j, j - offset, -(Nx_2 - N_d0 + 1), -(Nx_2 - N_d0 + 1) + 2, -1, 0, Nx_1 / 2 + 1 + N_d0 + 1, Nx_1 / 2 + 1 + N_d0 + 2,
			C3_1, C4, C2_1, C1_1, C3_1, C4, Coef, Coef_location);
			break;
		case 1:// 1 lei
			kind1_node(j, j - offset, -(Nx_2 - N_d0 + 1), -1, 0, 1, (Nx_1 - 1 + N_d0 - 2), (Nx_1 - 1 + N_d0 - 2) + 1, (Nx_1 - 1 + N_d0 - 2) + Nx_side, (Nx_1 - 1 + N_d0 - 2) + 1 + Nx_side,
				A4_1, A3_1, A5_1, A3_1, A1_1, A1_1, A2_1, A2_1, Coef, Coef_location);
			break;
		case 2:// 5 node 2dx
			inner_5node(j, j - offset, 2, -(Nx_2 - N_d0 + 1), Nx_1 - 2 + N_d0 - 1, Coef, Coef_location);
			break;
		case 3:// 1 lei que 2 dian
			kind2_node(j, j - offset, -(Nx_2 - N_d0 + 1), -1, 0, 1, Nx_1 - 3 + N_d0 - 1, Nx_1 - 3 + N_d0 - 1 + Nx_side,
				A4_1, A3_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		case (Nx_1 / 2 - 4) :// 1 lei que 2 dian
			kind2_node(j, j - offset, -(Nx_2 - N_d0 + 1), -1, 0, 1, N_d0 + Nx_1 / 2 + 4, N_d0 + Nx_1 / 2 + 4 + Nx_side,
			A4_1, A3_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		case (Nx_1 / 2 - 3) :// 5 node 2dx
			inner_5node(j, j - offset, 2, -(Nx_2 - N_d0 + 1), Nx_1 / 2 + 3 + N_d0, Coef, Coef_location);
			break;
		case (Nx_1 / 2 - 2) :// 1 lei
			kind1_node(j, j - offset, -(Nx_2 - N_d0 + 1), -1, 0, 1, (Nx_1 / 2 + 2 + N_d0), (Nx_1 / 2 + 2 + N_d0) + 1, (Nx_1 / 2 + 2 + N_d0) + Nx_side, (Nx_1 / 2 + 2 + N_d0) + 1 + Nx_side,
			A4_1, A3_1, A5_1, A3_1, A1_1, A1_1, A2_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 * 2 - Nx_3) % 2 == 0)
			{
				inner_4node(j, j - offset, -(Nx_2 - N_d0 + 1), -1, 0, 1, 1, 1, k*k * 4 * dx*dx - 4, 1, Coef, Coef_location);
			}
			else//1 lei
			{
				inner_4node(j, j - offset, -(Nx_2 - N_d0 + 1), -1, 0, 1, A4_1, A3_1, A5_1, A3_1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 * 2 + Nx_3 + Nx_1 / 2); i<(Nx_2 * 2 + Nx_3 + Nx_1); i++)//left 1 Nx_1: 26 node! part2
	{
		switch (i - Nx_2 * 2 - Nx_3 - Nx_1 / 2)
		{
		case 0:// 3 lei
			kind3_node(j, j - offset, -(Nx_2 - Node_11 + Nx_1) - 2, -(Nx_2 - Node_11 + Nx_1), 0, 1, Nx_1 / 2 + N_d0 + N_a - 3, Nx_1 / 2 + N_d0 + N_a - 2,
				C4, C3_1, C1_1, C2_1, C4, C3_1, Coef, Coef_location);
			break;
		case (Nx_1 / 2 - 1) ://node2: 3 lei
			kind3_node(j, j - offset, -(Nx_2 - Node_11 + Nx_1), -(Nx_2 - Node_11 + Nx_1) + 2, -1, 0, N_d0 + N_a + 1 + 1, N_d0 + N_a + 2 + 1,
			C3_1, C4, C2_1, C1_1, C3_1, C4, Coef, Coef_location);
			break;
		case 1:// 1 lei
			kind1_node(j, j - offset, -(Nx_2 - Node_11 + Nx_1), -1, 0, 1, (Nx_1 / 2 - 1 + N_d0 + N_a - 2), (Nx_1 / 2 - 1 + N_d0 + N_a - 2) + 1,
				(Nx_1 / 2 - 1 + N_d0 + N_a - 2) + Nx_side, (Nx_1 / 2 - 1 + N_d0 + N_a - 2) + 1 + Nx_side,
				A4_1, A3_1, A5_1, A3_1, A1_1, A1_1, A2_1, A2_1, Coef, Coef_location);
			break;
		case 2:// 5 node 2dx
			inner_5node(j, j - offset, 2, -(Nx_2 - Node_11 + Nx_1), Nx_1 / 2 - 2 + N_d0 + N_a - 1, Coef, Coef_location);
			break;
		case 3:// 1 lei que 2 dian
			kind2_node(j, j - offset, -(Nx_2 - Node_11 + Nx_1), -1, 0, 1, Nx_1 / 2 - 3 + N_d0 + N_a - 1, Nx_1 / 2 - 3 + N_d0 + N_a - 1 + Nx_side,
				A4_1, A3_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		case (Nx_1 / 2 - 4) :// 1 lei que 2 dian
			kind2_node(j, j - offset, -(Nx_2 - Node_11 + Nx_1), -1, 0, 1, N_d0 + N_a + 4, N_d0 + N_a + 4 + Nx_side,
			A4_1, A3_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		case (Nx_1 / 2 - 3) :// 5 node 2dx
			inner_5node(j, j - offset, 2, -(Nx_2 - Node_11 + Nx_1), 3 + N_d0 + N_a, Coef, Coef_location);
			break;
		case (Nx_1 / 2 - 2) :// 1 lei
			kind1_node(j, j - offset, -(Nx_2 - Node_11 + Nx_1), -1, 0, 1, (2 + N_d0 + N_a), (2 + N_d0 + N_a) + 1, (2 + N_d0 + N_a) + Nx_side, (2 + N_d0 + N_a) + 1 + Nx_side,
			A4_1, A3_1, A5_1, A3_1, A1_1, A1_1, A2_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 * 2 - Nx_3 - Nx_1 / 2) % 2 == 0)// 2dx
			{
				inner_4node(j, j - offset, -(Nx_2 - Node_11 + Nx_1), -1, 0, 1, 1, 1, k*k * 4 * dx*dx - 4, 1, Coef, Coef_location);
			}
			else//1 lei
			{
				inner_4node(j, j - offset, -(Nx_2 - Node_11 + Nx_1), -1, 0, 1, A4_1, A3_1, A5_1, A3_1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!!!!
	for (i = (Nx_2 * 2 + Nx_3 + Nx_1); i < (Nx_2 * 2 + Nx_3 + Nx_1 + Nx_side); i++)//middle left ,one line
	{
		switch (i - Nx_2 * 2 - Nx_3 - Nx_1)
		{
		case 0:
			boundary_node(j, j - offset, 1, -Nx_2 - Nx_1, Nx_side, Coef, Coef_location);
			break;
		case (Nx_side - 1) :
			boundary_node(j, j - offset, 2, -Nx_1 - Nx_side, Nx_side, Coef, Coef_location);
			break;
		case (N_d0 - 2) :// 2 lei
			kind2_node(j, j - offset, -(Nx_1 + N_d0 - 2), -1, 0, 1, Nx_side, Nx_side * 2, B2, B3, B1_1, B3, B4, B5, Coef, Coef_location);
			break;
		case (N_d0 - 1) :// 2 lei缺一点
			other_5_node(j, j - offset, -(Nx_1 - 2 + N_d0 - 1), -1, 0, Nx_side, 2 * Nx_side, B2, B3, B1_1, B4, B5, Coef, Coef_location);
			break;
		case N_d0:// 2 lei缺一点
			other_5_node(j, j - offset, -(Nx_1 / 2 + 3 + N_d0), 0, 1, Nx_side, 2 * Nx_side, B2, B1_1, B3, B4, B5, Coef, Coef_location);
			break;
		case (N_d0 + 1) :// 2 lei
			kind2_node(j, j - offset, -(Nx_1 / 2 + 1 + N_d0 + 1), -1, 0, 1, Nx_side, Nx_side * 2, B2, B3, B1_1, B3, B4, B5, Coef, Coef_location);
			break;
		case (N_d0 + N_a - 2) :// 2 lei
			kind2_node(j, j - offset, -(Nx_1 / 2 + N_a + N_d0 - 2), -1, 0, 1, Nx_side, Nx_side * 2, B2, B3, B1_1, B3, B4, B5, Coef, Coef_location);
			break;
		case (N_d0 + N_a - 1) :// 2 lei缺一点
			other_5_node(j, j - offset, -(Nx_1 / 2 - 2 + N_d0 + N_a - 1), -1, 0, Nx_side, 2 * Nx_side, B2, B3, B1_1, B4, B5, Coef, Coef_location);
			break;
		case (N_d0 + N_a) :// 2 lei缺一点
			other_5_node(j, j - offset, -(N_d0 + N_a + 3), 0, 1, Nx_side, 2 * Nx_side, B2, B1_1, B3, B4, B5, Coef, Coef_location);
			break;
		case (N_d0 + N_a + 1) :// 2 lei
			kind2_node(j, j - offset, -(1 + N_a + N_d0 + 1), -1, 0, 1, Nx_side, Nx_side * 2, B2, B3, B1_1, B3, B4, B5, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 * 2 - Nx_3 - Nx_1) < (N_d0 - 2))//up
			{
				inner_5node(j, j - offset, 4, -Nx_2 - Nx_1, Nx_side, Coef, Coef_location);
			}
			else if ((i - Nx_2 * 2 - Nx_3 - Nx_1) > (N_d0 + N_a + 1))//down
			{
				inner_5node(j, j - offset, 4, -Nx_1 - Nx_side, Nx_side, Coef, Coef_location);
			}
			else//middle
			{
				inner_5node(j, j - offset, 4, -(Nx_2 - Node_1 + Nx_1 + N_d0 + 1), Nx_side, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//**************************************************************
	cout << "middle_11区域节点数:" << j << endl;
	if (j == sizeM)
		cout << "middle_11 passed..." << endl;
	else
		cout << "middle_11 failed..." << endl;
	//修正矩阵，将不在本区域的节点刨除
	for (int i = 0; i != sizeM; i++)
		for (int j = 0; j != N_matrix; j++){
		if (Coef_location[i][j] < 0 || Coef_location[i][j] >= sizeN)
			Coef_location[i][j] = -1;
		}
	if (Debug)
	{
		ofstream outFile_Coef("middle_11.txt");
		if (!outFile_Coef)
		{
			cout << "Can't write the coefficient matrix into file !" << endl;
			exit(1);//Terminate with fault
		}
		else
		{
			int i, j;
			for (i = 0; i < sizeM; i++)
			{
				for (j = 0; j < N_matrix; j++)
				{
					if ((Coef[i][j].real != 0 || Coef[i][j].image != 0) && Coef_location[i][j] >= 0)
						outFile_Coef << Coef[i][j].real << "+j*" << Coef[i][j].image << " (" << i << "," << Coef_location[i][j] << ") ";
				}
				outFile_Coef << endl << endl;
			}
		}
		outFile_Coef.close();
	}

};
void scan_coef_Matrix_middle_7(complex** Coef, int num_total, int ** Coef_location, int sizeM, int sizeN, int offset)
{
	int i = 0, j = 0, m = 0;
	for (j = 0; j<num_total; j++)         //系数矩阵初始置零
	{
		for (i = 0; i<N_matrix; i++)
		{
			Coef[j][i].real = 0;
			Coef[j][i].image = 0;
			Coef_location[j][i] = -1;//初始化为负值
		}
	}
	j = 0;
	//系数矩阵扫描
	//******************************************************************************************
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9);
		i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Node_5 + 2); i++)// circle num13(5) :  80 node!  part1 !
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9)
		{
		case 0:
			boundary_node(j, j - offset, 1, -Num9 - Num8 - Num7 - Num6, Num5 + Num4 + Num3 + Num2, Coef, Coef_location);
			break;
		case (N_d0 - 3) ://2 lei shang 1j
			kind2_node(j, j - offset, -Num9 - Num8 - Num7 - Num6, -2, -1, 0, 1, Num5 + Num4 + Num3 + Num2,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (N_d0 + 1) :// 2 lei 2j shang
			kind2_node(j, j - offset, -(Num7 + Num6 + N_d0 - 2), -2, -1, 0, 1, Num5 - N_d0 + 2 + Num4,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node51:// inner 3 node
			inner_3node(j, j - offset, -1, 0, (Num5 - N_d0 - 2), 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node51 + 1) :// inner 3 node
			inner_3node(j, j - offset, 0, 1, (Num5 - N_d0 - 2) + N_c1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Node_5 - 3) ://2 lei xia 2j
			kind2_node(j, j - offset, -(Num7 / 2 + 1 + Num6 + Node_5), -1, 0, 1, 2, Num5 - Node_5 + Num4 + Num3 / 2 - 1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_5 + 1) ://2 lei xia 1j
			kind2_node(j, j - offset, -(Num9 - Node_9 + Num8 + Num7 + Num6 + Node_5), -1, 0, 1, 2, Num5 - Node_5 + Num4 + Num3 + Num2 + Node_1,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < (N_d0 - 3))// 4dx!
			{
				inner_5node(j, j - offset, 4, -Num9 - Num8 - Num7 - Num6, Num5 + Num4 + Num3 + Num2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) > (N_d0 - 3)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < (N_d0 + 1)))// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Num7 + Num6 + N_d0 - 2), Num5 - N_d0 + 2 + Num4, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) > (N_d0 + 1)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < Node51))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num6 + N_d0 + 2), (Num5 - N_d0 - 2), Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) > (Node51 + 1)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < (Node_5 - 3)))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num6 + N_d0 + 2) - 2, (Num5 - N_d0 - 2) + N_c1, Coef, Coef_location);
			}
			else// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Num7 / 2 + 1 + Num6 + Node_5), Num5 - Node_5 + Num4 + Num3 / 2 - 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Node_5 + 2);
		i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9); i++)// circle num13(5) :  80 node!  part2 !
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9)
		{
		case (Num5 - 1) :
			boundary_node(j, j - offset, 2, -Num8 - Num7 - Num6 - Num5, Num4 + Num3 + Num2 + Num1, Coef, Coef_location);
			break;
		case (Node_5 + N_a - 4) ://2 lei shang 1j
			kind2_node(j, j - offset, -(Num9 - Node_9 + Num8 + Num7 + Num6 + Node_5), -2, -1, 0, 1, Num5 - Node_5 + Num4 + Num3 + Num2 + Node_1,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_5 + N_a) :// 2 lei 2j shang
			kind2_node(j, j - offset, -(Num7 / 2 + Num6 + Node_5 + N_a - 3), -2, -1, 0, 1, Num5 - Node_5 - N_a + 3 + Num4 + Num3 / 2,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node52:// inner 3 node
			inner_3node(j, j - offset, -1, 0, (Num5 - Node_55 + 4 + Num4 - 1) - N_c1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node52 + 1) :// inner 3 node
			inner_3node(j, j - offset, 0, 1, (Num5 - Node_55 + 4 + Num4 - 1), k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Node_55 - 3) ://2 lei xia 2j
			kind2_node(j, j - offset, -(Num6 + Node_55 + 1), -1, 0, 1, 2, Num5 - Node_55 + Num4 + Num3 - 1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_55 + 1) ://2 lei xia 1j
			kind2_node(j, j - offset, -Num8 - Num7 - Num6 - Num5, -1, 0, 1, 2, Num4 + Num3 + Num2 + Num1,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < (Node_5 + N_a - 4))// 4dx!
			{
				inner_5node(j, j - offset, 4, -(Num9 - Node_9 + Num8 + Num7 + Num6 + Node_5), Num5 - Node_5 + Num4 + Num3 + Num2 + Node_1, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) > (Node_5 + N_a - 4)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < (Node_5 + N_a)))// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Num7 / 2 + Num6 + Node_5 + N_a - 3), Num5 - Node_5 - N_a + 3 + Num4 + Num3 / 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) > (Node_5 + N_a)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < Node52))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Node_55 - 4 + 1) + 2, (Num5 - Node_55 + 4 + Num4 - 1) - N_c1, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) > (Node52 + 1)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < (Node_55 - 3)))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Node_55 - 4 + 1), (Num5 - Node_55 + 4 + Num4 - 1), Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) > (Node_55 - 3)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < (Node_55 + 1)))// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num6 + Node_55 + 1), Num5 - Node_55 + Num4 + Num3 - 1, Coef, Coef_location);
			}
			else// 4dx!
			{
				inner_5node(j, j - offset, 4, -Num8 - Num7 - Num6 - Num5, Num4 + Num3 + Num2 + Num1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9);
		i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Num4 / 2); i++)// num 14(4) :  30 node! part1
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9)
		{
		case 0:// 1 lei 2j shang
			kind1_node(j, j - offset, -(Num5 - N_d0 - 2) - 2, -(Num5 - N_d0 - 2) - 1, -(Num5 - N_d0 - 2), 0, 1, Num4 + 4 - 2, Num4 + 4 - 1, Num4 + 4,
				A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case (Num4 / 2 - 1) :// 1 lei 2j xia
			kind1_node(j, j - offset, -(Num5 - N_d0 - 2) - N_c1, -(Num5 - N_d0 - 2) - N_c1 + 1, -(Num5 - N_d0 - 2) - N_c1 + 2, -1, 0,
			Num4 + 4, Num4 + 4 + 1, Num4 + 4 + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) <= Node41)// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num5 - N_d0 - 2), Num4 + 4, Coef, Coef_location);
			}
			else if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) >= Node42) // 1dx
			{
				inner_5node(j, j - offset, 1, -(Num5 - N_d0 - 2) - N_c1, Num4 + 4, Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_4node(j, j - offset, -1, 0, 1, Num4 + 4, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Num4 / 2);
		i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9); i++)// num14(4) :  30 node! part2
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9)
		{
		case (Num4 / 2) :// 1 lei 2j shang
			kind1_node(j, j - offset, -(Num5 - Node_55 + 4 + Num4 - 1) + N_c1 - 2, -(Num5 - Node_55 + 4 + Num4 - 1) + N_c1 - 1, -(Num5 - Node_55 + 4 + Num4 - 1) + N_c1, 0, 1,
			(Num3 - 5 + 1) - 2, (Num3 - 5 + 1) - 1, (Num3 - 5 + 1),
			A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case (Num4 - 1) :// 1 lei 2j xia
			kind1_node(j, j - offset, -(Num5 - Node_55 + 4 + Num4 - 1), -(Num5 - Node_55 + 4 + Num4 - 1) + 1, -(Num5 - Node_55 + 4 + Num4 - 1) + 2, -1, 0,
			(Num3 - 5 + 1), (Num3 - 5 + 1) + 1, (Num3 - 5 + 1) + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) <= Node43)// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num5 - Node_55 + 4 + Num4 - 1) + N_c1, (Num3 - 5 + 1), Coef, Coef_location);
			}
			else if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) >= Node44) // 1dx
			{
				inner_5node(j, j - offset, 1, -(Num5 - Node_55 + 4 + Num4 - 1), (Num3 - 5 + 1), Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_4node(j, j - offset, -1, 0, 1, (Num3 - 5 + 1), 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9);
		i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Num3 / 2); i++)// num 15(3) :  46 node! part1
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9)
		{
		case 0:// 1 lei 1j shang
			kind1_node(j, j - offset, -(Num5 - N_d0 + 2 + Num4) - 2, -(Num5 - N_d0 + 2 + Num4) - 1, -(Num5 - N_d0 + 2 + Num4), 0, 1,
				Num3 + Num2 + N_d0 - 2 - 2, Num3 + Num2 + N_d0 - 2 - 1, Num3 + Num2 + N_d0 - 2,
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case 3:// 2 lei 2j shang
			kind2_node(j, j - offset, -(Num5 - N_d0 + 2 + Num4), -2, -1, 0, 1, Num3 + Num2 + N_d0 - 2,
				B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case 4:// 3 lei 2j shang
			kind3_node(j, j - offset, -(Num4 + 4), -1, 0, 1, (Num3 + Num2 + N_d0 - 2) - 1, (Num3 + Num2 + N_d0 - 2),
				C2_2, C3_2, C1_2, C3_2, C4, C4, Coef, Coef_location);
			break;
		case (Num3 / 2 - 5) :// 3 lei 2j xia
			kind3_node(j, j - offset, -(Num4 + 4), -1, 0, 1, Num3 / 2 + 1 + Num2 + Node_1, Num3 / 2 + 1 + Num2 + Node_1 + 1,
			C2_2, C3_2, C1_2, C3_2, C4, C4, Coef, Coef_location);
			break;
		case (Num3 / 2 - 4) :// 2 lei 2j xia
			kind2_node(j, j - offset, -(Num5 - Node_5 + Num4 + Num3 / 2 - 1), -1, 0, 1, 2, Num3 / 2 + 1 + Num2 + Node_1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num3 / 2 - 1) :// 1 lei 1j xia
			kind1_node(j, j - offset, -(Num5 - Node_5 + Num4 + Num3 / 2 - 1), -(Num5 - Node_5 + Num4 + Num3 / 2 - 1) + 1, -(Num5 - Node_5 + Num4 + Num3 / 2 - 1) + 2, -1, 0,
			Num3 / 2 + 1 + Num2 + Node_1, Num3 / 2 + 1 + Num2 + Node_1 + 1, Num3 / 2 + 1 + Num2 + Node_1 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < 3)// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num5 - N_d0 + 2 + Num4), Num3 + Num2 + N_d0 - 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) > 4) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < (Num3 / 2 - 5)))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num4 + 4), Num3 - 5, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num5 - Node_5 + Num4 + Num3 / 2 - 1), Num3 / 2 + 1 + Num2 + Node_1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Num3 / 2);
		i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 * 2 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9); i++)// num 15(3) :  46 node! part2
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num3 / 2)
		{
		case 0:// 1 lei 1j shang
			kind1_node(j, j - offset, -(Num5 - Node_5 - N_a + 3 + Num4 + Num3 / 2) - 2, -(Num5 - Node_5 - N_a + 3 + Num4 + Num3 / 2) - 1, -(Num5 - Node_5 - N_a + 3 + Num4 + Num3 / 2), 0, 1,
				Num3 / 2 + Num2 + Node_1 + N_a - 3 - 2, Num3 / 2 + Num2 + Node_1 + N_a - 3 - 1, Num3 / 2 + Num2 + Node_1 + N_a - 3,
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case 3:// 2 lei 2j shang
			kind2_node(j, j - offset, -(Num5 - Node_5 - N_a + 3 + Num4 + Num3 / 2), -2, -1, 0, 1, Num3 / 2 + Num2 + Node_1 + N_a - 3,
				B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case 4:// 3 lei 2j shang
			kind3_node(j, j - offset, -(Num3 - 5 + 1), -1, 0, 1, (Num3 / 2 + Num2 + Node_1 + N_a - 3) - 1, (Num3 / 2 + Num2 + Node_1 + N_a - 3),
				C2_2, C3_2, C1_2, C3_2, C4, C4, Coef, Coef_location);
			break;
		case (Num3 / 2 - 5) :// 3 lei 2j xia
			kind3_node(j, j - offset, -(Num3 - 5 + 1), -1, 0, 1, 1 + Num2 + Node_11, 1 + Num2 + Node_11 + 1,
			C2_2, C3_2, C1_2, C3_2, C4, C4, Coef, Coef_location);
			break;
		case (Num3 / 2 - 4) :// 2 lei 2j xia
			kind2_node(j, j - offset, -(Num5 - Node_55 + Num4 + Num3 - 1), -1, 0, 1, 2, 1 + Num2 + Node_11,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num3 / 2 - 1) :// 1 lei 1j xia
			kind1_node(j, j - offset, -(Num5 - Node_55 + Num4 + Num3 - 1), -(Num5 - Node_55 + Num4 + Num3 - 1) + 1, -(Num5 - Node_55 + Num4 + Num3 - 1) + 2, -1, 0,
			1 + Num2 + Node_11, 1 + Num2 + Node_11 + 1, 1 + Num2 + Node_11 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num3 / 2) < 3)// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num5 - Node_5 - N_a + 3 + Num4 + Num3 / 2), Num3 / 2 + Num2 + Node_1 + N_a - 3, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num3 / 2) > 4) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num3 / 2) < (Num3 / 2 - 5)))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num3 - 5 + 1), 6 + Num2 - 1, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num5 - Node_55 + Num4 + Num3 - 1), 1 + Num2 + Node_11, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 * 2 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9);
		i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 * 2 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Num2 / 2); i++)// num 16(2) :  26 node! part1
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9)
		{
		case 0:// 3 lei 2j shang
			kind3_node(j, j - offset, -(Num3 - 3), -(Num3 - 3) + 2, 0, 1, Num2 + Node11 - 1, Num2 + Node11,
				C4, C3_2, C1_2, C2_2, C4, C3_2, Coef, Coef_location);
			break;
		case (Num2 / 2 - 1) ://node2: 3 lei 2j xia
			kind3_node(j, j - offset, -(Num3 - 5), -(Num3 - 5) + 2, -1, 0, Num2 / 2 + Node_1 - 4 + 1, Num2 / 2 + Node_1 - 4 + 1 + 1,
			C3_2, C4, C2_2, C1_2, C3_2, C4, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) % 2 == 0)// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num3 - 5),
					(Num2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) + Node11 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) / 2), Coef, Coef_location);
			}
			else//1 lei  2j left
			{
				kind1_node(j, j - offset, -(Num3 - 5), -1, 0, 1,
					(Num2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) + Node11 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) / 2),
					(Num2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) + Node11 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) / 2) + 1,
					(Num2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) + Node11 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) / 2) + (Num1 - N_d0 + 2),
					(Num2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) + Node11 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) / 2) + 1 + (Num1 - N_d0 + 2),
					A4_2, A3_2, A5_2, A3_2, A1_2, A1_2, A2_2, A2_2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 * 2 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Num2 / 2);
		i<(Nx_2 + Nx_3 + Num1 + Num2 * 2 + Num3 * 2 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9); i++)// num 16(2) :  26 node! part2
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2)
		{
		case 0:// 3 lei 2j shang
			kind3_node(j, j - offset, -(6 + Num2 - 1) - 2, -(6 + Num2 - 1), 0, 1, Num2 / 2 + Node13 - 1, Num2 / 2 + Node13,
				C4, C3_2, C1_2, C2_2, C4, C3_2, Coef, Coef_location);
			break;
		case (Num2 / 2 - 1) :// 3 lei 2j xia
			kind3_node(j, j - offset, -(6 + Num2 - 1), -(6 + Num2 - 1) + 2, -1, 0, Node_11 - 4 + 1, Node_11 - 4 + 1 + 1,
			C3_2, C4, C2_2, C1_2, C3_2, C4, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) % 2 == 0)// 1dx
			{
				inner_5node(j, j - offset, 1, -(6 + Num2 - 1), Num2 / 2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) + Node13 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) / 2,
					Coef, Coef_location);
			}
			else//1 lei 2j left
			{
				kind1_node(j, j - offset, -(6 + Num2 - 1), -1, 0, 1,
					(Num2 / 2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) + Node13 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) / 2),
					(Num2 / 2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) + Node13 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) / 2) + 1,
					(Num2 / 2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) + Node13 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) / 2) + (Num1 - Node_11 + Nx_3 - 1),
					(Num2 / 2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) + Node13 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) / 2) + 1 + (Num1 - Node_11 + Nx_3 - 1),
					A4_2, A3_2, A5_2, A3_2, A1_2, A1_2, A2_2, A2_2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (num_total - Nx_2 - Nx_3 - Num1); i<(num_total - Nx_2 - Nx_3 - Num1 + Node_1 + 2); i++)// circle num17(1) Nx_2:  74 node!  part1 !
	{
		switch (i - num_total + Nx_2 + Nx_3 + Num1)
		{
		case 0:
			boundary_node(j, j - offset, 1, -Num5 - Num4 - Num3 - Num2, Nx_2 + Nx_3, Coef, Coef_location);
			break;
		case (N_d0 - 3) ://2 lei shang 1j
			kind2_node(j, j - offset, -Num5 - Num4 - Num3 - Num2, -2, -1, 0, 1, Nx_2 + Nx_3,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_1 + 1) ://2 lei xia 1j
			kind2_node(j, j - offset, -(Num5 - Node_5 + Num4 + Num3 + Num2 + Node_1), -1, 0, 1, 2, Nx_2 + Nx_3,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - num_total + Nx_2 + Nx_3 + Num1) < (N_d0 - 3))// 4dx!
			{
				inner_5node(j, j - offset, 4, -Num5 - Num4 - Num3 - Num2, Nx_2 + Nx_3, Coef, Coef_location);
			}
			else if (((i - num_total + Nx_2 + Nx_3 + Num1) > (N_d0 - 3)) && ((i - num_total + Nx_2 + Nx_3 + Num1) < Node11))// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Num3 + Num2 + N_d0 - 2), Num1 - N_d0 + 2, Coef, Coef_location);
			}
			else if (((i - num_total + Nx_2 + Nx_3 + Num1) >= Node11) && ((i - num_total + Nx_2 + Nx_3 + Num1) <= Node12))// 2 lei right 2j !
			{
				kind2_node(j, j - offset, -((i - num_total + Nx_2 + Nx_3 + Num1) + Num2 - (i - num_total + Nx_2 + Nx_3 + Num1 - Node11) * 2),
					-1, 0, 1, Num1 - N_d0 + 2, Nx_2 + Nx_3, B2, B3, B1_2, B3, B4, B5, Coef, Coef_location);
			}
			else// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Num3 / 2 + 1 + Num2 + Node_1), Num1 - N_d0 + 2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (num_total - Nx_2 - Nx_3 - Num1 + Node_1 + 2); i<(num_total - Nx_2 - Nx_3); i++)// circle num17(1) Nx_2:  74 node!  part2 !
	{
		switch (i - num_total + Nx_2 + Nx_3 + Num1)
		{
		case (Num1 - 1) :
			boundary_node(j, j - offset, 2, -Num1 - Num2 - Num3 - Num4, Nx_2 + Nx_3, Coef, Coef_location);
			break;
		case (Node_1 + N_a - 4) ://2 lei shang 1j
			kind2_node(j, j - offset, -(Num5 - Node_5 + Num4 + Num3 + Num2 + Node_1), -2, -1, 0, 1, Nx_2 + Nx_3,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_11 + 1) ://2 lei xia 1j
			kind2_node(j, j - offset, -Num1 - Num2 - Num3 - Num4, -1, 0, 1, 2, Nx_2 + Nx_3,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - num_total + Nx_2 + Nx_3 + Num1) < (Node_1 + N_a - 4))// 4dx!
			{
				inner_5node(j, j - offset, 4, -(Num5 - Node_5 + Num4 + Num3 + Num2 + Node_1), Nx_2 + Nx_3, Coef, Coef_location);
			}
			else if (((i - num_total + Nx_2 + Nx_3 + Num1) > (Node_1 + N_a - 4)) && ((i - num_total + Nx_2 + Nx_3 + Num1) < Node13))// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Num3 / 2 + Num2 + Node_1 + N_a - 3), Num1 - Node_11 + Nx_3 - 1, Coef, Coef_location);
			}
			else if (((i - num_total + Nx_2 + Nx_3 + Num1) >= Node13) && ((i - num_total + Nx_2 + Nx_3 + Num1) <= Node14))// 2 lei left 2j !
			{
				kind2_node(j, j - offset, -((i - num_total + Nx_2 + Nx_3 + Num1) + Num2 / 2 - (i - num_total + Nx_2 + Nx_3 + Num1 - Node13) * 2),
					-1, 0, 1, Num1 - Node_11 + Nx_3 - 1, Nx_2 + Nx_3, B2, B3, B1_2, B3, B4, B5, Coef, Coef_location);
			}
			else if (((i - num_total + Nx_2 + Nx_3 + Num1) >= Node14) && ((i - num_total + Nx_2 + Nx_3 + Num1) <= Node_11))// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num2 + Node_11 + 1), Num1 - Node_11 + Nx_3 - 1, Coef, Coef_location);
			}
			else
			{
				inner_5node(j, j - offset, 4, -Num1 - Num2 - Num3 - Num4, Nx_2 + Nx_3, Coef, Coef_location);// 4dx!
			}
			break;
		}
		j = j + 1;
	}
	//cout << "offset:" << j << endl;
	//====================================================================================================
	//调用封装函数
	//left 1 --- right 1
	scan_coef_Matrix_middle_wrapper(j, Coef, Coef_location, num_total, sizeN - 2 * Num1, offset);
	//=======================================================================================================================
	//******************************************************************************************
	for (i = (Nx_2 + Nx_3); i<(Nx_2 + Nx_3 + Node_1_p2 + 2); i++)// circle num1 Nx_2:  74 node!  part1 !
	{
		switch (i - Nx_2 - Nx_3)
		{
		case 0:
			boundary_node(j, i, 1, -Nx_2 - Nx_3, Num1_p2 + Num2_p2 + Num3_p2 + Num4_p2, Coef, Coef_location);
			break;
		case (N_d0 - 3) ://2 lei shang 1j
			kind2_node(j, i, -Nx_2 - Nx_3, -2, -1, 0, 1, Num1_p2 + Num2_p2 + Num3_p2 + Num4_p2,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_1_p2 + 1) ://2 lei xia 1j
			kind2_node(j, i, -(Nx_2 + Nx_3), -1, 0, 1, 2, Num1_p2 - Node_1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + Node_5_p2,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3) < (N_d0 - 3))// 4dx!
			{
				inner_5node(j, i, 4, -Nx_3 - Nx_2, Num1_p2 + Num2_p2 + Num3_p2 + Num4_p2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3) > (N_d0 - 3)) && ((i - Nx_2 - Nx_3) < Node11_p2))// 2dx!
			{
				inner_5node(j, i, 2, -(Nx_3 + N_d0 - 2), Num1_p2 - N_d0 + 2 + Num2_p2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3) >= Node11_p2) && ((i - Nx_2 - Nx_3) <= Node12_p2))// 2 lei left 2j !
			{
				kind2_node(j, i, -Nx_3 - Nx_2, -(Nx_3 + N_d0 - 2), -1, 0, 1, Num1_p2 - (i - Nx_2 - Nx_3) + ((i - Nx_2 - Nx_3) - Node11_p2) * 2,
					B5, B4, B3, B1_2, B3, B2, Coef, Coef_location);
			}
			else// 2dx!
			{
				inner_5node(j, i, 2, -(Nx_3 + N_d0 - 2), Num1_p2 - Node_1_p2 + Num2_p2 + Num3_p2 / 2 - 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Node_1_p2 + 2); i<(Nx_2 + Nx_3 + Num1_p2); i++)// circle num1 Nx_2:  74 node!  part2 !
	{
		switch (i - Nx_2 - Nx_3)
		{
		case (Num1_p2 - 1) :
			boundary_node(j, i, 2, -Nx_2 - Nx_3, Num5_p2 + Num2_p2 + Num3_p2 + Num4_p2, Coef, Coef_location);
			break;
		case (Node_1_p2 + N_a - 4) ://2 lei shang 1j
			kind2_node(j, i, -Nx_2 - Nx_3, -2, -1, 0, 1, Num1_p2 - Node_1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + Node_5_p2,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_11_p2 + 1) ://2 lei xia 1j
			kind2_node(j, i, -(Nx_2 + Nx_3), -1, 0, 1, 2, Num5_p2 + Num2_p2 + Num3_p2 + Num4_p2,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3) < (Node_1_p2 + N_a - 4))// 4dx!
			{
				inner_5node(j, i, 4, -Nx_3 - Nx_2, Num1_p2 - Node_1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + Node_5_p2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3) > (Node_1_p2 + N_a - 4)) && ((i - Nx_2 - Nx_3) < Node13_p2))// 2dx!
			{
				inner_5node(j, i, 2, -(Node_11_p2 + 1), Num1_p2 - Node_1_p2 - N_a + 3 + Num2_p2 + Num3_p2 / 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3) >= Node13_p2) && ((i - Nx_2 - Nx_3) <= Node14_p2))// 2 lei left 2j !
			{
				kind2_node(j, i, -Nx_3 - Nx_2, -(Node_11_p2 + 1), -1, 0, 1, Num1_p2 - (i - Nx_2 - Nx_3) + Num2_p2 / 2 + ((i - Nx_2 - Nx_3) - Node13_p2) * 2,
					B5, B4, B3, B1_2, B3, B2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3) >= Node14_p2) && ((i - Nx_2 - Nx_3) <= Node_11_p2))// 2dx
			{
				inner_5node(j, i, 2, -(Node_11_p2 + 1), Num1_p2 - Node_11_p2 + Num2_p2 + Num3_p2 - 1, Coef, Coef_location);
			}
			else
			{
				inner_5node(j, i, 4, -(Nx_3 + Nx_2), Num5_p2 + Num2_p2 + Num3_p2 + Num4_p2, Coef, Coef_location);// 4dx!
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!!!!!!!!!
	for (i = (Nx_2 + Nx_3 + Num1_p2); i<(Nx_2 + Nx_3 + Num1_p2 + Num2_p2 / 2); i++)// num 2 :  26 node! part1
	{
		switch (i - Nx_2 - Nx_3 - Num1_p2)
		{
		case 0:// 3 lei 2j shang
			kind3_node(j, i, -(Num1_p2 - Node11_p2) - 1, -(Num1_p2 - Node11_p2), 0, 1, Num2_p2 + 5 - 2, Num2_p2 + 5,
				C4, C3_2, C1_2, C2_2, C4, C3_2, Coef, Coef_location);
			break;
		case (Num2_p2 / 2 - 1) ://node2: 3 lei 2j xia
			kind3_node(j, i, -(Num1_p2 - Node12_p2 + Num2_p2 / 2 - 1), -(Num1_p2 - Node12_p2 + Num2_p2 / 2 - 1) + 1, -1, 0,
			Num2_p2 + 5, Num2_p2 + 5 + 2, C3_2, C4, C2_2, C1_2, C3_2, C4, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1_p2) % 2 == 0)// 1dx
			{
				inner_5node(j, i, 1, -((i - Nx_2 - Nx_3 - Num1_p2) + Num1_p2 - Node11_p2 - (i - Nx_2 - Nx_3 - Num1_p2) / 2), Num2_p2 + 5, Coef, Coef_location);
			}
			else//1 lei  2j left
			{
				kind1_node(j, i, -((i - Nx_2 - Nx_3 - Num1_p2) + Num1_p2 - Node11_p2 - (i - Nx_2 - Nx_3 - Num1_p2) / 2) - (Nx_3 + N_d0 - 2),
					-((i - Nx_2 - Nx_3 - Num1_p2) + Num1_p2 - Node11_p2 - (i - Nx_2 - Nx_3 - Num1_p2) / 2) - (Nx_3 + N_d0 - 2) + 1,
					-((i - Nx_2 - Nx_3 - Num1_p2) + Num1_p2 - Node11_p2 - (i - Nx_2 - Nx_3 - Num1_p2) / 2),
					-((i - Nx_2 - Nx_3 - Num1_p2) + Num1_p2 - Node11_p2 - (i - Nx_2 - Nx_3 - Num1_p2) / 2) + 1,
					-1, 0, 1, 5 + Num2_p2, A2_2, A2_2, A1_2, A1_2, A3_2, A5_2, A3_2, A4_2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1_p2 + Num2_p2 / 2); i<(Nx_2 + Nx_3 + Num1_p2 + Num2_p2); i++)// num 2 :  26 node! part2
	{
		switch (i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 / 2)
		{
		case 0:// 3 lei 2j shang
			kind3_node(j, i, -(Num1_p2 - Node13_p2 + Num2_p2 / 2) - 1, -(Num1_p2 - Node13_p2 + Num2_p2 / 2), 0, 1,
				Num3_p2 - 5 - 2, Num3_p2 - 5, C4, C3_2, C1_2, C2_2, C4, C3_2, Coef, Coef_location);
			break;
		case (Num2_p2 / 2 - 1) :// 3 lei 2j xia
			kind3_node(j, i, -(Num1_p2 - Node14_p2 + Num2_p2 - 1), -(Num1_p2 - Node14_p2 + Num2_p2 - 1) + 1, -1, 0,
			Num3_p2 - 5, Num3_p2 - 5 + 2, C3_2, C4, C2_2, C1_2, C3_2, C4, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 / 2) % 2 == 0)// 1dx
			{
				inner_5node(j, i, 1, -((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 / 2) + Num1_p2 - Node13_p2 - (i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 / 2) / 2) - Num2_p2 / 2,
					Num3_p2 - 5, Coef, Coef_location);
			}
			else//1 lei 2j left
			{
				kind1_node(j, i, -((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 / 2) + Num1_p2 - Node13_p2 - (i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 / 2) / 2) - Num2_p2 / 2 - (Node_11_p2 + 1),
					-((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 / 2) + Num1_p2 - Node13_p2 - (i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 / 2) / 2) - Num2_p2 / 2 - (Node_11_p2 + 1) + 1,
					-((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 / 2) + Num1_p2 - Node13_p2 - (i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 / 2) / 2) - Num2_p2 / 2,
					-((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 / 2) + Num1_p2 - Node13_p2 - (i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 / 2) / 2) - Num2_p2 / 2 + 1,
					-1, 0, 1, Num3_p2 - 5, A2_2, A2_2, A1_2, A1_2, A3_2, A5_2, A3_2, A4_2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!!!
	for (i = (Nx_2 + Nx_3 + Num1_p2 + Num2_p2); i<(Nx_2 + Nx_3 + Num1_p2 + Num2_p2 + Num3_p2 / 2); i++)// num 3 :  46 node! part1
	{
		switch (i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2)
		{
		case 0:// 1 lei 1j shang
			kind1_node(j, i, -(Num1_p2 - N_d0 + 2 + Num2_p2) - 2, -(Num1_p2 - N_d0 + 2 + Num2_p2) - 1, -(Num1_p2 - N_d0 + 2 + Num2_p2),
				0, 1, Num3_p2 + Num4_p2 + N_d0 - 2 - 2, Num3_p2 + Num4_p2 + N_d0 - 2 - 1, Num3_p2 + Num4_p2 + N_d0 - 2,
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case 3:// 2 lei 2j shang
			kind2_node(j, i, -(Num1_p2 - N_d0 + 2 + Num2_p2), -2, -1, 0, 1, Num3_p2 + Num4_p2 + N_d0 - 2,
				B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case 4:// 3 lei 2j shang
			kind3_node(j, i, -(Num1_p2 - N_d0 + 2 + Num2_p2) - 1, -(Num1_p2 - N_d0 + 2 + Num2_p2), -1, 0, 1, Num3_p2 - 4,
				C4, C4, C3_2, C1_2, C3_2, C2_2, Coef, Coef_location);
			break;
		case (Num3_p2 / 2 - 5) :// 3 lei 2j xia
			kind3_node(j, i, -(Num1_p2 - Node_1_p2 + Num2_p2 + Num3_p2 / 2 - 1), -(Num1_p2 - Node_1_p2 + Num2_p2 + Num3_p2 / 2 - 1) + 1,
			-1, 0, 1, Num3_p2 - 4, C4, C4, C3_2, C1_2, C3_2, C2_2, Coef, Coef_location);
			break;
		case (Num3_p2 / 2 - 4) :// 2 lei 2j xia
			kind2_node(j, i, -(Num1_p2 - Node_1_p2 + Num2_p2 + Num3_p2 / 2 - 1), -1, 0, 1, 2, Num3_p2 / 2 + 1 + Num4_p2 + Node_5_p2,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num3_p2 / 2 - 1) :// 1 lei 1j xia
			kind1_node(j, i, -(Num1_p2 - Node_1_p2 + Num2_p2 + Num3_p2 / 2 - 1), -(Num1_p2 - Node_1_p2 + Num2_p2 + Num3_p2 / 2 - 1) + 1,
			-(Num1_p2 - Node_1_p2 + Num2_p2 + Num3_p2 / 2 - 1) + 2, -1, 0,
			Num3_p2 / 2 + 1 + Num4_p2 + Node_5_p2, Num3_p2 / 2 + 1 + Num4_p2 + Node_5_p2 + 1, Num3_p2 / 2 + 1 + Num4_p2 + Node_5_p2 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2) < 3)// 2dx
			{
				inner_5node(j, i, 2, -(Num1_p2 - N_d0 + 2 + Num2_p2), Num3_p2 + Num4_p2 + N_d0 - 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2) > 4) && ((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2) < (Num3_p2 / 2 - 5)))// 1dx
			{
				inner_5node(j, i, 1, -(Num2_p2 + 5), Num3_p2 - 4, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, i, 2, -(Num1_p2 - Node_1_p2 + Num2_p2 + Num3_p2 / 2 - 1), Num3_p2 / 2 + 1 + Num4_p2 + Node_5_p2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1_p2 + Num2_p2 + Num3_p2 / 2); i<(Nx_2 + Nx_3 + Num1_p2 + Num2_p2 + Num3_p2); i++)// num 3 :  46 node! part2
	{
		switch (i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 / 2)
		{
		case 0:// 1 lei 1j shang
			kind1_node(j, i, -(Num1_p2 - Node_1_p2 - N_a + 3 + Num2_p2 + Num3_p2 / 2) - 2, -(Num1_p2 - Node_1_p2 - N_a + 3 + Num2_p2 + Num3_p2 / 2) - 1,
				-(Num1_p2 - Node_1_p2 - N_a + 3 + Num2_p2 + Num3_p2 / 2), 0, 1, Num3_p2 / 2 + Num4_p2 + Node_5_p2 + N_a - 3 - 2,
				Num3_p2 / 2 + Num4_p2 + Node_5_p2 + N_a - 3 - 1, Num3_p2 / 2 + Num4_p2 + Node_5_p2 + N_a - 3,
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case 3:// 2 lei 2j shang
			kind2_node(j, i, -(Num1_p2 - Node_1_p2 - N_a + 3 + Num2_p2 + Num3_p2 / 2), -2, -1, 0, 1,
				Num3_p2 / 2 + Num4_p2 + Node_5_p2 + N_a - 3, B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case 4:// 3 lei 2j shang
			kind3_node(j, i, -(Num1_p2 - Node_1_p2 - N_a + 3 + Num2_p2 + Num3_p2 / 2) - 1, -(Num1_p2 - Node_1_p2 - N_a + 3 + Num2_p2 + Num3_p2 / 2),
				-1, 0, 1, 5 + Num4_p2 - 1, C4, C4, C3_2, C1_2, C3_2, C2_2, Coef, Coef_location);
			break;
		case (Num3_p2 / 2 - 5) :// 3 lei 2j xia
			kind3_node(j, i, -(Num1_p2 - Node_11_p2 + Num2_p2 + Num3_p2 - 1), -(Num1_p2 - Node_11_p2 + Num2_p2 + Num3_p2 - 1) + 1,
			-1, 0, 1, 5 + Num4_p2 - 1, C4, C4, C3_2, C1_2, C3_2, C2_2, Coef, Coef_location);
			break;
		case (Num3_p2 / 2 - 4) :// 2 lei 2j xia
			kind2_node(j, i, -(Num1_p2 - Node_11_p2 + Num2_p2 + Num3_p2 - 1), -1, 0, 1, 2, 1 + Num4_p2 + Node_55_p2,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num3_p2 / 2 - 1) :// 1 lei 1j xia
			kind1_node(j, i, -(Num1_p2 - Node_11_p2 + Num2_p2 + Num3_p2 - 1), -(Num1_p2 - Node_11_p2 + Num2_p2 + Num3_p2 - 1) + 1,
			-(Num1_p2 - Node_11_p2 + Num2_p2 + Num3_p2 - 1) + 2, -1, 0,
			1 + Num4_p2 + Node_55_p2, 1 + Num4_p2 + Node_55_p2 + 1, 1 + Num4_p2 + Node_55_p2 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 / 2) < 3)// 2dx
			{
				inner_5node(j, i, 2, -(Num1_p2 - Node_1_p2 - N_a + 3 + Num2_p2 + Num3_p2 / 2), Num3_p2 / 2 + Num4_p2 + Node_5_p2 + N_a - 3, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 / 2) > 4) && ((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 / 2) < (Num3_p2 / 2 - 5)))// 1dx
			{
				inner_5node(j, i, 1, -(Num3_p2 - 5), 5 + Num4_p2 - 1, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, i, 2, -(Num1_p2 - Node_11_p2 + Num2_p2 + Num3_p2 - 1), 1 + Num4_p2 + Node_55_p2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (Nx_2 + Nx_3 + Num1_p2 + Num2_p2 + Num3_p2);
		i<(Nx_2 + Nx_3 + Num1_p2 + Num2_p2 + Num3_p2 + Num4_p2 / 2); i++)// num 4 :  30 node! part1
	{
		switch (i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2)
		{
		case 0:// 1 lei 2j shang
			kind1_node(j, i, -(Num3_p2 - 4) - 2, -(Num3_p2 - 4) - 1, -(Num3_p2 - 4), 0, 1,
				Num4_p2 + N_d0 + 2 - 2, Num4_p2 + N_d0 + 2 - 1, Num4_p2 + N_d0 + 2,
				A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case (Num4_p2 / 2 - 1) :// 1 lei 2j xia
			kind1_node(j, i, -(Num3_p2 - 4), -(Num3_p2 - 4) + 1, -(Num3_p2 - 4) + 2, -1, 0,
			Num4_p2 + N_d0 + 2 - N_c1, Num4_p2 + N_d0 + 2 - N_c1 + 1, Num4_p2 + N_d0 + 2 - N_c1 + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2) <= Node41_p2)// 1dx
			{
				inner_5node(j, i, 1, -(Num3_p2 - 4), Num4_p2 + N_d0 + 2, Coef, Coef_location);
			}
			else if ((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2) >= Node42_p2) // 1dx
			{
				inner_5node(j, i, 1, -(Num3_p2 - 4), Num4_p2 + N_d0 + 2 - N_c1, Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_4node(j, i, -(Num3_p2 - 4), -1, 0, 1, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1_p2 + Num2_p2 + Num3_p2 + Num4_p2 / 2);
		i<(Nx_2 + Nx_3 + Num1_p2 + Num2_p2 + Num3_p2 + Num4_p2); i++)// num 4 :  30 node! part2
	{
		switch (i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2)
		{
		case (Num4_p2 / 2) :// 1 lei 2j shang
			kind1_node(j, i, -(5 + Num4_p2 - 1) - 2, -(5 + Num4_p2 - 1) - 1, -(5 + Num4_p2 - 1), 0, 1,
			(Node_55_p2 - 4 + 1 + N_c1) - 2, (Node_55_p2 - 4 + 1 + N_c1) - 1, (Node_55_p2 - 4 + 1 + N_c1),
			A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case (Num4_p2 - 1) :// 1 lei 2j xia
			kind1_node(j, i, -(5 + Num4_p2 - 1), -(5 + Num4_p2 - 1) + 1, -(5 + Num4_p2 - 1) + 2, -1, 0,
			(Node_55_p2 - 4 + 1), (Node_55_p2 - 4 + 1) + 1, (Node_55_p2 - 4 + 1) + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2) <= Node43_p2)// 1dx
			{
				inner_5node(j, i, 1, -(5 + Num4_p2 - 1), (Node_55_p2 - 4 + 1 + N_c1), Coef, Coef_location);
			}
			else if ((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2) >= Node44_p2) // 1dx
			{
				inner_5node(j, i, 1, -(5 + Num4_p2 - 1), (Node_55_p2 - 4 + 1), Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_4node(j, i, -(5 + Num4_p2 - 1), -1, 0, 1, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (Nx_2 + Nx_3 + Num1_p2 + Num2_p2 + Num3_p2 + Num4_p2);
		i<(Nx_2 + Nx_3 + Num1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + Node_5_p2 + 2); i++)// circle num5 :  80 node!  part1 !
	{
		switch (i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2)
		{
		case 0:
			boundary_node(j, i, 1, -Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2, Num5_p2 + Num6_p2 + Num7_p2 + Num8_p2, Coef, Coef_location);
			break;
		case (N_d0 - 3) ://2 lei shang 1j
			kind2_node(j, i, -Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2, -2, -1, 0, 1, Num5_p2 + Num6_p2 + Num7_p2 + Num8_p2,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (N_d0 + 1) :// 2 lei 2j shang
			kind2_node(j, i, -(Num3_p2 + Num4_p2 + N_d0 - 2), -2, -1, 0, 1, Num5_p2 - N_d0 + 2 + Num6_p2,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node51_p2:// inner 3 node
			inner_3node(j, i, -(Num4_p2 + N_d0 + 2), -1, 0, 1, 1, k*k*dx*dx - 4, Coef, Coef_location);
			break;
		case (Node51_p2 + 1) :// inner 3 node
			inner_3node(j, i, -(Num4_p2 + N_d0 + 2) + N_c1, 0, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node_5_p2 - 3) ://2 lei xia 2j
			kind2_node(j, i, -(Num3_p2 / 2 + 1 + Num4_p2 + Node_5_p2), -1, 0, 1, 2, Num5_p2 - Node_5_p2 + Num6_p2 + Node72_p2,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_5_p2 + 1) ://2 lei xia 1j
			kind2_node(j, i, -(Num1_p2 - Node_1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + Node_5_p2), -1, 0, 1, 2,
			Num5_p2 - Node_5_p2 + Num6_p2 + Num7_p2 + Num8_p2 + Node_9_p2,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2) < (N_d0 - 3))// 4dx!
			{
				inner_5node(j, i, 4, -Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2, Num5_p2 + Num6_p2 + Num7_p2 + Num8_p2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2) > (N_d0 - 3)) && ((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2) < (N_d0 + 1)))// 2dx!
			{
				inner_5node(j, i, 2, -(Num3_p2 + Num4_p2 + N_d0 - 2), Num5_p2 - N_d0 + 2 + Num6_p2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2) > (N_d0 + 1)) && ((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2) < Node51_p2))// 1dx
			{
				inner_5node(j, i, 1, -(Num4_p2 + N_d0 + 2), Num5_p2 - N_d0 - 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2) > (Node51_p2 + 1)) && ((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2) < (Node_5_p2 - 3)))// 1dx
			{
				inner_5node(j, i, 1, -(Num4_p2 + N_d0 + 2) + N_c1, Num5_p2 - N_d0 - 2 - 2, Coef, Coef_location);
			}
			else// 2dx!
			{
				inner_5node(j, i, 2, -(Num3_p2 / 2 + 1 + Num4_p2 + Node_5_p2), Num5_p2 - Node_5_p2 + Num6_p2 + Node72_p2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + Node_5_p2 + 2);
		i<(Nx_2 + Nx_3 + Num1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + Node53_p2 + 3); i++)// circle num5 :  80 node!  part2 !
	{
		switch (i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2)
		{
		default:
			if ((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2) < Node52_p2)// inner 5 node 4dx
			{
				inner_5node(j, i, 4, -(Num1_p2 - Node_1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + Node_5_p2), Num5_p2 - Node_5_p2 + Num6_p2 + Num7_p2 + Num8_p2 + Node_9_p2, Coef, Coef_location);
			}
			else if ((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2) > Node53_p2)// inner 5 node 4dx
			{
				inner_5node(j, i, 4, -(Num1_p2 - Node_1_p2 - 9 + Num2_p2 + Num3_p2 + Num4_p2 + Node53_p2),
					Num5_p2 - Node53_p2 + Num6_p2 + Num7_p2 + Num8_p2 + Node93_p2, Coef, Coef_location);
			}
			else// left  2lei 1j
			{
				kind2_node(j, i, -Nx_3 - Nx_2 - (Num1_p2 - Node_1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + Node_5_p2),
					-(Num1_p2 - Node_1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + Node_5_p2), -1, 0, 1,
					Num5_p2 - (i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2) + Num6_p2 + Node72_p2 + 1 + (i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Node52_p2) * 2,
					B5, B4, B3, B1_1, B3, B2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//**********************************************************************************!!!!!!!!!!!!
	for (i = (Nx_2 + Nx_3 + Num1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + Node53_p2 + 3);
		i<(Nx_2 + Nx_3 + Num1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + Num5_p2); i++)// circle num5 :  80 node!  part3 !
	{
		switch (i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2)
		{
		case (Num5_p2 - 1) :
			boundary_node(j, i, 2, -Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2, Num6_p2 + Num7_p2 + Num8_p2 + Num9_p2, Coef, Coef_location);
			break;
		case (Node53_p2 + 3) ://2 lei shang 1j
			kind2_node(j, i, -(Num1_p2 - Node_1_p2 - 9 + Num2_p2 + Num3_p2 + Num4_p2 + Node53_p2), -2, -1, 0, 1, Num5_p2 - Node53_p2 + Num6_p2 + Num7_p2 + Num8_p2 + Node93_p2,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node53_p2 + 7) :// 2 lei 2j shang
			kind2_node(j, i, -(Num3_p2 / 2 + Num4_p2 + Node53_p2 + 4), -2, -1, 0, 1, Num5_p2 - Node53_p2 - 4 + Num6_p2 + Node73_p2,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node54_p2:// inner 3 node
			inner_3node(j, i, -(Node_55_p2 - 4 + 1) - N_c1, -1, 0, 1, 1, k*k*dx*dx - 4, Coef, Coef_location);
			break;
		case (Node54_p2 + 1) :// inner 3 node
			inner_3node(j, i, -(Node_55_p2 - 4 + 1), 0, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node_55_p2 - 3) ://2 lei xia 2j
			kind2_node(j, i, -(Num4_p2 + Node_55_p2 + 1), -1, 0, 1, 2, Num5_p2 - Node_55_p2 + Num6_p2 + Num7_p2 - 1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_55_p2 + 1) ://2 lei xia 1j
			kind2_node(j, i, -Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2, -1, 0, 1, 2, Num6_p2 + Num7_p2 + Num8_p2 + Num9_p2,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2) < (Node53_p2 + 7))// 2dx!
			{
				inner_5node(j, i, 2, -(Num3_p2 / 2 + Num4_p2 + Node53_p2 + 4), Num5_p2 - Node53_p2 - 4 + Num6_p2 + Node73_p2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2) > (Node53_p2 + 7)) && ((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2) < Node54_p2))// 1dx
			{
				inner_5node(j, i, 1, -(Node_55_p2 - 4 + 1) - N_c1, (Num5_p2 - Node_55_p2 + 4 + Num6_p2 - 1) + 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2) > (Node54_p2 + 1)) && ((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2) < (Node_55_p2 - 3)))// 1dx
			{
				inner_5node(j, i, 1, -(Node_55_p2 - 4 + 1), (Num5_p2 - Node_55_p2 + 4 + Num6_p2 - 1), Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2) > (Node_55_p2 - 3)) && ((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2) < (Node_55_p2 + 1)))// 2dx
			{
				inner_5node(j, i, 2, -(Num4_p2 + Node_55_p2 + 1), Num5_p2 - Node_55_p2 + Num6_p2 + Num7_p2 - 1, Coef, Coef_location);
			}
			else// 4dx!
			{
				inner_5node(j, i, 4, -Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2, Num6_p2 + Num7_p2 + Num8_p2 + Num9_p2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (Nx_2 + Nx_3 + Num1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + Num5_p2);
		i<(Nx_2 + Nx_3 + Num1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + Num5_p2 + Num6_p2); i++)// num 6 :  16 node!
	{
		switch (i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2)
		{
		case 0:// 1 lei 2j shang
			kind1_node(j, i, -(Num5_p2 - N_d0 - 2) - 2, -(Num5_p2 - N_d0 - 2) - 1, -(Num5_p2 - N_d0 - 2), 0, 1, Num6_p2 + 4 - 2,
				Num6_p2 + 4 - 1, Num6_p2 + 4, A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case Node61_p2:// inner 3 node
			inner_3node(j, i, -(Num5_p2 - N_d0 - 2), -1, 0, 1, 1, k*k*dx*dx - 4, Coef, Coef_location);
			break;
		case (Node61_p2 + 1) :// inner 3 node
			inner_3node(j, i, -(Num5_p2 - N_d0 - 2 - 2), 0, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Num6_p2 / 2 - 1) :// 1 lei 2j xia
			kind1_node(j, i, -(Num5_p2 - N_d0 - 2 - 2), -(Num5_p2 - N_d0 - 2 - 2) + 1, -(Num5_p2 - N_d0 - 2 - 2) + 2, -1, 0,
			Num6_p2 + 4 - 2, Num6_p2 + 4 - 2 + 1, Num6_p2 + 4 - 2 + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		case (Num6_p2 / 2) :// 1 lei 2j shang
			kind1_node(j, i, -(Num5_p2 - Node_55_p2 + 4 + Num6_p2 - 1) - 2 - 2, -(Num5_p2 - Node_55_p2 + 4 + Num6_p2 - 1) - 2 - 1,
			-(Num5_p2 - Node_55_p2 + 4 + Num6_p2 - 1) - 2, 0, 1, Num7_p2 - 4 + 2 - 2, Num7_p2 - 4 + 2 - 1, Num7_p2 - 4 + 2,
			A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case Node62_p2:// inner 3 node
			inner_3node(j, i, -(Num5_p2 - Node_55_p2 + 4 + Num6_p2 - 1) - 2, -1, 0, 1, 1, k*k*dx*dx - 4, Coef, Coef_location);
			break;
		case (Node62_p2 + 1) :// inner 3 node
			inner_3node(j, i, -(Num5_p2 - Node_55_p2 + 4 + Num6_p2 - 1), 0, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Num6_p2 - 1) :// 1 lei 2j xia
			kind1_node(j, i, -(Num5_p2 - Node_55_p2 + 4 + Num6_p2 - 1), -(Num5_p2 - Node_55_p2 + 4 + Num6_p2 - 1) + 1,
			-(Num5_p2 - Node_55_p2 + 4 + Num6_p2 - 1) + 2, -1, 0, Num7_p2 - 4, Num7_p2 - 4 + 1, Num7_p2 - 4 + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2) < Node61_p2)// 1dx
			{
				inner_5node(j, i, 1, -(Num5_p2 - N_d0 - 2), Num6_p2 + 4, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2) > (Node61_p2 + 1)) && ((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2) < (Num6_p2 / 2 - 1)))// 2dx
			{
				inner_5node(j, i, 1, -(Num5_p2 - N_d0 - 2 - 2), Num6_p2 + 4 - 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2) > (Num6_p2 / 2)) && ((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2) < Node62_p2))// 2dx
			{
				inner_5node(j, i, 1, -(Num5_p2 - Node_55_p2 + 4 + Num6_p2 - 1) - 2, Num7_p2 - 4 + 2, Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_5node(j, i, 1, -(Num5_p2 - Node_55_p2 + 4 + Num6_p2 - 1), Num7_p2 - 4, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!!
	for (i = (Nx_2 + Nx_3 + Num1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + Num5_p2 + Num6_p2);
		i<(Nx_2 + Nx_3 + Num1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + Num5_p2 + Num6_p2 + Node72_p2 + 1); i++)// num 7 :  39 node! part1
	{
		switch (i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2 - Num6_p2)
		{
		case 0:// 1 lei 1j shang
			kind1_node(j, i, -(Num5_p2 - N_d0 + 2 + Num6_p2) - 2, -(Num5_p2 - N_d0 + 2 + Num6_p2) - 1, -(Num5_p2 - N_d0 + 2 + Num6_p2),
				0, 1, Num7_p2 + Num8_p2 + N_d0 - 2 - 2, Num7_p2 + Num8_p2 + N_d0 - 2 - 1, Num7_p2 + Num8_p2 + N_d0 - 2,
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case 3:// 2 lei 2j shang
			kind2_node(j, i, -(Num5_p2 - N_d0 + 2 + Num6_p2), -2, -1, 0, 1, Num7_p2 + Num8_p2 + N_d0 - 2,
				B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node71_p2://inner 4 node
			inner_4node(j, i, -(Num6_p2 + 4), -1, 0, Num7_p2 - 4, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node71_p2 + 1) ://inner 4 node
			inner_4node(j, i, -(Num6_p2 + 4) + 2, 0, 1, Num7_p2 - 4, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Node72_p2 - 3) :// 2 lei 2j xia
			kind2_node(j, i, -(Num5_p2 - Node_5_p2 + Num6_p2 + Node72_p2), -1, 0, 1, 2, Num7_p2 - Node72_p2 + Num8_p2 + Node_9_p2,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case Node72_p2:// 1 lei 1j xia
			kind1_node(j, i, -(Num5_p2 - Node_5_p2 + Num6_p2 + Node72_p2), -(Num5_p2 - Node_5_p2 + Num6_p2 + Node72_p2) + 1,
				-(Num5_p2 - Node_5_p2 + Num6_p2 + Node72_p2) + 2, -1, 0, Num7_p2 - Node72_p2 + Num8_p2 + Node_9_p2,
				Num7_p2 - Node72_p2 + Num8_p2 + Node_9_p2 + 1, Num7_p2 - Node72_p2 + Num8_p2 + Node_9_p2 + 2,
				A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2 - Num6_p2) < 3)// 2dx
			{
				inner_5node(j, i, 2, -(Num5_p2 - N_d0 + 2 + Num6_p2), Num7_p2 + Num8_p2 + N_d0 - 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2 - Num6_p2) > 3) && ((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2 - Num6_p2) < (Node71_p2)))// 1dx
			{
				inner_5node(j, i, 1, -(Num6_p2 + 4), Num7_p2 - 4, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2 - Num6_p2) > (Node71_p2 + 1)) && ((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2 - Num6_p2) < (Node72_p2 - 3)))// 1dx
			{
				inner_5node(j, i, 1, -(Num6_p2 + 4) + 2, Num7_p2 - 4, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, i, 2, -(Num5_p2 - Node_5_p2 + Num6_p2 + Node72_p2), Num7_p2 - Node72_p2 + Num8_p2 + Node_9_p2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + Num5_p2 + Num6_p2 + Node72_p2 + 1);
		i<(Nx_2 + Nx_3 + Num1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + Num5_p2 + Num6_p2 + Node73_p2); i++)// num 7 :  39 node! part2
	{
		switch (i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2 - Num6_p2 - Node72_p2 - 1)
		{
		case 0:// 3 lei 1j shang
			kind3_node(j, i, -(Num5_p2 - Node52_p2 + Num6_p2 + Node72_p2 + 1) - 1, -(Num5_p2 - Node52_p2 + Num6_p2 + Node72_p2 + 1),
				0, 1, (Num7_p2 - Node72_p2 + Num8_p2 + Node92_p2) - 2, (Num7_p2 - Node72_p2 + Num8_p2 + Node92_p2),
				C4, C3_1, C1_1, C2_1, C4, C3_1, Coef, Coef_location);
			break;
		case 10:// 3 lei 1j xia
			kind3_node(j, i, -(Num5_p2 - Node53_p2 + Num6_p2 + Node73_p2 - 1), -(Num5_p2 - Node53_p2 + Num6_p2 + Node73_p2 - 1) + 1,
				-1, 0, (Num7_p2 - Node72_p2 + Num8_p2 + Node92_p2), (Num7_p2 - Node72_p2 + Num8_p2 + Node92_p2) + 2,
				C3_1, C4, C2_1, C1_1, C3_1, C4, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2 - Num6_p2 - Node72_p2 - 1) % 2 == 0)// 1dx
			{
				inner_5node(j, i, 2, -((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2 - Num6_p2) + Num6_p2 + Num5_p2 - Node52_p2 - (i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2 - Num6_p2 - Node72_p2 - 1) / 2),
					(Num7_p2 - Node72_p2 + Num8_p2 + Node92_p2), Coef, Coef_location);
			}
			else//1 lei  1j left
			{
				kind1_node(j, i, -((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2 - Num6_p2) + Num6_p2 + Num5_p2 - Node52_p2 - (i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2 - Num6_p2 - Node72_p2 - 1) / 2) - (Num1_p2 - Node_1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + Node_5_p2),
					-((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2 - Num6_p2) + Num6_p2 + Num5_p2 - Node52_p2 - (i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2 - Num6_p2 - Node72_p2 - 1) / 2) - (Num1_p2 - Node_1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + Node_5_p2) + 1,
					-((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2 - Num6_p2) + Num6_p2 + Num5_p2 - Node52_p2 - (i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2 - Num6_p2 - Node72_p2 - 1) / 2),
					-((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2 - Num6_p2) + Num6_p2 + Num5_p2 - Node52_p2 - (i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2 - Num6_p2 - Node72_p2 - 1) / 2) + 1,
					-1, 0, 1, (Num7_p2 - Node72_p2 + Num8_p2 + Node92_p2), A2_1, A2_1, A1_1, A1_1, A3_1, A5_1, A3_1, A4_1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + Num5_p2 + Num6_p2 + Node73_p2);
		i<(Nx_2 + Nx_3 + Num1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + Num5_p2 + Num6_p2 + Num7_p2); i++)// num 7 :  39 node! part3
	{
		switch (i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2 - Num6_p2)
		{
		case Node73_p2:// 1 lei 1j shang
			kind1_node(j, i, -(Num5_p2 - Node53_p2 - 4 + Num6_p2 + Node73_p2) - 2, -(Num5_p2 - Node53_p2 - 4 + Num6_p2 + Node73_p2) - 1,
				-(Num5_p2 - Node53_p2 - 4 + Num6_p2 + Node73_p2), 0, 1, Num7_p2 - Node73_p2 + Num8_p2 + Node93_p2 + 4 - 2,
				Num7_p2 - Node73_p2 + Num8_p2 + Node93_p2 + 4 - 1, Num7_p2 - Node73_p2 + Num8_p2 + Node93_p2 + 4,
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case (Node73_p2 + 3) :// 2 lei 2j shang
			kind2_node(j, i, -(Num5_p2 - Node53_p2 - 4 + Num6_p2 + Node73_p2), -2, -1, 0, 1,
			Num7_p2 - Node73_p2 + Num8_p2 + Node93_p2 + 4, B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node74_p2://inner 4 node
			inner_4node(j, i, -(Num7_p2 - 5 + 1) - 2, -1, 0, 5 + Num8_p2 - 1, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node74_p2 + 1) ://inner 4 node
			inner_4node(j, i, -(Num7_p2 - 5 + 1), 0, 1, 5 + Num8_p2 - 1, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Num7_p2 - 4) :// 2 lei 2j xia
			kind2_node(j, i, -(Num5_p2 - Node_55_p2 + Num6_p2 + Num7_p2 - 1), -1, 0, 1, 2, 1 + Num8_p2 + Node_99_p2,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num7_p2 - 1) :// 1 lei 1j xia
			kind1_node(j, i, -(Num5_p2 - Node_55_p2 + Num6_p2 + Num7_p2 - 1), -(Num5_p2 - Node_55_p2 + Num6_p2 + Num7_p2 - 1) + 1,
			-(Num5_p2 - Node_55_p2 + Num6_p2 + Num7_p2 - 1) + 2, -1, 0,
			1 + Num8_p2 + Node_99_p2, 1 + Num8_p2 + Node_99_p2 + 1, 1 + Num8_p2 + Node_99_p2 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2 - Num6_p2) < (Node73_p2 + 3))// 2dx
			{
				inner_5node(j, i, 2, -(Num5_p2 - Node53_p2 - 4 + Num6_p2 + Node73_p2), Num7_p2 - Node73_p2 + Num8_p2 + Node93_p2 + 4, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2 - Num6_p2) > (Node73_p2 + 3)) && ((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2 - Num6_p2) < (Node74_p2)))// 1dx
			{
				inner_5node(j, i, 1, -(Num7_p2 - 5 + 1) - 2, 5 + Num8_p2 - 1, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2 - Num6_p2) > (Node74_p2 + 1)) && ((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2 - Num6_p2) < (Num7_p2 - 4)))// 1dx
			{
				inner_5node(j, i, 1, -(Num7_p2 - 5 + 1), 5 + Num8_p2 - 1, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, i, 2, -(Num5_p2 - Node_55_p2 + Num6_p2 + Num7_p2 - 1), 1 + Num8_p2 + Node_99_p2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//**************************************************************
	cout << "middle_7区域节点数:" << j << endl;
	if (j == sizeM)
		cout << "middle_7 passed..." << endl;
	else
		cout << "middle_7 failed..." << endl;
	//修正矩阵，将不在本区域的节点刨除
	for (int i = 0; i != sizeM; i++)
		for (int j = 0; j != N_matrix; j++){
		if (Coef_location[i][j] < 0 || Coef_location[i][j] >= sizeN)
			Coef_location[i][j] = -1;
		}
	if (Debug)
	{
		ofstream outFile_Coef("middle_7.txt");
		if (!outFile_Coef)
		{
			cout << "Can't write the coefficient matrix into file !" << endl;
			exit(1);//Terminate with fault
		}
		else
		{
			int i, j;
			for (i = 0; i < sizeM; i++)
			{
				for (j = 0; j < N_matrix; j++)
				{
					if ((Coef[i][j].real != 0 || Coef[i][j].image != 0) && Coef_location[i][j] >= 0)
						outFile_Coef << Coef[i][j].real << "+j*" << Coef[i][j].image << " (" << i << "," << Coef_location[i][j] << ") ";
				}
				outFile_Coef << endl << endl;
			}
		}
		outFile_Coef.close();
	}

};
void scan_coef_Matrix_middle_5(complex** Coef, int num_total, int ** Coef_location, int sizeM, int sizeN, int offset)
{
	int i = 0, j = 0, m = 0;
	for (j = 0; j<num_total; j++)         //系数矩阵初始置零
	{
		for (i = 0; i<N_matrix; i++)
		{
			Coef[j][i].real = 0;
			Coef[j][i].image = 0;
			Coef_location[j][i] = -1;//初始化为负值
		}
	}
	j = 0;
	//******************************************************************************************!!! num 31 !!!
	for (i = (M_post2_q3 + Num8_p2); i<(M_post2_q3 + Num8_p2 + Node72_p2 + 1); i++)// num 31(7) :  39 node! part1
	{
		switch (i - M_post2_q3 - Num8_p2)
		{
		case 0:// 1 lei 1j shang
			kind1_node(j, i, -(Num9_p2 - N_d0 + 2 + Num8_p2) - 2, -(Num9_p2 - N_d0 + 2 + Num8_p2) - 1, -(Num9_p2 - N_d0 + 2 + Num8_p2),
				0, 1, Num7_p2 + Num6_p2 + N_d0 - 2 - 2, Num7_p2 + Num6_p2 + N_d0 - 2 - 1, Num7_p2 + Num6_p2 + N_d0 - 2,
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case 3:// 2 lei 2j shang
			kind2_node(j, i, -(Num9_p2 - N_d0 + 2 + Num8_p2), -2, -1, 0, 1, Num7_p2 + Num6_p2 + N_d0 - 2,
				B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node71_p2://inner 4 node
			inner_4node(j, i, -(Num8_p2 + 4), -1, 0, Num7_p2 - 4, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node71_p2 + 1) ://inner 4 node
			inner_4node(j, i, -(Num8_p2 + 4), 0, 1, Num7_p2 - 4 + 2, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Node72_p2 - 3) :// 2 lei 2j xia
			kind2_node(j, i, -(Num9_p2 - Node_9_p2 + Num8_p2 + Node72_p2), -1, 0, 1, 2, (Num7_p2 - Node72_p2 + Num6_p2 + Node_5_p2),
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case Node72_p2:// 1 lei 1j xia
			kind1_node(j, i, -(Num9_p2 - Node_9_p2 + Num8_p2 + Node72_p2), -(Num9_p2 - Node_9_p2 + Num8_p2 + Node72_p2) + 1,
				-(Num9_p2 - Node_9_p2 + Num8_p2 + Node72_p2) + 2, -1, 0, (Num7_p2 - Node72_p2 + Num6_p2 + Node_5_p2),
				(Num7_p2 - Node72_p2 + Num6_p2 + Node_5_p2) + 1, (Num7_p2 - Node72_p2 + Num6_p2 + Node_5_p2) + 2,
				A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q3 - Num8_p2) < 3)// 2dx
			{
				inner_5node(j, i, 2, -(Num9_p2 - N_d0 + 2 + Num8_p2), Num7_p2 + Num6_p2 + N_d0 - 2, Coef, Coef_location);
			}
			else if (((i - M_post2_q3 - Num8_p2) > 3) && ((i - M_post2_q3 - Num8_p2) < (Node71_p2)))// 1dx
			{
				inner_5node(j, i, 1, -(Num8_p2 + 4), Num7_p2 - 4, Coef, Coef_location);
			}
			else if (((i - M_post2_q3 - Num8_p2) > (Node71_p2 + 1)) && ((i - M_post2_q3 - Num8_p2) < (Node72_p2 - 3)))// 1dx
			{
				inner_5node(j, i, 1, -(Num8_p2 + 4), Num7_p2 - 4 + 2, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, i, 2, -(Num9_p2 - Node_9_p2 + Num8_p2 + Node72_p2), (Num7_p2 - Node72_p2 + Num6_p2 + Node_5_p2), Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post2_q3 + Num8_p2 + Node72_p2 + 1); i<(M_post2_q3 + Num8_p2 + Node73_p2); i++)// num 31(7) :  39 node! part2
	{
		switch (i - M_post2_q3 - Num8_p2 - Node72_p2 - 1)
		{
		case 0:// 3 lei 1j shang
			kind3_node(j, i, -(Num9_p2 - Node92_p2 + Num8_p2 + Node72_p2) - 2, -(Num9_p2 - Node92_p2 + Num8_p2 + Node72_p2),
				0, 1, (Num7_p2 - Node72_p2 - 1 + Num6_p2 + Node52_p2) - 1, (Num7_p2 - Node72_p2 - 1 + Num6_p2 + Node52_p2),
				C4, C3_1, C1_1, C2_1, C4, C3_1, Coef, Coef_location);
			break;
		case 10:// 3 lei 1j xia
			kind3_node(j, i, -(Num9_p2 - Node92_p2 + Num8_p2 + Node72_p2), -(Num9_p2 - Node92_p2 + Num8_p2 + Node72_p2) + 2,
				-1, 0, (Num7_p2 - Node73_p2 + 1 + Num6_p2 + Node53_p2), (Num7_p2 - Node73_p2 + 1 + Num6_p2 + Node53_p2) + 1,
				C3_1, C4, C2_1, C1_1, C3_1, C4, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q3 - Num8_p2 - Node72_p2 - 1) % 2 == 0)// 1dx
			{
				inner_5node(j, i, 2, -(Num9_p2 - Node92_p2 + Num8_p2 + Node72_p2),
					Num7_p2 - (i - M_post2_q3 - Num8_p2) + Num6_p2 + Node52_p2 + (i - M_post2_q3 - Num8_p2 - Node72_p2 - 1) / 2,
					Coef, Coef_location);
			}
			else//1 lei  1j left
			{
				kind1_node(j, i, -(Num9_p2 - Node92_p2 + Num8_p2 + Node72_p2), -1, 0, 1,
					Num7_p2 - (i - M_post2_q3 - Num8_p2) + Num6_p2 + Node52_p2 + (i - M_post2_q3 - Num8_p2 - Node72_p2 - 1) / 2,
					Num7_p2 - (i - M_post2_q3 - Num8_p2) + Num6_p2 + Node52_p2 + (i - M_post2_q3 - Num8_p2 - Node72_p2 - 1) / 2 + 1,
					Num7_p2 - (i - M_post2_q3 - Num8_p2) + Num6_p2 + Node52_p2 + (i - M_post2_q3 - Num8_p2 - Node72_p2 - 1) / 2 + (Num5_p2 - Node_5_p2 + Num4_p2 + Num3_p2 + Num2_p2 + Node_1_p2),
					Num7_p2 - (i - M_post2_q3 - Num8_p2) + Num6_p2 + Node52_p2 + (i - M_post2_q3 - Num8_p2 - Node72_p2 - 1) / 2 + (Num5_p2 - Node_5_p2 + Num4_p2 + Num3_p2 + Num2_p2 + Node_1_p2) + 1,
					A4_1, A3_1, A5_1, A3_1, A1_1, A1_1, A2_1, A2_1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post2_q3 + Num8_p2 + Node73_p2); i<(M_post2_q3 + Num8_p2 + Num7_p2); i++)// num 31(7) :  39 node! part3
	{
		switch (i - M_post2_q3 - Num8_p2)
		{
		case Node73_p2:// 1 lei 1j shang
			kind1_node(j, i, -(Num9_p2 - Node93_p2 - 4 + Num8_p2 + Node73_p2) - 2, -(Num9_p2 - Node93_p2 - 4 + Num8_p2 + Node73_p2) - 1,
				-(Num9_p2 - Node93_p2 - 4 + Num8_p2 + Node73_p2), 0, 1, (Num7_p2 - Node73_p2 + Num6_p2 + Node53_p2 + 4) - 2,
				(Num7_p2 - Node73_p2 + Num6_p2 + Node53_p2 + 4) - 1, (Num7_p2 - Node73_p2 + Num6_p2 + Node53_p2 + 4),
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case (Node73_p2 + 3) :// 2 lei 2j shang
			kind2_node(j, i, -(Num9_p2 - Node93_p2 - 4 + Num8_p2 + Node73_p2), -2, -1, 0, 1,
			(Num7_p2 - Node73_p2 + Num6_p2 + Node53_p2 + 4), B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node74_p2://inner 4 node
			inner_4node(j, i, -(Num7_p2 - 5 + 1), -1, 0, (5 + Num6_p2 - 1) - 2, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node74_p2 + 1) ://inner 4 node
			inner_4node(j, i, -(Num7_p2 - 5 + 1), 0, 1, 5 + Num6_p2 - 1, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Num7_p2 - 4) :// 2 lei 2j xia
			kind2_node(j, i, -(Num9_p2 - Node_99_p2 + Num8_p2 + Num7_p2 - 1), -1, 0, 1, 2, (1 + Num6_p2 + Node_55_p2),
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num7_p2 - 1) :// 1 lei 1j xia
			kind1_node(j, i, -(Num9_p2 - Node_99_p2 + Num8_p2 + Num7_p2 - 1), -(Num9_p2 - Node_99_p2 + Num8_p2 + Num7_p2 - 1) + 1,
			-(Num9_p2 - Node_99_p2 + Num8_p2 + Num7_p2 - 1) + 2, -1, 0,
			(1 + Num6_p2 + Node_55_p2), (1 + Num6_p2 + Node_55_p2) + 1, (1 + Num6_p2 + Node_55_p2) + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q3 - Num8_p2) < (Node73_p2 + 3))// 2dx
			{
				inner_5node(j, i, 2, -(Num9_p2 - Node93_p2 - 4 + Num8_p2 + Node73_p2), (Num7_p2 - Node73_p2 + Num6_p2 + Node53_p2 + 4), Coef, Coef_location);
			}
			else if (((i - M_post2_q3 - Num8_p2) > (Node73_p2 + 3)) && ((i - M_post2_q3 - Num8_p2) < (Node74_p2)))// 1dx
			{
				inner_5node(j, i, 1, -(Num7_p2 - 5 + 1), (5 + Num6_p2 - 1) - 2, Coef, Coef_location);
			}
			else if (((i - M_post2_q3 - Num8_p2) > (Node74_p2 + 1)) && ((i - M_post2_q3 - Num8_p2) < (Num7_p2 - 4)))// 1dx
			{
				inner_5node(j, i, 1, -(Num7_p2 - 5 + 1), 5 + Num6_p2 - 1, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, i, 2, -(Num9_p2 - Node_99_p2 + Num8_p2 + Num7_p2 - 1), (1 + Num6_p2 + Node_55_p2), Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!! num 32 !!!
	for (i = (M_post2_q3 + Num8_p2 + Num7_p2); i<(M_post2_q3 + Num8_p2 + Num7_p2 + Num6_p2); i++)// num 32(6) :  16 node!
	{
		switch (i - M_post2_q3 - Num8_p2 - Num7_p2)
		{
		case 0:// 1 lei 2j shang
			kind1_node(j, i, -(Num7_p2 - 4) - 2, -(Num7_p2 - 4) - 1, -(Num7_p2 - 4), 0, 1, (Num6_p2 + N_d0 + 2) - 2,
				(Num6_p2 + N_d0 + 2) - 1, (Num6_p2 + N_d0 + 2), A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case Node61_p2:// inner 3 node
			inner_3node(j, i, -1, 0, (Num6_p2 + N_d0 + 2), 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node61_p2 + 1) :// inner 3 node
			inner_3node(j, i, 0, 1, (Num6_p2 + N_d0 + 2) + 2, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Num6_p2 / 2 - 1) :// 1 lei 2j xia
			kind1_node(j, i, -(Num7_p2 - 4) - 2, -(Num7_p2 - 4) - 2 + 1, -(Num7_p2 - 4) - 2 + 2, -1, 0,
			(Num6_p2 + N_d0 + 2) + 2, (Num6_p2 + N_d0 + 2) + 2 + 1, (Num6_p2 + N_d0 + 2) + 2 + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		case (Num6_p2 / 2) :// 1 lei 2j shang
			kind1_node(j, i, -(5 + Num6_p2 - 1) + 2 - 2, -(5 + Num6_p2 - 1) + 2 - 1, -(5 + Num6_p2 - 1) + 2,
			0, 1, (1 + Node_55_p2 - 4) - 2 - 2, (1 + Node_55_p2 - 4) - 2 - 1, (1 + Node_55_p2 - 4) - 2,
			A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case Node62_p2:// inner 3 node
			inner_3node(j, i, -1, 0, (1 + Node_55_p2 - 4) - 2, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node62_p2 + 1) :// inner 3 node
			inner_3node(j, i, 0, 1, (1 + Node_55_p2 - 4), k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Num6_p2 - 1) :// 1 lei 2j xia
			kind1_node(j, i, -(5 + Num6_p2 - 1), -(5 + Num6_p2 - 1) + 1, -(5 + Num6_p2 - 1) + 2, -1, 0,
			(1 + Node_55_p2 - 4), (1 + Node_55_p2 - 4) + 1, (1 + Node_55_p2 - 4) + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q3 - Num8_p2 - Num7_p2) < Node61_p2)// 1dx
			{
				inner_5node(j, i, 1, -(Num7_p2 - 4), (Num6_p2 + N_d0 + 2), Coef, Coef_location);
			}
			else if (((i - M_post2_q3 - Num8_p2 - Num7_p2) > (Node61_p2 + 1)) && ((i - M_post2_q3 - Num8_p2 - Num7_p2) < (Num6_p2 / 2 - 1)))// 2dx
			{
				inner_5node(j, i, 1, -(Num7_p2 - 4) - 2, (Num6_p2 + N_d0 + 2) + 2, Coef, Coef_location);
			}
			else if (((i - M_post2_q3 - Num8_p2 - Num7_p2) > (Num6_p2 / 2)) && ((i - M_post2_q3 - Num8_p2 - Num7_p2) < Node62_p2))// 2dx
			{
				inner_5node(j, i, 1, -(5 + Num6_p2 - 1) + 2, (1 + Node_55_p2 - 4) - 2, Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_5node(j, i, 1, -(5 + Num6_p2 - 1), (1 + Node_55_p2 - 4), Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!! num 33 !!!
	for (i = (M_post2_q3 + Num8_p2 + Num7_p2 + Num6_p2);
		i<(M_post2_q3 + Num8_p2 + Num7_p2 + Num6_p2 + Node_5_p2 + 2); i++)// circle num 33(5) :  80 node!  part1 !
	{
		switch (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2)
		{
		case 0:
			boundary_node(j, i, 1, -Num9_p2 - Num8_p2 - Num7_p2 - Num6_p2, Num5_p2 + Num4_p2 + Num3_p2 + Num2_p2, Coef, Coef_location);
			break;
		case (N_d0 - 3) ://2 lei shang 1j
			kind2_node(j, i, -Num9_p2 - Num8_p2 - Num7_p2 - Num6_p2, -2, -1, 0, 1, Num5_p2 + Num4_p2 + Num3_p2 + Num2_p2,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (N_d0 + 1) :// 2 lei 2j shang
			kind2_node(j, i, -(Num7_p2 + Num6_p2 + N_d0 - 2), -2, -1, 0, 1, Num5_p2 - N_d0 + 2 + Num4_p2,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node51_p2:// inner 3 node
			inner_3node(j, i, -1, 0, Num5_p2 - N_d0 - 2, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node51_p2 + 1) :// inner 3 node
			inner_3node(j, i, 0, 1, Num5_p2 - N_d0 - 2 + N_c1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Node_5_p2 - 3) ://2 lei xia 2j
			kind2_node(j, i, -(Num7_p2 - Node72_p2 + Num6_p2 + Node_5_p2), -1, 0, 1, 2, Num5_p2 - Node_5_p2 + Num4_p2 + Num3_p2 / 2 - 1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_5_p2 + 1) ://2 lei xia 1j
			kind2_node(j, i, -(Num9_p2 - Node_9_p2 + Num8_p2 + Num7_p2 + Num6_p2 + Node_5_p2), -1, 0, 1, 2,
			Num5_p2 - Node_5_p2 + Num4_p2 + Num3_p2 + Num2_p2 + Node_1_p2,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2) < (N_d0 - 3))// 4dx!
			{
				inner_5node(j, i, 4, -Num9_p2 - Num8_p2 - Num7_p2 - Num6_p2, Num5_p2 + Num4_p2 + Num3_p2 + Num2_p2, Coef, Coef_location);
			}
			else if (((i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2) > (N_d0 - 3)) && ((i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2) < (N_d0 + 1)))// 2dx!
			{
				inner_5node(j, i, 2, -(Num7_p2 + Num6_p2 + N_d0 - 2), Num5_p2 - N_d0 + 2 + Num4_p2, Coef, Coef_location);
			}
			else if (((i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2) > (N_d0 + 1)) && ((i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2) < Node51_p2))// 1dx
			{
				inner_5node(j, i, 1, -(Num6_p2 + N_d0 + 2), Num5_p2 - N_d0 - 2, Coef, Coef_location);
			}
			else if (((i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2) > (Node51_p2 + 1)) && ((i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2) < (Node_5_p2 - 3)))// 1dx
			{
				inner_5node(j, i, 1, -(Num6_p2 + N_d0 + 2) - 2, Num5_p2 - N_d0 - 2 + N_c1, Coef, Coef_location);
			}
			else// 2dx!
			{
				inner_5node(j, i, 2, -(Num7_p2 - Node72_p2 + Num6_p2 + Node_5_p2), Num5_p2 - Node_5_p2 + Num4_p2 + Num3_p2 / 2 - 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post2_q3 + Num8_p2 + Num7_p2 + Num6_p2 + Node_5_p2 + 2);
		i<(M_post2_q3 + Num8_p2 + Num7_p2 + Num6_p2 + Node53_p2 + 3); i++)// circle num 33(5) :  80 node!  part2 !
	{
		switch (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2)
		{
		default:
			if ((i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2) < Node52_p2)// inner 5 node 4dx
			{
				inner_5node(j, i, 4, -(Num9_p2 - Node_9_p2 + Num8_p2 + Num7_p2 + Num6_p2 + Node_5_p2), Num5_p2 - Node_5_p2 + Num4_p2 + Num3_p2 + Num2_p2 + Node_1_p2, Coef, Coef_location);
			}
			else if ((i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2) > Node53_p2)// inner 5 node 4dx
			{
				inner_5node(j, i, 4, -(Num9_p2 - Node93_p2 + Num8_p2 + Num7_p2 + Num6_p2 + Node53_p2),
					(Num5_p2 - Node_5_p2 + Num4_p2 + Num3_p2 + Num2_p2 + Node_1_p2), Coef, Coef_location);
			}
			else// left  2lei 1j
			{
				kind2_node(j, i, -(Num7_p2 - Node72_p2 - 1 - (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Node52_p2) * 2 + Num6_p2 + (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2)),
					-1, 0, 1, (Num5_p2 - Node_5_p2 + Num4_p2 + Num3_p2 + Num2_p2 + Node_1_p2),
					(Num5_p2 - Node_5_p2 + Num4_p2 + Num3_p2 + Num2_p2 + Node_1_p2) + (Nx_2 + Nx_3),
					B2, B3, B1_1, B3, B4, B5, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post2_q3 + Num8_p2 + Num7_p2 + Num6_p2 + Node53_p2 + 3);
		i<(M_post2_q3 + Num8_p2 + Num7_p2 + Num6_p2 + Num5_p2); i++)// circle num 33(5) :  80 node!  part3 !
	{
		switch (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2)
		{
		case (Num5_p2 - 1) :
			boundary_node(j, i, 2, -Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2, Num4_p2 + Num3_p2 + Num2_p2 + Num1_p2, Coef, Coef_location);
			break;
		case (Node53_p2 + 3) ://2 lei shang 1j
			kind2_node(j, i, -(Num9_p2 - Node93_p2 + Num8_p2 + Num7_p2 + Num6_p2 + Node53_p2), -2, -1, 0, 1, (Num5_p2 - Node_5_p2 + Num4_p2 + Num3_p2 + Num2_p2 + Node_1_p2),
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node53_p2 + 7) :// 2 lei 2j shang
			kind2_node(j, i, -(Num7_p2 - Node73_p2 + Num6_p2 + Node53_p2 + 4), -2, -1, 0, 1, Num5_p2 - Node53_p2 - 4 + Num4_p2 + Num3_p2 / 2,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node54_p2:// inner 3 node
			inner_3node(j, i, -1, 0, (Num5_p2 - Node_55_p2 + 4 + Num4_p2 - 1) - N_c1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node54_p2 + 1) :// inner 3 node
			inner_3node(j, i, 0, 1, (Num5_p2 - Node_55_p2 + 4 + Num4_p2 - 1), k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Node_55_p2 - 3) ://2 lei xia 2j
			kind2_node(j, i, -(Num6_p2 + Node_55_p2 + 1), -1, 0, 1, 2, Num5_p2 - Node_55_p2 + Num4_p2 + Num3_p2 - 1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_55_p2 + 1) ://2 lei xia 1j
			kind2_node(j, i, -Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2, -1, 0, 1, 2, Num4_p2 + Num3_p2 + Num2_p2 + Num1_p2,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2) < (Node53_p2 + 7))// 2dx!
			{
				inner_5node(j, i, 2, -(Num7_p2 - Node73_p2 + Num6_p2 + Node53_p2 + 4), Num5_p2 - Node53_p2 - 4 + Num4_p2 + Num3_p2 / 2, Coef, Coef_location);
			}
			else if (((i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2) > (Node53_p2 + 7)) && ((i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2) < Node54_p2))// 1dx
			{
				inner_5node(j, i, 1, -(Node_55_p2 - 4 + 1) + 2, (Num5_p2 - Node_55_p2 + 4 + Num4_p2 - 1) - N_c1, Coef, Coef_location);
			}
			else if (((i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2) > (Node54_p2 + 1)) && ((i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2) < (Node_55_p2 - 3)))// 1dx
			{
				inner_5node(j, i, 1, -(Node_55_p2 - 4 + 1), (Num5_p2 - Node_55_p2 + 4 + Num4_p2 - 1), Coef, Coef_location);
			}
			else if (((i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2) > (Node_55_p2 - 3)) && ((i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2) < (Node_55_p2 + 1)))// 2dx
			{
				inner_5node(j, i, 2, -(Num6_p2 + Node_55_p2 + 1), Num5_p2 - Node_55_p2 + Num4_p2 + Num3_p2 - 1, Coef, Coef_location);
			}
			else// 4dx!
			{
				inner_5node(j, i, 4, -Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2, Num4_p2 + Num3_p2 + Num2_p2 + Num1_p2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!! num 34 !!!
	for (i = (M_post2_q3 + Num8_p2 + Num7_p2 + Num6_p2 + Num5_p2);
		i<(M_post2_q3 + Num8_p2 + Num7_p2 + Num6_p2 + Num5_p2 + Num4_p2 / 2); i++)// num 34(4) :  30 node! part1
	{
		switch (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2)
		{
		case 0:// 1 lei 2j shang
			kind1_node(j, i, -(Num5_p2 - N_d0 - 2) - 2, -(Num5_p2 - N_d0 - 2) - 1, -(Num5_p2 - N_d0 - 2), 0, 1, Num4_p2 + 4 - 2,
				Num4_p2 + 4 - 1, Num4_p2 + 4, A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case (Num4_p2 / 2 - 1) :// 1 lei 2j xia
			kind1_node(j, i, -(Num5_p2 - N_d0 - 2) - N_c1, -(Num5_p2 - N_d0 - 2) - N_c1 + 1, -(Num5_p2 - N_d0 - 2) - N_c1 + 2,
			-1, 0, Num4_p2 + 4, Num4_p2 + 4 + 1, Num4_p2 + 4 + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2) <= Node41_p2)// 1dx
			{
				inner_5node(j, i, 1, -(Num5_p2 - N_d0 - 2), Num4_p2 + 4, Coef, Coef_location);
			}
			else if ((i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2) >= Node42_p2) // 1dx
			{
				inner_5node(j, i, 1, -(Num5_p2 - N_d0 - 2) - N_c1, Num4_p2 + 4, Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_4node(j, i, -1, 0, 1, Num4_p2 + 4, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post2_q3 + Num8_p2 + Num7_p2 + Num6_p2 + Num5_p2 + Num4_p2 / 2);
		i<(M_post2_q3 + Num8_p2 + Num7_p2 + Num6_p2 + Num5_p2 + Num4_p2); i++)// num 34(4) :  30 node! part2
	{
		switch (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2)
		{
		case (Num4_p2 / 2) :// 1 lei 2j shang
			kind1_node(j, i, -(Num5_p2 - Node_55_p2 + 4 + Num4_p2 - 1) + N_c1 - 2, -(Num5_p2 - Node_55_p2 + 4 + Num4_p2 - 1) + N_c1 - 1,
			-(Num5_p2 - Node_55_p2 + 4 + Num4_p2 - 1) + N_c1, 0, 1, (Num3_p2 - 5 + 1) - 2, (Num3_p2 - 5 + 1) - 1, (Num3_p2 - 5 + 1),
			A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case (Num4_p2 - 1) :// 1 lei 2j xia
			kind1_node(j, i, -(Num5_p2 - Node_55_p2 + 4 + Num4_p2 - 1), -(Num5_p2 - Node_55_p2 + 4 + Num4_p2 - 1) + 1,
			-(Num5_p2 - Node_55_p2 + 4 + Num4_p2 - 1) + 2, -1, 0, (Num3_p2 - 5 + 1), (Num3_p2 - 5 + 1) + 1, (Num3_p2 - 5 + 1) + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2) <= Node43_p2)// 1dx
			{
				inner_5node(j, i, 1, -(Num5_p2 - Node_55_p2 + 4 + Num4_p2 - 1) + N_c1, (Num3_p2 - 5 + 1), Coef, Coef_location);
			}
			else if ((i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2) >= Node44_p2) // 1dx
			{
				inner_5node(j, i, 1, -(Num5_p2 - Node_55_p2 + 4 + Num4_p2 - 1), (Num3_p2 - 5 + 1), Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_4node(j, i, -1, 0, 1, (Num3_p2 - 5 + 1), 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!! num 35 !!!
	for (i = (M_post2_q3 + Num8_p2 + Num7_p2 + Num6_p2 + Num5_p2 + Num4_p2);
		i<(M_post2_q3 + Num8_p2 + Num7_p2 + Num6_p2 + Num5_p2 + Num4_p2 + Num3_p2 / 2); i++)// num 35(3) :  46 node! part1
	{
		switch (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2)
		{
		case 0:// 1 lei 1j shang
			kind1_node(j, i, -(Num5_p2 - N_d0 + 2 + Num4_p2) - 2, -(Num5_p2 - N_d0 + 2 + Num4_p2) - 1, -(Num5_p2 - N_d0 + 2 + Num4_p2),
				0, 1, Num3_p2 + Num2_p2 + N_d0 - 2 - 2, Num3_p2 + Num2_p2 + N_d0 - 2 - 1, Num3_p2 + Num2_p2 + N_d0 - 2,
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case 3:// 2 lei 2j shang
			kind2_node(j, i, -(Num5_p2 - N_d0 + 2 + Num4_p2), -2, -1, 0, 1, Num3_p2 + Num2_p2 + N_d0 - 2,
				B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case 4:// 3 lei 2j shang
			kind3_node(j, i, -(Num4_p2 + 4), -1, 0, 1, (Num3_p2 + Num2_p2 + N_d0 - 2) - 1, (Num3_p2 + Num2_p2 + N_d0 - 2),
				C2_2, C3_2, C1_2, C3_2, C4, C4, Coef, Coef_location);
			break;
		case (Num3_p2 / 2 - 5) :// 3 lei 2j xia
			kind3_node(j, i, -(Num4_p2 + 4), -1, 0, 1, Num3_p2 / 2 + 1 + Num2_p2 + Node_1_p2, Num3_p2 / 2 + 1 + Num2_p2 + Node_1_p2 + 1,
			C2_2, C3_2, C1_2, C3_2, C4, C4, Coef, Coef_location);
			break;
		case (Num3_p2 / 2 - 4) :// 2 lei 2j xia
			kind2_node(j, i, -(Num5_p2 - Node_5_p2 + Num4_p2 + Num3_p2 / 2 - 1), -1, 0, 1, 2, Num3_p2 / 2 + 1 + Num2_p2 + Node_1_p2,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num3_p2 / 2 - 1) :// 1 lei 1j xia
			kind1_node(j, i, -(Num5_p2 - Node_5_p2 + Num4_p2 + Num3_p2 / 2 - 1), -(Num5_p2 - Node_5_p2 + Num4_p2 + Num3_p2 / 2 - 1) + 1, -(Num5_p2 - Node_5_p2 + Num4_p2 + Num3_p2 / 2 - 1) + 2,
			-1, 0, Num3_p2 / 2 + 1 + Num2_p2 + Node_1_p2, Num3_p2 / 2 + 1 + Num2_p2 + Node_1_p2 + 1, Num3_p2 / 2 + 1 + Num2_p2 + Node_1_p2 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2) < 3)// 2dx
			{
				inner_5node(j, i, 2, -(Num5_p2 - N_d0 + 2 + Num4_p2), Num3_p2 + Num2_p2 + N_d0 - 2, Coef, Coef_location);
			}
			else if (((i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2) > 4) && ((i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2) < (Num3_p2 / 2 - 5)))// 1dx
			{
				inner_5node(j, i, 1, -(Num4_p2 + 4), Num3_p2 - 5, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, i, 2, -(Num5_p2 - Node_5_p2 + Num4_p2 + Num3_p2 / 2 - 1), Num3_p2 / 2 + 1 + Num2_p2 + Node_1_p2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post2_q3 + Num8_p2 + Num7_p2 + Num6_p2 + Num5_p2 + Num4_p2 + Num3_p2 / 2);
		i<(M_post2_q3 + Num8_p2 + Num7_p2 + Num6_p2 + Num5_p2 + Num4_p2 + Num3_p2); i++)// num 35(3) :  46 node! part2
	{
		switch (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2 - Num3_p2 / 2)
		{
		case 0:// 1 lei 1j shang
			kind1_node(j, i, -(Num5_p2 - Node_5_p2 - N_a + 3 + Num4_p2 + Num3_p2 / 2) - 2, -(Num5_p2 - Node_5_p2 - N_a + 3 + Num4_p2 + Num3_p2 / 2) - 1,
				-(Num5_p2 - Node_5_p2 - N_a + 3 + Num4_p2 + Num3_p2 / 2), 0, 1,
				Num3_p2 / 2 + Num2_p2 + Node_1_p2 + N_a - 3 - 2, Num3_p2 / 2 + Num2_p2 + Node_1_p2 + N_a - 3 - 1, Num3_p2 / 2 + Num2_p2 + Node_1_p2 + N_a - 3,
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case 3:// 2 lei 2j shang
			kind2_node(j, i, -(Num5_p2 - Node_5_p2 - N_a + 3 + Num4_p2 + Num3_p2 / 2), -2, -1, 0, 1, Num3_p2 / 2 + Num2_p2 + Node_1_p2 + N_a - 3,
				B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case 4:// 3 lei 2j shang
			kind3_node(j, i, -(Num3_p2 - 5 + 1), -1, 0, 1, (Num3_p2 / 2 + Num2_p2 + Node_1_p2 + N_a - 3) - 1, (Num3_p2 / 2 + Num2_p2 + Node_1_p2 + N_a - 3),
				C2_2, C3_2, C1_2, C3_2, C4, C4, Coef, Coef_location);
			break;
		case (Num3_p2 / 2 - 5) :// 3 lei 2j xia
			kind3_node(j, i, -(Num3_p2 - 5 + 1), -1, 0, 1, 1 + Num2_p2 + Node_11_p2, 1 + Num2_p2 + Node_11_p2 + 1,
			C2_2, C3_2, C1_2, C3_2, C4, C4, Coef, Coef_location);
			break;
		case (Num3_p2 / 2 - 4) :// 2 lei 2j xia
			kind2_node(j, i, -(Num5_p2 - Node_55_p2 + Num4_p2 + Num3_p2 - 1), -1, 0, 1, 2, 1 + Num2_p2 + Node_11_p2,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num3_p2 / 2 - 1) :// 1 lei 1j xia
			kind1_node(j, i, -(Num5_p2 - Node_55_p2 + Num4_p2 + Num3_p2 - 1), -(Num5_p2 - Node_55_p2 + Num4_p2 + Num3_p2 - 1) + 1,
			-(Num5_p2 - Node_55_p2 + Num4_p2 + Num3_p2 - 1) + 2, -1, 0, 1 + Num2_p2 + Node_11_p2, 1 + Num2_p2 + Node_11_p2 + 1, 1 + Num2_p2 + Node_11_p2 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2 - Num3_p2 / 2) < 3)// 2dx
			{
				inner_5node(j, i, 2, -(Num5_p2 - Node_5_p2 - N_a + 3 + Num4_p2 + Num3_p2 / 2), Num3_p2 / 2 + Num2_p2 + Node_1_p2 + N_a - 3, Coef, Coef_location);
			}
			else if (((i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2 - Num3_p2 / 2) > 4) && ((i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2 - Num3_p2 / 2) < (Num3_p2 / 2 - 5)))// 1dx
			{
				inner_5node(j, i, 1, -(Num3_p2 - 5 + 1), 6 + Num2_p2 - 1, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, i, 2, -(Num5_p2 - Node_55_p2 + Num4_p2 + Num3_p2 - 1), 1 + Num2_p2 + Node_11_p2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!! num 36 !!!
	for (i = (M_post2_q3 + Num8_p2 + Num7_p2 + Num6_p2 + Num5_p2 + Num4_p2 + Num3_p2);
		i<(M_post2_q3 + Num8_p2 + Num7_p2 + Num6_p2 + Num5_p2 + Num4_p2 + Num3_p2 + Num2_p2 / 2); i++)// num 36(2) :  26 node! part1
	{
		switch (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2 - Num3_p2)
		{
		case 0:// 3 lei 2j shang
			kind3_node(j, i, -(Num3_p2 - 3), -(Num3_p2 - 3) + 2, 0, 1, Num2_p2 + Node11_p2 - 1, Num2_p2 + Node11_p2,
				C4, C3_2, C1_2, C2_2, C4, C3_2, Coef, Coef_location);
			break;
		case (Num2_p2 / 2 - 1) ://node2: 3 lei 2j xia
			kind3_node(j, i, -(Num3_p2 - 5), -(Num3_p2 - 5) + 2, -1, 0, Num2_p2 / 2 + Node_1_p2 - 4 + 1, Num2_p2 / 2 + Node_1_p2 - 4 + 1 + 1,
			C3_2, C4, C2_2, C1_2, C3_2, C4, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2 - Num3_p2) % 2 == 0)// 1dx
			{
				inner_5node(j, i, 1, -(Num3_p2 - 5),
					(Num2_p2 - (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2 - Num3_p2) + Node11_p2 + (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2 - Num3_p2) / 2), Coef, Coef_location);
			}
			else//1 lei  2j left
			{
				kind1_node(j, i, -(Num3_p2 - 5), -1, 0, 1,
					(Num2_p2 - (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2 - Num3_p2) + Node11_p2 + (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2 - Num3_p2) / 2),
					(Num2_p2 - (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2 - Num3_p2) + Node11_p2 + (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2 - Num3_p2) / 2) + 1,
					(Num2_p2 - (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2 - Num3_p2) + Node11_p2 + (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2 - Num3_p2) / 2) + (Num1_p2 - N_d0 + 2),
					(Num2_p2 - (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2 - Num3_p2) + Node11_p2 + (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2 - Num3_p2) / 2) + 1 + (Num1_p2 - N_d0 + 2),
					A4_2, A3_2, A5_2, A3_2, A1_2, A1_2, A2_2, A2_2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post2_q3 + Num8_p2 + Num7_p2 + Num6_p2 + Num5_p2 + Num4_p2 + Num3_p2 + Num2_p2 / 2);
		i<(M_post2_q3 + Num8_p2 + Num7_p2 + Num6_p2 + Num5_p2 + Num4_p2 + Num3_p2 + Num2_p2); i++)// num 36(2) :  26 node! part2
	{
		switch (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2 - Num3_p2 - Num2_p2 / 2)
		{
		case 0:// 3 lei 2j shang
			kind3_node(j, i, -(6 + Num2_p2 - 1) - 2, -(6 + Num2_p2 - 1), 0, 1, Num2_p2 / 2 + Node13_p2 - 1, Num2_p2 / 2 + Node13_p2,
				C4, C3_2, C1_2, C2_2, C4, C3_2, Coef, Coef_location);
			break;
		case (Num2_p2 / 2 - 1) :// 3 lei 2j xia
			kind3_node(j, i, -(6 + Num2_p2 - 1), -(6 + Num2_p2 - 1) + 2, -1, 0, Node_11_p2 - 4 + 1, Node_11_p2 - 4 + 1 + 1,
			C3_2, C4, C2_2, C1_2, C3_2, C4, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2 - Num3_p2 - Num2_p2 / 2) % 2 == 0)// 1dx
			{
				inner_5node(j, i, 1, -(6 + Num2_p2 - 1), Num2_p2 / 2 - (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2 - Num3_p2 - Num2_p2 / 2) + Node13_p2 + (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2 - Num3_p2 - Num2_p2 / 2) / 2,
					Coef, Coef_location);
			}
			else//1 lei 2j left
			{
				kind1_node(j, i, -(6 + Num2_p2 - 1), -1, 0, 1,
					(Num2_p2 / 2 - (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2 - Num3_p2 - Num2_p2 / 2) + Node13_p2 + (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2 - Num3_p2 - Num2_p2 / 2) / 2),
					(Num2_p2 / 2 - (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2 - Num3_p2 - Num2_p2 / 2) + Node13_p2 + (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2 - Num3_p2 - Num2_p2 / 2) / 2) + 1,
					(Num2_p2 / 2 - (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2 - Num3_p2 - Num2_p2 / 2) + Node13_p2 + (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2 - Num3_p2 - Num2_p2 / 2) / 2) + (Num1_p2 - Node_11_p2 + Nx_3 - 1),
					(Num2_p2 / 2 - (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2 - Num3_p2 - Num2_p2 / 2) + Node13_p2 + (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2 - Num3_p2 - Num2_p2 / 2) / 2) + 1 + (Num1_p2 - Node_11_p2 + Nx_3 - 1),
					A4_2, A3_2, A5_2, A3_2, A1_2, A1_2, A2_2, A2_2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!! num 37 !!!
	for (i = (num_total - Nx_2 - Nx_3 - Num1_p2); i<(num_total - Nx_2 - Nx_3 - Num1_p2 + Node_1_p2 + 2); i++)// circle num 37(1) Nx_2:  74 node!  part1 !
	{
		switch (i - num_total + Nx_2 + Nx_3 + Num1_p2)
		{
		case 0:
			boundary_node(j, i, 1, -Num5_p2 - Num4_p2 - Num3_p2 - Num2_p2, Nx_2 + Nx_3, Coef, Coef_location);
			break;
		case (N_d0 - 3) ://2 lei shang 1j
			kind2_node(j, i, -Num5_p2 - Num4_p2 - Num3_p2 - Num2_p2, -2, -1, 0, 1, Nx_2 + Nx_3,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_1_p2 + 1) ://2 lei xia 1j
			kind2_node(j, i, -(Num5_p2 - Node_5_p2 + Num4_p2 + Num3_p2 + Num2_p2 + Node_1_p2), -1, 0, 1, 2, Nx_2 + Nx_3,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - num_total + Nx_2 + Nx_3 + Num1_p2) < (N_d0 - 3))// 4dx!
			{
				inner_5node(j, i, 4, -Num5_p2 - Num4_p2 - Num3_p2 - Num2_p2, Nx_2 + Nx_3, Coef, Coef_location);
			}
			else if (((i - num_total + Nx_2 + Nx_3 + Num1_p2) > (N_d0 - 3)) && ((i - num_total + Nx_2 + Nx_3 + Num1_p2) < Node11_p2))// 2dx!
			{
				inner_5node(j, i, 2, -(Num3_p2 + Num2_p2 + N_d0 - 2), Num1_p2 - N_d0 + 2, Coef, Coef_location);
			}
			else if (((i - num_total + Nx_2 + Nx_3 + Num1_p2) >= Node11_p2) && ((i - num_total + Nx_2 + Nx_3 + Num1_p2) <= Node12_p2))// 2 lei right 2j !
			{
				kind2_node(j, i, -((i - num_total + Nx_2 + Nx_3 + Num1_p2) + Num2_p2 - (i - num_total + Nx_2 + Nx_3 + Num1_p2 - Node11_p2) * 2),
					-1, 0, 1, Num1_p2 - N_d0 + 2, Nx_2 + Nx_3, B2, B3, B1_2, B3, B4, B5, Coef, Coef_location);
			}
			else// 2dx!
			{
				inner_5node(j, i, 2, -(Num3_p2 / 2 + 1 + Num2_p2 + Node_1_p2), Num1_p2 - N_d0 + 2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (num_total - Nx_2 - Nx_3 - Num1_p2 + Node_1_p2 + 2); i<(num_total - Nx_2 - Nx_3); i++)// circle num 37(1) Nx_2:  74 node!  part2 !
	{
		switch (i - num_total + Nx_2 + Nx_3 + Num1_p2)
		{
		case (Num1_p2 - 1) :
			boundary_node(j, i, 2, -Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2, Nx_2 + Nx_3, Coef, Coef_location);
			break;
		case (Node_1_p2 + N_a - 4) ://2 lei shang 1j
			kind2_node(j, i, -(Num5_p2 - Node_5_p2 + Num4_p2 + Num3_p2 + Num2_p2 + Node_1_p2), -2, -1, 0, 1, Nx_2 + Nx_3,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_11_p2 + 1) ://2 lei xia 1j
			kind2_node(j, i, -Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2, -1, 0, 1, 2, Nx_2 + Nx_3,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - num_total + Nx_2 + Nx_3 + Num1_p2) < (Node_1_p2 + N_a - 4))// 4dx!
			{
				inner_5node(j, i, 4, -(Num5_p2 - Node_5_p2 + Num4_p2 + Num3_p2 + Num2_p2 + Node_1_p2), Nx_2 + Nx_3, Coef, Coef_location);
			}
			else if (((i - num_total + Nx_2 + Nx_3 + Num1_p2) > (Node_1_p2 + N_a - 4)) && ((i - num_total + Nx_2 + Nx_3 + Num1_p2) < Node13_p2))// 2dx!
			{
				inner_5node(j, i, 2, -(Num3_p2 / 2 + Num2_p2 + Node_1_p2 + N_a - 3), Num1_p2 - Node_11_p2 + Nx_3 - 1, Coef, Coef_location);
			}
			else if (((i - num_total + Nx_2 + Nx_3 + Num1_p2) >= Node13_p2) && ((i - num_total + Nx_2 + Nx_3 + Num1_p2) <= Node14_p2))// 2 lei left 2j !
			{
				kind2_node(j, i, -((i - num_total + Nx_2 + Nx_3 + Num1_p2) + Num2_p2 / 2 - (i - num_total + Nx_2 + Nx_3 + Num1_p2 - Node13_p2) * 2),
					-1, 0, 1, Num1_p2 - Node_11_p2 + Nx_3 - 1, Nx_2 + Nx_3, B2, B3, B1_2, B3, B4, B5, Coef, Coef_location);
			}
			else if (((i - num_total + Nx_2 + Nx_3 + Num1_p2) >= Node14_p2) && ((i - num_total + Nx_2 + Nx_3 + Num1_p2) <= Node_11_p2))// 2dx
			{
				inner_5node(j, i, 2, -(Num2_p2 + Node_11_p2 + 1), Num1_p2 - Node_11_p2 + Nx_3 - 1, Coef, Coef_location);
			}
			else
			{
				inner_5node(j, i, 4, -Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2, Nx_2 + Nx_3, Coef, Coef_location);// 4dx!
			}
			break;
		}
		j = j + 1;
	}

	//====================================================================================================
	//调用封装函数
	//left 1 --- right 1
	scan_coef_Matrix_middle_wrapper(j, Coef, Coef_location, num_total, sizeN - 2 * Num1, offset);
	//=======================================================================================================================
	//***************************************************************************************
	//cout << "circle num1:" << j << endl;
	for (i = (Nx_2 + Nx_3); i<(Nx_2 + Nx_3 + Node_1 + 2); i++)// circle num1 Nx_2:  74 node!  part1 !
	{
		switch (i - Nx_2 - Nx_3)
		{
		case 0:
			boundary_node(j, j - offset, 1, -Nx_2 - Nx_3, Num1 + Num2 + Num3 + Num4, Coef, Coef_location);
			break;
		case (N_d0 - 3) ://2 lei shang 1j
			kind2_node(j, j - offset, -Nx_2 - Nx_3, -2, -1, 0, 1, Num1 + Num2 + Num3 + Num4,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_1 + 1) ://2 lei xia 1j
			kind2_node(j, j - offset, -(Nx_2 + Nx_3), -1, 0, 1, 2, Num1 - Node_1 + Num2 + Num3 + Num4 + Node_5,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3) < (N_d0 - 3))// 4dx!
			{
				inner_5node(j, j - offset, 4, -Nx_3 - Nx_2, Num1 + Num2 + Num3 + Num4, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3) > (N_d0 - 3)) && ((i - Nx_2 - Nx_3) < Node11))// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Nx_3 + N_d0 - 2), Num1 - N_d0 + 2 + Num2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3) >= Node11) && ((i - Nx_2 - Nx_3) <= Node12))// 2 lei left 2j !
			{
				kind2_node(j, j - offset, -Nx_3 - Nx_2, -(Nx_3 + N_d0 - 2), -1, 0, 1, Num1 - (i - Nx_2 - Nx_3) + ((i - Nx_2 - Nx_3) - Node11) * 2,
					B5, B4, B3, B1_2, B3, B2, Coef, Coef_location);
			}
			else// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Nx_3 + N_d0 - 2), Num1 - Node_1 + Num2 + Num3 / 2 - 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Node_1 + 2); i<(Nx_2 + Nx_3 + Num1); i++)// circle num1 Nx_2:  74 node!  part2 !
	{
		switch (i - Nx_2 - Nx_3)
		{
		case (Num1 - 1) :
			boundary_node(j, j - offset, 2, -Nx_2 - Nx_3, Num5 + Num2 + Num3 + Num4, Coef, Coef_location);
			break;
		case (Node_1 + N_a - 4) ://2 lei shang 1j
			kind2_node(j, j - offset, -Nx_2 - Nx_3, -2, -1, 0, 1, Num1 - Node_1 + Num2 + Num3 + Num4 + Node_5,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_11 + 1) ://2 lei xia 1j
			kind2_node(j, j - offset, -(Nx_2 + Nx_3), -1, 0, 1, 2, Num5 + Num2 + Num3 + Num4,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3) < (Node_1 + N_a - 4))// 4dx!
			{
				inner_5node(j, j - offset, 4, -Nx_3 - Nx_2, Num1 - Node_1 + Num2 + Num3 + Num4 + Node_5, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3) > (Node_1 + N_a - 4)) && ((i - Nx_2 - Nx_3) < Node13))// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Node_11 + 1), Num1 - Node_1 - N_a + 3 + Num2 + Num3 / 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3) >= Node13) && ((i - Nx_2 - Nx_3) <= Node14))// 2 lei left 2j !
			{
				kind2_node(j, j - offset, -Nx_3 - Nx_2, -(Node_11 + 1), -1, 0, 1, Num1 - (i - Nx_2 - Nx_3) + Num2 / 2 + ((i - Nx_2 - Nx_3) - Node13) * 2,
					B5, B4, B3, B1_2, B3, B2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3) >= Node14) && ((i - Nx_2 - Nx_3) <= Node_11))// 2dx
			{
				inner_5node(j, j - offset, 2, -(Node_11 + 1), Num1 - Node_11 + Num2 + Num3 - 1, Coef, Coef_location);
			}
			else
			{
				inner_5node(j, j - offset, 4, -(Nx_3 + Nx_2), Num5 + Num2 + Num3 + Num4, Coef, Coef_location);// 4dx!
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!!!!!!!!!
	//cout << "circle num 2:" << j << endl;
	for (i = (Nx_2 + Nx_3 + Num1); i<(Nx_2 + Nx_3 + Num1 + Num2 / 2); i++)// num 2 :  26 node! part1
	{
		switch (i - Nx_2 - Nx_3 - Num1)
		{
		case 0:// 3 lei 2j shang
			kind3_node(j, j - offset, -(Num1 - Node11) - 1, -(Num1 - Node11), 0, 1, Num2 + 5 - 2, Num2 + 5,
				C4, C3_2, C1_2, C2_2, C4, C3_2, Coef, Coef_location);
			break;
		case (Num2 / 2 - 1) ://node2: 3 lei 2j xia
			kind3_node(j, j - offset, -(Num1 - Node12 + Num2 / 2 - 1), -(Num1 - Node12 + Num2 / 2 - 1) + 1, -1, 0, Num2 + 5, Num2 + 5 + 2,
			C3_2, C4, C2_2, C1_2, C3_2, C4, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1) % 2 == 0)// 1dx
			{
				inner_5node(j, j - offset, 1, -((i - Nx_2 - Nx_3 - Num1) + Num1 - Node11 - (i - Nx_2 - Nx_3 - Num1) / 2), Num2 + 5, Coef, Coef_location);
			}
			else//1 lei  2j left
			{
				kind1_node(j, j - offset, -((i - Nx_2 - Nx_3 - Num1) + Num1 - Node11 - (i - Nx_2 - Nx_3 - Num1) / 2) - (Nx_3 + N_d0 - 2),
					-((i - Nx_2 - Nx_3 - Num1) + Num1 - Node11 - (i - Nx_2 - Nx_3 - Num1) / 2) - (Nx_3 + N_d0 - 2) + 1,
					-((i - Nx_2 - Nx_3 - Num1) + Num1 - Node11 - (i - Nx_2 - Nx_3 - Num1) / 2),
					-((i - Nx_2 - Nx_3 - Num1) + Num1 - Node11 - (i - Nx_2 - Nx_3 - Num1) / 2) + 1,
					-1, 0, 1, 5 + Num2, A2_2, A2_2, A1_2, A1_2, A3_2, A5_2, A3_2, A4_2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 / 2); i<(Nx_2 + Nx_3 + Num1 + Num2); i++)// num 2 :  26 node! part2
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 / 2)
		{
		case 0:// 3 lei 2j shang
			kind3_node(j, j - offset, -(Num1 - Node13 + Num2 / 2) - 1, -(Num1 - Node13 + Num2 / 2), 0, 1, Num3 - 5 - 2, Num3 - 5,
				C4, C3_2, C1_2, C2_2, C4, C3_2, Coef, Coef_location);
			break;
		case (Num2 / 2 - 1) :// 3 lei 2j xia
			kind3_node(j, j - offset, -(Num1 - Node14 + Num2 - 1), -(Num1 - Node14 + Num2 - 1) + 1, -1, 0, Num3 - 5, Num3 - 5 + 2,
			C3_2, C4, C2_2, C1_2, C3_2, C4, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 / 2) % 2 == 0)// 1dx
			{
				inner_5node(j, j - offset, 1, -((i - Nx_2 - Nx_3 - Num1 - Num2 / 2) + Num1 - Node13 - (i - Nx_2 - Nx_3 - Num1 - Num2 / 2) / 2) - Num2 / 2,
					Num3 - 5, Coef, Coef_location);
			}
			else//1 lei 2j left
			{
				kind1_node(j, j - offset, -((i - Nx_2 - Nx_3 - Num1 - Num2 / 2) + Num1 - Node13 - (i - Nx_2 - Nx_3 - Num1 - Num2 / 2) / 2) - Num2 / 2 - (Node_11 + 1),
					-((i - Nx_2 - Nx_3 - Num1 - Num2 / 2) + Num1 - Node13 - (i - Nx_2 - Nx_3 - Num1 - Num2 / 2) / 2) - Num2 / 2 - (Node_11 + 1) + 1,
					-((i - Nx_2 - Nx_3 - Num1 - Num2 / 2) + Num1 - Node13 - (i - Nx_2 - Nx_3 - Num1 - Num2 / 2) / 2) - Num2 / 2,
					-((i - Nx_2 - Nx_3 - Num1 - Num2 / 2) + Num1 - Node13 - (i - Nx_2 - Nx_3 - Num1 - Num2 / 2) / 2) - Num2 / 2 + 1,
					-1, 0, 1, Num3 - 5, A2_2, A2_2, A1_2, A1_2, A3_2, A5_2, A3_2, A4_2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!!!
	//cout << "circle num 3:" << j << endl;
	for (i = (Nx_2 + Nx_3 + Num1 + Num2); i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 / 2); i++)// num 3 :  46 node! part1
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2)
		{
		case 0:// 1 lei 1j shang
			kind1_node(j, j - offset, -(Num1 - N_d0 + 2 + Num2) - 2, -(Num1 - N_d0 + 2 + Num2) - 1, -(Num1 - N_d0 + 2 + Num2), 0, 1, Num3 + Num4 + N_d0 - 2 - 2,
				Num3 + Num4 + N_d0 - 2 - 1, Num3 + Num4 + N_d0 - 2, A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case 3:// 2 lei 2j shang
			kind2_node(j, j - offset, -(Num1 - N_d0 + 2 + Num2), -2, -1, 0, 1, Num3 + Num4 + N_d0 - 2,
				B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case 4:// 3 lei 2j shang
			kind3_node(j, j - offset, -(Num1 - N_d0 + 2 + Num2) - 1, -(Num1 - N_d0 + 2 + Num2), -1, 0, 1, Num3 - 4,
				C4, C4, C3_2, C1_2, C3_2, C2_2, Coef, Coef_location);
			break;
		case (Num3 / 2 - 5) :// 3 lei 2j xia
			kind3_node(j, j - offset, -(Num1 - Node_1 + Num2 + Num3 / 2 - 1), -(Num1 - Node_1 + Num2 + Num3 / 2 - 1) + 1, -1, 0, 1, Num3 - 4,
			C4, C4, C3_2, C1_2, C3_2, C2_2, Coef, Coef_location);
			break;
		case (Num3 / 2 - 4) :// 2 lei 2j xia
			kind2_node(j, j - offset, -(Num1 - Node_1 + Num2 + Num3 / 2 - 1), -1, 0, 1, 2, Num3 / 2 + 1 + Num4 + Node_5,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num3 / 2 - 1) :// 1 lei 1j xia
			kind1_node(j, j - offset, -(Num1 - Node_1 + Num2 + Num3 / 2 - 1), -(Num1 - Node_1 + Num2 + Num3 / 2 - 1) + 1, -(Num1 - Node_1 + Num2 + Num3 / 2 - 1) + 2, -1, 0,
			Num3 / 2 + 1 + Num4 + Node_5, Num3 / 2 + 1 + Num4 + Node_5 + 1, Num3 / 2 + 1 + Num4 + Node_5 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2) < 3)// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num1 - N_d0 + 2 + Num2), Num3 + Num4 + N_d0 - 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2) > 4) && ((i - Nx_2 - Nx_3 - Num1 - Num2) < (Num3 / 2 - 5)))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num2 + 5), Num3 - 4, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num1 - Node_1 + Num2 + Num3 / 2 - 1), Num3 / 2 + 1 + Num4 + Node_5, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 / 2); i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3); i++)// num 3 :  46 node! part2
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 / 2)
		{
		case 0:// 1 lei 1j shang
			kind1_node(j, j - offset, -(Num1 - Node_1 - N_a + 3 + Num2 + Num3 / 2) - 2, -(Num1 - Node_1 - N_a + 3 + Num2 + Num3 / 2) - 1, -(Num1 - Node_1 - N_a + 3 + Num2 + Num3 / 2), 0, 1,
				Num3 / 2 + Num4 + Node_5 + N_a - 3 - 2, Num3 / 2 + Num4 + Node_5 + N_a - 3 - 1, Num3 / 2 + Num4 + Node_5 + N_a - 3,
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case 3:// 2 lei 2j shang
			kind2_node(j, j - offset, -(Num1 - Node_1 - N_a + 3 + Num2 + Num3 / 2), -2, -1, 0, 1, Num3 / 2 + Num4 + Node_5 + N_a - 3,
				B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case 4:// 3 lei 2j shang
			kind3_node(j, j - offset, -(Num1 - Node_1 - N_a + 3 + Num2 + Num3 / 2) - 1, -(Num1 - Node_1 - N_a + 3 + Num2 + Num3 / 2), -1, 0, 1, 5 + Num4 - 1,
				C4, C4, C3_2, C1_2, C3_2, C2_2, Coef, Coef_location);
			break;
		case (Num3 / 2 - 5) :// 3 lei 2j xia
			kind3_node(j, j - offset, -(Num1 - Node_11 + Num2 + Num3 - 1), -(Num1 - Node_11 + Num2 + Num3 - 1) + 1, -1, 0, 1, 5 + Num4 - 1,
			C4, C4, C3_2, C1_2, C3_2, C2_2, Coef, Coef_location);
			break;
		case (Num3 / 2 - 4) :// 2 lei 2j xia
			kind2_node(j, j - offset, -(Num1 - Node_11 + Num2 + Num3 - 1), -1, 0, 1, 2, 1 + Num4 + Node_55,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num3 / 2 - 1) :// 1 lei 1j xia
			kind1_node(j, j - offset, -(Num1 - Node_11 + Num2 + Num3 - 1), -(Num1 - Node_11 + Num2 + Num3 - 1) + 1, -(Num1 - Node_11 + Num2 + Num3 - 1) + 2, -1, 0,
			1 + Num4 + Node_55, 1 + Num4 + Node_55 + 1, 1 + Num4 + Node_55 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 / 2) < 3)// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num1 - Node_1 - N_a + 3 + Num2 + Num3 / 2), Num3 / 2 + Num4 + Node_5 + N_a - 3, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 / 2) > 4) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 / 2) < (Num3 / 2 - 5)))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num3 - 5), 5 + Num4 - 1, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num1 - Node_11 + Num2 + Num3 - 1), 1 + Num4 + Node_55, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	//cout << "circle num 4:" << j << endl;
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3); i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 / 2); i++)// num 4 :  30 node! part1
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3)
		{
		case 0:// 1 lei 2j shang
			kind1_node(j, j - offset, -(Num3 - 4) - 2, -(Num3 - 4) - 1, -(Num3 - 4), 0, 1, Num4 + N_d0 + 2 - 2, Num4 + N_d0 + 2 - 1, Num4 + N_d0 + 2,
				A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case (Num4 / 2 - 1) :// 1 lei 2j xia
			kind1_node(j, j - offset, -(Num3 - 4), -(Num3 - 4) + 1, -(Num3 - 4) + 2, -1, 0,
			Num4 + N_d0 + 2 - N_c1, Num4 + N_d0 + 2 - N_c1 + 1, Num4 + N_d0 + 2 - N_c1 + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3) <= Node41)// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num3 - 4), Num4 + N_d0 + 2, Coef, Coef_location);
			}
			else if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3) >= Node42) // 1dx
			{
				inner_5node(j, j - offset, 1, -(Num3 - 4), Num4 + N_d0 + 2 - N_c1, Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_4node(j, j - offset, -(Num3 - 4), -1, 0, 1, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 / 2); i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4); i++)// num 4 :  30 node! part2
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3)
		{
		case (Num4 / 2) :// 1 lei 2j shang
			kind1_node(j, j - offset, -(5 + Num4 - 1) - 2, -(5 + Num4 - 1) - 1, -(5 + Num4 - 1), 0, 1,
			(Node_55 - 4 + 1 + N_c1) - 2, (Node_55 - 4 + 1 + N_c1) - 1, (Node_55 - 4 + 1 + N_c1),
			A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case (Num4 - 1) :// 1 lei 2j xia
			kind1_node(j, j - offset, -(5 + Num4 - 1), -(5 + Num4 - 1) + 1, -(5 + Num4 - 1) + 2, -1, 0,
			(Node_55 - 4 + 1), (Node_55 - 4 + 1) + 1, (Node_55 - 4 + 1) + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3) <= Node43)// 1dx
			{
				inner_5node(j, j - offset, 1, -(5 + Num4 - 1), (Node_55 - 4 + 1 + N_c1), Coef, Coef_location);
			}
			else if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3) >= Node44) // 1dx
			{
				inner_5node(j, j - offset, 1, -(5 + Num4 - 1), (Node_55 - 4 + 1), Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_4node(j, j - offset, -(5 + Num4 - 1), -1, 0, 1, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	//cout << "circle num5:" << j << endl;
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4); i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Node_5 + 2); i++)// circle num5 :  80 node!  part1 !
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4)
		{
		case 0:
			boundary_node(j, j - offset, 1, -Num1 - Num2 - Num3 - Num4, Num5 + Num6 + Num7 + Num8, Coef, Coef_location);
			break;
		case (N_d0 - 3) ://2 lei shang 1j
			kind2_node(j, j - offset, -Num1 - Num2 - Num3 - Num4, -2, -1, 0, 1, Num5 + Num6 + Num7 + Num8,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (N_d0 + 1) :// 2 lei 2j shang
			kind2_node(j, j - offset, -(Num3 + Num4 + N_d0 - 2), -2, -1, 0, 1, Num5 - N_d0 + 2 + Num6,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node51:// inner 3 node
			inner_3node(j, j - offset, -(Num4 + N_d0 + 2), -1, 0, 1, 1, k*k*dx*dx - 4, Coef, Coef_location);
			break;
		case (Node51 + 1) :// inner 3 node
			inner_3node(j, j - offset, -(Num4 + N_d0 + 2) + N_c1, 0, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node_5 - 3) ://2 lei xia 2j
			kind2_node(j, j - offset, -(Num3 / 2 + 1 + Num4 + Node_5), -1, 0, 1, 2, Num5 - Node_5 + Num6 + Num7 / 2 - 1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_5 + 1) ://2 lei xia 1j
			kind2_node(j, j - offset, -(Num1 - Node_1 + Num2 + Num3 + Num4 + Node_5), -1, 0, 1, 2, Num5 - Node_5 + Num6 + Num7 + Num8 + Node_9,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < (N_d0 - 3))// 4dx!
			{
				inner_5node(j, j - offset, 4, -Num1 - Num2 - Num3 - Num4, Num5 + Num6 + Num7 + Num8, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) > (N_d0 - 3)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < (N_d0 + 1)))// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Num3 + Num4 + N_d0 - 2), Num5 - N_d0 + 2 + Num6, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) > (N_d0 + 1)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < Node51))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num4 + N_d0 + 2), Num5 - N_d0 - 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) > (Node51 + 1)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < (Node_5 - 3)))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num4 + N_d0 + 2) + N_c1, Num5 - N_d0 - 2 - 2, Coef, Coef_location);
			}
			else// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Num3 / 2 + 1 + Num4 + Node_5), Num5 - Node_5 + Num6 + Num7 / 2 - 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Node_5 + 2); i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5); i++)// circle num5 :  80 node!  part2 !
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4)
		{
		case (Num5 - 1) :
			boundary_node(j, j - offset, 2, -Num2 - Num3 - Num4 - Num5, Num6 + Num7 + Num8 + Num9, Coef, Coef_location);
			break;
		case (Node_5 + N_a - 4) ://2 lei shang 1j
			kind2_node(j, j - offset, -(Num1 - Node_1 + Num2 + Num3 + Num4 + Node_5), -2, -1, 0, 1, Num5 - Node_5 + Num6 + Num7 + Num8 + Node_9,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_5 + N_a) :// 2 lei 2j shang
			kind2_node(j, j - offset, -(Num3 / 2 + Num4 + Node_5 + N_a - 3), -2, -1, 0, 1, Num5 - Node_5 - N_a + 3 + Num6 + Num7 / 2,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node52:// inner 3 node
			inner_3node(j, j - offset, -(Node_55 - 4 + 1) - N_c1, -1, 0, 1, 1, k*k*dx*dx - 4, Coef, Coef_location);
			break;
		case (Node52 + 1) :// inner 3 node
			inner_3node(j, j - offset, -(Node_55 - 4 + 1), 0, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node_55 - 3) ://2 lei xia 2j
			kind2_node(j, j - offset, -(Num4 + Node_55 + 1), -1, 0, 1, 2, Num5 - Node_55 + Num6 + Num7 - 1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_55 + 1) ://2 lei xia 1j
			kind2_node(j, j - offset, -Num2 - Num3 - Num4 - Num5, -1, 0, 1, 2, Num6 + Num7 + Num8 + Num9,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < (Node_5 + N_a - 4))// 4dx!
			{
				inner_5node(j, j - offset, 4, -(Num1 - Node_1 + Num2 + Num3 + Num4 + Node_5), Num5 - Node_5 + Num6 + Num7 + Num8 + Node_9, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) > (Node_5 + N_a - 4)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < (Node_5 + N_a)))// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Num3 / 2 + Num4 + Node_5 + N_a - 3), Num5 - Node_5 - N_a + 3 + Num6 + Num7 / 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) > (Node_5 + N_a)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < Node52))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Node_55 - 4 + 1) - N_c1, (Num5 - Node_55 + 4 + Num6 - 1) + 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) > (Node52 + 1)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < (Node_55 - 3)))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Node_55 - 4 + 1), (Num5 - Node_55 + 4 + Num6 - 1), Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) > (Node_55 - 3)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < (Node_55 + 1)))// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num4 + Node_55 + 1), Num5 - Node_55 + Num6 + Num7 - 1, Coef, Coef_location);
			}
			else// 4dx!
			{
				inner_5node(j, j - offset, 4, -Num2 - Num3 - Num4 - Num5, Num6 + Num7 + Num8 + Num9, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//**************************************************************
	cout << "middle_7区域节点数:" << j << endl;
	if (j == sizeM)
		cout << "middle_5 passed..." << endl;
	else
		cout << "middle_5 failed..." << endl;
	//修正矩阵，将不在本区域的节点刨除
	for (int i = 0; i != sizeM; i++)
		for (int j = 0; j != N_matrix; j++){
		if (Coef_location[i][j] < 0 || Coef_location[i][j] >= sizeN)
			Coef_location[i][j] = -1;
		}
	if (Debug)
	{
		ofstream outFile_Coef("middle_5.txt");
		if (!outFile_Coef)
		{
			cout << "Can't write the coefficient matrix into file !" << endl;
			exit(1);//Terminate with fault
		}
		else
		{
			int i, j;
			for (i = 0; i < sizeM; i++)
			{
				for (j = 0; j < N_matrix; j++)
				{
					if ((Coef[i][j].real != 0 || Coef[i][j].image != 0) && Coef_location[i][j] >= 0)
						outFile_Coef << Coef[i][j].real << "+j*" << Coef[i][j].image << " (" << i << "," << Coef_location[i][j] << ") ";
				}
				outFile_Coef << endl << endl;
			}
		}
		outFile_Coef.close();
	}

};

//*********************************************************************************//
//post1小圆孔，长为两孔一缝隙：n代表重叠区域号，Ey[n][M]代表每次迭代所求的解
//*********************************************************************************//
void scan_coef_Matrix_post1(complex** Coef, int num_total, int ** Coef_location, int sizeM, int sizeN, int offset)
{
	int i = 0, j = 0, m = 0;
	for (j = 0; j<num_total; j++)         //系数矩阵初始置零
	{
		for (i = 0; i<N_matrix; i++)
		{
			Coef[j][i].real = 0;
			Coef[j][i].image = 0;
			Coef_location[j][i] = -1;//初始化为负值
		}
	}
	j = 0;
	//系数矩阵扫描
	//******************************************************************************************
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9);
		i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Node_5 + 2); i++)// circle num13(5) :  80 node!  part1 !
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9)
		{
		case 0:
			boundary_node(j, j - offset, 1, -Num9 - Num8 - Num7 - Num6, Num5 + Num4 + Num3 + Num2, Coef, Coef_location);
			break;
		case (N_d0 - 3) ://2 lei shang 1j
			kind2_node(j, j - offset, -Num9 - Num8 - Num7 - Num6, -2, -1, 0, 1, Num5 + Num4 + Num3 + Num2,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (N_d0 + 1) :// 2 lei 2j shang
			kind2_node(j, j - offset, -(Num7 + Num6 + N_d0 - 2), -2, -1, 0, 1, Num5 - N_d0 + 2 + Num4,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node51:// inner 3 node
			inner_3node(j, j - offset, -1, 0, (Num5 - N_d0 - 2), 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node51 + 1) :// inner 3 node
			inner_3node(j, j - offset, 0, 1, (Num5 - N_d0 - 2) + N_c1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Node_5 - 3) ://2 lei xia 2j
			kind2_node(j, j - offset, -(Num7 / 2 + 1 + Num6 + Node_5), -1, 0, 1, 2, Num5 - Node_5 + Num4 + Num3 / 2 - 1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_5 + 1) ://2 lei xia 1j
			kind2_node(j, j - offset, -(Num9 - Node_9 + Num8 + Num7 + Num6 + Node_5), -1, 0, 1, 2, Num5 - Node_5 + Num4 + Num3 + Num2 + Node_1,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < (N_d0 - 3))// 4dx!
			{
				inner_5node(j, j - offset, 4, -Num9 - Num8 - Num7 - Num6, Num5 + Num4 + Num3 + Num2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) > (N_d0 - 3)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < (N_d0 + 1)))// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Num7 + Num6 + N_d0 - 2), Num5 - N_d0 + 2 + Num4, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) > (N_d0 + 1)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < Node51))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num6 + N_d0 + 2), (Num5 - N_d0 - 2), Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) > (Node51 + 1)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < (Node_5 - 3)))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num6 + N_d0 + 2) - 2, (Num5 - N_d0 - 2) + N_c1, Coef, Coef_location);
			}
			else// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Num7 / 2 + 1 + Num6 + Node_5), Num5 - Node_5 + Num4 + Num3 / 2 - 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Node_5 + 2);
		i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9); i++)// circle num13(5) :  80 node!  part2 !
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9)
		{
		case (Num5 - 1) :
			boundary_node(j, j - offset, 2, -Num8 - Num7 - Num6 - Num5, Num4 + Num3 + Num2 + Num1, Coef, Coef_location);
			break;
		case (Node_5 + N_a - 4) ://2 lei shang 1j
			kind2_node(j, j - offset, -(Num9 - Node_9 + Num8 + Num7 + Num6 + Node_5), -2, -1, 0, 1, Num5 - Node_5 + Num4 + Num3 + Num2 + Node_1,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_5 + N_a) :// 2 lei 2j shang
			kind2_node(j, j - offset, -(Num7 / 2 + Num6 + Node_5 + N_a - 3), -2, -1, 0, 1, Num5 - Node_5 - N_a + 3 + Num4 + Num3 / 2,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node52:// inner 3 node
			inner_3node(j, j - offset, -1, 0, (Num5 - Node_55 + 4 + Num4 - 1) - N_c1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node52 + 1) :// inner 3 node
			inner_3node(j, j - offset, 0, 1, (Num5 - Node_55 + 4 + Num4 - 1), k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Node_55 - 3) ://2 lei xia 2j
			kind2_node(j, j - offset, -(Num6 + Node_55 + 1), -1, 0, 1, 2, Num5 - Node_55 + Num4 + Num3 - 1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_55 + 1) ://2 lei xia 1j
			kind2_node(j, j - offset, -Num8 - Num7 - Num6 - Num5, -1, 0, 1, 2, Num4 + Num3 + Num2 + Num1,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < (Node_5 + N_a - 4))// 4dx!
			{
				inner_5node(j, j - offset, 4, -(Num9 - Node_9 + Num8 + Num7 + Num6 + Node_5), Num5 - Node_5 + Num4 + Num3 + Num2 + Node_1, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) > (Node_5 + N_a - 4)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < (Node_5 + N_a)))// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Num7 / 2 + Num6 + Node_5 + N_a - 3), Num5 - Node_5 - N_a + 3 + Num4 + Num3 / 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) > (Node_5 + N_a)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < Node52))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Node_55 - 4 + 1) + 2, (Num5 - Node_55 + 4 + Num4 - 1) - N_c1, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) > (Node52 + 1)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < (Node_55 - 3)))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Node_55 - 4 + 1), (Num5 - Node_55 + 4 + Num4 - 1), Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) > (Node_55 - 3)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < (Node_55 + 1)))// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num6 + Node_55 + 1), Num5 - Node_55 + Num4 + Num3 - 1, Coef, Coef_location);
			}
			else// 4dx!
			{
				inner_5node(j, j - offset, 4, -Num8 - Num7 - Num6 - Num5, Num4 + Num3 + Num2 + Num1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9);
		i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Num4 / 2); i++)// num 14(4) :  30 node! part1
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9)
		{
		case 0:// 1 lei 2j shang
			kind1_node(j, j - offset, -(Num5 - N_d0 - 2) - 2, -(Num5 - N_d0 - 2) - 1, -(Num5 - N_d0 - 2), 0, 1, Num4 + 4 - 2, Num4 + 4 - 1, Num4 + 4,
				A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case (Num4 / 2 - 1) :// 1 lei 2j xia
			kind1_node(j, j - offset, -(Num5 - N_d0 - 2) - N_c1, -(Num5 - N_d0 - 2) - N_c1 + 1, -(Num5 - N_d0 - 2) - N_c1 + 2, -1, 0,
			Num4 + 4, Num4 + 4 + 1, Num4 + 4 + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) <= Node41)// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num5 - N_d0 - 2), Num4 + 4, Coef, Coef_location);
			}
			else if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) >= Node42) // 1dx
			{
				inner_5node(j, j - offset, 1, -(Num5 - N_d0 - 2) - N_c1, Num4 + 4, Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_4node(j, j - offset, -1, 0, 1, Num4 + 4, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Num4 / 2);
		i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9); i++)// num14(4) :  30 node! part2
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9)
		{
		case (Num4 / 2) :// 1 lei 2j shang
			kind1_node(j, j - offset, -(Num5 - Node_55 + 4 + Num4 - 1) + N_c1 - 2, -(Num5 - Node_55 + 4 + Num4 - 1) + N_c1 - 1, -(Num5 - Node_55 + 4 + Num4 - 1) + N_c1, 0, 1,
			(Num3 - 5 + 1) - 2, (Num3 - 5 + 1) - 1, (Num3 - 5 + 1),
			A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case (Num4 - 1) :// 1 lei 2j xia
			kind1_node(j, j - offset, -(Num5 - Node_55 + 4 + Num4 - 1), -(Num5 - Node_55 + 4 + Num4 - 1) + 1, -(Num5 - Node_55 + 4 + Num4 - 1) + 2, -1, 0,
			(Num3 - 5 + 1), (Num3 - 5 + 1) + 1, (Num3 - 5 + 1) + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) <= Node43)// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num5 - Node_55 + 4 + Num4 - 1) + N_c1, (Num3 - 5 + 1), Coef, Coef_location);
			}
			else if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) >= Node44) // 1dx
			{
				inner_5node(j, j - offset, 1, -(Num5 - Node_55 + 4 + Num4 - 1), (Num3 - 5 + 1), Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_4node(j, j - offset, -1, 0, 1, (Num3 - 5 + 1), 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9);
		i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Num3 / 2); i++)// num 15(3) :  46 node! part1
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9)
		{
		case 0:// 1 lei 1j shang
			kind1_node(j, j - offset, -(Num5 - N_d0 + 2 + Num4) - 2, -(Num5 - N_d0 + 2 + Num4) - 1, -(Num5 - N_d0 + 2 + Num4), 0, 1,
				Num3 + Num2 + N_d0 - 2 - 2, Num3 + Num2 + N_d0 - 2 - 1, Num3 + Num2 + N_d0 - 2,
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case 3:// 2 lei 2j shang
			kind2_node(j, j - offset, -(Num5 - N_d0 + 2 + Num4), -2, -1, 0, 1, Num3 + Num2 + N_d0 - 2,
				B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case 4:// 3 lei 2j shang
			kind3_node(j, j - offset, -(Num4 + 4), -1, 0, 1, (Num3 + Num2 + N_d0 - 2) - 1, (Num3 + Num2 + N_d0 - 2),
				C2_2, C3_2, C1_2, C3_2, C4, C4, Coef, Coef_location);
			break;
		case (Num3 / 2 - 5) :// 3 lei 2j xia
			kind3_node(j, j - offset, -(Num4 + 4), -1, 0, 1, Num3 / 2 + 1 + Num2 + Node_1, Num3 / 2 + 1 + Num2 + Node_1 + 1,
			C2_2, C3_2, C1_2, C3_2, C4, C4, Coef, Coef_location);
			break;
		case (Num3 / 2 - 4) :// 2 lei 2j xia
			kind2_node(j, j - offset, -(Num5 - Node_5 + Num4 + Num3 / 2 - 1), -1, 0, 1, 2, Num3 / 2 + 1 + Num2 + Node_1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num3 / 2 - 1) :// 1 lei 1j xia
			kind1_node(j, j - offset, -(Num5 - Node_5 + Num4 + Num3 / 2 - 1), -(Num5 - Node_5 + Num4 + Num3 / 2 - 1) + 1, -(Num5 - Node_5 + Num4 + Num3 / 2 - 1) + 2, -1, 0,
			Num3 / 2 + 1 + Num2 + Node_1, Num3 / 2 + 1 + Num2 + Node_1 + 1, Num3 / 2 + 1 + Num2 + Node_1 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < 3)// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num5 - N_d0 + 2 + Num4), Num3 + Num2 + N_d0 - 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) > 4) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < (Num3 / 2 - 5)))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num4 + 4), Num3 - 5, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num5 - Node_5 + Num4 + Num3 / 2 - 1), Num3 / 2 + 1 + Num2 + Node_1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Num3 / 2);
		i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 * 2 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9); i++)// num 15(3) :  46 node! part2
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num3 / 2)
		{
		case 0:// 1 lei 1j shang
			kind1_node(j, j - offset, -(Num5 - Node_5 - N_a + 3 + Num4 + Num3 / 2) - 2, -(Num5 - Node_5 - N_a + 3 + Num4 + Num3 / 2) - 1, -(Num5 - Node_5 - N_a + 3 + Num4 + Num3 / 2), 0, 1,
				Num3 / 2 + Num2 + Node_1 + N_a - 3 - 2, Num3 / 2 + Num2 + Node_1 + N_a - 3 - 1, Num3 / 2 + Num2 + Node_1 + N_a - 3,
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case 3:// 2 lei 2j shang
			kind2_node(j, j - offset, -(Num5 - Node_5 - N_a + 3 + Num4 + Num3 / 2), -2, -1, 0, 1, Num3 / 2 + Num2 + Node_1 + N_a - 3,
				B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case 4:// 3 lei 2j shang
			kind3_node(j, j - offset, -(Num3 - 5 + 1), -1, 0, 1, (Num3 / 2 + Num2 + Node_1 + N_a - 3) - 1, (Num3 / 2 + Num2 + Node_1 + N_a - 3),
				C2_2, C3_2, C1_2, C3_2, C4, C4, Coef, Coef_location);
			break;
		case (Num3 / 2 - 5) :// 3 lei 2j xia
			kind3_node(j, j - offset, -(Num3 - 5 + 1), -1, 0, 1, 1 + Num2 + Node_11, 1 + Num2 + Node_11 + 1,
			C2_2, C3_2, C1_2, C3_2, C4, C4, Coef, Coef_location);
			break;
		case (Num3 / 2 - 4) :// 2 lei 2j xia
			kind2_node(j, j - offset, -(Num5 - Node_55 + Num4 + Num3 - 1), -1, 0, 1, 2, 1 + Num2 + Node_11,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num3 / 2 - 1) :// 1 lei 1j xia
			kind1_node(j, j - offset, -(Num5 - Node_55 + Num4 + Num3 - 1), -(Num5 - Node_55 + Num4 + Num3 - 1) + 1, -(Num5 - Node_55 + Num4 + Num3 - 1) + 2, -1, 0,
			1 + Num2 + Node_11, 1 + Num2 + Node_11 + 1, 1 + Num2 + Node_11 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num3 / 2) < 3)// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num5 - Node_5 - N_a + 3 + Num4 + Num3 / 2), Num3 / 2 + Num2 + Node_1 + N_a - 3, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num3 / 2) > 4) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num3 / 2) < (Num3 / 2 - 5)))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num3 - 5 + 1), 6 + Num2 - 1, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num5 - Node_55 + Num4 + Num3 - 1), 1 + Num2 + Node_11, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 * 2 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9);
		i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 * 2 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Num2 / 2); i++)// num 16(2) :  26 node! part1
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9)
		{
		case 0:// 3 lei 2j shang
			kind3_node(j, j - offset, -(Num3 - 3), -(Num3 - 3) + 2, 0, 1, Num2 + Node11 - 1, Num2 + Node11,
				C4, C3_2, C1_2, C2_2, C4, C3_2, Coef, Coef_location);
			break;
		case (Num2 / 2 - 1) ://node2: 3 lei 2j xia
			kind3_node(j, j - offset, -(Num3 - 5), -(Num3 - 5) + 2, -1, 0, Num2 / 2 + Node_1 - 4 + 1, Num2 / 2 + Node_1 - 4 + 1 + 1,
			C3_2, C4, C2_2, C1_2, C3_2, C4, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) % 2 == 0)// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num3 - 5),
					(Num2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) + Node11 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) / 2), Coef, Coef_location);
			}
			else//1 lei  2j left
			{
				kind1_node(j, j - offset, -(Num3 - 5), -1, 0, 1,
					(Num2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) + Node11 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) / 2),
					(Num2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) + Node11 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) / 2) + 1,
					(Num2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) + Node11 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) / 2) + (Num1 - N_d0 + 2),
					(Num2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) + Node11 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) / 2) + 1 + (Num1 - N_d0 + 2),
					A4_2, A3_2, A5_2, A3_2, A1_2, A1_2, A2_2, A2_2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 * 2 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Num2 / 2);
		i<(Nx_2 + Nx_3 + Num1 + Num2 * 2 + Num3 * 2 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9); i++)// num 16(2) :  26 node! part2
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2)
		{
		case 0:// 3 lei 2j shang
			kind3_node(j, j - offset, -(6 + Num2 - 1) - 2, -(6 + Num2 - 1), 0, 1, Num2 / 2 + Node13 - 1, Num2 / 2 + Node13,
				C4, C3_2, C1_2, C2_2, C4, C3_2, Coef, Coef_location);
			break;
		case (Num2 / 2 - 1) :// 3 lei 2j xia
			kind3_node(j, j - offset, -(6 + Num2 - 1), -(6 + Num2 - 1) + 2, -1, 0, Node_11 - 4 + 1, Node_11 - 4 + 1 + 1,
			C3_2, C4, C2_2, C1_2, C3_2, C4, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) % 2 == 0)// 1dx
			{
				inner_5node(j, j - offset, 1, -(6 + Num2 - 1), Num2 / 2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) + Node13 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) / 2,
					Coef, Coef_location);
			}
			else//1 lei 2j left
			{
				kind1_node(j, j - offset, -(6 + Num2 - 1), -1, 0, 1,
					(Num2 / 2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) + Node13 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) / 2),
					(Num2 / 2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) + Node13 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) / 2) + 1,
					(Num2 / 2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) + Node13 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) / 2) + (Num1 - Node_11 + Nx_3 - 1),
					(Num2 / 2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) + Node13 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) / 2) + 1 + (Num1 - Node_11 + Nx_3 - 1),
					A4_2, A3_2, A5_2, A3_2, A1_2, A1_2, A2_2, A2_2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (num_total - Nx_2 - Nx_3 - Num1); i<(num_total - Nx_2 - Nx_3 - Num1 + Node_1 + 2); i++)// circle num17(1) Nx_2:  74 node!  part1 !
	{
		switch (i - num_total + Nx_2 + Nx_3 + Num1)
		{
		case 0:
			boundary_node(j, j - offset, 1, -Num5 - Num4 - Num3 - Num2, Nx_2 + Nx_3, Coef, Coef_location);
			break;
		case (N_d0 - 3) ://2 lei shang 1j
			kind2_node(j, j - offset, -Num5 - Num4 - Num3 - Num2, -2, -1, 0, 1, Nx_2 + Nx_3,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_1 + 1) ://2 lei xia 1j
			kind2_node(j, j - offset, -(Num5 - Node_5 + Num4 + Num3 + Num2 + Node_1), -1, 0, 1, 2, Nx_2 + Nx_3,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - num_total + Nx_2 + Nx_3 + Num1) < (N_d0 - 3))// 4dx!
			{
				inner_5node(j, j - offset, 4, -Num5 - Num4 - Num3 - Num2, Nx_2 + Nx_3, Coef, Coef_location);
			}
			else if (((i - num_total + Nx_2 + Nx_3 + Num1) > (N_d0 - 3)) && ((i - num_total + Nx_2 + Nx_3 + Num1) < Node11))// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Num3 + Num2 + N_d0 - 2), Num1 - N_d0 + 2, Coef, Coef_location);
			}
			else if (((i - num_total + Nx_2 + Nx_3 + Num1) >= Node11) && ((i - num_total + Nx_2 + Nx_3 + Num1) <= Node12))// 2 lei right 2j !
			{
				kind2_node(j, j - offset, -((i - num_total + Nx_2 + Nx_3 + Num1) + Num2 - (i - num_total + Nx_2 + Nx_3 + Num1 - Node11) * 2),
					-1, 0, 1, Num1 - N_d0 + 2, Nx_2 + Nx_3, B2, B3, B1_2, B3, B4, B5, Coef, Coef_location);
			}
			else// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Num3 / 2 + 1 + Num2 + Node_1), Num1 - N_d0 + 2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (num_total - Nx_2 - Nx_3 - Num1 + Node_1 + 2); i<(num_total - Nx_2 - Nx_3); i++)// circle num17(1) Nx_2:  74 node!  part2 !
	{
		switch (i - num_total + Nx_2 + Nx_3 + Num1)
		{
		case (Num1 - 1) :
			boundary_node(j, j - offset, 2, -Num1 - Num2 - Num3 - Num4, Nx_2 + Nx_3, Coef, Coef_location);
			break;
		case (Node_1 + N_a - 4) ://2 lei shang 1j
			kind2_node(j, j - offset, -(Num5 - Node_5 + Num4 + Num3 + Num2 + Node_1), -2, -1, 0, 1, Nx_2 + Nx_3,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_11 + 1) ://2 lei xia 1j
			kind2_node(j, j - offset, -Num1 - Num2 - Num3 - Num4, -1, 0, 1, 2, Nx_2 + Nx_3,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - num_total + Nx_2 + Nx_3 + Num1) < (Node_1 + N_a - 4))// 4dx!
			{
				inner_5node(j, j - offset, 4, -(Num5 - Node_5 + Num4 + Num3 + Num2 + Node_1), Nx_2 + Nx_3, Coef, Coef_location);
			}
			else if (((i - num_total + Nx_2 + Nx_3 + Num1) > (Node_1 + N_a - 4)) && ((i - num_total + Nx_2 + Nx_3 + Num1) < Node13))// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Num3 / 2 + Num2 + Node_1 + N_a - 3), Num1 - Node_11 + Nx_3 - 1, Coef, Coef_location);
			}
			else if (((i - num_total + Nx_2 + Nx_3 + Num1) >= Node13) && ((i - num_total + Nx_2 + Nx_3 + Num1) <= Node14))// 2 lei left 2j !
			{
				kind2_node(j, j - offset, -((i - num_total + Nx_2 + Nx_3 + Num1) + Num2 / 2 - (i - num_total + Nx_2 + Nx_3 + Num1 - Node13) * 2),
					-1, 0, 1, Num1 - Node_11 + Nx_3 - 1, Nx_2 + Nx_3, B2, B3, B1_2, B3, B4, B5, Coef, Coef_location);
			}
			else if (((i - num_total + Nx_2 + Nx_3 + Num1) >= Node14) && ((i - num_total + Nx_2 + Nx_3 + Num1) <= Node_11))// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num2 + Node_11 + 1), Num1 - Node_11 + Nx_3 - 1, Coef, Coef_location);
			}
			else
			{
				inner_5node(j, j - offset, 4, -Num1 - Num2 - Num3 - Num4, Nx_2 + Nx_3, Coef, Coef_location);// 4dx!
			}
			break;
		}
		j = j + 1;
	}
	//*****************************************************************************************
	for (i = (Nx_2); i<(Nx_2 + Nx_3); i++)//left 2 Nx_3:  30 node! 
	{
		switch (i - Nx_2)
		{
		case 0:// 1 lei shang
			kind1_node(j, i, -(Nx_2 - N_d0 + 2) - 2, -(Nx_2 - N_d0 + 2) - 1, -(Nx_2 - N_d0 + 2), 0, 1, Nx_3 + N_d0 - 2 - 2, Nx_3 + N_d0 - 2 - 1, Nx_3 + N_d0 - 2,
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case (Nx_3 / 2 - 1) :// 1 lei xia
			kind1_node(j, i, -(Nx_2 - N_d0 + 2), -(Nx_2 - N_d0 + 2) + 1, -(Nx_2 - N_d0 + 2) + 2, -1, 0,
			Nx_3 + N_d0 - 2, Nx_3 + N_d0 - 2 + 1, Nx_3 + N_d0 - 2 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		case (Nx_3 / 2) :// 1 lei shang
			kind1_node(j, i, -(Nx_2 - Node_11_p1 - 1 + Nx_3) - 2, -(Nx_2 - Node_11_p1 - 1 + Nx_3) - 1, -(Nx_2 - Node_11_p1 - 1 + Nx_3),
			0, 1, 1 + Node_11_p1 - 2, 1 + Node_11_p1 - 1, 1 + Node_11_p1,
			A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case (Nx_3 - 1) :// 1 lei xia
			kind1_node(j, i, -(Nx_2 - Node_11_p1 - 1 + Nx_3), -(Nx_2 - Node_11_p1 - 1 + Nx_3) + 1, -(Nx_2 - Node_11_p1 - 1 + Nx_3) + 2,
			-1, 0, 1 + Node_11_p1, 1 + Node_11_p1 + 1, 1 + Node_11_p1 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2) < (Nx_3 / 2 - 1))
			{
				inner_5node(j, i, 2, -(Nx_2 - N_d0 + 2), Nx_3 + N_d0 - 2, Coef, Coef_location);
			}
			else
			{
				inner_5node(j, i, 2, -(Nx_2 - Node_11_p1 - 1 + Nx_3), 1 + Node_11_p1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (Nx_2 + Nx_3); i<(Nx_2 + Nx_3 + Node_1_p1 + 2); i++)// circle num1 Nx_2:  74 node!  part1 !
	{
		switch (i - Nx_2 - Nx_3)
		{
		case 0:
			boundary_node(j, i, 1, -Nx_2 - Nx_3, Num1_p1 + Num2_p1 + Num3_p1 + Num4_p1, Coef, Coef_location);
			break;
		case (N_d0 - 3) ://2 lei shang 1j
			kind2_node(j, i, -Nx_2 - Nx_3, -2, -1, 0, 1, Num1_p1 + Num2_p1 + Num3_p1 + Num4_p1,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_1_p1 + 1) ://2 lei xia 1j
			kind2_node(j, i, -(Nx_2 + Nx_3), -1, 0, 1, 2, Num1_p1 - Node_1_p1 + Num2_p1 + Num3_p1 + Num4_p1 + Node_5_p1,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3) < (N_d0 - 3))// 4dx!
			{
				inner_5node(j, i, 4, -Nx_3 - Nx_2, Num1_p1 + Num2_p1 + Num3_p1 + Num4_p1, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3) > (N_d0 - 3)) && ((i - Nx_2 - Nx_3) < Node11_p1))// 2dx!
			{
				inner_5node(j, i, 2, -(Nx_3 + N_d0 - 2), Num1_p1 - N_d0 + 2 + Num2_p1, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3) >= Node11_p1) && ((i - Nx_2 - Nx_3) <= Node12_p1))// 2 lei left 2j !
			{
				kind2_node(j, i, -Nx_3 - Nx_2, -(Nx_3 + N_d0 - 2), -1, 0, 1, Num1_p1 - (i - Nx_2 - Nx_3) + ((i - Nx_2 - Nx_3) - Node11_p1) * 2,
					B5, B4, B3, B1_2, B3, B2, Coef, Coef_location);
			}
			else// 2dx!
			{
				inner_5node(j, i, 2, -(Nx_3 + N_d0 - 2), Num1_p1 - Node_1_p1 + Num2_p1 + Num3_p1 / 2 - 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Node_1_p1 + 2); i<(Nx_2 + Nx_3 + Num1_p1); i++)// circle num1 Nx_2:  74 node!  part2 !
	{
		switch (i - Nx_2 - Nx_3)
		{
		case (Num1_p1 - 1) :
			boundary_node(j, i, 2, -Nx_2 - Nx_3, Num5_p1 + Num2_p1 + Num3_p1 + Num4_p1, Coef, Coef_location);
			break;
		case (Node_1_p1 + N_a - 4) ://2 lei shang 1j
			kind2_node(j, i, -Nx_2 - Nx_3, -2, -1, 0, 1, Num1_p1 - Node_1_p1 + Num2_p1 + Num3_p1 + Num4_p1 + Node_5_p1,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_11_p1 + 1) ://2 lei xia 1j
			kind2_node(j, i, -(Nx_2 + Nx_3), -1, 0, 1, 2, Num5_p1 + Num2_p1 + Num3_p1 + Num4_p1,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3) < (Node_1_p1 + N_a - 4))// 4dx!
			{
				inner_5node(j, i, 4, -Nx_3 - Nx_2, Num1_p1 - Node_1_p1 + Num2_p1 + Num3_p1 + Num4_p1 + Node_5_p1, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3) > (Node_1_p1 + N_a - 4)) && ((i - Nx_2 - Nx_3) < Node13_p1))// 2dx!
			{
				inner_5node(j, i, 2, -(Node_11_p1 + 1), Num1_p1 - Node_1_p1 - N_a + 3 + Num2_p1 + Num3_p1 / 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3) >= Node13_p1) && ((i - Nx_2 - Nx_3) <= Node14_p1))// 2 lei left 2j !
			{
				kind2_node(j, i, -Nx_3 - Nx_2, -(Node_11_p1 + 1), -1, 0, 1, Num1_p1 - (i - Nx_2 - Nx_3) + Num2_p1 / 2 + ((i - Nx_2 - Nx_3) - Node13_p1) * 2,
					B5, B4, B3, B1_2, B3, B2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3) >= Node14_p1) && ((i - Nx_2 - Nx_3) <= Node_11_p1))// 2dx
			{
				inner_5node(j, i, 2, -(Node_11_p1 + 1), Num1_p1 - Node_11_p1 + Num2_p1 + Num3_p1 - 1, Coef, Coef_location);
			}
			else
			{
				inner_5node(j, i, 4, -(Nx_3 + Nx_2), Num5_p1 + Num2_p1 + Num3_p1 + Num4_p1, Coef, Coef_location);// 4dx!
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!!!!!!!!!
	for (i = (Nx_2 + Nx_3 + Num1_p1); i<(Nx_2 + Nx_3 + Num1_p1 + Num2_p1 / 2); i++)// num 2 :  26 node! part1
	{
		switch (i - Nx_2 - Nx_3 - Num1_p1)
		{
		case 0:// 3 lei 2j shang
			kind3_node(j, i, -(Num1_p1 - Node11_p1) - 1, -(Num1_p1 - Node11_p1), 0, 1, Num2_p1 + 5 - 2, Num2_p1 + 5,
				C4, C3_2, C1_2, C2_2, C4, C3_2, Coef, Coef_location);
			break;
		case (Num2_p1 / 2 - 1) ://node2: 3 lei 2j xia
			kind3_node(j, i, -(Num1_p1 - Node12_p1 + Num2_p1 / 2 - 1), -(Num1_p1 - Node12_p1 + Num2_p1 / 2 - 1) + 1, -1, 0,
			Num2_p1 + 5, Num2_p1 + 5 + 2, C3_2, C4, C2_2, C1_2, C3_2, C4, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1_p1) % 2 == 0)// 1dx
			{
				inner_5node(j, i, 1, -((i - Nx_2 - Nx_3 - Num1_p1) + Num1_p1 - Node11_p1 - (i - Nx_2 - Nx_3 - Num1_p1) / 2), Num2_p1 + 5, Coef, Coef_location);
			}
			else//1 lei  2j left
			{
				kind1_node(j, i, -((i - Nx_2 - Nx_3 - Num1_p1) + Num1_p1 - Node11_p1 - (i - Nx_2 - Nx_3 - Num1_p1) / 2) - (Nx_3 + N_d0 - 2),
					-((i - Nx_2 - Nx_3 - Num1_p1) + Num1_p1 - Node11_p1 - (i - Nx_2 - Nx_3 - Num1_p1) / 2) - (Nx_3 + N_d0 - 2) + 1,
					-((i - Nx_2 - Nx_3 - Num1_p1) + Num1_p1 - Node11_p1 - (i - Nx_2 - Nx_3 - Num1_p1) / 2),
					-((i - Nx_2 - Nx_3 - Num1_p1) + Num1_p1 - Node11_p1 - (i - Nx_2 - Nx_3 - Num1_p1) / 2) + 1,
					-1, 0, 1, 5 + Num2_p1, A2_2, A2_2, A1_2, A1_2, A3_2, A5_2, A3_2, A4_2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1_p1 + Num2_p1 / 2); i<(Nx_2 + Nx_3 + Num1_p1 + Num2_p1); i++)// num 2 :  26 node! part2
	{
		switch (i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1 / 2)
		{
		case 0:// 3 lei 2j shang
			kind3_node(j, i, -(Num1_p1 - Node13_p1 + Num2_p1 / 2) - 1, -(Num1_p1 - Node13_p1 + Num2_p1 / 2), 0, 1,
				Num3_p1 - 5 - 2, Num3_p1 - 5, C4, C3_2, C1_2, C2_2, C4, C3_2, Coef, Coef_location);
			break;
		case (Num2_p1 / 2 - 1) :// 3 lei 2j xia
			kind3_node(j, i, -(Num1_p1 - Node14_p1 + Num2_p1 - 1), -(Num1_p1 - Node14_p1 + Num2_p1 - 1) + 1, -1, 0,
			Num3_p1 - 5, Num3_p1 - 5 + 2, C3_2, C4, C2_2, C1_2, C3_2, C4, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1 / 2) % 2 == 0)// 1dx
			{
				inner_5node(j, i, 1, -((i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1 / 2) + Num1_p1 - Node13_p1 - (i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1 / 2) / 2) - Num2_p1 / 2,
					Num3_p1 - 5, Coef, Coef_location);
			}
			else//1 lei 2j left
			{
				kind1_node(j, i, -((i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1 / 2) + Num1_p1 - Node13_p1 - (i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1 / 2) / 2) - Num2_p1 / 2 - (Node_11_p1 + 1),
					-((i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1 / 2) + Num1_p1 - Node13_p1 - (i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1 / 2) / 2) - Num2_p1 / 2 - (Node_11_p1 + 1) + 1,
					-((i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1 / 2) + Num1_p1 - Node13_p1 - (i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1 / 2) / 2) - Num2_p1 / 2,
					-((i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1 / 2) + Num1_p1 - Node13_p1 - (i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1 / 2) / 2) - Num2_p1 / 2 + 1,
					-1, 0, 1, Num3_p1 - 5, A2_2, A2_2, A1_2, A1_2, A3_2, A5_2, A3_2, A4_2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!!!
	for (i = (Nx_2 + Nx_3 + Num1_p1 + Num2_p1); i<(Nx_2 + Nx_3 + Num1_p1 + Num2_p1 + Num3_p1 / 2); i++)// num 3 :  46 node! part1
	{
		switch (i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1)
		{
		case 0:// 1 lei 1j shang
			kind1_node(j, i, -(Num1_p1 - N_d0 + 2 + Num2_p1) - 2, -(Num1_p1 - N_d0 + 2 + Num2_p1) - 1, -(Num1_p1 - N_d0 + 2 + Num2_p1),
				0, 1, Num3_p1 + Num4_p1 + N_d0 - 2 - 2, Num3_p1 + Num4_p1 + N_d0 - 2 - 1, Num3_p1 + Num4_p1 + N_d0 - 2,
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case 3:// 2 lei 2j shang
			kind2_node(j, i, -(Num1_p1 - N_d0 + 2 + Num2_p1), -2, -1, 0, 1, Num3_p1 + Num4_p1 + N_d0 - 2,
				B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case 4:// 3 lei 2j shang
			kind3_node(j, i, -(Num1_p1 - N_d0 + 2 + Num2_p1) - 1, -(Num1_p1 - N_d0 + 2 + Num2_p1), -1, 0, 1, Num3_p1 - 4,
				C4, C4, C3_2, C1_2, C3_2, C2_2, Coef, Coef_location);
			break;
		case (Num3_p1 / 2 - 5) :// 3 lei 2j xia
			kind3_node(j, i, -(Num1_p1 - Node_1_p1 + Num2_p1 + Num3_p1 / 2 - 1), -(Num1_p1 - Node_1_p1 + Num2_p1 + Num3_p1 / 2 - 1) + 1,
			-1, 0, 1, Num3_p1 - 4, C4, C4, C3_2, C1_2, C3_2, C2_2, Coef, Coef_location);
			break;
		case (Num3_p1 / 2 - 4) :// 2 lei 2j xia
			kind2_node(j, i, -(Num1_p1 - Node_1_p1 + Num2_p1 + Num3_p1 / 2 - 1), -1, 0, 1, 2, Num3_p1 / 2 + 1 + Num4_p1 + Node_5_p1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num3_p1 / 2 - 1) :// 1 lei 1j xia
			kind1_node(j, i, -(Num1_p1 - Node_1_p1 + Num2_p1 + Num3_p1 / 2 - 1), -(Num1_p1 - Node_1_p1 + Num2_p1 + Num3_p1 / 2 - 1) + 1,
			-(Num1_p1 - Node_1_p1 + Num2_p1 + Num3_p1 / 2 - 1) + 2, -1, 0,
			Num3_p1 / 2 + 1 + Num4_p1 + Node_5_p1, Num3_p1 / 2 + 1 + Num4_p1 + Node_5_p1 + 1, Num3_p1 / 2 + 1 + Num4_p1 + Node_5_p1 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1) < 3)// 2dx
			{
				inner_5node(j, i, 2, -(Num1_p1 - N_d0 + 2 + Num2_p1), Num3_p1 + Num4_p1 + N_d0 - 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1) > 4) && ((i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1) < (Num3_p1 / 2 - 5)))// 1dx
			{
				inner_5node(j, i, 1, -(Num2_p1 + 5), Num3_p1 - 4, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, i, 2, -(Num1_p1 - Node_1_p1 + Num2_p1 + Num3_p1 / 2 - 1), Num3_p1 / 2 + 1 + Num4_p1 + Node_5_p1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1_p1 + Num2_p1 + Num3_p1 / 2); i<(Nx_2 + Nx_3 + Num1_p1 + Num2_p1 + Num3_p1); i++)// num 3 :  46 node! part2
	{
		switch (i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1 - Num3_p1 / 2)
		{
		case 0:// 1 lei 1j shang
			kind1_node(j, i, -(Num1_p1 - Node_1_p1 - N_a + 3 + Num2_p1 + Num3_p1 / 2) - 2, -(Num1_p1 - Node_1_p1 - N_a + 3 + Num2_p1 + Num3_p1 / 2) - 1,
				-(Num1_p1 - Node_1_p1 - N_a + 3 + Num2_p1 + Num3_p1 / 2), 0, 1, Num3_p1 / 2 + Num4_p1 + Node_5_p1 + N_a - 3 - 2,
				Num3_p1 / 2 + Num4_p1 + Node_5_p1 + N_a - 3 - 1, Num3_p1 / 2 + Num4_p1 + Node_5_p1 + N_a - 3,
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case 3:// 2 lei 2j shang
			kind2_node(j, i, -(Num1_p1 - Node_1_p1 - N_a + 3 + Num2_p1 + Num3_p1 / 2), -2, -1, 0, 1,
				Num3_p1 / 2 + Num4_p1 + Node_5_p1 + N_a - 3, B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case 4:// 3 lei 2j shang
			kind3_node(j, i, -(Num1_p1 - Node_1_p1 - N_a + 3 + Num2_p1 + Num3_p1 / 2) - 1, -(Num1_p1 - Node_1_p1 - N_a + 3 + Num2_p1 + Num3_p1 / 2),
				-1, 0, 1, 5 + Num4_p1 - 1, C4, C4, C3_2, C1_2, C3_2, C2_2, Coef, Coef_location);
			break;
		case (Num3_p1 / 2 - 5) :// 3 lei 2j xia
			kind3_node(j, i, -(Num1_p1 - Node_11_p1 + Num2_p1 + Num3_p1 - 1), -(Num1_p1 - Node_11_p1 + Num2_p1 + Num3_p1 - 1) + 1,
			-1, 0, 1, 5 + Num4_p1 - 1, C4, C4, C3_2, C1_2, C3_2, C2_2, Coef, Coef_location);
			break;
		case (Num3_p1 / 2 - 4) :// 2 lei 2j xia
			kind2_node(j, i, -(Num1_p1 - Node_11_p1 + Num2_p1 + Num3_p1 - 1), -1, 0, 1, 2, 1 + Num4_p1 + Node_55_p1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num3_p1 / 2 - 1) :// 1 lei 1j xia
			kind1_node(j, i, -(Num1_p1 - Node_11_p1 + Num2_p1 + Num3_p1 - 1), -(Num1_p1 - Node_11_p1 + Num2_p1 + Num3_p1 - 1) + 1,
			-(Num1_p1 - Node_11_p1 + Num2_p1 + Num3_p1 - 1) + 2, -1, 0,
			1 + Num4_p1 + Node_55_p1, 1 + Num4_p1 + Node_55_p1 + 1, 1 + Num4_p1 + Node_55_p1 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1 - Num3_p1 / 2) < 3)// 2dx
			{
				inner_5node(j, i, 2, -(Num1_p1 - Node_1_p1 - N_a + 3 + Num2_p1 + Num3_p1 / 2), Num3_p1 / 2 + Num4_p1 + Node_5_p1 + N_a - 3, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1 - Num3_p1 / 2) > 4) && ((i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1 - Num3_p1 / 2) < (Num3_p1 / 2 - 5)))// 1dx
			{
				inner_5node(j, i, 1, -(Num3_p1 - 5), 5 + Num4_p1 - 1, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, i, 2, -(Num1_p1 - Node_11_p1 + Num2_p1 + Num3_p1 - 1), 1 + Num4_p1 + Node_55_p1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (Nx_2 + Nx_3 + Num1_p1 + Num2_p1 + Num3_p1);
		i<(Nx_2 + Nx_3 + Num1_p1 + Num2_p1 + Num3_p1 + Num4_p1 / 2); i++)// num 4 :  30 node! part1
	{
		switch (i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1 - Num3_p1)
		{
		case 0:// 1 lei 2j shang
			kind1_node(j, i, -(Num3_p1 - 4) - 2, -(Num3_p1 - 4) - 1, -(Num3_p1 - 4), 0, 1,
				Num4_p1 + N_d0 + 2 - 2, Num4_p1 + N_d0 + 2 - 1, Num4_p1 + N_d0 + 2,
				A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case (Num4_p1 / 2 - 1) :// 1 lei 2j xia
			kind1_node(j, i, -(Num3_p1 - 4), -(Num3_p1 - 4) + 1, -(Num3_p1 - 4) + 2, -1, 0,
			Num4_p1 + N_d0 + 2 - N_c1, Num4_p1 + N_d0 + 2 - N_c1 + 1, Num4_p1 + N_d0 + 2 - N_c1 + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1 - Num3_p1) <= Node41_p1)// 1dx
			{
				inner_5node(j, i, 1, -(Num3_p1 - 4), Num4_p1 + N_d0 + 2, Coef, Coef_location);
			}
			else if ((i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1 - Num3_p1) >= Node42_p1) // 1dx
			{
				inner_5node(j, i, 1, -(Num3_p1 - 4), Num4_p1 + N_d0 + 2 - N_c1, Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_4node(j, i, -(Num3_p1 - 4), -1, 0, 1, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1_p1 + Num2_p1 + Num3_p1 + Num4_p1 / 2);
		i<(Nx_2 + Nx_3 + Num1_p1 + Num2_p1 + Num3_p1 + Num4_p1); i++)// num 4 :  30 node! part2
	{
		switch (i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1 - Num3_p1)
		{
		case (Num4_p1 / 2) :// 1 lei 2j shang
			kind1_node(j, i, -(5 + Num4_p1 - 1) - 2, -(5 + Num4_p1 - 1) - 1, -(5 + Num4_p1 - 1), 0, 1,
			(Node_55_p1 - 4 + 1 + N_c1) - 2, (Node_55_p1 - 4 + 1 + N_c1) - 1, (Node_55_p1 - 4 + 1 + N_c1),
			A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case (Num4_p1 - 1) :// 1 lei 2j xia
			kind1_node(j, i, -(5 + Num4_p1 - 1), -(5 + Num4_p1 - 1) + 1, -(5 + Num4_p1 - 1) + 2, -1, 0,
			(Node_55_p1 - 4 + 1), (Node_55_p1 - 4 + 1) + 1, (Node_55_p1 - 4 + 1) + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1 - Num3_p1) <= Node43_p1)// 1dx
			{
				inner_5node(j, i, 1, -(5 + Num4_p1 - 1), (Node_55_p1 - 4 + 1 + N_c1), Coef, Coef_location);
			}
			else if ((i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1 - Num3_p1) >= Node44_p1) // 1dx
			{
				inner_5node(j, i, 1, -(5 + Num4_p1 - 1), (Node_55_p1 - 4 + 1), Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_4node(j, i, -(5 + Num4_p1 - 1), -1, 0, 1, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (Nx_2 + Nx_3 + Num1_p1 + Num2_p1 + Num3_p1 + Num4_p1);
		i<(Nx_2 + Nx_3 + Num1_p1 + Num2_p1 + Num3_p1 + Num4_p1 + Node_5_p1 + 2); i++)// circle num5 :  80 node!  part1 !
	{
		switch (i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1 - Num3_p1 - Num4_p1)
		{
		case 0:
			boundary_node(j, i, 1, -Num1_p1 - Num2_p1 - Num3_p1 - Num4_p1, Num5_p1 + Num6_p1 + Num7_p1 + Num8_p1, Coef, Coef_location);
			break;
		case (N_d0 - 3) ://2 lei shang 1j
			kind2_node(j, i, -Num1_p1 - Num2_p1 - Num3_p1 - Num4_p1, -2, -1, 0, 1, Num5_p1 + Num6_p1 + Num7_p1 + Num8_p1,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (N_d0 + 1) :// 2 lei 2j shang
			kind2_node(j, i, -(Num3_p1 + Num4_p1 + N_d0 - 2), -2, -1, 0, 1, Num5_p1 - N_d0 + 2 + Num6_p1,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node51_p1:// inner 3 node
			inner_3node(j, i, -(Num4_p1 + N_d0 + 2), -1, 0, 1, 1, k*k*dx*dx - 4, Coef, Coef_location);
			break;
		case (Node51_p1 + 1) :// inner 3 node
			inner_3node(j, i, -(Num4_p1 + N_d0 + 2) + N_c1, 0, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node_5_p1 - 3) ://2 lei xia 2j
			kind2_node(j, i, -(Num3_p1 / 2 + 1 + Num4_p1 + Node_5_p1), -1, 0, 1, 2, Num5_p1 - Node_5_p1 + Num6_p1 + Num7_p1 / 2 - 1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_5_p1 + 1) ://2 lei xia 1j
			kind2_node(j, i, -(Num1_p1 - Node_1_p1 + Num2_p1 + Num3_p1 + Num4_p1 + Node_5_p1), -1, 0, 1, 2,
			Num5_p1 - Node_5_p1 + Num6_p1 + Num7_p1 + Num8_p1 + Node_9_p1, B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1 - Num3_p1 - Num4_p1) < (N_d0 - 3))// 4dx!
			{
				inner_5node(j, i, 4, -Num1_p1 - Num2_p1 - Num3_p1 - Num4_p1, Num5_p1 + Num6_p1 + Num7_p1 + Num8_p1, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1 - Num3_p1 - Num4_p1) > (N_d0 - 3)) && ((i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1 - Num3_p1 - Num4_p1) < (N_d0 + 1)))// 2dx!
			{
				inner_5node(j, i, 2, -(Num3_p1 + Num4_p1 + N_d0 - 2), Num5_p1 - N_d0 + 2 + Num6_p1, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1 - Num3_p1 - Num4_p1) > (N_d0 + 1)) && ((i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1 - Num3_p1 - Num4_p1) < Node51_p1))// 1dx
			{
				inner_5node(j, i, 1, -(Num4_p1 + N_d0 + 2), Num5_p1 - N_d0 - 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1 - Num3_p1 - Num4_p1) > (Node51_p1 + 1)) && ((i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1 - Num3_p1 - Num4_p1) < (Node_5_p1 - 3)))// 1dx
			{
				inner_5node(j, i, 1, -(Num4_p1 + N_d0 + 2) + N_c1, Num5_p1 - N_d0 - 2 - 2, Coef, Coef_location);
			}
			else// 2dx!
			{
				inner_5node(j, i, 2, -(Num3_p1 / 2 + 1 + Num4_p1 + Node_5_p1), Num5_p1 - Node_5_p1 + Num6_p1 + Num7_p1 / 2 - 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1_p1 + Num2_p1 + Num3_p1 + Num4_p1 + Node_5_p1 + 2);
		i<(Nx_2 + Nx_3 + Num1_p1 + Num2_p1 + Num3_p1 + Num4_p1 + Num5_p1); i++)// circle num5 :  83 node!  part2 !
	{
		switch (i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1 - Num3_p1 - Num4_p1)
		{
		case (Num5_p1 - 1) :
			boundary_node(j, i, 2, -Num2_p1 - Num3_p1 - Num4_p1 - Num5_p1, Num6_p1 + Num7_p1 + Num8_p1 + Num9_p1, Coef, Coef_location);
			break;
		case (Node_5_p1 + N_a - 4) ://2 lei shang 1j
			kind2_node(j, i, -(Num1_p1 - Node_1_p1 + Num2_p1 + Num3_p1 + Num4_p1 + Node_5_p1), -2, -1, 0, 1,
			Num5_p1 - Node_5_p1 + Num6_p1 + Num7_p1 + Num8_p1 + Node_9_p1, B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_5_p1 + N_a) :// 2 lei 2j shang
			kind2_node(j, i, -(Num3_p1 / 2 + Num4_p1 + Node_5_p1 + N_a - 3), -2, -1, 0, 1, Num5_p1 - Node_5_p1 - N_a + 3 + Num6_p1 + Num7_p1 / 2,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node52_p1:// inner 3 node
			inner_3node(j, i, -(Node_55_p1 - 4 + 1) - N_c1, -1, 0, 1, 1, k*k*dx*dx - 4, Coef, Coef_location);
			break;
		case (Node52_p1 + 1) :// inner 3 node
			inner_3node(j, i, -(Node_55_p1 - 4 + 1), 0, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node_55_p1 - 3) ://2 lei xia 2j
			kind2_node(j, i, -(Num4_p1 + Node_55_p1 + 1), -1, 0, 1, 2, Num5_p1 - Node_55_p1 + Num6_p1 + Num7_p1 - 1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_55_p1 + 1) ://2 lei xia 1j
			kind2_node(j, i, -Num2_p1 - Num3_p1 - Num4_p1 - Num5_p1, -1, 0, 1, 2, Num6_p1 + Num7_p1 + Num8_p1 + Num9_p1,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1 - Num3_p1 - Num4_p1) < (Node_5_p1 + N_a - 4))// 4dx!
			{
				inner_5node(j, i, 4, -(Num1_p1 - Node_1_p1 + Num2_p1 + Num3_p1 + Num4_p1 + Node_5_p1), Num5_p1 - Node_5_p1 + Num6_p1 + Num7_p1 + Num8_p1 + Node_9_p1, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1 - Num3_p1 - Num4_p1) > (Node_5_p1 + N_a - 4)) && ((i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1 - Num3_p1 - Num4_p1) < (Node_5_p1 + N_a)))// 2dx!
			{
				inner_5node(j, i, 2, -(Num3_p1 / 2 + Num4_p1 + Node_5_p1 + N_a - 3), Num5_p1 - Node_5_p1 - N_a + 3 + Num6_p1 + Num7_p1 / 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1 - Num3_p1 - Num4_p1) > (Node_5_p1 + N_a)) && ((i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1 - Num3_p1 - Num4_p1) < Node52_p1))// 1dx
			{
				inner_5node(j, i, 1, -(Node_55_p1 - 4 + 1) - N_c1, (Num5_p1 - Node_55_p1 + 4 + Num6_p1 - 1) + 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1 - Num3_p1 - Num4_p1) > (Node52_p1 + 1)) && ((i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1 - Num3_p1 - Num4_p1) < (Node_55_p1 - 3)))// 1dx
			{
				inner_5node(j, i, 1, -(Node_55_p1 - 4 + 1), (Num5_p1 - Node_55_p1 + 4 + Num6_p1 - 1), Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1 - Num3_p1 - Num4_p1) > (Node_55_p1 - 3)) && ((i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1 - Num3_p1 - Num4_p1) < (Node_55_p1 + 1)))// 2dx
			{
				inner_5node(j, i, 2, -(Num4_p1 + Node_55_p1 + 1), Num5_p1 - Node_55_p1 + Num6_p1 + Num7_p1 - 1, Coef, Coef_location);
			}
			else// 4dx!
			{
				inner_5node(j, i, 4, -Num2_p1 - Num3_p1 - Num4_p1 - Num5_p1, Num6_p1 + Num7_p1 + Num8_p1 + Num9_p1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (Nx_2 + Nx_3 + Num1_p1 + Num2_p1 + Num3_p1 + Num4_p1 + Num5_p1);
		i<(Nx_2 + Nx_3 + Num1_p1 + Num2_p1 + Num3_p1 + Num4_p1 + Num5_p1 + Num6_p1); i++)// num 6 :  16 node!
	{
		switch (i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1 - Num3_p1 - Num4_p1 - Num5_p1)
		{
		case 0:// 1 lei 2j shang
			kind1_node(j, i, -(Num5_p1 - N_d0 - 2) - 2, -(Num5_p1 - N_d0 - 2) - 1, -(Num5_p1 - N_d0 - 2), 0, 1,
				Num6_p1 + 4 - 2, Num6_p1 + 4 - 1, Num6_p1 + 4, A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case Node61_p1:// inner 3 node
			inner_3node(j, i, -(Num5_p1 - N_d0 - 2), -1, 0, 1, 1, k*k*dx*dx - 4, Coef, Coef_location);
			break;
		case (Node61_p1 + 1) :// inner 3 node
			inner_3node(j, i, -(Num5_p1 - N_d0 - 2 - 2), 0, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Num6_p1 / 2 - 1) :// 1 lei 2j xia
			kind1_node(j, i, -(Num5_p1 - N_d0 - 2 - 2), -(Num5_p1 - N_d0 - 2 - 2) + 1, -(Num5_p1 - N_d0 - 2 - 2) + 2, -1, 0,
			Num6_p1 + 4 - 2, Num6_p1 + 4 - 2 + 1, Num6_p1 + 4 - 2 + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		case (Num6_p1 / 2) :// 1 lei 2j shang
			kind1_node(j, i, -(Num5_p1 - Node_55_p1 + 4 + Num6_p1 - 1) - 2 - 2, -(Num5_p1 - Node_55_p1 + 4 + Num6_p1 - 1) - 2 - 1,
			-(Num5_p1 - Node_55_p1 + 4 + Num6_p1 - 1) - 2, 0, 1, Num7_p1 - 4 + 2 - 2, Num7_p1 - 4 + 2 - 1, Num7_p1 - 4 + 2,
			A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case Node62_p1:// inner 3 node
			inner_3node(j, i, -(Num5_p1 - Node_55_p1 + 4 + Num6_p1 - 1) - 2, -1, 0, 1, 1, k*k*dx*dx - 4, Coef, Coef_location);
			break;
		case (Node62_p1 + 1) :// inner 3 node
			inner_3node(j, i, -(Num5_p1 - Node_55_p1 + 4 + Num6_p1 - 1), 0, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Num6_p1 - 1) :// 1 lei 2j xia
			kind1_node(j, i, -(Num5_p1 - Node_55_p1 + 4 + Num6_p1 - 1), -(Num5_p1 - Node_55_p1 + 4 + Num6_p1 - 1) + 1,
			-(Num5_p1 - Node_55_p1 + 4 + Num6_p1 - 1) + 2, -1, 0, Num7_p1 - 4, Num7_p1 - 4 + 1, Num7_p1 - 4 + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1 - Num3_p1 - Num4_p1 - Num5_p1) < Node61_p1)// 1dx
			{
				inner_5node(j, i, 1, -(Num5_p1 - N_d0 - 2), Num6_p1 + 4, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1 - Num3_p1 - Num4_p1 - Num5_p1) > (Node61_p1 + 1)) && ((i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1 - Num3_p1 - Num4_p1 - Num5_p1) < (Num6_p1 / 2 - 1)))// 2dx
			{
				inner_5node(j, i, 1, -(Num5_p1 - N_d0 - 2 - 2), Num6_p1 + 4 - 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1 - Num3_p1 - Num4_p1 - Num5_p1) > (Num6_p1 / 2)) && ((i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1 - Num3_p1 - Num4_p1 - Num5_p1) < Node62_p1))// 2dx
			{
				inner_5node(j, i, 1, -(Num5_p1 - Node_55_p1 + 4 + Num6_p1 - 1) - 2, Num7_p1 - 4 + 2, Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_5node(j, i, 1, -(Num5_p1 - Node_55_p1 + 4 + Num6_p1 - 1), Num7_p1 - 4, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!!
	for (i = (Nx_2 + Nx_3 + Num1_p1 + Num2_p1 + Num3_p1 + Num4_p1 + Num5_p1 + Num6_p1);
		i<(Nx_2 + Nx_3 + Num1_p1 + Num2_p1 + Num3_p1 + Num4_p1 + Num5_p1 + Num6_p1 + Num7_p1 / 2); i++)// num 7 :  28 node! part1
	{
		switch (i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1 - Num3_p1 - Num4_p1 - Num5_p1 - Num6_p1)
		{
		case 0:// 1 lei 1j shang
			kind1_node(j, i, -(Num5_p1 - N_d0 + 2 + Num6_p1) - 2, -(Num5_p1 - N_d0 + 2 + Num6_p1) - 1, -(Num5_p1 - N_d0 + 2 + Num6_p1), 0, 1,
				Num7_p1 + Num8_p1 + N_d0 - 2 - 2, Num7_p1 + Num8_p1 + N_d0 - 2 - 1, Num7_p1 + Num8_p1 + N_d0 - 2,
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case 3:// 2 lei 2j shang
			kind2_node(j, i, -(Num5_p1 - N_d0 + 2 + Num6_p1), -2, -1, 0, 1, Num7_p1 + Num8_p1 + N_d0 - 2,
				B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node71_p1://inner 4 node
			inner_4node(j, i, -(Num6_p1 + 4), -1, 0, Num7_p1 - 4, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node71_p1 + 1) ://inner 4 node
			inner_4node(j, i, -(Num6_p1 + 4) + 2, 0, 1, Num7_p1 - 4, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Num7_p1 / 2 - 4) :// 2 lei 2j xia
			kind2_node(j, i, -(Num5_p1 - Node_5_p1 + Num6_p1 + Num7_p1 / 2 - 1), -1, 0, 1, 2, Num7_p1 / 2 + 1 + Num8_p1 + Node_9_p1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num7_p1 / 2 - 1) :// 1 lei 1j xia
			kind1_node(j, i, -(Num5_p1 - Node_5_p1 + Num6_p1 + Num7_p1 / 2 - 1), -(Num5_p1 - Node_5_p1 + Num6_p1 + Num7_p1 / 2 - 1) + 1,
			-(Num5_p1 - Node_5_p1 + Num6_p1 + Num7_p1 / 2 - 1) + 2, -1, 0, Num7_p1 / 2 + 1 + Num8_p1 + Node_9_p1,
			Num7_p1 / 2 + 1 + Num8_p1 + Node_9_p1 + 1, Num7_p1 / 2 + 1 + Num8_p1 + Node_9_p1 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1 - Num3_p1 - Num4_p1 - Num5_p1 - Num6_p1) < 3)// 2dx
			{
				inner_5node(j, i, 2, -(Num5_p1 - N_d0 + 2 + Num6_p1), Num7_p1 + Num8_p1 + N_d0 - 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1 - Num3_p1 - Num4_p1 - Num5_p1 - Num6_p1) > 3) && ((i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1 - Num3_p1 - Num4_p1 - Num5_p1 - Num6_p1) < (Node71_p1)))// 1dx
			{
				inner_5node(j, i, 1, -(Num6_p1 + 4), Num7_p1 - 4, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1 - Num3_p1 - Num4_p1 - Num5_p1 - Num6_p1) > (Node71_p1 + 1)) && ((i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1 - Num3_p1 - Num4_p1 - Num5_p1 - Num6_p1) < (Num7_p1 / 2 - 4)))// 1dx
			{
				inner_5node(j, i, 1, -(Num6_p1 + 4) + 2, Num7_p1 - 4, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, i, 2, -(Num5_p1 - Node_5_p1 + Num6_p1 + Num7_p1 / 2 - 1), Num7_p1 / 2 + 1 + Num8_p1 + Node_9_p1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1_p1 + Num2_p1 + Num3_p1 + Num4_p1 + Num5_p1 + Num6_p1 + Num7_p1 / 2);
		i<(Nx_2 + Nx_3 + Num1_p1 + Num2_p1 + Num3_p1 + Num4_p1 + Num5_p1 + Num6_p1 + Num7_p1); i++)// num 7 :  28 node! part2
	{
		switch (i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1 - Num3_p1 - Num4_p1 - Num5_p1 - Num6_p1)
		{
		case (Num7_p1 / 2) :// 1 lei 1j shang
			kind1_node(j, i, -(Num5_p1 - Node_5_p1 - N_a + 3 + Num6_p1 + Num7_p1 / 2) - 2, -(Num5_p1 - Node_5_p1 - N_a + 3 + Num6_p1 + Num7_p1 / 2) - 1,
			-(Num5_p1 - Node_5_p1 - N_a + 3 + Num6_p1 + Num7_p1 / 2), 0, 1, Num7_p1 / 2 + Num8_p1 + Node_9_p1 + N_a - 3 - 2,
			Num7_p1 / 2 + Num8_p1 + Node_9_p1 + N_a - 3 - 1, Num7_p1 / 2 + Num8_p1 + Node_9_p1 + N_a - 3,
			A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case (Num7_p1 / 2 + 3) :// 2 lei 2j shang
			kind2_node(j, i, -(Num5_p1 - Node_5_p1 - N_a + 3 + Num6_p1 + Num7_p1 / 2), -2, -1, 0, 1, Num7_p1 / 2 + Num8_p1 + Node_9_p1 + N_a - 3,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node72_p1://inner 4 node
			inner_4node(j, i, -(Num7_p1 - 5 + 1) - 2, -1, 0, 5 + Num8_p1 - 1, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node72_p1 + 1) ://inner 4 node
			inner_4node(j, i, -(Num7_p1 - 5 + 1), 0, 1, 5 + Num8_p1 - 1, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Num7_p1 - 4) :// 2 lei 2j xia
			kind2_node(j, i, -(Num5_p1 - Node_55_p1 + Num6_p1 + Num7_p1 - 1), -1, 0, 1, 2, 1 + Num8_p1 + Node_99_p1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num7_p1 - 1) :// 1 lei 1j xia
			kind1_node(j, i, -(Num5_p1 - Node_55_p1 + Num6_p1 + Num7_p1 - 1), -(Num5_p1 - Node_55_p1 + Num6_p1 + Num7_p1 - 1) + 1,
			-(Num5_p1 - Node_55_p1 + Num6_p1 + Num7_p1 - 1) + 2, -1, 0, 1 + Num8_p1 + Node_99_p1,
			1 + Num8_p1 + Node_99_p1 + 1, 1 + Num8_p1 + Node_99_p1 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1 - Num3_p1 - Num4_p1 - Num5_p1 - Num6_p1) < (Num7_p1 / 2 + 3))// 2dx
			{
				inner_5node(j, i, 2, -(Num5_p1 - Node_5_p1 - N_a + 3 + Num6_p1 + Num7_p1 / 2), Num7_p1 / 2 + Num8_p1 + Node_9_p1 + N_a - 3, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1 - Num3_p1 - Num4_p1 - Num5_p1 - Num6_p1) > (Num7_p1 / 2 + 3)) && ((i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1 - Num3_p1 - Num4_p1 - Num5_p1 - Num6_p1) < (Node72_p1)))// 1dx
			{
				inner_5node(j, i, 1, -(Num7_p1 - 5 + 1) - 2, 5 + Num8_p1 - 1, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1 - Num3_p1 - Num4_p1 - Num5_p1 - Num6_p1) > (Node72_p1 + 1)) && ((i - Nx_2 - Nx_3 - Num1_p1 - Num2_p1 - Num3_p1 - Num4_p1 - Num5_p1 - Num6_p1) < (Num7_p1 - 4)))// 1dx
			{
				inner_5node(j, i, 1, -(Num7_p1 - 5 + 1), 5 + Num8_p1 - 1, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, i, 2, -(Num5_p1 - Node_55_p1 + Num6_p1 + Num7_p1 - 1), 1 + Num8_p1 + Node_99_p1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (M_post1_q1 - Num9_p1 - Num8_p1);
		i<(M_post1_q1 - Num9_p1); i++)// num 8 :  12 node! 
	{
		switch (i - M_post1_q1 + Num9_p1 + Num8_p1)
		{
		case 0:// 1 lei shang 2j
			kind1_node(j, i, -(Num7_p1 - 4) - 2, -(Num7_p1 - 4) - 1, -(Num7_p1 - 4), 0, 1,
				Num8_p1 + N_d0 + 2 - 2, Num8_p1 + N_d0 + 2 - 1, Num8_p1 + N_d0 + 2,
				A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case 1:
			inner_5node(j, i, 1, -(Num7_p1 - 4), Num8_p1 + N_d0 + 2, Coef, Coef_location);
			break;
		case 2://inner 4 node
			inner_4node(j, i, -(Num7_p1 - 4), -1, 0, Num8_p1 + N_d0 + 2, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case 3://inner 4 node
			inner_4node(j, i, -(Num7_p1 - 4), 0, 1, Num8_p1 + N_d0 + 2, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case 4:
			inner_5node(j, i, 1, -(Num7_p1 - 4), Num8_p1 + N_d0 + 2, Coef, Coef_location);
			break;
		case (Num8_p1 / 2 - 1) :// 1 lei xia 2j
			kind1_node(j, i, -(Num7_p1 - 4), -(Num7_p1 - 4) + 1, -(Num7_p1 - 4) + 2, -1, 0, Num8_p1 + N_d0 + 2, Num8_p1 + N_d0 + 2 + 1, Num8_p1 + N_d0 + 2 + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		case (Num8_p1 / 2) :// 1 lei shang 2j
			kind1_node(j, i, -(5 + Num8_p1 - 1) - 2, -(5 + Num8_p1 - 1) - 1, -(5 + Num8_p1 - 1), 0, 1,
			Node_99_p1 - 4 + 1 - 2, Node_99_p1 - 4 + 1 - 1, Node_99_p1 - 4 + 1,
			A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case 7:
			inner_5node(j, i, 1, -(5 + Num8_p1 - 1), Node_99_p1 - 4 + 1, Coef, Coef_location);
			break;
		case 8://inner 4 node
			inner_4node(j, i, -(5 + Num8_p1 - 1), -1, 0, Node_99_p1 - 4 + 1, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case 9://inner 4 node
			inner_4node(j, i, -(5 + Num8_p1 - 1), 0, 1, Node_99_p1 - 4 + 1, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case 10:
			inner_5node(j, i, 1, -(5 + Num8_p1 - 1), Node_99_p1 - 4 + 1, Coef, Coef_location);
			break;
		case (Num8_p1 - 1) :// 1 lei xia 2j
			kind1_node(j, i, -(5 + Num8_p1 - 1), -(5 + Num8_p1 - 1) + 1, -(5 + Num8_p1 - 1) + 2, -1, 0,
			Node_99_p1 - 4 + 1, Node_99_p1 - 4 + 1 + 1, Node_99_p1 - 4 + 1 + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (M_post1_q1 - Num9_p1);
		i<(M_post1_q1 - Num9_p1 + Node_9_p1 + 2); i++)// circle num9 :  72 node!  part1 !
	{
		switch (i - M_post1_q1 + Num9_p1)
		{
		case 0:
			boundary_node(j, i, 1, -Num5_p1 - Num6_p1 - Num7_p1 - Num8_p1, Num10_p1 + Num11_p1 + Num12_p1 + Num9_p1, Coef, Coef_location);
			break;
		case (N_d0 - 3) ://2 lei shang 1j
			kind2_node(j, i, -Num5_p1 - Num6_p1 - Num7_p1 - Num8_p1, -2, -1, 0, 1, Num10_p1 + Num11_p1 + Num12_p1 + Num9_p1,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (N_d0 + 1) :// 2 lei 2j shang
			kind2_node(j, i, -(Num7_p1 + Num8_p1 + N_d0 - 2), -2, -1, 0, 1, Num9_p1 - N_d0 + 2 + Num10_p1,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node91_p1:// inner 4 node
			inner_4node(j, i, -(Num8_p1 + N_d0 + 2), -1, 0, Num9_p1 - N_d0 - 2, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node91_p1 + 1) :// inner 4 node
			inner_4node(j, i, -(Num8_p1 + N_d0 + 2), 0, 1, Num9_p1 - N_d0 - 2, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Node_9_p1 - 3) ://2 lei xia 2j
			kind2_node(j, i, -(Num7_p1 / 2 + 1 + Num8_p1 + Node_9_p1), -1, 0, 1, 2, Num9_p1 - Node_9_p1 + Num10_p1 + Node112_p1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_9_p1 + 1) ://2 lei xia 1j
			kind2_node(j, i, -(Num5_p1 - Node_5_p1 + Num6_p1 + Num7_p1 + Num8_p1 + Node_9_p1), -1, 0, 1, 2,
			Num9_p1 - Node_9_p1 + Num10_p1 + Num11_p1 + Num12_p1 + Node_13_p1, B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - M_post1_q1 + Num9_p1) < (N_d0 - 3))// 4dx!
			{
				inner_5node(j, i, 4, -Num5_p1 - Num6_p1 - Num7_p1 - Num8_p1, Num10_p1 + Num11_p1 + Num12_p1 + Num9_p1, Coef, Coef_location);
			}
			else if (((i - M_post1_q1 + Num9_p1) > (N_d0 - 3)) && ((i - M_post1_q1 + Num9_p1) < (N_d0 + 1)))// 2dx!
			{
				inner_5node(j, i, 2, -(Num7_p1 + Num8_p1 + N_d0 - 2), Num9_p1 - N_d0 + 2 + Num10_p1, Coef, Coef_location);
			}
			else if (((i - M_post1_q1 + Num9_p1) > (N_d0 + 1)) && ((i - M_post1_q1 + Num9_p1) < Node91_p1))// 1dx
			{
				inner_5node(j, i, 1, -(Num8_p1 + N_d0 + 2), Num9_p1 - N_d0 - 2, Coef, Coef_location);
			}
			else if (((i - M_post1_q1 + Num9_p1) > (Node91_p1 + 1)) && ((i - M_post1_q1 + Num9_p1) < (Node_9_p1 - 3)))// 1dx
			{
				inner_5node(j, i, 1, -(Num8_p1 + N_d0 + 2), Num9_p1 - N_d0 - 2, Coef, Coef_location);
			}
			else// 2dx!
			{
				inner_5node(j, i, 2, -(Num7_p1 / 2 + 1 + Num8_p1 + Node_9_p1), Num9_p1 - Node_9_p1 + Num10_p1 + Node112_p1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post1_q1 - Num9_p1 + Node_9_p1 + 2); i<(M_post1_q1 - Num9_p1 + Node93_p1 + 4); i++)// circle num9 :  72 node!  part2 !
	{
		switch (i - M_post1_q1 + Num9_p1 - Node_9_p1 - 2)
		{
		default:
			if ((i - M_post1_q1 + Num9_p1 - Node_9_p1 - 2) < 3)
			{
				inner_5node(j, i, 4, -(Num5_p1 - Node_5_p1 + Num6_p1 + Num7_p1 + Num8_p1 + Node_9_p1),
					Num9_p1 - Node_9_p1 + Num10_p1 + Num11_p1 + Num12_p1 + Node_13_p1, Coef, Coef_location);
			}
			else if (((i - M_post1_q1 + Num9_p1 - Node_9_p1 - 2) >= 3) && ((i - M_post1_q1 + Num9_p1 - Node_9_p1 - 2) <= 6))// 2lei zuo 1j
			{
				kind2_node(j, i, -(Num1_p1 - Node_1_p1 + Num2_p1 + Num3_p1 + Num4_p1 + Num5_p1 + Num6_p1 + Num7_p1 + Num8_p1 + Node_9_p1),
					-(Num5_p1 - Node_5_p1 + Num6_p1 + Num7_p1 + Num8_p1 + Node_9_p1), -1, 0, 1,
					Num9_p1 - Node92_p1 - (i - M_post1_q1 + Num9_p1 - Node_9_p1 - 2 - 3) + Num10_p1 + Node112_p1 + 1 + (i - M_post1_q1 + Num9_p1 - Node_9_p1 - 2 - 3) * 2,
					B5, B4, B3, B1_1, B3, B2, Coef, Coef_location);
			}
			else
			{
				inner_5node(j, i, 4, -(Num5_p1 - Node_5_p1 + Num6_p1 + Num7_p1 + Num8_p1 + Node_9_p1),
					Num9_p1 - Node93_p1 - 1 + Num10_p1 + Num11_p1 + Num12_p1 + Node133_p1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post1_q1 - Num9_p1 + Node93_p1 + 4); i<(M_post1_q1); i++)// circle num9 :  79 node!  part3 !
	{
		switch (i - M_post1_q1 + Num9_p1)
		{
		case (Num9_p1 - 1) :
			boundary_node(j, i, 2, -Num6_p1 - Num7_p1 - Num8_p1 - Num9_p1, Num10_p1 + Num11_p1 + Num12_p1 + Num13_p1, Coef, Coef_location);
			break;
		case (Node93_p1 + 4) ://2 lei shang 1j
			kind2_node(j, i, -(Num5_p1 - Node_5_p1 + Num6_p1 + Num7_p1 + Num8_p1 + Node_9_p1), -2, -1, 0, 1,
			Num9_p1 - Node93_p1 - 1 + Num10_p1 + Num11_p1 + Num12_p1 + Node133_p1, B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node94_p1 - 3) :// 2 lei 2j shang
			kind2_node(j, i, -(1 + Num8_p1 + Node_99_p1), -2, -1, 0, 1,
			Num9_p1 - Node_99_p1 + Num10_p1 + Num11_p1 - 1, B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node94_p1:// inner 4 node
			inner_4node(j, i, -(Node_99_p1 - 4 + 1), -1, 0, Num9_p1 - Node_99_p1 + 4 + Num10_p1 - 1, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node94_p1 + 1) :// inner 4 node
			inner_4node(j, i, -(Node_99_p1 - 4 + 1), 0, 1, Num9_p1 - Node_99_p1 + 4 + Num10_p1 - 1, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Node_99_p1 - 3) ://2 lei xia 2j
			kind2_node(j, i, -(Num8_p1 + Node_99_p1 + 1), -1, 0, 1, 2, Num9_p1 - Node_99_p1 + Num10_p1 + Num11_p1 - 1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_99_p1 + 1) ://2 lei xia 1j
			kind2_node(j, i, -Num6_p1 - Num7_p1 - Num8_p1 - Num9_p1, -1, 0, 1, 2, Num10_p1 + Num11_p1 + Num12_p1 + Num13_p1,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - M_post1_q1 + Num9_p1) < (Node94_p1 - 3))// 2dx!
			{
				inner_5node(j, i, 2, -(1 + Num8_p1 + Node_99_p1), Num9_p1 - Node_99_p1 + Num10_p1 + Num11_p1 - 1, Coef, Coef_location);
			}
			else if (((i - M_post1_q1 + Num9_p1) > (Node94_p1 - 3)) && ((i - M_post1_q1 + Num9_p1) < Node94_p1))// 1dx
			{
				inner_5node(j, i, 1, -(Node_99_p1 - 4 + 1), Num9_p1 - Node_99_p1 + 4 + Num10_p1 - 1, Coef, Coef_location);
			}
			else if (((i - M_post1_q1 + Num9_p1) > (Node94_p1 + 1)) && ((i - M_post1_q1 + Num9_p1) < (Node_99_p1 - 3)))// 1dx
			{
				inner_5node(j, i, 1, -(Node_99_p1 - 4 + 1), Num9_p1 - Node_99_p1 + 4 + Num10_p1 - 1, Coef, Coef_location);
			}
			else if (((i - M_post1_q1 + Num9_p1) > (Node_99_p1 - 3)) && ((i - M_post1_q1 + Num9_p1) < (Node_99_p1 + 1)))// 2dx
			{
				inner_5node(j, i, 2, -(1 + Num8_p1 + Node_99_p1), Num9_p1 - Node_99_p1 + Num10_p1 + Num11_p1 - 1, Coef, Coef_location);
			}
			else// 4dx!
			{
				inner_5node(j, i, 4, -Num6_p1 - Num7_p1 - Num8_p1 - Num9_p1, Num10_p1 + Num11_p1 + Num12_p1 + Num13_p1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (M_post1_q1); i<(M_post1_q1 + Num10_p1); i++)// num 10 :  12 node! 
	{
		switch (i - M_post1_q1)
		{
		case 0:// 1 lei shang 2j
			kind1_node(j, i, -(Num9_p1 - N_d0 - 2) - 2, -(Num9_p1 - N_d0 - 2) - 1, -(Num9_p1 - N_d0 - 2), 0, 1, Num10_p1 + 4 - 2,
				Num10_p1 + 4 - 1, Num10_p1 + 4, A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case 1:
			inner_5node(j, i, 1, -(Num9_p1 - N_d0 - 2), Num10_p1 + 4, Coef, Coef_location);
			break;
		case 2://inner 4 node
			inner_4node(j, i, -(Num9_p1 - N_d0 - 2), -1, 0, Num10_p1 + 4, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case 3://inner 4 node
			inner_4node(j, i, -(Num9_p1 - N_d0 - 2), 0, 1, Num10_p1 + 4, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case 4:
			inner_5node(j, i, 1, -(Num9_p1 - N_d0 - 2), Num10_p1 + 4, Coef, Coef_location);
			break;
		case (Num10_p1 / 2 - 1) :// 1 lei xia 2j
			kind1_node(j, i, -(Num9_p1 - N_d0 - 2), -(Num9_p1 - N_d0 - 2) + 1, -(Num9_p1 - N_d0 - 2) + 2, -1, 0, Num10_p1 + 4,
			Num10_p1 + 4 + 1, Num10_p1 + 4 + 2, A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		case (Num10_p1 / 2) :// 1 lei shang 2j
			kind1_node(j, i, -(Num9_p1 - Node_99_p1 + 4 + Num10_p1 - 1) - 2, -(Num9_p1 - Node_99_p1 + 4 + Num10_p1 - 1) - 1,
			-(Num9_p1 - Node_99_p1 + 4 + Num10_p1 - 1), 0, 1, Num11_p1 - 4 - 2, Num11_p1 - 4 - 1, Num11_p1 - 4,
			A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case 7:
			inner_5node(j, i, 1, -(Num9_p1 - Node_99_p1 + 4 + Num10_p1 - 1), Num11_p1 - 4, Coef, Coef_location);
			break;
		case 8://inner 4 node
			inner_4node(j, i, -(Num9_p1 - Node_99_p1 + 4 + Num10_p1 - 1), -1, 0, Num11_p1 - 4, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case 9://inner 4 node
			inner_4node(j, i, -(Num9_p1 - Node_99_p1 + 4 + Num10_p1 - 1), 0, 1, Num11_p1 - 4, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case 10:
			inner_5node(j, i, 1, -(Num9_p1 - Node_99_p1 + 4 + Num10_p1 - 1), Num11_p1 - 4, Coef, Coef_location);
			break;
		case (Num10_p1 - 1) :// 1 lei xia 2j
			kind1_node(j, i, -(Num9_p1 - Node_99_p1 + 4 + Num10_p1 - 1), -(Num9_p1 - Node_99_p1 + 4 + Num10_p1 - 1) + 1,
			-(Num9_p1 - Node_99_p1 + 4 + Num10_p1 - 1) + 2, -1, 0, Num11_p1 - 4, Num11_p1 - 4 + 1, Num11_p1 - 4 + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (M_post1_q1 + Num10_p1); i<(M_post1_q1 + Num10_p1 + Node112_p1 + 1); i++)// num 11 :  35 node! part1
	{
		switch (i - M_post1_q1 - Num10_p1)
		{
		case 0:// 1 lei 1j shang
			kind1_node(j, i, -(Num9_p1 - N_d0 + 2 + Num10_p1) - 2, -(Num9_p1 - N_d0 + 2 + Num10_p1) - 1, -(Num9_p1 - N_d0 + 2 + Num10_p1), 0, 1,
				Num11_p1 + Num12_p1 + N_d0 - 2 - 2, Num11_p1 + Num12_p1 + N_d0 - 2 - 1, Num11_p1 + Num12_p1 + N_d0 - 2,
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case 3:// 2 lei 2j shang
			kind2_node(j, i, -(Num9_p1 - N_d0 + 2 + Num10_p1), -2, -1, 0, 1, Num11_p1 + Num12_p1 + N_d0 - 2,
				B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node111_p1://inner 4 node
			inner_4node(j, i, -(Num10_p1 + 4), -1, 0, Num11_p1 - 4, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node111_p1 + 1) ://inner 4 node
			inner_4node(j, i, -(Num10_p1 + 4), 0, 1, Num11_p1 - 4 + 2, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Node112_p1 - 3) :// 2 lei 2j xia
			kind2_node(j, i, -(Num9_p1 - Node_9_p1 + Num10_p1 + Node112_p1), -1, 0, 1, 2, Num11_p1 - Node112_p1 + Num12_p1 + Node_13_p1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case Node112_p1:// 1 lei 1j xia
			kind1_node(j, i, -(Num9_p1 - Node_9_p1 + Num10_p1 + Node112_p1), -(Num9_p1 - Node_9_p1 + Num10_p1 + Node112_p1) + 1,
				-(Num9_p1 - Node_9_p1 + Num10_p1 + Node112_p1) + 2, -1, 0, Num11_p1 - Node112_p1 + Num12_p1 + Node_13_p1,
				Num11_p1 - Node112_p1 + Num12_p1 + Node_13_p1 + 1, Num11_p1 - Node112_p1 + Num12_p1 + Node_13_p1 + 2,
				A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - M_post1_q1 - Num10_p1) < 3)// 2dx
			{
				inner_5node(j, i, 2, -(Num9_p1 - N_d0 + 2 + Num10_p1), Num11_p1 + Num12_p1 + N_d0 - 2, Coef, Coef_location);
			}
			else if (((i - M_post1_q1 - Num10_p1) > 3) && ((i - M_post1_q1 - Num10_p1) < (Node111_p1)))// 1dx
			{
				inner_5node(j, i, 1, -(Num10_p1 + 4), Num11_p1 - 4, Coef, Coef_location);
			}
			else if (((i - M_post1_q1 - Num10_p1) > (Node111_p1 + 1)) && ((i - M_post1_q1 - Num10_p1) < (Node112_p1 - 3)))// 1dx
			{
				inner_5node(j, i, 1, -(Num10_p1 + 4), Num11_p1 - 4 + 2, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, i, 2, -(Num9_p1 - Node_9_p1 + Num10_p1 + Node112_p1), Num11_p1 - Node112_p1 + Num12_p1 + Node_13_p1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post1_q1 + Num10_p1 + Node112_p1 + 1); i<(M_post1_q1 + Num10_p1 + Node113_p1); i++)// num 11 :  35 node! part2
	{
		switch (i - M_post1_q1 - Num10_p1 - Node112_p1 - 1)
		{
		case 0:// 3 lei 1j shang
			kind3_node(j, i, -(Num9_p1 - Node92_p1 + Num10_p1 + Node112_p1 + 1) - 1, -(Num9_p1 - Node92_p1 + Num10_p1 + Node112_p1 + 1),
				0, 1, (Num11_p1 - Node112_p1 - 1 + Num12_p1 + Node132_p1), (Num11_p1 - Node112_p1 - 1 + Num12_p1 + Node132_p1) + 2,
				C4, C3_1, C1_1, C2_1, C4, C3_1, Coef, Coef_location);
			break;
		case 1:// 1 lei 1j zuo
			kind1_node(j, i, -(Num9_p1 - Node92_p1 + Num10_p1 + Node112_p1 + 1) - 1 - (Num5_p1 - Node_5_p1 + Num6_p1 + Num7_p1 + Num8_p1 + Node_9_p1),
				-(Num9_p1 - Node92_p1 + Num10_p1 + Node112_p1 + 1) - (Num5_p1 - Node_5_p1 + Num6_p1 + Num7_p1 + Num8_p1 + Node_9_p1),
				-(Num9_p1 - Node92_p1 + Num10_p1 + Node112_p1 + 1) - 1, -(Num9_p1 - Node92_p1 + Num10_p1 + Node112_p1 + 1),
				-1, 0, 1, (Num11_p1 - Node112_p1 - 1 + Num12_p1 + Node132_p1) + 2,
				A2_1, A2_1, A1_1, A1_1, A3_1, A5_1, A3_1, A4_1, Coef, Coef_location);
			break;
		case 2:// inner 5 node
			inner_5node(j, i, 2, -(Num9_p1 - Node92_p1 + Num10_p1 + Node112_p1 + 1) - 1,
				(Num11_p1 - Node112_p1 - 1 + Num12_p1 + Node132_p1) + 2, Coef, Coef_location);
			break;
		case 3:// 1 lei 1j zuo
			kind1_node(j, i, -(Num9_p1 - Node92_p1 + Num10_p1 + Node112_p1 + 1) - 2 - (Num5_p1 - Node_5_p1 + Num6_p1 + Num7_p1 + Num8_p1 + Node_9_p1),
				-(Num9_p1 - Node92_p1 + Num10_p1 + Node112_p1 + 1) - 1 - (Num5_p1 - Node_5_p1 + Num6_p1 + Num7_p1 + Num8_p1 + Node_9_p1),
				-(Num9_p1 - Node92_p1 + Num10_p1 + Node112_p1 + 1) - 2, -(Num9_p1 - Node92_p1 + Num10_p1 + Node112_p1 + 1) - 1,
				-1, 0, 1, (Num11_p1 - Node112_p1 - 1 + Num12_p1 + Node132_p1) + 2,
				A2_1, A2_1, A1_1, A1_1, A3_1, A5_1, A3_1, A4_1, Coef, Coef_location);
			break;
		case 4:// inner 5 node
			inner_5node(j, i, 2, -(Num9_p1 - Node92_p1 + Num10_p1 + Node112_p1 + 1) - 2,
				(Num11_p1 - Node112_p1 - 1 + Num12_p1 + Node132_p1) + 2, Coef, Coef_location);
			break;
		case 5:// 1 lei 1j zuo
			kind1_node(j, i, -(Num9_p1 - Node92_p1 + Num10_p1 + Node112_p1 + 1) - 3 - (Num5_p1 - Node_5_p1 + Num6_p1 + Num7_p1 + Num8_p1 + Node_9_p1),
				-(Num9_p1 - Node92_p1 + Num10_p1 + Node112_p1 + 1) - 2 - (Num5_p1 - Node_5_p1 + Num6_p1 + Num7_p1 + Num8_p1 + Node_9_p1),
				-(Num9_p1 - Node92_p1 + Num10_p1 + Node112_p1 + 1) - 3, -(Num9_p1 - Node92_p1 + Num10_p1 + Node112_p1 + 1) - 2,
				-1, 0, 1, (Num11_p1 - Node112_p1 - 1 + Num12_p1 + Node132_p1) + 2,
				A2_1, A2_1, A1_1, A1_1, A3_1, A5_1, A3_1, A4_1, Coef, Coef_location);
			break;
		case 6:// 3 lei 1j xia
			kind3_node(j, i, -(Num9_p1 - Node92_p1 + Num10_p1 + Node112_p1 + 1) - 3, -(Num9_p1 - Node92_p1 + Num10_p1 + Node112_p1 + 1) - 3 + 1,
				-1, 0, (Num11_p1 - Node112_p1 - 1 + Num12_p1 + Node132_p1) + 2, (Num11_p1 - Node112_p1 - 1 + Num12_p1 + Node132_p1) + 2 + 2,
				C3_1, C4, C2_1, C1_1, C3_1, C4, Coef, Coef_location);
			break;
		}
		j = j + 1;

	}
	for (i = (M_post1_q1 + Num10_p1 + Node113_p1); i<(M_post1_q1 + Num10_p1 + Num11_p1); i++)// num 11 :  35 node! part3
	{
		switch (i - M_post1_q1 - Num10_p1)
		{
		case Node113_p1:// 1 lei 1j shang
			kind1_node(j, i, -(Num9_p1 - Node93_p1 - 5 + Num10_p1 + Node113_p1) - 2, -(Num9_p1 - Node93_p1 - 5 + Num10_p1 + Node113_p1) - 1,
				-(Num9_p1 - Node93_p1 - 5 + Num10_p1 + Node113_p1), 0, 1, Num11_p1 - Node113_p1 + Num12_p1 + Node133_p1 + 4 - 2,
				Num11_p1 - Node113_p1 + Num12_p1 + Node133_p1 + 4 - 1, Num11_p1 - Node113_p1 + Num12_p1 + Node133_p1 + 4,
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case (Node113_p1 + 3) :// 2 lei 2j shang
			kind2_node(j, i, -(Num9_p1 - Node93_p1 - 5 + Num10_p1 + Node113_p1), -2, -1, 0, 1, Num11_p1 - Node113_p1 + Num12_p1 + Node133_p1 + 4,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node114_p1://inner 4 node
			inner_4node(j, i, -(Num11_p1 - 5 + 1), -1, 0, 5 + Num12_p1 - 1 - 2, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node114_p1 + 1) ://inner 4 node
			inner_4node(j, i, -(Num11_p1 - 5 + 1), 0, 1, 5 + Num12_p1 - 1, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Num11_p1 - 4) :// 2 lei 2j xia
			kind2_node(j, i, -(Num9_p1 - Node_99_p1 + Num10_p1 + Num11_p1 - 1), -1, 0, 1, 2, 1 + Num12_p1 + Node_1313_p1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num11_p1 - 1) :// 1 lei 1j xia
			kind1_node(j, i, -(Num9_p1 - Node_99_p1 + Num10_p1 + Num11_p1 - 1), -(Num9_p1 - Node_99_p1 + Num10_p1 + Num11_p1 - 1) + 1,
			-(Num9_p1 - Node_99_p1 + Num10_p1 + Num11_p1 - 1) + 2, -1, 0,
			1 + Num12_p1 + Node_1313_p1, 1 + Num12_p1 + Node_1313_p1 + 1, 1 + Num12_p1 + Node_1313_p1 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - M_post1_q1 - Num10_p1) < (Node113_p1 + 3))// 2dx
			{
				inner_5node(j, i, 2, -(Num9_p1 - Node93_p1 - 5 + Num10_p1 + Node113_p1), Num11_p1 - Node113_p1 + Num12_p1 + Node133_p1 + 4, Coef, Coef_location);
			}
			else if (((i - M_post1_q1 - Num10_p1) > (Node113_p1 + 3)) && ((i - M_post1_q1 - Num10_p1) < (Node114_p1)))// 1dx
			{
				inner_5node(j, i, 1, -(Num11_p1 - 5 + 1), 5 + Num12_p1 - 1 - 2, Coef, Coef_location);
			}
			else if (((i - M_post1_q1 - Num10_p1) > (Node114_p1 + 1)) && ((i - M_post1_q1 - Num10_p1) < (Num11_p1 - 4)))// 1dx
			{
				inner_5node(j, i, 1, -(Num11_p1 - 5 + 1), 5 + Num12_p1 - 1, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, i, 2, -(Num9_p1 - Node_99_p1 + Num10_p1 + Num11_p1 - 1), 1 + Num12_p1 + Node_1313_p1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (M_post1_q1 + Num10_p1 + Num11_p1);
		i<(M_post1_q1 + Num10_p1 + Num11_p1 + Num12_p1); i++)// num 12 :  16 node!  
	{
		switch (i - M_post1_q1 - Num10_p1 - Num11_p1)
		{
		case 0:// 1 lei 2j shang
			kind1_node(j, i, -(Num11_p1 - 4) - 2, -(Num11_p1 - 4) - 1, -(Num11_p1 - 4), 0, 1, Num12_p1 + N_d0 + 2 - 2,
				Num12_p1 + N_d0 + 2 - 1, Num12_p1 + N_d0 + 2, A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case Node121_p1:// inner 3 node
			inner_3node(j, i, -1, 0, (Num12_p1 + N_d0 + 2), 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node121_p1 + 1) :// inner 3 node
			inner_3node(j, i, 0, 1, (Num12_p1 + N_d0 + 2) + 2, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Num12_p1 / 2 - 1) :// 1 lei 2j xia
			kind1_node(j, i, (-(Num11_p1 - 4) - 2), (-(Num11_p1 - 4) - 2) + 1, (-(Num11_p1 - 4) - 2) + 2, -1, 0,
			(Num12_p1 + N_d0 + 2) + 2, (Num12_p1 + N_d0 + 2) + 2 + 1, (Num12_p1 + N_d0 + 2) + 2 + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		case (Num12_p1 / 2) :// 1 lei 2j shang
			kind1_node(j, i, -(5 + Num12_p1 - 1) + 2 - 2, -(5 + Num12_p1 - 1) + 2 - 1, -(5 + Num12_p1 - 1) + 2, 0, 1,
			(Node_1313_p1 - 4 + 1) - 2 - 2, (Node_1313_p1 - 4 + 1) - 2 - 1, (Node_1313_p1 - 4 + 1) - 2,
			A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case Node122_p1:// inner 3 node
			inner_3node(j, i, -1, 0, (Node_1313_p1 - 4 + 1) - 2, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node122_p1 + 1) :// inner 3 node
			inner_3node(j, i, 0, 1, (Node_1313_p1 - 4 + 1), k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Num12_p1 - 1) :// 1 lei 2j xia
			kind1_node(j, i, -(5 + Num12_p1 - 1), -(5 + Num12_p1 - 1) + 1, -(5 + Num12_p1 - 1) + 2, -1, 0,
			(Node_1313_p1 - 4 + 1), (Node_1313_p1 - 4 + 1) + 1, (Node_1313_p1 - 4 + 1) + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - M_post1_q1 - Num10_p1 - Num11_p1) < Node121_p1)// 1dx
			{
				inner_5node(j, i, 1, -(Num11_p1 - 4), Num12_p1 + N_d0 + 2, Coef, Coef_location);
			}
			else if (((i - M_post1_q1 - Num10_p1 - Num11_p1) > (Node121_p1 + 1)) && ((i - M_post1_q1 - Num10_p1 - Num11_p1) < (Num12_p1 / 2 - 1)))// 2dx
			{
				inner_5node(j, i, 1, -(Num11_p1 - 4) - 2, (Num12_p1 + N_d0 + 2) + 2, Coef, Coef_location);
			}
			else if (((i - M_post1_q1 - Num10_p1 - Num11_p1) > (Num12_p1 / 2)) && ((i - M_post1_q1 - Num10_p1 - Num11_p1) < Node122_p1))// 2dx
			{
				inner_5node(j, i, 1, -(5 + Num12_p1 - 1) + 2, (Node_1313_p1 - 4 + 1) - 2, Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_5node(j, i, 1, -(5 + Num12_p1 - 1), (Node_1313_p1 - 4 + 1), Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (M_post1_q1 + Num10_p1 + Num11_p1 + Num12_p1);
		i<(M_post1_q1 + Num10_p1 + Num11_p1 + Num12_p1 + Node_13_p1 + 2); i++)// circle num13 :  95 node!  part1 !
	{
		switch (i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1)
		{
		case 0:
			boundary_node(j, i, 1, -Num9_p1 - Num10_p1 - Num11_p1 - Num12_p1, Num13_p1 + Num14_p1 + Num15_p1 + Num16_p1, Coef, Coef_location);
			break;
		case (N_d0 - 3) ://2 lei shang 1j
			kind2_node(j, i, -Num9_p1 - Num10_p1 - Num11_p1 - Num12_p1, -2, -1, 0, 1, Num13_p1 + Num14_p1 + Num15_p1 + Num16_p1,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (N_d0 + 1) :// 2 lei 2j shang
			kind2_node(j, i, -(Num11_p1 + Num12_p1 + N_d0 - 2), -2, -1, 0, 1, Num13_p1 - N_d0 + 2 + Num14_p1,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node131_p1:// inner 3 node
			inner_3node(j, i, -1, 0, (Num13_p1 - N_d0 - 2), 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node131_p1 + 1) :// inner 3 node
			inner_3node(j, i, 0, 1, (Num13_p1 - N_d0 - 2) + N_c1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Node_13_p1 - 3) ://2 lei xia 2j
			kind2_node(j, i, -(Num11_p1 - Node112_p1 + Num12_p1 + Node_13_p1), -1, 0, 1, 2, Num13_p1 - Node_13_p1 + Num14_p1 + Node151_p1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_13_p1 + 1) ://2 lei xia 1j
			kind2_node(j, i, -(Num9_p1 - Node_9_p1 + Num10_p1 + Num11_p1 + Num12_p1 + Node_13_p1), -1, 0, 1, 2,
			Num13_p1 - Node_13_p1 + Num14_p1 + Num15_p1 + Num16_p1 + Node_17_p1, B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1) < (N_d0 - 3))// 4dx!
			{
				inner_5node(j, i, 4, -Num9_p1 - Num10_p1 - Num11_p1 - Num12_p1, Num13_p1 + Num14_p1 + Num15_p1 + Num16_p1, Coef, Coef_location);
			}
			else if (((i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1) > (N_d0 - 3)) && ((i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1) < (N_d0 + 1)))// 2dx!
			{
				inner_5node(j, i, 2, -(Num11_p1 + Num12_p1 + N_d0 - 2), Num13_p1 - N_d0 + 2 + Num14_p1, Coef, Coef_location);
			}
			else if (((i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1) > (N_d0 + 1)) && ((i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1) < Node131_p1))// 1dx
			{
				inner_5node(j, i, 1, -(Num12_p1 + N_d0 + 2), (Num13_p1 - N_d0 - 2), Coef, Coef_location);
			}
			else if (((i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1) > (Node131_p1 + 1)) && ((i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1) < (Node_13_p1 - 3)))// 1dx
			{
				inner_5node(j, i, 1, -(Num12_p1 + N_d0 + 2) - 2, (Num13_p1 - N_d0 - 2) + N_c1, Coef, Coef_location);
			}
			else// 2dx!
			{
				inner_5node(j, i, 2, -(Num11_p1 - Node112_p1 + Num12_p1 + Node_13_p1), Num13_p1 - Node_13_p1 + Num14_p1 + Node151_p1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post1_q1 + Num10_p1 + Num11_p1 + Num12_p1 + Node_13_p1 + 2);
		i<(M_post1_q1 + Num10_p1 + Num11_p1 + Num12_p1 + Node133_p1 + 3); i++)// circle num13 :  95 node!  part2 !
	{
		switch (i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1)
		{
		case (Node132_p1) ://2 lei shang
			kind2_node(j, i, -(Num9_p1 - Node_9_p1 + Num10_p1 + Num11_p1 + Num12_p1 + Node_13_p1), -2, -1, 0, 1,
			Num13_p1 - Node_13_p1 + Num14_p1 + Num15_p1 + Num16_p1 + Node_17_p1, B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node132_p1 + 1) ://3 lei shang
			kind3_node(j, i, -(Num9_p1 - Node_9_p1 + Num10_p1 + Num11_p1 + Num12_p1 + Node_13_p1) - 1, -(Num9_p1 - Node_9_p1 + Num10_p1 + Num11_p1 + Num12_p1 + Node_13_p1),
			-1, 0, 1, Num13_p1 - Node132_p1 + Num14_p1 + Node151_p1, C4, C4, C3_1, C1_1, C3_1, C2_1, Coef, Coef_location);
			break;
		case (Node133_p1 - 1) ://3 lei xia
			kind3_node(j, i, -(Num9_p1 - Node93_p1 - 1 + Num10_p1 + Num11_p1 + Num12_p1 + Node133_p1),
			-(Num9_p1 - Node93_p1 - 1 + Num10_p1 + Num11_p1 + Num12_p1 + Node133_p1) + 1, -1, 0, 1, Num13_p1 - Node132_p1 + Num14_p1 + Node151_p1,
			C4, C4, C3_1, C1_1, C3_1, C2_1, Coef, Coef_location);
			break;
		case (Node133_p1) ://2 lei xia
			kind2_node(j, i, -(Num9_p1 - Node93_p1 - 1 + Num10_p1 + Num11_p1 + Num12_p1 + Node133_p1), -1, 0, 1, 2,
			Num13_p1 - Node133_p1 + Num14_p1 + Num15_p1 + Num16_p1 + Node174_p1, B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1) < (Node132_p1))// 4dx
			{
				inner_5node(j, i, 4, -(Num9_p1 - Node_9_p1 + Num10_p1 + Num11_p1 + Num12_p1 + Node_13_p1),
					Num13_p1 - Node_13_p1 + Num14_p1 + Num15_p1 + Num16_p1 + Node_17_p1, Coef, Coef_location);
			}
			else if ((i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1) > (Node133_p1))// 4dx!
			{
				inner_5node(j, i, 4, -(Num9_p1 - Node93_p1 - 1 + Num10_p1 + Num11_p1 + Num12_p1 + Node133_p1),
					Num13_p1 - Node133_p1 + Num14_p1 + Num15_p1 + Num16_p1 + Node174_p1, Coef, Coef_location);
			}
			else if (((i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1) > (Node132_p1 + 1)) && ((i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1) < (Node133_p1 - 1)))
			{
				inner_5node(j, i, 2, -(Num11_p1 - Node112_p1 - 1 + Num12_p1 + Node132_p1 + 2),
					Num13_p1 - Node132_p1 + Num14_p1 + Node151_p1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post1_q1 + Num10_p1 + Num11_p1 + Num12_p1 + Node133_p1 + 3);
		i<(M_post1_q1 + Num10_p1 + Num11_p1 + Num12_p1 + Num13_p1); i++)// circle num13 :  95 node!  part3 !
	{
		switch (i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1)
		{
		case (Num13_p1 - 1) :
			boundary_node(j, i, 2, -Num10_p1 - Num11_p1 - Num12_p1 - Num13_p1, Num14_p1 + Num15_p1 + Num16_p1 + Num17_p1, Coef, Coef_location);
			break;
		case (Node133_p1 + 3) ://2 lei shang 1j
			kind2_node(j, i, -(Num9_p1 - Node93_p1 - 1 + Num10_p1 + Num11_p1 + Num12_p1 + Node133_p1), -2, -1, 0, 1,
			Num13_p1 - Node133_p1 + Num14_p1 + Num15_p1 + Num16_p1 + Node174_p1, B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node133_p1 + 7) :// 2 lei 2j shang
			kind2_node(j, i, -(Num11_p1 - Node113_p1 + Num12_p1 + Node133_p1 + 4), -2, -1, 0, 1, Num13_p1 - Node133_p1 - 4 + Num14_p1 + Node153_p1,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node134_p1:// inner 3 node
			inner_3node(j, i, -1, 0, (Num13_p1 - Node_1313_p1 + 4 + Num14_p1 - 1) - N_c1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node134_p1 + 1) :// inner 3 node
			inner_3node(j, i, 0, 1, (Num13_p1 - Node_1313_p1 + 4 + Num14_p1 - 1), k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Node_1313_p1 - 3) ://2 lei xia 2j
			kind2_node(j, i, -(Num12_p1 + Node_1313_p1 + 1), -1, 0, 1, 2, Num13_p1 - Node_1313_p1 + Num14_p1 + Num15_p1 - 1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_1313_p1 + 1) ://2 lei xia 1j
			kind2_node(j, i, -Num10_p1 - Num11_p1 - Num12_p1 - Num13_p1, -1, 0, 1, 2, Num14_p1 + Num15_p1 + Num16_p1 + Num17_p1,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1) < (Node133_p1 + 7))// 2dx!
			{
				inner_5node(j, i, 2, -(Num11_p1 - Node113_p1 + Num12_p1 + Node133_p1 + 4), Num13_p1 - Node133_p1 - 4 + Num14_p1 + Node153_p1, Coef, Coef_location);
			}
			else if (((i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1) > (Node133_p1 + 7)) && ((i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1) < Node134_p1))// 1dx
			{
				inner_5node(j, i, 1, -(Node_1313_p1 - 4 + 1) + 2, (Num13_p1 - Node_1313_p1 + 4 + Num14_p1 - 1) - N_c1, Coef, Coef_location);
			}
			else if (((i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1) > (Node134_p1 + 1)) && ((i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1) < (Node_1313_p1 - 3)))// 1dx
			{
				inner_5node(j, i, 1, -(Node_1313_p1 - 4 + 1), (Num13_p1 - Node_1313_p1 + 4 + Num14_p1 - 1), Coef, Coef_location);
			}
			else if (((i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1) > (Node_1313_p1 - 3)) && ((i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1) < (Node_1313_p1 + 1)))// 2dx
			{
				inner_5node(j, i, 2, -(Num12_p1 + Node_1313_p1 + 1), Num13_p1 - Node_1313_p1 + Num14_p1 + Num15_p1 - 1, Coef, Coef_location);
			}
			else// 4dx!
			{
				inner_5node(j, i, 4, -Num10_p1 - Num11_p1 - Num12_p1 - Num13_p1, Num14_p1 + Num15_p1 + Num16_p1 + Num17_p1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (M_post1_q1 + Num10_p1 + Num11_p1 + Num12_p1 + Num13_p1);
		i<(M_post1_q1 + Num10_p1 + Num11_p1 + Num12_p1 + Num13_p1 + Num14_p1 / 2); i++)// num 14 :  30 node! part1
	{
		switch (i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1 - Num13_p1)
		{
		case 0:// 1 lei 2j shang
			kind1_node(j, i, -(Num13_p1 - N_d0 - 2) - 2, -(Num13_p1 - N_d0 - 2) - 1, -(Num13_p1 - N_d0 - 2), 0, 1, Num14_p1 + 4 - 2,
				Num14_p1 + 4 - 1, Num14_p1 + 4, A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case (Num14_p1 / 2 - 1) :// 1 lei 2j xia
			kind1_node(j, i, -(Num13_p1 - N_d0 - 2) - N_c1, -(Num13_p1 - N_d0 - 2) - N_c1 + 1, -(Num13_p1 - N_d0 - 2) - N_c1 + 2, -1, 0,
			Num14_p1 + 4, Num14_p1 + 4 + 1, Num14_p1 + 4 + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1 - Num13_p1) <= Node141_p1)// 1dx
			{
				inner_5node(j, i, 1, -(Num13_p1 - N_d0 - 2), Num14_p1 + 4, Coef, Coef_location);
			}
			else if ((i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1 - Num13_p1) >= Node142_p1) // 1dx
			{
				inner_5node(j, i, 1, -(Num13_p1 - N_d0 - 2) - N_c1, Num14_p1 + 4, Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_4node(j, i, -1, 0, 1, Num14_p1 + 4, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post1_q1 + Num10_p1 + Num11_p1 + Num12_p1 + Num13_p1 + Num14_p1 / 2);
		i<(M_post1_q1 + Num10_p1 + Num11_p1 + Num12_p1 + Num13_p1 + Num14_p1); i++)// num14 :  30 node! part2
	{
		switch (i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1 - Num13_p1)
		{
		case (Num14_p1 / 2) :// 1 lei 2j shang
			kind1_node(j, i, -(Num13_p1 - Node_1313_p1 + 4 + Num14_p1 - 1) + N_c1 - 2, -(Num13_p1 - Node_1313_p1 + 4 + Num14_p1 - 1) + N_c1 - 1,
			-(Num13_p1 - Node_1313_p1 + 4 + Num14_p1 - 1) + N_c1, 0, 1, (Num15_p1 - 5 + 1) - 2, (Num15_p1 - 5 + 1) - 1, (Num15_p1 - 5 + 1),
			A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case (Num14_p1 - 1) :// 1 lei 2j xia
			kind1_node(j, i, -(Num13_p1 - Node_1313_p1 + 4 + Num14_p1 - 1), -(Num13_p1 - Node_1313_p1 + 4 + Num14_p1 - 1) + 1,
			-(Num13_p1 - Node_1313_p1 + 4 + Num14_p1 - 1) + 2, -1, 0, (Num15_p1 - 5 + 1), (Num15_p1 - 5 + 1) + 1, (Num15_p1 - 5 + 1) + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1 - Num13_p1) <= Node143_p1)// 1dx
			{
				inner_5node(j, i, 1, -(Num13_p1 - Node_1313_p1 + 4 + Num14_p1 - 1) + N_c1, (Num15_p1 - 5 + 1), Coef, Coef_location);
			}
			else if ((i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1 - Num13_p1) >= Node144_p1) // 1dx
			{
				inner_5node(j, i, 1, -(Num13_p1 - Node_1313_p1 + 4 + Num14_p1 - 1), (Num15_p1 - 5 + 1), Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_4node(j, i, -1, 0, 1, (Num15_p1 - 5 + 1), 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (M_post1_q1 + Num10_p1 + Num11_p1 + Num12_p1 + Num13_p1 + Num14_p1);
		i<(M_post1_q1 + Num10_p1 + Num11_p1 + Num12_p1 + Num13_p1 + Num14_p1 + Node151_p1 + 1); i++)// num 15 :  55 node! part1
	{
		switch (i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1 - Num13_p1 - Num14_p1)
		{
		case 0:// 1 lei 1j shang
			kind1_node(j, i, -(Num13_p1 - N_d0 + 2 + Num14_p1) - 2, -(Num13_p1 - N_d0 + 2 + Num14_p1) - 1, -(Num13_p1 - N_d0 + 2 + Num14_p1),
				0, 1, Num15_p1 + Num16_p1 + N_d0 - 2 - 2, Num15_p1 + Num16_p1 + N_d0 - 2 - 1, Num15_p1 + Num16_p1 + N_d0 - 2,
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case 3:// 2 lei 2j shang
			kind2_node(j, i, -(Num13_p1 - N_d0 + 2 + Num14_p1), -2, -1, 0, 1, Num15_p1 + Num16_p1 + N_d0 - 2,
				B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case 4:// 3 lei 2j shang
			kind3_node(j, i, -(Num14_p1 + 4), -1, 0, 1, Num15_p1 + Num16_p1 + N_d0 - 2 - 1, Num15_p1 + Num16_p1 + N_d0 - 2,
				C2_2, C3_2, C1_2, C3_2, C4, C4, Coef, Coef_location);
			break;
		case (Node151_p1 - 4) :// 3 lei 2j xia
			kind3_node(j, i, -(Num14_p1 + 4), -1, 0, 1, Num15_p1 - Node151_p1 + Num16_p1 + Node_17_p1, Num15_p1 - Node151_p1 + Num16_p1 + Node_17_p1 + 1,
			C2_2, C3_2, C1_2, C3_2, C4, C4, Coef, Coef_location);
			break;
		case (Node151_p1 - 3) :// 2 lei 2j xia
			kind2_node(j, i, -(Num13_p1 - Node_13_p1 + Num14_p1 + Node151_p1), -1, 0, 1, 2, Num15_p1 - Node151_p1 + Num16_p1 + Node_17_p1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case Node151_p1:// 1 lei 1j xia
			kind1_node(j, i, -(Num13_p1 - Node_13_p1 + Num14_p1 + Node151_p1), -(Num13_p1 - Node_13_p1 + Num14_p1 + Node151_p1) + 1,
				-(Num13_p1 - Node_13_p1 + Num14_p1 + Node151_p1) + 2, -1, 0, Num15_p1 - Node151_p1 + Num16_p1 + Node_17_p1,
				Num15_p1 - Node151_p1 + Num16_p1 + Node_17_p1 + 1, Num15_p1 - Node151_p1 + Num16_p1 + Node_17_p1 + 2,
				A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1 - Num13_p1 - Num14_p1) < 3)// 2dx
			{
				inner_5node(j, i, 2, -(Num13_p1 - N_d0 + 2 + Num14_p1), Num15_p1 + Num16_p1 + N_d0 - 2, Coef, Coef_location);
			}
			else if (((i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1 - Num13_p1 - Num14_p1) > 4) && ((i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1 - Num13_p1 - Num14_p1) < (Node151_p1 - 4)))// 1dx
			{
				inner_5node(j, i, 1, -(Num14_p1 + 4), Num15_p1 - 5, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, i, 2, -(Num13_p1 - Node_13_p1 + Num14_p1 + Node151_p1), Num15_p1 - Node151_p1 + Num16_p1 + Node_17_p1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post1_q1 + Num10_p1 + Num11_p1 + Num12_p1 + Num13_p1 + Num14_p1 + Node151_p1 + 1);
		i<(M_post1_q1 + Num10_p1 + Num11_p1 + Num12_p1 + Num13_p1 + Num14_p1 + Node153_p1); i++)// num 15 :  55 node! part2
	{
		switch (i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1 - Num13_p1 - Num14_p1)
		{
		case (Node151_p1 + 1) :// 1 lei 1j shang
			kind1_node(j, i, -(Num13_p1 - Node132_p1 + Num14_p1 + Node151_p1) - 2, -(Num13_p1 - Node132_p1 + Num14_p1 + Node151_p1) - 1,
			-(Num13_p1 - Node132_p1 + Num14_p1 + Node151_p1), 0, 1, (Num15_p1 - Node151_p1 + Num16_p1 + Node173_p1) - 2,
			(Num15_p1 - Node151_p1 + Num16_p1 + Node173_p1) - 1, (Num15_p1 - Node151_p1 + Num16_p1 + Node173_p1),
			A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case (Node152_p1) :// 2 lei 2j zuo
			kind2_node(j, i, -(Num13_p1 - Node132_p1 + Num14_p1 + Node151_p1) - (Num11_p1 - Node112_p1 - 1 + Num12_p1 + Node132_p1 + 2),
			-(Num13_p1 - Node132_p1 + Num14_p1 + Node151_p1), -1, 0, 1, (Num15_p1 - Node152_p1 + Node161_p1 + 1),
			B5, B4, B3, B1_2, B3, B2, Coef, Coef_location);
			break;
		case (Node152_p1 + 1) :// 2 lei 2j zuo
			kind2_node(j, i, -(Num13_p1 - Node132_p1 + Num14_p1 + Node151_p1) - (Num11_p1 - Node112_p1 - 1 + Num12_p1 + Node132_p1 + 2),
			-(Num13_p1 - Node132_p1 + Num14_p1 + Node151_p1), -1, 0, 1, (Num15_p1 - Node152_p1 + Node161_p1 + 1) + 1,
			B5, B4, B3, B1_2, B3, B2, Coef, Coef_location);
			break;
		case (Node152_p1 + 2) :// 2 lei 2j zuo
			kind2_node(j, i, -(Num13_p1 - Node132_p1 + Num14_p1 + Node151_p1) - (Num11_p1 - Node112_p1 - 1 + Num12_p1 + Node132_p1 + 2),
			-(Num13_p1 - Node132_p1 + Num14_p1 + Node151_p1), -1, 0, 1, (Num15_p1 - Node152_p1 + Node161_p1 + 1) + 2,
			B5, B4, B3, B1_2, B3, B2, Coef, Coef_location);
			break;
		case (Node153_p1 - 1) :// 1 lei 1j xia
			kind1_node(j, i, -(Num13_p1 - Node132_p1 + Num14_p1 + Node151_p1), -(Num13_p1 - Node132_p1 + Num14_p1 + Node151_p1) + 1,
			-(Num13_p1 - Node132_p1 + Num14_p1 + Node151_p1) + 2, -1, 0, Num15_p1 - Node153_p1 + Num16_p1 + Node174_p1,
			Num15_p1 - Node153_p1 + Num16_p1 + Node174_p1 + 1, Num15_p1 - Node153_p1 + Num16_p1 + Node174_p1 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1 - Num13_p1 - Num14_p1) < Node152_p1)// 1dx
			{
				inner_5node(j, i, 2, -(Num13_p1 - Node132_p1 + Num14_p1 + Node151_p1), (Num15_p1 - Node151_p1 + Num16_p1 + Node173_p1), Coef, Coef_location);
			}
			else if ((i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1 - Num13_p1 - Num14_p1) > (Node152_p1 + 2))// 1dx
			{
				inner_5node(j, i, 2, -(Num13_p1 - Node132_p1 + Num14_p1 + Node151_p1), Num15_p1 - Node153_p1 + Num16_p1 + Node174_p1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post1_q1 + Num10_p1 + Num11_p1 + Num12_p1 + Num13_p1 + Num14_p1 + Node153_p1);
		i<(M_post1_q1 + Num10_p1 + Num11_p1 + Num12_p1 + Num13_p1 + Num14_p1 + Num15_p1); i++)// num 15 :  55 node! part3
	{
		switch (i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1 - Num13_p1 - Num14_p1)
		{
		case Node153_p1:// 1 lei 1j shang
			kind1_node(j, i, -(Num13_p1 - Node133_p1 - 4 + Num14_p1 + Node153_p1) - 2, -(Num13_p1 - Node133_p1 - 4 + Num14_p1 + Node153_p1) - 1,
				-(Num13_p1 - Node133_p1 - 4 + Num14_p1 + Node153_p1), 0, 1, (Num15_p1 - Node153_p1 + Num16_p1 + Node174_p1 + 4) - 2,
				(Num15_p1 - Node153_p1 + Num16_p1 + Node174_p1 + 4) - 1, (Num15_p1 - Node153_p1 + Num16_p1 + Node174_p1 + 4),
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case (Node153_p1 + 3) :// 2 lei 2j shang
			kind2_node(j, i, -(Num13_p1 - Node133_p1 - 4 + Num14_p1 + Node153_p1), -2, -1, 0, 1, (Num15_p1 - Node153_p1 + Num16_p1 + Node174_p1 + 4),
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case (Node153_p1 + 4) :// 3 lei 2j shang
			kind3_node(j, i, -(Num15_p1 - 5 + 1), -1, 0, 1,
			(Num15_p1 - Node153_p1 + Num16_p1 + Node174_p1 + 4) - 1, (Num15_p1 - Node153_p1 + Num16_p1 + Node174_p1 + 4),
			C2_2, C3_2, C1_2, C3_2, C4, C4, Coef, Coef_location);
			break;
		case (Num15_p1 - 5) :// 3 lei 2j xia
			kind3_node(j, i, -(Num15_p1 - 5 + 1), -1, 0, 1, (1 + Num16_p1 + Node_1717_p1), (1 + Num16_p1 + Node_1717_p1) + 1,
			C2_2, C3_2, C1_2, C3_2, C4, C4, Coef, Coef_location);
			break;
		case (Num15_p1 - 4) :// 2 lei 2j xia
			kind2_node(j, i, -(Num13_p1 - Node_1313_p1 + Num14_p1 + Num15_p1 - 1), -1, 0, 1, 2, (1 + Num16_p1 + Node_1717_p1),
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num15_p1 - 1) :// 1 lei 1j xia
			kind1_node(j, i, -(Num13_p1 - Node_1313_p1 + Num14_p1 + Num15_p1 - 1), -(Num13_p1 - Node_1313_p1 + Num14_p1 + Num15_p1 - 1) + 1,
			-(Num13_p1 - Node_1313_p1 + Num14_p1 + Num15_p1 - 1) + 2, -1, 0, (1 + Num16_p1 + Node_1717_p1),
			(1 + Num16_p1 + Node_1717_p1) + 1, (1 + Num16_p1 + Node_1717_p1) + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1 - Num13_p1 - Num14_p1) < (Node153_p1 + 3))// 2dx
			{
				inner_5node(j, i, 2, -(Num13_p1 - Node133_p1 - 4 + Num14_p1 + Node153_p1), (Num15_p1 - Node153_p1 + Num16_p1 + Node174_p1 + 4), Coef, Coef_location);
			}
			else if (((i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1 - Num13_p1 - Num14_p1) > (Node153_p1 + 4)) && ((i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1 - Num13_p1 - Num14_p1) < (Num15_p1 - 4)))// 1dx
			{
				inner_5node(j, i, 1, -(Num15_p1 - 5 + 1), 6 + Num16_p1 - 1, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, i, 2, -(Num13_p1 - Node_1313_p1 + Num14_p1 + Num15_p1 - 1), (1 + Num16_p1 + Node_1717_p1), Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (M_post1_q1 + Num10_p1 + Num11_p1 + Num12_p1 + Num13_p1 + Num14_p1 + Num15_p1);
		i<(M_post1_q1 + Num10_p1 + Num11_p1 + Num12_p1 + Num13_p1 + Num14_p1 + Num15_p1 + Node161_p1 + 1); i++)// num 16 :  31 node! part1
	{
		switch (i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1 - Num13_p1 - Num14_p1 - Num15_p1)
		{
		case 0:// 3 lei 2j shang
			kind3_node(j, i, -(Num15_p1 - 3), -(Num15_p1 - 3) + 2, 0, 1, Num16_p1 + Node171_p1 - 1, Num16_p1 + Node171_p1,
				C4, C3_2, C1_2, C2_2, C4, C3_2, Coef, Coef_location);
			break;
		case Node161_p1://node2: 3 lei 2j xia
			kind3_node(j, i, -(Num15_p1 - 5), -(Num15_p1 - 5) + 2, -1, 0, (Num16_p1 - Node161_p1 + Node172_p1),
				(Num16_p1 - Node161_p1 + Node172_p1) + 1, C3_2, C4, C2_2, C1_2, C3_2, C4, Coef, Coef_location);
			break;
		default:
			if ((i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1 - Num13_p1 - Num14_p1 - Num15_p1) % 2 == 0)// 1dx
			{
				inner_5node(j, i, 1, -(Num15_p1 - 5),
					(Num16_p1 - (i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1 - Num13_p1 - Num14_p1 - Num15_p1) + Node171_p1 + (i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1 - Num13_p1 - Num14_p1 - Num15_p1) / 2), Coef, Coef_location);
			}
			else//1 lei  2j left
			{
				kind1_node(j, i, -(Num15_p1 - 5), -1, 0, 1,
					(Num16_p1 - (i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1 - Num13_p1 - Num14_p1 - Num15_p1) + Node171_p1 + (i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1 - Num13_p1 - Num14_p1 - Num15_p1) / 2),
					(Num16_p1 - (i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1 - Num13_p1 - Num14_p1 - Num15_p1) + Node171_p1 + (i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1 - Num13_p1 - Num14_p1 - Num15_p1) / 2) + 1,
					(Num16_p1 - (i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1 - Num13_p1 - Num14_p1 - Num15_p1) + Node171_p1 + (i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1 - Num13_p1 - Num14_p1 - Num15_p1) / 2) + (Num17_p1 - Node171_p1 + 4 + Num18_p1),
					(Num16_p1 - (i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1 - Num13_p1 - Num14_p1 - Num15_p1) + Node171_p1 + (i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1 - Num13_p1 - Num14_p1 - Num15_p1) / 2) + 1 + (Num17_p1 - Node171_p1 + 4 + Num18_p1),
					A4_2, A3_2, A5_2, A3_2, A1_2, A1_2, A2_2, A2_2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post1_q1 + Num10_p1 + Num11_p1 + Num12_p1 + Num13_p1 + Num14_p1 + Num15_p1 + Node161_p1 + 1);
		i<(M_post1_q1 + Num10_p1 + Num11_p1 + Num12_p1 + Num13_p1 + Num14_p1 + Num15_p1 + Node162_p1); i++)// num 16 :  31 node! part2
	{
		switch (i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1 - Num13_p1 - Num14_p1 - Num15_p1)
		{
		case (Node161_p1 + 1) :// 3 lei 2j shang !!! 1
			kind3_node(j, i, -(Num15_p1 - Node152_p1 + Node161_p1 + 1) - 1, -(Num15_p1 - Node152_p1 + Node161_p1 + 1), 0, 1,
			(Num16_p1 - Node161_p1 - 1 + Node173_p1 + 5) - 2, (Num16_p1 - Node161_p1 - 1 + Node173_p1 + 5),
			C4, C3_2, C1_2, C2_2, C4, C3_2, Coef, Coef_location);
			break;
		case (Node161_p1 + 2) :// 1 lei 2j zuo !!! 2
			kind1_node(j, i, -(Num15_p1 - Node152_p1 + Node161_p1 + 1) - 1 - (Num13_p1 - Node132_p1 + Num14_p1 + Node151_p1),
			-(Num15_p1 - Node152_p1 + Node161_p1 + 1) - (Num13_p1 - Node132_p1 + Num14_p1 + Node151_p1),
			-(Num15_p1 - Node152_p1 + Node161_p1 + 1) - 1, -(Num15_p1 - Node152_p1 + Node161_p1 + 1),
			-1, 0, 1, (Num16_p1 - Node161_p1 - 1 + Node173_p1 + 5), A2_2, A2_2, A1_2, A1_2, A3_2, A5_2, A3_2, A4_2, Coef, Coef_location);
			break;
		case (Node161_p1 + 3) :// inner 5 node !!! 3
			inner_5node(j, i, 1, -(Num15_p1 - Node152_p1 + Node161_p1 + 1) - 1, (Num16_p1 - Node161_p1 - 1 + Node173_p1 + 5), Coef, Coef_location);
			break;
		case (Node162_p1 - 2) :// 1 lei 2j zuo !!! 4
			kind1_node(j, i, -(Num15_p1 - Node152_p1 + Node161_p1 + 1) - 2 - (Num13_p1 - Node132_p1 + Num14_p1 + Node151_p1),
			-(Num15_p1 - Node152_p1 + Node161_p1 + 1) - 1 - (Num13_p1 - Node132_p1 + Num14_p1 + Node151_p1),
			-(Num15_p1 - Node152_p1 + Node161_p1 + 1) - 2, -(Num15_p1 - Node152_p1 + Node161_p1 + 1) - 1,
			-1, 0, 1, (Num16_p1 - Node161_p1 - 1 + Node173_p1 + 5), A2_2, A2_2, A1_2, A1_2, A3_2, A5_2, A3_2, A4_2, Coef, Coef_location);
			break;
		case (Node162_p1 - 1) :// 3 lei 2j xia !!! 5
			kind3_node(j, i, -(Num15_p1 - Node152_p1 + Node161_p1 + 1) - 2, -(Num15_p1 - Node152_p1 + Node161_p1 + 1) - 2 + 1, -1, 0,
			(Num16_p1 - Node161_p1 - 1 + Node173_p1 + 5), (Num16_p1 - Node161_p1 - 1 + Node173_p1 + 5) + 2,
			C3_2, C4, C2_2, C1_2, C3_2, C4, Coef, Coef_location);
			break;
		}
		j = j + 1;
	}
	for (i = (M_post1_q1 + Num10_p1 + Num11_p1 + Num12_p1 + Num13_p1 + Num14_p1 + Num15_p1 + Node162_p1);
		i<(M_post1_q1 + Num10_p1 + Num11_p1 + Num12_p1 + Num13_p1 + Num14_p1 + Num15_p1 + Num16_p1); i++)// num 16 :  31 node! part3
	{
		switch (i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1 - Num13_p1 - Num14_p1 - Num15_p1 - Node162_p1)
		{
		case 0:// 3 lei 2j shang
			kind3_node(j, i, -(6 + Num16_p1 - 1) - 2, -(6 + Num16_p1 - 1), 0, 1, Num16_p1 - Node162_p1 + Node175_p1 - 1,
				Num16_p1 - Node162_p1 + Node175_p1, C4, C3_2, C1_2, C2_2, C4, C3_2, Coef, Coef_location);
			break;
		case Node161_p1:// 3 lei 2j xia
			kind3_node(j, i, -(6 + Num16_p1 - 1), -(6 + Num16_p1 - 1) + 2, -1, 0, 1 + Node176_p1, 1 + Node176_p1 + 1,
				C3_2, C4, C2_2, C1_2, C3_2, C4, Coef, Coef_location);
			break;
		default:
			if ((i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1 - Num13_p1 - Num14_p1 - Num15_p1 - Node162_p1) % 2 == 0)// 1dx
			{
				inner_5node(j, i, 1, -(6 + Num16_p1 - 1),
					Num16_p1 - (i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1 - Num13_p1 - Num14_p1 - Num15_p1) + Node175_p1 + (i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1 - Num13_p1 - Num14_p1 - Num15_p1 - Node162_p1) / 2,
					Coef, Coef_location);
			}
			else//1 lei 2j left
			{
				kind1_node(j, i, -(6 + Num16_p1 - 1), -1, 0, 1,
					Num16_p1 - (i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1 - Num13_p1 - Num14_p1 - Num15_p1) + Node175_p1 + (i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1 - Num13_p1 - Num14_p1 - Num15_p1 - Node162_p1) / 2,
					Num16_p1 - (i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1 - Num13_p1 - Num14_p1 - Num15_p1) + Node175_p1 + (i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1 - Num13_p1 - Num14_p1 - Num15_p1 - Node162_p1) / 2 + 1,
					Num16_p1 - (i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1 - Num13_p1 - Num14_p1 - Num15_p1) + Node175_p1 + (i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1 - Num13_p1 - Num14_p1 - Num15_p1 - Node162_p1) / 2 + (Num17_p1 - Node_1717_p1 + Num18_p1 + Num19_p1 - 1),
					Num16_p1 - (i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1 - Num13_p1 - Num14_p1 - Num15_p1) + Node175_p1 + (i - M_post1_q1 - Num10_p1 - Num11_p1 - Num12_p1 - Num13_p1 - Num14_p1 - Num15_p1 - Node162_p1) / 2 + 1 + (Num17_p1 - Node_1717_p1 + Num18_p1 + Num19_p1 - 1),
					A4_2, A3_2, A5_2, A3_2, A1_2, A1_2, A2_2, A2_2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (M_post1_q2 - Num19_p1 - Num18_p1 - Num17_p1);
		i<(M_post1_q2 - Num19_p1 - Num18_p1 - Num17_p1 + Node_17_p1 + 2); i++)// circle num17 :  83 node!  part1 !
	{
		switch (i - M_post1_q2 + Num19_p1 + Num18_p1 + Num17_p1)
		{
		case 0:
			boundary_node(j, i, 1, -Num13_p1 - Num14_p1 - Num15_p1 - Num16_p1, Num17_p1 + Num18_p1 + Num19_p1 + Num18_p1, Coef, Coef_location);
			break;
		case (N_d0 - 3) ://2 lei shang 1j
			kind2_node(j, i, -Num13_p1 - Num14_p1 - Num15_p1 - Num16_p1, -2, -1, 0, 1, Num17_p1 + Num18_p1 + Num19_p1 + Num18_p1,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_17_p1 + 1) ://2 lei xia 1j
			kind2_node(j, i, -(Num13_p1 - Node_13_p1 + Num14_p1 + Num15_p1 + Num16_p1 + Node_17_p1), -1, 0, 1, 2,
			(Num17_p1 - Node_17_p1 + Num18_p1 + Num19_p1 + Num18_p1 + Node_17_p1), B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - M_post1_q2 + Num19_p1 + Num18_p1 + Num17_p1) < (N_d0 - 3))// 4dx!
			{
				inner_5node(j, i, 4, -Num13_p1 - Num14_p1 - Num15_p1 - Num16_p1, Num17_p1 + Num18_p1 + Num19_p1 + Num18_p1, Coef, Coef_location);
			}
			else if (((i - M_post1_q2 + Num19_p1 + Num18_p1 + Num17_p1) > (N_d0 - 3)) && ((i - M_post1_q2 + Num19_p1 + Num18_p1 + Num17_p1) < Node171_p1))// 2dx!
			{
				inner_5node(j, i, 2, -(Num15_p1 + Num16_p1 + N_d0 - 2), Num17_p1 - N_d0 + 2 + Num18_p1, Coef, Coef_location);
			}
			else if (((i - M_post1_q2 + Num19_p1 + Num18_p1 + Num17_p1) >= Node171_p1) && ((i - M_post1_q2 + Num19_p1 + Num18_p1 + Num17_p1) <= Node172_p1))// 2 lei right 2j !
			{
				kind2_node(j, i, -((i - M_post1_q2 + Num19_p1 + Num18_p1 + Num17_p1) + Num16_p1 - (i - M_post1_q2 + Num19_p1 + Num18_p1 + Num17_p1 - Node171_p1) * 2),
					-1, 0, 1, Num17_p1 - N_d0 + 2 + Num18_p1, Num17_p1 + Num18_p1 + Num19_p1 + Num18_p1, B2, B3, B1_2, B3, B4, B5, Coef, Coef_location);
			}
			else// 2dx!
			{
				inner_5node(j, i, 2, -(Num15_p1 - Node151_p1 + Num16_p1 + Node_17_p1), Num17_p1 - N_d0 + 2 + Num18_p1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post1_q2 - Num19_p1 - Num18_p1 - Num17_p1 + Node_17_p1 + 2);
		i<(M_post1_q2 - Num19_p1 - Num18_p1 - Num17_p1 + Node174_p1 + 3); i++)// circle num17 :  83 node!  part2 !
	{
		switch (i - M_post1_q2 + Num19_p1 + Num18_p1 + Num17_p1)
		{
		case Node173_p1:// 2lei shang 1j
			kind2_node(j, i, -(Num13_p1 - Node_13_p1 + Num14_p1 + Num15_p1 + Num16_p1 + Node_17_p1), -2, -1, 0, 1,
				(Num17_p1 - Node_17_p1 + Num18_p1 + Num19_p1 + Num18_p1 + Node_17_p1), B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node173_p1 + 3) :// 2lei shang 2j
			kind2_node(j, i, -(Num15_p1 - Node151_p1 + Num16_p1 + Node173_p1), -2, -1, 0, 1, (Num17_p1 - Node173_p1 + Num18_p1 + Node191_p1),
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case (Node173_p1 + 4) :// 3lei zuo 2j
			kind3_node(j, i, -(Num15_p1 - Node151_p1 + Num16_p1 + Node173_p1) - 1, -(Num15_p1 - Node151_p1 + Num16_p1 + Node173_p1),
			-1, 0, 1, (Num17_p1 - Node173_p1 - 4), C4, C4, C3_2, C1_2, C3_2, C2_2, Coef, Coef_location);
			break;
		case (Node173_p1 + 5) ://inner 5 node
			inner_5node(j, i, 1, -(Num16_p1 - Node161_p1 - 1 + Node173_p1 + 5), (Num17_p1 - Node173_p1 - 4), Coef, Coef_location);
			break;
			//************************************************
		case (Node174_p1 - 5) ://inner 5 node
			inner_5node(j, i, 1, -(Num16_p1 - Node161_p1 - 1 + Node173_p1 + 5), (Num17_p1 - Node174_p1 + 4 + 3), Coef, Coef_location);
			break;
		case (Node174_p1 - 4) :// 3lei zuo 2j
			kind3_node(j, i, -(Num15_p1 - Node153_p1 + Num16_p1 + Node174_p1), -(Num15_p1 - Node153_p1 + Num16_p1 + Node174_p1) + 1,
			-1, 0, 1, (Num17_p1 - Node174_p1 + 4 + 3), C4, C4, C3_2, C1_2, C3_2, C2_2, Coef, Coef_location);
			break;
		case (Node174_p1 - 3) :// 2lei xia 2j
			kind2_node(j, i, -(Num15_p1 - Node153_p1 + Num16_p1 + Node174_p1), -1, 0, 1, 2, Num17_p1 - Node174_p1 + Num18_p1 + Node193_p1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node174_p1) :// 2lei xia 1j
			kind2_node(j, i, -(Num13_p1 - Node133_p1 + Num14_p1 + Num15_p1 + Num16_p1 + Node174_p1), -1, 0, 1, 2,
			(Num17_p1 + Num18_p1 + Num19_p1 + Num18_p1), B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - M_post1_q2 + Num19_p1 + Num18_p1 + Num17_p1) < (Node173_p1))  // 2dx
			{
				inner_5node(j, i, 4, -(Num13_p1 - Node_13_p1 + Num14_p1 + Num15_p1 + Num16_p1 + Node_17_p1),
					(Num17_p1 - Node_17_p1 + Num18_p1 + Num19_p1 + Num18_p1 + Node_17_p1), Coef, Coef_location);
			}
			else if (((i - M_post1_q2 + Num19_p1 + Num18_p1 + Num17_p1) > (Node173_p1)) && ((i - M_post1_q2 + Num19_p1 + Num18_p1 + Num17_p1) < (Node173_p1 + 3)))  // 1dx
			{
				inner_5node(j, i, 2, -(Num15_p1 - Node151_p1 + Num16_p1 + Node173_p1),
					(Num17_p1 - Node173_p1 + Num18_p1 + Node191_p1), Coef, Coef_location);
			}
			else if (((i - M_post1_q2 + Num19_p1 + Num18_p1 + Num17_p1) > (Node173_p1 + 5)) && ((i - M_post1_q2 + Num19_p1 + Num18_p1 + Num17_p1) < (Node174_p1 - 5)))  // 1dx
			{
				inner_4node(j, i, -(Num16_p1 - Node161_p1 - 1 + Node173_p1 + 5), -1, 0, 1, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			}
			else if (((i - M_post1_q2 + Num19_p1 + Num18_p1 + Num17_p1) > (Node174_p1 - 3)) && ((i - M_post1_q2 + Num19_p1 + Num18_p1 + Num17_p1) < (Node174_p1)))  // 1dx
			{
				inner_5node(j, i, 2, -(Num15_p1 - Node153_p1 + Num16_p1 + Node174_p1),
					Num17_p1 - Node174_p1 + Num18_p1 + Node193_p1, Coef, Coef_location);
			}
			else
			{
				inner_5node(j, i, 4, -(Num13_p1 - Node133_p1 + Num14_p1 + Num15_p1 + Num16_p1 + Node174_p1),
					(Num17_p1 + Num18_p1 + Num19_p1 + Num18_p1), Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post1_q2 - Num19_p1 - Num18_p1 - Num17_p1 + Node174_p1 + 3);
		i<(M_post1_q2 - Num19_p1 - Num18_p1); i++)// circle num17 :  83 node!  part3 !
	{
		switch (i - M_post1_q2 + Num19_p1 + Num18_p1 + Num17_p1)
		{
		case (Num17_p1 - 1) :
			boundary_node(j, i, 2, -Num14_p1 - Num15_p1 - Num16_p1 - Num17_p1, Num18_p1 + Num19_p1 + Num18_p1 + Num17_p1, Coef, Coef_location);
			break;
		case (Node174_p1 + 3) ://2 lei shang 1j
			kind2_node(j, i, -(Num13_p1 - Node133_p1 + Num14_p1 + Num15_p1 + Num16_p1 + Node174_p1), -2, -1, 0, 1,
			(Num17_p1 + Num18_p1 + Num19_p1 + Num18_p1), B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_1717_p1 + 1) ://2 lei xia 1j
			kind2_node(j, i, -Num14_p1 - Num15_p1 - Num16_p1 - Num17_p1, -1, 0, 1, 2, Num18_p1 + Num19_p1 + Num18_p1 + Num17_p1,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - M_post1_q2 + Num19_p1 + Num18_p1 + Num17_p1) < Node175_p1)// 2dx!
			{
				inner_5node(j, i, 2, -(Num15_p1 - Node153_p1 + Num16_p1 + Node175_p1 - 4), Num17_p1 - Node_1717_p1 + Num18_p1 + Num19_p1 - 1, Coef, Coef_location);
			}
			else if (((i - M_post1_q2 + Num19_p1 + Num18_p1 + Num17_p1) >= Node175_p1) && ((i - M_post1_q2 + Num19_p1 + Num18_p1 + Num17_p1) <= Node176_p1))// 2 lei left 2j !
			{
				kind2_node(j, i, -((i - M_post1_q2 + Num19_p1 + Num18_p1 + Num17_p1) + Num16_p1 - Node162_p1 - (i - M_post1_q2 + Num19_p1 + Num18_p1 + Num17_p1 - Node175_p1) * 2),
					-1, 0, 1, Num17_p1 - Node_1717_p1 + Num18_p1 + Num19_p1 - 1, Num17_p1 + Num18_p1 + Num19_p1 + Num18_p1, B2, B3, B1_2, B3, B4, B5, Coef, Coef_location);
			}
			else if (((i - M_post1_q2 + Num19_p1 + Num18_p1 + Num17_p1) > Node176_p1) && ((i - M_post1_q2 + Num19_p1 + Num18_p1 + Num17_p1) <= Node_1717_p1))// 2dx
			{
				inner_5node(j, i, 2, -(Num16_p1 + Node_1717_p1 + 1), Num17_p1 - Node_1717_p1 + Num18_p1 + Num19_p1 - 1, Coef, Coef_location);
			}
			else
			{
				inner_5node(j, i, 4, -Num14_p1 - Num15_p1 - Num16_p1 - Num17_p1, Num18_p1 + Num19_p1 + Num18_p1 + Num17_p1, Coef, Coef_location);// 4dx!
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (M_post1_q2 - Num19_p1 - Num18_p1); i<(M_post1_q2 - Num19_p1); i++)// num 18 :  4 node! 
	{
		switch (i - M_post1_q2 + Num19_p1 + Num18_p1)
		{
		case 0:// 1 lei shang 2j
			kind1_node(j, i, -(Num17_p1 - Node173_p1 - 4) - 2, -(Num17_p1 - Node173_p1 - 4) - 1, -(Num17_p1 - Node173_p1 - 4),
				0, 1, (Num18_p1 + Node191_p1 + 4) - 2, (Num18_p1 + Node191_p1 + 4) - 1, (Num18_p1 + Node191_p1 + 4),
				A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case 1://inner 4 node
			inner_4node(j, i, -(Num17_p1 - Node173_p1 - 4), -1, 0, (Num18_p1 + Node191_p1 + 4), 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case 2://inner 4 node
			inner_4node(j, i, -(Num17_p1 - Node173_p1 - 4) + 3, 0, 1, (Num18_p1 + Node191_p1 + 4), 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Num18_p1 - 1) :// 1 lei xia 2j
			kind1_node(j, i, -(Num17_p1 - Node173_p1 - 4) + 3, -(Num17_p1 - Node173_p1 - 4) + 3 + 1, -(Num17_p1 - Node173_p1 - 4) + 3 + 2,
			-1, 0, (Num18_p1 + Node191_p1 + 4), (Num18_p1 + Node191_p1 + 4) + 1, (Num18_p1 + Node191_p1 + 4) + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (M_post1_q2 - Num19_p1); i<(M_post1_q2 - Num19_p1 + Node191_p1 + 1); i++)// num 19 :  40 node! part1
	{
		switch (i - M_post1_q2 + Num19_p1)
		{
		case 0:// 1 lei shang
			kind1_node(j, i, -(Num17_p1 - N_d0 + 2 + Num18_p1) - 2, -(Num17_p1 - N_d0 + 2 + Num18_p1) - 1, -(Num17_p1 - N_d0 + 2 + Num18_p1),
				0, 1, (Num19_p1 + Num18_p1 + N_d0 - 2) - 2, (Num19_p1 + Num18_p1 + N_d0 - 2) - 1, (Num19_p1 + Num18_p1 + N_d0 - 2),
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case Node191_p1:// 1 lei xia
			kind1_node(j, i, -(Num17_p1 - N_d0 + 2 + Num18_p1), -(Num17_p1 - N_d0 + 2 + Num18_p1) + 1, -(Num17_p1 - N_d0 + 2 + Num18_p1) + 2,
				-1, 0, (Num19_p1 + Num18_p1 + N_d0 - 2), (Num19_p1 + Num18_p1 + N_d0 - 2) + 1, (Num19_p1 + Num18_p1 + N_d0 - 2) + 2,
				A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			inner_5node(j, i, 2, -(Num17_p1 - N_d0 + 2 + Num18_p1), (Num19_p1 + Num18_p1 + N_d0 - 2), Coef, Coef_location);
			break;
		}
		j = j + 1;
	}
	for (i = (M_post1_q2 - Num19_p1 + Node191_p1 + 1); i<(M_post1_q2 - Num19_p1 + Node193_p1); i++)// num 19 :  40 node! part2
	{
		switch (i - M_post1_q2 + Num19_p1)
		{
		case (Node191_p1 + 1) :// 1 lei 1j shang
			kind1_node(j, i, -(Num17_p1 - Node173_p1 + Num18_p1 + Node191_p1) - 2, -(Num17_p1 - Node173_p1 + Num18_p1 + Node191_p1) - 1,
			-(Num17_p1 - Node173_p1 + Num18_p1 + Node191_p1), 0, 1, (Num19_p1 - Node191_p1 + Num18_p1 + Node173_p1) - 2,
			(Num19_p1 - Node191_p1 + Num18_p1 + Node173_p1) - 1, (Num19_p1 - Node191_p1 + Num18_p1 + Node173_p1),
			A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case (Node191_p1 + 2) :// 2dx
			inner_5node(j, i, 2, -(Num17_p1 - Node173_p1 + Num18_p1 + Node191_p1), (Num19_p1 - Node191_p1 + Num18_p1 + Node173_p1), Coef, Coef_location);
			break;
		case (Node191_p1 + 3) :// 2 lei 2j shang
			kind2_node(j, i, -(Num17_p1 - Node173_p1 + Num18_p1 + Node191_p1), -2, -1, 0, 1,
			(Num19_p1 - Node191_p1 + Num18_p1 + Node173_p1), B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node192_p1://inner 4 node
			inner_4node(j, i, -(3 + Node192_p1), -1, 0, Num19_p1 - Node192_p1 + 1, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node192_p1 + 1) ://inner 4 node
			inner_4node(j, i, -(3 + Node192_p1), 0, 1, Num19_p1 - Node192_p1 + 1, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Node193_p1 - 3) :// 2 lei 2j xia
			kind2_node(j, i, -(Num17_p1 - Node174_p1 + Num18_p1 + Node193_p1), -1, 0, 1, 2, Num19_p1 - Node193_p1 + Num18_p1 + Node174_p1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node193_p1 - 2) :// 2dx
			inner_5node(j, i, 2, -(Num17_p1 - Node174_p1 + Num18_p1 + Node193_p1), Num19_p1 - Node193_p1 + Num18_p1 + Node174_p1, Coef, Coef_location);
			break;
		case (Node193_p1 - 1) :// 1 lei 1j xia
			kind1_node(j, i, -(Num17_p1 - Node174_p1 + Num18_p1 + Node193_p1), -(Num17_p1 - Node174_p1 + Num18_p1 + Node193_p1) + 1,
			-(Num17_p1 - Node174_p1 + Num18_p1 + Node193_p1) + 2, -1, 0, Num19_p1 - Node193_p1 + Num18_p1 + Node174_p1,
			Num19_p1 - Node193_p1 + Num18_p1 + Node174_p1 + 1, Num19_p1 - Node193_p1 + Num18_p1 + Node174_p1 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			inner_5node(j, i, 1, -(3 + Node192_p1), Num19_p1 - Node192_p1 + 1, Coef, Coef_location);
			break;
		}
		j = j + 1;
	}
	for (i = (M_post1_q2 - Num19_p1 + Node193_p1); i<(M_post1_q2); i++)// num 19 :  40 node! part3
	{
		switch (i - M_post1_q2 + Num19_p1)
		{
		case Node193_p1:// 1 lei shang
			kind1_node(j, i, -(Num17_p1 - Node_1717_p1 + Num18_p1 - 1 + Num19_p1) - 2, -(Num17_p1 - Node_1717_p1 + Num18_p1 - 1 + Num19_p1) - 1,
				-(Num17_p1 - Node_1717_p1 + Num18_p1 - 1 + Num19_p1), 0, 1, (1 + Num18_p1 + Node_1717_p1) - 2,
				(1 + Num18_p1 + Node_1717_p1) - 1, (1 + Num18_p1 + Node_1717_p1),
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case (Num19_p1 - 1) :// 1 lei xia
			kind1_node(j, i, -(Num17_p1 - Node_1717_p1 + Num18_p1 - 1 + Num19_p1), -(Num17_p1 - Node_1717_p1 + Num18_p1 - 1 + Num19_p1) + 1,
			-(Num17_p1 - Node_1717_p1 + Num18_p1 - 1 + Num19_p1) + 2, -1, 0, (1 + Num18_p1 + Node_1717_p1),
			(1 + Num18_p1 + Node_1717_p1) + 1, (1 + Num18_p1 + Node_1717_p1) + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			inner_5node(j, i, 2, -(Num17_p1 - Node_1717_p1 + Num18_p1 - 1 + Num19_p1), (1 + Num18_p1 + Node_1717_p1), Coef, Coef_location);
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!!!
	for (i = (M_post1_q2); i<(M_post1_q2 + Num18_p1); i++)// num 20(18) :  4 node! 
	{
		switch (i - M_post1_q2)
		{
		case 0:// 1 lei shang 2j
			kind1_node(j, i, -(Num19_p1 - Node192_p1 + 1) - 2, -(Num19_p1 - Node192_p1 + 1) - 1, -(Num19_p1 - Node192_p1 + 1),
				0, 1, (Num18_p1 + Node173_p1 + 4) - 2, (Num18_p1 + Node173_p1 + 4) - 1, (Num18_p1 + Node173_p1 + 4),
				A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case 1://inner 4 node
			inner_4node(j, i, -(Num19_p1 - Node192_p1 + 1), -1, 0, (Num18_p1 + Node173_p1 + 4), 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case 2://inner 4 node
			inner_4node(j, i, -(Num19_p1 - Node192_p1 + 1), 0, 1, (Num18_p1 + Node173_p1 + 4) + 3, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Num18_p1 - 1) :// 1 lei xia 2j
			kind1_node(j, i, -(Num19_p1 - Node192_p1 + 1), -(Num19_p1 - Node192_p1 + 1) + 1, -(Num19_p1 - Node192_p1 + 1) + 2,
			-1, 0, (Num18_p1 + Node173_p1 + 4) + 3, (Num18_p1 + Node173_p1 + 4) + 3 + 1, (Num18_p1 + Node173_p1 + 4) + 3 + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		}
		j = j + 1;
	}
	cout << "之前无错5" << endl;
	//******************************************************************************************
	for (i = (M_post1_q2 + Num18_p1);
		i<(M_post1_q2 + Num18_p1 + Node_17_p1 + 2); i++)// circle num21(17) :  83 node!  part1 !
	{
		switch (i - M_post1_q2 - Num18_p1)
		{
		case 0:
			boundary_node(j, i, 1, -Num17_p1 - Num18_p1 - Num19_p1 - Num18_p1, Num17_p1 + Num16_p1 + Num15_p1 + Num14_p1, Coef, Coef_location);
			break;
		case (N_d0 - 3) ://2 lei shang 1j
			kind2_node(j, i, -Num17_p1 - Num18_p1 - Num19_p1 - Num18_p1, -2, -1, 0, 1, Num17_p1 + Num16_p1 + Num15_p1 + Num14_p1,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_17_p1 + 1) ://2 lei xia 1j
			kind2_node(j, i, -(Num17_p1 + Num18_p1 + Num19_p1 + Num18_p1), -1, 0, 1, 2,
			(Num17_p1 - Node_17_p1 + Num16_p1 + Num15_p1 + Num14_p1 + Node_13_p1), B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - M_post1_q2 - Num18_p1) < (N_d0 - 3))// 4dx!
			{
				inner_5node(j, i, 4, -Num17_p1 - Num18_p1 - Num19_p1 - Num18_p1, Num17_p1 + Num16_p1 + Num15_p1 + Num14_p1, Coef, Coef_location);
			}
			else if (((i - M_post1_q2 - Num18_p1) > (N_d0 - 3)) && ((i - M_post1_q2 - Num18_p1) < Node171_p1))// 2dx!
			{
				inner_5node(j, i, 2, -(Num19_p1 + Num18_p1 + N_d0 - 2), Num17_p1 - N_d0 + 2 + Num16_p1, Coef, Coef_location);
			}
			else if (((i - M_post1_q2 - Num18_p1) >= Node171_p1) && ((i - M_post1_q2 - Num18_p1) <= Node172_p1))// 2 lei right 2j !
			{
				kind2_node(j, i, -Num17_p1 - Num18_p1 - Num19_p1 - Num18_p1, -(Num19_p1 + Num18_p1 + N_d0 - 2), -1, 0, 1,
					Num17_p1 - (i - M_post1_q2 - Num18_p1) + (i - M_post1_q2 - Num18_p1 - Node171_p1) * 2,
					B5, B4, B3, B1_2, B3, B2, Coef, Coef_location);
			}
			else// 2dx!
			{
				inner_5node(j, i, 2, -(Num19_p1 + Num18_p1 + N_d0 - 2), Num17_p1 - Node_17_p1 + Num16_p1 + Node151_p1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post1_q2 + Num18_p1 + Node_17_p1 + 2);
		i<(M_post1_q2 + Num18_p1 + Node174_p1 + 3); i++)// circle num 21(17) :  83 node!  part2 !
	{
		switch (i - M_post1_q2 - Num18_p1)
		{
		case Node173_p1:// 2lei shang 1j
			kind2_node(j, i, -(Num17_p1 + Num18_p1 + Num19_p1 + Num18_p1), -2, -1, 0, 1,
				(Num17_p1 - Node_17_p1 + Num16_p1 + Num15_p1 + Num14_p1 + Node_13_p1), B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node173_p1 + 3) :// 2lei shang 2j
			kind2_node(j, i, -(Num19_p1 - Node191_p1 + Num18_p1 + Node173_p1), -2, -1, 0, 1,
			(Num17_p1 - Node173_p1 + Num16_p1 + Node151_p1), B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case (Node173_p1 + 4) :// 3lei zuo 2j
			kind3_node(j, i, -(Num18_p1 + Node173_p1 + 4), -1, 0, 1,
			(Num17_p1 - Node173_p1 + Num16_p1 + Node151_p1) - 1, (Num17_p1 - Node173_p1 + Num16_p1 + Node151_p1),
			C2_2, C3_2, C1_2, C3_2, C4, C4, Coef, Coef_location);
			break;
		case (Node173_p1 + 5) ://inner 5 node
			inner_5node(j, i, 1, -(Num18_p1 + Node173_p1 + 4), (Num17_p1 - Node173_p1 - 4 + Node161_p1), Coef, Coef_location);
			break;
			//************************************************
		case (Node174_p1 - 5) ://inner 5 node
			inner_5node(j, i, 1, -(1 + Node174_p1 - 4), (Num17_p1 - Node173_p1 - 4 + Node161_p1), Coef, Coef_location);
			break;
		case (Node174_p1 - 4) :// 3lei zuo 2j
			kind3_node(j, i, -(1 + Node174_p1 - 4), -1, 0, 1,
			(Num17_p1 - Node174_p1 + Num16_p1 + Node153_p1), (Num17_p1 - Node174_p1 + Num16_p1 + Node153_p1) + 1,
			C2_2, C3_2, C1_2, C3_2, C4, C4, Coef, Coef_location);
			break;
		case (Node174_p1 - 3) :// 2lei xia 2j
			kind2_node(j, i, -(Num19_p1 - Node193_p1 + Num18_p1 + Node174_p1), -1, 0, 1, 2,
			(Num17_p1 - Node174_p1 + Num16_p1 + Node153_p1), B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node174_p1) :// 2lei xia 1j
			kind2_node(j, i, -(Num17_p1 + Num18_p1 + Num19_p1 + Num18_p1), -1, 0, 1, 2,
			(Num17_p1 - Node174_p1 + Num16_p1 + Num15_p1 + Num14_p1 + Node133_p1), B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - M_post1_q2 - Num18_p1) < (Node173_p1))  // 2dx
			{
				inner_5node(j, i, 4, -(Num17_p1 + Num18_p1 + Num19_p1 + Num18_p1),
					(Num17_p1 - Node_17_p1 + Num16_p1 + Num15_p1 + Num14_p1 + Node_13_p1), Coef, Coef_location);
			}
			else if (((i - M_post1_q2 - Num18_p1) > (Node173_p1)) && ((i - M_post1_q2 - Num18_p1) < (Node173_p1 + 3)))  // 1dx
			{
				inner_5node(j, i, 2, -(Num19_p1 - Node191_p1 + Num18_p1 + Node173_p1),
					(Num17_p1 - Node173_p1 + Num16_p1 + Node151_p1), Coef, Coef_location);
			}
			else if (((i - M_post1_q2 - Num18_p1) > (Node173_p1 + 5)) && ((i - M_post1_q2 - Num18_p1) < (Node174_p1 - 5)))  // 1dx
			{
				inner_4node(j, i, -1, 0, 1, (Num17_p1 - Node173_p1 - 4 + Node161_p1), 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			}
			else if (((i - M_post1_q2 - Num18_p1) > (Node174_p1 - 3)) && ((i - M_post1_q2 - Num18_p1) < (Node174_p1)))  // 1dx
			{
				inner_5node(j, i, 2, -(Num19_p1 - Node193_p1 + Num18_p1 + Node174_p1),
					(Num17_p1 - Node174_p1 + Num16_p1 + Node153_p1), Coef, Coef_location);
			}
			else
			{
				inner_5node(j, i, 4, -(Num17_p1 + Num18_p1 + Num19_p1 + Num18_p1),
					(Num17_p1 - Node174_p1 + Num16_p1 + Num15_p1 + Num14_p1 + Node133_p1), Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post1_q2 + Num18_p1 + Node174_p1 + 3);
		i<(M_post1_q2 + Num18_p1 + Num17_p1); i++)// circle num 21(17) :  80 node!  part3 !
	{
		switch (i - M_post1_q2 - Num18_p1)
		{
		case (Num17_p1 - 1) :
			boundary_node(j, i, 2, -Num18_p1 - Num19_p1 - Num18_p1 - Num17_p1, Num16_p1 + Num15_p1 + Num14_p1 + Num13_p1, Coef, Coef_location);
			break;
		case (Node174_p1 + 3) ://2 lei shang 1j
			kind2_node(j, i, -(Num17_p1 + Num18_p1 + Num19_p1 + Num18_p1), -2, -1, 0, 1,
			(Num17_p1 - Node174_p1 + Num16_p1 + Num15_p1 + Num14_p1 + Node133_p1), B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_1717_p1 + 1) ://2 lei xia 1j
			kind2_node(j, i, -Num18_p1 - Num19_p1 - Num18_p1 - Num17_p1, -1, 0, 1, 2, Num16_p1 + Num15_p1 + Num14_p1 + Num13_p1,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - M_post1_q2 - Num18_p1) < Node175_p1)// 2dx!
			{
				inner_5node(j, i, 2, -(1 + Num18_p1 + Node_1717_p1), Num17_p1 - Node175_p1 + 4 + Num16_p1 + Node153_p1, Coef, Coef_location);
			}
			else if (((i - M_post1_q2 - Num18_p1) >= Node175_p1) && ((i - M_post1_q2 - Num18_p1) <= Node176_p1))// 2 lei left 2j !
			{
				kind2_node(j, i, -Num18_p1 - Num19_p1 - Num18_p1 - Num17_p1, -(1 + Num18_p1 + Node_1717_p1), -1, 0, 1,
					Num17_p1 - (i - M_post1_q2 - Num18_p1) + Node162_p1 + (i - M_post1_q2 - Num18_p1 - Node175_p1) * 2,
					B5, B4, B3, B1_2, B3, B2, Coef, Coef_location);
			}
			else if (((i - M_post1_q2 - Num18_p1) > Node176_p1) && ((i - M_post1_q2 - Num18_p1) <= Node_1717_p1))// 2dx
			{
				inner_5node(j, i, 2, -(1 + Num18_p1 + Node_1717_p1), Num17_p1 - Node_1717_p1 + Num16_p1 + Num15_p1 - 1, Coef, Coef_location);
			}
			else
			{
				inner_5node(j, i, 4, -Num18_p1 - Num19_p1 - Num18_p1 - Num17_p1, Num16_p1 + Num15_p1 + Num14_p1 + Num13_p1, Coef, Coef_location);// 4dx!
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (M_post1_q2 + Num18_p1 + Num17_p1);
		i<(M_post1_q2 + Num18_p1 + Num17_p1 + Node161_p1 + 1); i++)// num 22(16) :  31 node! part1
	{
		switch (i - M_post1_q2 - Num18_p1 - Num17_p1)
		{
		case 0:// 3 lei 2j shang
			kind3_node(j, i, -(Num17_p1 - Node171_p1) - 1, -(Num17_p1 - Node171_p1), 0, 1, (Num16_p1 + 5) - 2, (Num16_p1 + 5),
				C4, C3_2, C1_2, C2_2, C4, C3_2, Coef, Coef_location);
			break;
		case Node161_p1://node2: 3 lei 2j xia
			kind3_node(j, i, -(Num17_p1 - Node172_p1 + Node161_p1), -(Num17_p1 - Node172_p1 + Node161_p1) + 1, -1, 0, (Num16_p1 + 5),
				(Num16_p1 + 5) + 2, C3_2, C4, C2_2, C1_2, C3_2, C4, Coef, Coef_location);
			break;
		default:
			if ((i - M_post1_q2 - Num18_p1 - Num17_p1) % 2 == 0)// 1dx
			{
				inner_5node(j, i, 1,
					-((i - M_post1_q2 - Num18_p1 - Num17_p1) + Num17_p1 - Node171_p1 - (i - M_post1_q2 - Num18_p1 - Num17_p1) / 2),
					(Num16_p1 + 5), Coef, Coef_location);
			}
			else//1 lei  2j left
			{
				kind1_node(j, i, -((i - M_post1_q2 - Num18_p1 - Num17_p1) + Num17_p1 - Node171_p1 - (i - M_post1_q2 - Num18_p1 - Num17_p1) / 2) - (Num19_p1 + Num18_p1 + N_d0 - 2),
					-((i - M_post1_q2 - Num18_p1 - Num17_p1) + Num17_p1 - Node171_p1 - (i - M_post1_q2 - Num18_p1 - Num17_p1) / 2) + 1 - (Num19_p1 + Num18_p1 + N_d0 - 2),
					-((i - M_post1_q2 - Num18_p1 - Num17_p1) + Num17_p1 - Node171_p1 - (i - M_post1_q2 - Num18_p1 - Num17_p1) / 2),
					-((i - M_post1_q2 - Num18_p1 - Num17_p1) + Num17_p1 - Node171_p1 - (i - M_post1_q2 - Num18_p1 - Num17_p1) / 2) + 1,
					-1, 0, 1, (Num16_p1 + 5),
					A2_2, A2_2, A1_2, A1_2, A3_2, A5_2, A3_2, A4_2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post1_q2 + Num18_p1 + Num17_p1 + Node161_p1 + 1);
		i<(M_post1_q2 + Num18_p1 + Num17_p1 + Node162_p1); i++)// num 22(16) :  31 node! part2
	{
		switch (i - M_post1_q2 - Num18_p1 - Num17_p1)
		{
		case (Node161_p1 + 1) :// 3 lei 2j shang !!! 1
			kind3_node(j, i, -(Num17_p1 - Node173_p1 - 4 + Node161_p1) - 2, -(Num17_p1 - Node173_p1 - 4 + Node161_p1), 0, 1,
			(Num16_p1 - Node161_p1 - 1 + Node151_p1 + 4) - 1, (Num16_p1 - Node161_p1 - 1 + Node151_p1 + 4),
			C4, C3_2, C1_2, C2_2, C4, C3_2, Coef, Coef_location);
			break;
		case (Node161_p1 + 2) :// 1 lei 2j zuo !!! 2
			kind1_node(j, i, -(Num17_p1 - Node173_p1 - 4 + Node161_p1), -1, 0, 1,
			(Num16_p1 - Node161_p1 - 1 + Node151_p1 + 4) - 1, (Num16_p1 - Node161_p1 - 1 + Node151_p1 + 4),
			(Num16_p1 - Node161_p1 - 1 + Node151_p1 + 4) - 1 + (Num15_p1 - Node151_p1 + Num14_p1 + Node132_p1),
			(Num16_p1 - Node161_p1 - 1 + Node151_p1 + 4) + (Num15_p1 - Node151_p1 + Num14_p1 + Node132_p1),
			A4_2, A3_2, A5_2, A3_2, A1_2, A1_2, A2_2, A2_2, Coef, Coef_location);
			break;
		case (Node161_p1 + 3) :// inner 5 node !!! 3
			inner_5node(j, i, 1, -(Num17_p1 - Node173_p1 - 4 + Node161_p1), (Num16_p1 - Node161_p1 - 1 + Node151_p1 + 4) - 1, Coef, Coef_location);
			break;
		case (Node162_p1 - 2) :// 1 lei 2j zuo !!! 4
			kind1_node(j, i, -(Num17_p1 - Node173_p1 - 4 + Node161_p1), -1, 0, 1,
			(Num16_p1 - Node161_p1 - 1 + Node151_p1 + 4) - 2, (Num16_p1 - Node161_p1 - 1 + Node151_p1 + 4) - 1,
			(Num16_p1 - Node161_p1 - 1 + Node151_p1 + 4) - 2 + (Num15_p1 - Node151_p1 + Num14_p1 + Node132_p1),
			(Num16_p1 - Node161_p1 - 1 + Node151_p1 + 4) - 1 + (Num15_p1 - Node151_p1 + Num14_p1 + Node132_p1),
			A4_2, A3_2, A5_2, A3_2, A1_2, A1_2, A2_2, A2_2, Coef, Coef_location);
			break;
		case (Node162_p1 - 1) :// 3 lei 2j xia !!! 5
			kind3_node(j, i, -(Num17_p1 - Node173_p1 - 4 + Node161_p1), -(Num17_p1 - Node173_p1 - 4 + Node161_p1) + 2, -1, 0,
			(Num16_p1 - Node161_p1 - 1 + Node151_p1 + 4) - 2, (Num16_p1 - Node161_p1 - 1 + Node151_p1 + 4) - 2 + 1,
			C3_2, C4, C2_2, C1_2, C3_2, C4, Coef, Coef_location);
			break;
		}
		j = j + 1;
	}
	for (i = (M_post1_q2 + Num18_p1 + Num17_p1 + Node162_p1);
		i<(M_post1_q2 + Num18_p1 + Num17_p1 + Num16_p1); i++)// num 22(16) :  31 node! part3
	{
		switch (i - M_post1_q2 - Num18_p1 - Num17_p1 - Node162_p1)
		{
		case 0:// 3 lei 2j shang
			kind3_node(j, i, -(Num17_p1 - Node175_p1 + Node162_p1) - 1, -(Num17_p1 - Node175_p1 + Node162_p1), 0, 1,
				(Num16_p1 - Node162_p1 + Node153_p1 + 5) - 2, (Num16_p1 - Node162_p1 + Node153_p1 + 5),
				C4, C3_2, C1_2, C2_2, C4, C3_2, Coef, Coef_location);
			break;
		case Node161_p1:// 3 lei 2j xia
			kind3_node(j, i, -(Num17_p1 - Node176_p1 + Num16_p1 - 1), -(Num17_p1 - Node176_p1 + Num16_p1 - 1) + 1, -1, 0,
				(Num16_p1 - Node162_p1 + Node153_p1 + 5), (Num16_p1 - Node162_p1 + Node153_p1 + 5) + 2,
				C3_2, C4, C2_2, C1_2, C3_2, C4, Coef, Coef_location);
			break;
		default:
			if ((i - M_post1_q2 - Num18_p1 - Num17_p1 - Node162_p1) % 2 == 0)// 1dx
			{
				inner_5node(j, i, 1,
					-((i - M_post1_q2 - Num18_p1 - Num17_p1) + Num17_p1 - Node175_p1 - (i - M_post1_q2 - Num18_p1 - Num17_p1 - Node162_p1) / 2),
					(Num16_p1 - Node162_p1 + Node153_p1 + 5), Coef, Coef_location);
			}
			else//1 lei 2j left
			{
				kind1_node(j, i, -((i - M_post1_q2 - Num18_p1 - Num17_p1) + Num17_p1 - Node175_p1 - (i - M_post1_q2 - Num18_p1 - Num17_p1 - Node162_p1) / 2) - (1 + Num18_p1 + Node_1717_p1),
					-((i - M_post1_q2 - Num18_p1 - Num17_p1) + Num17_p1 - Node175_p1 - (i - M_post1_q2 - Num18_p1 - Num17_p1 - Node162_p1) / 2) + 1 - (1 + Num18_p1 + Node_1717_p1),
					-((i - M_post1_q2 - Num18_p1 - Num17_p1) + Num17_p1 - Node175_p1 - (i - M_post1_q2 - Num18_p1 - Num17_p1 - Node162_p1) / 2),
					-((i - M_post1_q2 - Num18_p1 - Num17_p1) + Num17_p1 - Node175_p1 - (i - M_post1_q2 - Num18_p1 - Num17_p1 - Node162_p1) / 2) + 1,
					-1, 0, 1, (Num16_p1 - Node162_p1 + Node153_p1 + 5),
					A2_2, A2_2, A1_2, A1_2, A3_2, A5_2, A3_2, A4_2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (M_post1_q2 + Num18_p1 + Num17_p1 + Num16_p1);
		i<(M_post1_q2 + Num18_p1 + Num17_p1 + Num16_p1 + Node151_p1 + 1); i++)// num 23(15) :  55 node! part1
	{
		switch (i - M_post1_q2 - Num18_p1 - Num17_p1 - Num16_p1)
		{
		case 0:// 1 lei 1j shang
			kind1_node(j, i, -(Num17_p1 - N_d0 + 2 + Num16_p1) - 2, -(Num17_p1 - N_d0 + 2 + Num16_p1) - 1, -(Num17_p1 - N_d0 + 2 + Num16_p1),
				0, 1, (Num15_p1 + Num14_p1 + N_d0 - 2) - 2, (Num15_p1 + Num14_p1 + N_d0 - 2) - 1, (Num15_p1 + Num14_p1 + N_d0 - 2),
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case 3:// 2 lei 2j shang
			kind2_node(j, i, -(Num17_p1 - N_d0 + 2 + Num16_p1), -2, -1, 0, 1, (Num15_p1 + Num14_p1 + N_d0 - 2),
				B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case 4:// 3 lei 2j zuo
			kind3_node(j, i, -(Num17_p1 - N_d0 + 2 + Num16_p1) - 1, -(Num17_p1 - N_d0 + 2 + Num16_p1), -1, 0, 1, Num15_p1 - 4,
				C4, C4, C3_2, C1_2, C3_2, C2_2, Coef, Coef_location);
			break;
		case (Node151_p1 - 4) :// 3 lei 2j zuo
			kind3_node(j, i, -(Num17_p1 - Node_17_p1 + Num16_p1 + Node151_p1), -(Num17_p1 - Node_17_p1 + Num16_p1 + Node151_p1) + 1,
			-1, 0, 1, Num15_p1 - 4,
			C4, C4, C3_2, C1_2, C3_2, C2_2, Coef, Coef_location);
			break;
		case (Node151_p1 - 3) :// 2 lei 2j xia
			kind2_node(j, i, -(Num17_p1 - Node_17_p1 + Num16_p1 + Node151_p1), -1, 0, 1, 2, (Num15_p1 - Node151_p1 + Num14_p1 + Node_13_p1),
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case Node151_p1:// 1 lei 1j xia
			kind1_node(j, i, -(Num17_p1 - Node_17_p1 + Num16_p1 + Node151_p1), -(Num17_p1 - Node_17_p1 + Num16_p1 + Node151_p1) + 1,
				-(Num17_p1 - Node_17_p1 + Num16_p1 + Node151_p1) + 2, -1, 0, (Num15_p1 - Node151_p1 + Num14_p1 + Node_13_p1),
				(Num15_p1 - Node151_p1 + Num14_p1 + Node_13_p1) + 1, (Num15_p1 - Node151_p1 + Num14_p1 + Node_13_p1) + 2,
				A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - M_post1_q2 - Num18_p1 - Num17_p1 - Num16_p1) < 3)// 2dx
			{
				inner_5node(j, i, 2, -(Num17_p1 - N_d0 + 2 + Num16_p1), (Num15_p1 + Num14_p1 + N_d0 - 2), Coef, Coef_location);
			}
			else if (((i - M_post1_q2 - Num18_p1 - Num17_p1 - Num16_p1) > 4) && ((i - M_post1_q2 - Num18_p1 - Num17_p1 - Num16_p1) < (Node151_p1 - 4)))// 1dx
			{
				inner_5node(j, i, 1, -(Num16_p1 + 5), Num15_p1 - 4, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, i, 2, -(Num17_p1 - Node_17_p1 + Num16_p1 + Node151_p1), (Num15_p1 - Node151_p1 + Num14_p1 + Node_13_p1), Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post1_q2 + Num18_p1 + Num17_p1 + Num16_p1 + Node151_p1 + 1);
		i<(M_post1_q2 + Num18_p1 + Num17_p1 + Num16_p1 + Node153_p1); i++)// num 23(15) :  62 node! part2
	{
		switch (i - M_post1_q2 - Num18_p1 - Num17_p1 - Num16_p1)
		{
		case (Node151_p1 + 1) :// 1 lei 1j shang
			kind1_node(j, i, -(Num17_p1 - Node173_p1 + Num16_p1 + Node151_p1) - 2, -(Num17_p1 - Node173_p1 + Num16_p1 + Node151_p1) - 1,
			-(Num17_p1 - Node173_p1 + Num16_p1 + Node151_p1), 0, 1, (Num15_p1 - Node151_p1 + Num14_p1 + Node132_p1) - 2,
			(Num15_p1 - Node151_p1 + Num14_p1 + Node132_p1) - 1, (Num15_p1 - Node151_p1 + Num14_p1 + Node132_p1),
			A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case (Node152_p1) :// 2 lei 2j zuo
			kind2_node(j, i, -(Num16_p1 - Node161_p1 - 1 + Node152_p1), -1, 0, 1, (Num15_p1 - Node151_p1 + Num14_p1 + Node132_p1),
			(Num15_p1 - Node151_p1 + Num14_p1 + Node132_p1) + (Num13_p1 - Node132_p1 - 1 + Num12_p1 + Node112_p1),
			B2, B3, B1_2, B3, B4, B5, Coef, Coef_location);
			break;
		case (Node152_p1 + 1) :// 2 lei 2j zuo
			kind2_node(j, i, -(Num16_p1 - Node161_p1 - 1 + Node152_p1) + 1, -1, 0, 1, (Num15_p1 - Node151_p1 + Num14_p1 + Node132_p1),
			(Num15_p1 - Node151_p1 + Num14_p1 + Node132_p1) + (Num13_p1 - Node132_p1 - 1 + Num12_p1 + Node112_p1),
			B2, B3, B1_2, B3, B4, B5, Coef, Coef_location);
			break;
		case (Node152_p1 + 2) :// 2 lei 2j zuo
			kind2_node(j, i, -(Num16_p1 - Node161_p1 - 1 + Node152_p1) + 2, -1, 0, 1, (Num15_p1 - Node151_p1 + Num14_p1 + Node132_p1),
			(Num15_p1 - Node151_p1 + Num14_p1 + Node132_p1) + (Num13_p1 - Node132_p1 - 1 + Num12_p1 + Node112_p1),
			B2, B3, B1_2, B3, B4, B5, Coef, Coef_location);
			break;
		case (Node153_p1 - 1) :// 1 lei 1j xia
			kind1_node(j, i, -(Num17_p1 - Node174_p1 + Num16_p1 + Node153_p1), -(Num17_p1 - Node174_p1 + Num16_p1 + Node153_p1) + 1,
			-(Num17_p1 - Node174_p1 + Num16_p1 + Node153_p1) + 2, -1, 0, (Num15_p1 - Node151_p1 + Num14_p1 + Node132_p1),
			(Num15_p1 - Node151_p1 + Num14_p1 + Node132_p1) + 1, (Num15_p1 - Node151_p1 + Num14_p1 + Node132_p1) + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - M_post1_q2 - Num18_p1 - Num17_p1 - Num16_p1) < Node152_p1)// 1dx
			{
				inner_5node(j, i, 2, -(Num17_p1 - Node173_p1 + Num16_p1 + Node151_p1), (Num15_p1 - Node151_p1 + Num14_p1 + Node132_p1), Coef, Coef_location);
			}
			else if ((i - M_post1_q2 - Num18_p1 - Num17_p1 - Num16_p1) > (Node152_p1 + 2))// 1dx
			{
				inner_5node(j, i, 2, -(Num17_p1 - Node174_p1 + Num16_p1 + Node153_p1), (Num15_p1 - Node151_p1 + Num14_p1 + Node132_p1), Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post1_q2 + Num18_p1 + Num17_p1 + Num16_p1 + Node153_p1);
		i<(M_post1_q2 + Num18_p1 + Num17_p1 + Num16_p1 + Num15_p1); i++)// num 23(15) :  55 node! part3
	{
		switch (i - M_post1_q2 - Num18_p1 - Num17_p1 - Num16_p1)
		{
		case Node153_p1:// 1 lei 1j shang
			kind1_node(j, i, -(Num17_p1 - Node174_p1 - 4 + Num16_p1 + Node153_p1) - 2, -(Num17_p1 - Node174_p1 - 4 + Num16_p1 + Node153_p1) - 1,
				-(Num17_p1 - Node174_p1 - 4 + Num16_p1 + Node153_p1), 0, 1, (Num15_p1 - Node153_p1 + Num14_p1 + Node133_p1 + 4) - 2,
				(Num15_p1 - Node153_p1 + Num14_p1 + Node133_p1 + 4) - 1, (Num15_p1 - Node153_p1 + Num14_p1 + Node133_p1 + 4),
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case (Node153_p1 + 3) :// 2 lei 2j shang
			kind2_node(j, i, -(Num17_p1 - Node174_p1 - 4 + Num16_p1 + Node153_p1), -2, -1, 0, 1,
			(Num15_p1 - Node153_p1 + Num14_p1 + Node133_p1 + 4), B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case (Node153_p1 + 4) :// 3 lei 2j shang
			kind3_node(j, i, -(Num17_p1 - Node174_p1 - 4 + Num16_p1 + Node153_p1) - 1, -(Num17_p1 - Node174_p1 - 4 + Num16_p1 + Node153_p1),
			-1, 0, 1, (5 + Num14_p1 - 1), C4, C4, C3_2, C1_2, C3_2, C2_2, Coef, Coef_location);
			break;
		case (Num15_p1 - 5) :// 3 lei 2j xia
			kind3_node(j, i, -(Num17_p1 - Node_1717_p1 + Num16_p1 + Num15_p1 - 1), -(Num17_p1 - Node_1717_p1 + Num16_p1 + Num15_p1 - 1) + 1,
			-1, 0, 1, (5 + Num14_p1 - 1),
			C4, C4, C3_2, C1_2, C3_2, C2_2, Coef, Coef_location);
			break;
		case (Num15_p1 - 4) :// 2 lei 2j xia
			kind2_node(j, i, -(Num17_p1 - Node_1717_p1 + Num16_p1 + Num15_p1 - 1), -1, 0, 1, 2, (1 + Num14_p1 + Node_1313_p1),
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num15_p1 - 1) :// 1 lei 1j xia
			kind1_node(j, i, -(Num17_p1 - Node_1717_p1 + Num16_p1 + Num15_p1 - 1), -(Num17_p1 - Node_1717_p1 + Num16_p1 + Num15_p1 - 1) + 1,
			-(Num17_p1 - Node_1717_p1 + Num16_p1 + Num15_p1 - 1) + 2, -1, 0, (1 + Num14_p1 + Node_1313_p1),
			(1 + Num14_p1 + Node_1313_p1) + 1, (1 + Num14_p1 + Node_1313_p1) + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - M_post1_q2 - Num18_p1 - Num17_p1 - Num16_p1) < (Node153_p1 + 3))// 2dx
			{
				inner_5node(j, i, 2, -(Num17_p1 - Node174_p1 - 4 + Num16_p1 + Node153_p1), (Num15_p1 - Node153_p1 + Num14_p1 + Node133_p1 + 4), Coef, Coef_location);
			}
			else if (((i - M_post1_q2 - Num18_p1 - Num17_p1 - Num16_p1) > (Node153_p1 + 4)) && ((i - M_post1_q2 - Num18_p1 - Num17_p1 - Num16_p1) < (Num15_p1 - 4)))// 1dx
			{
				inner_5node(j, i, 1, -(1 + Num15_p1 - 6), (5 + Num14_p1 - 1), Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, i, 2, -(Num17_p1 - Node_1717_p1 + Num16_p1 + Num15_p1 - 1), (1 + Num14_p1 + Node_1313_p1), Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!! num 24 !!!
	for (i = (M_post1_q2 + Num18_p1 + Num17_p1 + Num16_p1 + Num15_p1);
		i<(M_post1_q2 + Num18_p1 + Num17_p1 + Num16_p1 + Num15_p1 + Num14_p1 / 2); i++)// num 24(14) :  45 node! part1
	{
		switch (i - M_post1_q2 - Num18_p1 - Num17_p1 - Num16_p1 - Num15_p1)
		{
		case 0:// 1 lei 2j shang
			kind1_node(j, i, -(Num15_p1 - 4) - 2, -(Num15_p1 - 4) - 1, -(Num15_p1 - 4), 0, 1,
				Num14_p1 + N_d0 + 2 - 2, Num14_p1 + N_d0 + 2 - 1, Num14_p1 + N_d0 + 2,
				A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case (Num4_p1 / 2 - 1) :// 1 lei 2j xia
			kind1_node(j, i, -(Num15_p1 - 4), -(Num15_p1 - 4) + 1, -(Num15_p1 - 4) + 2, -1, 0,
			Num14_p1 + N_d0 + 2 - N_c1, Num14_p1 + N_d0 + 2 - N_c1 + 1, Num14_p1 + N_d0 + 2 - N_c1 + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - M_post1_q2 - Num18_p1 - Num17_p1 - Num16_p1 - Num15_p1) <= Node141_p1)// 1dx
			{
				inner_5node(j, i, 1, -(Num15_p1 - 4), Num14_p1 + N_d0 + 2, Coef, Coef_location);
			}
			else if ((i - M_post1_q2 - Num18_p1 - Num17_p1 - Num16_p1 - Num15_p1) >= Node142_p1) // 1dx
			{
				inner_5node(j, i, 1, -(Num15_p1 - 4), Num14_p1 + N_d0 + 2 - N_c1, Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_4node(j, i, -(Num15_p1 - 4), -1, 0, 1, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post1_q2 + Num18_p1 + Num17_p1 + Num16_p1 + Num15_p1 + Num14_p1 / 2);
		i<(M_post1_q2 + Num18_p1 + Num17_p1 + Num16_p1 + Num15_p1 + Num14_p1); i++)// num 4 :  30 node! part2
	{
		switch (i - M_post1_q2 - Num18_p1 - Num17_p1 - Num16_p1 - Num15_p1)
		{
		case (Num14_p1 / 2) :// 1 lei 2j shang
			kind1_node(j, i, -(5 + Num14_p1 - 1) - 2, -(5 + Num14_p1 - 1) - 1, -(5 + Num14_p1 - 1), 0, 1,
			(Node_1313_p1 - 4 + 1 + N_c1) - 2, (Node_1313_p1 - 4 + 1 + N_c1) - 1, (Node_1313_p1 - 4 + 1 + N_c1),
			A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case (Num14_p1 - 1) :// 1 lei 2j xia
			kind1_node(j, i, -(5 + Num14_p1 - 1), -(5 + Num14_p1 - 1) + 1, -(5 + Num14_p1 - 1) + 2, -1, 0,
			(Node_1313_p1 - 4 + 1), (Node_1313_p1 - 4 + 1) + 1, (Node_1313_p1 - 4 + 1) + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - M_post1_q2 - Num18_p1 - Num17_p1 - Num16_p1 - Num15_p1) <= Node143_p1)// 1dx
			{
				inner_5node(j, i, 1, -(5 + Num14_p1 - 1), (Node_1313_p1 - 4 + 1 + N_c1), Coef, Coef_location);
			}
			else if ((i - M_post1_q2 - Num18_p1 - Num17_p1 - Num16_p1 - Num15_p1) >= Node144_p1) // 1dx
			{
				inner_5node(j, i, 1, -(5 + Num14_p1 - 1), (Node_1313_p1 - 4 + 1), Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_4node(j, i, -(5 + Num14_p1 - 1), -1, 0, 1, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!! num 25 !!!
	for (i = (M_post1_q2 + Num18_p1 + Num17_p1 + Num16_p1 + Num15_p1 + Num14_p1);
		i<(M_post1_q2 + Num18_p1 + Num17_p1 + Num16_p1 + Num15_p1 + Num14_p1 + Node_13_p1 + 2); i++)// circle num 25(13) :  85 node!  part1 !
	{
		switch (i - M_post1_q2 - Num18_p1 - Num17_p1 - Num16_p1 - Num15_p1 - Num14_p1)
		{
		case 0:
			boundary_node(j, i, 1, -Num17_p1 - Num16_p1 - Num15_p1 - Num14_p1, Num13_p1 + Num12_p1 + Num11_p1 + Num10_p1, Coef, Coef_location);
			break;
		case (N_d0 - 3) ://2 lei shang 1j
			kind2_node(j, i, -Num17_p1 - Num16_p1 - Num15_p1 - Num14_p1, -2, -1, 0, 1, Num13_p1 + Num12_p1 + Num11_p1 + Num10_p1,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (N_d0 + 1) :// 2 lei 2j shang
			kind2_node(j, i, -(Num15_p1 + Num14_p1 + N_d0 - 2), -2, -1, 0, 1, Num13_p1 - N_d0 + 2 + Num12_p1,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node131_p1:// inner 3 node
			inner_3node(j, i, -(Num14_p1 + N_d0 + 2), -1, 0, 1, 1, k*k*dx*dx - 4, Coef, Coef_location);
			break;
		case (Node131_p1 + 1) :// inner 3 node
			inner_3node(j, i, -(Num14_p1 + N_d0 + 2) + N_c1, 0, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node_13_p1 - 3) ://2 lei xia 2j
			kind2_node(j, i, -(Num15_p1 - Node151_p1 + Num14_p1 + Node_13_p1), -1, 0, 1, 2, Num13_p1 - Node_13_p1 + Num12_p1 + Node112_p1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_13_p1 + 1) ://2 lei xia 1j
			kind2_node(j, i, -(Num17_p1 - Node_17_p1 + Num16_p1 + Num15_p1 + Num14_p1 + Node_13_p1), -1, 0, 1, 2,
			Num13_p1 - Node_13_p1 + Num12_p1 + Num11_p1 + Num10_p1 + Node_9_p1, B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - M_post1_q2 - Num18_p1 - Num17_p1 - Num16_p1 - Num15_p1 - Num14_p1) < (N_d0 - 3))// 4dx!
			{
				inner_5node(j, i, 4, -Num17_p1 - Num16_p1 - Num15_p1 - Num14_p1, Num13_p1 + Num12_p1 + Num11_p1 + Num10_p1, Coef, Coef_location);
			}
			else if (((i - M_post1_q2 - Num18_p1 - Num17_p1 - Num16_p1 - Num15_p1 - Num14_p1) > (N_d0 - 3)) && ((i - M_post1_q2 - Num18_p1 - Num17_p1 - Num16_p1 - Num15_p1 - Num14_p1) < (N_d0 + 1)))// 2dx!
			{
				inner_5node(j, i, 2, -(Num15_p1 + Num14_p1 + N_d0 - 2), Num13_p1 - N_d0 + 2 + Num12_p1, Coef, Coef_location);
			}
			else if (((i - M_post1_q2 - Num18_p1 - Num17_p1 - Num16_p1 - Num15_p1 - Num14_p1) > (N_d0 + 1)) && ((i - M_post1_q2 - Num18_p1 - Num17_p1 - Num16_p1 - Num15_p1 - Num14_p1) < Node131_p1))// 1dx
			{
				inner_5node(j, i, 1, -(Num14_p1 + N_d0 + 2), (Num13_p1 - N_d0 - 2), Coef, Coef_location);
			}
			else if (((i - M_post1_q2 - Num18_p1 - Num17_p1 - Num16_p1 - Num15_p1 - Num14_p1) > (Node131_p1 + 1)) && ((i - M_post1_q2 - Num18_p1 - Num17_p1 - Num16_p1 - Num15_p1 - Num14_p1) < (Node_13_p1 - 3)))// 1dx
			{
				inner_5node(j, i, 1, -(Num14_p1 + N_d0 + 2) + N_c1, (Num13_p1 - N_d0 - 2) - 2, Coef, Coef_location);
			}
			else// 2dx!
			{
				inner_5node(j, i, 2, -(Num15_p1 - Node151_p1 + Num14_p1 + Node_13_p1), Num13_p1 - Node_13_p1 + Num12_p1 + Node112_p1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post1_q2 + Num18_p1 + Num17_p1 + Num16_p1 + Num15_p1 + Num14_p1 + Node_13_p1 + 2);
		i<(M_post1_q2 + Num18_p1 + Num17_p1 + Num16_p1 + Num15_p1 + Num14_p1 + Node133_p1 + 3); i++)// circle num 25(13) :  85 node!  part2 !
	{
		switch (i - M_post1_q2 - Num18_p1 - Num17_p1 - Num16_p1 - Num15_p1 - Num14_p1)
		{
		case (Node132_p1) ://2 lei shang
			kind2_node(j, i, -(Num17_p1 - Node_17_p1 + Num16_p1 + Num15_p1 + Num14_p1 + Node_13_p1), -2, -1, 0, 1,
			(Num13_p1 - Node_13_p1 + Num12_p1 + Num11_p1 + Num10_p1 + Node_9_p1), B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node132_p1 + 1) ://3 lei shang
			kind3_node(j, i, -(Num15_p1 - Node151_p1 + Num14_p1 + Node132_p1), -1, 0, 1,
			(Num13_p1 - Node_13_p1 + Num12_p1 + Num11_p1 + Num10_p1 + Node_9_p1) - 1, (Num13_p1 - Node_13_p1 + Num12_p1 + Num11_p1 + Num10_p1 + Node_9_p1),
			C2_1, C3_1, C1_1, C3_1, C4, C4, Coef, Coef_location);
			break;
		case (Node133_p1 - 1) ://3 lei xia
			kind3_node(j, i, -(Num15_p1 - Node151_p1 + Num14_p1 + Node132_p1), -1, 0, 1,
			(Num13_p1 - Node133_p1 + Num12_p1 + Num11_p1 + Num10_p1 + Node93_p1 + 1),
			(Num13_p1 - Node133_p1 + Num12_p1 + Num11_p1 + Num10_p1 + Node93_p1 + 1) + 1, C2_1, C3_1, C1_1, C3_1, C4, C4, Coef, Coef_location);
			break;
		case (Node133_p1) ://2 lei xia
			kind2_node(j, i, -(Num17_p1 - Node174_p1 + Num16_p1 + Num15_p1 + Num14_p1 + Node133_p1), -1, 0, 1, 2,
			(Num13_p1 - Node133_p1 + Num12_p1 + Num11_p1 + Num10_p1 + Node93_p1 + 1), B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - M_post1_q2 - Num18_p1 - Num17_p1 - Num16_p1 - Num15_p1 - Num14_p1) < (Node132_p1))// 4dx
			{
				inner_5node(j, i, 4, -(Num17_p1 - Node_17_p1 + Num16_p1 + Num15_p1 + Num14_p1 + Node_13_p1),
					Num13_p1 - Node_13_p1 + Num12_p1 + Num11_p1 + Num10_p1 + Node_9_p1, Coef, Coef_location);
			}
			else if ((i - M_post1_q2 - Num18_p1 - Num17_p1 - Num16_p1 - Num15_p1 - Num14_p1) > (Node133_p1))// 4dx!
			{
				inner_5node(j, i, 4, -(Num17_p1 - Node174_p1 + Num16_p1 + Num15_p1 + Num14_p1 + Node133_p1),
					(Num13_p1 - Node133_p1 + Num12_p1 + Num11_p1 + Num10_p1 + Node93_p1 + 1), Coef, Coef_location);
			}
			else if (((i - M_post1_q2 - Num18_p1 - Num17_p1 - Num16_p1 - Num15_p1 - Num14_p1) > (Node132_p1 + 1)) && ((i - M_post1_q2 - Num18_p1 - Num17_p1 - Num16_p1 - Num15_p1 - Num14_p1) < (Node133_p1 - 1)))
			{
				inner_5node(j, i, 2, -(Num15_p1 - Node151_p1 + Num14_p1 + Node132_p1),
					Num13_p1 - Node132_p1 - 2 + Num12_p1 + Node112_p1 + 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post1_q2 + Num18_p1 + Num17_p1 + Num16_p1 + Num15_p1 + Num14_p1 + Node133_p1 + 3);
		i<(M_post1_q2 + Num18_p1 + Num17_p1 + Num16_p1 + Num15_p1 + Num14_p1 + Num13_p1); i++)// circle num 25(13) :  95 node!  part3 !
	{
		switch (i - M_post1_q2 - Num18_p1 - Num17_p1 - Num16_p1 - Num15_p1 - Num14_p1)
		{
		case (Num13_p1 - 1) :
			boundary_node(j, i, 2, -Num16_p1 - Num15_p1 - Num14_p1 - Num13_p1, Num12_p1 + Num11_p1 + Num10_p1 + Num9_p1, Coef, Coef_location);
			break;
		case (Node133_p1 + 3) ://2 lei shang 1j
			kind2_node(j, i, -(Num17_p1 - Node174_p1 + Num16_p1 + Num15_p1 + Num14_p1 + Node133_p1), -2, -1, 0, 1,
			(Num13_p1 - Node133_p1 + Num12_p1 + Num11_p1 + Num10_p1 + Node93_p1 + 1), B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node133_p1 + 7) :// 2 lei 2j shang
			kind2_node(j, i, -(Num15_p1 - Node153_p1 + Num14_p1 + Node133_p1 + 4), -2, -1, 0, 1, Num13_p1 - Node133_p1 - 4 + Num12_p1 + Node113_p1,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node134_p1:// inner 3 node
			inner_3node(j, i, -(1 + Node_1313_p1 - 4) - N_c1, -1, 0, 1, 1, k*k*dx*dx - 4, Coef, Coef_location);
			break;
		case (Node134_p1 + 1) :// inner 3 node
			inner_3node(j, i, -(1 + Node_1313_p1 - 4), 0, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node_1313_p1 - 3) ://2 lei xia 2j
			kind2_node(j, i, -(Num14_p1 + Node_1313_p1 + 1), -1, 0, 1, 2, Num13_p1 - Node_1313_p1 + Num12_p1 + Num11_p1 - 1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_1313_p1 + 1) ://2 lei xia 1j
			kind2_node(j, i, -Num16_p1 - Num15_p1 - Num14_p1 - Num13_p1, -1, 0, 1, 2, Num12_p1 + Num11_p1 + Num10_p1 + Num9_p1,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - M_post1_q2 - Num18_p1 - Num17_p1 - Num16_p1 - Num15_p1 - Num14_p1) < (Node133_p1 + 7))// 2dx!
			{
				inner_5node(j, i, 2, -(Num15_p1 - Node153_p1 + Num14_p1 + Node133_p1 + 4), Num13_p1 - Node133_p1 - 4 + Num12_p1 + Node113_p1, Coef, Coef_location);
			}
			else if (((i - M_post1_q2 - Num18_p1 - Num17_p1 - Num16_p1 - Num15_p1 - Num14_p1) > (Node133_p1 + 7)) && ((i - M_post1_q2 - Num18_p1 - Num17_p1 - Num16_p1 - Num15_p1 - Num14_p1) < Node134_p1))// 1dx
			{
				inner_5node(j, i, 1, -(1 + Node_1313_p1 - 4) - N_c1, (Num13_p1 - Node_1313_p1 + 4 + Num12_p1 - 1) + 2, Coef, Coef_location);
			}
			else if (((i - M_post1_q2 - Num18_p1 - Num17_p1 - Num16_p1 - Num15_p1 - Num14_p1) > (Node134_p1 + 1)) && ((i - M_post1_q2 - Num18_p1 - Num17_p1 - Num16_p1 - Num15_p1 - Num14_p1) < (Node_1313_p1 - 3)))// 1dx
			{
				inner_5node(j, i, 1, -(1 + Node_1313_p1 - 4), (Num13_p1 - Node_1313_p1 + 4 + Num12_p1 - 1), Coef, Coef_location);
			}
			else if (((i - M_post1_q2 - Num18_p1 - Num17_p1 - Num16_p1 - Num15_p1 - Num14_p1) > (Node_1313_p1 - 3)) && ((i - M_post1_q2 - Num18_p1 - Num17_p1 - Num16_p1 - Num15_p1 - Num14_p1) < (Node_1313_p1 + 1)))// 2dx
			{
				inner_5node(j, i, 2, -(Num14_p1 + Node_1313_p1 + 1), Num13_p1 - Node_1313_p1 + Num12_p1 + Num11_p1 - 1, Coef, Coef_location);
			}
			else// 4dx!
			{
				inner_5node(j, i, 4, -Num16_p1 - Num15_p1 - Num14_p1 - Num13_p1, Num12_p1 + Num11_p1 + Num10_p1 + Num9_p1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!! num 26 !!!
	for (i = (M_post1_q3 - Num9_p1 - Num10_p1 - Num11_p1 - Num12_p1);
		i<(M_post1_q3 - Num9_p1 - Num10_p1 - Num11_p1); i++)// num 26(12) :  29 node!  part1
	{
		switch (i - M_post1_q3 + Num9_p1 + Num10_p1 + Num11_p1 + Num12_p1)
		{
		case 0:// 1 lei 2j shang
			kind1_node(j, i, -(Num13_p1 - N_d0 - 2) - 2, -(Num13_p1 - N_d0 - 2) - 1, -(Num13_p1 - N_d0 - 2), 0, 1,
				Num12_p1 + 4 - 2, Num12_p1 + 4 - 1, Num12_p1 + 4, A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case Node121_p1:// inner 3 node
			inner_3node(j, i, -(Num13_p1 - N_d0 - 2), -1, 0, 1, 1, k*k*dx*dx - 4, Coef, Coef_location);
			break;
		case (Node121_p1 + 1) :// inner 3 node
			inner_3node(j, i, -(Num13_p1 - N_d0 - 2 - 2), 0, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Num12_p1 / 2 - 1) :// 1 lei 2j xia
			kind1_node(j, i, -(Num13_p1 - N_d0 - 2 - 2), -(Num13_p1 - N_d0 - 2 - 2) + 1, -(Num13_p1 - N_d0 - 2 - 2) + 2, -1, 0,
			Num12_p1 + 4 - 2, Num12_p1 + 4 - 2 + 1, Num12_p1 + 4 - 2 + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		case (Num12_p1 / 2) :// 1 lei 2j shang
			kind1_node(j, i, -(Num13_p1 - Node_1313_p1 + 4 + Num12_p1 - 1) - 2 - 2, -(Num13_p1 - Node_1313_p1 + 4 + Num12_p1 - 1) - 2 - 1,
			-(Num13_p1 - Node_1313_p1 + 4 + Num12_p1 - 1) - 2, 0, 1, Num11_p1 - 4 + 2 - 2, Num11_p1 - 4 + 2 - 1, Num11_p1 - 4 + 2,
			A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case Node122_p1:// inner 3 node
			inner_3node(j, i, -(Num13_p1 - Node_1313_p1 + 4 + Num12_p1 - 1) - 2, -1, 0, 1, 1, k*k*dx*dx - 4, Coef, Coef_location);
			break;
		case (Node122_p1 + 1) :// inner 3 node
			inner_3node(j, i, -(Num13_p1 - Node_1313_p1 + 4 + Num12_p1 - 1), 0, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Num12_p1 - 1) :// 1 lei 2j xia
			kind1_node(j, i, -(Num13_p1 - Node_1313_p1 + 4 + Num12_p1 - 1), -(Num13_p1 - Node_1313_p1 + 4 + Num12_p1 - 1) + 1,
			-(Num13_p1 - Node_1313_p1 + 4 + Num12_p1 - 1) + 2, -1, 0, Num11_p1 - 4, Num11_p1 - 4 + 1, Num11_p1 - 4 + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - M_post1_q3 + Num9_p1 + Num10_p1 + Num11_p1 + Num12_p1) < Node121_p1)// 1dx
			{
				inner_5node(j, i, 1, -(Num13_p1 - N_d0 - 2), Num12_p1 + 4, Coef, Coef_location);
			}
			else if (((i - M_post1_q3 + Num9_p1 + Num10_p1 + Num11_p1 + Num12_p1) > (Node121_p1 + 1)) && ((i - M_post1_q3 + Num9_p1 + Num10_p1 + Num11_p1 + Num12_p1) < (Num12_p1 / 2 - 1)))// 2dx
			{
				inner_5node(j, i, 1, -(Num13_p1 - N_d0 - 2 - 2), Num12_p1 + 4 - 2, Coef, Coef_location);
			}
			else if (((i - M_post1_q3 + Num9_p1 + Num10_p1 + Num11_p1 + Num12_p1) > (Num12_p1 / 2)) && ((i - M_post1_q3 + Num9_p1 + Num10_p1 + Num11_p1 + Num12_p1) < Node122_p1))// 2dx
			{
				inner_5node(j, i, 1, -(Num13_p1 - Node_1313_p1 + 4 + Num12_p1 - 1) - 2, Num11_p1 - 4 + 2, Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_5node(j, i, 1, -(Num13_p1 - Node_1313_p1 + 4 + Num12_p1 - 1), Num11_p1 - 4, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!! num 27 !!!
	for (i = (M_post1_q3 - Num9_p1 - Num10_p1 - Num11_p1); i<(M_post1_q3 - Num9_p1 - Num10_p1 - Num11_p1 + Node112_p1 + 1); i++)// num 27(11) :  41 node! part1
	{
		switch (i - M_post1_q3 + Num9_p1 + Num10_p1 + Num11_p1)
		{
		case 0:// 1 lei 1j shang
			kind1_node(j, i, -(Num13_p1 - N_d0 + 2 + Num12_p1) - 2, -(Num13_p1 - N_d0 + 2 + Num12_p1) - 1, -(Num13_p1 - N_d0 + 2 + Num12_p1),
				0, 1, (Num11_p1 + Num10_p1 + N_d0 - 2) - 2, (Num11_p1 + Num10_p1 + N_d0 - 2) - 1, (Num11_p1 + Num10_p1 + N_d0 - 2),
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case 3:// 2 lei 2j shang
			kind2_node(j, i, -(Num13_p1 - N_d0 + 2 + Num12_p1), -2, -1, 0, 1, (Num11_p1 + Num10_p1 + N_d0 - 2),
				B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node111_p1://inner 4 node
			inner_4node(j, i, -(Num12_p1 + 4), -1, 0, Num11_p1 - 4, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node111_p1 + 1) ://inner 4 node
			inner_4node(j, i, -(Num12_p1 + 4) + 2, 0, 1, Num11_p1 - 4, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Node112_p1 - 3) :// 2 lei 2j xia
			kind2_node(j, i, -(Num13_p1 - Node_13_p1 + Num12_p1 + Node112_p1), -1, 0, 1, 2, (Num11_p1 - Node112_p1 + Num10_p1 + Node_9_p1),
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case Node112_p1:// 1 lei 1j xia
			kind1_node(j, i, -(Num13_p1 - Node_13_p1 + Num12_p1 + Node112_p1), -(Num13_p1 - Node_13_p1 + Num12_p1 + Node112_p1) + 1,
				-(Num13_p1 - Node_13_p1 + Num12_p1 + Node112_p1) + 2, -1, 0, (Num11_p1 - Node112_p1 + Num10_p1 + Node_9_p1),
				(Num11_p1 - Node112_p1 + Num10_p1 + Node_9_p1) + 1, (Num11_p1 - Node112_p1 + Num10_p1 + Node_9_p1) + 2,
				A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - M_post1_q3 + Num9_p1 + Num10_p1 + Num11_p1) < 3)// 2dx
			{
				inner_5node(j, i, 2, -(Num13_p1 - N_d0 + 2 + Num12_p1), (Num11_p1 + Num10_p1 + N_d0 - 2), Coef, Coef_location);
			}
			else if (((i - M_post1_q3 + Num9_p1 + Num10_p1 + Num11_p1) > 3) && ((i - M_post1_q3 + Num9_p1 + Num10_p1 + Num11_p1) < (Node111_p1)))// 1dx
			{
				inner_5node(j, i, 1, -(Num12_p1 + 4), Num11_p1 - 4, Coef, Coef_location);
			}
			else if (((i - M_post1_q3 + Num9_p1 + Num10_p1 + Num11_p1) > (Node111_p1 + 1)) && ((i - M_post1_q3 + Num9_p1 + Num10_p1 + Num11_p1) < (Node112_p1 - 3)))// 1dx
			{
				inner_5node(j, i, 1, -(Num12_p1 + 4) + 2, Num11_p1 - 4, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, i, 2, -(Num13_p1 - Node_13_p1 + Num12_p1 + Node112_p1), (Num11_p1 - Node112_p1 + Num10_p1 + Node_9_p1), Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post1_q3 - Num9_p1 - Num10_p1 - Num11_p1 + Node112_p1 + 1); i<(M_post1_q3 - Num9_p1 - Num10_p1 - Num11_p1 + Node113_p1); i++)// num 27(11) :  41 node! part2
	{
		switch (i - M_post1_q3 + Num9_p1 + Num10_p1 + Num11_p1 - Node112_p1 - 1)
		{
		case 0:// 3 lei 1j shang
			kind3_node(j, i, -(Num13_p1 - Node132_p1 - 1 + Num12_p1 + Node112_p1) - 2, -(Num13_p1 - Node132_p1 - 1 + Num12_p1 + Node112_p1),
				0, 1, (Num11_p1 - Node112_p1 - 1 + Num10_p1 + Node92_p1) - 1, (Num11_p1 - Node112_p1 - 1 + Num10_p1 + Node92_p1),
				C4, C3_1, C1_1, C2_1, C4, C3_1, Coef, Coef_location);
			break;
		case 1:// 1 lei 1j you
			kind1_node(j, i, -(Num13_p1 - Node132_p1 - 1 + Num12_p1 + Node112_p1), -1, 0, 1,
				(Num11_p1 - Node112_p1 - 1 + Num10_p1 + Node92_p1) - 1, (Num11_p1 - Node112_p1 - 1 + Num10_p1 + Node92_p1),
				(Num11_p1 - Node112_p1 - 1 + Num10_p1 + Node92_p1) - 1 + (Num9_p1 - Node_9_p1 + Num8_p1 + Num7_p1 + Num6_p1 + Node_5_p1),
				(Num11_p1 - Node112_p1 - 1 + Num10_p1 + Node92_p1) + (Num9_p1 - Node_9_p1 + Num8_p1 + Num7_p1 + Num6_p1 + Node_5_p1),
				A4_1, A3_1, A5_1, A3_1, A1_1, A1_1, A2_1, A2_1, Coef, Coef_location);
			break;
		case 2:// inner 5 node
			inner_5node(j, i, 2, -(Num13_p1 - Node132_p1 - 1 + Num12_p1 + Node112_p1),
				(Num11_p1 - Node112_p1 - 1 + Num10_p1 + Node92_p1) - 1, Coef, Coef_location);
			break;
		case 3:// 1 lei 1j zuo
			kind1_node(j, i, -(Num13_p1 - Node132_p1 - 1 + Num12_p1 + Node112_p1), -1, 0, 1,
				(Num11_p1 - Node112_p1 - 1 + Num10_p1 + Node92_p1) - 2, (Num11_p1 - Node112_p1 - 1 + Num10_p1 + Node92_p1) - 1,
				(Num11_p1 - Node112_p1 - 1 + Num10_p1 + Node92_p1) - 2 + (Num9_p1 - Node_9_p1 + Num8_p1 + Num7_p1 + Num6_p1 + Node_5_p1),
				(Num11_p1 - Node112_p1 - 1 + Num10_p1 + Node92_p1) - 1 + (Num9_p1 - Node_9_p1 + Num8_p1 + Num7_p1 + Num6_p1 + Node_5_p1),
				A4_1, A3_1, A5_1, A3_1, A1_1, A1_1, A2_1, A2_1, Coef, Coef_location);
			break;
		case 4:// inner 5 node
			inner_5node(j, i, 2, -(Num13_p1 - Node132_p1 - 1 + Num12_p1 + Node112_p1),
				(Num11_p1 - Node112_p1 - 1 + Num10_p1 + Node92_p1) - 2, Coef, Coef_location);
			break;
		case 5:// 1 lei 1j zuo
			kind1_node(j, i, -(Num13_p1 - Node132_p1 - 1 + Num12_p1 + Node112_p1), -1, 0, 1,
				(Num11_p1 - Node112_p1 - 1 + Num10_p1 + Node92_p1) - 3, (Num11_p1 - Node112_p1 - 1 + Num10_p1 + Node92_p1) - 2,
				(Num11_p1 - Node112_p1 - 1 + Num10_p1 + Node92_p1) - 3 + (Num9_p1 - Node_9_p1 + Num8_p1 + Num7_p1 + Num6_p1 + Node_5_p1),
				(Num11_p1 - Node112_p1 - 1 + Num10_p1 + Node92_p1) - 2 + (Num9_p1 - Node_9_p1 + Num8_p1 + Num7_p1 + Num6_p1 + Node_5_p1),
				A4_1, A3_1, A5_1, A3_1, A1_1, A1_1, A2_1, A2_1, Coef, Coef_location);
			break;
		case 6:// 3 lei 1j xia
			kind3_node(j, i, -(Num13_p1 - Node132_p1 - 1 + Num12_p1 + Node112_p1), -(Num13_p1 - Node132_p1 - 1 + Num12_p1 + Node112_p1) + 2,
				-1, 0, (Num11_p1 - Node112_p1 - 1 + Num10_p1 + Node92_p1) - 3, (Num11_p1 - Node112_p1 - 1 + Num10_p1 + Node92_p1) - 3 + 1,
				C3_1, C4, C2_1, C1_1, C3_1, C4, Coef, Coef_location);
			break;
		}
		j = j + 1;
	}
	for (i = (M_post1_q3 - Num9_p1 - Num10_p1 - Num11_p1 + Node113_p1); i<(M_post1_q3 - Num9_p1 - Num10_p1); i++)// num 27(11) :  41 node! part3
	{
		switch (i - M_post1_q3 + Num9_p1 + Num10_p1 + Num11_p1)
		{
		case Node113_p1:// 1 lei 1j shang
			kind1_node(j, i, -(Num13_p1 - Node133_p1 - 4 + Num12_p1 + Node113_p1) - 2, -(Num13_p1 - Node133_p1 - 4 + Num12_p1 + Node113_p1) - 1,
				-(Num13_p1 - Node133_p1 - 4 + Num12_p1 + Node113_p1), 0, 1, (Num11_p1 - Node113_p1 + Num10_p1 + Node93_p1 + 5) - 2,
				(Num11_p1 - Node113_p1 + Num10_p1 + Node93_p1 + 5) - 1, (Num11_p1 - Node113_p1 + Num10_p1 + Node93_p1 + 5),
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case (Node113_p1 + 3) :// 2 lei 2j shang
			kind2_node(j, i, -(Num13_p1 - Node133_p1 - 4 + Num12_p1 + Node113_p1), -2, -1, 0, 1, (Num11_p1 - Node113_p1 + Num10_p1 + Node93_p1 + 5),
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node114_p1://inner 4 node
			inner_4node(j, i, -(Num11_p1 - 5 + 1) - 2, -1, 0, 5 + Num10_p1 - 1, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node114_p1 + 1) ://inner 4 node
			inner_4node(j, i, -(Num11_p1 - 5 + 1), 0, 1, 5 + Num10_p1 - 1, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Num11_p1 - 4) :// 2 lei 2j xia
			kind2_node(j, i, -(Num13_p1 - Node_1313_p1 + Num12_p1 + Num11_p1 - 1), -1, 0, 1, 2, (1 + Num10_p1 + Node_99_p1),
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num11_p1 - 1) :// 1 lei 1j xia
			kind1_node(j, i, -(Num13_p1 - Node_1313_p1 + Num12_p1 + Num11_p1 - 1), -(Num13_p1 - Node_1313_p1 + Num12_p1 + Num11_p1 - 1) + 1,
			-(Num13_p1 - Node_1313_p1 + Num12_p1 + Num11_p1 - 1) + 2, -1, 0,
			(1 + Num10_p1 + Node_99_p1), (1 + Num10_p1 + Node_99_p1) + 1, (1 + Num10_p1 + Node_99_p1) + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - M_post1_q3 + Num9_p1 + Num10_p1 + Num11_p1) < (Node113_p1 + 3))// 2dx
			{
				inner_5node(j, i, 2, -(Num13_p1 - Node133_p1 - 4 + Num12_p1 + Node113_p1), (Num11_p1 - Node113_p1 + Num10_p1 + Node93_p1 + 5), Coef, Coef_location);
			}
			else if (((i - M_post1_q3 + Num9_p1 + Num10_p1 + Num11_p1) > (Node113_p1 + 3)) && ((i - M_post1_q3 + Num9_p1 + Num10_p1 + Num11_p1) < (Node114_p1)))// 1dx
			{
				inner_5node(j, i, 1, -(Num11_p1 - 5 + 1) - 2, 5 + Num10_p1 - 1, Coef, Coef_location);
			}
			else if (((i - M_post1_q3 + Num9_p1 + Num10_p1 + Num11_p1) > (Node114_p1 + 1)) && ((i - M_post1_q3 + Num9_p1 + Num10_p1 + Num11_p1) < (Num11_p1 - 4)))// 1dx
			{
				inner_5node(j, i, 1, -(Num11_p1 - 5 + 1), 5 + Num10_p1 - 1, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, i, 2, -(Num13_p1 - Node_1313_p1 + Num12_p1 + Num11_p1 - 1), (1 + Num10_p1 + Node_99_p1), Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!! num 28 !!!
	for (i = (M_post1_q3 - Num9_p1 - Num10_p1); i<(M_post1_q3 - Num9_p1); i++)// num 28(10) :  12 node! 
	{
		switch (i - M_post1_q3 + Num9_p1 + Num10_p1)
		{
		case 0:// 1 lei shang 2j
			kind1_node(j, i, -(Num11_p1 - 4) - 2, -(Num11_p1 - 4) - 1, -(Num11_p1 - 4), 0, 1, Num10_p1 + N_d0 + 2 - 2, Num10_p1 + N_d0 + 2 - 1,
				Num10_p1 + N_d0 + 2, A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case 1:
			inner_5node(j, i, 1, -(Num11_p1 - 4), Num10_p1 + N_d0 + 2, Coef, Coef_location);
			break;
		case 2://inner 4 node
			inner_4node(j, i, -(Num11_p1 - 4), -1, 0, Num10_p1 + N_d0 + 2, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case 3://inner 4 node
			inner_4node(j, i, -(Num11_p1 - 4), 0, 1, Num10_p1 + N_d0 + 2, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case 4:
			inner_5node(j, i, 1, -(Num11_p1 - 4), Num10_p1 + N_d0 + 2, Coef, Coef_location);
			break;
		case (Num10_p1 / 2 - 1) :// 1 lei xia 2j
			kind1_node(j, i, -(Num11_p1 - 4), -(Num11_p1 - 4) + 1, -(Num11_p1 - 4) + 2, -1, 0, Num10_p1 + N_d0 + 2, Num10_p1 + N_d0 + 2 + 1,
			Num10_p1 + N_d0 + 2 + 2, A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		case (Num10_p1 / 2) :// 1 lei shang 2j
			kind1_node(j, i, -(5 + Num10_p1 - 1) - 2, -(5 + Num10_p1 - 1) - 1, -(5 + Num10_p1 - 1), 0, 1,
			Node_99_p1 - 4 + 1 - 2, Node_99_p1 - 4 + 1 - 1, Node_99_p1 - 4 + 1,
			A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case 7:
			inner_5node(j, i, 1, -(5 + Num10_p1 - 1), Node_99_p1 - 4 + 1, Coef, Coef_location);
			break;
		case 8://inner 4 node
			inner_4node(j, i, -(5 + Num10_p1 - 1), -1, 0, Node_99_p1 - 4 + 1, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case 9://inner 4 node
			inner_4node(j, i, -(5 + Num10_p1 - 1), 0, 1, Node_99_p1 - 4 + 1, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case 10:
			inner_5node(j, i, 1, -(5 + Num10_p1 - 1), Node_99_p1 - 4 + 1, Coef, Coef_location);
			break;
		case (Num10_p1 - 1) :// 1 lei xia 2j
			kind1_node(j, i, -(5 + Num10_p1 - 1), -(5 + Num10_p1 - 1) + 1, -(5 + Num10_p1 - 1) + 2, -1, 0,
			Node_99_p1 - 4 + 1, Node_99_p1 - 4 + 1 + 1, Node_99_p1 - 4 + 1 + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!! num 29 !!!
	for (i = (M_post1_q3 - Num9_p1); i<(M_post1_q3 - Num9_p1 + Node_9_p1 + 2); i++)// circle num 29(9) :  72 node!  part1 !
	{
		switch (i - M_post1_q3 + Num9_p1)
		{
		case 0:
			boundary_node(j, i, 1, -Num10_p1 - Num11_p1 - Num12_p1 - Num13_p1, Num9_p1 + Num8_p1 + Num7_p1 + Num6_p1, Coef, Coef_location);
			break;
		case (N_d0 - 3) ://2 lei shang 1j
			kind2_node(j, i, -Num10_p1 - Num11_p1 - Num12_p1 - Num13_p1, -2, -1, 0, 1, Num9_p1 + Num8_p1 + Num7_p1 + Num6_p1,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (N_d0 + 1) :// 2 lei 2j shang
			kind2_node(j, i, -(Num11_p1 + Num10_p1 + N_d0 - 2), -2, -1, 0, 1, Num9_p1 - N_d0 + 2 + Num8_p1,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node91_p1:// inner 4 node
			inner_4node(j, i, -(Num10_p1 + N_d0 + 2), -1, 0, Num9_p1 - N_d0 - 2, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node91_p1 + 1) :// inner 4 node
			inner_4node(j, i, -(Num10_p1 + N_d0 + 2), 0, 1, Num9_p1 - N_d0 - 2, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Node_9_p1 - 3) ://2 lei xia 2j
			kind2_node(j, i, -(Num11_p1 - Node112_p1 + Num10_p1 + Node_9_p1), -1, 0, 1, 2, Num9_p1 - Node_9_p1 + Num8_p1 + Num7_p1 / 2 - 1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_9_p1 + 1) ://2 lei xia 1j
			kind2_node(j, i, -(Num13_p1 - Node_13_p1 + Num12_p1 + Num11_p1 + Num10_p1 + Node_9_p1), -1, 0, 1, 2,
			Num9_p1 - Node_9_p1 + Num8_p1 + Num7_p1 + Num6_p1 + Node_5_p1, B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - M_post1_q3 + Num9_p1) < (N_d0 - 3))// 4dx!
			{
				inner_5node(j, i, 4, -Num10_p1 - Num11_p1 - Num12_p1 - Num13_p1, Num9_p1 + Num8_p1 + Num7_p1 + Num6_p1, Coef, Coef_location);
			}
			else if (((i - M_post1_q3 + Num9_p1) > (N_d0 - 3)) && ((i - M_post1_q3 + Num9_p1) < (N_d0 + 1)))// 2dx!
			{
				inner_5node(j, i, 2, -(Num11_p1 + Num10_p1 + N_d0 - 2), Num9_p1 - N_d0 + 2 + Num8_p1, Coef, Coef_location);
			}
			else if (((i - M_post1_q3 + Num9_p1) > (N_d0 + 1)) && ((i - M_post1_q3 + Num9_p1) < Node91_p1))// 1dx
			{
				inner_5node(j, i, 1, -(Num10_p1 + N_d0 + 2), Num9_p1 - N_d0 - 2, Coef, Coef_location);
			}
			else if (((i - M_post1_q3 + Num9_p1) > (Node91_p1 + 1)) && ((i - M_post1_q3 + Num9_p1) < (Node_9_p1 - 3)))// 1dx
			{
				inner_5node(j, i, 1, -(Num10_p1 + N_d0 + 2), Num9_p1 - N_d0 - 2, Coef, Coef_location);
			}
			else// 2dx!
			{
				inner_5node(j, i, 2, -(Num11_p1 - Node112_p1 + Num10_p1 + Node_9_p1), Num9_p1 - Node_9_p1 + Num8_p1 + Num7_p1 / 2 - 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post1_q3 - Num9_p1 + Node_9_p1 + 2); i<(M_post1_q3 - Num9_p1 + Node93_p1 + 4); i++)// circle num 29(9) :  72 node!  part2 !
	{
		switch (i - M_post1_q3 + Num9_p1 - Node_9_p1 - 2)
		{
		default:
			if ((i - M_post1_q3 + Num9_p1 - Node_9_p1 - 2) < 3)
			{
				inner_5node(j, i, 4, -(Num13_p1 - Node_13_p1 + Num12_p1 + Num11_p1 + Num10_p1 + Node_9_p1),
					Num9_p1 - Node_9_p1 + Num8_p1 + Num7_p1 + Num6_p1 + Node_5_p1, Coef, Coef_location);
			}
			else if (((i - M_post1_q3 + Num9_p1 - Node_9_p1 - 2) >= 3) && ((i - M_post1_q3 + Num9_p1 - Node_9_p1 - 2) <= 6))// 2lei you 1j
			{
				kind2_node(j, i, -(Num11_p1 - Node112_p1 - 1 - (i - M_post1_q3 + Num9_p1 - Node_9_p1 - 5) * 2 + Num10_p1 + (i - M_post1_q3 + Num9_p1)),
					-1, 0, 1, Num9_p1 - Node_9_p1 + Num8_p1 + Num7_p1 + Num6_p1 + Node_5_p1,
					Num9_p1 - Node_9_p1 + Num8_p1 + Num7_p1 + Num6_p1 + Num5_p1 + Num4_p1 + Num3_p1 + Num2_p1 + Node_1_p1,
					B2, B3, B1_1, B3, B4, B5, Coef, Coef_location);
			}
			else
			{
				inner_5node(j, i, 4, -(Num13_p1 - Node133_p1 + Num12_p1 + Num11_p1 + Num10_p1 + Node93_p1 + 1),
					Num9_p1 - Node_9_p1 + Num8_p1 + Num7_p1 + Num6_p1 + Node_5_p1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post1_q3 - Num9_p1 + Node93_p1 + 4); i<(M_post1_q3); i++)// circle num 29(9) :  79 node!  part3 !
	{
		switch (i - M_post1_q3 + Num9_p1)
		{
		case (Num9_p1 - 1) :
			boundary_node(j, i, 2, -Num12_p1 - Num11_p1 - Num10_p1 - Num9_p1, Num8_p1 + Num7_p1 + Num6_p1 + Num5_p1, Coef, Coef_location);
			break;
		case (Node93_p1 + 4) ://2 lei shang 1j
			kind2_node(j, i, -(Num13_p1 - Node133_p1 + Num12_p1 + Num11_p1 + Num10_p1 + Node93_p1 + 1), -2, -1, 0, 1,
			Num9_p1 - Node_9_p1 + Num8_p1 + Num7_p1 + Num6_p1 + Node_5_p1, B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node94_p1 - 3) :// 2 lei 2j shang
			kind2_node(j, i, -(Num11_p1 - Node113_p1 + Num10_p1 + Node93_p1 + 5), -2, -1, 0, 1,
			Num9_p1 - Node93_p1 - 5 + Num8_p1 + Num7_p1 / 2, B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node94_p1:// inner 4 node
			inner_4node(j, i, -(Node_99_p1 - 4 + 1), -1, 0, Num9_p1 - Node_99_p1 + 4 + Num8_p1 - 1, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node94_p1 + 1) :// inner 4 node
			inner_4node(j, i, -(Node_99_p1 - 4 + 1), 0, 1, Num9_p1 - Node_99_p1 + 4 + Num8_p1 - 1, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Node_99_p1 - 3) ://2 lei xia 2j
			kind2_node(j, i, -(Num10_p1 + Node_99_p1 + 1), -1, 0, 1, 2, Num9_p1 - Node_99_p1 + Num8_p1 + Num7_p1 - 1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_99_p1 + 1) ://2 lei xia 1j
			kind2_node(j, i, -Num12_p1 - Num11_p1 - Num10_p1 - Num9_p1, -1, 0, 1, 2, Num8_p1 + Num7_p1 + Num6_p1 + Num5_p1,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - M_post1_q3 + Num9_p1) < (Node94_p1 - 3))// 2dx!
			{
				inner_5node(j, i, 2, -(Num11_p1 - Node113_p1 + Num10_p1 + Node93_p1 + 5), Num9_p1 - Node93_p1 - 5 + Num8_p1 + Num7_p1 / 2, Coef, Coef_location);
			}
			else if (((i - M_post1_q3 + Num9_p1) > (Node94_p1 - 3)) && ((i - M_post1_q3 + Num9_p1) < Node94_p1))// 1dx
			{
				inner_5node(j, i, 1, -(Node_99_p1 - 4 + 1), Num9_p1 - Node_99_p1 + 4 + Num8_p1 - 1, Coef, Coef_location);
			}
			else if (((i - M_post1_q3 + Num9_p1) > (Node94_p1 + 1)) && ((i - M_post1_q3 + Num9_p1) < (Node_99_p1 - 3)))// 1dx
			{
				inner_5node(j, i, 1, -(Node_99_p1 - 4 + 1), Num9_p1 - Node_99_p1 + 4 + Num8_p1 - 1, Coef, Coef_location);
			}
			else if (((i - M_post1_q3 + Num9_p1) > (Node_99_p1 - 3)) && ((i - M_post1_q3 + Num9_p1) < (Node_99_p1 + 1)))// 2dx
			{
				inner_5node(j, i, 2, -(Num10_p1 + Node_99_p1 + 1), Num9_p1 - Node_99_p1 + Num8_p1 + Num7_p1 - 1, Coef, Coef_location);
			}
			else// 4dx!
			{
				inner_5node(j, i, 4, -Num12_p1 - Num11_p1 - Num10_p1 - Num9_p1, Num8_p1 + Num7_p1 + Num6_p1 + Num5_p1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!! num 30 !!!
	for (i = (M_post1_q3); i<(M_post1_q3 + Num8_p1); i++)// num 30(8) :  12 node! 
	{
		switch (i - M_post1_q3)
		{
		case 0:// 1 lei shang 2j
			kind1_node(j, i, -(Num9_p1 - N_d0 - 2) - 2, -(Num9_p1 - N_d0 - 2) - 1, -(Num9_p1 - N_d0 - 2), 0, 1, Num8_p1 + 4 - 2,
				Num8_p1 + 4 - 1, Num8_p1 + 4, A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case 1:
			inner_5node(j, i, 1, -(Num9_p1 - N_d0 - 2), Num8_p1 + 4, Coef, Coef_location);
			break;
		case 2://inner 4 node
			inner_4node(j, i, -(Num9_p1 - N_d0 - 2), -1, 0, Num8_p1 + 4, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case 3://inner 4 node
			inner_4node(j, i, -(Num9_p1 - N_d0 - 2), 0, 1, Num8_p1 + 4, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case 4:
			inner_5node(j, i, 1, -(Num9_p1 - N_d0 - 2), Num8_p1 + 4, Coef, Coef_location);
			break;
		case (Num8_p1 / 2 - 1) :// 1 lei xia 2j
			kind1_node(j, i, -(Num9_p1 - N_d0 - 2), -(Num9_p1 - N_d0 - 2) + 1, -(Num9_p1 - N_d0 - 2) + 2, -1, 0, Num8_p1 + 4,
			Num8_p1 + 4 + 1, Num8_p1 + 4 + 2, A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		case (Num8_p1 / 2) :// 1 lei shang 2j
			kind1_node(j, i, -(Num9_p1 - Node_99_p1 + 4 + Num8_p1 - 1) - 2, -(Num9_p1 - Node_99_p1 + 4 + Num8_p1 - 1) - 1,
			-(Num9_p1 - Node_99_p1 + 4 + Num8_p1 - 1), 0, 1, Num7_p1 - 4 - 2, Num7_p1 - 4 - 1, Num7_p1 - 4,
			A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case 7:
			inner_5node(j, i, 1, -(Num9_p1 - Node_99_p1 + 4 + Num8_p1 - 1), Num7_p1 - 4, Coef, Coef_location);
			break;
		case 8://inner 4 node
			inner_4node(j, i, -(Num9_p1 - Node_99_p1 + 4 + Num8_p1 - 1), -1, 0, Num7_p1 - 4, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case 9://inner 4 node
			inner_4node(j, i, -(Num9_p1 - Node_99_p1 + 4 + Num8_p1 - 1), 0, 1, Num7_p1 - 4, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case 10:
			inner_5node(j, i, 1, -(Num9_p1 - Node_99_p1 + 4 + Num8_p1 - 1), Num7_p1 - 4, Coef, Coef_location);
			break;
		case (Num8_p1 - 1) :// 1 lei xia 2j
			kind1_node(j, i, -(Num9_p1 - Node_99_p1 + 4 + Num8_p1 - 1), -(Num9_p1 - Node_99_p1 + 4 + Num8_p1 - 1) + 1,
			-(Num9_p1 - Node_99_p1 + 4 + Num8_p1 - 1) + 2, -1, 0, Num7_p1 - 4, Num7_p1 - 4 + 1, Num7_p1 - 4 + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!! num 31 !!!
	for (i = (M_post1_q3 + Num8_p1); i<(M_post1_q3 + Num8_p1 + Num7_p1 / 2); i++)// num 31(7) :  28 node! part1
	{
		switch (i - M_post1_q3 - Num8_p1)
		{
		case 0:// 1 lei 1j shang
			kind1_node(j, i, -(Num9_p1 - N_d0 + 2 + Num8_p1) - 2, -(Num9_p1 - N_d0 + 2 + Num8_p1) - 1, -(Num9_p1 - N_d0 + 2 + Num8_p1),
				0, 1, Num7_p1 + Num6_p1 + N_d0 - 2 - 2, Num7_p1 + Num6_p1 + N_d0 - 2 - 1, Num7_p1 + Num6_p1 + N_d0 - 2,
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case 3:// 2 lei 2j shang
			kind2_node(j, i, -(Num9_p1 - N_d0 + 2 + Num8_p1), -2, -1, 0, 1, Num7_p1 + Num6_p1 + N_d0 - 2,
				B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node71_p1://inner 4 node
			inner_4node(j, i, -(Num8_p1 + 4), -1, 0, Num7_p1 - 4, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node71_p1 + 1) ://inner 4 node
			inner_4node(j, i, -(Num8_p1 + 4), 0, 1, Num7_p1 - 4 + 2, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Num7_p1 / 2 - 4) :// 2 lei 2j xia
			kind2_node(j, i, -(Num9_p1 - Node_9_p1 + Num8_p1 + Num7_p1 / 2 - 1), -1, 0, 1, 2, Num7_p1 / 2 + 1 + Num6_p1 + Node_5_p1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num7_p1 / 2 - 1) :// 1 lei 1j xia
			kind1_node(j, i, -(Num9_p1 - Node_9_p1 + Num8_p1 + Num7_p1 / 2 - 1), -(Num9_p1 - Node_9_p1 + Num8_p1 + Num7_p1 / 2 - 1) + 1,
			-(Num9_p1 - Node_9_p1 + Num8_p1 + Num7_p1 / 2 - 1) + 2, -1, 0, Num7_p1 / 2 + 1 + Num6_p1 + Node_5_p1,
			Num7_p1 / 2 + 1 + Num6_p1 + Node_5_p1 + 1, Num7_p1 / 2 + 1 + Num6_p1 + Node_5_p1 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - M_post1_q3 - Num8_p1) < 3)// 2dx
			{
				inner_5node(j, i, 2, -(Num9_p1 - N_d0 + 2 + Num8_p1), Num7_p1 + Num6_p1 + N_d0 - 2, Coef, Coef_location);
			}
			else if (((i - M_post1_q3 - Num8_p1) > 3) && ((i - M_post1_q3 - Num8_p1) < (Node71_p1)))// 1dx
			{
				inner_5node(j, i, 1, -(Num8_p1 + 4), Num7_p1 - 4, Coef, Coef_location);
			}
			else if (((i - M_post1_q3 - Num8_p1) > (Node71_p1 + 1)) && ((i - M_post1_q3 - Num8_p1) < (Num7_p1 / 2 - 4)))// 1dx
			{
				inner_5node(j, i, 1, -(Num8_p1 + 4), Num7_p1 - 4 + 2, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, i, 2, -(Num9_p1 - Node_9_p1 + Num8_p1 + Num7_p1 / 2 - 1), Num7_p1 / 2 + 1 + Num6_p1 + Node_5_p1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post1_q3 + Num8_p1 + Num7_p1 / 2); i<(M_post1_q3 + Num8_p1 + Num7_p1); i++)// num 11(7) :  28 node! part2
	{
		switch (i - M_post1_q3 - Num8_p1)
		{
		case (Num7_p1 / 2) :// 1 lei 1j shang
			kind1_node(j, i, -(Num9_p1 - Node_9_p1 - N_a + 3 + Num8_p1 + Num7_p1 / 2) - 2, -(Num9_p1 - Node_9_p1 - N_a + 3 + Num8_p1 + Num7_p1 / 2) - 1,
			-(Num9_p1 - Node_9_p1 - N_a + 3 + Num8_p1 + Num7_p1 / 2), 0, 1, Num7_p1 / 2 + Num6_p1 + Node_5_p1 + N_a - 3 - 2,
			Num7_p1 / 2 + Num6_p1 + Node_5_p1 + N_a - 3 - 1, Num7_p1 / 2 + Num6_p1 + Node_5_p1 + N_a - 3,
			A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case (Num7_p1 / 2 + 3) :// 2 lei 2j shang
			kind2_node(j, i, -(Num9_p1 - Node_9_p1 - N_a + 3 + Num8_p1 + Num7_p1 / 2), -2, -1, 0, 1, Num7_p1 / 2 + Num6_p1 + Node_5_p1 + N_a - 3,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node72_p1://inner 4 node
			inner_4node(j, i, -(Num7_p1 - 5 + 1), -1, 0, 5 + Num6_p1 - 1 - 2, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node72_p1 + 1) ://inner 4 node
			inner_4node(j, i, -(Num7_p1 - 5 + 1), 0, 1, 5 + Num6_p1 - 1, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Num7_p1 - 4) :// 2 lei 2j xia
			kind2_node(j, i, -(Num9_p1 - Node_99_p1 + Num8_p1 + Num7_p1 - 1), -1, 0, 1, 2, 1 + Num6_p1 + Node_55_p1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num7_p1 - 1) :// 1 lei 1j xia
			kind1_node(j, i, -(Num9_p1 - Node_99_p1 + Num8_p1 + Num7_p1 - 1), -(Num9_p1 - Node_99_p1 + Num8_p1 + Num7_p1 - 1) + 1,
			-(Num9_p1 - Node_99_p1 + Num8_p1 + Num7_p1 - 1) + 2, -1, 0, 1 + Num6_p1 + Node_55_p1,
			1 + Num6_p1 + Node_55_p1 + 1, 1 + Num6_p1 + Node_55_p1 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - M_post1_q3 - Num8_p1) < (Num7_p1 / 2 + 3))// 2dx
			{
				inner_5node(j, i, 2, -(Num9_p1 - Node_9_p1 - N_a + 3 + Num8_p1 + Num7_p1 / 2), Num7_p1 / 2 + Num6_p1 + Node_5_p1 + N_a - 3, Coef, Coef_location);
			}
			else if (((i - M_post1_q3 - Num8_p1) > (Num7_p1 / 2 + 3)) && ((i - M_post1_q3 - Num8_p1) < (Node72_p1)))// 1dx
			{
				inner_5node(j, i, 1, -(Num7_p1 - 5 + 1), 5 + Num6_p1 - 1 - 2, Coef, Coef_location);
			}
			else if (((i - M_post1_q3 - Num8_p1) > (Node72_p1 + 1)) && ((i - M_post1_q3 - Num8_p1) < (Num7_p1 - 4)))// 1dx
			{
				inner_5node(j, i, 1, -(Num7_p1 - 5 + 1), 5 + Num6_p1 - 1, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, i, 2, -(Num9_p1 - Node_99_p1 + Num8_p1 + Num7_p1 - 1), 1 + Num6_p1 + Node_55_p1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (M_post1_q3 + Num8_p1 + Num7_p1); i<(M_post1_q3 + Num8_p1 + Num7_p1 + Num6_p1); i++)// num 12(6) :  16 node!
	{
		switch (i - M_post1_q3 - Num8_p1 - Num7_p1)
		{
		case 0:// 1 lei 2j shang
			kind1_node(j, i, -(Num7_p1 - 4) - 2, -(Num7_p1 - 4) - 1, -(Num7_p1 - 4), 0, 1,
				Num6_p1 + N_d0 + 2 - 2, Num6_p1 + N_d0 + 2 - 1, Num6_p1 + N_d0 + 2,
				A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case Node61_p1:// inner 3 node
			inner_3node(j, i, -1, 0, (Num6_p1 + N_d0 + 2), 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node61_p1 + 1) :// inner 3 node
			inner_3node(j, i, 0, 1, (Num6_p1 + N_d0 + 2) + 2, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Num6_p1 / 2 - 1) :// 1 lei 2j xia
			kind1_node(j, i, (-(Num7_p1 - 4) - 2), (-(Num7_p1 - 4) - 2) + 1, (-(Num7_p1 - 4) - 2) + 2, -1, 0,
			(Num6_p1 + N_d0 + 2) + 2, (Num6_p1 + N_d0 + 2) + 2 + 1, (Num6_p1 + N_d0 + 2) + 2 + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		case (Num6_p1 / 2) :// 1 lei 2j shang
			kind1_node(j, i, -(5 + Num6_p1 - 1) + 2 - 2, -(5 + Num6_p1 - 1) + 2 - 1, -(5 + Num6_p1 - 1) + 2, 0, 1,
			(Node_55_p1 - 4 + 1) - 2 - 2, (Node_55_p1 - 4 + 1) - 2 - 1, (Node_55_p1 - 4 + 1) - 2,
			A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case Node62_p1:// inner 3 node
			inner_3node(j, i, -1, 0, (Node_55_p1 - 4 + 1) - 2, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node62_p1 + 1) :// inner 3 node
			inner_3node(j, i, 0, 1, (Node_55_p1 - 4 + 1), k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Num6_p1 - 1) :// 1 lei 2j xia
			kind1_node(j, i, -(5 + Num6_p1 - 1), -(5 + Num6_p1 - 1) + 1, -(5 + Num6_p1 - 1) + 2, -1, 0,
			(Node_55_p1 - 4 + 1), (Node_55_p1 - 4 + 1) + 1, (Node_55_p1 - 4 + 1) + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - M_post1_q3 - Num8_p1 - Num7_p1) < Node61_p1)// 1dx
			{
				inner_5node(j, i, 1, -(Num7_p1 - 4), Num6_p1 + N_d0 + 2, Coef, Coef_location);
			}
			else if (((i - M_post1_q3 - Num8_p1 - Num7_p1) > (Node61_p1 + 1)) && ((i - M_post1_q3 - Num8_p1 - Num7_p1) < (Num6_p1 / 2 - 1)))// 2dx
			{
				inner_5node(j, i, 1, -(Num7_p1 - 4) - 2, (Num6_p1 + N_d0 + 2) + 2, Coef, Coef_location);
			}
			else if (((i - M_post1_q3 - Num8_p1 - Num7_p1) > (Num6_p1 / 2)) && ((i - M_post1_q3 - Num8_p1 - Num7_p1) < Node62_p1))// 2dx
			{
				inner_5node(j, i, 1, -(5 + Num6_p1 - 1) + 2, (Node_55_p1 - 4 + 1) - 2, Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_5node(j, i, 1, -(5 + Num6_p1 - 1), (Node_55_p1 - 4 + 1), Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (M_post1_q3 + Num8_p1 + Num7_p1 + Num6_p1); i<(M_post1_q3 + Num8_p1 + Num7_p1 + Num6_p1 + Node_5_p1 + 2); i++)// circle num13(5) :  83 node!  part1 !
	{
		switch (i - M_post1_q3 - Num8_p1 - Num7_p1 - Num6_p1)
		{
		case 0:
			boundary_node(j, i, 1, -Num9_p1 - Num8_p1 - Num7_p1 - Num6_p1, Num5_p1 + Num4_p1 + Num3_p1 + Num2_p1, Coef, Coef_location);
			break;
		case (N_d0 - 3) ://2 lei shang 1j
			kind2_node(j, i, -Num9_p1 - Num8_p1 - Num7_p1 - Num6_p1, -2, -1, 0, 1, Num5_p1 + Num4_p1 + Num3_p1 + Num2_p1,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (N_d0 + 1) :// 2 lei 2j shang
			kind2_node(j, i, -(Num7_p1 + Num6_p1 + N_d0 - 2), -2, -1, 0, 1, Num5_p1 - N_d0 + 2 + Num4_p1,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node51_p1:// inner 3 node
			inner_3node(j, i, -1, 0, (Num5_p1 - N_d0 - 2), 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node51_p1 + 1) :// inner 3 node
			inner_3node(j, i, 0, 1, (Num5_p1 - N_d0 - 2) + N_c1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Node_5_p1 - 3) ://2 lei xia 2j
			kind2_node(j, i, -(Num7_p1 / 2 + 1 + Num6_p1 + Node_5_p1), -1, 0, 1, 2, Num5_p1 - Node_5_p1 + Num4_p1 + Num3_p1 / 2 - 1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_5_p1 + 1) ://2 lei xia 1j
			kind2_node(j, i, -(Num9_p1 - Node_9_p1 + Num8_p1 + Num7_p1 + Num6_p1 + Node_5_p1), -1, 0, 1, 2,
			Num5_p1 - Node_5_p1 + Num4_p1 + Num3_p1 + Num2_p1 + Node_1_p1,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - M_post1_q3 - Num8_p1 - Num7_p1 - Num6_p1) < (N_d0 - 3))// 4dx!
			{
				inner_5node(j, i, 4, -Num9_p1 - Num8_p1 - Num7_p1 - Num6_p1, Num5_p1 + Num4_p1 + Num3_p1 + Num2_p1, Coef, Coef_location);
			}
			else if (((i - M_post1_q3 - Num8_p1 - Num7_p1 - Num6_p1) > (N_d0 - 3)) && ((i - M_post1_q3 - Num8_p1 - Num7_p1 - Num6_p1) < (N_d0 + 1)))// 2dx!
			{
				inner_5node(j, i, 2, -(Num7_p1 + Num6_p1 + N_d0 - 2), Num5_p1 - N_d0 + 2 + Num4_p1, Coef, Coef_location);
			}
			else if (((i - M_post1_q3 - Num8_p1 - Num7_p1 - Num6_p1) > (N_d0 + 1)) && ((i - M_post1_q3 - Num8_p1 - Num7_p1 - Num6_p1) < Node51_p1))// 1dx
			{
				inner_5node(j, i, 1, -(Num6_p1 + N_d0 + 2), (Num5_p1 - N_d0 - 2), Coef, Coef_location);
			}
			else if (((i - M_post1_q3 - Num8_p1 - Num7_p1 - Num6_p1) > (Node51_p1 + 1)) && ((i - M_post1_q3 - Num8_p1 - Num7_p1 - Num6_p1) < (Node_5_p1 - 3)))// 1dx
			{
				inner_5node(j, i, 1, -(Num6_p1 + N_d0 + 2) - 2, (Num5_p1 - N_d0 - 2) + N_c1, Coef, Coef_location);
			}
			else// 2dx!
			{
				inner_5node(j, i, 2, -(Num7_p1 / 2 + 1 + Num6_p1 + Node_5_p1), Num5_p1 - Node_5_p1 + Num4_p1 + Num3_p1 / 2 - 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post1_q3 + Num8_p1 + Num7_p1 + Num6_p1 + Node_5_p1 + 2);
		i<(M_post1_q3 + Num8_p1 + Num7_p1 + Num6_p1 + Num5_p1); i++)// circle num13(5) :  83 node!  part2 !
	{
		switch (i - M_post1_q3 - Num8_p1 - Num7_p1 - Num6_p1)
		{
		case (Num5_p1 - 1) :
			boundary_node(j, i, 2, -Num8_p1 - Num7_p1 - Num6_p1 - Num5_p1, Num4_p1 + Num3_p1 + Num2_p1 + Num1_p1, Coef, Coef_location);
			break;
		case (Node_5_p1 + N_a - 4) ://2 lei shang 1j
			kind2_node(j, i, -(Num9_p1 - Node_9_p1 + Num8_p1 + Num7_p1 + Num6_p1 + Node_5_p1), -2, -1, 0, 1,
			Num5_p1 - Node_5_p1 + Num4_p1 + Num3_p1 + Num2_p1 + Node_1_p1,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_5_p1 + N_a) :// 2 lei 2j shang
			kind2_node(j, i, -(Num7_p1 / 2 + Num6_p1 + Node_5_p1 + N_a - 3), -2, -1, 0, 1, Num5_p1 - Node_5_p1 - N_a + 3 + Num4_p1 + Num3_p1 / 2,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node52_p1:// inner 3 node
			inner_3node(j, i, -1, 0, (Num5_p1 - Node_55_p1 + 4 + Num4_p1 - 1) - N_c1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node52_p1 + 1) :// inner 3 node
			inner_3node(j, i, 0, 1, (Num5_p1 - Node_55_p1 + 4 + Num4_p1 - 1), k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Node_55_p1 - 3) ://2 lei xia 2j
			kind2_node(j, i, -(Num6_p1 + Node_55_p1 + 1), -1, 0, 1, 2, Num5_p1 - Node_55_p1 + Num4_p1 + Num3_p1 - 1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_55_p1 + 1) ://2 lei xia 1j
			kind2_node(j, i, -Num8_p1 - Num7_p1 - Num6_p1 - Num5_p1, -1, 0, 1, 2, Num4_p1 + Num3_p1 + Num2_p1 + Num1_p1,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - M_post1_q3 - Num8_p1 - Num7_p1 - Num6_p1) < (Node_5_p1 + N_a - 4))// 4dx!
			{
				inner_5node(j, i, 4, -(Num9_p1 - Node_9_p1 + Num8_p1 + Num7_p1 + Num6_p1 + Node_5_p1), Num5_p1 - Node_5_p1 + Num4_p1 + Num3_p1 + Num2_p1 + Node_1_p1, Coef, Coef_location);
			}
			else if (((i - M_post1_q3 - Num8_p1 - Num7_p1 - Num6_p1) > (Node_5_p1 + N_a - 4)) && ((i - M_post1_q3 - Num8_p1 - Num7_p1 - Num6_p1) < (Node_5_p1 + N_a)))// 2dx!
			{
				inner_5node(j, i, 2, -(Num7_p1 / 2 + Num6_p1 + Node_5_p1 + N_a - 3), Num5_p1 - Node_5_p1 - N_a + 3 + Num4_p1 + Num3_p1 / 2, Coef, Coef_location);
			}
			else if (((i - M_post1_q3 - Num8_p1 - Num7_p1 - Num6_p1) > (Node_5_p1 + N_a)) && ((i - M_post1_q3 - Num8_p1 - Num7_p1 - Num6_p1) < Node52_p1))// 1dx
			{
				inner_5node(j, i, 1, -(Node_55_p1 - 4 + 1) + 2, (Num5_p1 - Node_55_p1 + 4 + Num4_p1 - 1) - N_c1, Coef, Coef_location);
			}
			else if (((i - M_post1_q3 - Num8_p1 - Num7_p1 - Num6_p1) > (Node52_p1 + 1)) && ((i - M_post1_q3 - Num8_p1 - Num7_p1 - Num6_p1) < (Node_55_p1 - 3)))// 1dx
			{
				inner_5node(j, i, 1, -(Node_55_p1 - 4 + 1), (Num5_p1 - Node_55_p1 + 4 + Num4_p1 - 1), Coef, Coef_location);
			}
			else if (((i - M_post1_q3 - Num8_p1 - Num7_p1 - Num6_p1) > (Node_55_p1 - 3)) && ((i - M_post1_q3 - Num8_p1 - Num7_p1 - Num6_p1) < (Node_55_p1 + 1)))// 2dx
			{
				inner_5node(j, i, 2, -(Num6_p1 + Node_55_p1 + 1), Num5_p1 - Node_55_p1 + Num4_p1 + Num3_p1 - 1, Coef, Coef_location);
			}
			else// 4dx!
			{
				inner_5node(j, i, 4, -Num8_p1 - Num7_p1 - Num6_p1 - Num5_p1, Num4_p1 + Num3_p1 + Num2_p1 + Num1_p1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!! num 34 !!!
	for (i = (M_post1_q3 + Num8_p1 + Num7_p1 + Num6_p1 + Num5_p1);
		i<(M_post1_q3 + Num8_p1 + Num7_p1 + Num6_p1 + Num5_p1 + Num4_p1 / 2); i++)// num 34(4) :  30 node! part1
	{
		switch (i - M_post1_q3 - Num8_p1 - Num7_p1 - Num6_p1 - Num5_p1)
		{
		case 0:// 1 lei 2j shang
			kind1_node(j, i, -(Num5_p1 - N_d0 - 2) - 2, -(Num5_p1 - N_d0 - 2) - 1, -(Num5_p1 - N_d0 - 2), 0, 1, Num4_p1 + 4 - 2,
				Num4_p1 + 4 - 1, Num4_p1 + 4, A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case (Num4_p1 / 2 - 1) :// 1 lei 2j xia
			kind1_node(j, i, -(Num5_p1 - N_d0 - 2) - N_c1, -(Num5_p1 - N_d0 - 2) - N_c1 + 1, -(Num5_p1 - N_d0 - 2) - N_c1 + 2,
			-1, 0, Num4_p1 + 4, Num4_p1 + 4 + 1, Num4_p1 + 4 + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - M_post1_q3 - Num8_p1 - Num7_p1 - Num6_p1 - Num5_p1) <= Node41_p1)// 1dx
			{
				inner_5node(j, i, 1, -(Num5_p1 - N_d0 - 2), Num4_p1 + 4, Coef, Coef_location);
			}
			else if ((i - M_post1_q3 - Num8_p1 - Num7_p1 - Num6_p1 - Num5_p1) >= Node42_p1) // 1dx
			{
				inner_5node(j, i, 1, -(Num5_p1 - N_d0 - 2) - N_c1, Num4_p1 + 4, Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_4node(j, i, -1, 0, 1, Num4_p1 + 4, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post1_q3 + Num8_p1 + Num7_p1 + Num6_p1 + Num5_p1 + Num4_p1 / 2);
		i<(M_post1_q3 + Num8_p1 + Num7_p1 + Num6_p1 + Num5_p1 + Num4_p1); i++)// num 34(4) :  30 node! part2
	{
		switch (i - M_post1_q3 - Num8_p1 - Num7_p1 - Num6_p1 - Num5_p1)
		{
		case (Num4_p1 / 2) :// 1 lei 2j shang
			kind1_node(j, i, -(Num5_p1 - Node_55_p1 + 4 + Num4_p1 - 1) + N_c1 - 2, -(Num5_p1 - Node_55_p1 + 4 + Num4_p1 - 1) + N_c1 - 1,
			-(Num5_p1 - Node_55_p1 + 4 + Num4_p1 - 1) + N_c1, 0, 1, (Num3_p1 - 5 + 1) - 2, (Num3_p1 - 5 + 1) - 1, (Num3_p1 - 5 + 1),
			A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case (Num4_p1 - 1) :// 1 lei 2j xia
			kind1_node(j, i, -(Num5_p1 - Node_55_p1 + 4 + Num4_p1 - 1), -(Num5_p1 - Node_55_p1 + 4 + Num4_p1 - 1) + 1,
			-(Num5_p1 - Node_55_p1 + 4 + Num4_p1 - 1) + 2, -1, 0, (Num3_p1 - 5 + 1), (Num3_p1 - 5 + 1) + 1, (Num3_p1 - 5 + 1) + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - M_post1_q3 - Num8_p1 - Num7_p1 - Num6_p1 - Num5_p1) <= Node43_p1)// 1dx
			{
				inner_5node(j, i, 1, -(Num5_p1 - Node_55_p1 + 4 + Num4_p1 - 1) + N_c1, (Num3_p1 - 5 + 1), Coef, Coef_location);
			}
			else if ((i - M_post1_q3 - Num8_p1 - Num7_p1 - Num6_p1 - Num5_p1) >= Node44_p1) // 1dx
			{
				inner_5node(j, i, 1, -(Num5_p1 - Node_55_p1 + 4 + Num4_p1 - 1), (Num3_p1 - 5 + 1), Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_4node(j, i, -1, 0, 1, (Num3_p1 - 5 + 1), 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!! num 35 !!!
	for (i = (M_post1_q3 + Num8_p1 + Num7_p1 + Num6_p1 + Num5_p1 + Num4_p1);
		i<(M_post1_q3 + Num8_p1 + Num7_p1 + Num6_p1 + Num5_p1 + Num4_p1 + Num3_p1 / 2); i++)// num 35(3) :  46 node! part1
	{
		switch (i - M_post1_q3 - Num8_p1 - Num7_p1 - Num6_p1 - Num5_p1 - Num4_p1)
		{
		case 0:// 1 lei 1j shang
			kind1_node(j, i, -(Num5_p1 - N_d0 + 2 + Num4_p1) - 2, -(Num5_p1 - N_d0 + 2 + Num4_p1) - 1, -(Num5_p1 - N_d0 + 2 + Num4_p1),
				0, 1, Num3_p1 + Num2_p1 + N_d0 - 2 - 2, Num3_p1 + Num2_p1 + N_d0 - 2 - 1, Num3_p1 + Num2_p1 + N_d0 - 2,
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case 3:// 2 lei 2j shang
			kind2_node(j, i, -(Num5_p1 - N_d0 + 2 + Num4_p1), -2, -1, 0, 1, Num3_p1 + Num2_p1 + N_d0 - 2,
				B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case 4:// 3 lei 2j shang
			kind3_node(j, i, -(Num4_p1 + 4), -1, 0, 1, (Num3_p1 + Num2_p1 + N_d0 - 2) - 1, (Num3_p1 + Num2_p1 + N_d0 - 2),
				C2_2, C3_2, C1_2, C3_2, C4, C4, Coef, Coef_location);
			break;
		case (Num3_p1 / 2 - 5) :// 3 lei 2j xia
			kind3_node(j, i, -(Num4_p1 + 4), -1, 0, 1, Num3_p1 / 2 + 1 + Num2_p1 + Node_1_p1, Num3_p1 / 2 + 1 + Num2_p1 + Node_1_p1 + 1,
			C2_2, C3_2, C1_2, C3_2, C4, C4, Coef, Coef_location);
			break;
		case (Num3_p1 / 2 - 4) :// 2 lei 2j xia
			kind2_node(j, i, -(Num5_p1 - Node_5_p1 + Num4_p1 + Num3_p1 / 2 - 1), -1, 0, 1, 2, Num3_p1 / 2 + 1 + Num2_p1 + Node_1_p1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num3_p1 / 2 - 1) :// 1 lei 1j xia
			kind1_node(j, i, -(Num5_p1 - Node_5_p1 + Num4_p1 + Num3_p1 / 2 - 1), -(Num5_p1 - Node_5_p1 + Num4_p1 + Num3_p1 / 2 - 1) + 1,
			-(Num5_p1 - Node_5_p1 + Num4_p1 + Num3_p1 / 2 - 1) + 2, -1, 0,
			Num3_p1 / 2 + 1 + Num2_p1 + Node_1_p1, Num3_p1 / 2 + 1 + Num2_p1 + Node_1_p1 + 1, Num3_p1 / 2 + 1 + Num2_p1 + Node_1_p1 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - M_post1_q3 - Num8_p1 - Num7_p1 - Num6_p1 - Num5_p1 - Num4_p1) < 3)// 2dx
			{
				inner_5node(j, i, 2, -(Num5_p1 - N_d0 + 2 + Num4_p1), Num3_p1 + Num2_p1 + N_d0 - 2, Coef, Coef_location);
			}
			else if (((i - M_post1_q3 - Num8_p1 - Num7_p1 - Num6_p1 - Num5_p1 - Num4_p1) > 4) && ((i - M_post1_q3 - Num8_p1 - Num7_p1 - Num6_p1 - Num5_p1 - Num4_p1) < (Num3_p1 / 2 - 5)))// 1dx
			{
				inner_5node(j, i, 1, -(Num4_p1 + 4), Num3_p1 - 5, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, i, 2, -(Num5_p1 - Node_5_p1 + Num4_p1 + Num3_p1 / 2 - 1), Num3_p1 / 2 + 1 + Num2_p1 + Node_1_p1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post1_q3 + Num8_p1 + Num7_p1 + Num6_p1 + Num5_p1 + Num4_p1 + Num3_p1 / 2);
		i<(M_post1_q3 + Num8_p1 + Num7_p1 + Num6_p1 + Num5_p1 + Num4_p1 + Num3_p1); i++)// num 35(3) :  46 node! part2
	{
		switch (i - M_post1_q3 - Num8_p1 - Num7_p1 - Num6_p1 - Num5_p1 - Num4_p1 - Num3_p1 / 2)
		{
		case 0:// 1 lei 1j shang
			kind1_node(j, i, -(Num5_p1 - Node_5_p1 - N_a + 3 + Num4_p1 + Num3_p1 / 2) - 2, -(Num5_p1 - Node_5_p1 - N_a + 3 + Num4_p1 + Num3_p1 / 2) - 1,
				-(Num5_p1 - Node_5_p1 - N_a + 3 + Num4_p1 + Num3_p1 / 2), 0, 1,
				Num3_p1 / 2 + Num2_p1 + Node_1_p1 + N_a - 3 - 2, Num3_p1 / 2 + Num2_p1 + Node_1_p1 + N_a - 3 - 1, Num3_p1 / 2 + Num2_p1 + Node_1_p1 + N_a - 3,
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case 3:// 2 lei 2j shang
			kind2_node(j, i, -(Num5_p1 - Node_5_p1 - N_a + 3 + Num4_p1 + Num3_p1 / 2), -2, -1, 0, 1, Num3_p1 / 2 + Num2_p1 + Node_1_p1 + N_a - 3,
				B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case 4:// 3 lei 2j shang
			kind3_node(j, i, -(Num3_p1 - 5 + 1), -1, 0, 1, (Num3_p1 / 2 + Num2_p1 + Node_1_p1 + N_a - 3) - 1, (Num3_p1 / 2 + Num2_p1 + Node_1_p1 + N_a - 3),
				C2_2, C3_2, C1_2, C3_2, C4, C4, Coef, Coef_location);
			break;
		case (Num3_p1 / 2 - 5) :// 3 lei 2j xia
			kind3_node(j, i, -(Num3_p1 - 5 + 1), -1, 0, 1, 1 + Num2_p1 + Node_11_p1, 1 + Num2_p1 + Node_11_p1 + 1,
			C2_2, C3_2, C1_2, C3_2, C4, C4, Coef, Coef_location);
			break;
		case (Num3_p1 / 2 - 4) :// 2 lei 2j xia
			kind2_node(j, i, -(Num5_p1 - Node_55_p1 + Num4_p1 + Num3_p1 - 1), -1, 0, 1, 2, 1 + Num2_p1 + Node_11_p1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num3_p1 / 2 - 1) :// 1 lei 1j xia
			kind1_node(j, i, -(Num5_p1 - Node_55_p1 + Num4_p1 + Num3_p1 - 1), -(Num5_p1 - Node_55_p1 + Num4_p1 + Num3_p1 - 1) + 1,
			-(Num5_p1 - Node_55_p1 + Num4_p1 + Num3_p1 - 1) + 2, -1, 0, 1 + Num2_p1 + Node_11_p1, 1 + Num2_p1 + Node_11_p1 + 1, 1 + Num2_p1 + Node_11_p1 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - M_post1_q3 - Num8_p1 - Num7_p1 - Num6_p1 - Num5_p1 - Num4_p1 - Num3_p1 / 2) < 3)// 2dx
			{
				inner_5node(j, i, 2, -(Num5_p1 - Node_5_p1 - N_a + 3 + Num4_p1 + Num3_p1 / 2), Num3_p1 / 2 + Num2_p1 + Node_1_p1 + N_a - 3, Coef, Coef_location);
			}
			else if (((i - M_post1_q3 - Num8_p1 - Num7_p1 - Num6_p1 - Num5_p1 - Num4_p1 - Num3_p1 / 2) > 4) && ((i - M_post1_q3 - Num8_p1 - Num7_p1 - Num6_p1 - Num5_p1 - Num4_p1 - Num3_p1 / 2) < (Num3_p1 / 2 - 5)))// 1dx
			{
				inner_5node(j, i, 1, -(Num3_p1 - 5 + 1), 6 + Num2_p1 - 1, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, i, 2, -(Num5_p1 - Node_55_p1 + Num4_p1 + Num3_p1 - 1), 1 + Num2_p1 + Node_11_p1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!! num 36 !!!
	for (i = (M_post1_q3 + Num8_p1 + Num7_p1 + Num6_p1 + Num5_p1 + Num4_p1 + Num3_p1);
		i<(M_post1_q3 + Num8_p1 + Num7_p1 + Num6_p1 + Num5_p1 + Num4_p1 + Num3_p1 + Num2_p1 / 2); i++)// num 36(2) :  26 node! part1
	{
		switch (i - M_post1_q3 - Num8_p1 - Num7_p1 - Num6_p1 - Num5_p1 - Num4_p1 - Num3_p1)
		{
		case 0:// 3 lei 2j shang
			kind3_node(j, i, -(Num3_p1 - 3), -(Num3_p1 - 3) + 2, 0, 1, Num2_p1 + Node11_p1 - 1, Num2_p1 + Node11_p1,
				C4, C3_2, C1_2, C2_2, C4, C3_2, Coef, Coef_location);
			break;
		case (Num2_p1 / 2 - 1) ://node2: 3 lei 2j xia
			kind3_node(j, i, -(Num3_p1 - 5), -(Num3_p1 - 5) + 2, -1, 0, Num2_p1 / 2 + Node_1_p1 - 4 + 1, Num2_p1 / 2 + Node_1_p1 - 4 + 1 + 1,
			C3_2, C4, C2_2, C1_2, C3_2, C4, Coef, Coef_location);
			break;
		default:
			if ((i - M_post1_q3 - Num8_p1 - Num7_p1 - Num6_p1 - Num5_p1 - Num4_p1 - Num3_p1) % 2 == 0)// 1dx
			{
				inner_5node(j, i, 1, -(Num3_p1 - 5),
					(Num2_p1 - (i - M_post1_q3 - Num8_p1 - Num7_p1 - Num6_p1 - Num5_p1 - Num4_p1 - Num3_p1) + Node11_p1 + (i - M_post1_q3 - Num8_p1 - Num7_p1 - Num6_p1 - Num5_p1 - Num4_p1 - Num3_p1) / 2), Coef, Coef_location);
			}
			else//1 lei  2j left
			{
				kind1_node(j, i, -(Num3_p1 - 5), -1, 0, 1,
					(Num2_p1 - (i - M_post1_q3 - Num8_p1 - Num7_p1 - Num6_p1 - Num5_p1 - Num4_p1 - Num3_p1) + Node11_p1 + (i - M_post1_q3 - Num8_p1 - Num7_p1 - Num6_p1 - Num5_p1 - Num4_p1 - Num3_p1) / 2),
					(Num2_p1 - (i - M_post1_q3 - Num8_p1 - Num7_p1 - Num6_p1 - Num5_p1 - Num4_p1 - Num3_p1) + Node11_p1 + (i - M_post1_q3 - Num8_p1 - Num7_p1 - Num6_p1 - Num5_p1 - Num4_p1 - Num3_p1) / 2) + 1,
					(Num2_p1 - (i - M_post1_q3 - Num8_p1 - Num7_p1 - Num6_p1 - Num5_p1 - Num4_p1 - Num3_p1) + Node11_p1 + (i - M_post1_q3 - Num8_p1 - Num7_p1 - Num6_p1 - Num5_p1 - Num4_p1 - Num3_p1) / 2) + (Num1_p1 - N_d0 + 2),
					(Num2_p1 - (i - M_post1_q3 - Num8_p1 - Num7_p1 - Num6_p1 - Num5_p1 - Num4_p1 - Num3_p1) + Node11_p1 + (i - M_post1_q3 - Num8_p1 - Num7_p1 - Num6_p1 - Num5_p1 - Num4_p1 - Num3_p1) / 2) + 1 + (Num1_p1 - N_d0 + 2),
					A4_2, A3_2, A5_2, A3_2, A1_2, A1_2, A2_2, A2_2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post1_q3 + Num8_p1 + Num7_p1 + Num6_p1 + Num5_p1 + Num4_p1 + Num3_p1 + Num2_p1 / 2);
		i<(M_post1_q3 + Num8_p1 + Num7_p1 + Num6_p1 + Num5_p1 + Num4_p1 + Num3_p1 + Num2_p1); i++)// num 36(2) :  26 node! part2
	{
		switch (i - M_post1_q3 - Num8_p1 - Num7_p1 - Num6_p1 - Num5_p1 - Num4_p1 - Num3_p1 - Num2_p1 / 2)
		{
		case 0:// 3 lei 2j shang
			kind3_node(j, i, -(6 + Num2_p1 - 1) - 2, -(6 + Num2_p1 - 1), 0, 1, Num2_p1 / 2 + Node13_p1 - 1, Num2_p1 / 2 + Node13_p1,
				C4, C3_2, C1_2, C2_2, C4, C3_2, Coef, Coef_location);
			break;
		case (Num2_p1 / 2 - 1) :// 3 lei 2j xia
			kind3_node(j, i, -(6 + Num2_p1 - 1), -(6 + Num2_p1 - 1) + 2, -1, 0, Node_11_p1 - 4 + 1, Node_11_p1 - 4 + 1 + 1,
			C3_2, C4, C2_2, C1_2, C3_2, C4, Coef, Coef_location);
			break;
		default:
			if ((i - M_post1_q3 - Num8_p1 - Num7_p1 - Num6_p1 - Num5_p1 - Num4_p1 - Num3_p1 - Num2_p1 / 2) % 2 == 0)// 1dx
			{
				inner_5node(j, i, 1, -(6 + Num2_p1 - 1), Num2_p1 / 2 - (i - M_post1_q3 - Num8_p1 - Num7_p1 - Num6_p1 - Num5_p1 - Num4_p1 - Num3_p1 - Num2_p1 / 2) + Node13_p1 + (i - M_post1_q3 - Num8_p1 - Num7_p1 - Num6_p1 - Num5_p1 - Num4_p1 - Num3_p1 - Num2_p1 / 2) / 2,
					Coef, Coef_location);
			}
			else//1 lei 2j left
			{
				kind1_node(j, i, -(6 + Num2_p1 - 1), -1, 0, 1,
					(Num2_p1 / 2 - (i - M_post1_q3 - Num8_p1 - Num7_p1 - Num6_p1 - Num5_p1 - Num4_p1 - Num3_p1 - Num2_p1 / 2) + Node13_p1 + (i - M_post1_q3 - Num8_p1 - Num7_p1 - Num6_p1 - Num5_p1 - Num4_p1 - Num3_p1 - Num2_p1 / 2) / 2),
					(Num2_p1 / 2 - (i - M_post1_q3 - Num8_p1 - Num7_p1 - Num6_p1 - Num5_p1 - Num4_p1 - Num3_p1 - Num2_p1 / 2) + Node13_p1 + (i - M_post1_q3 - Num8_p1 - Num7_p1 - Num6_p1 - Num5_p1 - Num4_p1 - Num3_p1 - Num2_p1 / 2) / 2) + 1,
					(Num2_p1 / 2 - (i - M_post1_q3 - Num8_p1 - Num7_p1 - Num6_p1 - Num5_p1 - Num4_p1 - Num3_p1 - Num2_p1 / 2) + Node13_p1 + (i - M_post1_q3 - Num8_p1 - Num7_p1 - Num6_p1 - Num5_p1 - Num4_p1 - Num3_p1 - Num2_p1 / 2) / 2) + (Num1_p1 - Node_11_p1 + Nx_3 - 1),
					(Num2_p1 / 2 - (i - M_post1_q3 - Num8_p1 - Num7_p1 - Num6_p1 - Num5_p1 - Num4_p1 - Num3_p1 - Num2_p1 / 2) + Node13_p1 + (i - M_post1_q3 - Num8_p1 - Num7_p1 - Num6_p1 - Num5_p1 - Num4_p1 - Num3_p1 - Num2_p1 / 2) / 2) + 1 + (Num1_p1 - Node_11_p1 + Nx_3 - 1),
					A4_2, A3_2, A5_2, A3_2, A1_2, A1_2, A2_2, A2_2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!! num 37 !!!
	for (i = (num_total - Nx_2 - Nx_3 - Num1_p1); i<(num_total - Nx_2 - Nx_3 - Num1_p1 + Node_1_p1 + 2); i++)// circle num 37(1) Nx_2:  74 node!  part1 !
	{
		switch (i - num_total + Nx_2 + Nx_3 + Num1_p1)
		{
		case 0:
			boundary_node(j, i, 1, -Num5_p1 - Num4_p1 - Num3_p1 - Num2_p1, Nx_2 + Nx_3, Coef, Coef_location);
			break;
		case (N_d0 - 3) ://2 lei shang 1j
			kind2_node(j, i, -Num5_p1 - Num4_p1 - Num3_p1 - Num2_p1, -2, -1, 0, 1, Nx_2 + Nx_3,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_1_p1 + 1) ://2 lei xia 1j
			kind2_node(j, i, -(Num5_p1 - Node_5_p1 + Num4_p1 + Num3_p1 + Num2_p1 + Node_1_p1), -1, 0, 1, 2, Nx_2 + Nx_3,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - num_total + Nx_2 + Nx_3 + Num1_p1) < (N_d0 - 3))// 4dx!
			{
				inner_5node(j, i, 4, -Num5_p1 - Num4_p1 - Num3_p1 - Num2_p1, Nx_2 + Nx_3, Coef, Coef_location);
			}
			else if (((i - num_total + Nx_2 + Nx_3 + Num1_p1) > (N_d0 - 3)) && ((i - num_total + Nx_2 + Nx_3 + Num1_p1) < Node11_p1))// 2dx!
			{
				inner_5node(j, i, 2, -(Num3_p1 + Num2_p1 + N_d0 - 2), Num1_p1 - N_d0 + 2, Coef, Coef_location);
			}
			else if (((i - num_total + Nx_2 + Nx_3 + Num1_p1) >= Node11_p1) && ((i - num_total + Nx_2 + Nx_3 + Num1_p1) <= Node12_p1))// 2 lei right 2j !
			{
				kind2_node(j, i, -((i - num_total + Nx_2 + Nx_3 + Num1_p1) + Num2_p1 - (i - num_total + Nx_2 + Nx_3 + Num1_p1 - Node11_p1) * 2),
					-1, 0, 1, Num1_p1 - N_d0 + 2, Nx_2 + Nx_3, B2, B3, B1_2, B3, B4, B5, Coef, Coef_location);
			}
			else// 2dx!
			{
				inner_5node(j, i, 2, -(Num3_p1 / 2 + 1 + Num2_p1 + Node_1_p1), Num1_p1 - N_d0 + 2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (num_total - Nx_2 - Nx_3 - Num1_p1 + Node_1_p1 + 2); i<(num_total - Nx_2 - Nx_3); i++)// circle num 37(1) Nx_2:  74 node!  part2 !
	{
		switch (i - num_total + Nx_2 + Nx_3 + Num1_p1)
		{
		case (Num1_p1 - 1) :
			boundary_node(j, i, 2, -Num1_p1 - Num2_p1 - Num3_p1 - Num4_p1, Nx_2 + Nx_3, Coef, Coef_location);
			break;
		case (Node_1_p1 + N_a - 4) ://2 lei shang 1j
			kind2_node(j, i, -(Num5_p1 - Node_5_p1 + Num4_p1 + Num3_p1 + Num2_p1 + Node_1_p1), -2, -1, 0, 1, Nx_2 + Nx_3,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_11_p1 + 1) ://2 lei xia 1j
			kind2_node(j, i, -Num1_p1 - Num2_p1 - Num3_p1 - Num4_p1, -1, 0, 1, 2, Nx_2 + Nx_3,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - num_total + Nx_2 + Nx_3 + Num1_p1) < (Node_1_p1 + N_a - 4))// 4dx!
			{
				inner_5node(j, i, 4, -(Num5_p1 - Node_5_p1 + Num4_p1 + Num3_p1 + Num2_p1 + Node_1_p1), Nx_2 + Nx_3, Coef, Coef_location);
			}
			else if (((i - num_total + Nx_2 + Nx_3 + Num1_p1) > (Node_1_p1 + N_a - 4)) && ((i - num_total + Nx_2 + Nx_3 + Num1_p1) < Node13_p1))// 2dx!
			{
				inner_5node(j, i, 2, -(Num3_p1 / 2 + Num2_p1 + Node_1_p1 + N_a - 3), Num1_p1 - Node_11_p1 + Nx_3 - 1, Coef, Coef_location);
			}
			else if (((i - num_total + Nx_2 + Nx_3 + Num1_p1) >= Node13_p1) && ((i - num_total + Nx_2 + Nx_3 + Num1_p1) <= Node14_p1))// 2 lei left 2j !
			{
				kind2_node(j, i, -((i - num_total + Nx_2 + Nx_3 + Num1_p1) + Num2_p1 / 2 - (i - num_total + Nx_2 + Nx_3 + Num1_p1 - Node13_p1) * 2),
					-1, 0, 1, Num1_p1 - Node_11_p1 + Nx_3 - 1, Nx_2 + Nx_3, B2, B3, B1_2, B3, B4, B5, Coef, Coef_location);
			}
			else if (((i - num_total + Nx_2 + Nx_3 + Num1_p1) >= Node14_p1) && ((i - num_total + Nx_2 + Nx_3 + Num1_p1) <= Node_11_p1))// 2dx
			{
				inner_5node(j, i, 2, -(Num2_p1 + Node_11_p1 + 1), Num1_p1 - Node_11_p1 + Nx_3 - 1, Coef, Coef_location);
			}
			else
			{
				inner_5node(j, i, 4, -Num1_p1 - Num2_p1 - Num3_p1 - Num4_p1, Nx_2 + Nx_3, Coef, Coef_location);// 4dx!
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (num_total - Nx_2 - Nx_3); i<(num_total - Nx_2); i++)//right1 Nx_3:  30 node! 
	{
		switch (i - num_total + Nx_2 + Nx_3)
		{
		case 0:// 1 lei shang
			kind1_node(j, i, -(Nx_2 - N_d0 + 2) - 2, -(Nx_2 - N_d0 + 2) - 1, -(Nx_2 - N_d0 + 2), 0, 1, Nx_3 + N_d0 - 2 - 2, Nx_3 + N_d0 - 2 - 1, Nx_3 + N_d0 - 2,
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case (Nx_3 / 2 - 1) :// 1 lei xia
			kind1_node(j, i, -(Nx_2 - N_d0 + 2), -(Nx_2 - N_d0 + 2) + 1, -(Nx_2 - N_d0 + 2) + 2, -1, 0,
			Nx_3 + N_d0 - 2, Nx_3 + N_d0 - 2 + 1, Nx_3 + N_d0 - 2 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		case (Nx_3 / 2) :// 1 lei shang
			kind1_node(j, i, -(Nx_2 - Node_11_p1 - 1 + Nx_3) - 2, -(Nx_2 - Node_11_p1 - 1 + Nx_3) - 1, -(Nx_2 - Node_11_p1 - 1 + Nx_3), 0, 1,
			1 + Node_11_p1 - 2, 1 + Node_11_p1 - 1, 1 + Node_11_p1,
			A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case (Nx_3 - 1) :// 1 lei xia
			kind1_node(j, i, -(Nx_2 - Node_11_p1 - 1 + Nx_3), -(Nx_2 - Node_11_p1 - 1 + Nx_3) + 1, -(Nx_2 - Node_11_p1 - 1 + Nx_3) + 2, -1, 0,
			1 + Node_11_p1, 1 + Node_11_p1 + 1, 1 + Node_11_p1 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - num_total + Nx_2 + Nx_3) < (Nx_3 / 2 - 1))
			{
				inner_5node(j, i, 2, -(Nx_2 - N_d0 + 2), Nx_3 + N_d0 - 2, Coef, Coef_location);
			}
			else
			{
				inner_5node(j, i, 2, -(Nx_2 - Node_11_p1 - 1 + Nx_3), 1 + Node_11_p1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	//***************************************************************************************
	//最后一行 本节点向下的区域没有影响
	//cout << "circle num1:" << j << endl;
	for (i = (Nx_2 + Nx_3); i<(Nx_2 + Nx_3 + Node_1 + 2); i++)// circle num1 Nx_2:  74 node!  part1 !
	{
		switch (i - Nx_2 - Nx_3)
		{
		case 0:
			boundary_node(j, j - offset, 1, -Nx_2 - Nx_3, Num1 + Num2 + Num3 + Num4, Coef, Coef_location);
			break;
		case (N_d0 - 3) ://2 lei shang 1j
			kind2_node(j, j - offset, -Nx_2 - Nx_3, -2, -1, 0, 1, Num1 + Num2 + Num3 + Num4,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_1 + 1) ://2 lei xia 1j
			kind2_node(j, j - offset, -(Nx_2 + Nx_3), -1, 0, 1, 2, Num1 - Node_1 + Num2 + Num3 + Num4 + Node_5,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3) < (N_d0 - 3))// 4dx!
			{
				inner_5node(j, j - offset, 4, -Nx_3 - Nx_2, Num1 + Num2 + Num3 + Num4, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3) > (N_d0 - 3)) && ((i - Nx_2 - Nx_3) < Node11))// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Nx_3 + N_d0 - 2), Num1 - N_d0 + 2 + Num2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3) >= Node11) && ((i - Nx_2 - Nx_3) <= Node12))// 2 lei left 2j !
			{
				kind2_node(j, j - offset, -Nx_3 - Nx_2, -(Nx_3 + N_d0 - 2), -1, 0, 1, Num1 - (i - Nx_2 - Nx_3) + ((i - Nx_2 - Nx_3) - Node11) * 2,
					B5, B4, B3, B1_2, B3, B2, Coef, Coef_location);
			}
			else// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Nx_3 + N_d0 - 2), Num1 - Node_1 + Num2 + Num3 / 2 - 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Node_1 + 2); i<(Nx_2 + Nx_3 + Num1); i++)// circle num1 Nx_2:  74 node!  part2 !
	{
		switch (i - Nx_2 - Nx_3)
		{
		case (Num1 - 1) :
			boundary_node(j, j - offset, 2, -Nx_2 - Nx_3, Num5 + Num2 + Num3 + Num4, Coef, Coef_location);
			break;
		case (Node_1 + N_a - 4) ://2 lei shang 1j
			kind2_node(j, j - offset, -Nx_2 - Nx_3, -2, -1, 0, 1, Num1 - Node_1 + Num2 + Num3 + Num4 + Node_5,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_11 + 1) ://2 lei xia 1j
			kind2_node(j, j - offset, -(Nx_2 + Nx_3), -1, 0, 1, 2, Num5 + Num2 + Num3 + Num4,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3) < (Node_1 + N_a - 4))// 4dx!
			{
				inner_5node(j, j - offset, 4, -Nx_3 - Nx_2, Num1 - Node_1 + Num2 + Num3 + Num4 + Node_5, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3) > (Node_1 + N_a - 4)) && ((i - Nx_2 - Nx_3) < Node13))// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Node_11 + 1), Num1 - Node_1 - N_a + 3 + Num2 + Num3 / 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3) >= Node13) && ((i - Nx_2 - Nx_3) <= Node14))// 2 lei left 2j !
			{
				kind2_node(j, j - offset, -Nx_3 - Nx_2, -(Node_11 + 1), -1, 0, 1, Num1 - (i - Nx_2 - Nx_3) + Num2 / 2 + ((i - Nx_2 - Nx_3) - Node13) * 2,
					B5, B4, B3, B1_2, B3, B2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3) >= Node14) && ((i - Nx_2 - Nx_3) <= Node_11))// 2dx
			{
				inner_5node(j, j - offset, 2, -(Node_11 + 1), Num1 - Node_11 + Num2 + Num3 - 1, Coef, Coef_location);
			}
			else
			{
				inner_5node(j, j - offset, 4, -(Nx_3 + Nx_2), Num5 + Num2 + Num3 + Num4, Coef, Coef_location);// 4dx!
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!!!!!!!!!
	//cout << "circle num 2:" << j << endl;
	for (i = (Nx_2 + Nx_3 + Num1); i<(Nx_2 + Nx_3 + Num1 + Num2 / 2); i++)// num 2 :  26 node! part1
	{
		switch (i - Nx_2 - Nx_3 - Num1)
		{
		case 0:// 3 lei 2j shang
			kind3_node(j, j - offset, -(Num1 - Node11) - 1, -(Num1 - Node11), 0, 1, Num2 + 5 - 2, Num2 + 5,
				C4, C3_2, C1_2, C2_2, C4, C3_2, Coef, Coef_location);
			break;
		case (Num2 / 2 - 1) ://node2: 3 lei 2j xia
			kind3_node(j, j - offset, -(Num1 - Node12 + Num2 / 2 - 1), -(Num1 - Node12 + Num2 / 2 - 1) + 1, -1, 0, Num2 + 5, Num2 + 5 + 2,
			C3_2, C4, C2_2, C1_2, C3_2, C4, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1) % 2 == 0)// 1dx
			{
				inner_5node(j, j - offset, 1, -((i - Nx_2 - Nx_3 - Num1) + Num1 - Node11 - (i - Nx_2 - Nx_3 - Num1) / 2), Num2 + 5, Coef, Coef_location);
			}
			else//1 lei  2j left
			{
				kind1_node(j, j - offset, -((i - Nx_2 - Nx_3 - Num1) + Num1 - Node11 - (i - Nx_2 - Nx_3 - Num1) / 2) - (Nx_3 + N_d0 - 2),
					-((i - Nx_2 - Nx_3 - Num1) + Num1 - Node11 - (i - Nx_2 - Nx_3 - Num1) / 2) - (Nx_3 + N_d0 - 2) + 1,
					-((i - Nx_2 - Nx_3 - Num1) + Num1 - Node11 - (i - Nx_2 - Nx_3 - Num1) / 2),
					-((i - Nx_2 - Nx_3 - Num1) + Num1 - Node11 - (i - Nx_2 - Nx_3 - Num1) / 2) + 1,
					-1, 0, 1, 5 + Num2, A2_2, A2_2, A1_2, A1_2, A3_2, A5_2, A3_2, A4_2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 / 2); i<(Nx_2 + Nx_3 + Num1 + Num2); i++)// num 2 :  26 node! part2
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 / 2)
		{
		case 0:// 3 lei 2j shang
			kind3_node(j, j - offset, -(Num1 - Node13 + Num2 / 2) - 1, -(Num1 - Node13 + Num2 / 2), 0, 1, Num3 - 5 - 2, Num3 - 5,
				C4, C3_2, C1_2, C2_2, C4, C3_2, Coef, Coef_location);
			break;
		case (Num2 / 2 - 1) :// 3 lei 2j xia
			kind3_node(j, j - offset, -(Num1 - Node14 + Num2 - 1), -(Num1 - Node14 + Num2 - 1) + 1, -1, 0, Num3 - 5, Num3 - 5 + 2,
			C3_2, C4, C2_2, C1_2, C3_2, C4, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 / 2) % 2 == 0)// 1dx
			{
				inner_5node(j, j - offset, 1, -((i - Nx_2 - Nx_3 - Num1 - Num2 / 2) + Num1 - Node13 - (i - Nx_2 - Nx_3 - Num1 - Num2 / 2) / 2) - Num2 / 2,
					Num3 - 5, Coef, Coef_location);
			}
			else//1 lei 2j left
			{
				kind1_node(j, j - offset, -((i - Nx_2 - Nx_3 - Num1 - Num2 / 2) + Num1 - Node13 - (i - Nx_2 - Nx_3 - Num1 - Num2 / 2) / 2) - Num2 / 2 - (Node_11 + 1),
					-((i - Nx_2 - Nx_3 - Num1 - Num2 / 2) + Num1 - Node13 - (i - Nx_2 - Nx_3 - Num1 - Num2 / 2) / 2) - Num2 / 2 - (Node_11 + 1) + 1,
					-((i - Nx_2 - Nx_3 - Num1 - Num2 / 2) + Num1 - Node13 - (i - Nx_2 - Nx_3 - Num1 - Num2 / 2) / 2) - Num2 / 2,
					-((i - Nx_2 - Nx_3 - Num1 - Num2 / 2) + Num1 - Node13 - (i - Nx_2 - Nx_3 - Num1 - Num2 / 2) / 2) - Num2 / 2 + 1,
					-1, 0, 1, Num3 - 5, A2_2, A2_2, A1_2, A1_2, A3_2, A5_2, A3_2, A4_2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!!!
	//cout << "circle num 3:" << j << endl;
	for (i = (Nx_2 + Nx_3 + Num1 + Num2); i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 / 2); i++)// num 3 :  46 node! part1
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2)
		{
		case 0:// 1 lei 1j shang
			kind1_node(j, j - offset, -(Num1 - N_d0 + 2 + Num2) - 2, -(Num1 - N_d0 + 2 + Num2) - 1, -(Num1 - N_d0 + 2 + Num2), 0, 1, Num3 + Num4 + N_d0 - 2 - 2,
				Num3 + Num4 + N_d0 - 2 - 1, Num3 + Num4 + N_d0 - 2, A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case 3:// 2 lei 2j shang
			kind2_node(j, j - offset, -(Num1 - N_d0 + 2 + Num2), -2, -1, 0, 1, Num3 + Num4 + N_d0 - 2,
				B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case 4:// 3 lei 2j shang
			kind3_node(j, j - offset, -(Num1 - N_d0 + 2 + Num2) - 1, -(Num1 - N_d0 + 2 + Num2), -1, 0, 1, Num3 - 4,
				C4, C4, C3_2, C1_2, C3_2, C2_2, Coef, Coef_location);
			break;
		case (Num3 / 2 - 5) :// 3 lei 2j xia
			kind3_node(j, j - offset, -(Num1 - Node_1 + Num2 + Num3 / 2 - 1), -(Num1 - Node_1 + Num2 + Num3 / 2 - 1) + 1, -1, 0, 1, Num3 - 4,
			C4, C4, C3_2, C1_2, C3_2, C2_2, Coef, Coef_location);
			break;
		case (Num3 / 2 - 4) :// 2 lei 2j xia
			kind2_node(j, j - offset, -(Num1 - Node_1 + Num2 + Num3 / 2 - 1), -1, 0, 1, 2, Num3 / 2 + 1 + Num4 + Node_5,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num3 / 2 - 1) :// 1 lei 1j xia
			kind1_node(j, j - offset, -(Num1 - Node_1 + Num2 + Num3 / 2 - 1), -(Num1 - Node_1 + Num2 + Num3 / 2 - 1) + 1, -(Num1 - Node_1 + Num2 + Num3 / 2 - 1) + 2, -1, 0,
			Num3 / 2 + 1 + Num4 + Node_5, Num3 / 2 + 1 + Num4 + Node_5 + 1, Num3 / 2 + 1 + Num4 + Node_5 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2) < 3)// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num1 - N_d0 + 2 + Num2), Num3 + Num4 + N_d0 - 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2) > 4) && ((i - Nx_2 - Nx_3 - Num1 - Num2) < (Num3 / 2 - 5)))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num2 + 5), Num3 - 4, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num1 - Node_1 + Num2 + Num3 / 2 - 1), Num3 / 2 + 1 + Num4 + Node_5, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 / 2); i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3); i++)// num 3 :  46 node! part2
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 / 2)
		{
		case 0:// 1 lei 1j shang
			kind1_node(j, j - offset, -(Num1 - Node_1 - N_a + 3 + Num2 + Num3 / 2) - 2, -(Num1 - Node_1 - N_a + 3 + Num2 + Num3 / 2) - 1, -(Num1 - Node_1 - N_a + 3 + Num2 + Num3 / 2), 0, 1,
				Num3 / 2 + Num4 + Node_5 + N_a - 3 - 2, Num3 / 2 + Num4 + Node_5 + N_a - 3 - 1, Num3 / 2 + Num4 + Node_5 + N_a - 3,
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case 3:// 2 lei 2j shang
			kind2_node(j, j - offset, -(Num1 - Node_1 - N_a + 3 + Num2 + Num3 / 2), -2, -1, 0, 1, Num3 / 2 + Num4 + Node_5 + N_a - 3,
				B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case 4:// 3 lei 2j shang
			kind3_node(j, j - offset, -(Num1 - Node_1 - N_a + 3 + Num2 + Num3 / 2) - 1, -(Num1 - Node_1 - N_a + 3 + Num2 + Num3 / 2), -1, 0, 1, 5 + Num4 - 1,
				C4, C4, C3_2, C1_2, C3_2, C2_2, Coef, Coef_location);
			break;
		case (Num3 / 2 - 5) :// 3 lei 2j xia
			kind3_node(j, j - offset, -(Num1 - Node_11 + Num2 + Num3 - 1), -(Num1 - Node_11 + Num2 + Num3 - 1) + 1, -1, 0, 1, 5 + Num4 - 1,
			C4, C4, C3_2, C1_2, C3_2, C2_2, Coef, Coef_location);
			break;
		case (Num3 / 2 - 4) :// 2 lei 2j xia
			kind2_node(j, j - offset, -(Num1 - Node_11 + Num2 + Num3 - 1), -1, 0, 1, 2, 1 + Num4 + Node_55,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num3 / 2 - 1) :// 1 lei 1j xia
			kind1_node(j, j - offset, -(Num1 - Node_11 + Num2 + Num3 - 1), -(Num1 - Node_11 + Num2 + Num3 - 1) + 1, -(Num1 - Node_11 + Num2 + Num3 - 1) + 2, -1, 0,
			1 + Num4 + Node_55, 1 + Num4 + Node_55 + 1, 1 + Num4 + Node_55 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 / 2) < 3)// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num1 - Node_1 - N_a + 3 + Num2 + Num3 / 2), Num3 / 2 + Num4 + Node_5 + N_a - 3, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 / 2) > 4) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 / 2) < (Num3 / 2 - 5)))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num3 - 5), 5 + Num4 - 1, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num1 - Node_11 + Num2 + Num3 - 1), 1 + Num4 + Node_55, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	//cout << "circle num 4:" << j << endl;
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3); i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 / 2); i++)// num 4 :  30 node! part1
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3)
		{
		case 0:// 1 lei 2j shang
			kind1_node(j, j - offset, -(Num3 - 4) - 2, -(Num3 - 4) - 1, -(Num3 - 4), 0, 1, Num4 + N_d0 + 2 - 2, Num4 + N_d0 + 2 - 1, Num4 + N_d0 + 2,
				A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case (Num4 / 2 - 1) :// 1 lei 2j xia
			kind1_node(j, j - offset, -(Num3 - 4), -(Num3 - 4) + 1, -(Num3 - 4) + 2, -1, 0,
			Num4 + N_d0 + 2 - N_c1, Num4 + N_d0 + 2 - N_c1 + 1, Num4 + N_d0 + 2 - N_c1 + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3) <= Node41)// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num3 - 4), Num4 + N_d0 + 2, Coef, Coef_location);
			}
			else if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3) >= Node42) // 1dx
			{
				inner_5node(j, j - offset, 1, -(Num3 - 4), Num4 + N_d0 + 2 - N_c1, Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_4node(j, j - offset, -(Num3 - 4), -1, 0, 1, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 / 2); i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4); i++)// num 4 :  30 node! part2
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3)
		{
		case (Num4 / 2) :// 1 lei 2j shang
			kind1_node(j, j - offset, -(5 + Num4 - 1) - 2, -(5 + Num4 - 1) - 1, -(5 + Num4 - 1), 0, 1,
			(Node_55 - 4 + 1 + N_c1) - 2, (Node_55 - 4 + 1 + N_c1) - 1, (Node_55 - 4 + 1 + N_c1),
			A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case (Num4 - 1) :// 1 lei 2j xia
			kind1_node(j, j - offset, -(5 + Num4 - 1), -(5 + Num4 - 1) + 1, -(5 + Num4 - 1) + 2, -1, 0,
			(Node_55 - 4 + 1), (Node_55 - 4 + 1) + 1, (Node_55 - 4 + 1) + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3) <= Node43)// 1dx
			{
				inner_5node(j, j - offset, 1, -(5 + Num4 - 1), (Node_55 - 4 + 1 + N_c1), Coef, Coef_location);
			}
			else if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3) >= Node44) // 1dx
			{
				inner_5node(j, j - offset, 1, -(5 + Num4 - 1), (Node_55 - 4 + 1), Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_4node(j, j - offset, -(5 + Num4 - 1), -1, 0, 1, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	//cout << "circle num5:" << j << endl;
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4); i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Node_5 + 2); i++)// circle num5 :  80 node!  part1 !
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4)
		{
		case 0:
			boundary_node(j, j - offset, 1, -Num1 - Num2 - Num3 - Num4, Num5 + Num6 + Num7 + Num8, Coef, Coef_location);
			break;
		case (N_d0 - 3) ://2 lei shang 1j
			kind2_node(j, j - offset, -Num1 - Num2 - Num3 - Num4, -2, -1, 0, 1, Num5 + Num6 + Num7 + Num8,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (N_d0 + 1) :// 2 lei 2j shang
			kind2_node(j, j - offset, -(Num3 + Num4 + N_d0 - 2), -2, -1, 0, 1, Num5 - N_d0 + 2 + Num6,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node51:// inner 3 node
			inner_3node(j, j - offset, -(Num4 + N_d0 + 2), -1, 0, 1, 1, k*k*dx*dx - 4, Coef, Coef_location);
			break;
		case (Node51 + 1) :// inner 3 node
			inner_3node(j, j - offset, -(Num4 + N_d0 + 2) + N_c1, 0, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node_5 - 3) ://2 lei xia 2j
			kind2_node(j, j - offset, -(Num3 / 2 + 1 + Num4 + Node_5), -1, 0, 1, 2, Num5 - Node_5 + Num6 + Num7 / 2 - 1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_5 + 1) ://2 lei xia 1j
			kind2_node(j, j - offset, -(Num1 - Node_1 + Num2 + Num3 + Num4 + Node_5), -1, 0, 1, 2, Num5 - Node_5 + Num6 + Num7 + Num8 + Node_9,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < (N_d0 - 3))// 4dx!
			{
				inner_5node(j, j - offset, 4, -Num1 - Num2 - Num3 - Num4, Num5 + Num6 + Num7 + Num8, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) > (N_d0 - 3)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < (N_d0 + 1)))// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Num3 + Num4 + N_d0 - 2), Num5 - N_d0 + 2 + Num6, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) > (N_d0 + 1)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < Node51))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num4 + N_d0 + 2), Num5 - N_d0 - 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) > (Node51 + 1)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < (Node_5 - 3)))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num4 + N_d0 + 2) + N_c1, Num5 - N_d0 - 2 - 2, Coef, Coef_location);
			}
			else// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Num3 / 2 + 1 + Num4 + Node_5), Num5 - Node_5 + Num6 + Num7 / 2 - 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Node_5 + 2); i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5); i++)// circle num5 :  80 node!  part2 !
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4)
		{
		case (Num5 - 1) :
			boundary_node(j, j - offset, 2, -Num2 - Num3 - Num4 - Num5, Num6 + Num7 + Num8 + Num9, Coef, Coef_location);
			break;
		case (Node_5 + N_a - 4) ://2 lei shang 1j
			kind2_node(j, j - offset, -(Num1 - Node_1 + Num2 + Num3 + Num4 + Node_5), -2, -1, 0, 1, Num5 - Node_5 + Num6 + Num7 + Num8 + Node_9,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_5 + N_a) :// 2 lei 2j shang
			kind2_node(j, j - offset, -(Num3 / 2 + Num4 + Node_5 + N_a - 3), -2, -1, 0, 1, Num5 - Node_5 - N_a + 3 + Num6 + Num7 / 2,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node52:// inner 3 node
			inner_3node(j, j - offset, -(Node_55 - 4 + 1) - N_c1, -1, 0, 1, 1, k*k*dx*dx - 4, Coef, Coef_location);
			break;
		case (Node52 + 1) :// inner 3 node
			inner_3node(j, j - offset, -(Node_55 - 4 + 1), 0, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node_55 - 3) ://2 lei xia 2j
			kind2_node(j, j - offset, -(Num4 + Node_55 + 1), -1, 0, 1, 2, Num5 - Node_55 + Num6 + Num7 - 1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_55 + 1) ://2 lei xia 1j
			kind2_node(j, j - offset, -Num2 - Num3 - Num4 - Num5, -1, 0, 1, 2, Num6 + Num7 + Num8 + Num9,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < (Node_5 + N_a - 4))// 4dx!
			{
				inner_5node(j, j - offset, 4, -(Num1 - Node_1 + Num2 + Num3 + Num4 + Node_5), Num5 - Node_5 + Num6 + Num7 + Num8 + Node_9, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) > (Node_5 + N_a - 4)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < (Node_5 + N_a)))// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Num3 / 2 + Num4 + Node_5 + N_a - 3), Num5 - Node_5 - N_a + 3 + Num6 + Num7 / 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) > (Node_5 + N_a)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < Node52))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Node_55 - 4 + 1) - N_c1, (Num5 - Node_55 + 4 + Num6 - 1) + 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) > (Node52 + 1)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < (Node_55 - 3)))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Node_55 - 4 + 1), (Num5 - Node_55 + 4 + Num6 - 1), Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) > (Node_55 - 3)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < (Node_55 + 1)))// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num4 + Node_55 + 1), Num5 - Node_55 + Num6 + Num7 - 1, Coef, Coef_location);
			}
			else// 4dx!
			{
				inner_5node(j, j - offset, 4, -Num2 - Num3 - Num4 - Num5, Num6 + Num7 + Num8 + Num9, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//**************************************************************
	cout << "post1圆孔区域节点数：" << j << endl;
	if (j == sizeM)
		cout << "post1 passed..." << endl;
	else
		cout << "post1 failed..." << endl;
	//修正矩阵，将不在本区域的节点刨除
	for (int i = 0; i != sizeM; i++)
		for (int j = 0; j != N_matrix; j++){
		if (Coef_location[i][j] < 0 || Coef_location[i][j] >= sizeN)
			Coef_location[i][j] = -1;
		}
	if (Debug)
	{
		ofstream outFile_Coef("post1.txt");
		if (!outFile_Coef)
		{
			cout << "Can't write the coefficient matrix into file !" << endl;
			exit(1);//Terminate with fault
		}
		else
		{
			int i, j;
			for (i = 0; i < sizeM; i++)
			{
				for (j = 0; j < N_matrix; j++)
				{
					if ((Coef[i][j].real != 0 || Coef[i][j].image != 0) && Coef_location[i][j] >= 0)
						outFile_Coef << Coef[i][j].real << "+j*" << Coef[i][j].image << " (" << i << "," << Coef_location[i][j] << ") ";
				}
				outFile_Coef << endl << endl;
			}
		}
		outFile_Coef.close();
	}
};
//*********************************************************************************//
//post2大圆孔，长为两孔一缝隙：n代表重叠区域号，Ey[n][M]代表每次迭代所求的解
//*********************************************************************************//
void scan_coef_Matrix_post2(complex** Coef, int num_total, int ** Coef_location, int sizeM, int sizeN, int offset)
{
	int i = 0, j = 0, m = 0;
	for (j = 0; j<num_total; j++)         //系数矩阵初始置零
	{
		for (i = 0; i<N_matrix; i++)
		{
			Coef[j][i].real = 0;
			Coef[j][i].image = 0;
			Coef_location[j][i] = -1;//初始化为负值
		}
	}
	j = 0;
	//系数矩阵扫描
	//******************************************************************************************
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9);
		i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Node_5 + 2); i++)// circle num13(5) :  80 node!  part1 !
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9)
		{
		case 0:
			boundary_node(j, j - offset, 1, -Num9 - Num8 - Num7 - Num6, Num5 + Num4 + Num3 + Num2, Coef, Coef_location);
			break;
		case (N_d0 - 3) ://2 lei shang 1j
			kind2_node(j, j - offset, -Num9 - Num8 - Num7 - Num6, -2, -1, 0, 1, Num5 + Num4 + Num3 + Num2,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (N_d0 + 1) :// 2 lei 2j shang
			kind2_node(j, j - offset, -(Num7 + Num6 + N_d0 - 2), -2, -1, 0, 1, Num5 - N_d0 + 2 + Num4,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node51:// inner 3 node
			inner_3node(j, j - offset, -1, 0, (Num5 - N_d0 - 2), 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node51 + 1) :// inner 3 node
			inner_3node(j, j - offset, 0, 1, (Num5 - N_d0 - 2) + N_c1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Node_5 - 3) ://2 lei xia 2j
			kind2_node(j, j - offset, -(Num7 / 2 + 1 + Num6 + Node_5), -1, 0, 1, 2, Num5 - Node_5 + Num4 + Num3 / 2 - 1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_5 + 1) ://2 lei xia 1j
			kind2_node(j, j - offset, -(Num9 - Node_9 + Num8 + Num7 + Num6 + Node_5), -1, 0, 1, 2, Num5 - Node_5 + Num4 + Num3 + Num2 + Node_1,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < (N_d0 - 3))// 4dx!
			{
				inner_5node(j, j - offset, 4, -Num9 - Num8 - Num7 - Num6, Num5 + Num4 + Num3 + Num2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) > (N_d0 - 3)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < (N_d0 + 1)))// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Num7 + Num6 + N_d0 - 2), Num5 - N_d0 + 2 + Num4, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) > (N_d0 + 1)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < Node51))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num6 + N_d0 + 2), (Num5 - N_d0 - 2), Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) > (Node51 + 1)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < (Node_5 - 3)))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num6 + N_d0 + 2) - 2, (Num5 - N_d0 - 2) + N_c1, Coef, Coef_location);
			}
			else// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Num7 / 2 + 1 + Num6 + Node_5), Num5 - Node_5 + Num4 + Num3 / 2 - 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Node_5 + 2);
		i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9); i++)// circle num13(5) :  80 node!  part2 !
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9)
		{
		case (Num5 - 1) :
			boundary_node(j, j - offset, 2, -Num8 - Num7 - Num6 - Num5, Num4 + Num3 + Num2 + Num1, Coef, Coef_location);
			break;
		case (Node_5 + N_a - 4) ://2 lei shang 1j
			kind2_node(j, j - offset, -(Num9 - Node_9 + Num8 + Num7 + Num6 + Node_5), -2, -1, 0, 1, Num5 - Node_5 + Num4 + Num3 + Num2 + Node_1,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_5 + N_a) :// 2 lei 2j shang
			kind2_node(j, j - offset, -(Num7 / 2 + Num6 + Node_5 + N_a - 3), -2, -1, 0, 1, Num5 - Node_5 - N_a + 3 + Num4 + Num3 / 2,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node52:// inner 3 node
			inner_3node(j, j - offset, -1, 0, (Num5 - Node_55 + 4 + Num4 - 1) - N_c1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node52 + 1) :// inner 3 node
			inner_3node(j, j - offset, 0, 1, (Num5 - Node_55 + 4 + Num4 - 1), k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Node_55 - 3) ://2 lei xia 2j
			kind2_node(j, j - offset, -(Num6 + Node_55 + 1), -1, 0, 1, 2, Num5 - Node_55 + Num4 + Num3 - 1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_55 + 1) ://2 lei xia 1j
			kind2_node(j, j - offset, -Num8 - Num7 - Num6 - Num5, -1, 0, 1, 2, Num4 + Num3 + Num2 + Num1,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < (Node_5 + N_a - 4))// 4dx!
			{
				inner_5node(j, j - offset, 4, -(Num9 - Node_9 + Num8 + Num7 + Num6 + Node_5), Num5 - Node_5 + Num4 + Num3 + Num2 + Node_1, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) > (Node_5 + N_a - 4)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < (Node_5 + N_a)))// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Num7 / 2 + Num6 + Node_5 + N_a - 3), Num5 - Node_5 - N_a + 3 + Num4 + Num3 / 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) > (Node_5 + N_a)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < Node52))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Node_55 - 4 + 1) + 2, (Num5 - Node_55 + 4 + Num4 - 1) - N_c1, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) > (Node52 + 1)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < (Node_55 - 3)))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Node_55 - 4 + 1), (Num5 - Node_55 + 4 + Num4 - 1), Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) > (Node_55 - 3)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < (Node_55 + 1)))// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num6 + Node_55 + 1), Num5 - Node_55 + Num4 + Num3 - 1, Coef, Coef_location);
			}
			else// 4dx!
			{
				inner_5node(j, j - offset, 4, -Num8 - Num7 - Num6 - Num5, Num4 + Num3 + Num2 + Num1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9);
		i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Num4 / 2); i++)// num 14(4) :  30 node! part1
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9)
		{
		case 0:// 1 lei 2j shang
			kind1_node(j, j - offset, -(Num5 - N_d0 - 2) - 2, -(Num5 - N_d0 - 2) - 1, -(Num5 - N_d0 - 2), 0, 1, Num4 + 4 - 2, Num4 + 4 - 1, Num4 + 4,
				A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case (Num4 / 2 - 1) :// 1 lei 2j xia
			kind1_node(j, j - offset, -(Num5 - N_d0 - 2) - N_c1, -(Num5 - N_d0 - 2) - N_c1 + 1, -(Num5 - N_d0 - 2) - N_c1 + 2, -1, 0,
			Num4 + 4, Num4 + 4 + 1, Num4 + 4 + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) <= Node41)// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num5 - N_d0 - 2), Num4 + 4, Coef, Coef_location);
			}
			else if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) >= Node42) // 1dx
			{
				inner_5node(j, j - offset, 1, -(Num5 - N_d0 - 2) - N_c1, Num4 + 4, Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_4node(j, j - offset, -1, 0, 1, Num4 + 4, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Num4 / 2);
		i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9); i++)// num14(4) :  30 node! part2
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9)
		{
		case (Num4 / 2) :// 1 lei 2j shang
			kind1_node(j, j - offset, -(Num5 - Node_55 + 4 + Num4 - 1) + N_c1 - 2, -(Num5 - Node_55 + 4 + Num4 - 1) + N_c1 - 1, -(Num5 - Node_55 + 4 + Num4 - 1) + N_c1, 0, 1,
			(Num3 - 5 + 1) - 2, (Num3 - 5 + 1) - 1, (Num3 - 5 + 1),
			A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case (Num4 - 1) :// 1 lei 2j xia
			kind1_node(j, j - offset, -(Num5 - Node_55 + 4 + Num4 - 1), -(Num5 - Node_55 + 4 + Num4 - 1) + 1, -(Num5 - Node_55 + 4 + Num4 - 1) + 2, -1, 0,
			(Num3 - 5 + 1), (Num3 - 5 + 1) + 1, (Num3 - 5 + 1) + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) <= Node43)// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num5 - Node_55 + 4 + Num4 - 1) + N_c1, (Num3 - 5 + 1), Coef, Coef_location);
			}
			else if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) >= Node44) // 1dx
			{
				inner_5node(j, j - offset, 1, -(Num5 - Node_55 + 4 + Num4 - 1), (Num3 - 5 + 1), Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_4node(j, j - offset, -1, 0, 1, (Num3 - 5 + 1), 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9);
		i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Num3 / 2); i++)// num 15(3) :  46 node! part1
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9)
		{
		case 0:// 1 lei 1j shang
			kind1_node(j, j - offset, -(Num5 - N_d0 + 2 + Num4) - 2, -(Num5 - N_d0 + 2 + Num4) - 1, -(Num5 - N_d0 + 2 + Num4), 0, 1,
				Num3 + Num2 + N_d0 - 2 - 2, Num3 + Num2 + N_d0 - 2 - 1, Num3 + Num2 + N_d0 - 2,
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case 3:// 2 lei 2j shang
			kind2_node(j, j - offset, -(Num5 - N_d0 + 2 + Num4), -2, -1, 0, 1, Num3 + Num2 + N_d0 - 2,
				B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case 4:// 3 lei 2j shang
			kind3_node(j, j - offset, -(Num4 + 4), -1, 0, 1, (Num3 + Num2 + N_d0 - 2) - 1, (Num3 + Num2 + N_d0 - 2),
				C2_2, C3_2, C1_2, C3_2, C4, C4, Coef, Coef_location);
			break;
		case (Num3 / 2 - 5) :// 3 lei 2j xia
			kind3_node(j, j - offset, -(Num4 + 4), -1, 0, 1, Num3 / 2 + 1 + Num2 + Node_1, Num3 / 2 + 1 + Num2 + Node_1 + 1,
			C2_2, C3_2, C1_2, C3_2, C4, C4, Coef, Coef_location);
			break;
		case (Num3 / 2 - 4) :// 2 lei 2j xia
			kind2_node(j, j - offset, -(Num5 - Node_5 + Num4 + Num3 / 2 - 1), -1, 0, 1, 2, Num3 / 2 + 1 + Num2 + Node_1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num3 / 2 - 1) :// 1 lei 1j xia
			kind1_node(j, j - offset, -(Num5 - Node_5 + Num4 + Num3 / 2 - 1), -(Num5 - Node_5 + Num4 + Num3 / 2 - 1) + 1, -(Num5 - Node_5 + Num4 + Num3 / 2 - 1) + 2, -1, 0,
			Num3 / 2 + 1 + Num2 + Node_1, Num3 / 2 + 1 + Num2 + Node_1 + 1, Num3 / 2 + 1 + Num2 + Node_1 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < 3)// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num5 - N_d0 + 2 + Num4), Num3 + Num2 + N_d0 - 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) > 4) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) < (Num3 / 2 - 5)))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num4 + 4), Num3 - 5, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num5 - Node_5 + Num4 + Num3 / 2 - 1), Num3 / 2 + 1 + Num2 + Node_1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Num3 / 2);
		i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 * 2 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9); i++)// num 15(3) :  46 node! part2
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num3 / 2)
		{
		case 0:// 1 lei 1j shang
			kind1_node(j, j - offset, -(Num5 - Node_5 - N_a + 3 + Num4 + Num3 / 2) - 2, -(Num5 - Node_5 - N_a + 3 + Num4 + Num3 / 2) - 1, -(Num5 - Node_5 - N_a + 3 + Num4 + Num3 / 2), 0, 1,
				Num3 / 2 + Num2 + Node_1 + N_a - 3 - 2, Num3 / 2 + Num2 + Node_1 + N_a - 3 - 1, Num3 / 2 + Num2 + Node_1 + N_a - 3,
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case 3:// 2 lei 2j shang
			kind2_node(j, j - offset, -(Num5 - Node_5 - N_a + 3 + Num4 + Num3 / 2), -2, -1, 0, 1, Num3 / 2 + Num2 + Node_1 + N_a - 3,
				B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case 4:// 3 lei 2j shang
			kind3_node(j, j - offset, -(Num3 - 5 + 1), -1, 0, 1, (Num3 / 2 + Num2 + Node_1 + N_a - 3) - 1, (Num3 / 2 + Num2 + Node_1 + N_a - 3),
				C2_2, C3_2, C1_2, C3_2, C4, C4, Coef, Coef_location);
			break;
		case (Num3 / 2 - 5) :// 3 lei 2j xia
			kind3_node(j, j - offset, -(Num3 - 5 + 1), -1, 0, 1, 1 + Num2 + Node_11, 1 + Num2 + Node_11 + 1,
			C2_2, C3_2, C1_2, C3_2, C4, C4, Coef, Coef_location);
			break;
		case (Num3 / 2 - 4) :// 2 lei 2j xia
			kind2_node(j, j - offset, -(Num5 - Node_55 + Num4 + Num3 - 1), -1, 0, 1, 2, 1 + Num2 + Node_11,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num3 / 2 - 1) :// 1 lei 1j xia
			kind1_node(j, j - offset, -(Num5 - Node_55 + Num4 + Num3 - 1), -(Num5 - Node_55 + Num4 + Num3 - 1) + 1, -(Num5 - Node_55 + Num4 + Num3 - 1) + 2, -1, 0,
			1 + Num2 + Node_11, 1 + Num2 + Node_11 + 1, 1 + Num2 + Node_11 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num3 / 2) < 3)// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num5 - Node_5 - N_a + 3 + Num4 + Num3 / 2), Num3 / 2 + Num2 + Node_1 + N_a - 3, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num3 / 2) > 4) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num3 / 2) < (Num3 / 2 - 5)))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num3 - 5 + 1), 6 + Num2 - 1, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num5 - Node_55 + Num4 + Num3 - 1), 1 + Num2 + Node_11, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 * 2 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9);
		i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 * 2 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Num2 / 2); i++)// num 16(2) :  26 node! part1
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9)
		{
		case 0:// 3 lei 2j shang
			kind3_node(j, j - offset, -(Num3 - 3), -(Num3 - 3) + 2, 0, 1, Num2 + Node11 - 1, Num2 + Node11,
				C4, C3_2, C1_2, C2_2, C4, C3_2, Coef, Coef_location);
			break;
		case (Num2 / 2 - 1) ://node2: 3 lei 2j xia
			kind3_node(j, j - offset, -(Num3 - 5), -(Num3 - 5) + 2, -1, 0, Num2 / 2 + Node_1 - 4 + 1, Num2 / 2 + Node_1 - 4 + 1 + 1,
			C3_2, C4, C2_2, C1_2, C3_2, C4, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) % 2 == 0)// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num3 - 5),
					(Num2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) + Node11 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) / 2), Coef, Coef_location);
			}
			else//1 lei  2j left
			{
				kind1_node(j, j - offset, -(Num3 - 5), -1, 0, 1,
					(Num2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) + Node11 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) / 2),
					(Num2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) + Node11 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) / 2) + 1,
					(Num2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) + Node11 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) / 2) + (Num1 - N_d0 + 2),
					(Num2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) + Node11 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9) / 2) + 1 + (Num1 - N_d0 + 2),
					A4_2, A3_2, A5_2, A3_2, A1_2, A1_2, A2_2, A2_2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 * 2 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Num2 / 2);
		i<(Nx_2 + Nx_3 + Num1 + Num2 * 2 + Num3 * 2 + Num4 * 2 + Num5 * 2 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9); i++)// num 16(2) :  26 node! part2
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2)
		{
		case 0:// 3 lei 2j shang
			kind3_node(j, j - offset, -(6 + Num2 - 1) - 2, -(6 + Num2 - 1), 0, 1, Num2 / 2 + Node13 - 1, Num2 / 2 + Node13,
				C4, C3_2, C1_2, C2_2, C4, C3_2, Coef, Coef_location);
			break;
		case (Num2 / 2 - 1) :// 3 lei 2j xia
			kind3_node(j, j - offset, -(6 + Num2 - 1), -(6 + Num2 - 1) + 2, -1, 0, Node_11 - 4 + 1, Node_11 - 4 + 1 + 1,
			C3_2, C4, C2_2, C1_2, C3_2, C4, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) % 2 == 0)// 1dx
			{
				inner_5node(j, j - offset, 1, -(6 + Num2 - 1), Num2 / 2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) + Node13 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) / 2,
					Coef, Coef_location);
			}
			else//1 lei 2j left
			{
				kind1_node(j, j - offset, -(6 + Num2 - 1), -1, 0, 1,
					(Num2 / 2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) + Node13 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) / 2),
					(Num2 / 2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) + Node13 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) / 2) + 1,
					(Num2 / 2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) + Node13 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) / 2) + (Num1 - Node_11 + Nx_3 - 1),
					(Num2 / 2 - (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) + Node13 + (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 * 2 - Num4 * 2 - Num5 * 2 - Num6 * 2 - Num7 * 2 - Num8 * 2 - Num9 - Num2 / 2) / 2) + 1 + (Num1 - Node_11 + Nx_3 - 1),
					A4_2, A3_2, A5_2, A3_2, A1_2, A1_2, A2_2, A2_2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (num_total - Nx_2 - Nx_3 - Num1); i<(num_total - Nx_2 - Nx_3 - Num1 + Node_1 + 2); i++)// circle num17(1) Nx_2:  74 node!  part1 !
	{
		switch (i - num_total + Nx_2 + Nx_3 + Num1)
		{
		case 0:
			boundary_node(j, j - offset, 1, -Num5 - Num4 - Num3 - Num2, Nx_2 + Nx_3, Coef, Coef_location);
			break;
		case (N_d0 - 3) ://2 lei shang 1j
			kind2_node(j, j - offset, -Num5 - Num4 - Num3 - Num2, -2, -1, 0, 1, Nx_2 + Nx_3,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_1 + 1) ://2 lei xia 1j
			kind2_node(j, j - offset, -(Num5 - Node_5 + Num4 + Num3 + Num2 + Node_1), -1, 0, 1, 2, Nx_2 + Nx_3,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - num_total + Nx_2 + Nx_3 + Num1) < (N_d0 - 3))// 4dx!
			{
				inner_5node(j, j - offset, 4, -Num5 - Num4 - Num3 - Num2, Nx_2 + Nx_3, Coef, Coef_location);
			}
			else if (((i - num_total + Nx_2 + Nx_3 + Num1) > (N_d0 - 3)) && ((i - num_total + Nx_2 + Nx_3 + Num1) < Node11))// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Num3 + Num2 + N_d0 - 2), Num1 - N_d0 + 2, Coef, Coef_location);
			}
			else if (((i - num_total + Nx_2 + Nx_3 + Num1) >= Node11) && ((i - num_total + Nx_2 + Nx_3 + Num1) <= Node12))// 2 lei right 2j !
			{
				kind2_node(j, j - offset, -((i - num_total + Nx_2 + Nx_3 + Num1) + Num2 - (i - num_total + Nx_2 + Nx_3 + Num1 - Node11) * 2),
					-1, 0, 1, Num1 - N_d0 + 2, Nx_2 + Nx_3, B2, B3, B1_2, B3, B4, B5, Coef, Coef_location);
			}
			else// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Num3 / 2 + 1 + Num2 + Node_1), Num1 - N_d0 + 2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (num_total - Nx_2 - Nx_3 - Num1 + Node_1 + 2); i<(num_total - Nx_2 - Nx_3); i++)// circle num17(1) Nx_2:  74 node!  part2 !
	{
		switch (i - num_total + Nx_2 + Nx_3 + Num1)
		{
		case (Num1 - 1) :
			boundary_node(j, j - offset, 2, -Num1 - Num2 - Num3 - Num4, Nx_2 + Nx_3, Coef, Coef_location);
			break;
		case (Node_1 + N_a - 4) ://2 lei shang 1j
			kind2_node(j, j - offset, -(Num5 - Node_5 + Num4 + Num3 + Num2 + Node_1), -2, -1, 0, 1, Nx_2 + Nx_3,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_11 + 1) ://2 lei xia 1j
			kind2_node(j, j - offset, -Num1 - Num2 - Num3 - Num4, -1, 0, 1, 2, Nx_2 + Nx_3,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - num_total + Nx_2 + Nx_3 + Num1) < (Node_1 + N_a - 4))// 4dx!
			{
				inner_5node(j, j - offset, 4, -(Num5 - Node_5 + Num4 + Num3 + Num2 + Node_1), Nx_2 + Nx_3, Coef, Coef_location);
			}
			else if (((i - num_total + Nx_2 + Nx_3 + Num1) > (Node_1 + N_a - 4)) && ((i - num_total + Nx_2 + Nx_3 + Num1) < Node13))// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Num3 / 2 + Num2 + Node_1 + N_a - 3), Num1 - Node_11 + Nx_3 - 1, Coef, Coef_location);
			}
			else if (((i - num_total + Nx_2 + Nx_3 + Num1) >= Node13) && ((i - num_total + Nx_2 + Nx_3 + Num1) <= Node14))// 2 lei left 2j !
			{
				kind2_node(j, j - offset, -((i - num_total + Nx_2 + Nx_3 + Num1) + Num2 / 2 - (i - num_total + Nx_2 + Nx_3 + Num1 - Node13) * 2),
					-1, 0, 1, Num1 - Node_11 + Nx_3 - 1, Nx_2 + Nx_3, B2, B3, B1_2, B3, B4, B5, Coef, Coef_location);
			}
			else if (((i - num_total + Nx_2 + Nx_3 + Num1) >= Node14) && ((i - num_total + Nx_2 + Nx_3 + Num1) <= Node_11))// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num2 + Node_11 + 1), Num1 - Node_11 + Nx_3 - 1, Coef, Coef_location);
			}
			else
			{
				inner_5node(j, j - offset, 4, -Num1 - Num2 - Num3 - Num4, Nx_2 + Nx_3, Coef, Coef_location);// 4dx!
			}
			break;
		}
		j = j + 1;
	}
	//*****************************************************************************************
	for (i = (Nx_2); i<(Nx_2 + Nx_3); i++)//left 2 Nx_3:  30 node! 
	{
		switch (i - Nx_2)
		{
		case 0:// 1 lei shang
			kind1_node(j, i, -(Nx_2 - N_d0 + 2) - 2, -(Nx_2 - N_d0 + 2) - 1, -(Nx_2 - N_d0 + 2), 0, 1, Nx_3 + N_d0 - 2 - 2, Nx_3 + N_d0 - 2 - 1, Nx_3 + N_d0 - 2,
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case (Nx_3 / 2 - 1) :// 1 lei xia
			kind1_node(j, i, -(Nx_2 - N_d0 + 2), -(Nx_2 - N_d0 + 2) + 1, -(Nx_2 - N_d0 + 2) + 2, -1, 0,
			Nx_3 + N_d0 - 2, Nx_3 + N_d0 - 2 + 1, Nx_3 + N_d0 - 2 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		case (Nx_3 / 2) :// 1 lei shang
			kind1_node(j, i, -(Nx_2 - Node_11_p2 - 1 + Nx_3) - 2, -(Nx_2 - Node_11_p2 - 1 + Nx_3) - 1, -(Nx_2 - Node_11_p2 - 1 + Nx_3),
			0, 1, 1 + Node_11_p2 - 2, 1 + Node_11_p2 - 1, 1 + Node_11_p2,
			A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case (Nx_3 - 1) :// 1 lei xia
			kind1_node(j, i, -(Nx_2 - Node_11_p2 - 1 + Nx_3), -(Nx_2 - Node_11_p2 - 1 + Nx_3) + 1, -(Nx_2 - Node_11_p2 - 1 + Nx_3) + 2,
			-1, 0, 1 + Node_11_p2, 1 + Node_11_p2 + 1, 1 + Node_11_p2 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2) < (Nx_3 / 2 - 1))
			{
				inner_5node(j, i, 2, -(Nx_2 - N_d0 + 2), Nx_3 + N_d0 - 2, Coef, Coef_location);
			}
			else
			{
				inner_5node(j, i, 2, -(Nx_2 - Node_11_p2 - 1 + Nx_3), 1 + Node_11_p2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (Nx_2 + Nx_3); i<(Nx_2 + Nx_3 + Node_1_p2 + 2); i++)// circle num1 Nx_2:  74 node!  part1 !
	{
		switch (i - Nx_2 - Nx_3)
		{
		case 0:
			boundary_node(j, i, 1, -Nx_2 - Nx_3, Num1_p2 + Num2_p2 + Num3_p2 + Num4_p2, Coef, Coef_location);
			break;
		case (N_d0 - 3) ://2 lei shang 1j
			kind2_node(j, i, -Nx_2 - Nx_3, -2, -1, 0, 1, Num1_p2 + Num2_p2 + Num3_p2 + Num4_p2,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_1_p2 + 1) ://2 lei xia 1j
			kind2_node(j, i, -(Nx_2 + Nx_3), -1, 0, 1, 2, Num1_p2 - Node_1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + Node_5_p2,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3) < (N_d0 - 3))// 4dx!
			{
				inner_5node(j, i, 4, -Nx_3 - Nx_2, Num1_p2 + Num2_p2 + Num3_p2 + Num4_p2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3) > (N_d0 - 3)) && ((i - Nx_2 - Nx_3) < Node11_p2))// 2dx!
			{
				inner_5node(j, i, 2, -(Nx_3 + N_d0 - 2), Num1_p2 - N_d0 + 2 + Num2_p2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3) >= Node11_p2) && ((i - Nx_2 - Nx_3) <= Node12_p2))// 2 lei left 2j !
			{
				kind2_node(j, i, -Nx_3 - Nx_2, -(Nx_3 + N_d0 - 2), -1, 0, 1, Num1_p2 - (i - Nx_2 - Nx_3) + ((i - Nx_2 - Nx_3) - Node11_p2) * 2,
					B5, B4, B3, B1_2, B3, B2, Coef, Coef_location);
			}
			else// 2dx!
			{
				inner_5node(j, i, 2, -(Nx_3 + N_d0 - 2), Num1_p2 - Node_1_p2 + Num2_p2 + Num3_p2 / 2 - 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Node_1_p2 + 2); i<(Nx_2 + Nx_3 + Num1_p2); i++)// circle num1 Nx_2:  74 node!  part2 !
	{
		switch (i - Nx_2 - Nx_3)
		{
		case (Num1_p2 - 1) :
			boundary_node(j, i, 2, -Nx_2 - Nx_3, Num5_p2 + Num2_p2 + Num3_p2 + Num4_p2, Coef, Coef_location);
			break;
		case (Node_1_p2 + N_a - 4) ://2 lei shang 1j
			kind2_node(j, i, -Nx_2 - Nx_3, -2, -1, 0, 1, Num1_p2 - Node_1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + Node_5_p2,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_11_p2 + 1) ://2 lei xia 1j
			kind2_node(j, i, -(Nx_2 + Nx_3), -1, 0, 1, 2, Num5_p2 + Num2_p2 + Num3_p2 + Num4_p2,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3) < (Node_1_p2 + N_a - 4))// 4dx!
			{
				inner_5node(j, i, 4, -Nx_3 - Nx_2, Num1_p2 - Node_1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + Node_5_p2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3) > (Node_1_p2 + N_a - 4)) && ((i - Nx_2 - Nx_3) < Node13_p2))// 2dx!
			{
				inner_5node(j, i, 2, -(Node_11_p2 + 1), Num1_p2 - Node_1_p2 - N_a + 3 + Num2_p2 + Num3_p2 / 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3) >= Node13_p2) && ((i - Nx_2 - Nx_3) <= Node14_p2))// 2 lei left 2j !
			{
				kind2_node(j, i, -Nx_3 - Nx_2, -(Node_11_p2 + 1), -1, 0, 1, Num1_p2 - (i - Nx_2 - Nx_3) + Num2_p2 / 2 + ((i - Nx_2 - Nx_3) - Node13_p2) * 2,
					B5, B4, B3, B1_2, B3, B2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3) >= Node14_p2) && ((i - Nx_2 - Nx_3) <= Node_11_p2))// 2dx
			{
				inner_5node(j, i, 2, -(Node_11_p2 + 1), Num1_p2 - Node_11_p2 + Num2_p2 + Num3_p2 - 1, Coef, Coef_location);
			}
			else
			{
				inner_5node(j, i, 4, -(Nx_3 + Nx_2), Num5_p2 + Num2_p2 + Num3_p2 + Num4_p2, Coef, Coef_location);// 4dx!
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!!!!!!!!!
	for (i = (Nx_2 + Nx_3 + Num1_p2); i<(Nx_2 + Nx_3 + Num1_p2 + Num2_p2 / 2); i++)// num 2 :  26 node! part1
	{
		switch (i - Nx_2 - Nx_3 - Num1_p2)
		{
		case 0:// 3 lei 2j shang
			kind3_node(j, i, -(Num1_p2 - Node11_p2) - 1, -(Num1_p2 - Node11_p2), 0, 1, Num2_p2 + 5 - 2, Num2_p2 + 5,
				C4, C3_2, C1_2, C2_2, C4, C3_2, Coef, Coef_location);
			break;
		case (Num2_p2 / 2 - 1) ://node2: 3 lei 2j xia
			kind3_node(j, i, -(Num1_p2 - Node12_p2 + Num2_p2 / 2 - 1), -(Num1_p2 - Node12_p2 + Num2_p2 / 2 - 1) + 1, -1, 0,
			Num2_p2 + 5, Num2_p2 + 5 + 2, C3_2, C4, C2_2, C1_2, C3_2, C4, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1_p2) % 2 == 0)// 1dx
			{
				inner_5node(j, i, 1, -((i - Nx_2 - Nx_3 - Num1_p2) + Num1_p2 - Node11_p2 - (i - Nx_2 - Nx_3 - Num1_p2) / 2), Num2_p2 + 5, Coef, Coef_location);
			}
			else//1 lei  2j left
			{
				kind1_node(j, i, -((i - Nx_2 - Nx_3 - Num1_p2) + Num1_p2 - Node11_p2 - (i - Nx_2 - Nx_3 - Num1_p2) / 2) - (Nx_3 + N_d0 - 2),
					-((i - Nx_2 - Nx_3 - Num1_p2) + Num1_p2 - Node11_p2 - (i - Nx_2 - Nx_3 - Num1_p2) / 2) - (Nx_3 + N_d0 - 2) + 1,
					-((i - Nx_2 - Nx_3 - Num1_p2) + Num1_p2 - Node11_p2 - (i - Nx_2 - Nx_3 - Num1_p2) / 2),
					-((i - Nx_2 - Nx_3 - Num1_p2) + Num1_p2 - Node11_p2 - (i - Nx_2 - Nx_3 - Num1_p2) / 2) + 1,
					-1, 0, 1, 5 + Num2_p2, A2_2, A2_2, A1_2, A1_2, A3_2, A5_2, A3_2, A4_2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1_p2 + Num2_p2 / 2); i<(Nx_2 + Nx_3 + Num1_p2 + Num2_p2); i++)// num 2 :  26 node! part2
	{
		switch (i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 / 2)
		{
		case 0:// 3 lei 2j shang
			kind3_node(j, i, -(Num1_p2 - Node13_p2 + Num2_p2 / 2) - 1, -(Num1_p2 - Node13_p2 + Num2_p2 / 2), 0, 1,
				Num3_p2 - 5 - 2, Num3_p2 - 5, C4, C3_2, C1_2, C2_2, C4, C3_2, Coef, Coef_location);
			break;
		case (Num2_p2 / 2 - 1) :// 3 lei 2j xia
			kind3_node(j, i, -(Num1_p2 - Node14_p2 + Num2_p2 - 1), -(Num1_p2 - Node14_p2 + Num2_p2 - 1) + 1, -1, 0,
			Num3_p2 - 5, Num3_p2 - 5 + 2, C3_2, C4, C2_2, C1_2, C3_2, C4, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 / 2) % 2 == 0)// 1dx
			{
				inner_5node(j, i, 1, -((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 / 2) + Num1_p2 - Node13_p2 - (i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 / 2) / 2) - Num2_p2 / 2,
					Num3_p2 - 5, Coef, Coef_location);
			}
			else//1 lei 2j left
			{
				kind1_node(j, i, -((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 / 2) + Num1_p2 - Node13_p2 - (i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 / 2) / 2) - Num2_p2 / 2 - (Node_11_p2 + 1),
					-((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 / 2) + Num1_p2 - Node13_p2 - (i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 / 2) / 2) - Num2_p2 / 2 - (Node_11_p2 + 1) + 1,
					-((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 / 2) + Num1_p2 - Node13_p2 - (i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 / 2) / 2) - Num2_p2 / 2,
					-((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 / 2) + Num1_p2 - Node13_p2 - (i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 / 2) / 2) - Num2_p2 / 2 + 1,
					-1, 0, 1, Num3_p2 - 5, A2_2, A2_2, A1_2, A1_2, A3_2, A5_2, A3_2, A4_2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!!!
	for (i = (Nx_2 + Nx_3 + Num1_p2 + Num2_p2); i<(Nx_2 + Nx_3 + Num1_p2 + Num2_p2 + Num3_p2 / 2); i++)// num 3 :  46 node! part1
	{
		switch (i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2)
		{
		case 0:// 1 lei 1j shang
			kind1_node(j, i, -(Num1_p2 - N_d0 + 2 + Num2_p2) - 2, -(Num1_p2 - N_d0 + 2 + Num2_p2) - 1, -(Num1_p2 - N_d0 + 2 + Num2_p2),
				0, 1, Num3_p2 + Num4_p2 + N_d0 - 2 - 2, Num3_p2 + Num4_p2 + N_d0 - 2 - 1, Num3_p2 + Num4_p2 + N_d0 - 2,
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case 3:// 2 lei 2j shang
			kind2_node(j, i, -(Num1_p2 - N_d0 + 2 + Num2_p2), -2, -1, 0, 1, Num3_p2 + Num4_p2 + N_d0 - 2,
				B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case 4:// 3 lei 2j shang
			kind3_node(j, i, -(Num1_p2 - N_d0 + 2 + Num2_p2) - 1, -(Num1_p2 - N_d0 + 2 + Num2_p2), -1, 0, 1, Num3_p2 - 4,
				C4, C4, C3_2, C1_2, C3_2, C2_2, Coef, Coef_location);
			break;
		case (Num3_p2 / 2 - 5) :// 3 lei 2j xia
			kind3_node(j, i, -(Num1_p2 - Node_1_p2 + Num2_p2 + Num3_p2 / 2 - 1), -(Num1_p2 - Node_1_p2 + Num2_p2 + Num3_p2 / 2 - 1) + 1,
			-1, 0, 1, Num3_p2 - 4, C4, C4, C3_2, C1_2, C3_2, C2_2, Coef, Coef_location);
			break;
		case (Num3_p2 / 2 - 4) :// 2 lei 2j xia
			kind2_node(j, i, -(Num1_p2 - Node_1_p2 + Num2_p2 + Num3_p2 / 2 - 1), -1, 0, 1, 2, Num3_p2 / 2 + 1 + Num4_p2 + Node_5_p2,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num3_p2 / 2 - 1) :// 1 lei 1j xia
			kind1_node(j, i, -(Num1_p2 - Node_1_p2 + Num2_p2 + Num3_p2 / 2 - 1), -(Num1_p2 - Node_1_p2 + Num2_p2 + Num3_p2 / 2 - 1) + 1,
			-(Num1_p2 - Node_1_p2 + Num2_p2 + Num3_p2 / 2 - 1) + 2, -1, 0,
			Num3_p2 / 2 + 1 + Num4_p2 + Node_5_p2, Num3_p2 / 2 + 1 + Num4_p2 + Node_5_p2 + 1, Num3_p2 / 2 + 1 + Num4_p2 + Node_5_p2 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2) < 3)// 2dx
			{
				inner_5node(j, i, 2, -(Num1_p2 - N_d0 + 2 + Num2_p2), Num3_p2 + Num4_p2 + N_d0 - 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2) > 4) && ((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2) < (Num3_p2 / 2 - 5)))// 1dx
			{
				inner_5node(j, i, 1, -(Num2_p2 + 5), Num3_p2 - 4, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, i, 2, -(Num1_p2 - Node_1_p2 + Num2_p2 + Num3_p2 / 2 - 1), Num3_p2 / 2 + 1 + Num4_p2 + Node_5_p2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1_p2 + Num2_p2 + Num3_p2 / 2); i<(Nx_2 + Nx_3 + Num1_p2 + Num2_p2 + Num3_p2); i++)// num 3 :  46 node! part2
	{
		switch (i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 / 2)
		{
		case 0:// 1 lei 1j shang
			kind1_node(j, i, -(Num1_p2 - Node_1_p2 - N_a + 3 + Num2_p2 + Num3_p2 / 2) - 2, -(Num1_p2 - Node_1_p2 - N_a + 3 + Num2_p2 + Num3_p2 / 2) - 1,
				-(Num1_p2 - Node_1_p2 - N_a + 3 + Num2_p2 + Num3_p2 / 2), 0, 1, Num3_p2 / 2 + Num4_p2 + Node_5_p2 + N_a - 3 - 2,
				Num3_p2 / 2 + Num4_p2 + Node_5_p2 + N_a - 3 - 1, Num3_p2 / 2 + Num4_p2 + Node_5_p2 + N_a - 3,
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case 3:// 2 lei 2j shang
			kind2_node(j, i, -(Num1_p2 - Node_1_p2 - N_a + 3 + Num2_p2 + Num3_p2 / 2), -2, -1, 0, 1,
				Num3_p2 / 2 + Num4_p2 + Node_5_p2 + N_a - 3, B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case 4:// 3 lei 2j shang
			kind3_node(j, i, -(Num1_p2 - Node_1_p2 - N_a + 3 + Num2_p2 + Num3_p2 / 2) - 1, -(Num1_p2 - Node_1_p2 - N_a + 3 + Num2_p2 + Num3_p2 / 2),
				-1, 0, 1, 5 + Num4_p2 - 1, C4, C4, C3_2, C1_2, C3_2, C2_2, Coef, Coef_location);
			break;
		case (Num3_p2 / 2 - 5) :// 3 lei 2j xia
			kind3_node(j, i, -(Num1_p2 - Node_11_p2 + Num2_p2 + Num3_p2 - 1), -(Num1_p2 - Node_11_p2 + Num2_p2 + Num3_p2 - 1) + 1,
			-1, 0, 1, 5 + Num4_p2 - 1, C4, C4, C3_2, C1_2, C3_2, C2_2, Coef, Coef_location);
			break;
		case (Num3_p2 / 2 - 4) :// 2 lei 2j xia
			kind2_node(j, i, -(Num1_p2 - Node_11_p2 + Num2_p2 + Num3_p2 - 1), -1, 0, 1, 2, 1 + Num4_p2 + Node_55_p2,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num3_p2 / 2 - 1) :// 1 lei 1j xia
			kind1_node(j, i, -(Num1_p2 - Node_11_p2 + Num2_p2 + Num3_p2 - 1), -(Num1_p2 - Node_11_p2 + Num2_p2 + Num3_p2 - 1) + 1,
			-(Num1_p2 - Node_11_p2 + Num2_p2 + Num3_p2 - 1) + 2, -1, 0,
			1 + Num4_p2 + Node_55_p2, 1 + Num4_p2 + Node_55_p2 + 1, 1 + Num4_p2 + Node_55_p2 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 / 2) < 3)// 2dx
			{
				inner_5node(j, i, 2, -(Num1_p2 - Node_1_p2 - N_a + 3 + Num2_p2 + Num3_p2 / 2), Num3_p2 / 2 + Num4_p2 + Node_5_p2 + N_a - 3, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 / 2) > 4) && ((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 / 2) < (Num3_p2 / 2 - 5)))// 1dx
			{
				inner_5node(j, i, 1, -(Num3_p2 - 5), 5 + Num4_p2 - 1, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, i, 2, -(Num1_p2 - Node_11_p2 + Num2_p2 + Num3_p2 - 1), 1 + Num4_p2 + Node_55_p2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (Nx_2 + Nx_3 + Num1_p2 + Num2_p2 + Num3_p2);
		i<(Nx_2 + Nx_3 + Num1_p2 + Num2_p2 + Num3_p2 + Num4_p2 / 2); i++)// num 4 :  30 node! part1
	{
		switch (i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2)
		{
		case 0:// 1 lei 2j shang
			kind1_node(j, i, -(Num3_p2 - 4) - 2, -(Num3_p2 - 4) - 1, -(Num3_p2 - 4), 0, 1,
				Num4_p2 + N_d0 + 2 - 2, Num4_p2 + N_d0 + 2 - 1, Num4_p2 + N_d0 + 2,
				A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case (Num4_p2 / 2 - 1) :// 1 lei 2j xia
			kind1_node(j, i, -(Num3_p2 - 4), -(Num3_p2 - 4) + 1, -(Num3_p2 - 4) + 2, -1, 0,
			Num4_p2 + N_d0 + 2 - N_c1, Num4_p2 + N_d0 + 2 - N_c1 + 1, Num4_p2 + N_d0 + 2 - N_c1 + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2) <= Node41_p2)// 1dx
			{
				inner_5node(j, i, 1, -(Num3_p2 - 4), Num4_p2 + N_d0 + 2, Coef, Coef_location);
			}
			else if ((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2) >= Node42_p2) // 1dx
			{
				inner_5node(j, i, 1, -(Num3_p2 - 4), Num4_p2 + N_d0 + 2 - N_c1, Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_4node(j, i, -(Num3_p2 - 4), -1, 0, 1, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1_p2 + Num2_p2 + Num3_p2 + Num4_p2 / 2);
		i<(Nx_2 + Nx_3 + Num1_p2 + Num2_p2 + Num3_p2 + Num4_p2); i++)// num 4 :  30 node! part2
	{
		switch (i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2)
		{
		case (Num4_p2 / 2) :// 1 lei 2j shang
			kind1_node(j, i, -(5 + Num4_p2 - 1) - 2, -(5 + Num4_p2 - 1) - 1, -(5 + Num4_p2 - 1), 0, 1,
			(Node_55_p2 - 4 + 1 + N_c1) - 2, (Node_55_p2 - 4 + 1 + N_c1) - 1, (Node_55_p2 - 4 + 1 + N_c1),
			A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case (Num4_p2 - 1) :// 1 lei 2j xia
			kind1_node(j, i, -(5 + Num4_p2 - 1), -(5 + Num4_p2 - 1) + 1, -(5 + Num4_p2 - 1) + 2, -1, 0,
			(Node_55_p2 - 4 + 1), (Node_55_p2 - 4 + 1) + 1, (Node_55_p2 - 4 + 1) + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2) <= Node43_p2)// 1dx
			{
				inner_5node(j, i, 1, -(5 + Num4_p2 - 1), (Node_55_p2 - 4 + 1 + N_c1), Coef, Coef_location);
			}
			else if ((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2) >= Node44_p2) // 1dx
			{
				inner_5node(j, i, 1, -(5 + Num4_p2 - 1), (Node_55_p2 - 4 + 1), Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_4node(j, i, -(5 + Num4_p2 - 1), -1, 0, 1, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (Nx_2 + Nx_3 + Num1_p2 + Num2_p2 + Num3_p2 + Num4_p2);
		i<(Nx_2 + Nx_3 + Num1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + Node_5_p2 + 2); i++)// circle num5 :  80 node!  part1 !
	{
		switch (i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2)
		{
		case 0:
			boundary_node(j, i, 1, -Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2, Num5_p2 + Num6_p2 + Num7_p2 + Num8_p2, Coef, Coef_location);
			break;
		case (N_d0 - 3) ://2 lei shang 1j
			kind2_node(j, i, -Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2, -2, -1, 0, 1, Num5_p2 + Num6_p2 + Num7_p2 + Num8_p2,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (N_d0 + 1) :// 2 lei 2j shang
			kind2_node(j, i, -(Num3_p2 + Num4_p2 + N_d0 - 2), -2, -1, 0, 1, Num5_p2 - N_d0 + 2 + Num6_p2,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node51_p2:// inner 3 node
			inner_3node(j, i, -(Num4_p2 + N_d0 + 2), -1, 0, 1, 1, k*k*dx*dx - 4, Coef, Coef_location);
			break;
		case (Node51_p2 + 1) :// inner 3 node
			inner_3node(j, i, -(Num4_p2 + N_d0 + 2) + N_c1, 0, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node_5_p2 - 3) ://2 lei xia 2j
			kind2_node(j, i, -(Num3_p2 / 2 + 1 + Num4_p2 + Node_5_p2), -1, 0, 1, 2, Num5_p2 - Node_5_p2 + Num6_p2 + Node72_p2,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_5_p2 + 1) ://2 lei xia 1j
			kind2_node(j, i, -(Num1_p2 - Node_1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + Node_5_p2), -1, 0, 1, 2,
			Num5_p2 - Node_5_p2 + Num6_p2 + Num7_p2 + Num8_p2 + Node_9_p2,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2) < (N_d0 - 3))// 4dx!
			{
				inner_5node(j, i, 4, -Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2, Num5_p2 + Num6_p2 + Num7_p2 + Num8_p2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2) > (N_d0 - 3)) && ((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2) < (N_d0 + 1)))// 2dx!
			{
				inner_5node(j, i, 2, -(Num3_p2 + Num4_p2 + N_d0 - 2), Num5_p2 - N_d0 + 2 + Num6_p2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2) > (N_d0 + 1)) && ((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2) < Node51_p2))// 1dx
			{
				inner_5node(j, i, 1, -(Num4_p2 + N_d0 + 2), Num5_p2 - N_d0 - 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2) > (Node51_p2 + 1)) && ((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2) < (Node_5_p2 - 3)))// 1dx
			{
				inner_5node(j, i, 1, -(Num4_p2 + N_d0 + 2) + N_c1, Num5_p2 - N_d0 - 2 - 2, Coef, Coef_location);
			}
			else// 2dx!
			{
				inner_5node(j, i, 2, -(Num3_p2 / 2 + 1 + Num4_p2 + Node_5_p2), Num5_p2 - Node_5_p2 + Num6_p2 + Node72_p2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + Node_5_p2 + 2);
		i<(Nx_2 + Nx_3 + Num1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + Node53_p2 + 3); i++)// circle num5 :  80 node!  part2 !
	{
		switch (i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2)
		{
		default:
			if ((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2) < Node52_p2)// inner 5 node 4dx
			{
				inner_5node(j, i, 4, -(Num1_p2 - Node_1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + Node_5_p2), Num5_p2 - Node_5_p2 + Num6_p2 + Num7_p2 + Num8_p2 + Node_9_p2, Coef, Coef_location);
			}
			else if ((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2) > Node53_p2)// inner 5 node 4dx
			{
				inner_5node(j, i, 4, -(Num1_p2 - Node_1_p2 - 9 + Num2_p2 + Num3_p2 + Num4_p2 + Node53_p2),
					Num5_p2 - Node53_p2 + Num6_p2 + Num7_p2 + Num8_p2 + Node93_p2, Coef, Coef_location);
			}
			else// left  2lei 1j
			{
				kind2_node(j, i, -Nx_3 - Nx_2 - (Num1_p2 - Node_1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + Node_5_p2),
					-(Num1_p2 - Node_1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + Node_5_p2), -1, 0, 1,
					Num5_p2 - (i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2) + Num6_p2 + Node72_p2 + 1 + (i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Node52_p2) * 2,
					B5, B4, B3, B1_1, B3, B2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//**********************************************************************************!!!!!!!!!!!!
	for (i = (Nx_2 + Nx_3 + Num1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + Node53_p2 + 3);
		i<(Nx_2 + Nx_3 + Num1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + Num5_p2); i++)// circle num5 :  80 node!  part3 !
	{
		switch (i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2)
		{
		case (Num5_p2 - 1) :
			boundary_node(j, i, 2, -Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2, Num6_p2 + Num7_p2 + Num8_p2 + Num9_p2, Coef, Coef_location);
			break;
		case (Node53_p2 + 3) ://2 lei shang 1j
			kind2_node(j, i, -(Num1_p2 - Node_1_p2 - 9 + Num2_p2 + Num3_p2 + Num4_p2 + Node53_p2), -2, -1, 0, 1, Num5_p2 - Node53_p2 + Num6_p2 + Num7_p2 + Num8_p2 + Node93_p2,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node53_p2 + 7) :// 2 lei 2j shang
			kind2_node(j, i, -(Num3_p2 / 2 + Num4_p2 + Node53_p2 + 4), -2, -1, 0, 1, Num5_p2 - Node53_p2 - 4 + Num6_p2 + Node73_p2,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node54_p2:// inner 3 node
			inner_3node(j, i, -(Node_55_p2 - 4 + 1) - N_c1, -1, 0, 1, 1, k*k*dx*dx - 4, Coef, Coef_location);
			break;
		case (Node54_p2 + 1) :// inner 3 node
			inner_3node(j, i, -(Node_55_p2 - 4 + 1), 0, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node_55_p2 - 3) ://2 lei xia 2j
			kind2_node(j, i, -(Num4_p2 + Node_55_p2 + 1), -1, 0, 1, 2, Num5_p2 - Node_55_p2 + Num6_p2 + Num7_p2 - 1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_55_p2 + 1) ://2 lei xia 1j
			kind2_node(j, i, -Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2, -1, 0, 1, 2, Num6_p2 + Num7_p2 + Num8_p2 + Num9_p2,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2) < (Node53_p2 + 7))// 2dx!
			{
				inner_5node(j, i, 2, -(Num3_p2 / 2 + Num4_p2 + Node53_p2 + 4), Num5_p2 - Node53_p2 - 4 + Num6_p2 + Node73_p2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2) > (Node53_p2 + 7)) && ((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2) < Node54_p2))// 1dx
			{
				inner_5node(j, i, 1, -(Node_55_p2 - 4 + 1) - N_c1, (Num5_p2 - Node_55_p2 + 4 + Num6_p2 - 1) + 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2) > (Node54_p2 + 1)) && ((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2) < (Node_55_p2 - 3)))// 1dx
			{
				inner_5node(j, i, 1, -(Node_55_p2 - 4 + 1), (Num5_p2 - Node_55_p2 + 4 + Num6_p2 - 1), Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2) > (Node_55_p2 - 3)) && ((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2) < (Node_55_p2 + 1)))// 2dx
			{
				inner_5node(j, i, 2, -(Num4_p2 + Node_55_p2 + 1), Num5_p2 - Node_55_p2 + Num6_p2 + Num7_p2 - 1, Coef, Coef_location);
			}
			else// 4dx!
			{
				inner_5node(j, i, 4, -Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2, Num6_p2 + Num7_p2 + Num8_p2 + Num9_p2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (Nx_2 + Nx_3 + Num1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + Num5_p2);
		i<(Nx_2 + Nx_3 + Num1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + Num5_p2 + Num6_p2); i++)// num 6 :  16 node!
	{
		switch (i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2)
		{
		case 0:// 1 lei 2j shang
			kind1_node(j, i, -(Num5_p2 - N_d0 - 2) - 2, -(Num5_p2 - N_d0 - 2) - 1, -(Num5_p2 - N_d0 - 2), 0, 1, Num6_p2 + 4 - 2,
				Num6_p2 + 4 - 1, Num6_p2 + 4, A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case Node61_p2:// inner 3 node
			inner_3node(j, i, -(Num5_p2 - N_d0 - 2), -1, 0, 1, 1, k*k*dx*dx - 4, Coef, Coef_location);
			break;
		case (Node61_p2 + 1) :// inner 3 node
			inner_3node(j, i, -(Num5_p2 - N_d0 - 2 - 2), 0, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Num6_p2 / 2 - 1) :// 1 lei 2j xia
			kind1_node(j, i, -(Num5_p2 - N_d0 - 2 - 2), -(Num5_p2 - N_d0 - 2 - 2) + 1, -(Num5_p2 - N_d0 - 2 - 2) + 2, -1, 0,
			Num6_p2 + 4 - 2, Num6_p2 + 4 - 2 + 1, Num6_p2 + 4 - 2 + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		case (Num6_p2 / 2) :// 1 lei 2j shang
			kind1_node(j, i, -(Num5_p2 - Node_55_p2 + 4 + Num6_p2 - 1) - 2 - 2, -(Num5_p2 - Node_55_p2 + 4 + Num6_p2 - 1) - 2 - 1,
			-(Num5_p2 - Node_55_p2 + 4 + Num6_p2 - 1) - 2, 0, 1, Num7_p2 - 4 + 2 - 2, Num7_p2 - 4 + 2 - 1, Num7_p2 - 4 + 2,
			A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case Node62_p2:// inner 3 node
			inner_3node(j, i, -(Num5_p2 - Node_55_p2 + 4 + Num6_p2 - 1) - 2, -1, 0, 1, 1, k*k*dx*dx - 4, Coef, Coef_location);
			break;
		case (Node62_p2 + 1) :// inner 3 node
			inner_3node(j, i, -(Num5_p2 - Node_55_p2 + 4 + Num6_p2 - 1), 0, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Num6_p2 - 1) :// 1 lei 2j xia
			kind1_node(j, i, -(Num5_p2 - Node_55_p2 + 4 + Num6_p2 - 1), -(Num5_p2 - Node_55_p2 + 4 + Num6_p2 - 1) + 1,
			-(Num5_p2 - Node_55_p2 + 4 + Num6_p2 - 1) + 2, -1, 0, Num7_p2 - 4, Num7_p2 - 4 + 1, Num7_p2 - 4 + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2) < Node61_p2)// 1dx
			{
				inner_5node(j, i, 1, -(Num5_p2 - N_d0 - 2), Num6_p2 + 4, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2) > (Node61_p2 + 1)) && ((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2) < (Num6_p2 / 2 - 1)))// 2dx
			{
				inner_5node(j, i, 1, -(Num5_p2 - N_d0 - 2 - 2), Num6_p2 + 4 - 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2) > (Num6_p2 / 2)) && ((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2) < Node62_p2))// 2dx
			{
				inner_5node(j, i, 1, -(Num5_p2 - Node_55_p2 + 4 + Num6_p2 - 1) - 2, Num7_p2 - 4 + 2, Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_5node(j, i, 1, -(Num5_p2 - Node_55_p2 + 4 + Num6_p2 - 1), Num7_p2 - 4, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!!
	for (i = (Nx_2 + Nx_3 + Num1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + Num5_p2 + Num6_p2);
		i<(Nx_2 + Nx_3 + Num1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + Num5_p2 + Num6_p2 + Node72_p2 + 1); i++)// num 7 :  39 node! part1
	{
		switch (i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2 - Num6_p2)
		{
		case 0:// 1 lei 1j shang
			kind1_node(j, i, -(Num5_p2 - N_d0 + 2 + Num6_p2) - 2, -(Num5_p2 - N_d0 + 2 + Num6_p2) - 1, -(Num5_p2 - N_d0 + 2 + Num6_p2),
				0, 1, Num7_p2 + Num8_p2 + N_d0 - 2 - 2, Num7_p2 + Num8_p2 + N_d0 - 2 - 1, Num7_p2 + Num8_p2 + N_d0 - 2,
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case 3:// 2 lei 2j shang
			kind2_node(j, i, -(Num5_p2 - N_d0 + 2 + Num6_p2), -2, -1, 0, 1, Num7_p2 + Num8_p2 + N_d0 - 2,
				B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node71_p2://inner 4 node
			inner_4node(j, i, -(Num6_p2 + 4), -1, 0, Num7_p2 - 4, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node71_p2 + 1) ://inner 4 node
			inner_4node(j, i, -(Num6_p2 + 4) + 2, 0, 1, Num7_p2 - 4, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Node72_p2 - 3) :// 2 lei 2j xia
			kind2_node(j, i, -(Num5_p2 - Node_5_p2 + Num6_p2 + Node72_p2), -1, 0, 1, 2, Num7_p2 - Node72_p2 + Num8_p2 + Node_9_p2,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case Node72_p2:// 1 lei 1j xia
			kind1_node(j, i, -(Num5_p2 - Node_5_p2 + Num6_p2 + Node72_p2), -(Num5_p2 - Node_5_p2 + Num6_p2 + Node72_p2) + 1,
				-(Num5_p2 - Node_5_p2 + Num6_p2 + Node72_p2) + 2, -1, 0, Num7_p2 - Node72_p2 + Num8_p2 + Node_9_p2,
				Num7_p2 - Node72_p2 + Num8_p2 + Node_9_p2 + 1, Num7_p2 - Node72_p2 + Num8_p2 + Node_9_p2 + 2,
				A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2 - Num6_p2) < 3)// 2dx
			{
				inner_5node(j, i, 2, -(Num5_p2 - N_d0 + 2 + Num6_p2), Num7_p2 + Num8_p2 + N_d0 - 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2 - Num6_p2) > 3) && ((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2 - Num6_p2) < (Node71_p2)))// 1dx
			{
				inner_5node(j, i, 1, -(Num6_p2 + 4), Num7_p2 - 4, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2 - Num6_p2) > (Node71_p2 + 1)) && ((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2 - Num6_p2) < (Node72_p2 - 3)))// 1dx
			{
				inner_5node(j, i, 1, -(Num6_p2 + 4) + 2, Num7_p2 - 4, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, i, 2, -(Num5_p2 - Node_5_p2 + Num6_p2 + Node72_p2), Num7_p2 - Node72_p2 + Num8_p2 + Node_9_p2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + Num5_p2 + Num6_p2 + Node72_p2 + 1);
		i<(Nx_2 + Nx_3 + Num1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + Num5_p2 + Num6_p2 + Node73_p2); i++)// num 7 :  39 node! part2
	{
		switch (i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2 - Num6_p2 - Node72_p2 - 1)
		{
		case 0:// 3 lei 1j shang
			kind3_node(j, i, -(Num5_p2 - Node52_p2 + Num6_p2 + Node72_p2 + 1) - 1, -(Num5_p2 - Node52_p2 + Num6_p2 + Node72_p2 + 1),
				0, 1, (Num7_p2 - Node72_p2 + Num8_p2 + Node92_p2) - 2, (Num7_p2 - Node72_p2 + Num8_p2 + Node92_p2),
				C4, C3_1, C1_1, C2_1, C4, C3_1, Coef, Coef_location);
			break;
		case 10:// 3 lei 1j xia
			kind3_node(j, i, -(Num5_p2 - Node53_p2 + Num6_p2 + Node73_p2 - 1), -(Num5_p2 - Node53_p2 + Num6_p2 + Node73_p2 - 1) + 1,
				-1, 0, (Num7_p2 - Node72_p2 + Num8_p2 + Node92_p2), (Num7_p2 - Node72_p2 + Num8_p2 + Node92_p2) + 2,
				C3_1, C4, C2_1, C1_1, C3_1, C4, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2 - Num6_p2 - Node72_p2 - 1) % 2 == 0)// 1dx
			{
				inner_5node(j, i, 2, -((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2 - Num6_p2) + Num6_p2 + Num5_p2 - Node52_p2 - (i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2 - Num6_p2 - Node72_p2 - 1) / 2),
					(Num7_p2 - Node72_p2 + Num8_p2 + Node92_p2), Coef, Coef_location);
			}
			else//1 lei  1j left
			{
				kind1_node(j, i, -((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2 - Num6_p2) + Num6_p2 + Num5_p2 - Node52_p2 - (i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2 - Num6_p2 - Node72_p2 - 1) / 2) - (Num1_p2 - Node_1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + Node_5_p2),
					-((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2 - Num6_p2) + Num6_p2 + Num5_p2 - Node52_p2 - (i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2 - Num6_p2 - Node72_p2 - 1) / 2) - (Num1_p2 - Node_1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + Node_5_p2) + 1,
					-((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2 - Num6_p2) + Num6_p2 + Num5_p2 - Node52_p2 - (i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2 - Num6_p2 - Node72_p2 - 1) / 2),
					-((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2 - Num6_p2) + Num6_p2 + Num5_p2 - Node52_p2 - (i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2 - Num6_p2 - Node72_p2 - 1) / 2) + 1,
					-1, 0, 1, (Num7_p2 - Node72_p2 + Num8_p2 + Node92_p2), A2_1, A2_1, A1_1, A1_1, A3_1, A5_1, A3_1, A4_1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + Num5_p2 + Num6_p2 + Node73_p2);
		i<(Nx_2 + Nx_3 + Num1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + Num5_p2 + Num6_p2 + Num7_p2); i++)// num 7 :  39 node! part3
	{
		switch (i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2 - Num6_p2)
		{
		case Node73_p2:// 1 lei 1j shang
			kind1_node(j, i, -(Num5_p2 - Node53_p2 - 4 + Num6_p2 + Node73_p2) - 2, -(Num5_p2 - Node53_p2 - 4 + Num6_p2 + Node73_p2) - 1,
				-(Num5_p2 - Node53_p2 - 4 + Num6_p2 + Node73_p2), 0, 1, Num7_p2 - Node73_p2 + Num8_p2 + Node93_p2 + 4 - 2,
				Num7_p2 - Node73_p2 + Num8_p2 + Node93_p2 + 4 - 1, Num7_p2 - Node73_p2 + Num8_p2 + Node93_p2 + 4,
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case (Node73_p2 + 3) :// 2 lei 2j shang
			kind2_node(j, i, -(Num5_p2 - Node53_p2 - 4 + Num6_p2 + Node73_p2), -2, -1, 0, 1,
			Num7_p2 - Node73_p2 + Num8_p2 + Node93_p2 + 4, B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node74_p2://inner 4 node
			inner_4node(j, i, -(Num7_p2 - 5 + 1) - 2, -1, 0, 5 + Num8_p2 - 1, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node74_p2 + 1) ://inner 4 node
			inner_4node(j, i, -(Num7_p2 - 5 + 1), 0, 1, 5 + Num8_p2 - 1, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Num7_p2 - 4) :// 2 lei 2j xia
			kind2_node(j, i, -(Num5_p2 - Node_55_p2 + Num6_p2 + Num7_p2 - 1), -1, 0, 1, 2, 1 + Num8_p2 + Node_99_p2,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num7_p2 - 1) :// 1 lei 1j xia
			kind1_node(j, i, -(Num5_p2 - Node_55_p2 + Num6_p2 + Num7_p2 - 1), -(Num5_p2 - Node_55_p2 + Num6_p2 + Num7_p2 - 1) + 1,
			-(Num5_p2 - Node_55_p2 + Num6_p2 + Num7_p2 - 1) + 2, -1, 0,
			1 + Num8_p2 + Node_99_p2, 1 + Num8_p2 + Node_99_p2 + 1, 1 + Num8_p2 + Node_99_p2 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2 - Num6_p2) < (Node73_p2 + 3))// 2dx
			{
				inner_5node(j, i, 2, -(Num5_p2 - Node53_p2 - 4 + Num6_p2 + Node73_p2), Num7_p2 - Node73_p2 + Num8_p2 + Node93_p2 + 4, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2 - Num6_p2) > (Node73_p2 + 3)) && ((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2 - Num6_p2) < (Node74_p2)))// 1dx
			{
				inner_5node(j, i, 1, -(Num7_p2 - 5 + 1) - 2, 5 + Num8_p2 - 1, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2 - Num6_p2) > (Node74_p2 + 1)) && ((i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2 - Num6_p2) < (Num7_p2 - 4)))// 1dx
			{
				inner_5node(j, i, 1, -(Num7_p2 - 5 + 1), 5 + Num8_p2 - 1, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, i, 2, -(Num5_p2 - Node_55_p2 + Num6_p2 + Num7_p2 - 1), 1 + Num8_p2 + Node_99_p2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!!!!!!
	for (i = (Nx_2 + Nx_3 + Num1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + Num5_p2 + Num6_p2 + Num7_p2);
		i<(Nx_2 + Nx_3 + Num1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + Num5_p2 + Num6_p2 + Num7_p2 + Num8_p2); i++)// num 8 :  12 node! 
	{
		switch (i - Nx_2 - Nx_3 - Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2 - Num5_p2 - Num6_p2 - Num7_p2)
		{
		case 0:// 1 lei shang 2j
			kind1_node(j, i, -(Num7_p2 - 4) - 2, -(Num7_p2 - 4) - 1, -(Num7_p2 - 4), 0, 1, Num8_p2 + N_d0 + 2 - 2, Num8_p2 + N_d0 + 2 - 1,
				Num8_p2 + N_d0 + 2, A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case 1:
			inner_5node(j, i, 1, -(Num7_p2 - 4), Num8_p2 + N_d0 + 2, Coef, Coef_location);
			break;
		case 2://inner 4 node
			inner_4node(j, i, -(Num7_p2 - 4), -1, 0, Num8_p2 + N_d0 + 2, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case 3://inner 4 node
			inner_4node(j, i, -(Num7_p2 - 4), 0, 1, Num8_p2 + N_d0 + 2, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case 4:
			inner_5node(j, i, 1, -(Num7_p2 - 4), Num8_p2 + N_d0 + 2, Coef, Coef_location);
			break;
		case (Num8_p2 / 2 - 1) :// 1 lei xia 2j
			kind1_node(j, i, -(Num7_p2 - 4), -(Num7_p2 - 4) + 1, -(Num7_p2 - 4) + 2, -1, 0, Num8_p2 + N_d0 + 2, Num8_p2 + N_d0 + 2 + 1,
			Num8_p2 + N_d0 + 2 + 2, A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		case (Num8_p2 / 2) :// 1 lei shang 2j
			kind1_node(j, i, -(5 + Num8_p2 - 1) - 2, -(5 + Num8_p2 - 1) - 1, -(5 + Num8_p2 - 1), 0, 1,
			Node_99_p2 - 4 + 1 - 2, Node_99_p2 - 4 + 1 - 1, Node_99_p2 - 4 + 1,
			A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case 7:
			inner_5node(j, i, 1, -(5 + Num8_p2 - 1), Node_99_p2 - 4 + 1, Coef, Coef_location);
			break;
		case 8://inner 4 node
			inner_4node(j, i, -(5 + Num8_p2 - 1), -1, 0, Node_99_p2 - 4 + 1, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case 9://inner 4 node
			inner_4node(j, i, -(5 + Num8_p2 - 1), 0, 1, Node_99_p2 - 4 + 1, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case 10:
			inner_5node(j, i, 1, -(5 + Num8_p2 - 1), Node_99_p2 - 4 + 1, Coef, Coef_location);
			break;
		case (Num8_p2 - 1) :// 1 lei xia 2j
			kind1_node(j, i, -(5 + Num8_p2 - 1), -(5 + Num8_p2 - 1) + 1, -(5 + Num8_p2 - 1) + 2, -1, 0,
			Node_99_p2 - 4 + 1, Node_99_p2 - 4 + 1 + 1, Node_99_p2 - 4 + 1 + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (M_post2_q1 - Num9_p2);
		i<(M_post2_q1 - Num9_p2 + Node_9_p2 + 2); i++)// circle num9 :  79 node!  part1 !
	{
		switch (i - M_post2_q1 + Num9_p2)
		{
		case 0:
			boundary_node(j, i, 1, -Num5_p2 - Num6_p2 - Num7_p2 - Num8_p2, Num10_p2 + Num11_p2 + Num12_p2 + Num9_p2, Coef, Coef_location);
			break;
		case (N_d0 - 3) ://2 lei shang 1j
			kind2_node(j, i, -Num5_p2 - Num6_p2 - Num7_p2 - Num8_p2, -2, -1, 0, 1, Num9_p2 + Num10_p2 + Num11_p2 + Num12_p2,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (N_d0 + 1) :// 2 lei 2j shang
			kind2_node(j, i, -(Num7_p2 + Num8_p2 + N_d0 - 2), -2, -1, 0, 1, Num9_p2 - N_d0 + 2 + Num10_p2,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node91_p2:// inner 4 node
			inner_4node(j, i, -(Num8_p2 + N_d0 + 2), -1, 0, Num9_p2 - N_d0 - 2, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node91_p2 + 1) :// inner 4 node
			inner_4node(j, i, -(Num8_p2 + N_d0 + 2), 0, 1, Num9_p2 - N_d0 - 2, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Node_9_p2 - 3) ://2 lei xia 2j
			kind2_node(j, i, -(Num7_p2 - Node72_p2 + Num8_p2 + Node_9_p2), -1, 0, 1, 2, Num9_p2 - Node_9_p2 + Num10_p2 + Node112_p2,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_9_p2 + 1) ://2 lei xia 1j
			kind2_node(j, i, -(Num5_p2 - Node_5_p2 + Num6_p2 + Num7_p2 + Num8_p2 + Node_9_p2), -1, 0, 1, 2,
			Num9_p2 - Node_9_p2 + Num10_p2 + Num11_p2 + Num12_p2 + Node_13_p2, B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q1 + Num9_p2) < (N_d0 - 3))// 4dx!
			{
				inner_5node(j, i, 4, -Num5_p2 - Num6_p2 - Num7_p2 - Num8_p2, Num9_p2 + Num10_p2 + Num11_p2 + Num12_p2, Coef, Coef_location);
			}
			else if (((i - M_post2_q1 + Num9_p2) > (N_d0 - 3)) && ((i - M_post2_q1 + Num9_p2) < (N_d0 + 1)))// 2dx!
			{
				inner_5node(j, i, 2, -(Num7_p2 + Num8_p2 + N_d0 - 2), Num9_p2 - N_d0 + 2 + Num10_p2, Coef, Coef_location);
			}
			else if (((i - M_post2_q1 + Num9_p2) > (N_d0 + 1)) && ((i - M_post2_q1 + Num9_p2) < Node91_p2))// 1dx
			{
				inner_5node(j, i, 1, -(Num8_p2 + N_d0 + 2), Num9_p2 - N_d0 - 2, Coef, Coef_location);
			}
			else if (((i - M_post2_q1 + Num9_p2) > (Node91_p2 + 1)) && ((i - M_post2_q1 + Num9_p2) < (Node_9_p2 - 3)))// 1dx
			{
				inner_5node(j, i, 1, -(Num8_p2 + N_d0 + 2), Num9_p2 - N_d0 - 2, Coef, Coef_location);
			}
			else// 2dx!
			{
				inner_5node(j, i, 2, -(Num7_p2 - Node72_p2 + Num8_p2 + Node_9_p2), Num9_p2 - Node_9_p2 + Num10_p2 + Node112_p2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post2_q1 - Num9_p2 + Node_9_p2 + 2); i<(M_post2_q1 - Num9_p2 + Node93_p2 + 3); i++)// circle num9 :  79 node!  part2 !
	{
		switch (i - M_post2_q1 + Num9_p2)
		{
		case (Node92_p2 - 1) ://2 lei shang
			kind2_node(j, i, -(Num5_p2 - Node_5_p2 + Num6_p2 + Num7_p2 + Num8_p2 + Node_9_p2), -2, -1, 0, 1,
			Num9_p2 - Node_9_p2 + Num10_p2 + Num11_p2 + Num12_p2 + Node_13_p2, B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case Node92_p2://3 lei shang
			kind3_node(j, i, -(Num5_p2 - Node_5_p2 + Num6_p2 + Num7_p2 + Num8_p2 + Node_9_p2) - 1, -(Num5_p2 - Node_5_p2 + Num6_p2 + Num7_p2 + Num8_p2 + Node_9_p2),
				-1, 0, 1, Num9_p2 - Node92_p2 + Num10_p2 + Node112_p2 + 1, C4, C4, C3_1, C1_1, C3_1, C2_1, Coef, Coef_location);
			break;
		case Node93_p2://3 lei xia
			kind3_node(j, i, -(Num5_p2 - Node53_p2 + Num6_p2 + Num7_p2 + Num8_p2 + Node93_p2),
				-(Num5_p2 - Node53_p2 + Num6_p2 + Num7_p2 + Num8_p2 + Node93_p2) + 1, -1, 0, 1, Num9_p2 - Node92_p2 + Num10_p2 + Node112_p2 + 1,
				C4, C4, C3_1, C1_1, C3_1, C2_1, Coef, Coef_location);
			break;
		case (Node93_p2 + 1) ://2 lei xia
			kind2_node(j, i, -(Num5_p2 - Node53_p2 + Num6_p2 + Num7_p2 + Num8_p2 + Node93_p2), -1, 0, 1, 2,
			Num9_p2 - Node93_p2 + Num10_p2 + Num11_p2 + Num12_p2 + Node133_p2, B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q1 + Num9_p2) < (Node92_p2 - 1))// 4dx
			{
				inner_5node(j, i, 4, -(Num5_p2 - Node_5_p2 + Num6_p2 + Num7_p2 + Num8_p2 + Node_9_p2),
					Num9_p2 - Node_9_p2 + Num10_p2 + Num11_p2 + Num12_p2 + Node_13_p2, Coef, Coef_location);
			}
			else if ((i - M_post2_q1 + Num9_p2) > (Node93_p2 + 1))// 2dx!
			{
				inner_5node(j, i, 4, -(Num5_p2 - Node53_p2 + Num6_p2 + Num7_p2 + Num8_p2 + Node93_p2),
					Num9_p2 - Node93_p2 + Num10_p2 + Num11_p2 + Num12_p2 + Node133_p2, Coef, Coef_location);
			}
			else if (((i - M_post2_q1 + Num9_p2) > Node92_p2) && ((i - M_post2_q1 + Num9_p2) < Node93_p2))
			{
				inner_5node(j, i, 2, -(Num7_p2 - Node72_p2 + Num8_p2 + Node92_p2),
					Num9_p2 - Node92_p2 + Num10_p2 + Node112_p2 + 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post2_q1 - Num9_p2 + Node93_p2 + 3); i<(M_post2_q1); i++)// circle num9 :  79 node!  part3 !
	{
		switch (i - M_post2_q1 + Num9_p2)
		{
		case (Num9_p2 - 1) :
			boundary_node(j, i, 2, -Num6_p2 - Num7_p2 - Num8_p2 - Num9_p2, Num10_p2 + Num11_p2 + Num12_p2 + Num13_p2, Coef, Coef_location);
			break;
		case (Node93_p2 + 3) ://2 lei shang 1j
			kind2_node(j, i, -(Num5_p2 - Node53_p2 + Num6_p2 + Num7_p2 + Num8_p2 + Node93_p2), -2, -1, 0, 1,
			Num9_p2 - Node93_p2 + Num10_p2 + Num11_p2 + Num12_p2 + Node133_p2, B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node94_p2 - 3) :// 2 lei 2j shang
			kind2_node(j, i, -(Num7_p2 - Node73_p2 + Num8_p2 + Node93_p2 + 4), -2, -1, 0, 1,
			Num9_p2 - Node93_p2 - 4 + Num10_p2 + Node113_p2, B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node94_p2:// inner 4 node
			inner_4node(j, i, -(Node_99_p2 - 4 + 1), -1, 0, Num9_p2 - Node_99_p2 + 4 + Num10_p2 - 1, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node94_p2 + 1) :// inner 4 node
			inner_4node(j, i, -(Node_99_p2 - 4 + 1), 0, 1, Num9_p2 - Node_99_p2 + 4 + Num10_p2 - 1, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Node_99_p2 - 3) ://2 lei xia 2j
			kind2_node(j, i, -(Num8_p2 + Node_99_p2 + 1), -1, 0, 1, 2, Num9_p2 - Node_99_p2 + Num10_p2 + Num11_p2 - 1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_99_p2 + 1) ://2 lei xia 1j
			kind2_node(j, i, -Num6_p2 - Num7_p2 - Num8_p2 - Num9_p2, -1, 0, 1, 2, Num10_p2 + Num11_p2 + Num12_p2 + Num13_p2,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q1 + Num9_p2) < (Node94_p2 - 3))// 2dx!
			{
				inner_5node(j, i, 2, -(Num7_p2 - Node73_p2 + Num8_p2 + Node93_p2 + 4), Num9_p2 - Node93_p2 - 4 + Num10_p2 + Node113_p2, Coef, Coef_location);
			}
			else if (((i - M_post2_q1 + Num9_p2) > (Node94_p2 - 3)) && ((i - M_post2_q1 + Num9_p2) < Node94_p2))// 1dx
			{
				inner_5node(j, i, 1, -(Node_99_p2 - 4 + 1), Num9_p2 - Node_99_p2 + 4 + Num10_p2 - 1, Coef, Coef_location);
			}
			else if (((i - M_post2_q1 + Num9_p2) > (Node94_p2 + 1)) && ((i - M_post2_q1 + Num9_p2) < (Node_99_p2 - 3)))// 1dx
			{
				inner_5node(j, i, 1, -(Node_99_p2 - 4 + 1), Num9_p2 - Node_99_p2 + 4 + Num10_p2 - 1, Coef, Coef_location);
			}
			else if (((i - M_post2_q1 + Num9_p2) > (Node_99_p2 - 3)) && ((i - M_post2_q1 + Num9_p2) < (Node_99_p2 + 1)))// 2dx
			{
				inner_5node(j, i, 2, -(Num8_p2 + Node_99_p2 + 1), Num9_p2 - Node_99_p2 + Num10_p2 + Num11_p2 - 1, Coef, Coef_location);
			}
			else// 4dx!
			{
				inner_5node(j, i, 4, -Num6_p2 - Num7_p2 - Num8_p2 - Num9_p2, Num10_p2 + Num11_p2 + Num12_p2 + Num13_p2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (M_post2_q1); i<(M_post2_q1 + Num10_p2); i++)// num 10 :  12 node! 
	{
		switch (i - M_post2_q1)
		{
		case 0:// 1 lei shang 2j
			kind1_node(j, i, -(Num9_p2 - N_d0 - 2) - 2, -(Num9_p2 - N_d0 - 2) - 1, -(Num9_p2 - N_d0 - 2), 0, 1, Num10_p2 + 4 - 2,
				Num10_p2 + 4 - 1, Num10_p2 + 4, A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case 1:
			inner_5node(j, i, 1, -(Num9_p2 - N_d0 - 2), Num10_p2 + 4, Coef, Coef_location);
			break;
		case 2://inner 4 node
			inner_4node(j, i, -(Num9_p2 - N_d0 - 2), -1, 0, Num10_p2 + 4, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case 3://inner 4 node
			inner_4node(j, i, -(Num9_p2 - N_d0 - 2), 0, 1, Num10_p2 + 4, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case 4:
			inner_5node(j, i, 1, -(Num9_p2 - N_d0 - 2), Num10_p2 + 4, Coef, Coef_location);
			break;
		case (Num10_p2 / 2 - 1) :// 1 lei xia 2j
			kind1_node(j, i, -(Num9_p2 - N_d0 - 2), -(Num9_p2 - N_d0 - 2) + 1, -(Num9_p2 - N_d0 - 2) + 2, -1, 0, Num10_p2 + 4,
			Num10_p2 + 4 + 1, Num10_p2 + 4 + 2, A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		case (Num10_p2 / 2) :// 1 lei shang 2j
			kind1_node(j, i, -(Num9_p2 - Node_99_p2 + 4 + Num10_p2 - 1) - 2, -(Num9_p2 - Node_99_p2 + 4 + Num10_p2 - 1) - 1,
			-(Num9_p2 - Node_99_p2 + 4 + Num10_p2 - 1), 0, 1, Num11_p2 - 4 - 2, Num11_p2 - 4 - 1, Num11_p2 - 4,
			A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case 7:
			inner_5node(j, i, 1, -(Num9_p2 - Node_99_p2 + 4 + Num10_p2 - 1), Num11_p2 - 4, Coef, Coef_location);
			break;
		case 8://inner 4 node
			inner_4node(j, i, -(Num9_p2 - Node_99_p2 + 4 + Num10_p2 - 1), -1, 0, Num11_p2 - 4, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case 9://inner 4 node
			inner_4node(j, i, -(Num9_p2 - Node_99_p2 + 4 + Num10_p2 - 1), 0, 1, Num11_p2 - 4, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case 10:
			inner_5node(j, i, 1, -(Num9_p2 - Node_99_p2 + 4 + Num10_p2 - 1), Num11_p2 - 4, Coef, Coef_location);
			break;
		case (Num10_p2 - 1) :// 1 lei xia 2j
			kind1_node(j, i, -(Num9_p2 - Node_99_p2 + 4 + Num10_p2 - 1), -(Num9_p2 - Node_99_p2 + 4 + Num10_p2 - 1) + 1,
			-(Num9_p2 - Node_99_p2 + 4 + Num10_p2 - 1) + 2, -1, 0, Num11_p2 - 4, Num11_p2 - 4 + 1, Num11_p2 - 4 + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (M_post2_q1 + Num10_p2); i<(M_post2_q1 + Num10_p2 + Node112_p2 + 1); i++)// num 11 :  41 node! part1
	{
		switch (i - M_post2_q1 - Num10_p2)
		{
		case 0:// 1 lei 1j shang
			kind1_node(j, i, -(Num9_p2 - N_d0 + 2 + Num10_p2) - 2, -(Num9_p2 - N_d0 + 2 + Num10_p2) - 1, -(Num9_p2 - N_d0 + 2 + Num10_p2), 0, 1,
				Num11_p2 + Num12_p2 + N_d0 - 2 - 2, Num11_p2 + Num12_p2 + N_d0 - 2 - 1, Num11_p2 + Num12_p2 + N_d0 - 2,
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case 3:// 2 lei 2j shang
			kind2_node(j, i, -(Num9_p2 - N_d0 + 2 + Num10_p2), -2, -1, 0, 1, Num11_p2 + Num12_p2 + N_d0 - 2,
				B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node111_p2://inner 4 node
			inner_4node(j, i, -(Num10_p2 + 4), -1, 0, Num11_p2 - 4, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node111_p2 + 1) ://inner 4 node
			inner_4node(j, i, -(Num10_p2 + 4), 0, 1, Num11_p2 - 4 + 2, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Node112_p2 - 3) :// 2 lei 2j xia
			kind2_node(j, i, -(Num9_p2 - Node_9_p2 + Num10_p2 + Node112_p2), -1, 0, 1, 2, Num11_p2 - Node112_p2 + Num12_p2 + Node_13_p2,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case Node112_p2:// 1 lei 1j xia
			kind1_node(j, i, -(Num9_p2 - Node_9_p2 + Num10_p2 + Node112_p2), -(Num9_p2 - Node_9_p2 + Num10_p2 + Node112_p2) + 1,
				-(Num9_p2 - Node_9_p2 + Num10_p2 + Node112_p2) + 2, -1, 0, Num11_p2 - Node112_p2 + Num12_p2 + Node_13_p2,
				Num11_p2 - Node112_p2 + Num12_p2 + Node_13_p2 + 1, Num11_p2 - Node112_p2 + Num12_p2 + Node_13_p2 + 2,
				A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q1 - Num10_p2) < 3)// 2dx
			{
				inner_5node(j, i, 2, -(Num9_p2 - N_d0 + 2 + Num10_p2), Num11_p2 + Num12_p2 + N_d0 - 2, Coef, Coef_location);
			}
			else if (((i - M_post2_q1 - Num10_p2) > 3) && ((i - M_post2_q1 - Num10_p2) < (Node111_p2)))// 1dx
			{
				inner_5node(j, i, 1, -(Num10_p2 + 4), Num11_p2 - 4, Coef, Coef_location);
			}
			else if (((i - M_post2_q1 - Num10_p2) > (Node111_p2 + 1)) && ((i - M_post2_q1 - Num10_p2) < (Node112_p2 - 3)))// 1dx
			{
				inner_5node(j, i, 1, -(Num10_p2 + 4), Num11_p2 - 4 + 2, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, i, 2, -(Num9_p2 - Node_9_p2 + Num10_p2 + Node112_p2), Num11_p2 - Node112_p2 + Num12_p2 + Node_13_p2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post2_q1 + Num10_p2 + Node112_p2 + 1); i<(M_post2_q1 + Num10_p2 + Node113_p2); i++)// num 11 :  41 node! part2
	{
		switch (i - M_post2_q1 - Num10_p2)
		{
		case (Node112_p2 + 1) :// 1lei shang 1j
			kind1_node(j, i, -(Num9_p2 - Node92_p2 + Num10_p2 + Node112_p2 + 1) - 2, -(Num9_p2 - Node92_p2 + Num10_p2 + Node112_p2 + 1) - 1,
			-(Num9_p2 - Node92_p2 + Num10_p2 + Node112_p2 + 1), 0, 1, (Num11_p2 - Node112_p2 - 1 + Num12_p2 + Node132_p2) - 2,
			(Num11_p2 - Node112_p2 - 1 + Num12_p2 + Node132_p2) - 1, (Num11_p2 - Node112_p2 - 1 + Num12_p2 + Node132_p2),
			A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case (Node113_p2 - 1) :// 1lei xia 1j
			kind1_node(j, i, -(Num9_p2 - Node92_p2 + Num10_p2 + Node112_p2 + 1), -(Num9_p2 - Node92_p2 + Num10_p2 + Node112_p2 + 1) + 1,
			-(Num9_p2 - Node92_p2 + Num10_p2 + Node112_p2 + 1) + 2, -1, 0, (Num11_p2 - Node113_p2 + 1 + Num12_p2 + Node133_p2),
			(Num11_p2 - Node113_p2 + 1 + Num12_p2 + Node133_p2) + 1, (Num11_p2 - Node113_p2 + 1 + Num12_p2 + Node133_p2) + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q1 - Num10_p2) < (Node112_p2 + 4))// 4dx
			{
				inner_5node(j, i, 2, -(Num9_p2 - Node92_p2 + Num10_p2 + Node112_p2 + 1), (Num11_p2 - Node112_p2 - 1 + Num12_p2 + Node132_p2), Coef, Coef_location);
			}
			else if ((i - M_post2_q1 - Num10_p2) > (Node113_p2 - 4))// 4dx
			{
				inner_5node(j, i, 2, -(Num9_p2 - Node92_p2 + Num10_p2 + Node112_p2 + 1), (Num11_p2 - Node113_p2 + 1 + Num12_p2 + Node133_p2), Coef, Coef_location);
			}
			else// 2lei zuo 2j
			{
				kind2_node(j, i, -(Num9_p2 - Node92_p2 + Num10_p2 + Node112_p2 + 1) - (Num7_p2 - Node72_p2 + Num8_p2 + Node92_p2),
					-(Num9_p2 - Node92_p2 + Num10_p2 + Node112_p2 + 1), -1, 0, 1,
					Num11_p2 - (i - M_post2_q1 - Num10_p2) + Node122_p2 + 1 + (i - M_post2_q1 - Num10_p2 - Node112_p2 - 4) * 2,
					B5, B4, B3, B1_2, B3, B2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post2_q1 + Num10_p2 + Node113_p2); i<(M_post2_q1 + Num10_p2 + Num11_p2); i++)// num 11 :  41 node! part3
	{
		switch (i - M_post2_q1 - Num10_p2)
		{
		case Node113_p2:// 1 lei 1j shang
			kind1_node(j, i, -(Num9_p2 - Node93_p2 - 4 + Num10_p2 + Node113_p2) - 2, -(Num9_p2 - Node93_p2 - 4 + Num10_p2 + Node113_p2) - 1,
				-(Num9_p2 - Node93_p2 - 4 + Num10_p2 + Node113_p2), 0, 1, Num11_p2 - Node113_p2 + Num12_p2 + Node133_p2 + 4 - 2,
				Num11_p2 - Node113_p2 + Num12_p2 + Node133_p2 + 4 - 1, Num11_p2 - Node113_p2 + Num12_p2 + Node133_p2 + 4,
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case (Node113_p2 + 3) :// 2 lei 2j shang
			kind2_node(j, i, -(Num9_p2 - Node93_p2 - 4 + Num10_p2 + Node113_p2), -2, -1, 0, 1, Num11_p2 - Node113_p2 + Num12_p2 + Node133_p2 + 4,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node114_p2://inner 4 node
			inner_4node(j, i, -(Num11_p2 - 5 + 1), -1, 0, 5 + Num12_p2 - 1 - 2, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node114_p2 + 1) ://inner 4 node
			inner_4node(j, i, -(Num11_p2 - 5 + 1), 0, 1, 5 + Num12_p2 - 1, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Num11_p2 - 4) :// 2 lei 2j xia
			kind2_node(j, i, -(Num9_p2 - Node_99_p2 + Num10_p2 + Num11_p2 - 1), -1, 0, 1, 2, 1 + Num12_p2 + Node_1313_p2,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num11_p2 - 1) :// 1 lei 1j xia
			kind1_node(j, i, -(Num9_p2 - Node_99_p2 + Num10_p2 + Num11_p2 - 1), -(Num9_p2 - Node_99_p2 + Num10_p2 + Num11_p2 - 1) + 1,
			-(Num9_p2 - Node_99_p2 + Num10_p2 + Num11_p2 - 1) + 2, -1, 0,
			1 + Num12_p2 + Node_1313_p2, 1 + Num12_p2 + Node_1313_p2 + 1, 1 + Num12_p2 + Node_1313_p2 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q1 - Num10_p2) < (Node113_p2 + 3))// 2dx
			{
				inner_5node(j, i, 2, -(Num9_p2 - Node93_p2 - 4 + Num10_p2 + Node113_p2), Num11_p2 - Node113_p2 + Num12_p2 + Node133_p2 + 4, Coef, Coef_location);
			}
			else if (((i - M_post2_q1 - Num10_p2) > (Node113_p2 + 3)) && ((i - M_post2_q1 - Num10_p2) < (Node114_p2)))// 1dx
			{
				inner_5node(j, i, 1, -(Num11_p2 - 5 + 1), 5 + Num12_p2 - 1 - 2, Coef, Coef_location);
			}
			else if (((i - M_post2_q1 - Num10_p2) > (Node114_p2 + 1)) && ((i - M_post2_q1 - Num10_p2) < (Num11_p2 - 4)))// 1dx
			{
				inner_5node(j, i, 1, -(Num11_p2 - 5 + 1), 5 + Num12_p2 - 1, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, i, 2, -(Num9_p2 - Node_99_p2 + Num10_p2 + Num11_p2 - 1), 1 + Num12_p2 + Node_1313_p2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (M_post2_q1 + Num10_p2 + Num11_p2);
		i<(M_post2_q1 + Num10_p2 + Num11_p2 + Node122_p2 + 1); i++)// num 12 :  29 node!  part1
	{
		switch (i - M_post2_q1 - Num10_p2 - Num11_p2)
		{
		case 0:// 1 lei 2j shang
			kind1_node(j, i, -(Num11_p2 - 4) - 2, -(Num11_p2 - 4) - 1, -(Num11_p2 - 4), 0, 1, Num12_p2 + N_d0 + 2 - 2,
				Num12_p2 + N_d0 + 2 - 1, Num12_p2 + N_d0 + 2, A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case Node121_p2:// inner 3 node
			inner_3node(j, i, -1, 0, (Num12_p2 + N_d0 + 2), 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node121_p2 + 1) :// inner 3 node
			inner_3node(j, i, 0, 1, (Num12_p2 + N_d0 + 2) + 2, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case Node122_p2:// 1 lei 2j xia
			kind1_node(j, i, (-(Num11_p2 - 4) - 2), (-(Num11_p2 - 4) - 2) + 1, (-(Num11_p2 - 4) - 2) + 2, -1, 0,
				(Num12_p2 + N_d0 + 2) + 2, (Num12_p2 + N_d0 + 2) + 2 + 1, (Num12_p2 + N_d0 + 2) + 2 + 2,
				A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q1 - Num10_p2 - Num11_p2) < Node121_p2)// 1dx
			{
				inner_5node(j, i, 1, -(Num11_p2 - 4), Num12_p2 + N_d0 + 2, Coef, Coef_location);
			}
			else // 2dx
			{
				inner_5node(j, i, 1, -(Num11_p2 - 4) - 2, (Num12_p2 + N_d0 + 2) + 2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post2_q1 + Num10_p2 + Num11_p2 + Node122_p2 + 1); i<(M_post2_q1 + Num10_p2 + Num11_p2 + Node123_p2); i++)// num 12 :  29 node! part3
	{
		switch (i - M_post2_q1 - Num10_p2 - Num11_p2 - Node122_p2 - 1)
		{
		case (Node122_p2 + 1 - Node122_p2 - 1) :// 3 lei 2j shang !!! 1
			kind3_node(j, i, -(Num11_p2 - Node112_p2 - 3 + Node122_p2) - 1, -(Num11_p2 - Node112_p2 - 3 + Node122_p2), 0, 1,
			(Num12_p2 - Node122_p2 + Node132_p2 + 3) - 2, (Num12_p2 - Node122_p2 + Node132_p2 + 3),
			C4, C3_2, C1_2, C2_2, C4, C3_2, Coef, Coef_location);
			break;
		case (Node123_p2 - 1 - Node122_p2 - 1) :// 3 lei 2j xia !!! 5
			kind3_node(j, i, -(Num11_p2 - Node113_p2 + 3 + Node123_p2), -(Num11_p2 - Node113_p2 + 3 + Node123_p2) + 1, -1, 0,
			(Num12_p2 - Node122_p2 + Node132_p2 + 3), (Num12_p2 - Node122_p2 + Node132_p2 + 3) + 2,
			C3_2, C4, C2_2, C1_2, C3_2, C4, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q1 - Num10_p2 - Num11_p2 - Node122_p2 - 1) % 2 == 0)// 1dx
			{
				inner_5node(j, i, 1, -((i - M_post2_q1 - Num10_p2 - Num11_p2) + Num11_p2 - Node112_p2 - 4 - (i - M_post2_q1 - Num10_p2 - Num11_p2 - Node122_p2 - 1) / 2),
					(Num12_p2 - Node122_p2 + Node132_p2 + 3), Coef, Coef_location);
			}
			else// 1 lei 2j zuo
			{
				kind1_node(j, i, -((i - M_post2_q1 - Num10_p2 - Num11_p2) + Num11_p2 - Node112_p2 - 4 - (i - M_post2_q1 - Num10_p2 - Num11_p2 - Node122_p2 - 1) / 2) - (Num9_p2 - Node92_p2 + Num10_p2 + Node112_p2 + 1),
					-((i - M_post2_q1 - Num10_p2 - Num11_p2) + Num11_p2 - Node112_p2 - 4 - (i - M_post2_q1 - Num10_p2 - Num11_p2 - Node122_p2 - 1) / 2) - (Num9_p2 - Node92_p2 + Num10_p2 + Node112_p2 + 1) + 1,
					-((i - M_post2_q1 - Num10_p2 - Num11_p2) + Num11_p2 - Node112_p2 - 4 - (i - M_post2_q1 - Num10_p2 - Num11_p2 - Node122_p2 - 1) / 2),
					-((i - M_post2_q1 - Num10_p2 - Num11_p2) + Num11_p2 - Node112_p2 - 4 - (i - M_post2_q1 - Num10_p2 - Num11_p2 - Node122_p2 - 1) / 2) + 1,
					-1, 0, 1, (Num12_p2 - Node122_p2 + Node132_p2 + 3), A2_2, A2_2, A1_2, A1_2, A3_2, A5_2, A3_2, A4_2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post2_q1 + Num10_p2 + Num11_p2 + Node123_p2); i<(M_post2_q1 + Num10_p2 + Num11_p2 + Num12_p2); i++)// num 12 :  29 node! part3
	{
		switch (i - M_post2_q1 - Num10_p2 - Num11_p2)
		{
		case Node123_p2:// 1 lei 2j shang
			kind1_node(j, i, -(5 + Num12_p2 - 1) + 2 - 2, -(5 + Num12_p2 - 1) + 2 - 1, -(5 + Num12_p2 - 1) + 2, 0, 1,
				(Node_1313_p2 - 4 + 1) - 2 - 2, (Node_1313_p2 - 4 + 1) - 2 - 1, (Node_1313_p2 - 4 + 1) - 2,
				A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case Node124_p2:// inner 3 node
			inner_3node(j, i, -1, 0, (Node_1313_p2 - 4 + 1) - 2, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node124_p2 + 1) :// inner 3 node
			inner_3node(j, i, 0, 1, (Node_1313_p2 - 4 + 1), k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Num12_p2 - 1) :// 1 lei 2j xia
			kind1_node(j, i, -(5 + Num12_p2 - 1), -(5 + Num12_p2 - 1) + 1, -(5 + Num12_p2 - 1) + 2, -1, 0,
			(Node_1313_p2 - 4 + 1), (Node_1313_p2 - 4 + 1) + 1, (Node_1313_p2 - 4 + 1) + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q1 - Num10_p2 - Num11_p2) < Node124_p2)// 2dx
			{
				inner_5node(j, i, 1, -(5 + Num12_p2 - 1) + 2, (Node_1313_p2 - 4 + 1) - 2, Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_5node(j, i, 1, -(5 + Num12_p2 - 1), (Node_1313_p2 - 4 + 1), Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (M_post2_q1 + Num10_p2 + Num11_p2 + Num12_p2);
		i<(M_post2_q1 + Num10_p2 + Num11_p2 + Num12_p2 + Node_13_p2 + 2); i++)// circle num13 :  95 node!  part1 !
	{
		switch (i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2)
		{
		case 0:
			boundary_node(j, i, 1, -Num9_p2 - Num10_p2 - Num11_p2 - Num12_p2, Num13_p2 + Num14_p2 + Num15_p2 + Num16_p2, Coef, Coef_location);
			break;
		case (N_d0 - 3) ://2 lei shang 1j
			kind2_node(j, i, -Num9_p2 - Num10_p2 - Num11_p2 - Num12_p2, -2, -1, 0, 1, Num13_p2 + Num14_p2 + Num15_p2 + Num16_p2,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (N_d0 + 1) :// 2 lei 2j shang
			kind2_node(j, i, -(Num11_p2 + Num12_p2 + N_d0 - 2), -2, -1, 0, 1, Num13_p2 - N_d0 + 2 + Num14_p2,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node131_p2:// inner 3 node
			inner_3node(j, i, -1, 0, (Num13_p2 - N_d0 - 2), 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node131_p2 + 1) :// inner 3 node
			inner_3node(j, i, 0, 1, (Num13_p2 - N_d0 - 2) + N_c1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Node_13_p2 - 3) ://2 lei xia 2j
			kind2_node(j, i, -(Num11_p2 - Node112_p2 + Num12_p2 + Node_13_p2), -1, 0, 1, 2, Num13_p2 - Node_13_p2 + Num14_p2 + Node151_p2,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_13_p2 + 1) ://2 lei xia 1j
			kind2_node(j, i, -(Num9_p2 - Node_9_p2 + Num10_p2 + Num11_p2 + Num12_p2 + Node_13_p2), -1, 0, 1, 2,
			Num13_p2 - Node_13_p2 + Num14_p2 + Num15_p2 + Num16_p2 + Node_17_p2, B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2) < (N_d0 - 3))// 4dx!
			{
				inner_5node(j, i, 4, -Num9_p2 - Num10_p2 - Num11_p2 - Num12_p2, Num13_p2 + Num14_p2 + Num15_p2 + Num16_p2, Coef, Coef_location);
			}
			else if (((i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2) > (N_d0 - 3)) && ((i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2) < (N_d0 + 1)))// 2dx!
			{
				inner_5node(j, i, 2, -(Num11_p2 + Num12_p2 + N_d0 - 2), Num13_p2 - N_d0 + 2 + Num14_p2, Coef, Coef_location);
			}
			else if (((i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2) > (N_d0 + 1)) && ((i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2) < Node131_p2))// 1dx
			{
				inner_5node(j, i, 1, -(Num12_p2 + N_d0 + 2), (Num13_p2 - N_d0 - 2), Coef, Coef_location);
			}
			else if (((i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2) > (Node131_p2 + 1)) && ((i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2) < (Node_13_p2 - 3)))// 1dx
			{
				inner_5node(j, i, 1, -(Num12_p2 + N_d0 + 2) - 2, (Num13_p2 - N_d0 - 2) + N_c1, Coef, Coef_location);
			}
			else// 2dx!
			{
				inner_5node(j, i, 2, -(Num11_p2 - Node112_p2 + Num12_p2 + Node_13_p2), Num13_p2 - Node_13_p2 + Num14_p2 + Node151_p2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post2_q1 + Num10_p2 + Num11_p2 + Num12_p2 + Node_13_p2 + 2);
		i<(M_post2_q1 + Num10_p2 + Num11_p2 + Num12_p2 + Node133_p2 + 3); i++)// circle num13 :  95 node!  part2 !
	{
		switch (i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2)
		{
		case (Node_13_p2 + 2) :// 4dx
			inner_5node(j, i, 4, -(Num9_p2 - Node_9_p2 + Num10_p2 + Num11_p2 + Num12_p2 + Node_13_p2), Num13_p2 - Node_13_p2 + Num14_p2 + Num15_p2 + Num16_p2 + Node_17_p2, Coef, Coef_location);
			break;
		case (Node_13_p2 + 3) :// 2lei shang 1j
			kind2_node(j, i, -(Num9_p2 - Node_9_p2 + Num10_p2 + Num11_p2 + Num12_p2 + Node_13_p2), -2, -1, 0, 1, Num13_p2 - Node_13_p2 + Num14_p2 + Num15_p2 + Num16_p2 + Node_17_p2,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node132_p2 + 2) :// 2lei shang 2j
			kind2_node(j, i, -(Num11_p2 - Node112_p2 - 1 + Num12_p2 + Node132_p2), -2, -1, 0, 1, Num13_p2 - Node132_p2 + Num14_p2 + Node151_p2 + 1,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case (Node132_p2 + 3) :// 3lei zuo 1j
			kind3_node(j, i, -(Num11_p2 - Node112_p2 - 1 + Num12_p2 + Node132_p2) - 1, -(Num11_p2 - Node112_p2 - 1 + Num12_p2 + Node132_p2),
			-1, 0, 1, Num13_p2 - Node132_p2 - 3 + Node143_p2 - 4, C4, C4, C3_2, C1_2, C3_2, C2_2, Coef, Coef_location);
			break;
		case (Node133_p2 - 3) :// 3lei zuo 1j
			kind3_node(j, i, -(Num11_p2 - Node113_p2 + 1 + Num12_p2 + Node133_p2), -(Num11_p2 - Node113_p2 + 1 + Num12_p2 + Node133_p2) + 1,
			-1, 0, 1, Num13_p2 - Node132_p2 - 3 + Node143_p2 - 4, C4, C4, C3_2, C1_2, C3_2, C2_2, Coef, Coef_location);
			break;
		case (Node133_p2 - 2) :// 2lei xia 2j
			kind2_node(j, i, -(Num11_p2 - Node113_p2 + 1 + Num12_p2 + Node133_p2), -1, 0, 1, 2, Num13_p2 - Node133_p2 + Num14_p2 + Node153_p2 - 1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node133_p2 + 1) :// 2lei xia 1j
			kind2_node(j, i, -(Num9_p2 - Node93_p2 + Num10_p2 + Num11_p2 + Num12_p2 + Node133_p2), -1, 0, 1, 2,
			(Num13_p2 - Node133_p2 + Num14_p2 + Num15_p2 + Num16_p2 + Node174_p2), B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node133_p2 + 2) :// 4dx
			inner_5node(j, i, 4, -(Num9_p2 - Node93_p2 + Num10_p2 + Num11_p2 + Num12_p2 + Node133_p2),
			(Num13_p2 - Node133_p2 + Num14_p2 + Num15_p2 + Num16_p2 + Node174_p2), Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2) < (Node132_p2 + 2))  // 2dx
			{
				inner_5node(j, i, 2, -(Num11_p2 - Node112_p2 - 1 + Num12_p2 + Node132_p2), Num13_p2 - Node132_p2 + Num14_p2 + Node151_p2 + 1, Coef, Coef_location);
			}
			else if ((i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2) > (Node133_p2 - 2))  // 2dx
			{
				inner_5node(j, i, 2, -(Num11_p2 - Node113_p2 + 1 + Num12_p2 + Node133_p2), Num13_p2 - Node133_p2 + Num14_p2 + Node153_p2 - 1, Coef, Coef_location);
			}
			else  // 1dx
			{
				inner_5node(j, i, 1, -(Num12_p2 - Node122_p2 - 1 + Node132_p2 + 4), Num13_p2 - Node132_p2 - 3 + Node143_p2 - 4, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post2_q1 + Num10_p2 + Num11_p2 + Num12_p2 + Node133_p2 + 3);
		i<(M_post2_q1 + Num10_p2 + Num11_p2 + Num12_p2 + Num13_p2); i++)// circle num13 :  95 node!  part3 !
	{
		switch (i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2)
		{
		case (Num13_p2 - 1) :
			boundary_node(j, i, 2, -Num10_p2 - Num11_p2 - Num12_p2 - Num13_p2, Num14_p2 + Num15_p2 + Num16_p2 + Num17_p2, Coef, Coef_location);
			break;
		case (Node133_p2 + 3) ://2 lei shang 1j
			kind2_node(j, i, -(Num9_p2 - Node93_p2 + Num10_p2 + Num11_p2 + Num12_p2 + Node133_p2), -2, -1, 0, 1,
			(Num13_p2 - Node133_p2 + Num14_p2 + Num15_p2 + Num16_p2 + Node174_p2), B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node133_p2 + 7) :// 2 lei 2j shang
			kind2_node(j, i, -(Num11_p2 - Node113_p2 + Num12_p2 + Node133_p2 + 4), -2, -1, 0, 1, Num13_p2 - Node133_p2 - 4 + Num14_p2 + Node153_p2,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node134_p2:// inner 3 node
			inner_3node(j, i, -1, 0, (Num13_p2 - Node_1313_p2 + 4 + Num14_p2 - 1) - N_c1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node134_p2 + 1) :// inner 3 node
			inner_3node(j, i, 0, 1, (Num13_p2 - Node_1313_p2 + 4 + Num14_p2 - 1), k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Node_1313_p2 - 3) ://2 lei xia 2j
			kind2_node(j, i, -(Num12_p2 + Node_1313_p2 + 1), -1, 0, 1, 2, Num13_p2 - Node_1313_p2 + Num14_p2 + Num15_p2 - 1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_1313_p2 + 1) ://2 lei xia 1j
			kind2_node(j, i, -Num10_p2 - Num11_p2 - Num12_p2 - Num13_p2, -1, 0, 1, 2, Num14_p2 + Num15_p2 + Num16_p2 + Num17_p2,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2) < (Node133_p2 + 7))// 2dx!
			{
				inner_5node(j, i, 2, -(Num11_p2 - Node113_p2 + Num12_p2 + Node133_p2 + 4), Num13_p2 - Node133_p2 - 4 + Num14_p2 + Node153_p2, Coef, Coef_location);
			}
			else if (((i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2) > (Node133_p2 + 7)) && ((i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2) < Node134_p2))// 1dx
			{
				inner_5node(j, i, 1, -(Node_1313_p2 - 4 + 1) + 2, (Num13_p2 - Node_1313_p2 + 4 + Num14_p2 - 1) - N_c1, Coef, Coef_location);
			}
			else if (((i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2) > (Node134_p2 + 1)) && ((i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2) < (Node_1313_p2 - 3)))// 1dx
			{
				inner_5node(j, i, 1, -(Node_1313_p2 - 4 + 1), (Num13_p2 - Node_1313_p2 + 4 + Num14_p2 - 1), Coef, Coef_location);
			}
			else if (((i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2) > (Node_1313_p2 - 3)) && ((i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2) < (Node_1313_p2 + 1)))// 2dx
			{
				inner_5node(j, i, 2, -(Num12_p2 + Node_1313_p2 + 1), Num13_p2 - Node_1313_p2 + Num14_p2 + Num15_p2 - 1, Coef, Coef_location);
			}
			else// 4dx!
			{
				inner_5node(j, i, 4, -Num10_p2 - Num11_p2 - Num12_p2 - Num13_p2, Num14_p2 + Num15_p2 + Num16_p2 + Num17_p2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (M_post2_q1 + Num10_p2 + Num11_p2 + Num12_p2 + Num13_p2);
		i<(M_post2_q1 + Num10_p2 + Num11_p2 + Num12_p2 + Num13_p2 + Node142_p2 + 5); i++)// num 14 :  45 node! part1
	{
		switch (i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2 - Num13_p2)
		{
		case 0:// 1 lei 2j shang
			kind1_node(j, i, -(Num13_p2 - N_d0 - 2) - 2, -(Num13_p2 - N_d0 - 2) - 1, -(Num13_p2 - N_d0 - 2), 0, 1, Num14_p2 + 4 - 2,
				Num14_p2 + 4 - 1, Num14_p2 + 4, A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case (Node142_p2 + 4) :// 1 lei 2j xia
			kind1_node(j, i, -(Num13_p2 - N_d0 - 2) - N_c1, -(Num13_p2 - N_d0 - 2) - N_c1 + 1, -(Num13_p2 - N_d0 - 2) - N_c1 + 2, -1, 0,
			Num14_p2 + 4, Num14_p2 + 4 + 1, Num14_p2 + 4 + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2 - Num13_p2) <= Node141_p2)// 1dx
			{
				inner_5node(j, i, 1, -(Num13_p2 - N_d0 - 2), Num14_p2 + 4, Coef, Coef_location);
			}
			else if ((i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2 - Num13_p2) >= Node142_p2) // 1dx
			{
				inner_5node(j, i, 1, -(Num13_p2 - N_d0 - 2) - N_c1, Num14_p2 + 4, Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_4node(j, i, -1, 0, 1, Num14_p2 + 4, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post2_q1 + Num10_p2 + Num11_p2 + Num12_p2 + Num13_p2 + Node142_p2 + 5);
		i<(M_post2_q1 + Num10_p2 + Num11_p2 + Num12_p2 + Num13_p2 + Node145_p2 - 4); i++)// num 14 :  45 node! part2
	{
		switch (i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2 - Num13_p2)
		{
		case (Node142_p2 + 5) :// 1 lei 2j shang
			kind1_node(j, i, -(Num13_p2 - Node132_p2 - 3 + Node142_p2 + 5) - 2, -(Num13_p2 - Node132_p2 - 3 + Node142_p2 + 5) - 1, -(Num13_p2 - Node132_p2 - 3 + Node142_p2 + 5),
			0, 1, (Num14_p2 - Node142_p2 - 5 + Node151_p2 + 4) - 2, (Num14_p2 - Node142_p2 - 5 + Node151_p2 + 4) - 1, (Num14_p2 - Node142_p2 - 5 + Node151_p2 + 4),
			A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case (Node145_p2 - 5) :// 1 lei 2j xia
			kind1_node(j, i, -(Num13_p2 - Node132_p2 - 3 + Node142_p2 + 5), -(Num13_p2 - Node132_p2 - 3 + Node142_p2 + 5), -(Num13_p2 - Node132_p2 - 3 + Node142_p2 + 5) + 2,
			-1, 0, (Num14_p2 - Node142_p2 - 5 + Node151_p2 + 4) - N_c1, (Num14_p2 - Node142_p2 - 5 + Node151_p2 + 4) - N_c1 + 1, (Num14_p2 - Node142_p2 - 5 + Node151_p2 + 4) - N_c1 + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2 - Num13_p2) <= Node143_p2)// 1dx
			{
				inner_5node(j, i, 1, -(Num13_p2 - Node132_p2 - 3 + Node142_p2 + 5), (Num14_p2 - Node142_p2 - 5 + Node151_p2 + 4), Coef, Coef_location);
			}
			else if ((i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2 - Num13_p2) >= Node144_p2) // 1dx
			{
				inner_5node(j, i, 1, -(Num13_p2 - Node132_p2 - 3 + Node142_p2 + 5), (Num14_p2 - Node142_p2 - 5 + Node151_p2 + 4) - N_c1, Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_4node(j, i, -(Num13_p2 - Node132_p2 - 3 + Node142_p2 + 5), -1, 0, 1, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post2_q1 + Num10_p2 + Num11_p2 + Num12_p2 + Num13_p2 + Node145_p2 - 4);
		i<(M_post2_q1 + Num10_p2 + Num11_p2 + Num12_p2 + Num13_p2 + Num14_p2); i++)// num14 :  45 node! part3
	{
		switch (i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2 - Num13_p2)
		{
		case (Node145_p2 - 4) :// 1 lei 2j shang
			kind1_node(j, i, -(Num13_p2 - Node_1313_p2 + 4 + Num14_p2 - 1) + N_c1 - 2, -(Num13_p2 - Node_1313_p2 + 4 + Num14_p2 - 1) + N_c1 - 1,
			-(Num13_p2 - Node_1313_p2 + 4 + Num14_p2 - 1) + N_c1, 0, 1,
			(Num15_p2 - 5 + 1) - 2, (Num15_p2 - 5 + 1) - 1, (Num15_p2 - 5 + 1),
			A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case (Num14_p2 - 1) :// 1 lei 2j xia
			kind1_node(j, i, -(Num13_p2 - Node_1313_p2 + 4 + Num14_p2 - 1), -(Num13_p2 - Node_1313_p2 + 4 + Num14_p2 - 1) + 1, -(Num13_p2 - Node_1313_p2 + 4 + Num14_p2 - 1) + 2,
			-1, 0, (Num15_p2 - 5 + 1), (Num15_p2 - 5 + 1) + 1, (Num15_p2 - 5 + 1) + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2 - Num13_p2) <= Node145_p2)// 1dx
			{
				inner_5node(j, i, 1, -(Num13_p2 - Node_1313_p2 + 4 + Num14_p2 - 1) + N_c1, (Num15_p2 - 5 + 1), Coef, Coef_location);
			}
			else if ((i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2 - Num13_p2) >= Node146_p2) // 1dx
			{
				inner_5node(j, i, 1, -(Num13_p2 - Node_1313_p2 + 4 + Num14_p2 - 1), (Num15_p2 - 5 + 1), Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_4node(j, i, -1, 0, 1, (Num15_p2 - 5 + 1), 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (M_post2_q1 + Num10_p2 + Num11_p2 + Num12_p2 + Num13_p2 + Num14_p2);
		i<(M_post2_q1 + Num10_p2 + Num11_p2 + Num12_p2 + Num13_p2 + Num14_p2 + Node151_p2 + 1); i++)// num 15 :  62 node! part1
	{
		switch (i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2 - Num13_p2 - Num14_p2)
		{
		case 0:// 1 lei 1j shang
			kind1_node(j, i, -(Num13_p2 - N_d0 + 2 + Num14_p2) - 2, -(Num13_p2 - N_d0 + 2 + Num14_p2) - 1, -(Num13_p2 - N_d0 + 2 + Num14_p2),
				0, 1, Num15_p2 + Num16_p2 + N_d0 - 2 - 2, Num15_p2 + Num16_p2 + N_d0 - 2 - 1, Num15_p2 + Num16_p2 + N_d0 - 2,
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case 3:// 2 lei 2j shang
			kind2_node(j, i, -(Num13_p2 - N_d0 + 2 + Num14_p2), -2, -1, 0, 1, Num15_p2 + Num16_p2 + N_d0 - 2,
				B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case 4:// 3 lei 2j shang
			kind3_node(j, i, -(Num14_p2 + 4), -1, 0, 1, Num15_p2 + Num16_p2 + N_d0 - 2 - 1, Num15_p2 + Num16_p2 + N_d0 - 2,
				C2_2, C3_2, C1_2, C3_2, C4, C4, Coef, Coef_location);
			break;
		case (Node151_p2 - 4) :// 3 lei 2j xia
			kind3_node(j, i, -(Num14_p2 + 4), -1, 0, 1, Num15_p2 - Node151_p2 + Num16_p2 + Node_17_p2, Num15_p2 - Node151_p2 + Num16_p2 + Node_17_p2 + 1,
			C2_2, C3_2, C1_2, C3_2, C4, C4, Coef, Coef_location);
			break;
		case (Node151_p2 - 3) :// 2 lei 2j xia
			kind2_node(j, i, -(Num13_p2 - Node_13_p2 + Num14_p2 + Node151_p2), -1, 0, 1, 2, Num15_p2 - Node151_p2 + Num16_p2 + Node_17_p2,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case Node151_p2:// 1 lei 1j xia
			kind1_node(j, i, -(Num13_p2 - Node_13_p2 + Num14_p2 + Node151_p2), -(Num13_p2 - Node_13_p2 + Num14_p2 + Node151_p2) + 1, -(Num13_p2 - Node_13_p2 + Num14_p2 + Node151_p2) + 2,
				-1, 0, Num15_p2 - Node151_p2 + Num16_p2 + Node_17_p2, Num15_p2 - Node151_p2 + Num16_p2 + Node_17_p2 + 1, Num15_p2 - Node151_p2 + Num16_p2 + Node_17_p2 + 2,
				A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2 - Num13_p2 - Num14_p2) < 3)// 2dx
			{
				inner_5node(j, i, 2, -(Num13_p2 - N_d0 + 2 + Num14_p2), Num15_p2 + Num16_p2 + N_d0 - 2, Coef, Coef_location);
			}
			else if (((i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2 - Num13_p2 - Num14_p2) > 4) && ((i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2 - Num13_p2 - Num14_p2) < (Node151_p2 - 4)))// 1dx
			{
				inner_5node(j, i, 1, -(Num14_p2 + 4), Num15_p2 - 5, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, i, 2, -(Num13_p2 - Node_13_p2 + Num14_p2 + Node151_p2), Num15_p2 - Node151_p2 + Num16_p2 + Node_17_p2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post2_q1 + Num10_p2 + Num11_p2 + Num12_p2 + Num13_p2 + Num14_p2 + Node151_p2 + 1);
		i<(M_post2_q1 + Num10_p2 + Num11_p2 + Num12_p2 + Num13_p2 + Num14_p2 + Node153_p2); i++)// num 15 :  62 node! part2
	{
		switch (i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2 - Num13_p2 - Num14_p2)
		{
		case (Node151_p2 + 1) :// 1 lei 1j shang
			kind1_node(j, i, -(Num13_p2 - Node132_p2 + Num14_p2 + Node151_p2 + 1) - 2, -(Num13_p2 - Node132_p2 + Num14_p2 + Node151_p2 + 1) - 1,
			-(Num13_p2 - Node132_p2 + Num14_p2 + Node151_p2 + 1), 0, 1, (Num15_p2 - Node151_p2 - 1 + Num16_p2 + Node173_p2) - 2,
			(Num15_p2 - Node151_p2 - 1 + Num16_p2 + Node173_p2) - 1, (Num15_p2 - Node151_p2 - 1 + Num16_p2 + Node173_p2),
			A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case (Node151_p2 + 2) ://  2dx
			inner_5node(j, i, 2, -(Num13_p2 - Node132_p2 + Num14_p2 + Node151_p2 + 1), (Num15_p2 - Node151_p2 - 1 + Num16_p2 + Node173_p2), Coef, Coef_location);
			break;
		case (Node151_p2 + 3) :// 2 lei 2j shang
			kind2_node(j, i, -(Num13_p2 - Node132_p2 + Num14_p2 + Node151_p2 + 1), -2, -1, 0, 1, (Num15_p2 - Node151_p2 - 1 + Num16_p2 + Node173_p2),
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node152_p2:// inner 3 node
			inner_3node(j, i, -(Num14_p2 - Node143_p2 + Node152_p2), -1, 0, 1, 1, k*k*dx*dx - 4, Coef, Coef_location);
			break;
		case (Node152_p2 + 1) :// inner 3 node
			inner_3node(j, i, -(Num14_p2 - Node143_p2 + Node152_p2) + N_c1, 0, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node153_p2 - 3) ://2 lei xia 2j
			kind2_node(j, i, -(Num13_p2 - Node133_p2 + Num14_p2 + Node153_p2 - 1), -1, 0, 1, 2, (Num15_p2 - Node153_p2 + 1 + Num16_p2 + Node174_p2),
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node153_p2 - 2) :// 2dx
			inner_5node(j, i, 2, -(Num13_p2 - Node133_p2 + Num14_p2 + Node153_p2 - 1), (Num15_p2 - Node153_p2 + 1 + Num16_p2 + Node174_p2), Coef, Coef_location);
			break;
		case (Node153_p2 - 1) :// 1 lei 1j xia
			kind1_node(j, i, -(Num13_p2 - Node133_p2 + Num14_p2 + Node153_p2 - 1), -(Num13_p2 - Node133_p2 + Num14_p2 + Node153_p2 - 1) + 1,
			-(Num13_p2 - Node133_p2 + Num14_p2 + Node153_p2 - 1) + 2, -1, 0, (Num15_p2 - Node153_p2 + 1 + Num16_p2 + Node174_p2),
			(Num15_p2 - Node153_p2 + 1 + Num16_p2 + Node174_p2) + 1, (Num15_p2 - Node153_p2 + 1 + Num16_p2 + Node174_p2) + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2 - Num13_p2 - Num14_p2) < Node152_p2)// 1dx
			{
				inner_5node(j, i, 1, -(Num14_p2 - Node143_p2 + Node152_p2), Num15_p2 - Node152_p2 + 1 + Node162_p2, Coef, Coef_location);
			}
			else if ((i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2 - Num13_p2 - Num14_p2) > (Node152_p2 + 1))// 1dx
			{
				inner_5node(j, i, 1, -(Num14_p2 - Node143_p2 + Node152_p2) + N_c1, Num15_p2 - Node152_p2 + 1 + Node162_p2 - 2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post2_q1 + Num10_p2 + Num11_p2 + Num12_p2 + Num13_p2 + Num14_p2 + Node153_p2);
		i<(M_post2_q1 + Num10_p2 + Num11_p2 + Num12_p2 + Num13_p2 + Num14_p2 + Num15_p2); i++)// num 15 :  62 node! part3
	{
		switch (i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2 - Num13_p2 - Num14_p2)
		{
		case Node153_p2:// 1 lei 1j shang
			kind1_node(j, i, -(Num13_p2 - Node133_p2 - 4 + Num14_p2 + Node153_p2) - 2, -(Num13_p2 - Node133_p2 - 4 + Num14_p2 + Node153_p2) - 1,
				-(Num13_p2 - Node133_p2 - 4 + Num14_p2 + Node153_p2), 0, 1, (Num15_p2 - Node153_p2 + Num16_p2 + Node174_p2 + 4) - 2,
				(Num15_p2 - Node153_p2 + Num16_p2 + Node174_p2 + 4) - 1, (Num15_p2 - Node153_p2 + Num16_p2 + Node174_p2 + 4),
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case (Node153_p2 + 3) :// 2 lei 2j shang
			kind2_node(j, i, -(Num13_p2 - Node133_p2 - 4 + Num14_p2 + Node153_p2), -2, -1, 0, 1, (Num15_p2 - Node153_p2 + Num16_p2 + Node174_p2 + 4),
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case (Node153_p2 + 4) :// 3 lei 2j shang
			kind3_node(j, i, -(Num15_p2 - 5 + 1), -1, 0, 1,
			(Num15_p2 - Node153_p2 + Num16_p2 + Node174_p2 + 4) - 1, (Num15_p2 - Node153_p2 + Num16_p2 + Node174_p2 + 4),
			C2_2, C3_2, C1_2, C3_2, C4, C4, Coef, Coef_location);
			break;
		case (Num15_p2 - 5) :// 3 lei 2j xia
			kind3_node(j, i, -(Num15_p2 - 5 + 1), -1, 0, 1, (1 + Num16_p2 + Node_1717_p2), (1 + Num16_p2 + Node_1717_p2) + 1,
			C2_2, C3_2, C1_2, C3_2, C4, C4, Coef, Coef_location);
			break;
		case (Num15_p2 - 4) :// 2 lei 2j xia
			kind2_node(j, i, -(Num13_p2 - Node_1313_p2 + Num14_p2 + Num15_p2 - 1), -1, 0, 1, 2, (1 + Num16_p2 + Node_1717_p2),
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num15_p2 - 1) :// 1 lei 1j xia
			kind1_node(j, i, -(Num13_p2 - Node_1313_p2 + Num14_p2 + Num15_p2 - 1), -(Num13_p2 - Node_1313_p2 + Num14_p2 + Num15_p2 - 1) + 1,
			-(Num13_p2 - Node_1313_p2 + Num14_p2 + Num15_p2 - 1) + 2, -1, 0, (1 + Num16_p2 + Node_1717_p2),
			(1 + Num16_p2 + Node_1717_p2) + 1, (1 + Num16_p2 + Node_1717_p2) + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2 - Num13_p2 - Num14_p2) < (Node153_p2 + 3))// 2dx
			{
				inner_5node(j, i, 2, -(Num13_p2 - Node133_p2 - 4 + Num14_p2 + Node153_p2), (Num15_p2 - Node153_p2 + Num16_p2 + Node174_p2 + 4), Coef, Coef_location);
			}
			else if (((i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2 - Num13_p2 - Num14_p2) > (Node153_p2 + 4)) && ((i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2 - Num13_p2 - Num14_p2) < (Num15_p2 - 4)))// 1dx
			{
				inner_5node(j, i, 1, -(Num15_p2 - 5 + 1), 6 + Num16_p2 - 1, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, i, 2, -(Num13_p2 - Node_1313_p2 + Num14_p2 + Num15_p2 - 1), (1 + Num16_p2 + Node_1717_p2), Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (M_post2_q1 + Num10_p2 + Num11_p2 + Num12_p2 + Num13_p2 + Num14_p2 + Num15_p2);
		i<(M_post2_q1 + Num10_p2 + Num11_p2 + Num12_p2 + Num13_p2 + Num14_p2 + Num15_p2 + Node161_p2 + 1); i++)// num 16 :  34 node! part1
	{
		switch (i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2 - Num13_p2 - Num14_p2 - Num15_p2)
		{
		case 0:// 3 lei 2j shang
			kind3_node(j, i, -(Num15_p2 - 3), -(Num15_p2 - 3) + 2, 0, 1, Num16_p2 + Node171_p2 - 1, Num16_p2 + Node171_p2,
				C4, C3_2, C1_2, C2_2, C4, C3_2, Coef, Coef_location);
			break;
		case Node161_p2://node2: 3 lei 2j xia
			kind3_node(j, i, -(Num15_p2 - 5), -(Num15_p2 - 5) + 2, -1, 0, (Num16_p2 - Node161_p2 + Node172_p2),
				(Num16_p2 - Node161_p2 + Node172_p2) + 1, C3_2, C4, C2_2, C1_2, C3_2, C4, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2 - Num13_p2 - Num14_p2 - Num15_p2) % 2 == 0)// 1dx
			{
				inner_5node(j, i, 1, -(Num15_p2 - 5),
					(Num16_p2 - (i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2 - Num13_p2 - Num14_p2 - Num15_p2) + Node171_p2 + (i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2 - Num13_p2 - Num14_p2 - Num15_p2) / 2), Coef, Coef_location);
			}
			else//1 lei  2j left
			{
				kind1_node(j, i, -(Num15_p2 - 5), -1, 0, 1,
					(Num16_p2 - (i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2 - Num13_p2 - Num14_p2 - Num15_p2) + Node171_p2 + (i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2 - Num13_p2 - Num14_p2 - Num15_p2) / 2),
					(Num16_p2 - (i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2 - Num13_p2 - Num14_p2 - Num15_p2) + Node171_p2 + (i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2 - Num13_p2 - Num14_p2 - Num15_p2) / 2) + 1,
					(Num16_p2 - (i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2 - Num13_p2 - Num14_p2 - Num15_p2) + Node171_p2 + (i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2 - Num13_p2 - Num14_p2 - Num15_p2) / 2) + (Num17_p2 - Node171_p2 + 4 + Num18_p2),
					(Num16_p2 - (i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2 - Num13_p2 - Num14_p2 - Num15_p2) + Node171_p2 + (i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2 - Num13_p2 - Num14_p2 - Num15_p2) / 2) + 1 + (Num17_p2 - Node171_p2 + 4 + Num18_p2),
					A4_2, A3_2, A5_2, A3_2, A1_2, A1_2, A2_2, A2_2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post2_q1 + Num10_p2 + Num11_p2 + Num12_p2 + Num13_p2 + Num14_p2 + Num15_p2 + Node161_p2 + 1);
		i<(M_post2_q1 + Num10_p2 + Num11_p2 + Num12_p2 + Num13_p2 + Num14_p2 + Num15_p2 + Node163_p2); i++)// num 16 :  34 node! part2
	{
		switch (i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2 - Num13_p2 - Num14_p2 - Num15_p2)
		{
		case (Node161_p2 + 1) :// 1 lei 2j shang
			kind1_node(j, i, -(Num15_p2 - Node151_p2 - 3 + Node161_p2) - 2, -(Num15_p2 - Node151_p2 - 3 + Node161_p2) - 1,
			-(Num15_p2 - Node151_p2 - 3 + Node161_p2), 0, 1, (Num16_p2 - Node161_p2 + Node173_p2 + 2) - 2, (Num16_p2 - Node161_p2 + Node173_p2 + 2) - 1,
			(Num16_p2 - Node161_p2 + Node173_p2 + 2), A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case Node162_p2:// inner 3 node
			inner_3node(j, i, -(Num15_p2 - Node151_p2 - 3 + Node161_p2), -1, 0, 1, 1, k*k*dx*dx - 4, Coef, Coef_location);
			break;
		case (Node162_p2 + 1) :// inner 3 node
			inner_3node(j, i, -(Num15_p2 - Node151_p2 - 3 + Node161_p2) + 2, 0, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node163_p2 - 1) :// 1 lei 2j xia
			kind1_node(j, i, -(Num15_p2 - Node151_p2 - 3 + Node161_p2) + 2, -(Num15_p2 - Node151_p2 - 3 + Node161_p2) + 2 + 1,
			-(Num15_p2 - Node151_p2 - 3 + Node161_p2) + 2 + 2, -1, 0, (Num16_p2 - Node161_p2 + Node173_p2 + 2) - 2, (Num16_p2 - Node161_p2 + Node173_p2 + 2) - 2 + 1,
			(Num16_p2 - Node161_p2 + Node173_p2 + 2) - 2 + 2, A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2 - Num13_p2 - Num14_p2 - Num15_p2) < Node162_p2)// 1dx
			{
				inner_5node(j, i, 1, -(Num15_p2 - Node151_p2 - 3 + Node161_p2), (Num16_p2 - Node161_p2 + Node173_p2 + 2), Coef, Coef_location);
			}
			else // 2dx
			{
				inner_5node(j, i, 1, -(Num15_p2 - Node151_p2 - 3 + Node161_p2) + 2, (Num16_p2 - Node161_p2 + Node173_p2 + 2) - 2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post2_q1 + Num10_p2 + Num11_p2 + Num12_p2 + Num13_p2 + Num14_p2 + Num15_p2 + Node163_p2);
		i<(M_post2_q1 + Num10_p2 + Num11_p2 + Num12_p2 + Num13_p2 + Num14_p2 + Num15_p2 + Num16_p2); i++)// num 16 :  34 node! part3
	{
		switch (i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2 - Num13_p2 - Num14_p2 - Num15_p2 - Node163_p2)
		{
		case 0:// 3 lei 2j shang
			kind3_node(j, i, -(6 + Num16_p2 - 1) - 2, -(6 + Num16_p2 - 1), 0, 1, Num16_p2 - Node163_p2 + Node175_p2 - 1, Num16_p2 - Node163_p2 + Node175_p2,
				C4, C3_2, C1_2, C2_2, C4, C3_2, Coef, Coef_location);
			break;
		case Node161_p2:// 3 lei 2j xia
			kind3_node(j, i, -(6 + Num16_p2 - 1), -(6 + Num16_p2 - 1) + 2, -1, 0, 1 + Node176_p2, 1 + Node176_p2 + 1,
				C3_2, C4, C2_2, C1_2, C3_2, C4, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2 - Num13_p2 - Num14_p2 - Num15_p2 - Node163_p2) % 2 == 0)// 1dx
			{
				inner_5node(j, i, 1, -(6 + Num16_p2 - 1),
					Num16_p2 - (i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2 - Num13_p2 - Num14_p2 - Num15_p2) + Node175_p2 + (i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2 - Num13_p2 - Num14_p2 - Num15_p2 - Node163_p2) / 2,
					Coef, Coef_location);
			}
			else//1 lei 2j left
			{
				kind1_node(j, i, -(6 + Num16_p2 - 1), -1, 0, 1,
					Num16_p2 - (i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2 - Num13_p2 - Num14_p2 - Num15_p2) + Node175_p2 + (i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2 - Num13_p2 - Num14_p2 - Num15_p2 - Node163_p2) / 2,
					Num16_p2 - (i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2 - Num13_p2 - Num14_p2 - Num15_p2) + Node175_p2 + (i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2 - Num13_p2 - Num14_p2 - Num15_p2 - Node163_p2) / 2 + 1,
					Num16_p2 - (i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2 - Num13_p2 - Num14_p2 - Num15_p2) + Node175_p2 + (i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2 - Num13_p2 - Num14_p2 - Num15_p2 - Node163_p2) / 2 + (Num17_p2 - Node_1717_p2 + Num18_p2 + Num19_p2 - 1),
					Num16_p2 - (i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2 - Num13_p2 - Num14_p2 - Num15_p2) + Node175_p2 + (i - M_post2_q1 - Num10_p2 - Num11_p2 - Num12_p2 - Num13_p2 - Num14_p2 - Num15_p2 - Node163_p2) / 2 + 1 + (Num17_p2 - Node_1717_p2 + Num18_p2 + Num19_p2 - 1),
					A4_2, A3_2, A5_2, A3_2, A1_2, A1_2, A2_2, A2_2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (M_post2_q2 - Num19_p2 - Num18_p2 - Num17_p2);
		i<(M_post2_q2 - Num19_p2 - Num18_p2 - Num17_p2 + Node_17_p2 + 2); i++)// circle num17 :  80 node!  part1 !
	{
		switch (i - M_post2_q2 + Num19_p2 + Num18_p2 + Num17_p2)
		{
		case 0:
			boundary_node(j, i, 1, -Num13_p2 - Num14_p2 - Num15_p2 - Num16_p2, Num17_p2 + Num18_p2 + Num19_p2 + Num18_p2, Coef, Coef_location);
			break;
		case (N_d0 - 3) ://2 lei shang 1j
			kind2_node(j, i, -Num13_p2 - Num14_p2 - Num15_p2 - Num16_p2, -2, -1, 0, 1, Num17_p2 + Num18_p2 + Num19_p2 + Num18_p2,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_17_p2 + 1) ://2 lei xia 1j
			kind2_node(j, i, -(Num13_p2 - Node_13_p2 + Num14_p2 + Num15_p2 + Num16_p2 + Node_17_p2), -1, 0, 1, 2,
			(Num17_p2 - Node_17_p2 + Num18_p2 + Num19_p2 + Num18_p2 + Node_17_p2), B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q2 + Num19_p2 + Num18_p2 + Num17_p2) < (N_d0 - 3))// 4dx!
			{
				inner_5node(j, i, 4, -Num13_p2 - Num14_p2 - Num15_p2 - Num16_p2, Num17_p2 + Num18_p2 + Num19_p2 + Num18_p2, Coef, Coef_location);
			}
			else if (((i - M_post2_q2 + Num19_p2 + Num18_p2 + Num17_p2) > (N_d0 - 3)) && ((i - M_post2_q2 + Num19_p2 + Num18_p2 + Num17_p2) < Node171_p2))// 2dx!
			{
				inner_5node(j, i, 2, -(Num15_p2 + Num16_p2 + N_d0 - 2), Num17_p2 - N_d0 + 2 + Num18_p2, Coef, Coef_location);
			}
			else if (((i - M_post2_q2 + Num19_p2 + Num18_p2 + Num17_p2) >= Node171_p2) && ((i - M_post2_q2 + Num19_p2 + Num18_p2 + Num17_p2) <= Node172_p2))// 2 lei right 2j !
			{
				kind2_node(j, i, -((i - M_post2_q2 + Num19_p2 + Num18_p2 + Num17_p2) + Num16_p2 - (i - M_post2_q2 + Num19_p2 + Num18_p2 + Num17_p2 - Node171_p2) * 2),
					-1, 0, 1, Num17_p2 - N_d0 + 2 + Num18_p2, Num17_p2 + Num18_p2 + Num19_p2 + Num18_p2, B2, B3, B1_2, B3, B4, B5, Coef, Coef_location);
			}
			else// 2dx!
			{
				inner_5node(j, i, 2, -(Num15_p2 - Node151_p2 + Num16_p2 + Node_17_p2), Num17_p2 - N_d0 + 2 + Num18_p2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post2_q2 - Num19_p2 - Num18_p2 - Num17_p2 + Node_17_p2 + 2);
		i<(M_post2_q2 - Num19_p2 - Num18_p2 - Num17_p2 + Node174_p2 + 3); i++)// circle num17 :  80 node!  part2 !
	{
		switch (i - M_post2_q2 + Num19_p2 + Num18_p2 + Num17_p2)
		{
		case (Node_17_p2 + 2) :// 4dx
			inner_5node(j, i, 4, -(Num13_p2 - Node_13_p2 + Num14_p2 + Num15_p2 + Num16_p2 + Node_17_p2), (Num17_p2 - Node_17_p2 + Num18_p2 + Num19_p2 + Num18_p2 + Node_17_p2), Coef, Coef_location);
			break;
		case (Node_17_p2 + 3) :// 2lei shang 1j
			kind2_node(j, i, -(Num13_p2 - Node_13_p2 + Num14_p2 + Num15_p2 + Num16_p2 + Node_17_p2), -2, -1, 0, 1, (Num17_p2 - Node_17_p2 + Num18_p2 + Num19_p2 + Num18_p2 + Node_17_p2),
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node173_p2 + 2) :// 2lei shang 2j
			kind2_node(j, i, -(Num15_p2 - Node151_p2 - 1 + Num16_p2 + Node173_p2), -2, -1, 0, 1, Num17_p2 - Node173_p2 + Num18_p2 + Node191_p2 + 1,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case (Node173_p2 + 5) ://inner 4 node
			inner_4node(j, i, -(Num16_p2 - Node161_p2 - 1 + Node173_p2 + 3), -1, 0, Num17_p2 - Node173_p2 - 3, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node173_p2 + 6) ://inner 4 node
			inner_4node(j, i, -(Num16_p2 - Node161_p2 - 1 + Node173_p2 + 3) + 2, 0, 1, Num17_p2 - Node173_p2 - 3, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Node174_p2 - 2) :// 2lei xia 2j
			kind2_node(j, i, -(Num15_p2 - Node153_p2 + 1 + Num16_p2 + Node174_p2), -1, 0, 1, 2, Num17_p2 - Node174_p2 + Num18_p2 + Node193_p2 - 1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node174_p2 + 1) :// 2lei xia 1j
			kind2_node(j, i, -(Num13_p2 - Node133_p2 + Num14_p2 + Num15_p2 + Num16_p2 + Node174_p2), -1, 0, 1, 2,
			(Num17_p2 - Node174_p2 + Num18_p2 + Num19_p2 + Num18_p2 + Node174_p2), B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node174_p2 + 2) :// 4dx
			inner_5node(j, i, 4, -(Num13_p2 - Node133_p2 + Num14_p2 + Num15_p2 + Num16_p2 + Node174_p2),
			(Num17_p2 - Node174_p2 + Num18_p2 + Num19_p2 + Num18_p2 + Node174_p2), Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q2 + Num19_p2 + Num18_p2 + Num17_p2) < (Node173_p2 + 2))  // 2dx
			{
				inner_5node(j, i, 2, -(Num15_p2 - Node151_p2 - 1 + Num16_p2 + Node173_p2), Num17_p2 - Node173_p2 + Num18_p2 + Node191_p2 + 1, Coef, Coef_location);
			}
			else if (((i - M_post2_q2 + Num19_p2 + Num18_p2 + Num17_p2) > (Node173_p2 + 2)) && ((i - M_post2_q2 + Num19_p2 + Num18_p2 + Num17_p2) < (Node173_p2 + 5)))  // 1dx
			{
				inner_5node(j, i, 1, -(Num16_p2 - Node161_p2 - 1 + Node173_p2 + 3), Num17_p2 - Node173_p2 - 3, Coef, Coef_location);
			}
			else if (((i - M_post2_q2 + Num19_p2 + Num18_p2 + Num17_p2) > (Node173_p2 + 6)) && ((i - M_post2_q2 + Num19_p2 + Num18_p2 + Num17_p2) < (Node174_p2 - 2)))  // 1dx
			{
				inner_5node(j, i, 1, -(Num16_p2 - Node161_p2 - 1 + Node173_p2 + 3) + 2, Num17_p2 - Node173_p2 - 3, Coef, Coef_location);
			}
			else
			{
				inner_5node(j, i, 2, -(Num15_p2 - Node153_p2 + 1 + Num16_p2 + Node174_p2), Num17_p2 - Node174_p2 + Num18_p2 + Node193_p2 - 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post2_q2 - Num19_p2 - Num18_p2 - Num17_p2 + Node174_p2 + 3);
		i<(M_post2_q2 - Num19_p2 - Num18_p2); i++)// circle num17 :  80 node!  part3 !
	{
		switch (i - M_post2_q2 + Num19_p2 + Num18_p2 + Num17_p2)
		{
		case (Num17_p2 - 1) :
			boundary_node(j, i, 2, -Num14_p2 - Num15_p2 - Num16_p2 - Num17_p2, Num18_p2 + Num19_p2 + Num18_p2 + Num17_p2, Coef, Coef_location);
			break;
		case (Node174_p2 + 3) ://2 lei shang 1j
			kind2_node(j, i, -(Num13_p2 - Node133_p2 + Num14_p2 + Num15_p2 + Num16_p2 + Node174_p2), -2, -1, 0, 1, (Num17_p2 - Node174_p2 + Num18_p2 + Num19_p2 + Num18_p2 + Node174_p2),
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_1717_p2 + 1) ://2 lei xia 1j
			kind2_node(j, i, -Num14_p2 - Num15_p2 - Num16_p2 - Num17_p2, -1, 0, 1, 2, Num18_p2 + Num19_p2 + Num18_p2 + Num17_p2,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q2 + Num19_p2 + Num18_p2 + Num17_p2) < Node175_p2)// 2dx!
			{
				inner_5node(j, i, 2, -(Num15_p2 - Node153_p2 + Num16_p2 + Node175_p2 - 4), Num17_p2 - Node_1717_p2 + Num18_p2 + Num19_p2 - 1, Coef, Coef_location);
			}
			else if (((i - M_post2_q2 + Num19_p2 + Num18_p2 + Num17_p2) >= Node175_p2) && ((i - M_post2_q2 + Num19_p2 + Num18_p2 + Num17_p2) <= Node176_p2))// 2 lei left 2j !
			{
				kind2_node(j, i, -((i - M_post2_q2 + Num19_p2 + Num18_p2 + Num17_p2) + Num16_p2 - Node163_p2 - (i - M_post2_q2 + Num19_p2 + Num18_p2 + Num17_p2 - Node175_p2) * 2),
					-1, 0, 1, Num17_p2 - Node_1717_p2 + Num18_p2 + Num19_p2 - 1, Num17_p2 + Num18_p2 + Num19_p2 + Num18_p2, B2, B3, B1_2, B3, B4, B5, Coef, Coef_location);
			}
			else if (((i - M_post2_q2 + Num19_p2 + Num18_p2 + Num17_p2) > Node176_p2) && ((i - M_post2_q2 + Num19_p2 + Num18_p2 + Num17_p2) <= Node_1717_p2))// 2dx
			{
				inner_5node(j, i, 2, -(Num16_p2 + Node_1717_p2 + 1), Num17_p2 - Node_1717_p2 + Num18_p2 + Num19_p2 - 1, Coef, Coef_location);
			}
			else
			{
				inner_5node(j, i, 4, -Num14_p2 - Num15_p2 - Num16_p2 - Num17_p2, Num18_p2 + Num19_p2 + Num18_p2 + Num17_p2, Coef, Coef_location);// 4dx!
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (M_post2_q2 - Num19_p2 - Num18_p2); i<(M_post2_q2 - Num19_p2); i++)// num 18 :  6 node! 
	{
		switch (i - M_post2_q2 + Num19_p2 + Num18_p2)
		{
		case 0:// 1 lei shang 2j
			kind1_node(j, i, -(Num17_p2 - Node173_p2 - 3) - 2, -(Num17_p2 - Node173_p2 - 3) - 1, -(Num17_p2 - Node173_p2 - 3),
				0, 1, (Num18_p2 + Node192_p2 - 2) - 2, (Num18_p2 + Node192_p2 - 2) - 1, (Num18_p2 + Node192_p2 - 2),
				A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case 1:
			inner_5node(j, i, 1, -(Num17_p2 - Node173_p2 - 3), (Num18_p2 + Node192_p2 - 2), Coef, Coef_location);
			break;
		case 2://inner 4 node
			inner_4node(j, i, -(Num17_p2 - Node173_p2 - 3), -1, 0, (Num18_p2 + Node192_p2 - 2), 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case 3://inner 4 node
			inner_4node(j, i, -(Num17_p2 - Node173_p2 - 3), 0, 1, (Num18_p2 + Node192_p2 - 2), 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case 4:
			inner_5node(j, i, 1, -(Num17_p2 - Node173_p2 - 3), (Num18_p2 + Node192_p2 - 2), Coef, Coef_location);
			break;
		case (Num18_p2 - 1) :// 1 lei xia 2j
			kind1_node(j, i, -(Num17_p2 - Node173_p2 - 3), -(Num17_p2 - Node173_p2 - 3) + 1, -(Num17_p2 - Node173_p2 - 3) + 2,
			-1, 0, (Num18_p2 + Node192_p2 - 2), (Num18_p2 + Node192_p2 - 2) + 1, (Num18_p2 + Node192_p2 - 2) + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (M_post2_q2 - Num19_p2); i<(M_post2_q2 - Num19_p2 + Node191_p2 + 1); i++)// num 19 :  42 node! part1
	{
		switch (i - M_post2_q2 + Num19_p2)
		{
		case 0:// 1 lei shang
			kind1_node(j, i, -(Num17_p2 - N_d0 + 2 + Num18_p2) - 2, -(Num17_p2 - N_d0 + 2 + Num18_p2) - 1, -(Num17_p2 - N_d0 + 2 + Num18_p2),
				0, 1, (Num19_p2 + Num18_p2 + N_d0 - 2) - 2, (Num19_p2 + Num18_p2 + N_d0 - 2) - 1, (Num19_p2 + Num18_p2 + N_d0 - 2),
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case Node191_p2:// 1 lei xia
			kind1_node(j, i, -(Num17_p2 - N_d0 + 2 + Num18_p2), -(Num17_p2 - N_d0 + 2 + Num18_p2) + 1, -(Num17_p2 - N_d0 + 2 + Num18_p2) + 2,
				-1, 0, (Num19_p2 + Num18_p2 + N_d0 - 2), (Num19_p2 + Num18_p2 + N_d0 - 2) + 1, (Num19_p2 + Num18_p2 + N_d0 - 2) + 2,
				A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			inner_5node(j, i, 2, -(Num17_p2 - N_d0 + 2 + Num18_p2), (Num19_p2 + Num18_p2 + N_d0 - 2), Coef, Coef_location);
			break;
		}
		j = j + 1;
	}
	for (i = (M_post2_q2 - Num19_p2 + Node191_p2 + 1); i<(M_post2_q2 - Num19_p2 + Node193_p2); i++)// num 19 :  42 node! part2
	{
		switch (i - M_post2_q2 + Num19_p2)
		{
		case (Node191_p2 + 1) :// 1 lei 1j shang
			kind1_node(j, i, -(Num17_p2 - Node173_p2 + Num18_p2 + Node191_p2 + 1) - 2, -(Num17_p2 - Node173_p2 + Num18_p2 + Node191_p2 + 1) - 1,
			-(Num17_p2 - Node173_p2 + Num18_p2 + Node191_p2 + 1), 0, 1, (Num19_p2 - Node191_p2 - 1 + Num18_p2 + Node173_p2) - 2,
			(Num19_p2 - Node191_p2 - 1 + Num18_p2 + Node173_p2) - 1, (Num19_p2 - Node191_p2 - 1 + Num18_p2 + Node173_p2),
			A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case (Node191_p2 + 2) :// 2dx
			inner_5node(j, i, 2, -(Num17_p2 - Node173_p2 + Num18_p2 + Node191_p2 + 1), (Num19_p2 - Node191_p2 - 1 + Num18_p2 + Node173_p2), Coef, Coef_location);
			break;
		case (Node191_p2 + 3) :// 2 lei 2j shang
			kind2_node(j, i, -(Num17_p2 - Node173_p2 + Num18_p2 + Node191_p2 + 1), -2, -1, 0, 1, (Num19_p2 - Node191_p2 - 1 + Num18_p2 + Node173_p2),
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node192_p2://inner 4 node
			inner_4node(j, i, -(4 + Node192_p2), -1, 0, Num19_p2 - Node192_p2 + 2, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node192_p2 + 1) ://inner 4 node
			inner_4node(j, i, -(4 + Node192_p2), 0, 1, Num19_p2 - Node192_p2 + 2, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Node193_p2 - 3) :// 2 lei 2j xia
			kind2_node(j, i, -(Num17_p2 - Node174_p2 + Num18_p2 + Node193_p2 - 1), -1, 0, 1, 2, Num19_p2 - Node193_p2 + 1 + Num18_p2 + Node174_p2,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node193_p2 - 2) :// 2dx
			inner_5node(j, i, 2, -(Num17_p2 - Node174_p2 + Num18_p2 + Node193_p2 - 1), Num19_p2 - Node193_p2 + 1 + Num18_p2 + Node174_p2, Coef, Coef_location);
			break;
		case (Node193_p2 - 1) :// 1 lei 1j xia
			kind1_node(j, i, -(Num17_p2 - Node174_p2 + Num18_p2 + Node193_p2 - 1), -(Num17_p2 - Node174_p2 + Num18_p2 + Node193_p2 - 1) + 1,
			-(Num17_p2 - Node174_p2 + Num18_p2 + Node193_p2 - 1) + 2, -1, 0, Num19_p2 - Node193_p2 + 1 + Num18_p2 + Node174_p2,
			Num19_p2 - Node193_p2 + 1 + Num18_p2 + Node174_p2 + 1, Num19_p2 - Node193_p2 + 1 + Num18_p2 + Node174_p2 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			inner_5node(j, i, 1, -(4 + Node192_p2), Num19_p2 - Node192_p2 + 2, Coef, Coef_location);
			break;
		}
		j = j + 1;
	}
	for (i = (M_post2_q2 - Num19_p2 + Node193_p2); i<(M_post2_q2); i++)// num 19 :  42 node! part3
	{
		switch (i - M_post2_q2 + Num19_p2)
		{
		case Node193_p2:// 1 lei shang
			kind1_node(j, i, -(Num17_p2 - Node_1717_p2 + Num18_p2 - 1 + Num19_p2) - 2, -(Num17_p2 - Node_1717_p2 + Num18_p2 - 1 + Num19_p2) - 1,
				-(Num17_p2 - Node_1717_p2 + Num18_p2 - 1 + Num19_p2), 0, 1, (1 + Num18_p2 + Node_1717_p2) - 2,
				(1 + Num18_p2 + Node_1717_p2) - 1, (1 + Num18_p2 + Node_1717_p2),
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case (Num19_p2 - 1) :// 1 lei xia
			kind1_node(j, i, -(Num17_p2 - Node_1717_p2 + Num18_p2 - 1 + Num19_p2), -(Num17_p2 - Node_1717_p2 + Num18_p2 - 1 + Num19_p2) + 1,
			-(Num17_p2 - Node_1717_p2 + Num18_p2 - 1 + Num19_p2) + 2, -1, 0, (1 + Num18_p2 + Node_1717_p2),
			(1 + Num18_p2 + Node_1717_p2) + 1, (1 + Num18_p2 + Node_1717_p2) + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			inner_5node(j, i, 2, -(Num17_p2 - Node_1717_p2 + Num18_p2 - 1 + Num19_p2), (1 + Num18_p2 + Node_1717_p2), Coef, Coef_location);
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!!!
	for (i = (M_post2_q2); i<(M_post2_q2 + Num18_p2); i++)// num 20(18) :  6 node! 
	{
		switch (i - M_post2_q2)
		{
		case 0:// 1 lei shang 2j
			kind1_node(j, i, -(Num19_p2 - Node192_p2 + 2) - 2, -(Num19_p2 - Node192_p2 + 2) - 1, -(Num19_p2 - Node192_p2 + 2),
				0, 1, (Num18_p2 + Node173_p2 + 3) - 2, (Num18_p2 + Node173_p2 + 3) - 1, (Num18_p2 + Node173_p2 + 3),
				A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case 1:
			inner_5node(j, i, 1, -(Num19_p2 - Node192_p2 + 2), (Num18_p2 + Node173_p2 + 3), Coef, Coef_location);
			break;
		case 2://inner 4 node
			inner_4node(j, i, -(Num19_p2 - Node192_p2 + 2), -1, 0, (Num18_p2 + Node173_p2 + 3), 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case 3://inner 4 node
			inner_4node(j, i, -(Num19_p2 - Node192_p2 + 2), 0, 1, (Num18_p2 + Node173_p2 + 3), 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case 4:
			inner_5node(j, i, 1, -(Num19_p2 - Node192_p2 + 2), (Num18_p2 + Node173_p2 + 3), Coef, Coef_location);
			break;
		case (Num18_p2 - 1) :// 1 lei xia 2j
			kind1_node(j, i, -(Num19_p2 - Node192_p2 + 2), -(Num19_p2 - Node192_p2 + 2) + 1, -(Num19_p2 - Node192_p2 + 2) + 2,
			-1, 0, (Num18_p2 + Node173_p2 + 3), (Num18_p2 + Node173_p2 + 3) + 1, (Num18_p2 + Node173_p2 + 3) + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		}
		j = j + 1;
	}
	cout << "之前无错5" << endl;
	//******************************************************************************************
	for (i = (M_post2_q2 + Num18_p2);
		i<(M_post2_q2 + Num18_p2 + Node_17_p2 + 2); i++)// circle num21(17) :  80 node!  part1 !
	{
		switch (i - M_post2_q2 - Num18_p2)
		{
		case 0:
			boundary_node(j, i, 1, -Num17_p2 - Num18_p2 - Num19_p2 - Num18_p2, Num17_p2 + Num16_p2 + Num15_p2 + Num14_p2, Coef, Coef_location);
			break;
		case (N_d0 - 3) ://2 lei shang 1j
			kind2_node(j, i, -Num17_p2 - Num18_p2 - Num19_p2 - Num18_p2, -2, -1, 0, 1, Num17_p2 + Num16_p2 + Num15_p2 + Num14_p2,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_17_p2 + 1) ://2 lei xia 1j
			kind2_node(j, i, -(Num17_p2 - Node_17_p2 + Num18_p2 + Num19_p2 + Num18_p2 + Node_17_p2), -1, 0, 1, 2,
			(Num17_p2 - Node_17_p2 + Num16_p2 + Num15_p2 + Num14_p2 + Node_13_p2), B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q2 - Num18_p2) < (N_d0 - 3))// 4dx!
			{
				inner_5node(j, i, 4, -Num17_p2 - Num18_p2 - Num19_p2 - Num18_p2, Num17_p2 + Num16_p2 + Num15_p2 + Num14_p2, Coef, Coef_location);
			}
			else if (((i - M_post2_q2 - Num18_p2) > (N_d0 - 3)) && ((i - M_post2_q2 - Num18_p2) < Node171_p2))// 2dx!
			{
				inner_5node(j, i, 2, -(Num19_p2 + Num18_p2 + N_d0 - 2), Num17_p2 - N_d0 + 2 + Num16_p2, Coef, Coef_location);
			}
			else if (((i - M_post2_q2 - Num18_p2) >= Node171_p2) && ((i - M_post2_q2 - Num18_p2) <= Node172_p2))// 2 lei right 2j !
			{
				kind2_node(j, i, -Num17_p2 - Num18_p2 - Num19_p2 - Num18_p2, -(Num19_p2 + Num18_p2 + N_d0 - 2), -1, 0, 1,
					Num17_p2 - (i - M_post2_q2 - Num18_p2) + (i - M_post2_q2 - Num18_p2 - Node171_p2) * 2,
					B5, B4, B3, B1_2, B3, B2, Coef, Coef_location);
			}
			else// 2dx!
			{
				inner_5node(j, i, 2, -(Num19_p2 + Num18_p2 + N_d0 - 2), Num17_p2 - Node_17_p2 + Num16_p2 + Node151_p2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post2_q2 + Num18_p2 + Node_17_p2 + 2);
		i<(M_post2_q2 + Num18_p2 + Node174_p2 + 3); i++)// circle num 21(17) :  80 node!  part2 !
	{
		switch (i - M_post2_q2 - Num18_p2)
		{
		case (Node_17_p2 + 2) :// 4dx
			inner_5node(j, i, 4, -(Num17_p2 - Node_17_p2 + Num18_p2 + Num19_p2 + Num18_p2 + Node_17_p2),
			(Num17_p2 - Node_17_p2 + Num16_p2 + Num15_p2 + Num14_p2 + Node_13_p2), Coef, Coef_location);
			break;
		case (Node_17_p2 + 3) :// 2lei shang 1j
			kind2_node(j, i, -(Num17_p2 - Node_17_p2 + Num18_p2 + Num19_p2 + Num18_p2 + Node_17_p2),
			-2, -1, 0, 1, (Num17_p2 - Node_17_p2 + Num16_p2 + Num15_p2 + Num14_p2 + Node_13_p2),
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node173_p2 + 2) :// 2lei shang 2j
			kind2_node(j, i, -(Num19_p2 - Node191_p2 - 1 + Num18_p2 + Node173_p2), -2, -1, 0, 1, Num17_p2 - Node173_p2 + Num16_p2 + Node151_p2 + 1,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case (Node173_p2 + 5) ://inner 4 node
			inner_4node(j, i, -(Num18_p2 + Node173_p2 + 3), -1, 0, (Num17_p2 - Node173_p2 - 2 + Node161_p2), 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node173_p2 + 6) ://inner 4 node
			inner_4node(j, i, -(Num18_p2 + Node173_p2 + 3), 0, 1, (Num17_p2 - Node173_p2 - 2 + Node161_p2) + 2, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Node174_p2 - 2) :// 2lei xia 2j
			kind2_node(j, i, -(Num19_p2 - Node193_p2 + 1 + Num18_p2 + Node174_p2), -1, 0, 1, 2, Num17_p2 - Node174_p2 + Num16_p2 + Node153_p2 - 1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node174_p2 + 1) :// 2lei xia 1j
			kind2_node(j, i, -(Num17_p2 + Num18_p2 + Num19_p2 + Num18_p2), -1, 0, 1, 2,
			(Num17_p2 - Node174_p2 + Num16_p2 + Num15_p2 + Num14_p2 + Node133_p2), B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node174_p2 + 2) :// 4dx
			inner_5node(j, i, 4, -(Num17_p2 + Num18_p2 + Num19_p2 + Num18_p2),
			(Num17_p2 - Node174_p2 + Num16_p2 + Num15_p2 + Num14_p2 + Node133_p2), Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q2 - Num18_p2) < (Node173_p2 + 2))  // 2dx
			{
				inner_5node(j, i, 2, -(Num19_p2 - Node191_p2 - 1 + Num18_p2 + Node173_p2), Num17_p2 - Node173_p2 + Num16_p2 + Node151_p2 + 1, Coef, Coef_location);
			}
			else if (((i - M_post2_q2 - Num18_p2) > (Node173_p2 + 2)) && ((i - M_post2_q2 - Num18_p2) < (Node173_p2 + 5)))  // 1dx
			{
				inner_5node(j, i, 1, -(Num18_p2 + Node173_p2 + 3), (Num17_p2 - Node173_p2 - 2 + Node161_p2), Coef, Coef_location);
			}
			else if (((i - M_post2_q2 - Num18_p2) > (Node173_p2 + 6)) && ((i - M_post2_q2 - Num18_p2) < (Node174_p2 - 2)))  // 1dx
			{
				inner_5node(j, i, 1, -(Num18_p2 + Node173_p2 + 3), (Num17_p2 - Node173_p2 - 2 + Node161_p2) + 2, Coef, Coef_location);
			}
			else
			{
				inner_5node(j, i, 2, -(Num19_p2 - Node193_p2 + 1 + Num18_p2 + Node174_p2), Num17_p2 - Node174_p2 + Num16_p2 + Node153_p2 - 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post2_q2 + Num18_p2 + Node174_p2 + 3);
		i<(M_post2_q2 + Num18_p2 + Num17_p2); i++)// circle num 21(17) :  80 node!  part3 !
	{
		switch (i - M_post2_q2 - Num18_p2)
		{
		case (Num17_p2 - 1) :
			boundary_node(j, i, 2, -Num18_p2 - Num19_p2 - Num18_p2 - Num17_p2, Num16_p2 + Num15_p2 + Num14_p2 + Num13_p2, Coef, Coef_location);
			break;
		case (Node174_p2 + 3) ://2 lei shang 1j
			kind2_node(j, i, -(Num17_p2 + Num18_p2 + Num19_p2 + Num18_p2), -2, -1, 0, 1,
			(Num17_p2 - Node174_p2 + Num16_p2 + Num15_p2 + Num14_p2 + Node133_p2), B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_1717_p2 + 1) ://2 lei xia 1j
			kind2_node(j, i, -Num18_p2 - Num19_p2 - Num18_p2 - Num17_p2, -1, 0, 1, 2, Num16_p2 + Num15_p2 + Num14_p2 + Num13_p2,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q2 - Num18_p2) < Node175_p2)// 2dx!
			{
				inner_5node(j, i, 2, -(1 + Num18_p2 + Node_1717_p2), Num17_p2 - Node175_p2 + 4 + Num16_p2 + Node153_p2, Coef, Coef_location);
			}
			else if (((i - M_post2_q2 - Num18_p2) >= Node175_p2) && ((i - M_post2_q2 - Num18_p2) <= Node176_p2))// 2 lei left 2j !
			{
				kind2_node(j, i, -Num18_p2 - Num19_p2 - Num18_p2 - Num17_p2, -(1 + Num18_p2 + Node_1717_p2), -1, 0, 1,
					Num17_p2 - (i - M_post2_q2 - Num18_p2) + Node163_p2 + (i - M_post2_q2 - Num18_p2 - Node175_p2) * 2,
					B5, B4, B3, B1_2, B3, B2, Coef, Coef_location);
			}
			else if (((i - M_post2_q2 - Num18_p2) > Node176_p2) && ((i - M_post2_q2 - Num18_p2) <= Node_1717_p2))// 2dx
			{
				inner_5node(j, i, 2, -(1 + Num18_p2 + Node_1717_p2), Num17_p2 - Node_1717_p2 + Num16_p2 + Num15_p2 - 1, Coef, Coef_location);
			}
			else
			{
				inner_5node(j, i, 4, -Num18_p2 - Num19_p2 - Num18_p2 - Num17_p2, Num16_p2 + Num15_p2 + Num14_p2 + Num13_p2, Coef, Coef_location);// 4dx!
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (M_post2_q2 + Num18_p2 + Num17_p2);
		i<(M_post2_q2 + Num18_p2 + Num17_p2 + Node161_p2 + 1); i++)// num 22(16) :  34 node! part1
	{
		switch (i - M_post2_q2 - Num18_p2 - Num17_p2)
		{
		case 0:// 3 lei 2j shang
			kind3_node(j, i, -(Num17_p2 - Node171_p2) - 1, -(Num17_p2 - Node171_p2), 0, 1, (Num16_p2 + 5) - 2, (Num16_p2 + 5),
				C4, C3_2, C1_2, C2_2, C4, C3_2, Coef, Coef_location);
			break;
		case Node161_p2://node2: 3 lei 2j xia
			kind3_node(j, i, -(Num17_p2 - Node172_p2 + Node161_p2), -(Num17_p2 - Node172_p2 + Node161_p2) + 1, -1, 0, (Num16_p2 + 5),
				(Num16_p2 + 5) + 2, C3_2, C4, C2_2, C1_2, C3_2, C4, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q2 - Num18_p2 - Num17_p2) % 2 == 0)// 1dx
			{
				inner_5node(j, i, 1,
					-((i - M_post2_q2 - Num18_p2 - Num17_p2) + Num17_p2 - Node171_p2 - (i - M_post2_q2 - Num18_p2 - Num17_p2) / 2),
					(Num16_p2 + 5), Coef, Coef_location);
			}
			else//1 lei  2j left
			{
				kind1_node(j, i, -((i - M_post2_q2 - Num18_p2 - Num17_p2) + Num17_p2 - Node171_p2 - (i - M_post2_q2 - Num18_p2 - Num17_p2) / 2) - (Num19_p2 + Num18_p2 + N_d0 - 2),
					-((i - M_post2_q2 - Num18_p2 - Num17_p2) + Num17_p2 - Node171_p2 - (i - M_post2_q2 - Num18_p2 - Num17_p2) / 2) + 1 - (Num19_p2 + Num18_p2 + N_d0 - 2),
					-((i - M_post2_q2 - Num18_p2 - Num17_p2) + Num17_p2 - Node171_p2 - (i - M_post2_q2 - Num18_p2 - Num17_p2) / 2),
					-((i - M_post2_q2 - Num18_p2 - Num17_p2) + Num17_p2 - Node171_p2 - (i - M_post2_q2 - Num18_p2 - Num17_p2) / 2) + 1,
					-1, 0, 1, (Num16_p2 + 5),
					A2_2, A2_2, A1_2, A1_2, A3_2, A5_2, A3_2, A4_2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post2_q2 + Num18_p2 + Num17_p2 + Node161_p2 + 1);
		i<(M_post2_q2 + Num18_p2 + Num17_p2 + Node163_p2); i++)// num 22(16) :  34 node! part2
	{
		switch (i - M_post2_q2 - Num18_p2 - Num17_p2)
		{
		case (Node161_p2 + 1) :// 1 lei 2j shang
			kind1_node(j, i, -(Num17_p2 - Node173_p2 - 2 + Node161_p2) - 2, -(Num17_p2 - Node173_p2 - 2 + Node161_p2) - 1,
			-(Num17_p2 - Node173_p2 - 2 + Node161_p2), 0, 1, (Num16_p2 - Node161_p2 + Node151_p2 + 3) - 2, (Num16_p2 - Node161_p2 + Node151_p2 + 3) - 1,
			(Num16_p2 - Node161_p2 + Node151_p2 + 3), A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case Node162_p2:// inner 3 node
			inner_3node(j, i, -1, 0, (Num16_p2 - Node161_p2 + Node151_p2 + 3), 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node162_p2 + 1) :// inner 3 node
			inner_3node(j, i, 0, 1, (Num16_p2 - Node161_p2 + Node151_p2 + 3) + 2, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Node163_p2 - 1) :// 1 lei 2j xia
			kind1_node(j, i, -(Num17_p2 - Node173_p2 - 2 + Node161_p2) - 2, -(Num17_p2 - Node173_p2 - 2 + Node161_p2) - 2 + 1,
			-(Num17_p2 - Node173_p2 - 2 + Node161_p2) - 2 + 2, -1, 0, (Num16_p2 - Node161_p2 + Node151_p2 + 3) + 2, (Num16_p2 - Node161_p2 + Node151_p2 + 3) + 2 + 1,
			(Num16_p2 - Node161_p2 + Node151_p2 + 3) + 2 + 2, A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q2 - Num18_p2 - Num17_p2) < Node162_p2)// 1dx
			{
				inner_5node(j, i, 1, -(Num17_p2 - Node173_p2 - 2 + Node161_p2), (Num16_p2 - Node161_p2 + Node151_p2 + 3), Coef, Coef_location);
			}
			else // 2dx
			{
				inner_5node(j, i, 1, -(Num17_p2 - Node173_p2 - 2 + Node161_p2) - 2, (Num16_p2 - Node161_p2 + Node151_p2 + 3) + 2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post2_q2 + Num18_p2 + Num17_p2 + Node163_p2);
		i<(M_post2_q2 + Num18_p2 + Num17_p2 + Num16_p2); i++)// num 22(16) :  34 node! part3
	{
		switch (i - M_post2_q2 - Num18_p2 - Num17_p2 - Node163_p2)
		{
		case 0:// 3 lei 2j shang
			kind3_node(j, i, -(Num17_p2 - Node175_p2 + Node163_p2) - 1, -(Num17_p2 - Node175_p2 + Node163_p2), 0, 1,
				(Num16_p2 - Node163_p2 + Node153_p2 + 5) - 2, (Num16_p2 - Node163_p2 + Node153_p2 + 5),
				C4, C3_2, C1_2, C2_2, C4, C3_2, Coef, Coef_location);
			break;
		case Node161_p2:// 3 lei 2j xia
			kind3_node(j, i, -(Num17_p2 - Node176_p2 + Num16_p2 - 1), -(Num17_p2 - Node176_p2 + Num16_p2 - 1) + 1, -1, 0,
				(Num16_p2 - Node163_p2 + Node153_p2 + 5), (Num16_p2 - Node163_p2 + Node153_p2 + 5) + 2,
				C3_2, C4, C2_2, C1_2, C3_2, C4, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q2 - Num18_p2 - Num17_p2 - Node163_p2) % 2 == 0)// 1dx
			{
				inner_5node(j, i, 1,
					-((i - M_post2_q2 - Num18_p2 - Num17_p2) + Num17_p2 - Node175_p2 - (i - M_post2_q2 - Num18_p2 - Num17_p2 - Node163_p2) / 2),
					(Num16_p2 - Node163_p2 + Node153_p2 + 5), Coef, Coef_location);
			}
			else//1 lei 2j left
			{
				kind1_node(j, i, -((i - M_post2_q2 - Num18_p2 - Num17_p2) + Num17_p2 - Node175_p2 - (i - M_post2_q2 - Num18_p2 - Num17_p2 - Node163_p2) / 2) - (1 + Num18_p2 + Node_1717_p2),
					-((i - M_post2_q2 - Num18_p2 - Num17_p2) + Num17_p2 - Node175_p2 - (i - M_post2_q2 - Num18_p2 - Num17_p2 - Node163_p2) / 2) + 1 - (1 + Num18_p2 + Node_1717_p2),
					-((i - M_post2_q2 - Num18_p2 - Num17_p2) + Num17_p2 - Node175_p2 - (i - M_post2_q2 - Num18_p2 - Num17_p2 - Node163_p2) / 2),
					-((i - M_post2_q2 - Num18_p2 - Num17_p2) + Num17_p2 - Node175_p2 - (i - M_post2_q2 - Num18_p2 - Num17_p2 - Node163_p2) / 2) + 1,
					-1, 0, 1, (Num16_p2 - Node163_p2 + Node153_p2 + 5),
					A2_2, A2_2, A1_2, A1_2, A3_2, A5_2, A3_2, A4_2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (M_post2_q2 + Num18_p2 + Num17_p2 + Num16_p2);
		i<(M_post2_q2 + Num18_p2 + Num17_p2 + Num16_p2 + Node151_p2 + 1); i++)// num 23(15) :  62 node! part1
	{
		switch (i - M_post2_q2 - Num18_p2 - Num17_p2 - Num16_p2)
		{
		case 0:// 1 lei 1j shang
			kind1_node(j, i, -(Num17_p2 - N_d0 + 2 + Num16_p2) - 2, -(Num17_p2 - N_d0 + 2 + Num16_p2) - 1, -(Num17_p2 - N_d0 + 2 + Num16_p2),
				0, 1, (Num15_p2 + Num14_p2 + N_d0 - 2) - 2, (Num15_p2 + Num14_p2 + N_d0 - 2) - 1, (Num15_p2 + Num14_p2 + N_d0 - 2),
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case 3:// 2 lei 2j shang
			kind2_node(j, i, -(Num17_p2 - N_d0 + 2 + Num16_p2), -2, -1, 0, 1, (Num15_p2 + Num14_p2 + N_d0 - 2),
				B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case 4:// 3 lei 2j zuo
			kind3_node(j, i, -(Num17_p2 - N_d0 + 2 + Num16_p2) - 1, -(Num17_p2 - N_d0 + 2 + Num16_p2), -1, 0, 1, Num15_p2 - 4,
				C4, C4, C3_2, C1_2, C3_2, C2_2, Coef, Coef_location);
			break;
		case (Node151_p2 - 4) :// 3 lei 2j zuo
			kind3_node(j, i, -(Num17_p2 - Node_17_p2 + Num16_p2 + Node151_p2), -(Num17_p2 - Node_17_p2 + Num16_p2 + Node151_p2) + 1,
			-1, 0, 1, Num15_p2 - 4,
			C4, C4, C3_2, C1_2, C3_2, C2_2, Coef, Coef_location);
			break;
		case (Node151_p2 - 3) :// 2 lei 2j xia
			kind2_node(j, i, -(Num17_p2 - Node_17_p2 + Num16_p2 + Node151_p2), -1, 0, 1, 2, Num15_p2 - Node151_p2 + Num14_p2 + Node_13_p2,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case Node151_p2:// 1 lei 1j xia
			kind1_node(j, i, -(Num17_p2 - Node_17_p2 + Num16_p2 + Node151_p2), -(Num17_p2 - Node_17_p2 + Num16_p2 + Node151_p2) + 1,
				-(Num17_p2 - Node_17_p2 + Num16_p2 + Node151_p2) + 2, -1, 0, Num15_p2 - Node151_p2 + Num14_p2 + Node_13_p2,
				Num15_p2 - Node151_p2 + Num14_p2 + Node_13_p2 + 1, Num15_p2 - Node151_p2 + Num14_p2 + Node_13_p2 + 2,
				A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q2 - Num18_p2 - Num17_p2 - Num16_p2) < 3)// 2dx
			{
				inner_5node(j, i, 2, -(Num17_p2 - N_d0 + 2 + Num16_p2), (Num15_p2 + Num14_p2 + N_d0 - 2), Coef, Coef_location);
			}
			else if (((i - M_post2_q2 - Num18_p2 - Num17_p2 - Num16_p2) > 4) && ((i - M_post2_q2 - Num18_p2 - Num17_p2 - Num16_p2) < (Node151_p2 - 4)))// 1dx
			{
				inner_5node(j, i, 1, -(Num16_p2 + 5), Num15_p2 - 4, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, i, 2, -(Num17_p2 - Node_17_p2 + Num16_p2 + Node151_p2), Num15_p2 - Node151_p2 + Num14_p2 + Node_13_p2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post2_q2 + Num18_p2 + Num17_p2 + Num16_p2 + Node151_p2 + 1);
		i<(M_post2_q2 + Num18_p2 + Num17_p2 + Num16_p2 + Node153_p2); i++)// num 23(15) :  62 node! part2
	{
		switch (i - M_post2_q2 - Num18_p2 - Num17_p2 - Num16_p2)
		{
		case (Node151_p2 + 1) :// 1 lei 1j shang
			kind1_node(j, i, -(Num17_p2 - Node173_p2 + Num16_p2 + Node151_p2 + 1) - 2, -(Num17_p2 - Node173_p2 + Num16_p2 + Node151_p2 + 1) - 1,
			-(Num17_p2 - Node173_p2 + Num16_p2 + Node151_p2 + 1), 0, 1, (Num15_p2 - Node151_p2 - 1 + Num14_p2 + Node132_p2) - 2,
			(Num15_p2 - Node151_p2 - 1 + Num14_p2 + Node132_p2) - 1, (Num15_p2 - Node151_p2 - 1 + Num14_p2 + Node132_p2),
			A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case (Node151_p2 + 2) ://  2dx
			inner_5node(j, i, 2, -(Num17_p2 - Node173_p2 + Num16_p2 + Node151_p2 + 1), (Num15_p2 - Node151_p2 - 1 + Num14_p2 + Node132_p2), Coef, Coef_location);
			break;
		case (Node151_p2 + 3) :// 2 lei 2j shang
			kind2_node(j, i, -(Num17_p2 - Node173_p2 + Num16_p2 + Node151_p2 + 1), -2, -1, 0, 1, (Num15_p2 - Node151_p2 - 1 + Num14_p2 + Node132_p2),
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node152_p2:// inner 3 node
			inner_3node(j, i, -1, 0, (Num15_p2 - Node152_p2 + Node143_p2), 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node152_p2 + 1) :// inner 3 node
			inner_3node(j, i, 0, 1, (Num15_p2 - Node152_p2 + Node143_p2) + N_c1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Node153_p2 - 3) ://2 lei xia 2j
			kind2_node(j, i, -(Num17_p2 - Node174_p2 + Num16_p2 + Node153_p2 - 1), -1, 0, 1, 2, (Num15_p2 - Node153_p2 + 1 + Num14_p2 + Node133_p2),
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node153_p2 - 2) :// 2dx
			inner_5node(j, i, 2, -(Num17_p2 - Node174_p2 + Num16_p2 + Node153_p2 - 1), (Num15_p2 - Node153_p2 + 1 + Num14_p2 + Node133_p2), Coef, Coef_location);
			break;
		case (Node153_p2 - 1) :// 1 lei 1j xia
			kind1_node(j, i, -(Num17_p2 - Node174_p2 + Num16_p2 + Node153_p2 - 1), -(Num17_p2 - Node174_p2 + Num16_p2 + Node153_p2 - 1) + 1,
			-(Num17_p2 - Node174_p2 + Num16_p2 + Node153_p2 - 1) + 2, -1, 0, (Num15_p2 - Node153_p2 + 1 + Num14_p2 + Node133_p2),
			(Num15_p2 - Node153_p2 + 1 + Num14_p2 + Node133_p2) + 1, (Num15_p2 - Node153_p2 + 1 + Num14_p2 + Node133_p2) + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q2 - Num18_p2 - Num17_p2 - Num16_p2) < Node152_p2)// 1dx
			{
				inner_5node(j, i, 1, -(Num16_p2 - Node162_p2 + Node152_p2 - 1), (Num15_p2 - Node152_p2 + Node143_p2), Coef, Coef_location);
			}
			else if ((i - M_post2_q2 - Num18_p2 - Num17_p2 - Num16_p2) > (Node152_p2 + 1))// 1dx
			{
				inner_5node(j, i, 1, -(Num16_p2 - Node162_p2 + Node152_p2 - 1) - 2, (Num15_p2 - Node152_p2 + Node143_p2) + N_c1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post2_q2 + Num18_p2 + Num17_p2 + Num16_p2 + Node153_p2);
		i<(M_post2_q2 + Num18_p2 + Num17_p2 + Num16_p2 + Num15_p2); i++)// num 23(15) :  62 node! part3
	{
		switch (i - M_post2_q2 - Num18_p2 - Num17_p2 - Num16_p2)
		{
		case Node153_p2:// 1 lei 1j shang
			kind1_node(j, i, -(Num17_p2 - Node174_p2 - 4 + Num16_p2 + Node153_p2) - 2, -(Num17_p2 - Node174_p2 - 4 + Num16_p2 + Node153_p2) - 1,
				-(Num17_p2 - Node174_p2 - 4 + Num16_p2 + Node153_p2), 0, 1, (Num15_p2 - Node153_p2 + Num14_p2 + Node133_p2 + 4) - 2,
				(Num15_p2 - Node153_p2 + Num14_p2 + Node133_p2 + 4) - 1, (Num15_p2 - Node153_p2 + Num14_p2 + Node133_p2 + 4),
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case (Node153_p2 + 3) :// 2 lei 2j shang
			kind2_node(j, i, -(Num17_p2 - Node174_p2 - 4 + Num16_p2 + Node153_p2), -2, -1, 0, 1, (Num15_p2 - Node153_p2 + Num14_p2 + Node133_p2 + 4),
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case (Node153_p2 + 4) :// 3 lei 2j shang
			kind3_node(j, i, -(Num17_p2 - Node174_p2 - 4 + Num16_p2 + Node153_p2) - 1, -(Num17_p2 - Node174_p2 - 4 + Num16_p2 + Node153_p2),
			-1, 0, 1, (5 + Num14_p2 - 1),
			C4, C4, C3_2, C1_2, C3_2, C2_2, Coef, Coef_location);
			break;
		case (Num15_p2 - 5) :// 3 lei 2j xia
			kind3_node(j, i, -(Num17_p2 - Node_1717_p2 + Num16_p2 + Num15_p2 - 1), -(Num17_p2 - Node_1717_p2 + Num16_p2 + Num15_p2 - 1) + 1,
			-1, 0, 1, (5 + Num14_p2 - 1),
			C4, C4, C3_2, C1_2, C3_2, C2_2, Coef, Coef_location);
			break;
		case (Num15_p2 - 4) :// 2 lei 2j xia
			kind2_node(j, i, -(Num17_p2 - Node_1717_p2 + Num16_p2 + Num15_p2 - 1), -1, 0, 1, 2, (1 + Num14_p2 + Node_1313_p2),
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num15_p2 - 1) :// 1 lei 1j xia
			kind1_node(j, i, -(Num17_p2 - Node_1717_p2 + Num16_p2 + Num15_p2 - 1), -(Num17_p2 - Node_1717_p2 + Num16_p2 + Num15_p2 - 1) + 1,
			-(Num17_p2 - Node_1717_p2 + Num16_p2 + Num15_p2 - 1) + 2, -1, 0, (1 + Num14_p2 + Node_1313_p2),
			(1 + Num14_p2 + Node_1313_p2) + 1, (1 + Num14_p2 + Node_1313_p2) + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q2 - Num18_p2 - Num17_p2 - Num16_p2) < (Node153_p2 + 3))// 2dx
			{
				inner_5node(j, i, 2, -(Num17_p2 - Node174_p2 - 4 + Num16_p2 + Node153_p2), (Num15_p2 - Node153_p2 + Num14_p2 + Node133_p2 + 4), Coef, Coef_location);
			}
			else if (((i - M_post2_q2 - Num18_p2 - Num17_p2 - Num16_p2) > (Node153_p2 + 4)) && ((i - M_post2_q2 - Num18_p2 - Num17_p2 - Num16_p2) < (Num15_p2 - 4)))// 1dx
			{
				inner_5node(j, i, 1, -(1 + Num15_p2 - 6), (5 + Num14_p2 - 1), Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, i, 2, -(Num17_p2 - Node_1717_p2 + Num16_p2 + Num15_p2 - 1), (1 + Num14_p2 + Node_1313_p2), Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!! num 24 !!!
	for (i = (M_post2_q2 + Num18_p2 + Num17_p2 + Num16_p2 + Num15_p2);
		i<(M_post2_q2 + Num18_p2 + Num17_p2 + Num16_p2 + Num15_p2 + Node142_p2 + 5); i++)// num 24(14) :  45 node! part1
	{
		switch (i - M_post2_q2 - Num18_p2 - Num17_p2 - Num16_p2 - Num15_p2)
		{
		case 0:// 1 lei 2j shang
			kind1_node(j, i, -(Num15_p2 - 4) - 2, -(Num15_p2 - 4) - 1, -(Num15_p2 - 4), 0, 1, (Num14_p2 + N_d0 + 2) - 2,
				(Num14_p2 + N_d0 + 2) - 1, (Num14_p2 + N_d0 + 2), A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case (Node142_p2 + 4) :// 1 lei 2j xia
			kind1_node(j, i, -(Num15_p2 - 4), -(Num15_p2 - 4) + 1, -(Num15_p2 - 4) + 2, -1, 0,
			(Num14_p2 + N_d0 + 2) - N_c1, (Num14_p2 + N_d0 + 2) - N_c1 + 1, (Num14_p2 + N_d0 + 2) - N_c1 + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q2 - Num18_p2 - Num17_p2 - Num16_p2 - Num15_p2) <= Node141_p2)// 1dx
			{
				inner_5node(j, i, 1, -(Num15_p2 - 4), (Num14_p2 + N_d0 + 2), Coef, Coef_location);
			}
			else if ((i - M_post2_q2 - Num18_p2 - Num17_p2 - Num16_p2 - Num15_p2) >= Node142_p2) // 1dx
			{
				inner_5node(j, i, 1, -(Num15_p2 - 4), (Num14_p2 + N_d0 + 2) - N_c1, Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_4node(j, i, -(Num15_p2 - 4), -1, 0, 1, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post2_q2 + Num18_p2 + Num17_p2 + Num16_p2 + Num15_p2 + Node142_p2 + 5);
		i<(M_post2_q2 + Num18_p2 + Num17_p2 + Num16_p2 + Num15_p2 + Node145_p2 - 4); i++)// num 24(14) :  45 node! part2
	{
		switch (i - M_post2_q2 - Num18_p2 - Num17_p2 - Num16_p2 - Num15_p2)
		{
		case (Node142_p2 + 5) :// 1 lei 2j shang
			kind1_node(j, i, -(Num15_p2 - Node152_p2 + Node143_p2) - 2, -(Num15_p2 - Node152_p2 + Node143_p2) - 1,
			-(Num15_p2 - Node152_p2 + Node143_p2), 0, 1, (Num14_p2 - Node143_p2 + Node132_p2 + 7) - 2,
			(Num14_p2 - Node143_p2 + Node132_p2 + 7) - 1, (Num14_p2 - Node143_p2 + Node132_p2 + 7),
			A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case (Node145_p2 - 5) :// 1 lei 2j xia
			kind1_node(j, i, -(Num15_p2 - Node152_p2 + Node143_p2 + N_c1), -(Num15_p2 - Node152_p2 + Node143_p2 + N_c1) + 1,
			-(Num15_p2 - Node152_p2 + Node143_p2 + N_c1) + 2, -1, 0, (Num14_p2 - Node143_p2 + Node132_p2 + 7),
			(Num14_p2 - Node143_p2 + Node132_p2 + 7) + 1, (Num14_p2 - Node143_p2 + Node132_p2 + 7) + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q2 - Num18_p2 - Num17_p2 - Num16_p2 - Num15_p2) <= Node143_p2)// 1dx
			{
				inner_5node(j, i, 1, -(Num15_p2 - Node152_p2 + Node143_p2), (Num14_p2 - Node143_p2 + Node132_p2 + 7), Coef, Coef_location);
			}
			else if ((i - M_post2_q2 - Num18_p2 - Num17_p2 - Num16_p2 - Num15_p2) >= Node144_p2) // 1dx
			{
				inner_5node(j, i, 1, -(Num15_p2 - Node152_p2 + Node143_p2 + N_c1), (Num14_p2 - Node143_p2 + Node132_p2 + 7), Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_4node(j, i, -1, 0, 1, (Num14_p2 - Node143_p2 + Node132_p2 + 7), 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post2_q2 + Num18_p2 + Num17_p2 + Num16_p2 + Num15_p2 + Node145_p2 - 4);
		i<(M_post2_q2 + Num18_p2 + Num17_p2 + Num16_p2 + Num15_p2 + Num14_p2); i++)// num 24(14) :  45 node! part3
	{
		switch (i - M_post2_q2 - Num18_p2 - Num17_p2 - Num16_p2 - Num15_p2)
		{
		case (Node145_p2 - 4) :// 1 lei 2j shang
			kind1_node(j, i, -(5 + Num14_p2 - 1) - 2, -(5 + Num14_p2 - 1) - 1, -(5 + Num14_p2 - 1), 0, 1,
			(1 + Node_1313_p2 - 4) + N_c1 - 2, (1 + Node_1313_p2 - 4) + N_c1 - 1, (1 + Node_1313_p2 - 4) + N_c1,
			A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case (Num14_p2 - 1) :// 1 lei 2j xia
			kind1_node(j, i, -(5 + Num14_p2 - 1), -(5 + Num14_p2 - 1) + 1, -(5 + Num14_p2 - 1) + 2,
			-1, 0, (1 + Node_1313_p2 - 4), (1 + Node_1313_p2 - 4) + 1, (1 + Node_1313_p2 - 4) + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q2 - Num18_p2 - Num17_p2 - Num16_p2 - Num15_p2) <= Node145_p2)// 1dx
			{
				inner_5node(j, i, 1, -(5 + Num14_p2 - 1), (1 + Node_1313_p2 - 4) + N_c1, Coef, Coef_location);
			}
			else if ((i - M_post2_q2 - Num18_p2 - Num17_p2 - Num16_p2 - Num15_p2) >= Node146_p2) // 1dx
			{
				inner_5node(j, i, 1, -(5 + Num14_p2 - 1), (1 + Node_1313_p2 - 4), Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_4node(j, i, -(5 + Num14_p2 - 1), -1, 0, 1, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!! num 25 !!!
	for (i = (M_post2_q2 + Num18_p2 + Num17_p2 + Num16_p2 + Num15_p2 + Num14_p2);
		i<(M_post2_q2 + Num18_p2 + Num17_p2 + Num16_p2 + Num15_p2 + Num14_p2 + Node_13_p2 + 2); i++)// circle num 25(13) :  95 node!  part1 !
	{
		switch (i - M_post2_q2 - Num18_p2 - Num17_p2 - Num16_p2 - Num15_p2 - Num14_p2)
		{
		case 0:
			boundary_node(j, i, 1, -Num17_p2 - Num16_p2 - Num15_p2 - Num14_p2, Num13_p2 + Num12_p2 + Num11_p2 + Num10_p2, Coef, Coef_location);
			break;
		case (N_d0 - 3) ://2 lei shang 1j
			kind2_node(j, i, -Num17_p2 - Num16_p2 - Num15_p2 - Num14_p2, -2, -1, 0, 1, Num13_p2 + Num12_p2 + Num11_p2 + Num10_p2,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (N_d0 + 1) :// 2 lei 2j shang
			kind2_node(j, i, -(Num15_p2 + Num14_p2 + N_d0 - 2), -2, -1, 0, 1, Num13_p2 - N_d0 + 2 + Num12_p2,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node131_p2:// inner 3 node
			inner_3node(j, i, -(Num14_p2 + N_d0 + 2), -1, 0, 1, 1, k*k*dx*dx - 4, Coef, Coef_location);
			break;
		case (Node131_p2 + 1) :// inner 3 node
			inner_3node(j, i, -(Num14_p2 + N_d0 + 2) + N_c1, 0, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node_13_p2 - 3) ://2 lei xia 2j
			kind2_node(j, i, -(Num15_p2 - Node151_p2 + Num14_p2 + Node_13_p2), -1, 0, 1, 2, Num13_p2 - Node_13_p2 + Num12_p2 + Node112_p2,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_13_p2 + 1) ://2 lei xia 1j
			kind2_node(j, i, -(Num17_p2 - Node_17_p2 + Num16_p2 + Num15_p2 + Num14_p2 + Node_13_p2), -1, 0, 1, 2,
			Num13_p2 - Node_13_p2 + Num12_p2 + Num11_p2 + Num10_p2 + Node_9_p2, B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q2 - Num18_p2 - Num17_p2 - Num16_p2 - Num15_p2 - Num14_p2) < (N_d0 - 3))// 4dx!
			{
				inner_5node(j, i, 4, -Num17_p2 - Num16_p2 - Num15_p2 - Num14_p2, Num13_p2 + Num12_p2 + Num11_p2 + Num10_p2, Coef, Coef_location);
			}
			else if (((i - M_post2_q2 - Num18_p2 - Num17_p2 - Num16_p2 - Num15_p2 - Num14_p2) > (N_d0 - 3)) && ((i - M_post2_q2 - Num18_p2 - Num17_p2 - Num16_p2 - Num15_p2 - Num14_p2) < (N_d0 + 1)))// 2dx!
			{
				inner_5node(j, i, 2, -(Num15_p2 + Num14_p2 + N_d0 - 2), Num13_p2 - N_d0 + 2 + Num12_p2, Coef, Coef_location);
			}
			else if (((i - M_post2_q2 - Num18_p2 - Num17_p2 - Num16_p2 - Num15_p2 - Num14_p2) > (N_d0 + 1)) && ((i - M_post2_q2 - Num18_p2 - Num17_p2 - Num16_p2 - Num15_p2 - Num14_p2) < Node131_p2))// 1dx
			{
				inner_5node(j, i, 1, -(Num14_p2 + N_d0 + 2), (Num13_p2 - N_d0 - 2), Coef, Coef_location);
			}
			else if (((i - M_post2_q2 - Num18_p2 - Num17_p2 - Num16_p2 - Num15_p2 - Num14_p2) > (Node131_p2 + 1)) && ((i - M_post2_q2 - Num18_p2 - Num17_p2 - Num16_p2 - Num15_p2 - Num14_p2) < (Node_13_p2 - 3)))// 1dx
			{
				inner_5node(j, i, 1, -(Num14_p2 + N_d0 + 2) + N_c1, (Num13_p2 - N_d0 - 2) - 2, Coef, Coef_location);
			}
			else// 2dx!
			{
				inner_5node(j, i, 2, -(Num15_p2 - Node151_p2 + Num14_p2 + Node_13_p2), Num13_p2 - Node_13_p2 + Num12_p2 + Node112_p2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post2_q2 + Num18_p2 + Num17_p2 + Num16_p2 + Num15_p2 + Num14_p2 + Node_13_p2 + 2);
		i<(M_post2_q2 + Num18_p2 + Num17_p2 + Num16_p2 + Num15_p2 + Num14_p2 + Node133_p2 + 3); i++)// circle num 25(13) :  95 node!  part2 !
	{
		switch (i - M_post2_q2 - Num18_p2 - Num17_p2 - Num16_p2 - Num15_p2 - Num14_p2)
		{
		case (Node_13_p2 + 2) :// 4dx
			inner_5node(j, i, 4, -(Num17_p2 - Node_17_p2 + Num16_p2 + Num15_p2 + Num14_p2 + Node_13_p2), Num13_p2 - Node_13_p2 + Num12_p2 + Num11_p2 + Num10_p2 + Node_9_p2, Coef, Coef_location);
			break;
		case (Node_13_p2 + 3) :// 2lei shang 1j
			kind2_node(j, i, -(Num17_p2 - Node_17_p2 + Num16_p2 + Num15_p2 + Num14_p2 + Node_13_p2), -2, -1, 0, 1, Num13_p2 - Node_13_p2 + Num12_p2 + Num11_p2 + Num10_p2 + Node_9_p2,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node132_p2 + 2) :// 2lei shang 2j
			kind2_node(j, i, -(Num15_p2 - Node151_p2 - 1 + Num14_p2 + Node132_p2), -2, -1, 0, 1, (Num13_p2 - Node132_p2 + Num12_p2 + Node112_p2 + 1),
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case (Node132_p2 + 3) :// 3lei zuo 1j
			kind3_node(j, i, -(Num14_p2 - Node143_p2 + Node132_p2 + 7), -1, 0, 1,
			(Num13_p2 - Node132_p2 + Num12_p2 + Node112_p2 + 1) - 1, (Num13_p2 - Node132_p2 + Num12_p2 + Node112_p2 + 1),
			C2_2, C3_2, C1_2, C3_2, C4, C4, Coef, Coef_location);
			break;
		case (Node133_p2 - 3) :// 3lei zuo 1j
			kind3_node(j, i, -(Num14_p2 - Node143_p2 + Node132_p2 + 7), -1, 0, 1,
			(Num13_p2 - Node133_p2 + Num12_p2 + Node113_p2 - 1), (Num13_p2 - Node133_p2 + Num12_p2 + Node113_p2 - 1) + 1,
			C2_2, C3_2, C1_2, C3_2, C4, C4, Coef, Coef_location);
			break;
		case (Node133_p2 - 2) :// 2lei xia 2j
			kind2_node(j, i, -(Num15_p2 - Node153_p2 + 1 + Num14_p2 + Node133_p2), -1, 0, 1, 2, (Num13_p2 - Node133_p2 + Num12_p2 + Node113_p2 - 1),
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node133_p2 + 1) :// 2lei xia 1j
			kind2_node(j, i, -(Num17_p2 - Node174_p2 + Num16_p2 + Num15_p2 + Num14_p2 + Node133_p2), -1, 0, 1, 2,
			(Num13_p2 - Node133_p2 + Num12_p2 + Num11_p2 + Num10_p2 + Node93_p2), B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node133_p2 + 2) :// 4dx
			inner_5node(j, i, 4, -(Num17_p2 - Node174_p2 + Num16_p2 + Num15_p2 + Num14_p2 + Node133_p2),
			(Num13_p2 - Node133_p2 + Num12_p2 + Num11_p2 + Num10_p2 + Node93_p2), Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q2 - Num18_p2 - Num17_p2 - Num16_p2 - Num15_p2 - Num14_p2) < (Node132_p2 + 2))  // 2dx
			{
				inner_5node(j, i, 2, -(Num15_p2 - Node151_p2 - 1 + Num14_p2 + Node132_p2), Num13_p2 - Node132_p2 + Num12_p2 + Node112_p2 + 1, Coef, Coef_location);
			}
			else if ((i - M_post2_q2 - Num18_p2 - Num17_p2 - Num16_p2 - Num15_p2 - Num14_p2) > (Node133_p2 - 2))  // 2dx
			{
				inner_5node(j, i, 2, -(Num15_p2 - Node153_p2 + 1 + Num14_p2 + Node133_p2), Num13_p2 - Node133_p2 + Num12_p2 + Node113_p2 - 1, Coef, Coef_location);
			}
			else  // 1dx
			{
				inner_5node(j, i, 1, -(Num14_p2 - Node143_p2 + Node132_p2 + 7), Num13_p2 - Node132_p2 - 3 + Node122_p2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post2_q2 + Num18_p2 + Num17_p2 + Num16_p2 + Num15_p2 + Num14_p2 + Node133_p2 + 3);
		i<(M_post2_q2 + Num18_p2 + Num17_p2 + Num16_p2 + Num15_p2 + Num14_p2 + Num13_p2); i++)// circle num 25(13) :  95 node!  part3 !
	{
		switch (i - M_post2_q2 - Num18_p2 - Num17_p2 - Num16_p2 - Num15_p2 - Num14_p2)
		{
		case (Num13_p2 - 1) :
			boundary_node(j, i, 2, -Num16_p2 - Num15_p2 - Num14_p2 - Num13_p2, Num12_p2 + Num11_p2 + Num10_p2 + Num9_p2, Coef, Coef_location);
			break;
		case (Node133_p2 + 3) ://2 lei shang 1j
			kind2_node(j, i, -(Num17_p2 - Node174_p2 + Num16_p2 + Num15_p2 + Num14_p2 + Node133_p2), -2, -1, 0, 1,
			(Num13_p2 - Node133_p2 + Num12_p2 + Num11_p2 + Num10_p2 + Node93_p2), B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node133_p2 + 7) :// 2 lei 2j shang
			kind2_node(j, i, -(Num15_p2 - Node153_p2 + Num14_p2 + Node133_p2 + 4), -2, -1, 0, 1, Num13_p2 - Node133_p2 - 4 + Num12_p2 + Node113_p2,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node134_p2:// inner 3 node
			inner_3node(j, i, -(1 + Node_1313_p2 - 4) - N_c1, -1, 0, 1, 1, k*k*dx*dx - 4, Coef, Coef_location);
			break;
		case (Node134_p2 + 1) :// inner 3 node
			inner_3node(j, i, -(1 + Node_1313_p2 - 4), 0, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node_1313_p2 - 3) ://2 lei xia 2j
			kind2_node(j, i, -(Num14_p2 + Node_1313_p2 + 1), -1, 0, 1, 2, Num13_p2 - Node_1313_p2 + Num12_p2 + Num11_p2 - 1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_1313_p2 + 1) ://2 lei xia 1j
			kind2_node(j, i, -Num16_p2 - Num15_p2 - Num14_p2 - Num13_p2, -1, 0, 1, 2, Num12_p2 + Num11_p2 + Num10_p2 + Num9_p2,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q2 - Num18_p2 - Num17_p2 - Num16_p2 - Num15_p2 - Num14_p2) < (Node133_p2 + 7))// 2dx!
			{
				inner_5node(j, i, 2, -(Num15_p2 - Node153_p2 + Num14_p2 + Node133_p2 + 4), Num13_p2 - Node133_p2 - 4 + Num12_p2 + Node113_p2, Coef, Coef_location);
			}
			else if (((i - M_post2_q2 - Num18_p2 - Num17_p2 - Num16_p2 - Num15_p2 - Num14_p2) > (Node133_p2 + 7)) && ((i - M_post2_q2 - Num18_p2 - Num17_p2 - Num16_p2 - Num15_p2 - Num14_p2) < Node134_p2))// 1dx
			{
				inner_5node(j, i, 1, -(1 + Node_1313_p2 - 4) - N_c1, (Num13_p2 - Node_1313_p2 + 4 + Num12_p2 - 1) + 2, Coef, Coef_location);
			}
			else if (((i - M_post2_q2 - Num18_p2 - Num17_p2 - Num16_p2 - Num15_p2 - Num14_p2) > (Node134_p2 + 1)) && ((i - M_post2_q2 - Num18_p2 - Num17_p2 - Num16_p2 - Num15_p2 - Num14_p2) < (Node_1313_p2 - 3)))// 1dx
			{
				inner_5node(j, i, 1, -(1 + Node_1313_p2 - 4), (Num13_p2 - Node_1313_p2 + 4 + Num12_p2 - 1), Coef, Coef_location);
			}
			else if (((i - M_post2_q2 - Num18_p2 - Num17_p2 - Num16_p2 - Num15_p2 - Num14_p2) > (Node_1313_p2 - 3)) && ((i - M_post2_q2 - Num18_p2 - Num17_p2 - Num16_p2 - Num15_p2 - Num14_p2) < (Node_1313_p2 + 1)))// 2dx
			{
				inner_5node(j, i, 2, -(Num14_p2 + Node_1313_p2 + 1), Num13_p2 - Node_1313_p2 + Num12_p2 + Num11_p2 - 1, Coef, Coef_location);
			}
			else// 4dx!
			{
				inner_5node(j, i, 4, -Num16_p2 - Num15_p2 - Num14_p2 - Num13_p2, Num12_p2 + Num11_p2 + Num10_p2 + Num9_p2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!! num 26 !!!
	for (i = (M_post2_q3 - Num9_p2 - Num10_p2 - Num11_p2 - Num12_p2);
		i<(M_post2_q3 - Num9_p2 - Num10_p2 - Num11_p2 - Num12_p2 + Node122_p2 + 1); i++)// num 26(12) :  29 node!  part1
	{
		switch (i - M_post2_q3 + Num9_p2 + Num10_p2 + Num11_p2 + Num12_p2)
		{
		case 0:// 1 lei 2j shang
			kind1_node(j, i, -(Num13_p2 - N_d0 - 2) - 2, -(Num13_p2 - N_d0 - 2) - 1, -(Num13_p2 - N_d0 - 2), 0, 1, (Num12_p2 + 4) - 2,
				(Num12_p2 + 4) - 1, (Num12_p2 + 4), A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case Node121_p2:// inner 3 node
			inner_3node(j, i, -(Num13_p2 - N_d0 - 2), -1, 0, 1, 1, k*k*dx*dx - 4, Coef, Coef_location);
			break;
		case (Node121_p2 + 1) :// inner 3 node
			inner_3node(j, i, -(Num13_p2 - N_d0 - 2) + 2, 0, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case Node122_p2:// 1 lei 2j xia
			kind1_node(j, i, -(Num13_p2 - N_d0 - 2) + 2, -(Num13_p2 - N_d0 - 2) + 2 + 1, -(Num13_p2 - N_d0 - 2) + 2 + 2, -1, 0,
				(Num12_p2 + 4) - 2, (Num12_p2 + 4) - 2 + 1, (Num12_p2 + 4) - 2 + 2,
				A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q3 + Num9_p2 + Num10_p2 + Num11_p2 + Num12_p2) < Node121_p2)// 1dx
			{
				inner_5node(j, i, 1, -(Num13_p2 - N_d0 - 2), (Num12_p2 + 4), Coef, Coef_location);
			}
			else // 2dx
			{
				inner_5node(j, i, 1, -(Num13_p2 - N_d0 - 2) + 2, (Num12_p2 + 4) - 2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post2_q3 - Num9_p2 - Num10_p2 - Num11_p2 - Num12_p2 + Node122_p2 + 1);
		i<(M_post2_q3 - Num9_p2 - Num10_p2 - Num11_p2 - Num12_p2 + Node123_p2); i++)// num 26(12) :  29 node! part3
	{
		switch (i - M_post2_q3 + Num9_p2 + Num10_p2 + Num11_p2 + Num12_p2 - Node122_p2 - 1)
		{
		case (Node122_p2 + 1 - Node122_p2 - 1) :// 3 lei 2j shang !!! 1
			kind3_node(j, i, -(Num13_p2 - Node132_p2 - 3 + Node122_p2) - 2, -(Num13_p2 - Node132_p2 - 3 + Node122_p2), 0, 1,
			(Num12_p2 - Node122_p2 + Node112_p2 + 3) - 1, (Num12_p2 - Node122_p2 + Node112_p2 + 3),
			C4, C3_2, C1_2, C2_2, C4, C3_2, Coef, Coef_location);
			break;
		case (Node123_p2 - 1 - Node122_p2 - 1) :// 3 lei 2j xia !!! 5
			kind3_node(j, i, -(Num13_p2 - Node132_p2 - 3 + Node122_p2), -(Num13_p2 - Node132_p2 - 3 + Node122_p2) + 2, -1, 0,
			(Num12_p2 - Node123_p2 + Node113_p2 - 3), (Num12_p2 - Node123_p2 + Node113_p2 - 3) + 1,
			C3_2, C4, C2_2, C1_2, C3_2, C4, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q3 + Num9_p2 + Num10_p2 + Num11_p2 + Num12_p2 - Node122_p2 - 1) % 2 == 0)// 1dx
			{
				inner_5node(j, i, 1, -(Num13_p2 - Node132_p2 - 3 + Node122_p2),
					Num12_p2 - (i - M_post2_q3 + Num9_p2 + Num10_p2 + Num11_p2 + Num12_p2) + Node112_p2 + 4 + (i - M_post2_q3 + Num9_p2 + Num10_p2 + Num11_p2 + Num12_p2 - Node122_p2 - 1) / 2,
					Coef, Coef_location);
			}
			else// 1 lei 2j zuo
			{
				kind1_node(j, i, -(Num13_p2 - Node132_p2 - 3 + Node122_p2), -1, 0, 1,
					Num12_p2 - (i - M_post2_q3 + Num9_p2 + Num10_p2 + Num11_p2 + Num12_p2) + Node112_p2 + 4 + (i - M_post2_q3 + Num9_p2 + Num10_p2 + Num11_p2 + Num12_p2 - Node122_p2 - 1) / 2,
					Num12_p2 - (i - M_post2_q3 + Num9_p2 + Num10_p2 + Num11_p2 + Num12_p2) + Node112_p2 + 4 + (i - M_post2_q3 + Num9_p2 + Num10_p2 + Num11_p2 + Num12_p2 - Node122_p2 - 1) / 2 + 1,
					Num12_p2 - (i - M_post2_q3 + Num9_p2 + Num10_p2 + Num11_p2 + Num12_p2) + Node112_p2 + 4 + (i - M_post2_q3 + Num9_p2 + Num10_p2 + Num11_p2 + Num12_p2 - Node122_p2 - 1) / 2 + (Num11_p2 - Node112_p2 + Num10_p2 + Node_9_p2 + 3),
					Num12_p2 - (i - M_post2_q3 + Num9_p2 + Num10_p2 + Num11_p2 + Num12_p2) + Node112_p2 + 4 + (i - M_post2_q3 + Num9_p2 + Num10_p2 + Num11_p2 + Num12_p2 - Node122_p2 - 1) / 2 + (Num11_p2 - Node112_p2 + Num10_p2 + Node_9_p2 + 3) + 1,
					A4_2, A3_2, A5_2, A3_2, A1_2, A1_2, A2_2, A2_2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post2_q3 - Num9_p2 - Num10_p2 - Num11_p2 - Num12_p2 + Node123_p2); i<(M_post2_q3 - Num9_p2 - Num10_p2 - Num11_p2); i++)// num 26(12) :  29 node! part3
	{
		switch (i - M_post2_q3 + Num9_p2 + Num10_p2 + Num11_p2 + Num12_p2)
		{
		case Node123_p2:// 1 lei 2j shang
			kind1_node(j, i, -(Num13_p2 - Node_1313_p2 + 4 + Num12_p2 - 1) - 2 - 2, -(Num13_p2 - Node_1313_p2 + 4 + Num12_p2 - 1) - 2 - 1,
				-(Num13_p2 - Node_1313_p2 + 4 + Num12_p2 - 1) - 2, 0, 1,
				(1 + Num11_p2 - 5) + 2 - 2, (1 + Num11_p2 - 5) + 2 - 1, (1 + Num11_p2 - 5) + 2,
				A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case Node124_p2:// inner 3 node
			inner_3node(j, i, -(Num13_p2 - Node_1313_p2 + 4 + Num12_p2 - 1) - 2, -1, 0, 1, 1, k*k*dx*dx - 4, Coef, Coef_location);
			break;
		case (Node124_p2 + 1) :// inner 3 node
			inner_3node(j, i, -(Num13_p2 - Node_1313_p2 + 4 + Num12_p2 - 1), 0, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Num12_p2 - 1) :// 1 lei 2j xia
			kind1_node(j, i, -(Num13_p2 - Node_1313_p2 + 4 + Num12_p2 - 1), -(Num13_p2 - Node_1313_p2 + 4 + Num12_p2 - 1) + 1,
			-(Num13_p2 - Node_1313_p2 + 4 + Num12_p2 - 1) + 2, -1, 0,
			(1 + Num11_p2 - 5), (1 + Num11_p2 - 5) + 1, (1 + Num11_p2 - 5) + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q3 + Num9_p2 + Num10_p2 + Num11_p2 + Num12_p2) < Node124_p2)// 2dx
			{
				inner_5node(j, i, 1, -(Num13_p2 - Node_1313_p2 + 4 + Num12_p2 - 1) - 2, (1 + Num11_p2 - 5) + 2, Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_5node(j, i, 1, -(Num13_p2 - Node_1313_p2 + 4 + Num12_p2 - 1), (1 + Num11_p2 - 5), Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!! num 27 !!!
	for (i = (M_post2_q3 - Num9_p2 - Num10_p2 - Num11_p2); i<(M_post2_q3 - Num9_p2 - Num10_p2 - Num11_p2 + Node112_p2 + 1); i++)// num 27(11) :  41 node! part1
	{
		switch (i - M_post2_q3 + Num9_p2 + Num10_p2 + Num11_p2)
		{
		case 0:// 1 lei 1j shang
			kind1_node(j, i, -(Num13_p2 - N_d0 + 2 + Num12_p2) - 2, -(Num13_p2 - N_d0 + 2 + Num12_p2) - 1, -(Num13_p2 - N_d0 + 2 + Num12_p2),
				0, 1, (Num11_p2 + Num10_p2 + N_d0 - 2) - 2, (Num11_p2 + Num10_p2 + N_d0 - 2) - 1, (Num11_p2 + Num10_p2 + N_d0 - 2),
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case 3:// 2 lei 2j shang
			kind2_node(j, i, -(Num13_p2 - N_d0 + 2 + Num12_p2), -2, -1, 0, 1, (Num11_p2 + Num10_p2 + N_d0 - 2),
				B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node111_p2://inner 4 node
			inner_4node(j, i, -(Num12_p2 + 4), -1, 0, Num11_p2 - 4, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node111_p2 + 1) ://inner 4 node
			inner_4node(j, i, -(Num12_p2 + 4) + 2, 0, 1, Num11_p2 - 4, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Node112_p2 - 3) :// 2 lei 2j xia
			kind2_node(j, i, -(Num13_p2 - Node_13_p2 + Num12_p2 + Node112_p2), -1, 0, 1, 2, Num11_p2 - Node112_p2 + Num10_p2 + Node_9_p2,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case Node112_p2:// 1 lei 1j xia
			kind1_node(j, i, -(Num13_p2 - Node_13_p2 + Num12_p2 + Node112_p2), -(Num13_p2 - Node_13_p2 + Num12_p2 + Node112_p2) + 1,
				-(Num13_p2 - Node_13_p2 + Num12_p2 + Node112_p2) + 2, -1, 0, Num11_p2 - Node112_p2 + Num10_p2 + Node_9_p2,
				Num11_p2 - Node112_p2 + Num10_p2 + Node_9_p2 + 1, Num11_p2 - Node112_p2 + Num10_p2 + Node_9_p2 + 2,
				A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q3 + Num9_p2 + Num10_p2 + Num11_p2) < 3)// 2dx
			{
				inner_5node(j, i, 2, -(Num13_p2 - N_d0 + 2 + Num12_p2), (Num11_p2 + Num10_p2 + N_d0 - 2), Coef, Coef_location);
			}
			else if (((i - M_post2_q3 + Num9_p2 + Num10_p2 + Num11_p2) > 3) && ((i - M_post2_q3 + Num9_p2 + Num10_p2 + Num11_p2) < (Node111_p2)))// 1dx
			{
				inner_5node(j, i, 1, -(Num12_p2 + 4), Num11_p2 - 4, Coef, Coef_location);
			}
			else if (((i - M_post2_q3 + Num9_p2 + Num10_p2 + Num11_p2) > (Node111_p2 + 1)) && ((i - M_post2_q3 + Num9_p2 + Num10_p2 + Num11_p2) < (Node112_p2 - 3)))// 1dx
			{
				inner_5node(j, i, 1, -(Num12_p2 + 4) + 2, Num11_p2 - 4, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, i, 2, -(Num13_p2 - Node_13_p2 + Num12_p2 + Node112_p2), Num11_p2 - Node112_p2 + Num10_p2 + Node_9_p2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post2_q3 - Num9_p2 - Num10_p2 - Num11_p2 + Node112_p2 + 1); i<(M_post2_q3 - Num9_p2 - Num10_p2 - Num11_p2 + Node113_p2); i++)// num 27(11) :  41 node! part2
	{
		switch (i - M_post2_q3 + Num9_p2 + Num10_p2 + Num11_p2)
		{
		case (Node112_p2 + 1) :// 1lei shang 1j
			kind1_node(j, i, -(Num13_p2 - Node132_p2 + Num12_p2 + Node112_p2 + 1) - 2, -(Num13_p2 - Node132_p2 + Num12_p2 + Node112_p2 + 1) - 1,
			-(Num13_p2 - Node132_p2 + Num12_p2 + Node112_p2 + 1), 0, 1, (Num11_p2 - Node112_p2 - 1 + Num10_p2 + Node92_p2) - 2,
			(Num11_p2 - Node112_p2 - 1 + Num10_p2 + Node92_p2) - 1, (Num11_p2 - Node112_p2 - 1 + Num10_p2 + Node92_p2),
			A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case (Node113_p2 - 1) :// 1lei xia 1j
			kind1_node(j, i, -(Num13_p2 - Node133_p2 + Num12_p2 + Node113_p2 - 1), -(Num13_p2 - Node133_p2 + Num12_p2 + Node113_p2 - 1) + 1,
			-(Num13_p2 - Node133_p2 + Num12_p2 + Node113_p2 - 1) + 2, -1, 0, (Num11_p2 - Node112_p2 - 1 + Num10_p2 + Node92_p2),
			(Num11_p2 - Node112_p2 - 1 + Num10_p2 + Node92_p2) + 1, (Num11_p2 - Node112_p2 - 1 + Num10_p2 + Node92_p2) + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q3 + Num9_p2 + Num10_p2 + Num11_p2) < (Node112_p2 + 4))// 4dx
			{
				inner_5node(j, i, 2, -(Num13_p2 - Node132_p2 + Num12_p2 + Node112_p2 + 1), (Num11_p2 - Node112_p2 - 1 + Num10_p2 + Node92_p2), Coef, Coef_location);
			}
			else if ((i - M_post2_q3 + Num9_p2 + Num10_p2 + Num11_p2) > (Node113_p2 - 4))// 4dx
			{
				inner_5node(j, i, 2, -(Num13_p2 - Node133_p2 + Num12_p2 + Node113_p2 - 1), (Num11_p2 - Node112_p2 - 1 + Num10_p2 + Node92_p2), Coef, Coef_location);
			}
			else// 2lei zuo 2j
			{
				kind2_node(j, i, -(Num12_p2 - Node122_p2 - 1 - (i - M_post2_q3 + Num9_p2 + Num10_p2 + Num11_p2 - Node112_p2 - 4) * 2 + (i - M_post2_q3 + Num9_p2 + Num10_p2 + Num11_p2)),
					-1, 0, 1, (Num11_p2 - Node112_p2 - 1 + Num10_p2 + Node92_p2),
					(Num11_p2 - Node112_p2 - 1 + Num10_p2 + Node92_p2) + (Num9_p2 - Node92_p2 + Num8_p2 + Node72_p2),
					B2, B3, B1_2, B3, B4, B5, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post2_q3 - Num9_p2 - Num10_p2 - Num11_p2 + Node113_p2); i<(M_post2_q3 - Num9_p2 - Num10_p2); i++)// num 27(11) :  41 node! part3
	{
		switch (i - M_post2_q3 + Num9_p2 + Num10_p2 + Num11_p2)
		{
		case Node113_p2:// 1 lei 1j shang
			kind1_node(j, i, -(Num13_p2 - Node133_p2 - 4 + Num12_p2 + Node113_p2) - 2, -(Num13_p2 - Node133_p2 - 4 + Num12_p2 + Node113_p2) - 1,
				-(Num13_p2 - Node133_p2 - 4 + Num12_p2 + Node113_p2), 0, 1, (Num11_p2 - Node113_p2 + Num10_p2 + Node93_p2 + 4) - 2,
				(Num11_p2 - Node113_p2 + Num10_p2 + Node93_p2 + 4) - 1, (Num11_p2 - Node113_p2 + Num10_p2 + Node93_p2 + 4),
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case (Node113_p2 + 3) :// 2 lei 2j shang
			kind2_node(j, i, -(Num13_p2 - Node133_p2 - 4 + Num12_p2 + Node113_p2), -2, -1, 0, 1, (Num11_p2 - Node113_p2 + Num10_p2 + Node93_p2 + 4),
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node114_p2://inner 4 node
			inner_4node(j, i, -(Num11_p2 - 5 + 1) - 2, -1, 0, 5 + Num10_p2 - 1, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node114_p2 + 1) ://inner 4 node
			inner_4node(j, i, -(Num11_p2 - 5 + 1), 0, 1, 5 + Num10_p2 - 1, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Num11_p2 - 4) :// 2 lei 2j xia
			kind2_node(j, i, -(Num13_p2 - Node_1313_p2 + Num12_p2 + Num11_p2 - 1), -1, 0, 1, 2, (1 + Num10_p2 + Node_99_p2),
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num11_p2 - 1) :// 1 lei 1j xia
			kind1_node(j, i, -(Num13_p2 - Node_1313_p2 + Num12_p2 + Num11_p2 - 1), -(Num13_p2 - Node_1313_p2 + Num12_p2 + Num11_p2 - 1) + 1,
			-(Num13_p2 - Node_1313_p2 + Num12_p2 + Num11_p2 - 1) + 2, -1, 0,
			(1 + Num10_p2 + Node_99_p2), (1 + Num10_p2 + Node_99_p2) + 1, (1 + Num10_p2 + Node_99_p2) + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q3 + Num9_p2 + Num10_p2 + Num11_p2) < (Node113_p2 + 3))// 2dx
			{
				inner_5node(j, i, 2, -(Num13_p2 - Node133_p2 - 4 + Num12_p2 + Node113_p2), (Num11_p2 - Node113_p2 + Num10_p2 + Node93_p2 + 4), Coef, Coef_location);
			}
			else if (((i - M_post2_q3 + Num9_p2 + Num10_p2 + Num11_p2) > (Node113_p2 + 3)) && ((i - M_post2_q3 + Num9_p2 + Num10_p2 + Num11_p2) < (Node114_p2)))// 1dx
			{
				inner_5node(j, i, 1, -(Num11_p2 - 5 + 1) - 2, 5 + Num10_p2 - 1, Coef, Coef_location);
			}
			else if (((i - M_post2_q3 + Num9_p2 + Num10_p2 + Num11_p2) > (Node114_p2 + 1)) && ((i - M_post2_q3 + Num9_p2 + Num10_p2 + Num11_p2) < (Num11_p2 - 4)))// 1dx
			{
				inner_5node(j, i, 1, -(Num11_p2 - 5 + 1), 5 + Num10_p2 - 1, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, i, 2, -(Num13_p2 - Node_1313_p2 + Num12_p2 + Num11_p2 - 1), (1 + Num10_p2 + Node_99_p2), Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!! num 28 !!!
	for (i = (M_post2_q3 - Num9_p2 - Num10_p2); i<(M_post2_q3 - Num9_p2); i++)// num 28(10) :  12 node! 
	{
		switch (i - M_post2_q3 + Num9_p2 + Num10_p2)
		{
		case 0:// 1 lei shang 2j
			kind1_node(j, i, -(Num11_p2 - 4) - 2, -(Num11_p2 - 4) - 1, -(Num11_p2 - 4), 0, 1, Num10_p2 + N_d0 + 2 - 2, Num10_p2 + N_d0 + 2 - 1,
				Num10_p2 + N_d0 + 2, A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case 1:
			inner_5node(j, i, 1, -(Num11_p2 - 4), Num10_p2 + N_d0 + 2, Coef, Coef_location);
			break;
		case 2://inner 4 node
			inner_4node(j, i, -(Num11_p2 - 4), -1, 0, Num10_p2 + N_d0 + 2, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case 3://inner 4 node
			inner_4node(j, i, -(Num11_p2 - 4), 0, 1, Num10_p2 + N_d0 + 2, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case 4:
			inner_5node(j, i, 1, -(Num11_p2 - 4), Num10_p2 + N_d0 + 2, Coef, Coef_location);
			break;
		case (Num10_p2 / 2 - 1) :// 1 lei xia 2j
			kind1_node(j, i, -(Num11_p2 - 4), -(Num11_p2 - 4) + 1, -(Num11_p2 - 4) + 2, -1, 0, Num10_p2 + N_d0 + 2, Num10_p2 + N_d0 + 2 + 1,
			Num10_p2 + N_d0 + 2 + 2, A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		case (Num10_p2 / 2) :// 1 lei shang 2j
			kind1_node(j, i, -(5 + Num10_p2 - 1) - 2, -(5 + Num10_p2 - 1) - 1, -(5 + Num10_p2 - 1), 0, 1,
			Node_99_p2 - 4 + 1 - 2, Node_99_p2 - 4 + 1 - 1, Node_99_p2 - 4 + 1,
			A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case 7:
			inner_5node(j, i, 1, -(5 + Num10_p2 - 1), Node_99_p2 - 4 + 1, Coef, Coef_location);
			break;
		case 8://inner 4 node
			inner_4node(j, i, -(5 + Num10_p2 - 1), -1, 0, Node_99_p2 - 4 + 1, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case 9://inner 4 node
			inner_4node(j, i, -(5 + Num10_p2 - 1), 0, 1, Node_99_p2 - 4 + 1, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case 10:
			inner_5node(j, i, 1, -(5 + Num10_p2 - 1), Node_99_p2 - 4 + 1, Coef, Coef_location);
			break;
		case (Num10_p2 - 1) :// 1 lei xia 2j
			kind1_node(j, i, -(5 + Num10_p2 - 1), -(5 + Num10_p2 - 1) + 1, -(5 + Num10_p2 - 1) + 2, -1, 0,
			Node_99_p2 - 4 + 1, Node_99_p2 - 4 + 1 + 1, Node_99_p2 - 4 + 1 + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!! num 29 !!!
	for (i = (M_post2_q3 - Num9_p2); i<(M_post2_q3 - Num9_p2 + Node_9_p2 + 2); i++)// circle num 29(9) :  79 node!  part1 !
	{
		switch (i - M_post2_q3 + Num9_p2)
		{
		case 0:
			boundary_node(j, i, 1, -Num10_p2 - Num11_p2 - Num12_p2 - Num13_p2, Num9_p2 + Num8_p2 + Num7_p2 + Num6_p2, Coef, Coef_location);
			break;
		case (N_d0 - 3) ://2 lei shang 1j
			kind2_node(j, i, -Num10_p2 - Num11_p2 - Num12_p2 - Num13_p2, -2, -1, 0, 1, Num9_p2 + Num8_p2 + Num7_p2 + Num6_p2,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (N_d0 + 1) :// 2 lei 2j shang
			kind2_node(j, i, -(Num11_p2 + Num10_p2 + N_d0 - 2), -2, -1, 0, 1, Num9_p2 - N_d0 + 2 + Num8_p2,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node91_p2:// inner 4 node
			inner_4node(j, i, -(Num10_p2 + N_d0 + 2), -1, 0, Num9_p2 - N_d0 - 2, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node91_p2 + 1) :// inner 4 node
			inner_4node(j, i, -(Num10_p2 + N_d0 + 2), 0, 1, Num9_p2 - N_d0 - 2, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Node_9_p2 - 3) ://2 lei xia 2j
			kind2_node(j, i, -(Num11_p2 - Node112_p2 + Num10_p2 + Node_9_p2), -1, 0, 1, 2, Num9_p2 - Node_9_p2 + Num8_p2 + Node72_p2,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_9_p2 + 1) ://2 lei xia 1j
			kind2_node(j, i, -(Num13_p2 - Node_13_p2 + Num12_p2 + Num11_p2 + Num10_p2 + Node_9_p2), -1, 0, 1, 2,
			Num9_p2 - Node_9_p2 + Num8_p2 + Num7_p2 + Num6_p2 + Node_5_p2, B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q3 + Num9_p2) < (N_d0 - 3))// 4dx!
			{
				inner_5node(j, i, 4, -Num10_p2 - Num11_p2 - Num12_p2 - Num13_p2, Num9_p2 + Num8_p2 + Num7_p2 + Num6_p2, Coef, Coef_location);
			}
			else if (((i - M_post2_q3 + Num9_p2) > (N_d0 - 3)) && ((i - M_post2_q3 + Num9_p2) < (N_d0 + 1)))// 2dx!
			{
				inner_5node(j, i, 2, -(Num11_p2 + Num10_p2 + N_d0 - 2), Num9_p2 - N_d0 + 2 + Num8_p2, Coef, Coef_location);
			}
			else if (((i - M_post2_q3 + Num9_p2) > (N_d0 + 1)) && ((i - M_post2_q3 + Num9_p2) < Node91_p2))// 1dx
			{
				inner_5node(j, i, 1, -(Num10_p2 + N_d0 + 2), Num9_p2 - N_d0 - 2, Coef, Coef_location);
			}
			else if (((i - M_post2_q3 + Num9_p2) > (Node91_p2 + 1)) && ((i - M_post2_q3 + Num9_p2) < (Node_9_p2 - 3)))// 1dx
			{
				inner_5node(j, i, 1, -(Num10_p2 + N_d0 + 2), Num9_p2 - N_d0 - 2, Coef, Coef_location);
			}
			else// 2dx!
			{
				inner_5node(j, i, 2, -(Num11_p2 - Node112_p2 + Num10_p2 + Node_9_p2), Num9_p2 - Node_9_p2 + Num8_p2 + Node72_p2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post2_q3 - Num9_p2 + Node_9_p2 + 2); i<(M_post2_q3 - Num9_p2 + Node93_p2 + 3); i++)// circle num 29(9) :  79 node!  part2 !
	{
		switch (i - M_post2_q3 + Num9_p2)
		{
		case (Node92_p2 - 1) ://2 lei shang
			kind2_node(j, i, -(Num13_p2 - Node_13_p2 + Num12_p2 + Num11_p2 + Num10_p2 + Node_9_p2), -2, -1, 0, 1,
			Num9_p2 - Node_9_p2 + Num8_p2 + Num7_p2 + Num6_p2 + Node_5_p2, B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case Node92_p2://3 lei you
			kind3_node(j, i, -(Num11_p2 - Node112_p2 - 1 + Num10_p2 + Node92_p2), -1, 0, 1,
				(Num9_p2 - Node_9_p2 + Num8_p2 + Num7_p2 + Num6_p2 + Node_5_p2) - 1,
				Num9_p2 - Node_9_p2 + Num8_p2 + Num7_p2 + Num6_p2 + Node_5_p2, C2_1, C3_1, C1_1, C3_1, C4, C4, Coef, Coef_location);
			break;
		case Node93_p2://3 lei you
			kind3_node(j, i, -(Num11_p2 - Node112_p2 - 1 + Num10_p2 + Node92_p2), -1, 0, 1,
				(Num9_p2 - Node93_p2 + Num8_p2 + Num7_p2 + Num6_p2 + Node53_p2),
				(Num9_p2 - Node93_p2 + Num8_p2 + Num7_p2 + Num6_p2 + Node53_p2) + 1, C2_1, C3_1, C1_1, C3_1, C4, C4, Coef, Coef_location);
			break;
		case (Node93_p2 + 1) ://2 lei xia
			kind2_node(j, i, -(Num13_p2 - Node133_p2 + Num12_p2 + Num11_p2 + Num10_p2 + Node93_p2), -1, 0, 1, 2,
			(Num9_p2 - Node93_p2 + Num8_p2 + Num7_p2 + Num6_p2 + Node53_p2), B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q3 + Num9_p2) < (Node92_p2 - 1))// 4dx
			{
				inner_5node(j, i, 4, -(Num13_p2 - Node_13_p2 + Num12_p2 + Num11_p2 + Num10_p2 + Node_9_p2),
					Num9_p2 - Node_9_p2 + Num8_p2 + Num7_p2 + Num6_p2 + Node_5_p2, Coef, Coef_location);
			}
			else if ((i - M_post2_q3 + Num9_p2) > (Node93_p2 + 1))// 2dx!
			{
				inner_5node(j, i, 4, -(Num13_p2 - Node133_p2 + Num12_p2 + Num11_p2 + Num10_p2 + Node93_p2),
					(Num9_p2 - Node93_p2 + Num8_p2 + Num7_p2 + Num6_p2 + Node53_p2), Coef, Coef_location);
			}
			else if (((i - M_post2_q3 + Num9_p2) > Node92_p2) && ((i - M_post2_q3 + Num9_p2) < Node93_p2))
			{
				inner_5node(j, i, 2, -(Num11_p2 - Node112_p2 - 1 + Num10_p2 + Node92_p2),
					Num9_p2 - Node92_p2 + Num8_p2 + Node72_p2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post2_q3 - Num9_p2 + Node93_p2 + 3); i<(M_post2_q3); i++)// circle num 29(9) :  79 node!  part3 !
	{
		switch (i - M_post2_q3 + Num9_p2)
		{
		case (Num9_p2 - 1) :
			boundary_node(j, i, 2, -Num12_p2 - Num11_p2 - Num10_p2 - Num9_p2, Num8_p2 + Num7_p2 + Num6_p2 + Num5_p2, Coef, Coef_location);
			break;
		case (Node93_p2 + 3) ://2 lei shang 1j
			kind2_node(j, i, -(Num13_p2 - Node133_p2 + Num12_p2 + Num11_p2 + Num10_p2 + Node93_p2), -2, -1, 0, 1,
			(Num9_p2 - Node93_p2 + Num8_p2 + Num7_p2 + Num6_p2 + Node53_p2), B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node94_p2 - 3) :// 2 lei 2j shang
			kind2_node(j, i, -(Num11_p2 - Node113_p2 + Num10_p2 + Node93_p2 + 4), -2, -1, 0, 1,
			Num9_p2 - Node93_p2 - 4 + Num8_p2 + Node73_p2, B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node94_p2:// inner 4 node
			inner_4node(j, i, -(Node_99_p2 - 4 + 1), -1, 0, Num9_p2 - Node_99_p2 + 4 + Num8_p2 - 1, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node94_p2 + 1) :// inner 4 node
			inner_4node(j, i, -(Node_99_p2 - 4 + 1), 0, 1, Num9_p2 - Node_99_p2 + 4 + Num8_p2 - 1, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Node_99_p2 - 3) ://2 lei xia 2j
			kind2_node(j, i, -(Num10_p2 + Node_99_p2 + 1), -1, 0, 1, 2, Num9_p2 - Node_99_p2 + Num8_p2 + Num7_p2 - 1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_99_p2 + 1) ://2 lei xia 1j
			kind2_node(j, i, -Num12_p2 - Num11_p2 - Num10_p2 - Num9_p2, -1, 0, 1, 2, Num8_p2 + Num7_p2 + Num6_p2 + Num5_p2,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q3 + Num9_p2) < (Node94_p2 - 3))// 2dx!
			{
				inner_5node(j, i, 2, -(Num11_p2 - Node113_p2 + Num10_p2 + Node93_p2 + 4), Num9_p2 - Node93_p2 - 4 + Num8_p2 + Node73_p2, Coef, Coef_location);
			}
			else if (((i - M_post2_q3 + Num9_p2) > (Node94_p2 - 3)) && ((i - M_post2_q3 + Num9_p2) < Node94_p2))// 1dx
			{
				inner_5node(j, i, 1, -(Node_99_p2 - 4 + 1), Num9_p2 - Node_99_p2 + 4 + Num8_p2 - 1, Coef, Coef_location);
			}
			else if (((i - M_post2_q3 + Num9_p2) > (Node94_p2 + 1)) && ((i - M_post2_q3 + Num9_p2) < (Node_99_p2 - 3)))// 1dx
			{
				inner_5node(j, i, 1, -(Node_99_p2 - 4 + 1), Num9_p2 - Node_99_p2 + 4 + Num8_p2 - 1, Coef, Coef_location);
			}
			else if (((i - M_post2_q3 + Num9_p2) > (Node_99_p2 - 3)) && ((i - M_post2_q3 + Num9_p2) < (Node_99_p2 + 1)))// 2dx
			{
				inner_5node(j, i, 2, -(Num10_p2 + Node_99_p2 + 1), Num9_p2 - Node_99_p2 + Num8_p2 + Num7_p2 - 1, Coef, Coef_location);
			}
			else// 4dx!
			{
				inner_5node(j, i, 4, -Num12_p2 - Num11_p2 - Num10_p2 - Num9_p2, Num8_p2 + Num7_p2 + Num6_p2 + Num5_p2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!! num 30 !!!
	for (i = (M_post2_q3); i<(M_post2_q3 + Num8_p2); i++)// num 30(8) :  12 node! 
	{
		switch (i - M_post2_q3)
		{
		case 0:// 1 lei shang 2j
			kind1_node(j, i, -(Num9_p2 - N_d0 - 2) - 2, -(Num9_p2 - N_d0 - 2) - 1, -(Num9_p2 - N_d0 - 2), 0, 1, Num8_p2 + 4 - 2,
				Num8_p2 + 4 - 1, Num8_p2 + 4, A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case 1:
			inner_5node(j, i, 1, -(Num9_p2 - N_d0 - 2), Num8_p2 + 4, Coef, Coef_location);
			break;
		case 2://inner 4 node
			inner_4node(j, i, -(Num9_p2 - N_d0 - 2), -1, 0, Num8_p2 + 4, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case 3://inner 4 node
			inner_4node(j, i, -(Num9_p2 - N_d0 - 2), 0, 1, Num8_p2 + 4, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case 4:
			inner_5node(j, i, 1, -(Num9_p2 - N_d0 - 2), Num8_p2 + 4, Coef, Coef_location);
			break;
		case (Num8_p2 / 2 - 1) :// 1 lei xia 2j
			kind1_node(j, i, -(Num9_p2 - N_d0 - 2), -(Num9_p2 - N_d0 - 2) + 1, -(Num9_p2 - N_d0 - 2) + 2, -1, 0, Num8_p2 + 4,
			Num8_p2 + 4 + 1, Num8_p2 + 4 + 2, A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		case (Num8_p2 / 2) :// 1 lei shang 2j
			kind1_node(j, i, -(Num9_p2 - Node_99_p2 + 4 + Num8_p2 - 1) - 2, -(Num9_p2 - Node_99_p2 + 4 + Num8_p2 - 1) - 1,
			-(Num9_p2 - Node_99_p2 + 4 + Num8_p2 - 1), 0, 1, Num7_p2 - 4 - 2, Num7_p2 - 4 - 1, Num7_p2 - 4,
			A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case 7:
			inner_5node(j, i, 1, -(Num9_p2 - Node_99_p2 + 4 + Num8_p2 - 1), Num7_p2 - 4, Coef, Coef_location);
			break;
		case 8://inner 4 node
			inner_4node(j, i, -(Num9_p2 - Node_99_p2 + 4 + Num8_p2 - 1), -1, 0, Num7_p2 - 4, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case 9://inner 4 node
			inner_4node(j, i, -(Num9_p2 - Node_99_p2 + 4 + Num8_p2 - 1), 0, 1, Num7_p2 - 4, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case 10:
			inner_5node(j, i, 1, -(Num9_p2 - Node_99_p2 + 4 + Num8_p2 - 1), Num7_p2 - 4, Coef, Coef_location);
			break;
		case (Num8_p2 - 1) :// 1 lei xia 2j
			kind1_node(j, i, -(Num9_p2 - Node_99_p2 + 4 + Num8_p2 - 1), -(Num9_p2 - Node_99_p2 + 4 + Num8_p2 - 1) + 1,
			-(Num9_p2 - Node_99_p2 + 4 + Num8_p2 - 1) + 2, -1, 0, Num7_p2 - 4, Num7_p2 - 4 + 1, Num7_p2 - 4 + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!! num 31 !!!
	for (i = (M_post2_q3 + Num8_p2); i<(M_post2_q3 + Num8_p2 + Node72_p2 + 1); i++)// num 31(7) :  39 node! part1
	{
		switch (i - M_post2_q3 - Num8_p2)
		{
		case 0:// 1 lei 1j shang
			kind1_node(j, i, -(Num9_p2 - N_d0 + 2 + Num8_p2) - 2, -(Num9_p2 - N_d0 + 2 + Num8_p2) - 1, -(Num9_p2 - N_d0 + 2 + Num8_p2),
				0, 1, Num7_p2 + Num6_p2 + N_d0 - 2 - 2, Num7_p2 + Num6_p2 + N_d0 - 2 - 1, Num7_p2 + Num6_p2 + N_d0 - 2,
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case 3:// 2 lei 2j shang
			kind2_node(j, i, -(Num9_p2 - N_d0 + 2 + Num8_p2), -2, -1, 0, 1, Num7_p2 + Num6_p2 + N_d0 - 2,
				B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node71_p2://inner 4 node
			inner_4node(j, i, -(Num8_p2 + 4), -1, 0, Num7_p2 - 4, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node71_p2 + 1) ://inner 4 node
			inner_4node(j, i, -(Num8_p2 + 4), 0, 1, Num7_p2 - 4 + 2, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Node72_p2 - 3) :// 2 lei 2j xia
			kind2_node(j, i, -(Num9_p2 - Node_9_p2 + Num8_p2 + Node72_p2), -1, 0, 1, 2, (Num7_p2 - Node72_p2 + Num6_p2 + Node_5_p2),
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case Node72_p2:// 1 lei 1j xia
			kind1_node(j, i, -(Num9_p2 - Node_9_p2 + Num8_p2 + Node72_p2), -(Num9_p2 - Node_9_p2 + Num8_p2 + Node72_p2) + 1,
				-(Num9_p2 - Node_9_p2 + Num8_p2 + Node72_p2) + 2, -1, 0, (Num7_p2 - Node72_p2 + Num6_p2 + Node_5_p2),
				(Num7_p2 - Node72_p2 + Num6_p2 + Node_5_p2) + 1, (Num7_p2 - Node72_p2 + Num6_p2 + Node_5_p2) + 2,
				A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q3 - Num8_p2) < 3)// 2dx
			{
				inner_5node(j, i, 2, -(Num9_p2 - N_d0 + 2 + Num8_p2), Num7_p2 + Num6_p2 + N_d0 - 2, Coef, Coef_location);
			}
			else if (((i - M_post2_q3 - Num8_p2) > 3) && ((i - M_post2_q3 - Num8_p2) < (Node71_p2)))// 1dx
			{
				inner_5node(j, i, 1, -(Num8_p2 + 4), Num7_p2 - 4, Coef, Coef_location);
			}
			else if (((i - M_post2_q3 - Num8_p2) > (Node71_p2 + 1)) && ((i - M_post2_q3 - Num8_p2) < (Node72_p2 - 3)))// 1dx
			{
				inner_5node(j, i, 1, -(Num8_p2 + 4), Num7_p2 - 4 + 2, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, i, 2, -(Num9_p2 - Node_9_p2 + Num8_p2 + Node72_p2), (Num7_p2 - Node72_p2 + Num6_p2 + Node_5_p2), Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post2_q3 + Num8_p2 + Node72_p2 + 1); i<(M_post2_q3 + Num8_p2 + Node73_p2); i++)// num 31(7) :  39 node! part2
	{
		switch (i - M_post2_q3 - Num8_p2 - Node72_p2 - 1)
		{
		case 0:// 3 lei 1j shang
			kind3_node(j, i, -(Num9_p2 - Node92_p2 + Num8_p2 + Node72_p2) - 2, -(Num9_p2 - Node92_p2 + Num8_p2 + Node72_p2),
				0, 1, (Num7_p2 - Node72_p2 - 1 + Num6_p2 + Node52_p2) - 1, (Num7_p2 - Node72_p2 - 1 + Num6_p2 + Node52_p2),
				C4, C3_1, C1_1, C2_1, C4, C3_1, Coef, Coef_location);
			break;
		case 10:// 3 lei 1j xia
			kind3_node(j, i, -(Num9_p2 - Node92_p2 + Num8_p2 + Node72_p2), -(Num9_p2 - Node92_p2 + Num8_p2 + Node72_p2) + 2,
				-1, 0, (Num7_p2 - Node73_p2 + 1 + Num6_p2 + Node53_p2), (Num7_p2 - Node73_p2 + 1 + Num6_p2 + Node53_p2) + 1,
				C3_1, C4, C2_1, C1_1, C3_1, C4, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q3 - Num8_p2 - Node72_p2 - 1) % 2 == 0)// 1dx
			{
				inner_5node(j, i, 2, -(Num9_p2 - Node92_p2 + Num8_p2 + Node72_p2),
					Num7_p2 - (i - M_post2_q3 - Num8_p2) + Num6_p2 + Node52_p2 + (i - M_post2_q3 - Num8_p2 - Node72_p2 - 1) / 2,
					Coef, Coef_location);
			}
			else//1 lei  1j left
			{
				kind1_node(j, i, -(Num9_p2 - Node92_p2 + Num8_p2 + Node72_p2), -1, 0, 1,
					Num7_p2 - (i - M_post2_q3 - Num8_p2) + Num6_p2 + Node52_p2 + (i - M_post2_q3 - Num8_p2 - Node72_p2 - 1) / 2,
					Num7_p2 - (i - M_post2_q3 - Num8_p2) + Num6_p2 + Node52_p2 + (i - M_post2_q3 - Num8_p2 - Node72_p2 - 1) / 2 + 1,
					Num7_p2 - (i - M_post2_q3 - Num8_p2) + Num6_p2 + Node52_p2 + (i - M_post2_q3 - Num8_p2 - Node72_p2 - 1) / 2 + (Num5_p2 - Node_5_p2 + Num4_p2 + Num3_p2 + Num2_p2 + Node_1_p2),
					Num7_p2 - (i - M_post2_q3 - Num8_p2) + Num6_p2 + Node52_p2 + (i - M_post2_q3 - Num8_p2 - Node72_p2 - 1) / 2 + (Num5_p2 - Node_5_p2 + Num4_p2 + Num3_p2 + Num2_p2 + Node_1_p2) + 1,
					A4_1, A3_1, A5_1, A3_1, A1_1, A1_1, A2_1, A2_1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post2_q3 + Num8_p2 + Node73_p2); i<(M_post2_q3 + Num8_p2 + Num7_p2); i++)// num 31(7) :  39 node! part3
	{
		switch (i - M_post2_q3 - Num8_p2)
		{
		case Node73_p2:// 1 lei 1j shang
			kind1_node(j, i, -(Num9_p2 - Node93_p2 - 4 + Num8_p2 + Node73_p2) - 2, -(Num9_p2 - Node93_p2 - 4 + Num8_p2 + Node73_p2) - 1,
				-(Num9_p2 - Node93_p2 - 4 + Num8_p2 + Node73_p2), 0, 1, (Num7_p2 - Node73_p2 + Num6_p2 + Node53_p2 + 4) - 2,
				(Num7_p2 - Node73_p2 + Num6_p2 + Node53_p2 + 4) - 1, (Num7_p2 - Node73_p2 + Num6_p2 + Node53_p2 + 4),
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case (Node73_p2 + 3) :// 2 lei 2j shang
			kind2_node(j, i, -(Num9_p2 - Node93_p2 - 4 + Num8_p2 + Node73_p2), -2, -1, 0, 1,
			(Num7_p2 - Node73_p2 + Num6_p2 + Node53_p2 + 4), B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node74_p2://inner 4 node
			inner_4node(j, i, -(Num7_p2 - 5 + 1), -1, 0, (5 + Num6_p2 - 1) - 2, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node74_p2 + 1) ://inner 4 node
			inner_4node(j, i, -(Num7_p2 - 5 + 1), 0, 1, 5 + Num6_p2 - 1, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Num7_p2 - 4) :// 2 lei 2j xia
			kind2_node(j, i, -(Num9_p2 - Node_99_p2 + Num8_p2 + Num7_p2 - 1), -1, 0, 1, 2, (1 + Num6_p2 + Node_55_p2),
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num7_p2 - 1) :// 1 lei 1j xia
			kind1_node(j, i, -(Num9_p2 - Node_99_p2 + Num8_p2 + Num7_p2 - 1), -(Num9_p2 - Node_99_p2 + Num8_p2 + Num7_p2 - 1) + 1,
			-(Num9_p2 - Node_99_p2 + Num8_p2 + Num7_p2 - 1) + 2, -1, 0,
			(1 + Num6_p2 + Node_55_p2), (1 + Num6_p2 + Node_55_p2) + 1, (1 + Num6_p2 + Node_55_p2) + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q3 - Num8_p2) < (Node73_p2 + 3))// 2dx
			{
				inner_5node(j, i, 2, -(Num9_p2 - Node93_p2 - 4 + Num8_p2 + Node73_p2), (Num7_p2 - Node73_p2 + Num6_p2 + Node53_p2 + 4), Coef, Coef_location);
			}
			else if (((i - M_post2_q3 - Num8_p2) > (Node73_p2 + 3)) && ((i - M_post2_q3 - Num8_p2) < (Node74_p2)))// 1dx
			{
				inner_5node(j, i, 1, -(Num7_p2 - 5 + 1), (5 + Num6_p2 - 1) - 2, Coef, Coef_location);
			}
			else if (((i - M_post2_q3 - Num8_p2) > (Node74_p2 + 1)) && ((i - M_post2_q3 - Num8_p2) < (Num7_p2 - 4)))// 1dx
			{
				inner_5node(j, i, 1, -(Num7_p2 - 5 + 1), 5 + Num6_p2 - 1, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, i, 2, -(Num9_p2 - Node_99_p2 + Num8_p2 + Num7_p2 - 1), (1 + Num6_p2 + Node_55_p2), Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!! num 32 !!!
	for (i = (M_post2_q3 + Num8_p2 + Num7_p2); i<(M_post2_q3 + Num8_p2 + Num7_p2 + Num6_p2); i++)// num 32(6) :  16 node!
	{
		switch (i - M_post2_q3 - Num8_p2 - Num7_p2)
		{
		case 0:// 1 lei 2j shang
			kind1_node(j, i, -(Num7_p2 - 4) - 2, -(Num7_p2 - 4) - 1, -(Num7_p2 - 4), 0, 1, (Num6_p2 + N_d0 + 2) - 2,
				(Num6_p2 + N_d0 + 2) - 1, (Num6_p2 + N_d0 + 2), A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case Node61_p2:// inner 3 node
			inner_3node(j, i, -1, 0, (Num6_p2 + N_d0 + 2), 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node61_p2 + 1) :// inner 3 node
			inner_3node(j, i, 0, 1, (Num6_p2 + N_d0 + 2) + 2, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Num6_p2 / 2 - 1) :// 1 lei 2j xia
			kind1_node(j, i, -(Num7_p2 - 4) - 2, -(Num7_p2 - 4) - 2 + 1, -(Num7_p2 - 4) - 2 + 2, -1, 0,
			(Num6_p2 + N_d0 + 2) + 2, (Num6_p2 + N_d0 + 2) + 2 + 1, (Num6_p2 + N_d0 + 2) + 2 + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		case (Num6_p2 / 2) :// 1 lei 2j shang
			kind1_node(j, i, -(5 + Num6_p2 - 1) + 2 - 2, -(5 + Num6_p2 - 1) + 2 - 1, -(5 + Num6_p2 - 1) + 2,
			0, 1, (1 + Node_55_p2 - 4) - 2 - 2, (1 + Node_55_p2 - 4) - 2 - 1, (1 + Node_55_p2 - 4) - 2,
			A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case Node62_p2:// inner 3 node
			inner_3node(j, i, -1, 0, (1 + Node_55_p2 - 4) - 2, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node62_p2 + 1) :// inner 3 node
			inner_3node(j, i, 0, 1, (1 + Node_55_p2 - 4), k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Num6_p2 - 1) :// 1 lei 2j xia
			kind1_node(j, i, -(5 + Num6_p2 - 1), -(5 + Num6_p2 - 1) + 1, -(5 + Num6_p2 - 1) + 2, -1, 0,
			(1 + Node_55_p2 - 4), (1 + Node_55_p2 - 4) + 1, (1 + Node_55_p2 - 4) + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q3 - Num8_p2 - Num7_p2) < Node61_p2)// 1dx
			{
				inner_5node(j, i, 1, -(Num7_p2 - 4), (Num6_p2 + N_d0 + 2), Coef, Coef_location);
			}
			else if (((i - M_post2_q3 - Num8_p2 - Num7_p2) > (Node61_p2 + 1)) && ((i - M_post2_q3 - Num8_p2 - Num7_p2) < (Num6_p2 / 2 - 1)))// 2dx
			{
				inner_5node(j, i, 1, -(Num7_p2 - 4) - 2, (Num6_p2 + N_d0 + 2) + 2, Coef, Coef_location);
			}
			else if (((i - M_post2_q3 - Num8_p2 - Num7_p2) > (Num6_p2 / 2)) && ((i - M_post2_q3 - Num8_p2 - Num7_p2) < Node62_p2))// 2dx
			{
				inner_5node(j, i, 1, -(5 + Num6_p2 - 1) + 2, (1 + Node_55_p2 - 4) - 2, Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_5node(j, i, 1, -(5 + Num6_p2 - 1), (1 + Node_55_p2 - 4), Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!! num 33 !!!
	for (i = (M_post2_q3 + Num8_p2 + Num7_p2 + Num6_p2);
		i<(M_post2_q3 + Num8_p2 + Num7_p2 + Num6_p2 + Node_5_p2 + 2); i++)// circle num 33(5) :  80 node!  part1 !
	{
		switch (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2)
		{
		case 0:
			boundary_node(j, i, 1, -Num9_p2 - Num8_p2 - Num7_p2 - Num6_p2, Num5_p2 + Num4_p2 + Num3_p2 + Num2_p2, Coef, Coef_location);
			break;
		case (N_d0 - 3) ://2 lei shang 1j
			kind2_node(j, i, -Num9_p2 - Num8_p2 - Num7_p2 - Num6_p2, -2, -1, 0, 1, Num5_p2 + Num4_p2 + Num3_p2 + Num2_p2,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (N_d0 + 1) :// 2 lei 2j shang
			kind2_node(j, i, -(Num7_p2 + Num6_p2 + N_d0 - 2), -2, -1, 0, 1, Num5_p2 - N_d0 + 2 + Num4_p2,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node51_p2:// inner 3 node
			inner_3node(j, i, -1, 0, Num5_p2 - N_d0 - 2, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node51_p2 + 1) :// inner 3 node
			inner_3node(j, i, 0, 1, Num5_p2 - N_d0 - 2 + N_c1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Node_5_p2 - 3) ://2 lei xia 2j
			kind2_node(j, i, -(Num7_p2 - Node72_p2 + Num6_p2 + Node_5_p2), -1, 0, 1, 2, Num5_p2 - Node_5_p2 + Num4_p2 + Num3_p2 / 2 - 1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_5_p2 + 1) ://2 lei xia 1j
			kind2_node(j, i, -(Num9_p2 - Node_9_p2 + Num8_p2 + Num7_p2 + Num6_p2 + Node_5_p2), -1, 0, 1, 2,
			Num5_p2 - Node_5_p2 + Num4_p2 + Num3_p2 + Num2_p2 + Node_1_p2,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2) < (N_d0 - 3))// 4dx!
			{
				inner_5node(j, i, 4, -Num9_p2 - Num8_p2 - Num7_p2 - Num6_p2, Num5_p2 + Num4_p2 + Num3_p2 + Num2_p2, Coef, Coef_location);
			}
			else if (((i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2) > (N_d0 - 3)) && ((i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2) < (N_d0 + 1)))// 2dx!
			{
				inner_5node(j, i, 2, -(Num7_p2 + Num6_p2 + N_d0 - 2), Num5_p2 - N_d0 + 2 + Num4_p2, Coef, Coef_location);
			}
			else if (((i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2) > (N_d0 + 1)) && ((i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2) < Node51_p2))// 1dx
			{
				inner_5node(j, i, 1, -(Num6_p2 + N_d0 + 2), Num5_p2 - N_d0 - 2, Coef, Coef_location);
			}
			else if (((i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2) > (Node51_p2 + 1)) && ((i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2) < (Node_5_p2 - 3)))// 1dx
			{
				inner_5node(j, i, 1, -(Num6_p2 + N_d0 + 2) - 2, Num5_p2 - N_d0 - 2 + N_c1, Coef, Coef_location);
			}
			else// 2dx!
			{
				inner_5node(j, i, 2, -(Num7_p2 - Node72_p2 + Num6_p2 + Node_5_p2), Num5_p2 - Node_5_p2 + Num4_p2 + Num3_p2 / 2 - 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post2_q3 + Num8_p2 + Num7_p2 + Num6_p2 + Node_5_p2 + 2);
		i<(M_post2_q3 + Num8_p2 + Num7_p2 + Num6_p2 + Node53_p2 + 3); i++)// circle num 33(5) :  80 node!  part2 !
	{
		switch (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2)
		{
		default:
			if ((i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2) < Node52_p2)// inner 5 node 4dx
			{
				inner_5node(j, i, 4, -(Num9_p2 - Node_9_p2 + Num8_p2 + Num7_p2 + Num6_p2 + Node_5_p2), Num5_p2 - Node_5_p2 + Num4_p2 + Num3_p2 + Num2_p2 + Node_1_p2, Coef, Coef_location);
			}
			else if ((i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2) > Node53_p2)// inner 5 node 4dx
			{
				inner_5node(j, i, 4, -(Num9_p2 - Node93_p2 + Num8_p2 + Num7_p2 + Num6_p2 + Node53_p2),
					(Num5_p2 - Node_5_p2 + Num4_p2 + Num3_p2 + Num2_p2 + Node_1_p2), Coef, Coef_location);
			}
			else// left  2lei 1j
			{
				kind2_node(j, i, -(Num7_p2 - Node72_p2 - 1 - (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Node52_p2) * 2 + Num6_p2 + (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2)),
					-1, 0, 1, (Num5_p2 - Node_5_p2 + Num4_p2 + Num3_p2 + Num2_p2 + Node_1_p2),
					(Num5_p2 - Node_5_p2 + Num4_p2 + Num3_p2 + Num2_p2 + Node_1_p2) + (Nx_2 + Nx_3),
					B2, B3, B1_1, B3, B4, B5, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post2_q3 + Num8_p2 + Num7_p2 + Num6_p2 + Node53_p2 + 3);
		i<(M_post2_q3 + Num8_p2 + Num7_p2 + Num6_p2 + Num5_p2); i++)// circle num 33(5) :  80 node!  part3 !
	{
		switch (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2)
		{
		case (Num5_p2 - 1) :
			boundary_node(j, i, 2, -Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2, Num4_p2 + Num3_p2 + Num2_p2 + Num1_p2, Coef, Coef_location);
			break;
		case (Node53_p2 + 3) ://2 lei shang 1j
			kind2_node(j, i, -(Num9_p2 - Node93_p2 + Num8_p2 + Num7_p2 + Num6_p2 + Node53_p2), -2, -1, 0, 1, (Num5_p2 - Node_5_p2 + Num4_p2 + Num3_p2 + Num2_p2 + Node_1_p2),
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node53_p2 + 7) :// 2 lei 2j shang
			kind2_node(j, i, -(Num7_p2 - Node73_p2 + Num6_p2 + Node53_p2 + 4), -2, -1, 0, 1, Num5_p2 - Node53_p2 - 4 + Num4_p2 + Num3_p2 / 2,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node54_p2:// inner 3 node
			inner_3node(j, i, -1, 0, (Num5_p2 - Node_55_p2 + 4 + Num4_p2 - 1) - N_c1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node54_p2 + 1) :// inner 3 node
			inner_3node(j, i, 0, 1, (Num5_p2 - Node_55_p2 + 4 + Num4_p2 - 1), k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			break;
		case (Node_55_p2 - 3) ://2 lei xia 2j
			kind2_node(j, i, -(Num6_p2 + Node_55_p2 + 1), -1, 0, 1, 2, Num5_p2 - Node_55_p2 + Num4_p2 + Num3_p2 - 1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_55_p2 + 1) ://2 lei xia 1j
			kind2_node(j, i, -Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2, -1, 0, 1, 2, Num4_p2 + Num3_p2 + Num2_p2 + Num1_p2,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2) < (Node53_p2 + 7))// 2dx!
			{
				inner_5node(j, i, 2, -(Num7_p2 - Node73_p2 + Num6_p2 + Node53_p2 + 4), Num5_p2 - Node53_p2 - 4 + Num4_p2 + Num3_p2 / 2, Coef, Coef_location);
			}
			else if (((i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2) > (Node53_p2 + 7)) && ((i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2) < Node54_p2))// 1dx
			{
				inner_5node(j, i, 1, -(Node_55_p2 - 4 + 1) + 2, (Num5_p2 - Node_55_p2 + 4 + Num4_p2 - 1) - N_c1, Coef, Coef_location);
			}
			else if (((i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2) > (Node54_p2 + 1)) && ((i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2) < (Node_55_p2 - 3)))// 1dx
			{
				inner_5node(j, i, 1, -(Node_55_p2 - 4 + 1), (Num5_p2 - Node_55_p2 + 4 + Num4_p2 - 1), Coef, Coef_location);
			}
			else if (((i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2) > (Node_55_p2 - 3)) && ((i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2) < (Node_55_p2 + 1)))// 2dx
			{
				inner_5node(j, i, 2, -(Num6_p2 + Node_55_p2 + 1), Num5_p2 - Node_55_p2 + Num4_p2 + Num3_p2 - 1, Coef, Coef_location);
			}
			else// 4dx!
			{
				inner_5node(j, i, 4, -Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2, Num4_p2 + Num3_p2 + Num2_p2 + Num1_p2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!! num 34 !!!
	for (i = (M_post2_q3 + Num8_p2 + Num7_p2 + Num6_p2 + Num5_p2);
		i<(M_post2_q3 + Num8_p2 + Num7_p2 + Num6_p2 + Num5_p2 + Num4_p2 / 2); i++)// num 34(4) :  30 node! part1
	{
		switch (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2)
		{
		case 0:// 1 lei 2j shang
			kind1_node(j, i, -(Num5_p2 - N_d0 - 2) - 2, -(Num5_p2 - N_d0 - 2) - 1, -(Num5_p2 - N_d0 - 2), 0, 1, Num4_p2 + 4 - 2,
				Num4_p2 + 4 - 1, Num4_p2 + 4, A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case (Num4_p2 / 2 - 1) :// 1 lei 2j xia
			kind1_node(j, i, -(Num5_p2 - N_d0 - 2) - N_c1, -(Num5_p2 - N_d0 - 2) - N_c1 + 1, -(Num5_p2 - N_d0 - 2) - N_c1 + 2,
			-1, 0, Num4_p2 + 4, Num4_p2 + 4 + 1, Num4_p2 + 4 + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2) <= Node41_p2)// 1dx
			{
				inner_5node(j, i, 1, -(Num5_p2 - N_d0 - 2), Num4_p2 + 4, Coef, Coef_location);
			}
			else if ((i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2) >= Node42_p2) // 1dx
			{
				inner_5node(j, i, 1, -(Num5_p2 - N_d0 - 2) - N_c1, Num4_p2 + 4, Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_4node(j, i, -1, 0, 1, Num4_p2 + 4, 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post2_q3 + Num8_p2 + Num7_p2 + Num6_p2 + Num5_p2 + Num4_p2 / 2);
		i<(M_post2_q3 + Num8_p2 + Num7_p2 + Num6_p2 + Num5_p2 + Num4_p2); i++)// num 34(4) :  30 node! part2
	{
		switch (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2)
		{
		case (Num4_p2 / 2) :// 1 lei 2j shang
			kind1_node(j, i, -(Num5_p2 - Node_55_p2 + 4 + Num4_p2 - 1) + N_c1 - 2, -(Num5_p2 - Node_55_p2 + 4 + Num4_p2 - 1) + N_c1 - 1,
			-(Num5_p2 - Node_55_p2 + 4 + Num4_p2 - 1) + N_c1, 0, 1, (Num3_p2 - 5 + 1) - 2, (Num3_p2 - 5 + 1) - 1, (Num3_p2 - 5 + 1),
			A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case (Num4_p2 - 1) :// 1 lei 2j xia
			kind1_node(j, i, -(Num5_p2 - Node_55_p2 + 4 + Num4_p2 - 1), -(Num5_p2 - Node_55_p2 + 4 + Num4_p2 - 1) + 1,
			-(Num5_p2 - Node_55_p2 + 4 + Num4_p2 - 1) + 2, -1, 0, (Num3_p2 - 5 + 1), (Num3_p2 - 5 + 1) + 1, (Num3_p2 - 5 + 1) + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2) <= Node43_p2)// 1dx
			{
				inner_5node(j, i, 1, -(Num5_p2 - Node_55_p2 + 4 + Num4_p2 - 1) + N_c1, (Num3_p2 - 5 + 1), Coef, Coef_location);
			}
			else if ((i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2) >= Node44_p2) // 1dx
			{
				inner_5node(j, i, 1, -(Num5_p2 - Node_55_p2 + 4 + Num4_p2 - 1), (Num3_p2 - 5 + 1), Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_4node(j, i, -1, 0, 1, (Num3_p2 - 5 + 1), 1, k*k*dx*dx - 4, 1, 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!! num 35 !!!
	for (i = (M_post2_q3 + Num8_p2 + Num7_p2 + Num6_p2 + Num5_p2 + Num4_p2);
		i<(M_post2_q3 + Num8_p2 + Num7_p2 + Num6_p2 + Num5_p2 + Num4_p2 + Num3_p2 / 2); i++)// num 35(3) :  46 node! part1
	{
		switch (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2)
		{
		case 0:// 1 lei 1j shang
			kind1_node(j, i, -(Num5_p2 - N_d0 + 2 + Num4_p2) - 2, -(Num5_p2 - N_d0 + 2 + Num4_p2) - 1, -(Num5_p2 - N_d0 + 2 + Num4_p2),
				0, 1, Num3_p2 + Num2_p2 + N_d0 - 2 - 2, Num3_p2 + Num2_p2 + N_d0 - 2 - 1, Num3_p2 + Num2_p2 + N_d0 - 2,
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case 3:// 2 lei 2j shang
			kind2_node(j, i, -(Num5_p2 - N_d0 + 2 + Num4_p2), -2, -1, 0, 1, Num3_p2 + Num2_p2 + N_d0 - 2,
				B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case 4:// 3 lei 2j shang
			kind3_node(j, i, -(Num4_p2 + 4), -1, 0, 1, (Num3_p2 + Num2_p2 + N_d0 - 2) - 1, (Num3_p2 + Num2_p2 + N_d0 - 2),
				C2_2, C3_2, C1_2, C3_2, C4, C4, Coef, Coef_location);
			break;
		case (Num3_p2 / 2 - 5) :// 3 lei 2j xia
			kind3_node(j, i, -(Num4_p2 + 4), -1, 0, 1, Num3_p2 / 2 + 1 + Num2_p2 + Node_1_p2, Num3_p2 / 2 + 1 + Num2_p2 + Node_1_p2 + 1,
			C2_2, C3_2, C1_2, C3_2, C4, C4, Coef, Coef_location);
			break;
		case (Num3_p2 / 2 - 4) :// 2 lei 2j xia
			kind2_node(j, i, -(Num5_p2 - Node_5_p2 + Num4_p2 + Num3_p2 / 2 - 1), -1, 0, 1, 2, Num3_p2 / 2 + 1 + Num2_p2 + Node_1_p2,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num3_p2 / 2 - 1) :// 1 lei 1j xia
			kind1_node(j, i, -(Num5_p2 - Node_5_p2 + Num4_p2 + Num3_p2 / 2 - 1), -(Num5_p2 - Node_5_p2 + Num4_p2 + Num3_p2 / 2 - 1) + 1, -(Num5_p2 - Node_5_p2 + Num4_p2 + Num3_p2 / 2 - 1) + 2,
			-1, 0, Num3_p2 / 2 + 1 + Num2_p2 + Node_1_p2, Num3_p2 / 2 + 1 + Num2_p2 + Node_1_p2 + 1, Num3_p2 / 2 + 1 + Num2_p2 + Node_1_p2 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2) < 3)// 2dx
			{
				inner_5node(j, i, 2, -(Num5_p2 - N_d0 + 2 + Num4_p2), Num3_p2 + Num2_p2 + N_d0 - 2, Coef, Coef_location);
			}
			else if (((i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2) > 4) && ((i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2) < (Num3_p2 / 2 - 5)))// 1dx
			{
				inner_5node(j, i, 1, -(Num4_p2 + 4), Num3_p2 - 5, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, i, 2, -(Num5_p2 - Node_5_p2 + Num4_p2 + Num3_p2 / 2 - 1), Num3_p2 / 2 + 1 + Num2_p2 + Node_1_p2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post2_q3 + Num8_p2 + Num7_p2 + Num6_p2 + Num5_p2 + Num4_p2 + Num3_p2 / 2);
		i<(M_post2_q3 + Num8_p2 + Num7_p2 + Num6_p2 + Num5_p2 + Num4_p2 + Num3_p2); i++)// num 35(3) :  46 node! part2
	{
		switch (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2 - Num3_p2 / 2)
		{
		case 0:// 1 lei 1j shang
			kind1_node(j, i, -(Num5_p2 - Node_5_p2 - N_a + 3 + Num4_p2 + Num3_p2 / 2) - 2, -(Num5_p2 - Node_5_p2 - N_a + 3 + Num4_p2 + Num3_p2 / 2) - 1,
				-(Num5_p2 - Node_5_p2 - N_a + 3 + Num4_p2 + Num3_p2 / 2), 0, 1,
				Num3_p2 / 2 + Num2_p2 + Node_1_p2 + N_a - 3 - 2, Num3_p2 / 2 + Num2_p2 + Node_1_p2 + N_a - 3 - 1, Num3_p2 / 2 + Num2_p2 + Node_1_p2 + N_a - 3,
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case 3:// 2 lei 2j shang
			kind2_node(j, i, -(Num5_p2 - Node_5_p2 - N_a + 3 + Num4_p2 + Num3_p2 / 2), -2, -1, 0, 1, Num3_p2 / 2 + Num2_p2 + Node_1_p2 + N_a - 3,
				B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case 4:// 3 lei 2j shang
			kind3_node(j, i, -(Num3_p2 - 5 + 1), -1, 0, 1, (Num3_p2 / 2 + Num2_p2 + Node_1_p2 + N_a - 3) - 1, (Num3_p2 / 2 + Num2_p2 + Node_1_p2 + N_a - 3),
				C2_2, C3_2, C1_2, C3_2, C4, C4, Coef, Coef_location);
			break;
		case (Num3_p2 / 2 - 5) :// 3 lei 2j xia
			kind3_node(j, i, -(Num3_p2 - 5 + 1), -1, 0, 1, 1 + Num2_p2 + Node_11_p2, 1 + Num2_p2 + Node_11_p2 + 1,
			C2_2, C3_2, C1_2, C3_2, C4, C4, Coef, Coef_location);
			break;
		case (Num3_p2 / 2 - 4) :// 2 lei 2j xia
			kind2_node(j, i, -(Num5_p2 - Node_55_p2 + Num4_p2 + Num3_p2 - 1), -1, 0, 1, 2, 1 + Num2_p2 + Node_11_p2,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num3_p2 / 2 - 1) :// 1 lei 1j xia
			kind1_node(j, i, -(Num5_p2 - Node_55_p2 + Num4_p2 + Num3_p2 - 1), -(Num5_p2 - Node_55_p2 + Num4_p2 + Num3_p2 - 1) + 1,
			-(Num5_p2 - Node_55_p2 + Num4_p2 + Num3_p2 - 1) + 2, -1, 0, 1 + Num2_p2 + Node_11_p2, 1 + Num2_p2 + Node_11_p2 + 1, 1 + Num2_p2 + Node_11_p2 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2 - Num3_p2 / 2) < 3)// 2dx
			{
				inner_5node(j, i, 2, -(Num5_p2 - Node_5_p2 - N_a + 3 + Num4_p2 + Num3_p2 / 2), Num3_p2 / 2 + Num2_p2 + Node_1_p2 + N_a - 3, Coef, Coef_location);
			}
			else if (((i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2 - Num3_p2 / 2) > 4) && ((i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2 - Num3_p2 / 2) < (Num3_p2 / 2 - 5)))// 1dx
			{
				inner_5node(j, i, 1, -(Num3_p2 - 5 + 1), 6 + Num2_p2 - 1, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, i, 2, -(Num5_p2 - Node_55_p2 + Num4_p2 + Num3_p2 - 1), 1 + Num2_p2 + Node_11_p2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!! num 36 !!!
	for (i = (M_post2_q3 + Num8_p2 + Num7_p2 + Num6_p2 + Num5_p2 + Num4_p2 + Num3_p2);
		i<(M_post2_q3 + Num8_p2 + Num7_p2 + Num6_p2 + Num5_p2 + Num4_p2 + Num3_p2 + Num2_p2 / 2); i++)// num 36(2) :  26 node! part1
	{
		switch (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2 - Num3_p2)
		{
		case 0:// 3 lei 2j shang
			kind3_node(j, i, -(Num3_p2 - 3), -(Num3_p2 - 3) + 2, 0, 1, Num2_p2 + Node11_p2 - 1, Num2_p2 + Node11_p2,
				C4, C3_2, C1_2, C2_2, C4, C3_2, Coef, Coef_location);
			break;
		case (Num2_p2 / 2 - 1) ://node2: 3 lei 2j xia
			kind3_node(j, i, -(Num3_p2 - 5), -(Num3_p2 - 5) + 2, -1, 0, Num2_p2 / 2 + Node_1_p2 - 4 + 1, Num2_p2 / 2 + Node_1_p2 - 4 + 1 + 1,
			C3_2, C4, C2_2, C1_2, C3_2, C4, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2 - Num3_p2) % 2 == 0)// 1dx
			{
				inner_5node(j, i, 1, -(Num3_p2 - 5),
					(Num2_p2 - (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2 - Num3_p2) + Node11_p2 + (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2 - Num3_p2) / 2), Coef, Coef_location);
			}
			else//1 lei  2j left
			{
				kind1_node(j, i, -(Num3_p2 - 5), -1, 0, 1,
					(Num2_p2 - (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2 - Num3_p2) + Node11_p2 + (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2 - Num3_p2) / 2),
					(Num2_p2 - (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2 - Num3_p2) + Node11_p2 + (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2 - Num3_p2) / 2) + 1,
					(Num2_p2 - (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2 - Num3_p2) + Node11_p2 + (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2 - Num3_p2) / 2) + (Num1_p2 - N_d0 + 2),
					(Num2_p2 - (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2 - Num3_p2) + Node11_p2 + (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2 - Num3_p2) / 2) + 1 + (Num1_p2 - N_d0 + 2),
					A4_2, A3_2, A5_2, A3_2, A1_2, A1_2, A2_2, A2_2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (M_post2_q3 + Num8_p2 + Num7_p2 + Num6_p2 + Num5_p2 + Num4_p2 + Num3_p2 + Num2_p2 / 2);
		i<(M_post2_q3 + Num8_p2 + Num7_p2 + Num6_p2 + Num5_p2 + Num4_p2 + Num3_p2 + Num2_p2); i++)// num 36(2) :  26 node! part2
	{
		switch (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2 - Num3_p2 - Num2_p2 / 2)
		{
		case 0:// 3 lei 2j shang
			kind3_node(j, i, -(6 + Num2_p2 - 1) - 2, -(6 + Num2_p2 - 1), 0, 1, Num2_p2 / 2 + Node13_p2 - 1, Num2_p2 / 2 + Node13_p2,
				C4, C3_2, C1_2, C2_2, C4, C3_2, Coef, Coef_location);
			break;
		case (Num2_p2 / 2 - 1) :// 3 lei 2j xia
			kind3_node(j, i, -(6 + Num2_p2 - 1), -(6 + Num2_p2 - 1) + 2, -1, 0, Node_11_p2 - 4 + 1, Node_11_p2 - 4 + 1 + 1,
			C3_2, C4, C2_2, C1_2, C3_2, C4, Coef, Coef_location);
			break;
		default:
			if ((i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2 - Num3_p2 - Num2_p2 / 2) % 2 == 0)// 1dx
			{
				inner_5node(j, i, 1, -(6 + Num2_p2 - 1), Num2_p2 / 2 - (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2 - Num3_p2 - Num2_p2 / 2) + Node13_p2 + (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2 - Num3_p2 - Num2_p2 / 2) / 2,
					Coef, Coef_location);
			}
			else//1 lei 2j left
			{
				kind1_node(j, i, -(6 + Num2_p2 - 1), -1, 0, 1,
					(Num2_p2 / 2 - (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2 - Num3_p2 - Num2_p2 / 2) + Node13_p2 + (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2 - Num3_p2 - Num2_p2 / 2) / 2),
					(Num2_p2 / 2 - (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2 - Num3_p2 - Num2_p2 / 2) + Node13_p2 + (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2 - Num3_p2 - Num2_p2 / 2) / 2) + 1,
					(Num2_p2 / 2 - (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2 - Num3_p2 - Num2_p2 / 2) + Node13_p2 + (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2 - Num3_p2 - Num2_p2 / 2) / 2) + (Num1_p2 - Node_11_p2 + Nx_3 - 1),
					(Num2_p2 / 2 - (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2 - Num3_p2 - Num2_p2 / 2) + Node13_p2 + (i - M_post2_q3 - Num8_p2 - Num7_p2 - Num6_p2 - Num5_p2 - Num4_p2 - Num3_p2 - Num2_p2 / 2) / 2) + 1 + (Num1_p2 - Node_11_p2 + Nx_3 - 1),
					A4_2, A3_2, A5_2, A3_2, A1_2, A1_2, A2_2, A2_2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!! num 37 !!!
	for (i = (num_total - Nx_2 - Nx_3 - Num1_p2); i<(num_total - Nx_2 - Nx_3 - Num1_p2 + Node_1_p2 + 2); i++)// circle num 37(1) Nx_2:  74 node!  part1 !
	{
		switch (i - num_total + Nx_2 + Nx_3 + Num1_p2)
		{
		case 0:
			boundary_node(j, i, 1, -Num5_p2 - Num4_p2 - Num3_p2 - Num2_p2, Nx_2 + Nx_3, Coef, Coef_location);
			break;
		case (N_d0 - 3) ://2 lei shang 1j
			kind2_node(j, i, -Num5_p2 - Num4_p2 - Num3_p2 - Num2_p2, -2, -1, 0, 1, Nx_2 + Nx_3,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_1_p2 + 1) ://2 lei xia 1j
			kind2_node(j, i, -(Num5_p2 - Node_5_p2 + Num4_p2 + Num3_p2 + Num2_p2 + Node_1_p2), -1, 0, 1, 2, Nx_2 + Nx_3,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - num_total + Nx_2 + Nx_3 + Num1_p2) < (N_d0 - 3))// 4dx!
			{
				inner_5node(j, i, 4, -Num5_p2 - Num4_p2 - Num3_p2 - Num2_p2, Nx_2 + Nx_3, Coef, Coef_location);
			}
			else if (((i - num_total + Nx_2 + Nx_3 + Num1_p2) > (N_d0 - 3)) && ((i - num_total + Nx_2 + Nx_3 + Num1_p2) < Node11_p2))// 2dx!
			{
				inner_5node(j, i, 2, -(Num3_p2 + Num2_p2 + N_d0 - 2), Num1_p2 - N_d0 + 2, Coef, Coef_location);
			}
			else if (((i - num_total + Nx_2 + Nx_3 + Num1_p2) >= Node11_p2) && ((i - num_total + Nx_2 + Nx_3 + Num1_p2) <= Node12_p2))// 2 lei right 2j !
			{
				kind2_node(j, i, -((i - num_total + Nx_2 + Nx_3 + Num1_p2) + Num2_p2 - (i - num_total + Nx_2 + Nx_3 + Num1_p2 - Node11_p2) * 2),
					-1, 0, 1, Num1_p2 - N_d0 + 2, Nx_2 + Nx_3, B2, B3, B1_2, B3, B4, B5, Coef, Coef_location);
			}
			else// 2dx!
			{
				inner_5node(j, i, 2, -(Num3_p2 / 2 + 1 + Num2_p2 + Node_1_p2), Num1_p2 - N_d0 + 2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (num_total - Nx_2 - Nx_3 - Num1_p2 + Node_1_p2 + 2); i<(num_total - Nx_2 - Nx_3); i++)// circle num 37(1) Nx_2:  74 node!  part2 !
	{
		switch (i - num_total + Nx_2 + Nx_3 + Num1_p2)
		{
		case (Num1_p2 - 1) :
			boundary_node(j, i, 2, -Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2, Nx_2 + Nx_3, Coef, Coef_location);
			break;
		case (Node_1_p2 + N_a - 4) ://2 lei shang 1j
			kind2_node(j, i, -(Num5_p2 - Node_5_p2 + Num4_p2 + Num3_p2 + Num2_p2 + Node_1_p2), -2, -1, 0, 1, Nx_2 + Nx_3,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_11_p2 + 1) ://2 lei xia 1j
			kind2_node(j, i, -Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2, -1, 0, 1, 2, Nx_2 + Nx_3,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - num_total + Nx_2 + Nx_3 + Num1_p2) < (Node_1_p2 + N_a - 4))// 4dx!
			{
				inner_5node(j, i, 4, -(Num5_p2 - Node_5_p2 + Num4_p2 + Num3_p2 + Num2_p2 + Node_1_p2), Nx_2 + Nx_3, Coef, Coef_location);
			}
			else if (((i - num_total + Nx_2 + Nx_3 + Num1_p2) > (Node_1_p2 + N_a - 4)) && ((i - num_total + Nx_2 + Nx_3 + Num1_p2) < Node13_p2))// 2dx!
			{
				inner_5node(j, i, 2, -(Num3_p2 / 2 + Num2_p2 + Node_1_p2 + N_a - 3), Num1_p2 - Node_11_p2 + Nx_3 - 1, Coef, Coef_location);
			}
			else if (((i - num_total + Nx_2 + Nx_3 + Num1_p2) >= Node13_p2) && ((i - num_total + Nx_2 + Nx_3 + Num1_p2) <= Node14_p2))// 2 lei left 2j !
			{
				kind2_node(j, i, -((i - num_total + Nx_2 + Nx_3 + Num1_p2) + Num2_p2 / 2 - (i - num_total + Nx_2 + Nx_3 + Num1_p2 - Node13_p2) * 2),
					-1, 0, 1, Num1_p2 - Node_11_p2 + Nx_3 - 1, Nx_2 + Nx_3, B2, B3, B1_2, B3, B4, B5, Coef, Coef_location);
			}
			else if (((i - num_total + Nx_2 + Nx_3 + Num1_p2) >= Node14_p2) && ((i - num_total + Nx_2 + Nx_3 + Num1_p2) <= Node_11_p2))// 2dx
			{
				inner_5node(j, i, 2, -(Num2_p2 + Node_11_p2 + 1), Num1_p2 - Node_11_p2 + Nx_3 - 1, Coef, Coef_location);
			}
			else
			{
				inner_5node(j, i, 4, -Num1_p2 - Num2_p2 - Num3_p2 - Num4_p2, Nx_2 + Nx_3, Coef, Coef_location);// 4dx!
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	for (i = (num_total - Nx_2 - Nx_3); i<(num_total - Nx_2); i++)//right1 Nx_3:  30 node! 
	{
		switch (i - num_total + Nx_2 + Nx_3)
		{
		case 0:// 1 lei shang
			kind1_node(j, i, -(Nx_2 - N_d0 + 2) - 2, -(Nx_2 - N_d0 + 2) - 1, -(Nx_2 - N_d0 + 2), 0, 1, Nx_3 + N_d0 - 2 - 2, Nx_3 + N_d0 - 2 - 1, Nx_3 + N_d0 - 2,
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case (Nx_3 / 2 - 1) :// 1 lei xia
			kind1_node(j, i, -(Nx_2 - N_d0 + 2), -(Nx_2 - N_d0 + 2) + 1, -(Nx_2 - N_d0 + 2) + 2, -1, 0,
			Nx_3 + N_d0 - 2, Nx_3 + N_d0 - 2 + 1, Nx_3 + N_d0 - 2 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		case (Nx_3 / 2) :// 1 lei shang
			kind1_node(j, i, -(Nx_2 - Node_11_p2 - 1 + Nx_3) - 2, -(Nx_2 - Node_11_p2 - 1 + Nx_3) - 1, -(Nx_2 - Node_11_p2 - 1 + Nx_3), 0, 1,
			1 + Node_11_p2 - 2, 1 + Node_11_p2 - 1, 1 + Node_11_p2,
			A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case (Nx_3 - 1) :// 1 lei xia
			kind1_node(j, i, -(Nx_2 - Node_11_p2 - 1 + Nx_3), -(Nx_2 - Node_11_p2 - 1 + Nx_3) + 1, -(Nx_2 - Node_11_p2 - 1 + Nx_3) + 2, -1, 0,
			1 + Node_11_p2, 1 + Node_11_p2 + 1, 1 + Node_11_p2 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - num_total + Nx_2 + Nx_3) < (Nx_3 / 2 - 1))
			{
				inner_5node(j, i, 2, -(Nx_2 - N_d0 + 2), Nx_3 + N_d0 - 2, Coef, Coef_location);
			}
			else
			{
				inner_5node(j, i, 2, -(Nx_2 - Node_11_p2 - 1 + Nx_3), 1 + Node_11_p2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//***************************************************************************************
	//最后一行 本节点向下的区域没有影响
	//cout << "circle num1:" << j << endl;
	for (i = (Nx_2 + Nx_3); i<(Nx_2 + Nx_3 + Node_1 + 2); i++)// circle num1 Nx_2:  74 node!  part1 !
	{
		switch (i - Nx_2 - Nx_3)
		{
		case 0:
			boundary_node(j, j - offset, 1, -Nx_2 - Nx_3, Num1 + Num2 + Num3 + Num4, Coef, Coef_location);
			break;
		case (N_d0 - 3) ://2 lei shang 1j
			kind2_node(j, j - offset, -Nx_2 - Nx_3, -2, -1, 0, 1, Num1 + Num2 + Num3 + Num4,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_1 + 1) ://2 lei xia 1j
			kind2_node(j, j - offset, -(Nx_2 + Nx_3), -1, 0, 1, 2, Num1 - Node_1 + Num2 + Num3 + Num4 + Node_5,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3) < (N_d0 - 3))// 4dx!
			{
				inner_5node(j, j - offset, 4, -Nx_3 - Nx_2, Num1 + Num2 + Num3 + Num4, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3) > (N_d0 - 3)) && ((i - Nx_2 - Nx_3) < Node11))// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Nx_3 + N_d0 - 2), Num1 - N_d0 + 2 + Num2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3) >= Node11) && ((i - Nx_2 - Nx_3) <= Node12))// 2 lei left 2j !
			{
				kind2_node(j, j - offset, -Nx_3 - Nx_2, -(Nx_3 + N_d0 - 2), -1, 0, 1, Num1 - (i - Nx_2 - Nx_3) + ((i - Nx_2 - Nx_3) - Node11) * 2,
					B5, B4, B3, B1_2, B3, B2, Coef, Coef_location);
			}
			else// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Nx_3 + N_d0 - 2), Num1 - Node_1 + Num2 + Num3 / 2 - 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Node_1 + 2); i<(Nx_2 + Nx_3 + Num1); i++)// circle num1 Nx_2:  74 node!  part2 !
	{
		switch (i - Nx_2 - Nx_3)
		{
		case (Num1 - 1) :
			boundary_node(j, j - offset, 2, -Nx_2 - Nx_3, Num5 + Num2 + Num3 + Num4, Coef, Coef_location);
			break;
		case (Node_1 + N_a - 4) ://2 lei shang 1j
			kind2_node(j, j - offset, -Nx_2 - Nx_3, -2, -1, 0, 1, Num1 - Node_1 + Num2 + Num3 + Num4 + Node_5,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_11 + 1) ://2 lei xia 1j
			kind2_node(j, j - offset, -(Nx_2 + Nx_3), -1, 0, 1, 2, Num5 + Num2 + Num3 + Num4,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3) < (Node_1 + N_a - 4))// 4dx!
			{
				inner_5node(j, j - offset, 4, -Nx_3 - Nx_2, Num1 - Node_1 + Num2 + Num3 + Num4 + Node_5, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3) > (Node_1 + N_a - 4)) && ((i - Nx_2 - Nx_3) < Node13))// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Node_11 + 1), Num1 - Node_1 - N_a + 3 + Num2 + Num3 / 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3) >= Node13) && ((i - Nx_2 - Nx_3) <= Node14))// 2 lei left 2j !
			{
				kind2_node(j, j - offset, -Nx_3 - Nx_2, -(Node_11 + 1), -1, 0, 1, Num1 - (i - Nx_2 - Nx_3) + Num2 / 2 + ((i - Nx_2 - Nx_3) - Node13) * 2,
					B5, B4, B3, B1_2, B3, B2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3) >= Node14) && ((i - Nx_2 - Nx_3) <= Node_11))// 2dx
			{
				inner_5node(j, j - offset, 2, -(Node_11 + 1), Num1 - Node_11 + Num2 + Num3 - 1, Coef, Coef_location);
			}
			else
			{
				inner_5node(j, j - offset, 4, -(Nx_3 + Nx_2), Num5 + Num2 + Num3 + Num4, Coef, Coef_location);// 4dx!
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!!!!!!!!!
	//cout << "circle num 2:" << j << endl;
	for (i = (Nx_2 + Nx_3 + Num1); i<(Nx_2 + Nx_3 + Num1 + Num2 / 2); i++)// num 2 :  26 node! part1
	{
		switch (i - Nx_2 - Nx_3 - Num1)
		{
		case 0:// 3 lei 2j shang
			kind3_node(j, j - offset, -(Num1 - Node11) - 1, -(Num1 - Node11), 0, 1, Num2 + 5 - 2, Num2 + 5,
				C4, C3_2, C1_2, C2_2, C4, C3_2, Coef, Coef_location);
			break;
		case (Num2 / 2 - 1) ://node2: 3 lei 2j xia
			kind3_node(j, j - offset, -(Num1 - Node12 + Num2 / 2 - 1), -(Num1 - Node12 + Num2 / 2 - 1) + 1, -1, 0, Num2 + 5, Num2 + 5 + 2,
			C3_2, C4, C2_2, C1_2, C3_2, C4, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1) % 2 == 0)// 1dx
			{
				inner_5node(j, j - offset, 1, -((i - Nx_2 - Nx_3 - Num1) + Num1 - Node11 - (i - Nx_2 - Nx_3 - Num1) / 2), Num2 + 5, Coef, Coef_location);
			}
			else//1 lei  2j left
			{
				kind1_node(j, j - offset, -((i - Nx_2 - Nx_3 - Num1) + Num1 - Node11 - (i - Nx_2 - Nx_3 - Num1) / 2) - (Nx_3 + N_d0 - 2),
					-((i - Nx_2 - Nx_3 - Num1) + Num1 - Node11 - (i - Nx_2 - Nx_3 - Num1) / 2) - (Nx_3 + N_d0 - 2) + 1,
					-((i - Nx_2 - Nx_3 - Num1) + Num1 - Node11 - (i - Nx_2 - Nx_3 - Num1) / 2),
					-((i - Nx_2 - Nx_3 - Num1) + Num1 - Node11 - (i - Nx_2 - Nx_3 - Num1) / 2) + 1,
					-1, 0, 1, 5 + Num2, A2_2, A2_2, A1_2, A1_2, A3_2, A5_2, A3_2, A4_2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 / 2); i<(Nx_2 + Nx_3 + Num1 + Num2); i++)// num 2 :  26 node! part2
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 / 2)
		{
		case 0:// 3 lei 2j shang
			kind3_node(j, j - offset, -(Num1 - Node13 + Num2 / 2) - 1, -(Num1 - Node13 + Num2 / 2), 0, 1, Num3 - 5 - 2, Num3 - 5,
				C4, C3_2, C1_2, C2_2, C4, C3_2, Coef, Coef_location);
			break;
		case (Num2 / 2 - 1) :// 3 lei 2j xia
			kind3_node(j, j - offset, -(Num1 - Node14 + Num2 - 1), -(Num1 - Node14 + Num2 - 1) + 1, -1, 0, Num3 - 5, Num3 - 5 + 2,
			C3_2, C4, C2_2, C1_2, C3_2, C4, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 / 2) % 2 == 0)// 1dx
			{
				inner_5node(j, j - offset, 1, -((i - Nx_2 - Nx_3 - Num1 - Num2 / 2) + Num1 - Node13 - (i - Nx_2 - Nx_3 - Num1 - Num2 / 2) / 2) - Num2 / 2,
					Num3 - 5, Coef, Coef_location);
			}
			else//1 lei 2j left
			{
				kind1_node(j, j - offset, -((i - Nx_2 - Nx_3 - Num1 - Num2 / 2) + Num1 - Node13 - (i - Nx_2 - Nx_3 - Num1 - Num2 / 2) / 2) - Num2 / 2 - (Node_11 + 1),
					-((i - Nx_2 - Nx_3 - Num1 - Num2 / 2) + Num1 - Node13 - (i - Nx_2 - Nx_3 - Num1 - Num2 / 2) / 2) - Num2 / 2 - (Node_11 + 1) + 1,
					-((i - Nx_2 - Nx_3 - Num1 - Num2 / 2) + Num1 - Node13 - (i - Nx_2 - Nx_3 - Num1 - Num2 / 2) / 2) - Num2 / 2,
					-((i - Nx_2 - Nx_3 - Num1 - Num2 / 2) + Num1 - Node13 - (i - Nx_2 - Nx_3 - Num1 - Num2 / 2) / 2) - Num2 / 2 + 1,
					-1, 0, 1, Num3 - 5, A2_2, A2_2, A1_2, A1_2, A3_2, A5_2, A3_2, A4_2, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************!!!!
	//cout << "circle num 3:" << j << endl;
	for (i = (Nx_2 + Nx_3 + Num1 + Num2); i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 / 2); i++)// num 3 :  46 node! part1
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2)
		{
		case 0:// 1 lei 1j shang
			kind1_node(j, j - offset, -(Num1 - N_d0 + 2 + Num2) - 2, -(Num1 - N_d0 + 2 + Num2) - 1, -(Num1 - N_d0 + 2 + Num2), 0, 1, Num3 + Num4 + N_d0 - 2 - 2,
				Num3 + Num4 + N_d0 - 2 - 1, Num3 + Num4 + N_d0 - 2, A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case 3:// 2 lei 2j shang
			kind2_node(j, j - offset, -(Num1 - N_d0 + 2 + Num2), -2, -1, 0, 1, Num3 + Num4 + N_d0 - 2,
				B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case 4:// 3 lei 2j shang
			kind3_node(j, j - offset, -(Num1 - N_d0 + 2 + Num2) - 1, -(Num1 - N_d0 + 2 + Num2), -1, 0, 1, Num3 - 4,
				C4, C4, C3_2, C1_2, C3_2, C2_2, Coef, Coef_location);
			break;
		case (Num3 / 2 - 5) :// 3 lei 2j xia
			kind3_node(j, j - offset, -(Num1 - Node_1 + Num2 + Num3 / 2 - 1), -(Num1 - Node_1 + Num2 + Num3 / 2 - 1) + 1, -1, 0, 1, Num3 - 4,
			C4, C4, C3_2, C1_2, C3_2, C2_2, Coef, Coef_location);
			break;
		case (Num3 / 2 - 4) :// 2 lei 2j xia
			kind2_node(j, j - offset, -(Num1 - Node_1 + Num2 + Num3 / 2 - 1), -1, 0, 1, 2, Num3 / 2 + 1 + Num4 + Node_5,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num3 / 2 - 1) :// 1 lei 1j xia
			kind1_node(j, j - offset, -(Num1 - Node_1 + Num2 + Num3 / 2 - 1), -(Num1 - Node_1 + Num2 + Num3 / 2 - 1) + 1, -(Num1 - Node_1 + Num2 + Num3 / 2 - 1) + 2, -1, 0,
			Num3 / 2 + 1 + Num4 + Node_5, Num3 / 2 + 1 + Num4 + Node_5 + 1, Num3 / 2 + 1 + Num4 + Node_5 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2) < 3)// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num1 - N_d0 + 2 + Num2), Num3 + Num4 + N_d0 - 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2) > 4) && ((i - Nx_2 - Nx_3 - Num1 - Num2) < (Num3 / 2 - 5)))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num2 + 5), Num3 - 4, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num1 - Node_1 + Num2 + Num3 / 2 - 1), Num3 / 2 + 1 + Num4 + Node_5, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 / 2); i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3); i++)// num 3 :  46 node! part2
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 / 2)
		{
		case 0:// 1 lei 1j shang
			kind1_node(j, j - offset, -(Num1 - Node_1 - N_a + 3 + Num2 + Num3 / 2) - 2, -(Num1 - Node_1 - N_a + 3 + Num2 + Num3 / 2) - 1, -(Num1 - Node_1 - N_a + 3 + Num2 + Num3 / 2), 0, 1,
				Num3 / 2 + Num4 + Node_5 + N_a - 3 - 2, Num3 / 2 + Num4 + Node_5 + N_a - 3 - 1, Num3 / 2 + Num4 + Node_5 + N_a - 3,
				A2_1, A1_1, A3_1, A5_1, A4_1, A2_1, A1_1, A3_1, Coef, Coef_location);
			break;
		case 3:// 2 lei 2j shang
			kind2_node(j, j - offset, -(Num1 - Node_1 - N_a + 3 + Num2 + Num3 / 2), -2, -1, 0, 1, Num3 / 2 + Num4 + Node_5 + N_a - 3,
				B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case 4:// 3 lei 2j shang
			kind3_node(j, j - offset, -(Num1 - Node_1 - N_a + 3 + Num2 + Num3 / 2) - 1, -(Num1 - Node_1 - N_a + 3 + Num2 + Num3 / 2), -1, 0, 1, 5 + Num4 - 1,
				C4, C4, C3_2, C1_2, C3_2, C2_2, Coef, Coef_location);
			break;
		case (Num3 / 2 - 5) :// 3 lei 2j xia
			kind3_node(j, j - offset, -(Num1 - Node_11 + Num2 + Num3 - 1), -(Num1 - Node_11 + Num2 + Num3 - 1) + 1, -1, 0, 1, 5 + Num4 - 1,
			C4, C4, C3_2, C1_2, C3_2, C2_2, Coef, Coef_location);
			break;
		case (Num3 / 2 - 4) :// 2 lei 2j xia
			kind2_node(j, j - offset, -(Num1 - Node_11 + Num2 + Num3 - 1), -1, 0, 1, 2, 1 + Num4 + Node_55,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Num3 / 2 - 1) :// 1 lei 1j xia
			kind1_node(j, j - offset, -(Num1 - Node_11 + Num2 + Num3 - 1), -(Num1 - Node_11 + Num2 + Num3 - 1) + 1, -(Num1 - Node_11 + Num2 + Num3 - 1) + 2, -1, 0,
			1 + Num4 + Node_55, 1 + Num4 + Node_55 + 1, 1 + Num4 + Node_55 + 2,
			A3_1, A1_1, A2_1, A4_1, A5_1, A3_1, A1_1, A2_1, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 / 2) < 3)// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num1 - Node_1 - N_a + 3 + Num2 + Num3 / 2), Num3 / 2 + Num4 + Node_5 + N_a - 3, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 / 2) > 4) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 / 2) < (Num3 / 2 - 5)))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num3 - 5), 5 + Num4 - 1, Coef, Coef_location);
			}
			else// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num1 - Node_11 + Num2 + Num3 - 1), 1 + Num4 + Node_55, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	//cout << "circle num 4:" << j << endl;
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3); i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 / 2); i++)// num 4 :  30 node! part1
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3)
		{
		case 0:// 1 lei 2j shang
			kind1_node(j, j - offset, -(Num3 - 4) - 2, -(Num3 - 4) - 1, -(Num3 - 4), 0, 1, Num4 + N_d0 + 2 - 2, Num4 + N_d0 + 2 - 1, Num4 + N_d0 + 2,
				A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case (Num4 / 2 - 1) :// 1 lei 2j xia
			kind1_node(j, j - offset, -(Num3 - 4), -(Num3 - 4) + 1, -(Num3 - 4) + 2, -1, 0,
			Num4 + N_d0 + 2 - N_c1, Num4 + N_d0 + 2 - N_c1 + 1, Num4 + N_d0 + 2 - N_c1 + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3) <= Node41)// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num3 - 4), Num4 + N_d0 + 2, Coef, Coef_location);
			}
			else if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3) >= Node42) // 1dx
			{
				inner_5node(j, j - offset, 1, -(Num3 - 4), Num4 + N_d0 + 2 - N_c1, Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_4node(j, j - offset, -(Num3 - 4), -1, 0, 1, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 / 2); i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4); i++)// num 4 :  30 node! part2
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3)
		{
		case (Num4 / 2) :// 1 lei 2j shang
			kind1_node(j, j - offset, -(5 + Num4 - 1) - 2, -(5 + Num4 - 1) - 1, -(5 + Num4 - 1), 0, 1,
			(Node_55 - 4 + 1 + N_c1) - 2, (Node_55 - 4 + 1 + N_c1) - 1, (Node_55 - 4 + 1 + N_c1),
			A2_2, A1_2, A3_2, A5_2, A4_2, A2_2, A1_2, A3_2, Coef, Coef_location);
			break;
		case (Num4 - 1) :// 1 lei 2j xia
			kind1_node(j, j - offset, -(5 + Num4 - 1), -(5 + Num4 - 1) + 1, -(5 + Num4 - 1) + 2, -1, 0,
			(Node_55 - 4 + 1), (Node_55 - 4 + 1) + 1, (Node_55 - 4 + 1) + 2,
			A3_2, A1_2, A2_2, A4_2, A5_2, A3_2, A1_2, A2_2, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3) <= Node43)// 1dx
			{
				inner_5node(j, j - offset, 1, -(5 + Num4 - 1), (Node_55 - 4 + 1 + N_c1), Coef, Coef_location);
			}
			else if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3) >= Node44) // 1dx
			{
				inner_5node(j, j - offset, 1, -(5 + Num4 - 1), (Node_55 - 4 + 1), Coef, Coef_location);
			}
			else// 1dx 4_node
			{
				inner_4node(j, j - offset, -(5 + Num4 - 1), -1, 0, 1, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//******************************************************************************************
	//cout << "circle num5:" << j << endl;
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4); i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Node_5 + 2); i++)// circle num5 :  80 node!  part1 !
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4)
		{
		case 0:
			boundary_node(j, j - offset, 1, -Num1 - Num2 - Num3 - Num4, Num5 + Num6 + Num7 + Num8, Coef, Coef_location);
			break;
		case (N_d0 - 3) ://2 lei shang 1j
			kind2_node(j, j - offset, -Num1 - Num2 - Num3 - Num4, -2, -1, 0, 1, Num5 + Num6 + Num7 + Num8,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (N_d0 + 1) :// 2 lei 2j shang
			kind2_node(j, j - offset, -(Num3 + Num4 + N_d0 - 2), -2, -1, 0, 1, Num5 - N_d0 + 2 + Num6,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node51:// inner 3 node
			inner_3node(j, j - offset, -(Num4 + N_d0 + 2), -1, 0, 1, 1, k*k*dx*dx - 4, Coef, Coef_location);
			break;
		case (Node51 + 1) :// inner 3 node
			inner_3node(j, j - offset, -(Num4 + N_d0 + 2) + N_c1, 0, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node_5 - 3) ://2 lei xia 2j
			kind2_node(j, j - offset, -(Num3 / 2 + 1 + Num4 + Node_5), -1, 0, 1, 2, Num5 - Node_5 + Num6 + Num7 / 2 - 1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_5 + 1) ://2 lei xia 1j
			kind2_node(j, j - offset, -(Num1 - Node_1 + Num2 + Num3 + Num4 + Node_5), -1, 0, 1, 2, Num5 - Node_5 + Num6 + Num7 + Num8 + Node_9,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < (N_d0 - 3))// 4dx!
			{
				inner_5node(j, j - offset, 4, -Num1 - Num2 - Num3 - Num4, Num5 + Num6 + Num7 + Num8, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) > (N_d0 - 3)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < (N_d0 + 1)))// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Num3 + Num4 + N_d0 - 2), Num5 - N_d0 + 2 + Num6, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) > (N_d0 + 1)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < Node51))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num4 + N_d0 + 2), Num5 - N_d0 - 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) > (Node51 + 1)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < (Node_5 - 3)))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Num4 + N_d0 + 2) + N_c1, Num5 - N_d0 - 2 - 2, Coef, Coef_location);
			}
			else// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Num3 / 2 + 1 + Num4 + Node_5), Num5 - Node_5 + Num6 + Num7 / 2 - 1, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	for (i = (Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Node_5 + 2); i<(Nx_2 + Nx_3 + Num1 + Num2 + Num3 + Num4 + Num5); i++)// circle num5 :  80 node!  part2 !
	{
		switch (i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4)
		{
		case (Num5 - 1) :
			boundary_node(j, j - offset, 2, -Num2 - Num3 - Num4 - Num5, Num6 + Num7 + Num8 + Num9, Coef, Coef_location);
			break;
		case (Node_5 + N_a - 4) ://2 lei shang 1j
			kind2_node(j, j - offset, -(Num1 - Node_1 + Num2 + Num3 + Num4 + Node_5), -2, -1, 0, 1, Num5 - Node_5 + Num6 + Num7 + Num8 + Node_9,
			B3, B5, B4, B1_1, B2, B3, Coef, Coef_location);
			break;
		case (Node_5 + N_a) :// 2 lei 2j shang
			kind2_node(j, j - offset, -(Num3 / 2 + Num4 + Node_5 + N_a - 3), -2, -1, 0, 1, Num5 - Node_5 - N_a + 3 + Num6 + Num7 / 2,
			B3, B5, B4, B1_2, B2, B3, Coef, Coef_location);
			break;
		case Node52:// inner 3 node
			inner_3node(j, j - offset, -(Node_55 - 4 + 1) - N_c1, -1, 0, 1, 1, k*k*dx*dx - 4, Coef, Coef_location);
			break;
		case (Node52 + 1) :// inner 3 node
			inner_3node(j, j - offset, -(Node_55 - 4 + 1), 0, 1, 1, k*k*dx*dx - 4, 1, Coef, Coef_location);
			break;
		case (Node_55 - 3) ://2 lei xia 2j
			kind2_node(j, j - offset, -(Num4 + Node_55 + 1), -1, 0, 1, 2, Num5 - Node_55 + Num6 + Num7 - 1,
			B3, B2, B1_2, B4, B5, B3, Coef, Coef_location);
			break;
		case (Node_55 + 1) ://2 lei xia 1j
			kind2_node(j, j - offset, -Num2 - Num3 - Num4 - Num5, -1, 0, 1, 2, Num6 + Num7 + Num8 + Num9,
			B3, B2, B1_1, B4, B5, B3, Coef, Coef_location);
			break;
		default:
			if ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < (Node_5 + N_a - 4))// 4dx!
			{
				inner_5node(j, j - offset, 4, -(Num1 - Node_1 + Num2 + Num3 + Num4 + Node_5), Num5 - Node_5 + Num6 + Num7 + Num8 + Node_9, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) > (Node_5 + N_a - 4)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < (Node_5 + N_a)))// 2dx!
			{
				inner_5node(j, j - offset, 2, -(Num3 / 2 + Num4 + Node_5 + N_a - 3), Num5 - Node_5 - N_a + 3 + Num6 + Num7 / 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) > (Node_5 + N_a)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < Node52))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Node_55 - 4 + 1) - N_c1, (Num5 - Node_55 + 4 + Num6 - 1) + 2, Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) > (Node52 + 1)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < (Node_55 - 3)))// 1dx
			{
				inner_5node(j, j - offset, 1, -(Node_55 - 4 + 1), (Num5 - Node_55 + 4 + Num6 - 1), Coef, Coef_location);
			}
			else if (((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) > (Node_55 - 3)) && ((i - Nx_2 - Nx_3 - Num1 - Num2 - Num3 - Num4) < (Node_55 + 1)))// 2dx
			{
				inner_5node(j, j - offset, 2, -(Num4 + Node_55 + 1), Num5 - Node_55 + Num6 + Num7 - 1, Coef, Coef_location);
			}
			else// 4dx!
			{
				inner_5node(j, j - offset, 4, -Num2 - Num3 - Num4 - Num5, Num6 + Num7 + Num8 + Num9, Coef, Coef_location);
			}
			break;
		}
		j = j + 1;
	}
	//**************************************************************
	cout << "post2圆孔区域节点数：" << j << endl;
	if (j == sizeM)
		cout << "post2 passed..." << endl;
	else
		cout << "post2 failed..." << endl;
	//修正矩阵，将不在本区域的节点刨除
	for (int i = 0; i != sizeM; i++)
		for (int j = 0; j != N_matrix; j++){
		if (Coef_location[i][j] < 0 || Coef_location[i][j] >= sizeN)
			Coef_location[i][j] = -1;
		}
	if (Debug)
	{
		ofstream outFile_Coef("post2.txt");
		if (!outFile_Coef)
		{
			cout << "Can't write the coefficient matrix into file !" << endl;
			exit(1);//Terminate with fault
		}
		else
		{
			int i, j;
			for (i = 0; i < sizeM; i++)
			{
				for (j = 0; j < N_matrix; j++)
				{
					if ((Coef[i][j].real != 0 || Coef[i][j].image != 0) && Coef_location[i][j] >= 0)
						outFile_Coef << Coef[i][j].real << "+j*" << Coef[i][j].image << " (" << i << "," << Coef_location[i][j] << ") ";
				}
				outFile_Coef << endl << endl;
			}
		}
		outFile_Coef.close();
	}
};

//共轭梯度子函数 求解A*x=b
//说明：A为系数矩阵
//      b为右端向量
//      size为矩阵阶数
//       输出为解 x
//*******************************************************************
void cgmethod(complex** A, int** Coef_Location, complex* right, complex* x, int sizeM, int sizeN, int N_matrix);
//入射激励设置
void scan_binc(complex* b, int num_total)	//4dx!!!!
{
	for (int i = 0; i != num_total; i++)
	{
		b[i].real = 0;
		b[i].image = 0;
	}
	//参数扫描
	int j = 0;
	for (int i = 0; i<Nx_side; i++)//left截断左边界 1阶mur吸收边界条件 均匀网格：4dx
	{
		switch (i)
		{
		case 0:          //角点处理：场平均吸收边条:kz/k!?!??!!?
			break;
		case (Nx_side - 1) :
			break;
		case (N_d0 - 1) :
			break;
		case (N_d0 + N_a) :
			break;
		default:
			if ((i >= N_d0) && (i < (N_d0 + N_a)))// add source
			{
				b[j].real = sin(pi*(i - N_d0 + 1)*(4 * dx) / a) * (cos(2 * kz*(4 * dx)) - 1);
				b[j].image = -sin(pi*(i - N_d0 + 1)*(4 * dx) / a) * sin(2 * kz*(4 * dx));
			}
			else
			{
			}
			break;
		}
		j = j + 1;
	}

	if (Debug)
	{
		ofstream out_b("binc.txt");
		for (int i = 0; i != num_total; i++)
		{
			if (b[i].real != 0 || b[i].image != 0)
			{
				out_b << b[i].real << "+j*" << b[i].image << "    " << i;
				out_b << endl;
			}
		}
	}
}
//更新B
void getB(complex* right, complex**A, int **Coef_Location, complex *x, int sizeM, int sizeN, int startpos);



//投影过程的封装函数
void projective_process(complex **A, int ** b, complex *B, complex *X, complex* Binc, complex* x, int startpos, int extra_nodes, double sinta, int M, int N){
	for (int i = 0; i != M; i++)
		B[i] = Binc[i + startpos];
	cgmethod(A, b, B, X, M, N, N_matrix);
	for (int i = 0; i != N; i++)
	{
		X[i] = X[i] * sinta;
		x[i + startpos + extra_nodes] = x[i + startpos + extra_nodes] + X[i];
	}
	getB(Binc, A, b, X, M, N, startpos);
}

//********************************************************************************//
//数据处理
//********************************************************************************//
void get_s_parameter(complex **Ey_total2_end){
	double s11 = 0, s21 = 0;//s11,s12
	complex x1[N_a + 2], x2, x3;//x1为左端反射波，x2为右端透射波
	x1[0].real = 0;
	x1[0].image = 0;
	x1[N_a + 1].real = 0;
	x1[N_a + 1].image = 0;
	int j;
	for (j = N_d0; j<(N_d0 + N_a); j++)
	{
		x1[j - N_d0 + 1].real = Ey_total2_end[0][j].real;
		x1[j - N_d0 + 1].image = Ey_total2_end[0][j].image;
	}
	//	complex Einc2[N_a+2];
	//	for( j=0; j<(N_a+2); j++ )
	//	{
	//		Einc2[j].real=sin( pi*j*4*dx/a );
	//		Einc2[j].image=0;
	//		x1[j].real=x1[j].real-Einc2[j].real;
	//		x1[j].image=x1[j].image-Einc2[j].image;
	//	}
	x2.real = 0;
	x2.image = 0;
	x3.real = 0;
	x3.image = 0;
	for (j = 0; j<(N_a + 1); j++)
	{
		x2.real = (x1[j].real + x1[j + 1].real) / 2;
		x2.image = (x1[j].image + x1[j + 1].image) / 2;
		x3.real = x3.real + x2.real * 4 * dx;
		x3.image = x3.image + x2.image * 4 * dx;
		//	cout<<x2.real<<"  "<<x2.image<<endl;
	}
	s11 = sqrt(x3.real*x3.real + x3.image*x3.image) / ((a / pi) * (1 - cos(pi)));
	cout << "s11 = " << s11 << " -- " << (20 * log10(s11)) << " dB" << endl;
	//求 s21
	for (j = N_d0; j<(N_d0 + N_a); j++)
	{
		x1[j - N_d0 + 1].real = Ey_total2_end[1][j + M_end - Nx_side].real;
		x1[j - N_d0 + 1].image = Ey_total2_end[1][j + M_end - Nx_side].image;
	}
	x2.real = 0;
	x2.image = 0;
	x3.real = 0;
	x3.image = 0;
	for (j = 0; j<(N_a + 1); j++)
	{
		x2.real = (x1[j].real + x1[j + 1].real) / 2;
		x2.image = (x1[j].image + x1[j + 1].image) / 2;
		x3.real = x3.real + x2.real * 4 * dx;
		x3.image = x3.image + x2.image * 4 * dx;
	}
	s21 = sqrt(x3.real*x3.real + x3.image*x3.image) / ((a / pi) * (1 - cos(pi)));
	cout << "s21 = " << s21 << " -- " << (20 * log10(s21)) << " dB" << endl;
	cout << "s11*s11+s21*s21 = " << (s11*s11 + s21*s21) << endl;
	//输出频率点及S参数
	ofstream outfile3("result1.txt", ios::app);
	outfile3 << "频率：" << (1200 / 1.483 / N_dx) << endl;
	outfile3 << "s21 = " << (20 * log10(s21)) << " dB" << endl;
	outfile3 << "s11 = " << (20 * log10(s11)) << " dB" << endl;
	outfile3.close();
}

//全部场值图的辅助函数
void show_middle(int N_p, int N_q, double** Ey_who, double** Ey_mid)
{
	int j = 0, k = 0;
	int Nx_22 = Nx_side + 2 * (N_w + 1);
	// line 1-2
	for (k = 0; k<(N_slot + 1); k++)//left
	{
		for (j = 0; j<(N_d0 - 2); j++)
			Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + k][j] = Ey_mid[N_q][(Nx_2 + Nx_3)*k + j];//1
		for (j = (N_d0 - 2); j <= (N_d0 + N_w + 2); j++)
			Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + k][j] = Ey_mid[N_q][(Nx_2 + Nx_3)*k + N_d0 - 1 + (j - N_d0 + 2) * 2];//2
		for (j = (N_d0 + N_w + 3); j<(N_d0 + N_w + N_a - 1); j++)
			Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + k][j] = Ey_mid[N_q][(Nx_2 + Nx_3)*k + Node_1 + 1 + (j - N_d0 - N_w - 3)];//3
		for (j = (N_d0 + N_w + N_a - 1); j <= (Nx_22 - N_d0 + 1); j++)
			Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + k][j] = Ey_mid[N_q][(Nx_2 + Nx_3)*k + Node_1 + N_a - 2 + (j - N_d0 - N_w - N_a + 1) * 2];//4
		for (j = (Nx_22 - N_d0 + 2); j<(Nx_22); j++)
			Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + k][j] = Ey_mid[N_q][(Nx_2 + Nx_3)*k + Node_11 + 1 + (j - Nx_22 + N_d0 - 2)];//5
	}
	// line 3--Num 5
	for (j = 0; j<(N_d0 - 2); j++)
		Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + N_slot + 1][j] = Ey_mid[N_q][(Nx_2 + Nx_3) * 1 + Num1 + Num2 + Num3 + Num4 + j];//1
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + N_slot + 1][N_d0 - 2] = Ey_mid[N_q][(Nx_2 + Nx_3) * 1 + Num1 + Num2 + Num3 + Num4 + N_d0 - 1];//2
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + N_slot + 1][N_d0 - 1] = Ey_mid[N_q][(Nx_2 + Nx_3) * 1 + Num1 + Num2 + Num3 + Num4 + N_d0 + 1];//2
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + N_slot + 1][N_d0] = Ey_mid[N_q][(Nx_2 + Nx_3) * 1 + Num1 + Num2 + Num3 + Num4 + N_d0 + 5];//2
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + N_slot + 1][N_d0 + N_w] = Ey_mid[N_q][(Nx_2 + Nx_3) * 1 + Num1 + Num2 + Num3 + Num4 + N_d0 + 8];//2
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + N_slot + 1][N_d0 + N_w + 1] = Ey_mid[N_q][(Nx_2 + Nx_3) * 1 + Num1 + Num2 + Num3 + Num4 + Node_5 - 3];//2
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + N_slot + 1][N_d0 + N_w + 2] = Ey_mid[N_q][(Nx_2 + Nx_3) * 1 + Num1 + Num2 + Num3 + Num4 + Node_5 - 1];//2
	for (j = (N_d0 + N_w + 3); j<(N_d0 + N_w + N_a - 1); j++)
		Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + N_slot + 1][j] = Ey_mid[N_q][(Nx_2 + Nx_3) * 1 + Num1 + Num2 + Num3 + Num4 + Node_5 + 1 + (j - N_d0 - N_w - 3)];//3
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + N_slot + 1][N_d0 + N_w + N_a - 1] = Ey_mid[N_q][(Nx_2 + Nx_3) * 1 + Num1 + Num2 + Num3 + Num4 + Node_5 + N_a - 2];//4
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + N_slot + 1][N_d0 + N_w + N_a] = Ey_mid[N_q][(Nx_2 + Nx_3) * 1 + Num1 + Num2 + Num3 + Num4 + Node_5 + N_a];//4
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + N_slot + 1][N_d0 + N_w + N_a + 1] = Ey_mid[N_q][(Nx_2 + Nx_3) * 1 + Num1 + Num2 + Num3 + Num4 + Node52 - 1];//4
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + N_slot + 1][Nx_22 - N_d0 - 1] = Ey_mid[N_q][(Nx_2 + Nx_3) * 1 + Num1 + Num2 + Num3 + Num4 + Node52 + 2];//4
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + N_slot + 1][Nx_22 - N_d0] = Ey_mid[N_q][(Nx_2 + Nx_3) * 1 + Num1 + Num2 + Num3 + Num4 + Node_55 - 3];//4
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + N_slot + 1][Nx_22 - N_d0 + 1] = Ey_mid[N_q][(Nx_2 + Nx_3) * 1 + Num1 + Num2 + Num3 + Num4 + Node_55 - 1];//4
	for (j = (Nx_22 - N_d0 + 2); j<(Nx_22); j++)
		Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + N_slot + 1][j] = Ey_mid[N_q][(Nx_2 + Nx_3) * 1 + Num1 + Num2 + Num3 + Num4 + Node_55 + 1 + (j - Nx_22 + N_d0 - 2)];//5

	// line 4--Num 9
	for (j = 0; j<(N_d0 - 2); j++)
		Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + N_slot + 2][j] = Ey_mid[N_q][(Nx_2 + Nx_3) * 1 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 + Num7 + Num8 + j];//1
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + N_slot + 2][N_d0 - 2] = Ey_mid[N_q][(Nx_2 + Nx_3) * 1 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 + Num7 + Num8 + N_d0 - 1];//2
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + N_slot + 2][N_d0 - 1] = Ey_mid[N_q][(Nx_2 + Nx_3) * 1 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 + Num7 + Num8 + N_d0 + 1];//2
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + N_slot + 2][N_d0 + N_w + 1] = Ey_mid[N_q][(Nx_2 + Nx_3) * 1 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 + Num7 + Num8 + Node_9 - 3];//2
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + N_slot + 2][N_d0 + N_w + 2] = Ey_mid[N_q][(Nx_2 + Nx_3) * 1 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 + Num7 + Num8 + Node_9 - 1];//2
	for (j = (N_d0 + N_w + 3); j<(N_d0 + N_w + N_a - 1); j++)
		Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + N_slot + 2][j] = Ey_mid[N_q][(Nx_2 + Nx_3) * 1 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 + Num7 + Num8 + Node_9 + 1 + (j - N_d0 - N_w - 3)];//3
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + N_slot + 2][N_d0 + N_w + N_a - 1] = Ey_mid[N_q][(Nx_2 + Nx_3) * 1 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 + Num7 + Num8 + Node_9 + N_a - 2];//4
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + N_slot + 2][N_d0 + N_w + N_a] = Ey_mid[N_q][(Nx_2 + Nx_3) * 1 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 + Num7 + Num8 + Node_9 + N_a];//4
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + N_slot + 2][Nx_22 - N_d0] = Ey_mid[N_q][(Nx_2 + Nx_3) * 1 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 + Num7 + Num8 + Node_99 - 3];//4
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + N_slot + 2][Nx_22 - N_d0 + 1] = Ey_mid[N_q][(Nx_2 + Nx_3) * 1 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 + Num7 + Num8 + Node_99 - 1];//4
	for (j = (Nx_22 - N_d0 + 2); j<(Nx_22); j++)
		Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + N_slot + 2][j] = Ey_mid[N_q][(Nx_2 + Nx_3) * 1 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 + Num7 + Num8 + Node_99 + 1 + (j - Nx_22 + N_d0 - 2)];//5
	// line 5--Num 13
	for (j = 0; j<(N_d0 - 2); j++)
		Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + N_slot + 3][j] = Ey_mid[N_q][(Nx_2 + Nx_3) * 1 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + j];//1
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + N_slot + 3][N_d0 - 2] = Ey_mid[N_q][(Nx_2 + Nx_3) * 1 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + N_d0 - 1];//2
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + N_slot + 3][N_d0 - 1] = Ey_mid[N_q][(Nx_2 + Nx_3) * 1 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + N_d0 + 1];//2
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + N_slot + 3][N_d0] = Ey_mid[N_q][(Nx_2 + Nx_3) * 1 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + N_d0 + 5];//2
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + N_slot + 3][N_d0 + N_w] = Ey_mid[N_q][(Nx_2 + Nx_3) * 1 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + N_d0 + 8];//2
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + N_slot + 3][N_d0 + N_w + 1] = Ey_mid[N_q][(Nx_2 + Nx_3) * 1 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Node_5 - 3];//2
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + N_slot + 3][N_d0 + N_w + 2] = Ey_mid[N_q][(Nx_2 + Nx_3) * 1 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Node_5 - 1];//2
	for (j = (N_d0 + N_w + 3); j<(N_d0 + N_w + N_a - 1); j++)
		Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + N_slot + 3][j] = Ey_mid[N_q][(Nx_2 + Nx_3) * 1 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Node_5 + 1 + (j - N_d0 - N_w - 3)];//3
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + N_slot + 3][N_d0 + N_w + N_a - 1] = Ey_mid[N_q][(Nx_2 + Nx_3) * 1 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Node_5 + N_a - 2];//4
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + N_slot + 3][N_d0 + N_w + N_a] = Ey_mid[N_q][(Nx_2 + Nx_3) * 1 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Node_5 + N_a];//4
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + N_slot + 3][N_d0 + N_w + N_a + 1] = Ey_mid[N_q][(Nx_2 + Nx_3) * 1 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Node52 - 1];//4
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + N_slot + 3][Nx_22 - N_d0 - 1] = Ey_mid[N_q][(Nx_2 + Nx_3) * 1 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Node52 + 2];//4
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + N_slot + 3][Nx_22 - N_d0] = Ey_mid[N_q][(Nx_2 + Nx_3) * 1 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Node_55 - 3];//4
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + N_slot + 3][Nx_22 - N_d0 + 1] = Ey_mid[N_q][(Nx_2 + Nx_3) * 1 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Node_55 - 1];//4
	for (j = (Nx_22 - N_d0 + 2); j<(Nx_22); j++)
		Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + N_slot + 3][j] = Ey_mid[N_q][(Nx_2 + Nx_3) * 1 + Num1 + Num2 + Num3 + Num4 + Num5 + Num6 * 2 + Num7 * 2 + Num8 * 2 + Num9 + Node_55 + 1 + (j - Nx_22 + N_d0 - 2)];//5

}
void show_post1(int N_p, int N_q, double** Ey_who, double** Ey_mid)
{
	int j = 0, k = 0;
	int Nx_22 = Nx_side + 2 * (N_w + 1);
	//line 1、2
	for (k = 0; k<(N_slot + 1); k++)//left
	{
		for (j = 0; j<(N_d0 - 2); j++)
			Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + k][j] = Ey_mid[N_q][(Nx_2 + Nx_3)*k + j];//1
		for (j = (N_d0 - 2); j <= (N_d0 + N_w + 2); j++)
			Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + k][j] = Ey_mid[N_q][(Nx_2 + Nx_3)*k + N_d0 - 1 + (j - N_d0 + 2) * 2];//2
		for (j = (N_d0 + N_w + 3); j<(N_d0 + N_w + N_a - 1); j++)
			Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + k][j] = Ey_mid[N_q][(Nx_2 + Nx_3)*k + Node_1_p2 + 1 + (j - N_d0 - N_w - 3)];//3
		for (j = (N_d0 + N_w + N_a - 1); j <= (Nx_22 - N_d0 + 1); j++)
			Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + k][j] = Ey_mid[N_q][(Nx_2 + Nx_3)*k + Node_1_p2 + N_a - 2 + (j - N_d0 - N_w - N_a + 1) * 2];//4
		for (j = (Nx_22 - N_d0 + 2); j<(Nx_22); j++)
			Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + k][j] = Ey_mid[N_q][(Nx_2 + Nx_3)*k + Node_11_p2 + 1 + (j - Nx_22 + N_d0 - 2)];//5
	}
	// line 3--Num 5_p2
	for (j = 0; j<(N_d0 - 2); j++)
		Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 2][j] = Ey_mid[N_q][(Nx_2 + Nx_3) * 1 + Num1_p1 + Num2_p1 + Num3_p1 + Num4_p1 + j];//1
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 2][N_d0 - 2] = Ey_mid[N_q][(Nx_2 + Nx_3) * 1 + Num1_p1 + Num2_p1 + Num3_p1 + Num4_p1 + N_d0 - 1];//2
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 2][N_d0 - 1] = Ey_mid[N_q][(Nx_2 + Nx_3) * 1 + Num1_p1 + Num2_p1 + Num3_p1 + Num4_p1 + N_d0 + 1];//2
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 2][N_d0 + 0] = Ey_mid[N_q][(Nx_2 + Nx_3) * 1 + Num1_p1 + Num2_p1 + Num3_p1 + Num4_p1 + Node51_p1 - 1];//2
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 2][N_d0 + N_w] = Ey_mid[N_q][(Nx_2 + Nx_3) * 1 + Num1_p1 + Num2_p1 + Num3_p1 + Num4_p1 + Node51_p1 + 2];//2
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 2][N_d0 + N_w + 1] = Ey_mid[N_q][(Nx_2 + Nx_3) * 1 + Num1_p1 + Num2_p1 + Num3_p1 + Num4_p1 + Node_5_p1 - 3];//2
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 2][N_d0 + N_w + 2] = Ey_mid[N_q][(Nx_2 + Nx_3) * 1 + Num1_p1 + Num2_p1 + Num3_p1 + Num4_p1 + Node_5_p1 - 1];//2
	for (j = (N_d0 + N_w + 3); j<(N_d0 + N_w + N_a - 1); j++)
		Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 2][j] = Ey_mid[N_q][(Nx_2 + Nx_3) * 1 + Num1_p1 + Num2_p1 + Num3_p1 + Num4_p1 + Node_5_p1 + 1 + (j - N_d0 - N_w - 3)];//3
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 2][N_d0 + N_w + N_a - 1] = Ey_mid[N_q][(Nx_2 + Nx_3) * 1 + Num1_p1 + Num2_p1 + Num3_p1 + Num4_p1 + Node_5_p1 + N_a - 2];//4
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 2][N_d0 + N_w + N_a] = Ey_mid[N_q][(Nx_2 + Nx_3) * 1 + Num1_p1 + Num2_p1 + Num3_p1 + Num4_p1 + Node_5_p1 + N_a];//4
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 2][N_d0 + N_w + N_a + 1] = Ey_mid[N_q][(Nx_2 + Nx_3) * 1 + Num1_p1 + Num2_p1 + Num3_p1 + Num4_p1 + Node52_p1 - 1];//4
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 2][Nx_22 - N_d0 - 1] = Ey_mid[N_q][(Nx_2 + Nx_3) * 1 + Num1_p1 + Num2_p1 + Num3_p1 + Num4_p1 + Node52_p1 + 2];//4
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 2][Nx_22 - N_d0] = Ey_mid[N_q][(Nx_2 + Nx_3) * 1 + Num1_p1 + Num2_p1 + Num3_p1 + Num4_p1 + Node_55_p1 - 3];//4
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 2][Nx_22 - N_d0 + 1] = Ey_mid[N_q][(Nx_2 + Nx_3) * 1 + Num1_p1 + Num2_p1 + Num3_p1 + Num4_p1 + Node_55_p1 - 1];//4
	for (j = (Nx_22 - N_d0 + 2); j<(Nx_22); j++)
		Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 2][j] = Ey_mid[N_q][(Nx_2 + Nx_3) * 1 + Num1_p1 + Num2_p1 + Num3_p1 + Num4_p1 + Node_55_p1 + 1 + (j - Nx_22 + N_d0 - 2)];//5

	// line 4--Num 9_p2
	for (j = 0; j<(N_d0 - 2); j++)
		Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 3][j] = Ey_mid[N_q][M_post1_q1 - Num9_p1 + j];//1
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 3][N_d0 - 2] = Ey_mid[N_q][M_post1_q1 - Num9_p1 + N_d0 - 1];//2
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 3][N_d0 - 1] = Ey_mid[N_q][M_post1_q1 - Num9_p1 + N_d0 + 1];//2
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 3][N_d0 + N_w + 1] = Ey_mid[N_q][M_post1_q1 - Num9_p1 + Node_9_p1 - 3];//2
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 3][N_d0 + N_w + 2] = Ey_mid[N_q][M_post1_q1 - Num9_p1 + Node_9_p1 - 1];//2
	for (j = (N_d0 + N_w + 3); j<(N_d0 + N_w + N_a - 1); j++)
		Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 3][j] = Ey_mid[N_q][M_post1_q1 - Num9_p1 + Node_9_p1 + 1 + (j - N_d0 - N_w - 3)];//3
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 3][N_d0 + N_w + N_a - 1] = Ey_mid[N_q][M_post1_q1 - Num9_p1 + Node_9_p1 + N_a - 2];//4
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 3][N_d0 + N_w + N_a] = Ey_mid[N_q][M_post1_q1 - Num9_p1 + Node_9_p1 + N_a];//4
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 3][Nx_22 - N_d0] = Ey_mid[N_q][M_post1_q1 - Num9_p1 + Node_99_p1 - 3];//4
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 3][Nx_22 - N_d0 + 1] = Ey_mid[N_q][M_post1_q1 - Num9_p1 + Node_99_p1 - 1];//4
	for (j = (Nx_22 - N_d0 + 2); j<(Nx_22); j++)
		Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 3][j] = Ey_mid[N_q][M_post1_q1 - Num9_p1 + Node_99_p1 + 1 + (j - Nx_22 + N_d0 - 2)];//5

	// line 5--Num 13_p2
	for (j = 0; j<(N_d0 - 2); j++)
		Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 4][j] = Ey_mid[N_q][M_post1_q1 + Num10_p1 + Num11_p1 + Num12_p1 + j];//1
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 4][N_d0 - 2] = Ey_mid[N_q][M_post1_q1 + Num10_p1 + Num11_p1 + Num12_p1 + N_d0 - 1];//2
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 4][N_d0 - 1] = Ey_mid[N_q][M_post1_q1 + Num10_p1 + Num11_p1 + Num12_p1 + N_d0 + 1];//2
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 4][N_d0 + 0] = Ey_mid[N_q][M_post1_q1 + Num10_p1 + Num11_p1 + Num12_p1 + Node131_p1 - 1];//2
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 4][N_d0 + N_w] = Ey_mid[N_q][M_post1_q1 + Num10_p1 + Num11_p1 + Num12_p1 + Node131_p1 + 2];//2
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 4][N_d0 + N_w + 1] = Ey_mid[N_q][M_post1_q1 + Num10_p1 + Num11_p1 + Num12_p1 + Node_13_p2 - 3];//2
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 4][N_d0 + N_w + 2] = Ey_mid[N_q][M_post1_q1 + Num10_p1 + Num11_p1 + Num12_p1 + Node_13_p2 - 1];//2

	for (j = (N_d0 + N_w + 3); j<(N_d0 + N_w + 7); j++)
		Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 4][j] = Ey_mid[N_q][M_post1_q1 + Num10_p1 + Num11_p1 + Num12_p1 + Node_13_p1 + 1 + (j - N_d0 - N_w - 3)];//3
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 4][N_d0 + N_w + 7] = Ey_mid[N_q][M_post1_q1 + Num10_p1 + Num11_p1 + Num12_p1 + Node132_p1 + 2];//3
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 4][N_d0 + N_w + 8] = Ey_mid[N_q][M_post1_q1 + Num10_p1 + Num11_p1 + Num12_p1 + Node132_p1 + 4];//3
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 4][N_d0 + N_w + 9] = Ey_mid[N_q][M_post1_q1 + Num10_p1 + Num11_p1 + Num12_p1 + Node132_p1 + 6];//3
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 4][N_d0 + N_w + 10] = Ey_mid[N_q][M_post1_q1 + Num10_p1 + Num11_p1 + Num12_p1 + Node132_p1 + 8];//3
	for (j = (N_d0 + N_w + 11); j<(N_d0 + N_w + N_a - 1); j++)
		Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 4][j] = Ey_mid[N_q][M_post1_q1 + Num10_p1 + Num11_p1 + Num12_p1 + Node133_p1 + (j - N_d0 - N_w - 11)];//3

	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 4][N_d0 + N_w + N_a - 1] = Ey_mid[N_q][M_post1_q1 + Num10_p1 + Num11_p1 + Num12_p1 + Node133_p1 + 5];//4
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 4][N_d0 + N_w + N_a] = Ey_mid[N_q][M_post1_q1 + Num10_p1 + Num11_p1 + Num12_p1 + Node133_p1 + 7];//4
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 4][N_d0 + N_w + N_a + 1] = Ey_mid[N_q][M_post1_q1 + Num10_p1 + Num11_p1 + Num12_p1 + Node134_p1 - 1];//4
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 4][Nx_22 - N_d0 - 1] = Ey_mid[N_q][M_post1_q1 + Num10_p1 + Num11_p1 + Num12_p1 + Node134_p1 + 2];//4
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 4][Nx_22 - N_d0] = Ey_mid[N_q][M_post1_q1 + Num10_p1 + Num11_p1 + Num12_p1 + Node_1313_p1 - 3];//4
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 4][Nx_22 - N_d0 + 1] = Ey_mid[N_q][M_post1_q1 + Num10_p1 + Num11_p1 + Num12_p1 + Node_1313_p1 - 1];//4
	for (j = (Nx_22 - N_d0 + 2); j<(Nx_22); j++)
		Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 4][j] = Ey_mid[N_q][M_post1_q1 + Num10_p1 + Num11_p1 + Num12_p1 + Node_1313_p1 + 1 + (j - Nx_22 + N_d0 - 2)];//5

	// line 6--Num 17_p2   
	for (j = 0; j<(N_d0 - 2); j++)
		Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 5][j] = Ey_mid[N_q][M_post1_q2 - Num17_p1 - Num18_p1 - Num19_p1 + j];//1
	for (j = (N_d0 - 2); j <= (N_d0 + N_w + 2); j++)
		Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 5][j] = Ey_mid[N_q][M_post1_q2 - Num17_p1 - Num18_p1 - Num19_p1 + N_d0 - 1 + (j - N_d0 + 2) * 2];//2

	for (j = (N_d0 + N_w + 3); j<(N_d0 + N_w + 7); j++)
		Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 5][j] = Ey_mid[N_q][M_post1_q2 - Num17_p1 - Num18_p1 - Num19_p1 + Node_17_p1 + 1 + (j - N_d0 - N_w - 3)];//3
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 5][N_d0 + N_w + 7] = Ey_mid[N_q][M_post1_q2 - Num17_p1 - Num18_p1 - Num19_p1 + Node173_p1 + 2];//3
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 5][N_d0 + N_w + 8] = Ey_mid[N_q][M_post1_q2 - Num17_p1 - Num18_p1 - Num19_p1 + Node173_p1 + 5];//3
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 5][N_d0 + N_w + 9] = Ey_mid[N_q][M_post1_q2 - Num17_p1 - Num18_p1 - Num19_p1 + Node174_p1 - 5];//3
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 5][N_d0 + N_w + 10] = Ey_mid[N_q][M_post1_q2 - Num17_p1 - Num18_p1 - Num19_p1 + Node174_p1 - 2];//3
	for (j = (N_d0 + N_w + 11); j<(N_d0 + N_w + N_a - 1); j++)
		Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 5][j] = Ey_mid[N_q][M_post1_q2 - Num17_p1 - Num18_p1 - Num19_p1 + Node174_p1 + (j - N_d0 - N_w - 11)];//3

	for (j = (N_d0 + N_w + N_a - 1); j <= (Nx_22 - N_d0 + 1); j++)
		Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 5][j] = Ey_mid[N_q][M_post1_q2 - Num17_p1 - Num18_p1 - Num19_p1 + Node174_p1 + 5 + (j - N_d0 - N_w - N_a + 1) * 2];//4
	for (j = (Nx_22 - N_d0 + 2); j<(Nx_22); j++)
		Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 5][j] = Ey_mid[N_q][M_post1_q2 - Num17_p1 - Num18_p1 - Num19_p1 + Node_1717_p1 + 1 + (j - Nx_22 + N_d0 - 2)];//5

	// line 7--Num 21_p2  
	for (j = 0; j<(N_d0 - 2); j++)
		Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 6][j] = Ey_mid[N_q][M_post1_q2 + Num18_p1 + j];//1
	for (j = (N_d0 - 2); j <= (N_d0 + N_w + 2); j++)
		Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 6][j] = Ey_mid[N_q][M_post1_q2 + Num18_p1 + N_d0 - 1 + (j - N_d0 + 2) * 2];//2

	for (j = (N_d0 + N_w + 3); j<(N_d0 + N_w + 7); j++)
		Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 6][j] = Ey_mid[N_q][M_post1_q2 + Num18_p1 + Node_17_p1 + 1 + (j - N_d0 - N_w - 3)];//3
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 6][N_d0 + N_w + 7] = Ey_mid[N_q][M_post1_q2 + Num18_p1 + Node173_p1 + 2];//3
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 6][N_d0 + N_w + 8] = Ey_mid[N_q][M_post1_q2 + Num18_p1 + Node173_p1 + 5];//3
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 6][N_d0 + N_w + 9] = Ey_mid[N_q][M_post1_q2 + Num18_p1 + Node174_p1 - 5];//3
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 6][N_d0 + N_w + 10] = Ey_mid[N_q][M_post1_q2 + Num18_p1 + Node174_p1 - 2];//3
	for (j = (N_d0 + N_w + 11); j<(N_d0 + N_w + N_a - 1); j++)
		Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 6][j] = Ey_mid[N_q][M_post1_q2 + Num18_p1 + Node174_p1 + (j - N_d0 - N_w - 11)];//3

	for (j = (N_d0 + N_w + N_a - 1); j <= (Nx_22 - N_d0 + 1); j++)
		Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 6][j] = Ey_mid[N_q][M_post1_q2 + Num18_p1 + Node174_p1 + 5 + (j - N_d0 - N_w - N_a + 1) * 2];//4
	for (j = (Nx_22 - N_d0 + 2); j<(Nx_22); j++)
		Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 6][j] = Ey_mid[N_q][M_post1_q2 + Num18_p1 + Node_1717_p1 + 1 + (j - Nx_22 + N_d0 - 2)];//5

	// line 8--Num 25_p2 (Num 13_p2)
	for (j = 0; j<(N_d0 - 2); j++)
		Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 7][j] = Ey_mid[N_q][M_post1_q2 + Num18_p1 + Num17_p1 + Num16_p1 + Num15_p1 + Num14_p1 + j];//1
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 7][N_d0 - 2] = Ey_mid[N_q][M_post1_q2 + Num18_p1 + Num17_p1 + Num16_p1 + Num15_p1 + Num14_p1 + N_d0 - 1];//2
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 7][N_d0 - 1] = Ey_mid[N_q][M_post1_q2 + Num18_p1 + Num17_p1 + Num16_p1 + Num15_p1 + Num14_p1 + N_d0 + 1];//2
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 7][N_d0 + 0] = Ey_mid[N_q][M_post1_q2 + Num18_p1 + Num17_p1 + Num16_p1 + Num15_p1 + Num14_p1 + Node131_p1 - 1];//2
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 7][N_d0 + N_w] = Ey_mid[N_q][M_post1_q2 + Num18_p1 + Num17_p1 + Num16_p1 + Num15_p1 + Num14_p1 + Node131_p1 + 2];//2
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 7][N_d0 + N_w + 1] = Ey_mid[N_q][M_post1_q2 + Num18_p1 + Num17_p1 + Num16_p1 + Num15_p1 + Num14_p1 + Node_13_p2 - 3];//2
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 7][N_d0 + N_w + 2] = Ey_mid[N_q][M_post1_q2 + Num18_p1 + Num17_p1 + Num16_p1 + Num15_p1 + Num14_p1 + Node_13_p2 - 1];//2

	for (j = (N_d0 + N_w + 3); j<(N_d0 + N_w + 7); j++)
		Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 7][j] = Ey_mid[N_q][M_post1_q2 + Num18_p1 + Num17_p1 + Num16_p1 + Num15_p1 + Num14_p1 + Node_13_p1 + 1 + (j - N_d0 - N_w - 3)];//3
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 7][N_d0 + N_w + 7] = Ey_mid[N_q][M_post1_q2 + Num18_p1 + Num17_p1 + Num16_p1 + Num15_p1 + Num14_p1 + Node132_p1 + 2];//3
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 7][N_d0 + N_w + 8] = Ey_mid[N_q][M_post1_q2 + Num18_p1 + Num17_p1 + Num16_p1 + Num15_p1 + Num14_p1 + Node132_p1 + 4];//3
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 7][N_d0 + N_w + 9] = Ey_mid[N_q][M_post1_q2 + Num18_p1 + Num17_p1 + Num16_p1 + Num15_p1 + Num14_p1 + Node132_p1 + 6];//3
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 7][N_d0 + N_w + 10] = Ey_mid[N_q][M_post1_q2 + Num18_p1 + Num17_p1 + Num16_p1 + Num15_p1 + Num14_p1 + Node132_p1 + 8];//3
	for (j = (N_d0 + N_w + 11); j<(N_d0 + N_w + N_a - 1); j++)
		Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 7][j] = Ey_mid[N_q][M_post1_q2 + Num18_p1 + Num17_p1 + Num16_p1 + Num15_p1 + Num14_p1 + Node133_p1 + (j - N_d0 - N_w - 11)];//3

	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 7][N_d0 + N_w + N_a - 1] = Ey_mid[N_q][M_post1_q2 + Num18_p1 + Num17_p1 + Num16_p1 + Num15_p1 + Num14_p1 + Node133_p1 + 5];//4
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 7][N_d0 + N_w + N_a] = Ey_mid[N_q][M_post1_q2 + Num18_p1 + Num17_p1 + Num16_p1 + Num15_p1 + Num14_p1 + Node133_p1 + 7];//4
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 7][N_d0 + N_w + N_a + 1] = Ey_mid[N_q][M_post1_q2 + Num18_p1 + Num17_p1 + Num16_p1 + Num15_p1 + Num14_p1 + Node134_p1 - 1];//4
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 7][Nx_22 - N_d0 - 1] = Ey_mid[N_q][M_post1_q2 + Num18_p1 + Num17_p1 + Num16_p1 + Num15_p1 + Num14_p1 + Node134_p1 + 2];//4
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 7][Nx_22 - N_d0] = Ey_mid[N_q][M_post1_q2 + Num18_p1 + Num17_p1 + Num16_p1 + Num15_p1 + Num14_p1 + Node_1313_p1 - 3];//4
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 7][Nx_22 - N_d0 + 1] = Ey_mid[N_q][M_post1_q2 + Num18_p1 + Num17_p1 + Num16_p1 + Num15_p1 + Num14_p1 + Node_1313_p1 - 1];//4
	for (j = (Nx_22 - N_d0 + 2); j<(Nx_22); j++)
		Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 7][j] = Ey_mid[N_q][M_post1_q2 + Num18_p1 + Num17_p1 + Num16_p1 + Num15_p1 + Num14_p1 + Node_1313_p1 + 1 + (j - Nx_22 + N_d0 - 2)];//5

	// line 9--Num 29_p2 (4--Num 9_p2)
	for (j = 0; j<(N_d0 - 2); j++)
		Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 8][j] = Ey_mid[N_q][M_post1_q3 - Num9_p1 + j];//1
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 8][N_d0 - 2] = Ey_mid[N_q][M_post1_q3 - Num9_p1 + N_d0 - 1];//2
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 8][N_d0 - 1] = Ey_mid[N_q][M_post1_q3 - Num9_p1 + N_d0 + 1];//2
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 8][N_d0 + N_w + 1] = Ey_mid[N_q][M_post1_q3 - Num9_p1 + Node_9_p1 - 3];//2
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 8][N_d0 + N_w + 2] = Ey_mid[N_q][M_post1_q3 - Num9_p1 + Node_9_p1 - 1];//2
	for (j = (N_d0 + N_w + 3); j<(N_d0 + N_w + N_a - 1); j++)
		Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 8][j] = Ey_mid[N_q][M_post1_q3 - Num9_p1 + Node_9_p1 + 1 + (j - N_d0 - N_w - 3)];//3
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 8][N_d0 + N_w + N_a - 1] = Ey_mid[N_q][M_post1_q3 - Num9_p1 + Node_9_p1 + N_a - 2];//4
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 8][N_d0 + N_w + N_a] = Ey_mid[N_q][M_post1_q3 - Num9_p1 + Node_9_p1 + N_a];//4
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 8][Nx_22 - N_d0] = Ey_mid[N_q][M_post1_q3 - Num9_p1 + Node_99_p1 - 3];//4
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 8][Nx_22 - N_d0 + 1] = Ey_mid[N_q][M_post1_q3 - Num9_p1 + Node_99_p1 - 1];//4
	for (j = (Nx_22 - N_d0 + 2); j<(Nx_22); j++)
		Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 8][j] = Ey_mid[N_q][M_post1_q3 - Num9_p1 + Node_99_p1 + 1 + (j - Nx_22 + N_d0 - 2)];//5

	// line 10--Num 33_p2 (3--Num 5_p2)
	for (j = 0; j<(N_d0 - 2); j++)
		Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 9][j] = Ey_mid[N_q][M_post1_q3 + Num8_p1 + Num7_p1 + Num6_p1 + j];//1
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 9][N_d0 - 2] = Ey_mid[N_q][M_post1_q3 + Num8_p1 + Num7_p1 + Num6_p1 + N_d0 - 1];//2
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 9][N_d0 - 1] = Ey_mid[N_q][M_post1_q3 + Num8_p1 + Num7_p1 + Num6_p1 + N_d0 + 1];//2
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 9][N_d0 + 0] = Ey_mid[N_q][M_post1_q3 + Num8_p1 + Num7_p1 + Num6_p1 + Node51_p1 - 1];//2
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 9][N_d0 + N_w] = Ey_mid[N_q][M_post1_q3 + Num8_p1 + Num7_p1 + Num6_p1 + Node51_p1 + 2];//2
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 9][N_d0 + N_w + 1] = Ey_mid[N_q][M_post1_q3 + Num8_p1 + Num7_p1 + Num6_p1 + Node_5_p1 - 3];//2
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 9][N_d0 + N_w + 2] = Ey_mid[N_q][M_post1_q3 + Num8_p1 + Num7_p1 + Num6_p1 + Node_5_p1 - 1];//2
	for (j = (N_d0 + N_w + 3); j<(N_d0 + N_w + N_a - 1); j++)
		Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 9][j] = Ey_mid[N_q][M_post1_q3 + Num8_p1 + Num7_p1 + Num6_p1 + Node_5_p1 + 1 + (j - N_d0 - N_w - 3)];//3
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 9][N_d0 + N_w + N_a - 1] = Ey_mid[N_q][M_post1_q3 + Num8_p1 + Num7_p1 + Num6_p1 + Node_5_p1 + N_a - 2];//4
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 9][N_d0 + N_w + N_a] = Ey_mid[N_q][M_post1_q3 + Num8_p1 + Num7_p1 + Num6_p1 + Node_5_p1 + N_a];//4
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 9][N_d0 + N_w + N_a + 1] = Ey_mid[N_q][M_post1_q3 + Num8_p1 + Num7_p1 + Num6_p1 + Node52_p1 - 1];//4
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 9][Nx_22 - N_d0 - 1] = Ey_mid[N_q][M_post1_q3 + Num8_p1 + Num7_p1 + Num6_p1 + Node52_p1 + 2];//4
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 9][Nx_22 - N_d0] = Ey_mid[N_q][M_post1_q3 + Num8_p1 + Num7_p1 + Num6_p1 + Node_55_p1 - 3];//4
	Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 9][Nx_22 - N_d0 + 1] = Ey_mid[N_q][M_post1_q3 + Num8_p1 + Num7_p1 + Num6_p1 + Node_55_p1 - 1];//4
	for (j = (Nx_22 - N_d0 + 2); j<(Nx_22); j++)
		Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 9][j] = Ey_mid[N_q][M_post1_q3 + Num8_p1 + Num7_p1 + Num6_p1 + Node_55_p1 + 1 + (j - Nx_22 + N_d0 - 2)];//5

	for (k = 0; k<(N_slot + 1); k++)//left
	{
		for (j = 0; j<(N_d0 - 2); j++)
			Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 10 + k][j] = Ey_mid[N_q][M_post1 - Num1_p1 - (Nx_2 + Nx_3) + (Nx_2 + Nx_3)*k + j];//1
		for (j = (N_d0 - 2); j <= (N_d0 + N_w + 2); j++)
			Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 10 + k][j] = Ey_mid[N_q][M_post1 - Num1_p1 - (Nx_2 + Nx_3) + (Nx_2 + Nx_3)*k + N_d0 - 1 + (j - N_d0 + 2) * 2];//2
		for (j = (N_d0 + N_w + 3); j<(N_d0 + N_w + N_a - 1); j++)
			Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 10 + k][j] = Ey_mid[N_q][M_post1 - Num1_p1 - (Nx_2 + Nx_3) + (Nx_2 + Nx_3)*k + Node_1_p2 + 1 + (j - N_d0 - N_w - 3)];//3
		for (j = (N_d0 + N_w + N_a - 1); j <= (Nx_22 - N_d0 + 1); j++)
			Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 10 + k][j] = Ey_mid[N_q][M_post1 - Num1_p1 - (Nx_2 + Nx_3) + (Nx_2 + Nx_3)*k + Node_1_p2 + N_a - 2 + (j - N_d0 - N_w - N_a + 1) * 2];//4
		for (j = (Nx_22 - N_d0 + 2); j<(Nx_22); j++)
			Ey_who[N_d1 + 1 + (N_s + N_w)*N_p + 10 + k][j] = Ey_mid[N_q][M_post1 - Num1_p1 - (Nx_2 + Nx_3) + (Nx_2 + Nx_3)*k + Node_11_p2 + 1 + (j - Nx_22 + N_d0 - 2)];//5
	}


}
//显示完整场图
void show_value(int M){
	vector<double> x(M);

	ifstream in_e("Evalue.txt");
	for (int i = 0; i != M; i++)
		in_e >> x[i];
	in_e.close();
	const int overlapping_nodes = 2 * Num2 + Num3;
	int i, j;
	const int Nx_22 = 58;//Nx_side+2*(N_w+1);
	int k = 0;
	double **Ey_whole = new double*[(N_d1 + 1) * 2 + (N_circle + 2)*(N_w + N_s) + N_s - 1];//括号内代表行数
	for (j = 0; j<((N_d1 + 1) * 2 + (N_circle + 2)*(N_w + N_s) + N_s - 1); j++)
	{
		Ey_whole[j] = new double[Nx_22];
	}
	//将总场值按照分区存储
	double **Ey_total2_end = new double*[N_end];//
	for (j = 0; j<N_end; j++)
	{
		Ey_total2_end[j] = new double[M_end];
	}
	double **Ey_total2_middle = new double*[N_middle];//
	for (j = 0; j<N_middle; j++)
	{
		Ey_total2_middle[j] = new double[M_middle];
	}
	double **Ey_total2_post1 = new double*[N_post1];//
	for (j = 0; j<N_post1; j++)
	{
		Ey_total2_post1[j] = new double[M_post1];
	}
	double **Ey_total2_post2 = new double*[N_post2];//
	for (j = 0; j<N_post2; j++)
	{
		Ey_total2_post2[j] = new double[M_post2];
	}
	int startpos[13];
	startpos[0] = 0;												//end0
	startpos[1] = startpos[0] + M_end - overlapping_nodes;
	startpos[2] = startpos[1] + M_middle - overlapping_nodes;
	startpos[3] = startpos[2] + M_middle - overlapping_nodes;		//post1
	startpos[4] = startpos[3] + M_post1 - overlapping_nodes;
	startpos[5] = startpos[4] + M_middle - overlapping_nodes;
	startpos[6] = startpos[5] + M_middle - overlapping_nodes;		//post2
	startpos[7] = startpos[6] + M_post2 - overlapping_nodes;
	startpos[8] = startpos[7] + M_middle - overlapping_nodes;
	startpos[9] = startpos[8] + M_middle - overlapping_nodes;		//post1
	startpos[10] = startpos[9] + M_post1 - overlapping_nodes;
	startpos[11] = startpos[10] + M_middle - overlapping_nodes;
	startpos[12] = startpos[11] + M_middle - overlapping_nodes;		//end1
	cout << "startpos of zone:" << endl;
	for (int i = 0; i != 13; i++)
		cout << setw(2) << i << "/" << setw(5) << startpos[i] << endl;
	cout << M - M_end << endl;

	//将x转化为Ey_total2_*
	for (int i = 0; i != M_end; i++)									//end0
		Ey_total2_end[0][i] = x[i + startpos[0]];
	for (int i = 0; i != M_middle; i++)
		Ey_total2_middle[0][i] = x[i + startpos[1]];
	for (int i = 0; i != M_middle; i++)
		Ey_total2_middle[1][i] = x[i + startpos[2]];
	for (int i = 0; i != M_post1; i++)
		Ey_total2_post1[0][i] = x[i + startpos[3]];						//post1
	for (int i = 0; i != M_middle; i++)
		Ey_total2_middle[2][i] = x[i + startpos[4]];
	for (int i = 0; i != M_middle; i++)
		Ey_total2_middle[3][i] = x[i + startpos[5]];
	for (int i = 0; i != M_middle; i++)
		Ey_total2_post2[0][i] = x[i + startpos[6]];						//post2
	for (int i = 0; i != M_middle; i++)
		Ey_total2_middle[4][i] = x[i + startpos[7]];
	for (int i = 0; i != M_middle; i++)
		Ey_total2_middle[5][i] = x[i + startpos[8]];
	for (int i = 0; i != M_post1; i++)
		Ey_total2_post1[0][i] = x[i + startpos[9]];						//post1
	for (int i = 0; i != M_end; i++)
		for (int i = 0; i != M_middle; i++)
			Ey_total2_middle[4][i] = x[i + startpos[10]];
	for (int i = 0; i != M_middle; i++)
		Ey_total2_end[1][i] = x[i + startpos[11]];						//end1

	//初始化，最大值
	for (i = 0; i<((N_d1 + 1) * 2 + (N_circle + 2)*(N_w + N_s) + N_s - 1); i++)
	{
		for (j = 0; j<Nx_22; j++)
		{
			Ey_whole[i][j] = 0.5;//!!
		}
	}

	//****************************************************
	//********************region_side_in******************
	//****************************************************
	for (i = 0; i<(N_d1 + 1); i++)
	{
		for (j = 0; j<N_d0; j++)
			Ey_whole[i][j] = Ey_total2_end[0][j + Nx_side*i];
		for (j = N_d0; j<(N_d0 + N_a); j++)
			Ey_whole[i][j + N_w + 1] = Ey_total2_end[0][j + Nx_side*i];
		for (j = (N_d0 + N_a); j<Nx_side; j++)
			Ey_whole[i][j + (N_w + 1) * 2] = Ey_total2_end[0][j + Nx_side*i];
	}

	//****************************************************
	//********************region_middle*******************
	//****************************************************

	// 普通圆孔0-1
	show_middle(0, 0, Ey_whole, Ey_total2_middle);
	show_middle(1, 1, Ey_whole, Ey_total2_middle);


	// 第一个post1 小圆孔 !!
	show_post1(2, 0, Ey_whole, Ey_total2_post1);
	// 普通圆孔3-4
	show_middle(4, 2, Ey_whole, Ey_total2_middle);
	show_middle(5, 3, Ey_whole, Ey_total2_middle);
	// post2 大圆孔 !!!  ( 6&7 )  10 line
	//line 1、2
	for (k = 0; k<(N_slot + 1); k++)//left
	{
		for (j = 0; j<(N_d0 - 2); j++)
			Ey_whole[N_d1 + 1 + (N_s + N_w) * 6 + k][j] = Ey_total2_post2[0][(Nx_2 + Nx_3)*k + j];//1
		for (j = (N_d0 - 2); j <= (N_d0 + N_w + 2); j++)
			Ey_whole[N_d1 + 1 + (N_s + N_w) * 6 + k][j] = Ey_total2_post2[0][(Nx_2 + Nx_3)*k + N_d0 - 1 + (j - N_d0 + 2) * 2];//2
		for (j = (N_d0 + N_w + 3); j<(N_d0 + N_w + N_a - 1); j++)
			Ey_whole[N_d1 + 1 + (N_s + N_w) * 6 + k][j] = Ey_total2_post2[0][(Nx_2 + Nx_3)*k + Node_1_p2 + 1 + (j - N_d0 - N_w - 3)];//3
		for (j = (N_d0 + N_w + N_a - 1); j <= (Nx_22 - N_d0 + 1); j++)
			Ey_whole[N_d1 + 1 + (N_s + N_w) * 6 + k][j] = Ey_total2_post2[0][(Nx_2 + Nx_3)*k + Node_1_p2 + N_a - 2 + (j - N_d0 - N_w - N_a + 1) * 2];//4
		for (j = (Nx_22 - N_d0 + 2); j<(Nx_22); j++)
			Ey_whole[N_d1 + 1 + (N_s + N_w) * 6 + k][j] = Ey_total2_post2[0][(Nx_2 + Nx_3)*k + Node_11_p2 + 1 + (j - Nx_22 + N_d0 - 2)];//5
	}
	// line 3--Num 5_p2
	for (j = 0; j<(N_d0 - 2); j++)
		Ey_whole[N_d1 + 1 + (N_s + N_w) * 6 + N_slot + 1][j] = Ey_total2_post2[0][(Nx_2 + Nx_3) * 1 + Num1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + j];//1
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 6 + N_slot + 1][N_d0 - 2] = Ey_total2_post2[0][(Nx_2 + Nx_3) * 1 + Num1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + N_d0 - 1];//2
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 6 + N_slot + 1][N_d0 - 1] = Ey_total2_post2[0][(Nx_2 + Nx_3) * 1 + Num1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + N_d0 + 1];//2
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 6 + N_slot + 1][N_d0] = Ey_total2_post2[0][(Nx_2 + Nx_3) * 1 + Num1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + N_d0 + 5];//2
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 6 + N_slot + 1][N_d0 + N_w] = Ey_total2_post2[0][(Nx_2 + Nx_3) * 1 + Num1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + N_d0 + 8];//2
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 6 + N_slot + 1][N_d0 + N_w + 1] = Ey_total2_post2[0][(Nx_2 + Nx_3) * 1 + Num1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + Node_5_p2 - 3];//2
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 6 + N_slot + 1][N_d0 + N_w + 2] = Ey_total2_post2[0][(Nx_2 + Nx_3) * 1 + Num1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + Node_5_p2 - 1];//2
	for (j = (N_d0 + N_w + 3); j<(N_d0 + N_w + N_a - 1); j++)
		Ey_whole[N_d1 + 1 + (N_s + N_w) * 6 + N_slot + 1][j] = Ey_total2_post2[0][(Nx_2 + Nx_3) * 1 + Num1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + Node_5_p2 + 1 + (j - N_d0 - N_w - 3)];//3
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 6 + N_slot + 1][N_d0 + N_w + N_a - 1] = Ey_total2_post2[0][(Nx_2 + Nx_3) * 1 + Num1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + Node_5_p2 + N_a - 2];//4
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 6 + N_slot + 1][N_d0 + N_w + N_a] = Ey_total2_post2[0][(Nx_2 + Nx_3) * 1 + Num1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + Node_5_p2 + N_a];//4
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 6 + N_slot + 1][N_d0 + N_w + N_a + 1] = Ey_total2_post2[0][(Nx_2 + Nx_3) * 1 + Num1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + Node54_p2 - 1];//4
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 6 + N_slot + 1][Nx_22 - N_d0 - 1] = Ey_total2_post2[0][(Nx_2 + Nx_3) * 1 + Num1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + Node54_p2 + 2];//4
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 6 + N_slot + 1][Nx_22 - N_d0] = Ey_total2_post2[0][(Nx_2 + Nx_3) * 1 + Num1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + Node_55_p2 - 3];//4
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 6 + N_slot + 1][Nx_22 - N_d0 + 1] = Ey_total2_post2[0][(Nx_2 + Nx_3) * 1 + Num1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + Node_55_p2 - 1];//4
	for (j = (Nx_22 - N_d0 + 2); j<(Nx_22); j++)
		Ey_whole[N_d1 + 1 + (N_s + N_w) * 6 + N_slot + 1][j] = Ey_total2_post2[0][(Nx_2 + Nx_3) * 1 + Num1_p2 + Num2_p2 + Num3_p2 + Num4_p2 + Node_55_p2 + 1 + (j - Nx_22 + N_d0 - 2)];//5

	// line 4--Num 9_p2
	for (j = 0; j<(N_d0 - 2); j++)
		Ey_whole[N_d1 + 1 + (N_s + N_w) * 6 + N_slot + 2][j] = Ey_total2_post2[0][M_post2_q1 - Num9_p2 + j];//1
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 6 + N_slot + 2][N_d0 - 2] = Ey_total2_post2[0][M_post2_q1 - Num9_p2 + N_d0 - 1];//2
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 6 + N_slot + 2][N_d0 - 1] = Ey_total2_post2[0][M_post2_q1 - Num9_p2 + N_d0 + 1];//2
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 6 + N_slot + 2][N_d0 + N_w + 1] = Ey_total2_post2[0][M_post2_q1 - Num9_p2 + Node_9_p2 - 3];//2
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 6 + N_slot + 2][N_d0 + N_w + 2] = Ey_total2_post2[0][M_post2_q1 - Num9_p2 + Node_9_p2 - 1];//2

	for (j = (N_d0 + N_w + 3); j<(N_d0 + N_w + 6); j++)
		Ey_whole[N_d1 + 1 + (N_s + N_w) * 6 + N_slot + 2][j] = Ey_total2_post2[0][M_post2_q1 - Num9_p2 + Node_9_p2 + 1 + (j - N_d0 - N_w - 3)];//3
	for (j = (N_d0 + N_w + 6); j<(N_d0 + N_w + 12); j++)
		Ey_whole[N_d1 + 1 + (N_s + N_w) * 6 + N_slot + 2][j] = Ey_total2_post2[0][M_post2_q1 - Num9_p2 + Node_9_p2 + 5 + (j - N_d0 - N_w - 6) * 2];//3
	for (j = (N_d0 + N_w + 12); j<(N_d0 + N_w + N_a - 1); j++)
		Ey_whole[N_d1 + 1 + (N_s + N_w) * 6 + N_slot + 2][j] = Ey_total2_post2[0][M_post2_q1 - Num9_p2 + Node93_p2 + 1 + (j - N_d0 - N_w - 12)];//3

	Ey_whole[N_d1 + 1 + (N_s + N_w) * 6 + N_slot + 2][N_d0 + N_w + N_a - 1] = Ey_total2_post2[0][M_post2_q1 - Num9_p2 + Node93_p2 + 5];//4
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 6 + N_slot + 2][N_d0 + N_w + N_a] = Ey_total2_post2[0][M_post2_q1 - Num9_p2 + Node93_p2 + 7];//4
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 6 + N_slot + 2][Nx_22 - N_d0] = Ey_total2_post2[0][M_post2_q1 - Num9_p2 + Node_99_p2 - 3];//4
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 6 + N_slot + 2][Nx_22 - N_d0 + 1] = Ey_total2_post2[0][M_post2_q1 - Num9_p2 + Node_99_p2 - 1];//4
	for (j = (Nx_22 - N_d0 + 2); j<(Nx_22); j++)
		Ey_whole[N_d1 + 1 + (N_s + N_w) * 6 + N_slot + 2][j] = Ey_total2_post2[0][M_post2_q1 - Num9_p2 + Node_99_p2 + 1 + (j - Nx_22 + N_d0 - 2)];//5

	// line 5--Num 13_p2
	for (j = 0; j<(N_d0 - 2); j++)
		Ey_whole[N_d1 + 1 + (N_s + N_w) * 6 + N_slot + 3][j] = Ey_total2_post2[0][M_post2_q1 + Num10_p2 + Num11_p2 + Num12_p2 + j];//1
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 6 + N_slot + 3][N_d0 - 2] = Ey_total2_post2[0][M_post2_q1 + Num10_p2 + Num11_p2 + Num12_p2 + N_d0 - 1];//2
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 6 + N_slot + 3][N_d0 - 1] = Ey_total2_post2[0][M_post2_q1 + Num10_p2 + Num11_p2 + Num12_p2 + N_d0 + 1];//2
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 6 + N_slot + 3][N_d0] = Ey_total2_post2[0][M_post2_q1 + Num10_p2 + Num11_p2 + Num12_p2 + N_d0 + 5];//2
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 6 + N_slot + 3][N_d0 + N_w] = Ey_total2_post2[0][M_post2_q1 + Num10_p2 + Num11_p2 + Num12_p2 + N_d0 + 8];//2
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 6 + N_slot + 3][N_d0 + N_w + 1] = Ey_total2_post2[0][M_post2_q1 + Num10_p2 + Num11_p2 + Num12_p2 + Node_13_p2 - 3];//2
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 6 + N_slot + 3][N_d0 + N_w + 2] = Ey_total2_post2[0][M_post2_q1 + Num10_p2 + Num11_p2 + Num12_p2 + Node_13_p2 - 1];//2

	for (j = (N_d0 + N_w + 3); j<(N_d0 + N_w + 6); j++)
		Ey_whole[N_d1 + 1 + (N_s + N_w) * 6 + N_slot + 3][j] = Ey_total2_post2[0][M_post2_q1 + Num10_p2 + Num11_p2 + Num12_p2 + Node_13_p2 + 1 + (j - N_d0 - N_w - 3)];//3
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 6 + N_slot + 3][N_d0 + N_w + 6] = Ey_total2_post2[0][M_post2_q1 + Num10_p2 + Num11_p2 + Num12_p2 + Node_13_p2 + 5];//3
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 6 + N_slot + 3][N_d0 + N_w + 7] = Ey_total2_post2[0][M_post2_q1 + Num10_p2 + Num11_p2 + Num12_p2 + Node_13_p2 + 8];//3
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 6 + N_slot + 3][N_d0 + N_w + 8] = Ey_total2_post2[0][M_post2_q1 + Num10_p2 + Num11_p2 + Num12_p2 + Node_13_p2 + 12];//3
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 6 + N_slot + 3][N_d0 + N_w + 9] = Ey_total2_post2[0][M_post2_q1 + Num10_p2 + Num11_p2 + Num12_p2 + Node_13_p2 + 16];//3
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 6 + N_slot + 3][N_d0 + N_w + 10] = Ey_total2_post2[0][M_post2_q1 + Num10_p2 + Num11_p2 + Num12_p2 + Node133_p2 - 4];//3
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 6 + N_slot + 3][N_d0 + N_w + 11] = Ey_total2_post2[0][M_post2_q1 + Num10_p2 + Num11_p2 + Num12_p2 + Node133_p2 - 1];//3
	for (j = (N_d0 + N_w + 12); j<(N_d0 + N_w + N_a - 1); j++)
		Ey_whole[N_d1 + 1 + (N_s + N_w) * 6 + N_slot + 3][j] = Ey_total2_post2[0][M_post2_q1 + Num10_p2 + Num11_p2 + Num12_p2 + Node133_p2 + 1 + (j - N_d0 - N_w - 12)];//3

	Ey_whole[N_d1 + 1 + (N_s + N_w) * 6 + N_slot + 3][N_d0 + N_w + N_a - 1] = Ey_total2_post2[0][M_post2_q1 + Num10_p2 + Num11_p2 + Num12_p2 + Node133_p2 + 5];//4
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 6 + N_slot + 3][N_d0 + N_w + N_a] = Ey_total2_post2[0][M_post2_q1 + Num10_p2 + Num11_p2 + Num12_p2 + Node133_p2 + 7];//4
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 6 + N_slot + 3][N_d0 + N_w + N_a + 1] = Ey_total2_post2[0][M_post2_q1 + Num10_p2 + Num11_p2 + Num12_p2 + Node134_p2 - 1];//4
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 6 + N_slot + 3][Nx_22 - N_d0 - 1] = Ey_total2_post2[0][M_post2_q1 + Num10_p2 + Num11_p2 + Num12_p2 + Node134_p2 + 2];//4
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 6 + N_slot + 3][Nx_22 - N_d0] = Ey_total2_post2[0][M_post2_q1 + Num10_p2 + Num11_p2 + Num12_p2 + Node_1313_p2 - 3];//4
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 6 + N_slot + 3][Nx_22 - N_d0 + 1] = Ey_total2_post2[0][M_post2_q1 + Num10_p2 + Num11_p2 + Num12_p2 + Node_1313_p2 - 1];//4
	for (j = (Nx_22 - N_d0 + 2); j<(Nx_22); j++)
		Ey_whole[N_d1 + 1 + (N_s + N_w) * 6 + N_slot + 3][j] = Ey_total2_post2[0][M_post2_q1 + Num10_p2 + Num11_p2 + Num12_p2 + Node_1313_p2 + 1 + (j - Nx_22 + N_d0 - 2)];//5

	// line 6--Num 17_p2
	for (j = 0; j<(N_d0 - 2); j++)
		Ey_whole[N_d1 + 1 + (N_s + N_w) * 7][j] = Ey_total2_post2[0][M_post2_q2 - Num17_p2 - Num18_p2 - Num19_p2 + j];//1
	for (j = (N_d0 - 2); j <= (N_d0 + N_w + 2); j++)
		Ey_whole[N_d1 + 1 + (N_s + N_w) * 7][j] = Ey_total2_post2[0][M_post2_q2 - Num17_p2 - Num18_p2 - Num19_p2 + N_d0 - 1 + (j - N_d0 + 2) * 2];//2

	for (j = (N_d0 + N_w + 3); j<(N_d0 + N_w + 6); j++)
		Ey_whole[N_d1 + 1 + (N_s + N_w) * 7][j] = Ey_total2_post2[0][M_post2_q2 - Num17_p2 - Num18_p2 - Num19_p2 + Node_17_p2 + 1 + (j - N_d0 - N_w - 3)];//3
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 7][N_d0 + N_w + 6] = Ey_total2_post2[0][M_post2_q2 - Num17_p2 - Num18_p2 - Num19_p2 + Node173_p2 + 1];//3
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 7][N_d0 + N_w + 7] = Ey_total2_post2[0][M_post2_q2 - Num17_p2 - Num18_p2 - Num19_p2 + Node173_p2 + 4];//3
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 7][N_d0 + N_w + 10] = Ey_total2_post2[0][M_post2_q2 - Num17_p2 - Num18_p2 - Num19_p2 + Node174_p2 - 4];//3
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 7][N_d0 + N_w + 11] = Ey_total2_post2[0][M_post2_q2 - Num17_p2 - Num18_p2 - Num19_p2 + Node174_p2 - 1];//3
	for (j = (N_d0 + N_w + 12); j<(N_d0 + N_w + N_a - 1); j++)
		Ey_whole[N_d1 + 1 + (N_s + N_w) * 7][j] = Ey_total2_post2[0][M_post2_q2 - Num17_p2 - Num18_p2 - Num19_p2 + Node174_p2 + 1 + (j - N_d0 - N_w - 12)];//3

	for (j = (N_d0 + N_w + N_a - 1); j <= (Nx_22 - N_d0 + 1); j++)
		Ey_whole[N_d1 + 1 + (N_s + N_w) * 7][j] = Ey_total2_post2[0][M_post2_q2 - Num17_p2 - Num18_p2 - Num19_p2 + Node174_p2 + 5 + (j - N_d0 - N_w - N_a + 1) * 2];//4
	for (j = (Nx_22 - N_d0 + 2); j<(Nx_22); j++)
		Ey_whole[N_d1 + 1 + (N_s + N_w) * 7][j] = Ey_total2_post2[0][M_post2_q2 - Num17_p2 - Num18_p2 - Num19_p2 + Node_1717_p2 + 1 + (j - Nx_22 + N_d0 - 2)];//5

	// line 7--Num 21_p2
	for (j = 0; j<(N_d0 - 2); j++)
		Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 1][j] = Ey_total2_post2[0][M_post2_q2 + Num18_p2 + j];//1
	for (j = (N_d0 - 2); j <= (N_d0 + N_w + 2); j++)
		Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 1][j] = Ey_total2_post2[0][M_post2_q2 + Num18_p2 + N_d0 - 1 + (j - N_d0 + 2) * 2];//2

	for (j = (N_d0 + N_w + 3); j<(N_d0 + N_w + 6); j++)
		Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 1][j] = Ey_total2_post2[0][M_post2_q2 + Num18_p2 + Node_17_p2 + 1 + (j - N_d0 - N_w - 3)];//3
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 1][N_d0 + N_w + 6] = Ey_total2_post2[0][M_post2_q2 + Num18_p2 + Node173_p2 + 1];//3
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 1][N_d0 + N_w + 7] = Ey_total2_post2[0][M_post2_q2 + Num18_p2 + Node173_p2 + 4];//3
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 1][N_d0 + N_w + 10] = Ey_total2_post2[0][M_post2_q2 + Num18_p2 + Node174_p2 - 4];//3
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 1][N_d0 + N_w + 11] = Ey_total2_post2[0][M_post2_q2 + Num18_p2 + Node174_p2 - 1];//3
	for (j = (N_d0 + N_w + 12); j<(N_d0 + N_w + N_a - 1); j++)
		Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 1][j] = Ey_total2_post2[0][M_post2_q2 + Num18_p2 + Node174_p2 + 1 + (j - N_d0 - N_w - 12)];//3

	for (j = (N_d0 + N_w + N_a - 1); j <= (Nx_22 - N_d0 + 1); j++)
		Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 1][j] = Ey_total2_post2[0][M_post2_q2 + Num18_p2 + Node174_p2 + 5 + (j - N_d0 - N_w - N_a + 1) * 2];//4
	for (j = (Nx_22 - N_d0 + 2); j<(Nx_22); j++)
		Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 1][j] = Ey_total2_post2[0][M_post2_q2 + Num18_p2 + Node_1717_p2 + 1 + (j - Nx_22 + N_d0 - 2)];//5

	// line 8--Num 25_p2 (Num 13_p2)
	for (j = 0; j<(N_d0 - 2); j++)
		Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 2][j] = Ey_total2_post2[0][M_post2_q2 + Num14_p2 + Num15_p2 + Num16_p2 + Num17_p2 + Num18_p2 + j];//1
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 2][N_d0 - 2] = Ey_total2_post2[0][M_post2_q2 + Num14_p2 + Num15_p2 + Num16_p2 + Num17_p2 + Num18_p2 + N_d0 - 1];//2
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 2][N_d0 - 1] = Ey_total2_post2[0][M_post2_q2 + Num14_p2 + Num15_p2 + Num16_p2 + Num17_p2 + Num18_p2 + N_d0 + 1];//2
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 2][N_d0] = Ey_total2_post2[0][M_post2_q2 + Num14_p2 + Num15_p2 + Num16_p2 + Num17_p2 + Num18_p2 + N_d0 + 5];//2
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 2][N_d0 + N_w] = Ey_total2_post2[0][M_post2_q2 + Num14_p2 + Num15_p2 + Num16_p2 + Num17_p2 + Num18_p2 + N_d0 + 8];//2
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 2][N_d0 + N_w + 1] = Ey_total2_post2[0][M_post2_q2 + Num14_p2 + Num15_p2 + Num16_p2 + Num17_p2 + Num18_p2 + Node_13_p2 - 3];//2
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 2][N_d0 + N_w + 2] = Ey_total2_post2[0][M_post2_q2 + Num14_p2 + Num15_p2 + Num16_p2 + Num17_p2 + Num18_p2 + Node_13_p2 - 1];//2

	for (j = (N_d0 + N_w + 3); j<(N_d0 + N_w + 6); j++)
		Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 2][j] = Ey_total2_post2[0][M_post2_q2 + Num14_p2 + Num15_p2 + Num16_p2 + Num17_p2 + Num18_p2 + Node_13_p2 + 1 + (j - N_d0 - N_w - 3)];//3
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 2][N_d0 + N_w + 6] = Ey_total2_post2[0][M_post2_q2 + Num14_p2 + Num15_p2 + Num16_p2 + Num17_p2 + Num18_p2 + Node_13_p2 + 5];//3
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 2][N_d0 + N_w + 7] = Ey_total2_post2[0][M_post2_q2 + Num14_p2 + Num15_p2 + Num16_p2 + Num17_p2 + Num18_p2 + Node_13_p2 + 8];//3
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 2][N_d0 + N_w + 8] = Ey_total2_post2[0][M_post2_q2 + Num14_p2 + Num15_p2 + Num16_p2 + Num17_p2 + Num18_p2 + Node_13_p2 + 12];//3
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 2][N_d0 + N_w + 9] = Ey_total2_post2[0][M_post2_q2 + Num14_p2 + Num15_p2 + Num16_p2 + Num17_p2 + Num18_p2 + Node_13_p2 + 16];//3
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 2][N_d0 + N_w + 10] = Ey_total2_post2[0][M_post2_q2 + Num14_p2 + Num15_p2 + Num16_p2 + Num17_p2 + Num18_p2 + Node133_p2 - 4];//3
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 2][N_d0 + N_w + 11] = Ey_total2_post2[0][M_post2_q2 + Num14_p2 + Num15_p2 + Num16_p2 + Num17_p2 + Num18_p2 + Node133_p2 - 1];//3
	for (j = (N_d0 + N_w + 12); j<(N_d0 + N_w + N_a - 1); j++)
		Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 2][j] = Ey_total2_post2[0][M_post2_q2 + Num14_p2 + Num15_p2 + Num16_p2 + Num17_p2 + Num18_p2 + Node133_p2 + 1 + (j - N_d0 - N_w - 12)];//3

	Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 2][N_d0 + N_w + N_a - 1] = Ey_total2_post2[0][M_post2_q2 + Num14_p2 + Num15_p2 + Num16_p2 + Num17_p2 + Num18_p2 + Node133_p2 + 5];//4
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 2][N_d0 + N_w + N_a] = Ey_total2_post2[0][M_post2_q2 + Num14_p2 + Num15_p2 + Num16_p2 + Num17_p2 + Num18_p2 + Node133_p2 + 7];//4
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 2][N_d0 + N_w + N_a + 1] = Ey_total2_post2[0][M_post2_q2 + Num14_p2 + Num15_p2 + Num16_p2 + Num17_p2 + Num18_p2 + Node134_p2 - 1];//4
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 2][Nx_22 - N_d0 - 1] = Ey_total2_post2[0][M_post2_q2 + Num14_p2 + Num15_p2 + Num16_p2 + Num17_p2 + Num18_p2 + Node134_p2 + 2];//4
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 2][Nx_22 - N_d0] = Ey_total2_post2[0][M_post2_q2 + Num14_p2 + Num15_p2 + Num16_p2 + Num17_p2 + Num18_p2 + Node_1313_p2 - 3];//4
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 2][Nx_22 - N_d0 + 1] = Ey_total2_post2[0][M_post2_q2 + Num14_p2 + Num15_p2 + Num16_p2 + Num17_p2 + Num18_p2 + Node_1313_p2 - 1];//4
	for (j = (Nx_22 - N_d0 + 2); j<(Nx_22); j++)
		Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 2][j] = Ey_total2_post2[0][M_post2_q2 + Num14_p2 + Num15_p2 + Num16_p2 + Num17_p2 + Num18_p2 + Node_1313_p2 + 1 + (j - Nx_22 + N_d0 - 2)];//5

	// line 9--Num 29_p2 (4--Num 9_p2)
	for (j = 0; j<(N_d0 - 2); j++)
		Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 3][j] = Ey_total2_post2[0][M_post2_q3 - Num9_p2 + j];//1
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 3][N_d0 - 2] = Ey_total2_post2[0][M_post2_q3 - Num9_p2 + N_d0 - 1];//2
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 3][N_d0 - 1] = Ey_total2_post2[0][M_post2_q3 - Num9_p2 + N_d0 + 1];//2
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 3][N_d0 + N_w + 1] = Ey_total2_post2[0][M_post2_q3 - Num9_p2 + Node_9_p2 - 3];//2
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 3][N_d0 + N_w + 2] = Ey_total2_post2[0][M_post2_q3 - Num9_p2 + Node_9_p2 - 1];//2

	for (j = (N_d0 + N_w + 3); j<(N_d0 + N_w + 6); j++)
		Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 3][j] = Ey_total2_post2[0][M_post2_q3 - Num9_p2 + Node_9_p2 + 1 + (j - N_d0 - N_w - 3)];//3
	for (j = (N_d0 + N_w + 6); j<(N_d0 + N_w + 12); j++)
		Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 3][j] = Ey_total2_post2[0][M_post2_q3 - Num9_p2 + Node_9_p2 + 5 + (j - N_d0 - N_w - 6) * 2];//3
	for (j = (N_d0 + N_w + 12); j<(N_d0 + N_w + N_a - 1); j++)
		Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 3][j] = Ey_total2_post2[0][M_post2_q3 - Num9_p2 + Node93_p2 + 1 + (j - N_d0 - N_w - 12)];//3

	Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 3][N_d0 + N_w + N_a - 1] = Ey_total2_post2[0][M_post2_q3 - Num9_p2 + Node93_p2 + 5];//4
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 3][N_d0 + N_w + N_a] = Ey_total2_post2[0][M_post2_q3 - Num9_p2 + Node93_p2 + 7];//4
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 3][Nx_22 - N_d0] = Ey_total2_post2[0][M_post2_q3 - Num9_p2 + Node_99_p2 - 3];//4
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 3][Nx_22 - N_d0 + 1] = Ey_total2_post2[0][M_post2_q3 - Num9_p2 + Node_99_p2 - 1];//4
	for (j = (Nx_22 - N_d0 + 2); j<(Nx_22); j++)
		Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 3][j] = Ey_total2_post2[0][M_post2_q3 - Num9_p2 + Node_99_p2 + 1 + (j - Nx_22 + N_d0 - 2)];//5

	// line 10--Num 33_p2 (3--Num 5_p2)
	for (j = 0; j<(N_d0 - 2); j++)
		Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 4][j] = Ey_total2_post2[0][M_post2_q3 + Num8_p2 + Num7_p2 + Num6_p2 + j];//1
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 4][N_d0 - 2] = Ey_total2_post2[0][M_post2_q3 + Num8_p2 + Num7_p2 + Num6_p2 + N_d0 - 1];//2
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 4][N_d0 - 1] = Ey_total2_post2[0][M_post2_q3 + Num8_p2 + Num7_p2 + Num6_p2 + N_d0 + 1];//2
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 4][N_d0] = Ey_total2_post2[0][M_post2_q3 + Num8_p2 + Num7_p2 + Num6_p2 + N_d0 + 5];//2
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 4][N_d0 + N_w] = Ey_total2_post2[0][M_post2_q3 + Num8_p2 + Num7_p2 + Num6_p2 + N_d0 + 8];//2
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 4][N_d0 + N_w + 1] = Ey_total2_post2[0][M_post2_q3 + Num8_p2 + Num7_p2 + Num6_p2 + Node_5_p2 - 3];//2
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 4][N_d0 + N_w + 2] = Ey_total2_post2[0][M_post2_q3 + Num8_p2 + Num7_p2 + Num6_p2 + Node_5_p2 - 1];//2
	for (j = (N_d0 + N_w + 3); j<(N_d0 + N_w + N_a - 1); j++)
		Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 4][j] = Ey_total2_post2[0][M_post2_q3 + Num8_p2 + Num7_p2 + Num6_p2 + Node_5_p2 + 1 + (j - N_d0 - N_w - 3)];//3
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 4][N_d0 + N_w + N_a - 1] = Ey_total2_post2[0][M_post2_q3 + Num8_p2 + Num7_p2 + Num6_p2 + Node_5_p2 + N_a - 2];//4
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 4][N_d0 + N_w + N_a] = Ey_total2_post2[0][M_post2_q3 + Num8_p2 + Num7_p2 + Num6_p2 + Node_5_p2 + N_a];//4
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 4][N_d0 + N_w + N_a + 1] = Ey_total2_post2[0][M_post2_q3 + Num8_p2 + Num7_p2 + Num6_p2 + Node54_p2 - 1];//4
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 4][Nx_22 - N_d0 - 1] = Ey_total2_post2[0][M_post2_q3 + Num8_p2 + Num7_p2 + Num6_p2 + Node54_p2 + 2];//4
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 4][Nx_22 - N_d0] = Ey_total2_post2[0][M_post2_q3 + Num8_p2 + Num7_p2 + Num6_p2 + Node_55_p2 - 3];//4
	Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 4][Nx_22 - N_d0 + 1] = Ey_total2_post2[0][M_post2_q3 + Num8_p2 + Num7_p2 + Num6_p2 + Node_55_p2 - 1];//4
	for (j = (Nx_22 - N_d0 + 2); j<(Nx_22); j++)
		Ey_whole[N_d1 + 1 + (N_s + N_w) * 7 + 4][j] = Ey_total2_post2[0][M_post2_q3 + Num8_p2 + Num7_p2 + Num6_p2 + Node_55_p2 + 1 + (j - Nx_22 + N_d0 - 2)];//5

	// 普通圆孔5-6
	show_middle(8, 4, Ey_whole, Ey_total2_middle);
	show_middle(9, 5, Ey_whole, Ey_total2_middle);
	// 第二个post1 小圆孔 !!
	show_post1(10, 1, Ey_whole, Ey_total2_post1);
	// 普通圆孔7-9
	show_middle(12, 6, Ey_whole, Ey_total2_middle);
	show_middle(13, 7, Ey_whole, Ey_total2_middle);

	//****************************************************
	//********************region_side_out*****************
	//****************************************************
	for (k = 0; k<(N_slot + 1); k++)//right
	{
		for (j = 0; j<(N_d0 - 2); j++)
			Ey_whole[N_d1 + 1 + N_slot2 + (N_s + N_w)*(N_circle + 2) + k][j] = Ey_total2_end[1][(Nx_2 + Nx_3)*k + j];//1
		for (j = (N_d0 - 2); j <= (N_d0 + N_w + 2); j++)
			Ey_whole[N_d1 + 1 + N_slot2 + (N_s + N_w)*(N_circle + 2) + k][j] = Ey_total2_end[1][(Nx_2 + Nx_3)*k + N_d0 - 1 + (j - N_d0 + 2) * 2];//2
		for (j = (N_d0 + N_w + 3); j<(N_d0 + N_w + N_a - 1); j++)
			Ey_whole[N_d1 + 1 + N_slot2 + (N_s + N_w)*(N_circle + 2) + k][j] = Ey_total2_end[1][(Nx_2 + Nx_3)*k + Node_1 + 1 + (j - N_d0 - N_w - 3)];//3
		for (j = (N_d0 + N_w + N_a - 1); j <= (Nx_22 - N_d0 + 1); j++)
			Ey_whole[N_d1 + 1 + N_slot2 + (N_s + N_w)*(N_circle + 2) + k][j] = Ey_total2_end[1][(Nx_2 + Nx_3)*k + Node_1 + N_a - 2 + (j - N_d0 - N_w - N_a + 1) * 2];//4
		for (j = (Nx_22 - N_d0 + 2); j<(Nx_22); j++)
			Ey_whole[N_d1 + 1 + N_slot2 + (N_s + N_w)*(N_circle + 2) + k][j] = Ey_total2_end[1][(Nx_2 + Nx_3)*k + Node_11 + 1 + (j - Nx_22 + N_d0 - 2)];//5
	}
	for (i = 0; i<(N_d1 + 1); i++)
	{
		for (j = 0; j<N_d0; j++)
			Ey_whole[N_d1 + 1 + N_s - 1 + (N_s + N_w)*(N_circle + 2) + i][j] = Ey_total2_end[1][j + Nx_side*i + Nx_2 * 2 + Nx_3 * 1 + Nx_1];
		for (j = N_d0; j<(N_d0 + N_a); j++)
			Ey_whole[N_d1 + 1 + N_s - 1 + (N_s + N_w)*(N_circle + 2) + i][j + N_w + 1] = Ey_total2_end[1][j + Nx_side*i + Nx_2 * 2 + Nx_3 * 1 + Nx_1];
		for (j = (N_d0 + N_a); j<Nx_side; j++)
			Ey_whole[N_d1 + 1 + N_s - 1 + (N_s + N_w)*(N_circle + 2) + i][j + (N_w + 1) * 2] = Ey_total2_end[1][j + Nx_side*i + Nx_2 * 2 + Nx_3 * 1 + Nx_1];
	}
	//将矩阵导出********
	ofstream outfile("result0.txt", ios::app);
	for (i = 0; i<((N_d1 + 1) * 2 + (N_circle + 2)*(N_w + N_s) + N_s - 1); i++)
	{
		for (j = 0; j<Nx_22; j++)
		{
			outfile << Ey_whole[i][j] << endl;
		}
	}
	outfile.close();

	cout << "whole_Ey执行完成！" << endl;
	//显示完整场图!

	////delete Ey_whole
	for (j = 0; j<((N_d1 + 1) * 2 + (N_circle + 2)*(N_w + N_s) + N_s - 1); j++)
	{
		delete[] Ey_whole[j];
		Ey_whole[j] = 0;
	}
	delete[] Ey_whole;
	Ey_whole = 0;

}
//用s参数验证程序正确性
//输入全部场值
//导出场值
void toEvalue(complex *x, int num_total);
//导出坐标，绘图使用的xy坐标
void get_xy_positon(){
	//输出坐标信息
	ofstream out_xy("xy.txt");
	for (int i = 0; i != ((N_d1 + 1) * 2 + (N_circle + 2)*(N_w + N_s) + N_s - 1); i++)
		for (int j = 0; j != 58; j++)
		{
		out_xy << j << " " << i << endl;
		}
	out_xy.close();
}
void show_binc(int k, complex* binc, int num_total);
void main()// DDM-10 region
{
	time_t begin;
	begin = time(&begin);
    int i=0,j=0;
	
	cout << "SIW filter 模型计算：" << endl;

	const int overlapping_nodes = 2 * Num2 + Num3;		//重叠区域
	const int extra_nodes0 = Nx_side + Nx_1;
	const int extra_nodes1 = Num2 + Num3 + Num4 + Num5;	//有影响的节点数
	const int extra_nodes2 = extra_nodes1 + Num6 + Num7_p2;	//有影响的节点数（包含大孔的复杂区域）

	int N = N_matrix;	//压缩矩阵维数
	int M = (M_end - overlapping_nodes) + (M_middle - overlapping_nodes)* N_middle + (M_post1 - overlapping_nodes)*N_post1
		+ (M_post2 - overlapping_nodes)*N_post2 + M_end;	//节点总数

	//各个包含重叠分区矩阵的维数
	int M0 = M_end + extra_nodes1;
	int N0 = M_end;
	int M1 = M_middle + 2 * extra_nodes1;
	int N1 = M_middle;
	int M1_1 = M_middle + extra_nodes0 + extra_nodes1;
	int N1_1 = M_middle;
	int M1_5 = M_middle + extra_nodes1 + extra_nodes2;
	int N1_5 = M_middle;
	int M1_7 = M1_5;
	int N1_7 = N1_5;
	int M1_11 = M_middle + extra_nodes0 + extra_nodes1;
	int N1_11 = M_middle;
	int M1_post1 = M_post1 + 2 * extra_nodes1;
	int N1_post1 = M_post1;
	int M1_post2 = M_post2 + 2 * extra_nodes1;
	int N1_post2 = M_post2;
	int M2 = M_end + extra_nodes1;
	int N2 = M_end;

	//每个单独的分区矩阵的节点数
	int zone_N0 = M_end - overlapping_nodes;
	int zone_N1 = M_middle - 2 * overlapping_nodes;
	int zone_N1_post1 = M_post1 - 2 * overlapping_nodes;
	int zone_N1_post2 = M_post2 - 2 * overlapping_nodes;
	int zone_N2 = M_end - overlapping_nodes;

	//计算的起始点
	int startpos[13];	
	startpos[0] = 0;																				//end
	startpos[1] = startpos[0] + zone_N0 - extra_nodes0;												
	startpos[2] = startpos[1] + extra_nodes0 + overlapping_nodes + zone_N1 - extra_nodes1;
	startpos[3] = startpos[2] + extra_nodes1 + overlapping_nodes + zone_N1 - extra_nodes1;			//post1
	startpos[4] = startpos[3] + extra_nodes1 + overlapping_nodes + zone_N1_post1 - extra_nodes1;
	startpos[5] = startpos[4] + extra_nodes1 + overlapping_nodes + zone_N1  - extra_nodes1;
	startpos[6] = startpos[5] + extra_nodes1 + overlapping_nodes + zone_N1  - extra_nodes1;			//post2
	startpos[7] = startpos[6] + extra_nodes1 + overlapping_nodes + zone_N1_post2  - extra_nodes2;
	startpos[8] = startpos[7] + extra_nodes2 + overlapping_nodes + zone_N1 - extra_nodes1;
	startpos[9] = startpos[8] + extra_nodes1 + overlapping_nodes + zone_N1  - extra_nodes1;			//post1
	startpos[10] = startpos[9] + extra_nodes1 + overlapping_nodes + zone_N1_post1  - extra_nodes1;
	startpos[11] = startpos[10] + extra_nodes1 + overlapping_nodes + zone_N1  - extra_nodes1;
	startpos[12] = startpos[11] + extra_nodes1 + overlapping_nodes + zone_N1 - extra_nodes1;		//end
	cout << "startpos of computing:" << endl;
	for (int i = 0; i != 13; i++)
		cout << setw(2) << i << "/" << setw(5) << startpos[i] << endl;


	//显示完整场值图 从文件中读取
	show_value(M);
	complex zero(0, 0);
	complex* Binc = new complex[M];
	for (int i = 0; i != M; i++)
		Binc[i] = zero;

	//添加激励源
	scan_binc(Binc, M);
	//初试并设置分区矩阵参数 未考虑每个单独的区域！！！！
	//A_* 表示系数，b_* 表示位置
	complex **A_end0 = new complex*[M0];
	int** b_end0 = new int*[M0];
	for (int i = 0; i != M0; i++)
	{
		A_end0[i] = new complex[N];
		b_end0[i] = new int[N];
	}

	complex **A_middle = new complex*[M1];
	int** b_middle = new int*[M1];
	for (int i = 0; i != M1; i++)
	{
		A_middle[i] = new complex[N];
		b_middle[i] = new int[N];
	}

	complex **A_middle_1 = new complex*[M1_1];
	int** b_middle_1 = new int*[M1_1];
	for (int i = 0; i != M1_1; i++)
	{
		A_middle_1[i] = new complex[N];
		b_middle_1[i] = new int[N];
	}

	complex **A_middle_5 = new complex*[M1_5];
	int** b_middle_5 = new int*[M1_5];
	for (int i = 0; i != M1_5; i++)
	{
		A_middle_5[i] = new complex[N];
		b_middle_5[i] = new int[N];
	}

	complex **A_middle_7 = new complex*[M1_7];
	int** b_middle_7 = new int*[M1_7];
	for (int i = 0; i != M1_7; i++)
	{
		A_middle_7[i] = new complex[N];
		b_middle_7[i] = new int[N];
	}

	complex **A_middle_11 = new complex*[M1_11];
	int** b_middle_11 = new int*[M1_11];
	for (int i = 0; i != M1_11; i++)
	{
		A_middle_11[i] = new complex[N];
		b_middle_11[i] = new int[N];
	}

	complex **A_middle_post1 = new complex*[M1_post1];
	int** b_middle_post1 = new int*[M1_post1];
	for (int i = 0; i != M1_post1; i++)
	{
		A_middle_post1[i] = new complex[N];
		b_middle_post1[i] = new int[N];
	}
	
	complex **A_middle_post2 = new complex*[M1_post2];
	int** b_middle_post2 = new int*[M1_post2];
	for (int i = 0; i != M1_post2; i++)
	{
		A_middle_post2[i] = new complex[N];
		b_middle_post2[i] = new int[N];
	}
	
	complex **A_end1 = new complex*[M2];
	int** b_end1 = new int*[M2];
	for (int i = 0; i != M2; i++)
	{
		A_end1[i] = new complex[N];
		b_end1[i] = new int[N];
	}

	//获取系数矩阵(post1的左右两边没有处理）
	scan_coef_Matrix_end1(A_end0, M_end, b_end0, M0, N0,0);										//end
	scan_coef_Matrix_middle_1(A_middle_1, M_middle, b_middle_1, M1_1, N1_1, extra_nodes0);
	scan_coef_Matrix_middle(A_middle, M_middle, b_middle, M1, N1, extra_nodes1);
	scan_coef_Matrix_post1(A_middle_post1, M_post1, b_middle_post1, M1_post1, N1_post1, extra_nodes1);			//post1
	scan_coef_Matrix_middle_5(A_middle_5, M_middle, b_middle_5, M1_5, N1_5, extra_nodes1);
	scan_coef_Matrix_post2(A_middle_post2, M_post2, b_middle_post2, M1_post2, N1_post2, extra_nodes1);			//post2
	scan_coef_Matrix_middle_7(A_middle_7, M_middle, b_middle_7, M1_7, N1_7, extra_nodes2);
	scan_coef_Matrix_middle_11(A_middle_11, M_middle, b_middle_11, M1_11, N1_11, extra_nodes1);
	scan_coef_Matrix_end2(A_end1, M_end, b_end1, M2, N2, extra_nodes1);										//end

	//中间存储变量 X_* 表示最小二乘法的结果
	complex* X_end0 = new complex[N0];
	complex* X_middle = new complex[N1];
	complex* X_middle_post1 = new complex[N1_post1];
	complex* X_middle_post2 = new complex[N1_post2];
	complex* X_end1 = new complex[N2];

	// B_* 表示最小二乘法的右端项
	complex* B_end0 = new complex[M0];
	complex* B_middle = new complex[M1];
	complex* B_middle_1 = new complex[M1_1];			//1,11 使用
	complex* B_middle_5 = new complex[M1_5];			//5,7 使用
	complex* B_middle_post1 = new complex[M1_post1];
	complex* B_middle_post2 = new complex[M1_post2];
	complex* B_end1 = new complex[M2];

	//存储的最终结果 x 每个节点的数值结果
	complex* x = new complex[M];

	//收敛条件设置
	double bnorm = 0.0;
	bnorm = VectorNorm(M, Binc);
	double eps = bnorm*1e-3;
	int times = 0;
	double sinta;
	double sinta1 = 1.75;
	//sinta1 = 1.0; 

	cout << "收敛目标：" << eps << endl;
	cout << "开始投影:" << endl;
	//迭代过程 未完成！！！！！
	while (bnorm > eps)
	{
		//快速投影系数
		if (times == 0)
			sinta = 1;
		else
			sinta = sinta1;

		int k = 0;
		//调用投影函数
		projective_process(A_end0, b_end0, B_end0, X_end0, Binc, x, startpos[0], 0, sinta, M0, N0);		//end
		//show_binc(k++, Binc, M);
		projective_process(A_middle_1, b_middle_1, B_middle_1, X_middle, Binc, x, startpos[1], extra_nodes0, sinta, M1_1, N1_1);
		//show_binc(k++, Binc, M);
		projective_process(A_middle, b_middle, B_middle, X_middle, Binc, x, startpos[2], extra_nodes1, sinta, M1, N1);
		//show_binc(k++, Binc, M);
		projective_process(A_middle_post1, b_middle_post1, B_middle_post1, X_middle_post1, Binc, x, startpos[3], extra_nodes1, sinta, M1_post1, N1_post1);		//post1
		//show_binc(k++, Binc, M);
		projective_process(A_middle, b_middle, B_middle, X_middle, Binc, x, startpos[4], extra_nodes1, sinta, M1, N1);
		//show_binc(k++, Binc, M);
		projective_process(A_middle_5, b_middle_5, B_middle_5, X_middle, Binc, x, startpos[5], extra_nodes1, sinta, M1_5, N1_5);
		//show_binc(k++, Binc, M);
		projective_process(A_middle_post2, b_middle_post2, B_middle_post2, X_middle_post2, Binc, x, startpos[6], extra_nodes1, sinta, M1_post2, N1_post2);		//post2
		//show_binc(k++, Binc, M);
		projective_process(A_middle_7, b_middle_7, B_middle_5, X_middle, Binc, x, startpos[7], extra_nodes2, sinta, M1_7, N1_7);
		//show_binc(k++, Binc, M);
		projective_process(A_middle, b_middle, B_middle, X_middle, Binc, x, startpos[8], extra_nodes1, sinta, M1, N1);
		//show_binc(k++, Binc, M);
		projective_process(A_middle_post1, b_middle_post1, B_middle_post1, X_middle_post1, Binc, x, startpos[9], extra_nodes1, sinta, M1_post1, N1_post1);		//post1
		//show_binc(k++, Binc, M);
		projective_process(A_middle, b_middle, B_middle, X_middle, Binc, x, startpos[10], extra_nodes1, sinta, M1, N1);
		//show_binc(k++, Binc, M);
		projective_process(A_middle_11, b_middle_11, B_middle_1, X_middle, Binc, x, startpos[11], extra_nodes1, sinta, M1_11, N1_11);
		//show_binc(k++, Binc, M);
		projective_process(A_end1, b_end1, B_end1, X_end1, Binc, x, startpos[12], extra_nodes1, sinta, M2, N2);		//end
		//show_binc(k++, Binc, M);

		bnorm = 0.0;
		bnorm = VectorNorm(M, Binc);

		cout << ++times << " cond:" << bnorm << endl;
	}

	ofstream out_e("Evalue.txt");
	for (int i = 0; i != M; i++)
		out_e << x[i].amplitude() << endl;
	out_e.close();
	ofstream out_all_x("x_all.txt");
	for (int i = 0; i != M; i++)
		out_all_x << x[i].real << " " << x[i].image << endl;
	out_all_x.close();
	
	//用s参数验证程序正确性***************************
	double s11 = 0, s21 = 0;//s11,s12
	complex x1[N_a + 2], x2, x3;//x1为左端反射波，x2为右端透射波
	x1[0].real = 0;
	x1[0].image = 0;
	x1[N_a + 1].real = 0;
	x1[N_a + 1].image = 0;
	for (j = N_d0; j<(N_d0 + N_a); j++)
	{
		x1[j - N_d0 + 1].real = x[j].real;
		x1[j - N_d0 + 1].image = x[j].image;
	}
	complex Einc2[N_a+2];
	for( j=0; j<(N_a+2); j++ )
	{
		Einc2[j].real=sin( pi*j*4*dx/a );
		Einc2[j].image=0;
		x1[j].real=x1[j].real-Einc2[j].real;
		x1[j].image=x1[j].image-Einc2[j].image;
	}
	x2.real = 0;
	x2.image = 0;
	x3.real = 0;
	x3.image = 0;
	for (j = 0; j<(N_a + 1); j++)
	{
		x2.real = (x1[j].real + x1[j + 1].real) / 2;
		x2.image = (x1[j].image + x1[j + 1].image) / 2;
		x3.real = x3.real + x2.real * 4 * dx;
		x3.image = x3.image + x2.image * 4 * dx;
		//	cout<<x2.real<<"  "<<x2.image<<endl;
	}
	s11 = sqrt(x3.real*x3.real + x3.image*x3.image) / ((a / pi) * (1 - cos(pi)));
	cout << "s11 = " << s11 << " -- " << (20 * log10(s11)) << " dB" << endl;
	//求 s21
	for (j = N_d0; j<(N_d0 + N_a); j++)
	{
		x1[j - N_d0 + 1].real = x[j + M - Nx_side].real;
		x1[j - N_d0 + 1].image = x[j + M - Nx_side].image;
	}
	x2.real = 0;
	x2.image = 0;
	x3.real = 0;
	x3.image = 0;
	for (j = 0; j<(N_a + 1); j++)
	{
		x2.real = (x1[j].real + x1[j + 1].real) / 2;
		x2.image = (x1[j].image + x1[j + 1].image) / 2;
		x3.real = x3.real + x2.real * 4 * dx;
		x3.image = x3.image + x2.image * 4 * dx;
	}
	s21 = sqrt(x3.real*x3.real + x3.image*x3.image) / ((a / pi) * (1 - cos(pi)));
	cout << "s21 = " << s21 << " -- " << (20 * log10(s21)) << " dB" << endl;
	cout << "s11*s11+s21*s21 = " << (s11*s11 + s21*s21) << endl;
	//输出频率点及S参数
	ofstream outfile3("result1.txt", ios::app);
	outfile3 << "频率：" << (1200 / 1.483 / N_dx) << endl;
	outfile3 << "s21 = " << (20 * log10(s21)) << " dB" << endl;
	outfile3 << "s11 = " << (20 * log10(s11)) << " dB" << endl;
	outfile3 << "s21=" << s21 << " s11=" << s11 << endl;
	outfile3.close();
	
	

	//释放内存
	delete Binc;
	for (int i = 0; i != M0; i++)
		delete A_end0[i];
	delete A_end0;
	for (int i = 0; i != N0; i++)
		delete b_end0[i];
	delete b_end0;
	for (int i = 0; i != M1; i++)
		delete A_middle[i];
	delete A_middle;
	for (int i = 0; i != M1; i++)
		delete b_middle[i];
	delete b_middle;
	for (int i = 0; i != M1_1; i++)
		delete A_middle_1[i];
	delete A_middle_1;
	for (int i = 0; i != M1_1; i++)
		delete b_middle_1[i];
	delete b_middle_1;
	for (int i = 0; i != M1_5; i++)
		delete A_middle_5[i];
	delete A_middle_5;
	for (int i = 0; i != M1_5; i++)
		delete b_middle_5[i];
	delete b_middle_5;
	for (int i = 0; i != M1_7; i++)
		delete A_middle_7[i];
	delete A_middle_7;
	for (int i = 0; i != M1_7; i++)
		delete b_middle_7[i];
	delete b_middle_7;
	for (int i = 0; i != M1_11; i++)
		delete A_middle_11[i];
	delete A_middle_11;
	for (int i = 0; i != M1_11; i++)
		delete b_middle_11[i];
	delete b_middle_11;
	for (int i = 0; i != M1_post1; i++)
		delete A_middle_post1[i];
	delete A_middle_post1;
	for (int i = 0; i != M1_post1; i++)
		delete b_middle_post1[i];
	delete b_middle_post1;
	for (int i = 0; i != M1_post2; i++)
		delete A_middle_post2[i];
	delete A_middle_post2;
	for (int i = 0; i != M1_post2; i++)
		delete b_middle_post2[i];
	delete b_middle_post2;
	for (int i = 0; i != M2; i++)
		delete A_end1[i];
	delete A_end1;
	for (int i = 0; i != M2; i++)
		delete b_end1[i];
	delete b_end1;
	delete X_end0;
	delete X_middle;
	delete X_middle_post1;
	delete X_middle_post2;
	delete X_end1;
	delete B_end0;
	delete B_middle;
	delete B_middle_1;
	delete B_middle_5;
	delete B_end1;
	delete B_middle_post1;
	delete B_middle_post2;
	delete x;

	//计算所用时间并输出
	time_t end;
	end = time(&end);
    cout <<"程序所用时间："<<( end - begin )<<"(s)"<<endl;
	system("pause");
}




