#include "Tracer.h"
#include "stdio.h"
#include "cmath"
#include <iostream>
void main(int argc, char** argv)
{

  CTracer tracer;
  CScene scene;
  double pi = 3.14159265358979323;
  int xRes = 1024;  // Default resolution
  int yRes = 768;
  double alpha, beta;
  double Xpos, Ypos, Zpos, Xup, Yup, Zup;
  double Xr, Yr, Zr, Xv, Yv, Zv, M1, G1, c1, k1, flag;
	 xRes = 512;
	 yRes = 512;
	alpha = 1.57;
	Xpos = 1e+11;
	Ypos = 0;
	Zpos = 1e+10;
	Xup = -1;
	Yup = 0;
	Zup = 6;
	Xr = 0;
	Yr = 1;
	Zr = 0;
	Xv = -6;
	Yv = 0;
	Zv = 1;
	M1 = 8.57e+36; //mass
	G1 = 6.674e-11;
	c1 = 3e+8; //3e+8
	k1 = 3; //rad of ring = k * R_hole
	flag = 0;
  if(argc == 2) // There is input file in parameters
  {
    FILE* file = fopen(argv[1], "r");
    if(file)
    {
      int xResFromFile = 0;
      int yResFromFile = 0;
	  double alphaFile = 0, XposFile = 0, YposFile = 0, ZposFile = 0, XupFile = 0, YupFile = 0, ZupFile = 0;
	  double XrFile = 0, YrFile = 0, ZrFile = 0, XvFile = 0, YvFile = 0, ZvFile = 0, MFile = 0, GFile = 0, cFile = 0, kFile = 0, flagF = 0;
	  
      if(fscanf(file, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &xResFromFile, &yResFromFile, &alphaFile, 
		  &XposFile, &YposFile, &ZposFile, &XupFile, 
		  &YupFile, &ZupFile, &XrFile, &YrFile, &ZrFile, &XvFile, &YvFile, &ZvFile, &MFile, &GFile, &cFile, &kFile, &flagF) == 20)
      {
        xRes = xResFromFile;
        yRes = yResFromFile;
		alpha = alphaFile;
		Xpos = XposFile;
		Ypos = YposFile;
		Zpos = ZposFile;
		Xup = XupFile;
		Yup = YupFile;
		Zup = ZupFile;
		Xr = XrFile;
		Yr = YrFile;
		Zr = ZrFile;
		Xv = XvFile;
		Yv = YvFile;
		Zv = ZvFile;
		M1 = MFile; //mass
		G1 = GFile;
		c1 = cFile; //3e+8
		k1 = kFile; //rad of ring = k * R_hole
		flag = flagF;
      }
      else
        printf("Invalid config format! Using default parameters.\r\n");

      fclose(file);
    }
    else
      printf("Invalid config path! Using default parameters.\r\n");
  }
  else
    printf("No config! Using default parameters.\r\n");
  beta =  2 * atan(((double(yRes)) / (double(xRes))) * tan(alpha / 2.0));
  tracer.m_pScene = &scene;
  tracer.RenderImage(xRes, yRes, alpha, beta, Xpos, Ypos, Zpos, Xup, Yup, Zup, Xr, Yr, Zr, Xv, Yv, Zv, M1, G1, c1, k1, flag);
  tracer.SaveImageToFile("Result.png");
}