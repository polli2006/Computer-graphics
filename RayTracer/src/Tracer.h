#pragma once
#include <iostream>
#include "glm/glm.hpp"
#include "Types.h"
#include "Scene.h"

#include "string"
#include "atlimage.h"

class CTracer
{
public:
  double k, R, M, c, G;
  CImage *pImage, *pBack;
  double VectorMult(glm::vec3 a, glm::vec3 b);
  double RayDisk(SRay ray, double RSphere, glm::vec3 Spos);
  double RaySphere(SRay ray, glm::vec3 spos, double r);
  SRay MakeRay(glm::uvec2 pixelPos);  // Create ray for specified pixel
  glm::vec3 TraceRay(SRay ray); // Trace ray, compute its color
  void RenderImage(int xRes, int yRes, double alpha, double beta, double Xpos, double Ypos, double Zpos, double Xup, 
	double Yup, double Zup, double Xr, double Yr, double Zr, double Xv, double Yv, double Zv, double M, double G, double c, double k, double flag);
  void SaveImageToFile(std::string fileName);
  CImage* LoadImageFromFile(std::string fileName);

public:
  SCamera m_camera;
  CScene* m_pScene;
};