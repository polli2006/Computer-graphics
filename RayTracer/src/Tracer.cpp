#include "Tracer.h"
//using namespace std;
using namespace glm;

double CTracer::VectorMult(vec3 a, vec3 b)
{
	vec3 v = vec3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
	return length(v);
}
double minim(double a, double b)
{
	if (a <= b) return a;
	else return b;
}
double maxim(double a, double b)
{	
	if (a >= b) return a;
	else return b;
}
double CTracer::RayDisk(SRay ray, double RSphere, vec3 Spos)
{
	double Rdisk = k * RSphere;
	ray.m_dir = normalize(ray.m_dir);
	double lx = ray.m_dir.x, ly = ray.m_dir.y, lz = ray.m_dir.z;
	if (fabs(lz) < 0.00000001)
	{
		return -1; //parallel
	}
	double t0 = -ray.m_start.z / lz;
	if (t0 < 0) return t0;
	vec3 coord = ray.m_start + float(t0) * vec3(lx, ly, lz);
	double scal = dot(coord, coord);
	if (!(sqrt(scal) >= RSphere && sqrt(scal) <= Rdisk)) t0 = -1;
	return t0;
}

double CTracer::RaySphere(SRay ray,  vec3 Spos, double RSphere)
{
	ray.m_dir = normalize(ray.m_dir);
	double b = dot(ray.m_start - Spos, ray.m_dir);
	double c0 = dot(ray.m_start - Spos, ray.m_start - Spos) - R * R;
	double d = b * b - c0;
	double t = -1;
	if (d >= 0)
	{
		double t1 = -b + sqrt(d);
		double t2 = -b - sqrt(d);
		double tmin  = minim(t1,t2);
		double tmax = maxim(t1,t2);
		t = (tmin >= 0) ? tmin : tmax;
	}
	return t;
}
SRay CTracer::MakeRay(glm::uvec2 pixelPos)
{
	double W, H;
	W =  m_camera.m_resolution.x ;
	H =  m_camera.m_resolution.y ;
	double k1 = (pixelPos.x + 0.5) / W - 0.5;
	double k2 = (pixelPos.y + 0.5) / H - 0.5;
	vec3 newRay = m_camera.m_forward + (float(k1)) * m_camera.m_right + (float(k2)) * m_camera.m_up;
	SRay r;
	r.m_start = m_camera.m_pos;
	r.m_dir = newRay;
    return r; 
}

glm::vec3 CTracer::TraceRay(SRay ray)
{
	double RSphere = R;
	vec3 Spos = vec3(0, 0, 0);
	double td, ts, deltaT = 2, rabs, leng;
	vec3 r0, a0, color = vec3(0,0,255), dif, normdir, disk, new_start, new_dir;
	int iter = 0, limit = round ((length(m_camera.m_pos) + R * k) / (c * deltaT));
	double dist, distcam, xd, yd, k0;
	dif = m_camera.m_pos - Spos;
	int color0 = 2, sw, sh, i, j;
	double newX, newY;
	unsigned char *pixel_addr;
	distcam = sqrt(dot(dif, dif));
	do
	{
		dif = ray.m_start - Spos;
		dist = length(dif);
		//if ((dist >= distcam + 1e11) && (iter > 0)) break;
		//if (iter > 500) break;
		ray.m_dir = normalize(ray.m_dir) * float(c);	
		r0 = Spos - ray.m_start; 
		rabs = length(r0);
		k0 = G * M / (rabs * rabs * rabs);
		a0 = r0 * float(k0);
		new_start = ray.m_start + ray.m_dir * float(deltaT) + a0 * float(deltaT * deltaT / 2.0); 
		new_dir = ray.m_dir + a0 * float(deltaT);
		normdir = normalize(ray.m_dir);
		ts = RaySphere(ray, Spos, R);
		td = RayDisk(ray, RSphere, Spos);

		if ((VectorMult(new_dir, ray.m_dir) / (length(new_dir) * length(ray.m_dir)) < 1e-3) &&(ts < 0)&&(td < 0)) 
		{
			++iter;
		}
		else iter = 0;
		
		if (iter >= limit) break;
		if ((dot(ray.m_start, ray.m_start) - R * R) * (dot(new_start, new_start) - R * R) < 0) //in sphere
		{
			/*ts = RaySphere(ray, Spos, R);
			if (ts >= 0)*/ color = vec3(0, 0, 0);
			color0 = 0;
			break;
		}
		if ((ray.m_start.z * new_start.z <= 0) && (ray.m_start.x * ray.m_start.x + ray.m_start.y * ray.m_start.y <= R * R * k * k ||
			new_start.x * new_start.x + new_start.y * new_start.y <= R * k * R * k))
		{
			//color = vec3(0, 255, 0);
			if (td >= 0)
			{
				disk = ray.m_start + float(td) * normdir;
				xd = disk.x;
				yd = disk.y;
				newX = xd / (2 * R * k) + 0.5;
				newY = yd / (2 * R * k) + 0.5;			
				sw = pImage->GetWidth();
				sh = pImage->GetHeight();
				i = floor(newX * sw) ;
				j = floor(newY * sh);
				if (i < 0 ) i = 0;
				if (i > pImage->GetWidth()) 
					i = pImage->GetWidth();
				if (j < 0 ) j = 0;
				if (j >  pImage->GetHeight()) 
					j =  pImage->GetHeight();
				pixel_addr = (unsigned char *)pImage->GetPixelAddress(i, j);
				if (pixel_addr[3] < 1e-6) 
				{
					td = -1;  //next
				}
				else 	
				{
					color0 = 1;
					color = vec3(float(*(pixel_addr + 2)) / 255.0f, float(*(pixel_addr + 1)) / 255.0f, float(*(pixel_addr)) / 255.0f);
					break;
				}
			}
		}
		ray.m_start += ray.m_dir * float(deltaT) + a0 * float(deltaT * deltaT / 2.0); 
		ray.m_dir += a0 * float(deltaT);
	}
	while (/*(ts < 0) && (td < 0)*/1);
	double psi, hi, pi = 3.141593;
	sw = pBack->GetWidth();
	sh = pBack->GetHeight();
	double H = minim(sw, sh);
	if (color0 == 2)
	{
		normdir = normalize(ray.m_dir);
		psi = atan2(normdir.x, normdir.y) + pi;
		hi = asin(normdir.z) + pi / 2.0;
		psi = psi * H / pi;
		hi = hi * H / pi;
		if (psi > 2 * H) {psi = 2 * H - 1e-6;}
		if (psi < 0) psi = 0;
		if (hi > H) {hi = H - 1e-6;}
		if (hi < 0) hi = 0;
		i = floor(psi);
		j = floor(hi);
		pixel_addr = (unsigned char *)pBack->GetPixelAddress(i, j);
		color = vec3(float(*(pixel_addr + 2)) / 255.0f, float(*(pixel_addr + 1)) / 255.0f, float(*(pixel_addr)) / 255.0f);
	}

    return color; 
}

void CTracer::RenderImage(int xRes, int yRes, double alpha, double beta, double Xpos, double Ypos, double Zpos, double Xup, 
	double Yup, double Zup, double Xr, double Yr, double Zr, double Xv, double Yv, double Zv, double M1, double G1, double c1, double k1, double flag)
{

  // Rendering
  k = k1;
  M = M1;
  G = G1;
  c = c1;
  R = 2 * G1 * M1 / (c1 * c1);
  double dist = yRes / (2.0 * tan(beta / 2.0));
  m_camera.m_resolution = uvec2(xRes, yRes);
  m_camera.m_pixels.resize(xRes * yRes);
  m_camera.m_viewAngle.x = alpha; 
  m_camera.m_viewAngle.y = beta; 
  m_camera.m_pos = vec3(Xpos, Ypos, Zpos);
  m_camera.m_up = vec3(Xup, Yup, Zup);
  m_camera.m_forward = vec3(Xv, Yv, Zv);
  m_camera.m_right = vec3(Xr, Yr, Zr);
  m_camera.m_up *= float(1 / length(m_camera.m_up)) * float(yRes);
  m_camera.m_forward *= float(1 / length(m_camera.m_forward)) * float(dist);
  m_camera.m_right *= float(1 / length(m_camera.m_right)) * float(xRes);

  pImage = LoadImageFromFile("data/disk_32.png");
  pBack = LoadImageFromFile("data/stars.jpg");
  int i1, j1;
  if (flag == 1)
  {
	 #pragma omp parallel for
	 for( i1 = 0; i1 < yRes; i1++)
     #pragma omp parallel for
	    for( j1 = 0; j1 < xRes; j1++)
		{
		  SRay ray = MakeRay(uvec2(j1, i1));
		  m_camera.m_pixels[i1 * xRes + j1] = TraceRay(ray);
		}
  }
  else
  {
	 for( i1 = 0; i1 < yRes; i1++)
		 for( j1 = 0; j1 < xRes; j1++)
    {
      SRay ray = MakeRay(uvec2(j1, i1));
      m_camera.m_pixels[i1 * xRes + j1] = TraceRay(ray);
	 }
  }
 /* CImage* pImage = LoadImageFromFile("data/disk_32.png");
  if(pImage->GetBPP() == 32)

  {
    auto pData = (unsigned char*)pImage->GetBits();
    auto pCurrentLine = pData;
    int pitch = pImage->GetPitch();

    for(int i = 0; i < pImage->GetHeight(); i++) // Image lines
    {
      for(int j = 0; j < pImage->GetWidth(); j++) // Pixels in line
      {
        unsigned char b = pCurrentLine[j * 4];
        unsigned char g = pCurrentLine[j * 4 + 1];
        unsigned char r = pCurrentLine[j * 4 + 2];
        unsigned char alpha = pCurrentLine[j * 4 + 3];
      }
	  pCurrentLine += pitch;
    }
  } */
}

void CTracer::SaveImageToFile(std::string fileName)
{
  CImage image;

  int width = m_camera.m_resolution.x;
  int height = m_camera.m_resolution.y;

  image.Create(width, height, 24);
    
	int pitch = image.GetPitch();
	unsigned char* imageBuffer = (unsigned char*)image.GetBits();

	if (pitch < 0)
	{
		imageBuffer += pitch * (height - 1);
		pitch =- pitch;
	}

	int i, j;
	int imageDisplacement = 0;
	int textureDisplacement = 0;

	for (i = 0; i < height; i++)
	{
    for (j = 0; j < width; j++)
    {
      vec3 color = m_camera.m_pixels[textureDisplacement + j];

      imageBuffer[imageDisplacement + j * 3] = clamp(color.b, 0.0f, 1.0f) * 255.0f;
      imageBuffer[imageDisplacement + j * 3 + 1] = clamp(color.g, 0.0f, 1.0f) * 255.0f;
      imageBuffer[imageDisplacement + j * 3 + 2] = clamp(color.r, 0.0f, 1.0f) * 255.0f;
    }

		imageDisplacement += pitch;
		textureDisplacement += width;
	}

  image.Save(fileName.c_str());
	image.Destroy();
}

CImage* CTracer::LoadImageFromFile(std::string fileName)
{
  CImage* pImage = new CImage;

  if(SUCCEEDED(pImage->Load(fileName.c_str())))
  return pImage;
  else
  {
    delete pImage;
    return NULL;
  }
}