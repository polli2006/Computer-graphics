#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <initializer_list>
#include <limits>
#include <vector>
#include <tuple>
#include <memory>
#include <queue> 

using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::ofstream;
using std::vector;
using std::tuple;                                     
using std::tie;
using std::make_tuple;
using std::shared_ptr;
using std::make_shared;
using std::queue;
using std::stoi;

#include "io.h"
#include "matrix.h"
#include "MyObject.h"

Matrix <int> binar(const Image&);
Matrix <int> components_search(Matrix <int> m);
Matrix <int> black_noise (Matrix <int> m);
tuple < vector <int>, vector <int>, vector <int>>
    center_square (Matrix <int> m);
tuple <vector <int>, vector <int>>
    radius(vector <int>, vector <int>, Matrix <int>);
int sqr(int a);
tuple <int, Matrix<int>, int, int, int> start_suit_gear(const Matrix<int> &, const Image&);
tuple <int, Matrix <int>> suit_gear(const Matrix <int> &, const Matrix<int> &);
void distance_rows(Matrix <double>);
void distance_cols(Matrix <double>);
int cogs(int, int, int, const Matrix <int> &, int, int);

double border = 127.5;
double inf = 6000000;   
int mm, nn, mm1, nn1;
int count; // 0 = fon; 1 = comp1; ... ; n = delete;
int gear_cntr_x, gear_cntr_y, ax_cntr_x, ax_cntr_y, ax1_cntr_x, ax1_cntr_y, ax_numb, ax1_numb;
int image_name, red, green, blue;

tuple<int, vector<shared_ptr<IObject>>, Image>
repair_mechanism(const Image& in, const Image& im1, const Image& im2, const Image& im3)
{

    //first image
    mm = int(in.n_rows);
    nn = int(in.n_cols);
    Matrix <int> comp = binar(in);
    comp = components_search(comp);
    comp = black_noise(comp);
    vector <int> sq;
    vector <int> x;
    vector<int> y;
    tie(sq, x, y)= center_square(comp);

    vector <int> max_rad;
    vector <int> min_rad;
    tie(max_rad, min_rad) = radius(x, y, comp);

    int result_idx = 0;
    auto object_array = vector<shared_ptr<IObject>>();
    Image out(mm, nn);
    int f = 0;
    for (int i = 1; i < count; ++ i)
    	{
    		if (max_rad[i] - min_rad[i] <= 2)
    		{
    			f = 1;
    		}
    	} 

    if (f) //it will be a base part
    {
        ax_cntr_x = ax1_cntr_x;
        ax_cntr_y = ax1_cntr_y;
        ax_numb = ax1_numb;

        //new images, 0001_1.bmp etc

        mm1 = mm; //permanent
        nn1 = nn;
        int count1 = count;
        int flag[4], squ[4], maxr[4], minr[4];
        Matrix <int> copy[4];
        
        tie(flag[1], copy[1], squ[1], maxr[1], minr[1]) = start_suit_gear(comp, im1); //center
        tie(flag[2], copy[2], squ[2], maxr[2], minr[2]) = start_suit_gear(comp, im2);
        tie(flag[3], copy[3], squ[3], maxr[3], minr[3]) = start_suit_gear(comp, im3);
        int maxsq = 0, maxnum = 0;
        for (int l = 1; l < 4; ++ l)
            if ((flag[l] == 1) && (squ[l] > maxsq))
            {
                maxsq = squ[l];
                maxnum = l;
            }
        comp = copy[maxnum];
        result_idx = maxnum;
        
        //out
        
        mm = mm1; nn = nn1;
        count = count1;
        int numb;
        for (int i = 1; i < count; ++ i)
        {
            numb = cogs(i, x[i], y[i], comp, max_rad[i], min_rad[i]);
            if (i != ax_numb)
            {
                Gear g;
                g.location = make_tuple(y[i], x[i]);
                g.min_r = min_rad[i];
                g.max_r = max_rad[i];
                g.is_broken = false;
                g.num_cogs = numb;
                shared_ptr <IObject> pointer = make_shared <Gear> (g);
                object_array.push_back(pointer);
            }
            else
            {
                Axis a;
                a.location = make_tuple (y[i], x[i]);
                shared_ptr <IObject> pointer = make_shared <Axis>(a);
                object_array.push_back(pointer);
            }
        }
    } 
    else //it will be a 2nd bonus part
    {  
    	int mnr, mxr;

        double max;
        Matrix <double> mt(mm, nn);
        vector <int> dist_center_x(count), dist_center_y(count);
        for (int i = 1; i < count; ++ i)
        {
            for (int j = 0; j < mm; ++ j)
                for (int l = 0; l < nn; ++l)
                {
                    if (comp(j, l) != i)
                        mt(j, l) = 0;
                    else
                        mt(j, l) = inf;
                }

            distance_rows(mt);

            distance_cols(mt);

            max = 0;
            dist_center_x[i] = 1;
            dist_center_y[i] = 1;
            for (int j = 0; j < mm; ++ j)
                for (int l = 0; l < nn; ++l)
                    if (mt(j, l) > max)
                    {
                        max = mt(j, l);
                        dist_center_x[i] = j;
                        dist_center_y[i] = l;
                    } 
        } 
        int broken = 0;
        vector <int> cogs_num (count);
        for (int i = 1; i < count; ++ i)
        {
            if (sqr(dist_center_x[i] - x[i]) + sqr(dist_center_y[i] - y[i]) >= 3)
            {
                x[i] = dist_center_x[i];
                y[i] = dist_center_y[i];
                ax_cntr_x = x[i];
                ax_cntr_y = y[i];
                broken = i;
            }            
        }

        Matrix <int> comp2 = comp.deep_copy();
        int count1;
        int flag[4], squ[4], maxr[4], minr[4];
        Matrix <int> copy[4];
        if (broken == 0)
        {
        	
            double l;
            vector <int> mas(360);
            for (int j = 1; j < count; ++ j)
            {
            	l = (max_rad[j] + min_rad[j]) / 2;
	            for (int i = 0; i < 360; ++ i)
	            {
	            	mas[i] = comp(x[j] + l * cos (i), y[j] + l * sin(i)); //polar
	            }
	            int first, leng, max1, min1;
	            first = 0;
	            leng = 0;
	            max1 = 0;
	            min1 = 360;
	            for (int i = 0; i < 360; ++ i)
	            	if (j != mas[i]) ++ first; //first cog
	            	else break;
	            for (int i = first; i < 360; ++ i)
	            {
	            	if (mas[i] == j)
	            	{
	            		if (leng > max1) max1 = leng;
	            		if ((leng != 0) && (leng < min1)) min1 = leng;
	            		leng = 0;
	            	}
	            	else
	            	{
	            		++ leng;
	            	}
	            }
	            leng += first;
	            if (leng > max1) max1 = leng;
	            if ((leng != 0) && (leng < min1)) min1 = leng;
	            

	            if (max1 / min1 >= 2) broken = j;	
        	}
    	
        }
        	
            for (int i = 0; i < mm; ++ i) 
                for (int j = 0; j < nn; ++j)
                    if (comp(i, j) == broken) comp(i, j) = 0;
            
            mm1 = mm; //permanent
            nn1 = nn;
            count1= count;
            ax_numb = broken;
            tie(flag[1], copy[1], squ[1], maxr[1], minr[1]) = start_suit_gear(comp, im1); //center
            ax_numb = broken;
            tie(flag[2], copy[2], squ[2], maxr[2], minr[2]) = start_suit_gear(comp, im2);
            ax_numb = broken;
            tie(flag[3], copy[3], squ[3], maxr[3], minr[3]) = start_suit_gear(comp, im3);
            ax_numb = broken;
            int maxsq = 0, maxnum = 0;
            for (int l = 1; l < 4; ++ l)
            {
                if ((flag[l] == 1) && (squ[l] > maxsq))
                {
                    maxsq = squ[l];
                    maxnum = l;
                }
            }
            comp = copy[maxnum];
            result_idx = maxnum; 
            mxr = maxr[maxnum];
            mnr = minr[maxnum];
            mm = mm1;
            nn = nn1;
            count = count1;
        

        for (int i = 1; i < count; ++ i)
        {
            if (i == broken)
            {
                min_rad[i] = mnr;
                max_rad[i] = mxr;
                ax_numb = i;
            }
            cogs_num[i] = cogs(i, x[i], y[i], comp2, max_rad[i], min_rad[i]);
        }

        for (int i = 1; i < count; ++ i)
        {
                Gear g;
                g.location = make_tuple(y[i], x[i]);

                if (i != broken)
                {
                    g.is_broken = false;
                }   
                else
                {
                    g.is_broken = true;
                } 
                g.min_r = min_rad[i];
                g.max_r = max_rad[i];   
                g.num_cogs = cogs_num[i];
                shared_ptr <IObject> pointer = make_shared <Gear> (g);
                object_array.push_back(pointer);
        }
    	
	}
    for (int i = 0; i < mm; i++)
        for (int j = 0; j < nn; j++)
        {
            if (comp(i, j) == 0)
            {
                out(i, j) = make_tuple(0, 0, 0);
            }
            else
            {
                out(i, j) = make_tuple(red, green, blue);
            }
        }
    return make_tuple(result_idx, object_array, out);
}

Matrix <int> binar(const Image& in)
{
    int r, g, b;
    Matrix <int> m(mm, nn);
    for (int i = 0; i < mm; ++i)
    {
        for (int j = 0; j < nn; ++j)
        {
            tie(r, g, b) = in(i, j);

            if ((r < border) && (g < border) && (b < border)) 
            {
                m(i, j) = 0; //it is black
            }
            else {
                m(i, j) = -1;  //it is white
                tie(red, green, blue) = in(i, j);
            } 
        }
    }
    return m;
}

Matrix <int> components_search(Matrix <int> m)  //~void 
{
    count = 1;
    int a, b, l,  noise[3][2], ind;   // white noise
    queue <int> elem_que;
    for (int j = 0; j < nn; ++j)
    {
        for (int i = 0; i < mm; ++i)
        {
            for (int k = 0; k < 3; ++k)
                for (int q = 0; q < 2; ++q)
                    noise[k][q] = -1;
            l = 0;
            ind = 0;
            if (m(i, j) == -1)
            {
                elem_que.push(i);
                elem_que.push(j);
                while (elem_que.empty() == false)
                {
                    a = elem_que.front();
                    elem_que.pop();
                    b = elem_que.front();
                    elem_que.pop();
                    ++ind;
                    noise[l][0] = a;
                    noise[l][1] = b;
                    if (l <= 1) ++l;
                    
                    m(a, b) = count;
                    
                    //Look at the neighbours

                    if ((a - 1 >= 0) && (b - 1 >= 0 ) && (m(a - 1, b - 1) == -1))
                    {        
                        elem_que.push(a - 1);
                        elem_que.push(b - 1);
                        m(a - 1, b - 1) = count;
                    }
                    if ((b - 1 >= 0 ) && (m(a, b - 1) == -1))
                    {
                        elem_que.push(a);
                        elem_que.push(b - 1);
                        m(a, b - 1) = count;
                    }
                    if ((a - 1 >= 0) && (m(a - 1, b) == -1))
                    {
                        elem_que.push(a - 1);
                        elem_que.push(b);
                        m(a - 1, b) = count;
                    }
                    if ((a + 1 < mm) && (b + 1 < nn ) && (m(a + 1, b + 1) == -1))
                    {
                        elem_que.push(a + 1);
                        elem_que.push(b + 1);
                        m(a + 1, b + 1) = count;
                    }
                    if ((a + 1 < mm) && (m(a + 1, b) == -1))
                    {
                        elem_que.push(a + 1);
                        elem_que.push(b);
                        m(a + 1, b) = count;
                    }
                    if ((b + 1 < nn) && (m(a, b + 1) == -1))
                    {
                        elem_que.push(a);
                        elem_que.push(b + 1);
                        m(a, b + 1) = count;
                    }
                    if ((a - 1 >= 0 ) && (b + 1 < nn) && (m(a - 1, b + 1) == -1))
                    {
                        elem_que.push(a - 1);
                        elem_que.push(b + 1);
                        m(a - 1, b + 1) = count;
                    }
                    if ((b - 1 >= 0 ) && (a + 1 < mm) && (m(a + 1, b - 1) == -1))
                    {
                        elem_que.push(a + 1);
                        elem_que.push(b - 1);
                        m(a + 1, b - 1) = count;
                    }

                }
                if (ind <= 3) //delete white noise
                {
                     for (int k = 0; k < 3; ++k)
                        if ((noise[k][0] >= 0) && (noise[k][1] >= 0)) 
                            m(noise[k][0], noise[k][1]) = 0; 
                }
                else
                {
                    count ++;
                }
            }
        }
    }
    return m;
}

Matrix <int> black_noise (Matrix <int> m)
{
    int white_sum = 0, black_sum = 0, comp_num[8];
    for (int a = 0; a < mm; ++ a)
        for (int b = 0; b < nn; ++ b)
        {
            white_sum = 0;
            black_sum = 0;
            for (int k = 0; k < 8; ++ k)
                comp_num[k] = 0;   
            if (m(a, b) == 0)
            {
                if ((a - 1 >= 0) && (b - 1 >= 0 ))
                {
                    if ((m(a - 1, b - 1) == 0))
                        ++ black_sum;
                    else
                    { 
                        ++ white_sum;
                        comp_num[white_sum - 1] = m(a, b);
                    }
                }
                if ((b - 1 >= 0 ))
                {
                    if (m(a, b - 1) == 0)
                        ++ black_sum;
                    else
                    {
                        ++ white_sum;
                        comp_num[white_sum - 1] = m(a, b);
                    }
                }
                if ((a - 1 >= 0))
                {
                    if (m(a - 1, b) == 0)
                        ++ black_sum;   
                    else
                    {
                        ++ white_sum;
                        comp_num[white_sum - 1] = m(a, b);
                    }
                }
                if ((a + 1 < mm) && (b + 1 < nn ))
                {
                    if (m(a + 1, b + 1) == 0)
                         ++ black_sum;
                    else
                    {
                        ++ white_sum;
                        comp_num[white_sum - 1] = m(a, b);
                    }
                }
                if ((a + 1 < mm))
                {
                    if (m(a + 1, b) == 0)
                        ++ black_sum;
                    else
                    {
                         ++ white_sum;
                        comp_num[white_sum - 1] = m(a, b);
                    }
                }
                if ((b + 1 < nn))
                {
                    if (m(a, b + 1) == 0)
                        ++ black_sum;
                    else
                    {
                        ++ white_sum;
                        comp_num[white_sum - 1] = m(a, b);
                    }
                }
                if ((a - 1 >= 0 ) && (b + 1 < nn))
                {
                    if (m(a - 1, b + 1) == 0)
                        ++ black_sum;
                    else
                    {
                        ++ white_sum;
                        comp_num[white_sum - 1] = m(a, b);
                    }
                }
                if ((b - 1 >= 0 ) && (a + 1 < mm))
                {
                    if (m(a + 1, b - 1) == 0)
                        ++ black_sum;
                    else
                    {
                        ++ white_sum;
                        comp_num[white_sum - 1] = m(a, b);
                    }
                }
                if (white_sum > black_sum)
                {
                    int flag = 0;
                    for (int p = 0; p < 8; ++ p)
                        for (int q = p + 1; q < 8; ++ q)
                            if (comp_num[p] != comp_num[q]) flag = 1; 
                    if (!flag) m(a, b) = comp_num[0];
                }
            }
        }
    return m;
}

tuple < vector <int>, vector <int>, vector <int>>
    center_square (Matrix <int> m)
{
    vector <int> square(count);
    vector <int> center_x (count);
    vector <int> center_y(count);

    for (int i = 1; i < count; ++i)
    {
        square[i] = 0;
        center_x[i] = 0;
        center_y[i] = 0;
    }
    for (int i = 0; i < mm; ++i)
    {
        for (int j = 0; j < nn; ++j)
        {
            if (m(i, j) != 0)
            {
                ++ square[m(i, j)];
                center_x [m(i, j)] += i;
                center_y [m(i, j)] += j;
 
            }
        }
    }
    int min = square[1], ind = 1, n;
    for (int i = 1; i < count; ++ i)  
    {
        if (square[i] < min) 
        {
            min = square[i];
            ind = i;
        }
        n = square[i];
        center_x[i] = int (round(double(center_x[i]) / n)); 
        center_y[i] = int (round(double(center_y[i]) / n)); 
    }
    ax1_cntr_x = center_x[ind];
    ax1_cntr_y = center_y[ind];
    ax1_numb = ind;
    return make_tuple(square, center_x, center_y);   
}

tuple <vector <int>, vector <int>>
    radius(vector <int> center_x, vector <int> center_y, Matrix<int> m)
{
    vector <int> min_distance_x (count);
    vector <int> min_distance_y (count);
    vector <int> min_rad (count);
    vector <int> max_distance_x (count);
    vector <int> max_distance_y (count);
    vector <int> max_rad (count);

    for (int i = 0; i < count; ++i)
    {
        min_distance_x[i] = mm;
        min_distance_y[i] = nn;
        max_distance_x[i] = center_x[i];
        max_distance_y[i] = center_y[i];
    }
    for (int i = 0; i < mm; ++ i)
        for (int j = 0; j < nn; ++ j)
        {
            if (m(i, j) == 0)
            {
                for (int k = 0; k < count; ++ k)
                {
                    double rad1 = sqrt(sqr(i - center_x[k]) + sqr(j - center_y[k]));
                    double rad2 = sqrt(sqr(min_distance_x[k] - center_x[k]) 
                                        + sqr(min_distance_y[k] - center_y[k]));
                    if (rad1 < rad2) 
                    {
                        min_distance_x[k] = i;
                        min_distance_y[k] = j;
                    }
                }
            }
            if (m(i, j) != 0)
            {
                double rad1 = sqrt(sqr(i - center_x[m(i, j)]) + sqr(j - center_y[m(i, j)]));
                double rad2 = sqrt(sqr(max_distance_x[m(i, j)] - center_x[m(i, j)]) 
                                    + sqr(max_distance_y[m(i, j)] - center_y[m(i, j)]));
                if (rad1 > rad2) 
                {
                    max_distance_x[m(i, j)] = i;
                    max_distance_y[m(i, j)] = j;
                }
            }
        }
    for (int i = 1; i < count; ++ i)
    {
        min_rad[i] = int(round(sqrt(sqr(min_distance_x[i] - center_x[i]) 
                                + sqr(min_distance_y[i] - center_y[i]))));

        max_rad[i] = int(round(sqrt(sqr(max_distance_x[i] - center_x[i]) 
                                + sqr(max_distance_y[i] - center_y[i]))));

    }
    return make_tuple(max_rad, min_rad);
} 
int sqr(int a)
{
    return a * a;
}
tuple <int, Matrix<int>, int, int, int> start_suit_gear(const Matrix<int> &copy_comp, const Image& im)
{
    mm = int(im.n_rows);
    nn = int(im.n_cols);
    Matrix <int> m1 = binar(im);
    m1 = components_search(m1);
    m1 = black_noise(m1);
    vector <int> sq1;
    vector <int> x1;
    vector <int> y1;
    vector <int> minr;
    vector <int> maxr;
    tie(sq1, x1, y1) = center_square(m1);
    tie(maxr, minr) = radius(x1, y1, m1);
    gear_cntr_x = x1[1];
    gear_cntr_y = y1[1];
    int flag;
    Matrix<int> comp1;
    tie(flag, comp1) = suit_gear(copy_comp, m1);
    return make_tuple(flag, comp1, sq1[1],maxr[1], minr[1]);
}

tuple<int, Matrix<int>> suit_gear(const Matrix<int> &copy_comp, const Matrix<int> &m1)
{
    int delta_x = ax_cntr_x - gear_cntr_x;
    int delta_y = ax_cntr_y - gear_cntr_y;
    Matrix<int> ans = copy_comp.deep_copy();
    for (int i = 0; i < mm; ++ i) 
    {
        for (int j = 0; j < nn; ++j)
        {
            if (m1(i, j) != 0)
            {
                if ((i + delta_x >= mm1) || (j + delta_y >= nn1) ||
                    (i + delta_x < 0) || (j + delta_y <0))
                    return make_tuple(0, ans);
                else
                {
                    if ((copy_comp(i + delta_x, j + delta_y) != 0)
                        && (copy_comp(i + delta_x, j + delta_y) != ax_numb))
                        return make_tuple(0, ans);
                    else
                        ans(i + delta_x, j + delta_y) = ax_numb;
                }

            }
        }
    }
    return make_tuple(1, ans);
}

int is_bord(Matrix <int> im, int i, int j)
{
    if (((i + 1 <= mm) && (im(i + 1, j) == 0)) || 
        ((i - 1 >= 0) && (im(i - 1, j) == 0)) ||
        ((j + 1 <= nn) && (im(i, j + 1) == 0)) || 
        ((j - 1 >= 0) && (im(i, j - 1) == 0)))
        return 1;
    else return 0;
}
void distance_rows(Matrix <double> f)
{
    double s;
    vector <double> D(nn);
    vector <int> v(nn);
    vector <double> z(nn);
    int k; 
    for (int ind = 0; ind < mm; ++ ind)
    {
        k = 0;
        v[0] = 0;
        z[0] = - inf;
        z[1] = inf;

        for (int q = 1; q < nn - 1; ++ q)
        {
            do
            {
                s =(f(ind, q) + q * q - (f(ind, v[k]) + v[k] * v[k])) / 
                    (2 * q - 2 * v[k]);
                if (s <= z[k]) -- k;
                else
                {
                    ++ k;
                    v[k] = q;
                    z[k] = s;
                    z[k + 1] = inf;
                    break;
                }
            }
            while (1);
        }
        k = 0;
        for (int q = 0; q < nn; ++ q)
        {
            while (z[k + 1] < q) ++ k;
            D[q] = (q - v[k]) * (q - v[k]) + f(ind, v[k]);
        }
        for (int q = 0; q < nn; ++q)
            f(ind, q) = D[q];
    }
    return;    
} 
void distance_cols(Matrix <double> f)
{
    double s;
    vector <double> D(nn);
    vector <int> v(nn);
    vector <double> z(nn);
    int k; 
    for (int ind = 0; ind < nn; ++ ind)
    {
        k = 0;
        v[0] = 0;
        z[0] = - inf;
        z[1] = inf;

        for (int q = 1; q < mm - 1; ++ q)
        {
            do
            {
                s =(f(q, ind) + q * q - (f(v[k], ind) + v[k] * v[k])) / 
                    (2 * q - 2 * v[k]);
                if (s <= z[k]) -- k;
                else
                {
                    ++ k;
                    v[k] = q;
                    z[k] = s;
                    z[k + 1] = inf;
                    break;
                }
            }
            while (1);
        }
        k = 0;
        for (int q = 0; q < mm; ++ q)
        {
            while (z[k + 1] < q) ++ k;
            D[q] = (q - v[k]) * (q - v[k]) + f(v[k], ind);
        }
        for (int q = 0; q < mm; ++q)
            f(q, ind) = D[q];
    }
    return;    
} 
int cogs(int com, int x, int y, const Matrix <int> &m1, int mxr, int mnr)
{
    int count1 = count;
    int n = 0; 
    Matrix <int> m = m1.deep_copy(); 
    mm = int(m1.n_rows);
    nn = int(m1.n_cols);
    for (int i = 0; i < mm; ++i)
        for (int j = 0; j < nn; ++j)
        {
            if (m(i, j) != com)
                m(i, j) = 0;
            if ((m(i, j) == com)&&(sqrt(sqr(i - x) + sqr(j - y)) <= double(mnr + mxr)/2.0))
                m(i, j) = 0;
            if (m(i, j) != 0) m(i, j) = -1;
        } 
        m = components_search(m);
        n = count - 1;
        count = count1;
        return (n);
}

int main(int argc, char **argv)
{
    if (argc != 4)
    {
        cout << "Usage: " << endl << argv[0]
             << " <in_image.bmp> <out_image.bmp> <out_result.txt>" << endl;
        return 0;
    }

    try 
    {
        Image src_image = load_image(argv[1]);
        string str1, str01, str02, str03, name = argv[1];
        int l = name.length();
        for (int k = 0; k < l - 4; ++ k)
            str1 += name[k];
        str01 = str1 + "_1.bmp";
        str02 = str1 + "_2.bmp";
        str03 = str1 + "_3.bmp";
        Image im1 = load_image(str01.c_str());
        Image im2 = load_image(str02.c_str());
        Image im3 = load_image(str03.c_str());
        ofstream fout(argv[3]); //ôîðìèðóåò òåêñòîâûé ôàéë

        vector<shared_ptr<IObject>> object_array;
        Image dst_image;
        int result_idx;

        tie(result_idx, object_array, dst_image) = repair_mechanism(src_image, im1, im2,im3); //ðàñïàêîâûâàåò tuple â ïåðåìåííûå ïî ïîðÿäêó
        save_image(dst_image, argv[2]);

        fout << result_idx << endl;
        fout << object_array.size() << endl;
        for (const auto &obj : object_array)
            obj->Write(fout);

    } catch (const string &s) {
        cerr << "Error: " << s << endl;
        return 1;
    }
}
