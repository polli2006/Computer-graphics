#include <string>
#include <vector>
#include <fstream>
#include <cassert>
#include <iostream>
#include <cmath>
#include "matrix.h"
#include "classifier.h"
#include "EasyBMP.h"
#include "linear.h"
#include "argvparser.h"

using std::string;
using std::vector;
using std::ifstream;
using std::ofstream;
using std::pair;
using std::make_pair;
using std::cout;
using std::cerr;
using std::endl;

using CommandLineProcessing::ArgvParser;

typedef vector<pair<BMP*, int> > TDataSet;
typedef vector<pair<string, int> > TFileList;
typedef vector<pair<vector<float>, int> > TFeatures;

double pi = 3.14159265;

Matrix <double> Sobel(const Matrix <double> &, int, int, int, int);
vector <double> hog(const Matrix <double> &, const Matrix<double> &, int, int, int, int);
Matrix <double> matr(const Matrix <double> &, int, int, int, int, int, int);
vector <double> color(const TDataSet&, int, int, int, int);
vector <double> pyramid(const Matrix <double> &, int, int, int);
vector <double> start(const Matrix <double> &, int, int, int);
vector <double> kernel (vector <double> &);
// Load list of files and its labels from 'data_file' and
// stores it in 'file_list'
void LoadFileList(const string& data_file, TFileList* file_list) 
{
    ifstream stream(data_file.c_str());

    string filename;
    int label;
    
    int char_idx = data_file.size() - 1;
    for (; char_idx >= 0; --char_idx)
        if (data_file[char_idx] == '/' || data_file[char_idx] == '\\')
            break;
    string data_path = data_file.substr(0,char_idx+1);
    
    while(!stream.eof() && !stream.fail()) {
        stream >> filename >> label;
        if (filename.size())
            file_list->push_back(make_pair(data_path + filename, label));
    }

    stream.close();
}

// Load images by list of files 'file_list' and store them in 'data_set'
void LoadImages(const TFileList& file_list, TDataSet* data_set) {
    for (size_t img_idx = 0; img_idx < file_list.size(); ++img_idx) {
            // Create image
        BMP* image = new BMP();
            // Read image from file
        image->ReadFromFile(file_list[img_idx].first.c_str());
            // Add image and it's label to dataset
        data_set->push_back(make_pair(image, file_list[img_idx].second));
    }
}

// Save result of prediction to file
void SavePredictions(const TFileList& file_list,
                     const TLabels& labels, 
                     const string& prediction_file) {
        // Check that list of files and list of labels has equal size 
    assert(file_list.size() == labels.size());
        // Open 'prediction_file' for writing
    ofstream stream(prediction_file.c_str());

        // Write file names and labels to stream
    for (size_t image_idx = 0; image_idx < file_list.size(); ++image_idx)
        stream << file_list[image_idx].first << " " << labels[image_idx] << endl;
    stream.close();
}

// Exatract features from dataset.
// You should implement this function by yourself =)
void ExtractFeatures(const TDataSet& data_set, TFeatures* features) 
{
    RGBApixel pixel;
    int i, j;
    int cell_dim;
    int m, n;
    for (size_t image_idx = 0; image_idx < data_set.size(); ++ image_idx) 
    {
        n = data_set[image_idx].first->TellHeight();
        m = data_set[image_idx].first->TellWidth();
        Matrix <double> image(m, n);
        //Grayscale
        for (i = 0; i < m; ++ i)
            for (j = 0; j < n; ++ j)
            {
                pixel = data_set[image_idx].first->GetPixel(i, j);
                image(i, j) = 0.299 * pixel.Red + 0.587 * pixel.Green + 0.114 * pixel.Blue;
            }
        cell_dim = 8;
        //BASE PART
        vector <double> histogram  = start(image, m, n, cell_dim);
        //PYRAMID
        cell_dim = 4;
        vector <double> histogram1 = pyramid(image, m, n, cell_dim);

        //hog0+hog1+hog2+hog3+hog4
        histogram.insert(histogram.end(), histogram1.begin(), histogram1.end());
        //colors
        cell_dim = 8;
        vector <double> color_hog = color(data_set, image_idx, m, n, cell_dim);
        //hog+hog_color
        histogram.insert(histogram.end(), color_hog.begin(), color_hog.end());
        //Nonlinear kernel SVM
        vector<double> nonlinear = kernel(histogram);
        histogram.insert(histogram.end(), nonlinear.begin(), nonlinear.end());

        //end
        vector <float> desc;
        for (i = 0; i < int(histogram.size()); ++ i)
        {
            desc.push_back(float(histogram[i]));
        }
        features->push_back(make_pair(desc, data_set[image_idx].second));

    }
}
vector <double> kernel (vector <double> &histogram)
{
    int n = 1;
    vector <double> nonlinear;
    double L = 0.25;
    double lyambda = -n*L;
    double eps = 0.01;
    int i;
    double re1, im1,re2, im2, re0, im0;
    for (i = 0; i < int(histogram.size()); ++ i)
    {
        if ((histogram[i] >= -eps)&&(histogram[i] <= eps))
        {
            re1 = 0;
            im1 = 0;
            re2 = 0;
            im2 = 0;
            re0 = 0;
            im0 = 0;
        }
        else
        {
            if (histogram[i] < 0)
                histogram[i] = -histogram[i];
            //lyambda = 0
            re0 = cos(0 * log(histogram[i])) * sqrt(histogram[i] / cosh(0));
            im0 = sin(0);
            //+lyambda
            re1 = cos(-lyambda * log(histogram[i])) * sqrt(histogram[i] / cosh(pi * lyambda));
            im1 = sin(-lyambda * log(histogram[i])) * sqrt(histogram[i] /cosh(pi * lyambda));

            //-lyambda
            re2 = cos(lyambda * log(histogram[i])) * sqrt(histogram[i] /cosh(-pi * lyambda));
            im2 = sin(lyambda * log(histogram[i])) * sqrt(histogram[i] / cosh(-pi * lyambda));
        }       
        nonlinear.push_back(re0);
        nonlinear.push_back(im0);
        nonlinear.push_back(re1);
        nonlinear.push_back(im1);
        nonlinear.push_back(re2);
        nonlinear.push_back(im2);        
    }
    return nonlinear;   
}
vector <double> pyramid(const Matrix <double> &im, int m, int n, int cell_dim)
{

    int i0, m0, n0, j0;
    //hog1
    i0 = 0;
    j0 = 0;
    m0 = m/2;
    n0 = n/2;
    Matrix <double> part1 = matr(im, i0, j0, m0, n0, m/2, n/2);
    vector <double> hog1 =start(part1, m/2, n/2, cell_dim);
    //hog2
    i0 = 0;
    j0 = n/2;
    m0 = m/2;
    n0 = (n/2)*2;
    Matrix <double> part2 = matr(im, i0, j0, m0, n0, m/2, n/2);
    vector <double> hog2 =start(part2, m/2, n/2, cell_dim);
    //hog3
    i0 = m/2;
    j0 = 0;
    m0 = (m/2)*2;
    n0 = n/2;
    Matrix <double> part3 = matr(im, i0, j0, m0, n0, m/2, n/2);
    vector <double> hog3 =start(part3, m/2, n/2, cell_dim);
    //hog4
    i0 = m/2;
    j0 = n/2;
    m0 = (m/2)*2;
    n0 = (n/2)*2;
    Matrix <double> part4 = matr(im, i0, j0, m0, n0, m/2, n/2);
    vector <double> hog4 =start(part4, m/2, n/2, cell_dim);
    hog1.insert(hog1.end(), hog2.begin(), hog2.end());
    hog1.insert(hog1.end(), hog3.begin(), hog3.end());
    hog1.insert(hog1.end(), hog4.begin(), hog4.end());
    return hog1;
}
vector <double> start(const Matrix <double> &im, int m, int n, int cell_dim)
{
    int i, j, m1, n1;
    int x1 = -1;
    int x3 = 1;
    //horizont
    Matrix <double> horiz = Sobel(im, m, n - 2, x1, x3);
    //verical
    x1 = 1;
    x3 = -1;
    Matrix <double> vertic = Sobel(im, m - 2, n, x1, x3);
    //gradient
    m1 = m - 2;
    n1 = n - 2;
    Matrix <double> grad_abs(m1, n1);
    Matrix <double> grad_direct(m1, n1);
    for (i = 0; i < m1; ++ i)
        for (j = 0; j < n1; ++ j)
        {
            grad_abs(i, j) = sqrt(horiz(i + 1, j) * horiz(i + 1, j) 
                + vertic(i, j + 1) * vertic(i, j + 1));
            grad_direct(i, j) = atan2(-vertic(i, j + 1), horiz(i + 1, j));
        }
    //histogram
    int sector_numb = 8;
    vector <double> histogram  = hog(grad_abs, grad_direct, cell_dim, sector_numb, m1, n1); 
    return histogram;    
}
Matrix <double> Sobel(const Matrix <double>& image, int a, int b, int x1, int x3)
{
    Matrix <double> filter (a, b);
    int i0, j0;
    if (x1 < 0)
    {
        i0 = 0;
        j0 = 1;
    }
    else
    {
        i0 = 1;
        j0 = 0;
    }
    for (int i = i0; i < a + i0; ++ i)
        for (int j = j0; j < b + j0; ++ j)
        {

            filter(i - i0, j - j0) = image(i - i0, j - j0) * x1 + image(i + i0, j + j0) * x3;
        }
    return filter;
}
//create new matrix from old
Matrix <double> matr(const Matrix <double> & mtrx, int i0, int j0, int m0, int n0, int m, int n)
{
    Matrix <double> temp(m, n);
    for (int i = i0; i < m0; ++ i)
        for (int j = j0; j< n0; ++ j)
            temp(i - i0, j - j0) = mtrx(i, j);

    return temp;
}
vector <double> hog(const Matrix <double> &grad_abs, const Matrix <double> & grad_direct, 
    int cell_dim, int sector_numb, int m1, int n1)
{
    int cell_cols, cell_rows, numb_of_sect, numb_of_cell, a, b;
    int size1, size2, i, j;
    double a1, sum;
    int cell_numb = cell_dim * cell_dim;
    vector <double> histogram (cell_dim * cell_dim * sector_numb, 0);
    a1 = 2 * pi / sector_numb;
    cell_cols = (n1) / cell_dim; //size in width of one cell
    cell_rows = (m1) / cell_dim; //size in height
    size1 = cell_rows * cell_dim; 
    size2 = cell_cols * cell_dim;
    for (i = 0; i < size1; ++ i)
        for (j = 0; j < size2; ++j)
        {
            a = i / cell_rows;
            b = j / cell_cols;
            numb_of_cell = a * cell_dim + b;
            numb_of_sect = (grad_direct(i, j) + pi) / a1;
            if (numb_of_sect < 0) numb_of_sect = 0;
            if (numb_of_sect >= sector_numb) numb_of_sect = sector_numb - 1;
            histogram[numb_of_cell * sector_numb + numb_of_sect] += grad_abs(i, j);   
        }
    //norma
    a = 0;
    sum = 0;
    vector <double> mas_sum(cell_numb);
    for (i = 0; i < cell_numb * sector_numb; ++ i)
    {
        ++a;
        sum += histogram[i] * histogram[i];
        if (a == sector_numb)
        {
            if (sum > 0)
                mas_sum[i / sector_numb] = sqrt(sum);
            else mas_sum[i / sector_numb] = 1;
            a = 0;
            sum = 0;
        }
    }
    for (i = 0; i < cell_numb * sector_numb; ++ i)
    {

        histogram[i] /=(mas_sum[i / sector_numb]);
    }

    return histogram;
}
vector <double> color(const TDataSet& im, int image_idx, int m1, int n1, int cell_dim)
{
    int size1, size2, i, j;
    int cell_cols, cell_rows, numb_of_cell, a, b;
    vector <double> red (cell_dim * cell_dim, 0);
    vector <double> green (cell_dim * cell_dim, 0);
    vector <double> blue (cell_dim * cell_dim, 0);
    cell_cols = (n1) / cell_dim; //size in width of one cell
    cell_rows = (m1) / cell_dim; //size in height
    size1 = cell_rows * cell_dim; 
    size2 = cell_cols * cell_dim;
    int pix_in_block = cell_cols * cell_rows;
    RGBApixel pixel;
    for (i = 0; i < size1; ++ i)
        for (j = 0; j < size2; ++j)
        {
            a = i / cell_rows;
            b = j / cell_cols;
            numb_of_cell = a * cell_dim + b;

            pixel = im[image_idx].first->GetPixel(i, j);
            red[numb_of_cell] += pixel.Red;
            green[numb_of_cell] += pixel.Green;
            blue[numb_of_cell] += pixel.Blue;   
        }
    //norma
    for (i = 0; i < cell_dim * cell_dim; ++ i)
    {
        if (pix_in_block != 0)
        {
            red[i] /= (pix_in_block * 255);
            green[i] /= (pix_in_block * 255);
            blue[i] /= (pix_in_block * 255);
        }
    }
    red.insert(red.end(), green.begin(), green.end());
    red.insert(red.end(), blue.begin(), blue.end());
    return red;

}
// Clear dataset structure
void ClearDataset(TDataSet* data_set) {
        // Delete all images from dataset
    for (size_t image_idx = 0; image_idx < data_set->size(); ++image_idx)
        delete (*data_set)[image_idx].first;
        // Clear dataset
    data_set->clear();
}

// Train SVM classifier using data from 'data_file' and save trained model
// to 'model_file'
void TrainClassifier(const string& data_file, const string& model_file) {
        // List of image file names and its labels
    TFileList file_list;
        // Structure of images and its labels
    TDataSet data_set;
        // Structure of features of images and its labels
    TFeatures features;
        // Model which would be trained
    TModel model;
        // Parameters of classifier
    TClassifierParams params;
    
        // Load list of image file names and its labels
    LoadFileList(data_file, &file_list);
        // Load images
    LoadImages(file_list, &data_set);
        // Extract features from images
    ExtractFeatures(data_set, &features);
        // PLACE YOUR CODE HERE
        // You can change parameters of classifier here
    params.C = 0.01;
    TClassifier classifier(params);
        // Train classifier
    classifier.Train(features, &model);
        // Save model to file
    model.Save(model_file);
        // Clear dataset structure
    ClearDataset(&data_set);
}

// Predict data from 'data_file' using model from 'model_file' and
// save predictions to 'prediction_file'
void PredictData(const string& data_file,
                 const string& model_file,
                 const string& prediction_file) {
        // List of image file names and its labels
    TFileList file_list;
        // Structure of images and its labels
    TDataSet data_set;
        // Structure of features of images and its labels
    TFeatures features;
        // List of image labels
    TLabels labels;

        // Load list of image file names and its labels
    LoadFileList(data_file, &file_list);
        // Load images
    LoadImages(file_list, &data_set);
        // Extract features from images
    ExtractFeatures(data_set, &features);

        // Classifier 
    TClassifier classifier = TClassifier(TClassifierParams());
        // Trained model
    TModel model;
        // Load model from file
    model.Load(model_file);
        // Predict images by its features using 'model' and store predictions
        // to 'labels'
    classifier.Predict(features, model, &labels);

        // Save predictions
    SavePredictions(file_list, labels, prediction_file);
        // Clear dataset structure
    ClearDataset(&data_set);
}

int main(int argc, char** argv) {
    // Command line options parser
    ArgvParser cmd;
        // Description of program
    cmd.setIntroductoryDescription("Machine graphics course, task 2. CMC MSU, 2014.");
        // Add help option
    cmd.setHelpOption("h", "help", "Print this help message");
        // Add other options
    cmd.defineOption("data_set", "File with dataset",
        ArgvParser::OptionRequiresValue | ArgvParser::OptionRequired);
    cmd.defineOption("model", "Path to file to save or load model",
        ArgvParser::OptionRequiresValue | ArgvParser::OptionRequired);
    cmd.defineOption("predicted_labels", "Path to file to save prediction results",
        ArgvParser::OptionRequiresValue);
    cmd.defineOption("train", "Train classifier");
    cmd.defineOption("predict", "Predict dataset");
        
        // Add options aliases
    cmd.defineOptionAlternative("data_set", "d");
    cmd.defineOptionAlternative("model", "m");
    cmd.defineOptionAlternative("predicted_labels", "l");
    cmd.defineOptionAlternative("train", "t");
    cmd.defineOptionAlternative("predict", "p");

        // Parse options
    int result = cmd.parse(argc, argv);

        // Check for errors or help option
    if (result) {
        cout << cmd.parseErrorDescription(result) << endl;
        return result;
    }

        // Get values 
    string data_file = cmd.optionValue("data_set");
    string model_file = cmd.optionValue("model");
    bool train = cmd.foundOption("train");
    bool predict = cmd.foundOption("predict");

        // If we need to train classifier
    if (train)
        TrainClassifier(data_file, model_file);
        // If we need to predict data
    if (predict) {
            // You must declare file to save images
        if (!cmd.foundOption("predicted_labels")) {
            cerr << "Error! Option --predicted_labels not found!" << endl;
            return 1;
        }
            // File to save predictions
        string prediction_file = cmd.optionValue("predicted_labels");
            // Predict data
        PredictData(data_file, model_file, prediction_file);
    }
}