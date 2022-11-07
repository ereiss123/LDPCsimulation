/* errtopng
   Program to analyze error patterns and convert them to
   PNG images.
*/


#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <cstdio>
#include <malloc.h>
#include <png.h>

using namespace std;

// This takes the float value 'val', converts it to red, green & blue values, then 
// sets those values into the image memory buffer location pointed to by 'ptr'
inline void setRGB(png_byte *ptr, float val);

// This function actually writes out the PNG image file. The string 'title' is
// also written into the image file

int writeImage(const char * filename, int width, int height, vector<vector<float> > & buffer, char* title);

void printMatrix(vector<vector<float> > & matrix);
int readMatrix(char * fname, vector<vector<float> > & matrix);
void shiftMatrix(vector<vector<float> > & matrix, float shiftval, float scaleVal);
void mergeMatrices(vector<vector<float> > & matrix1, vector<vector<float> > & matrix2);
void countErrorTrace(vector<vector<float> > & matrix ,vector<float> & errTrace);
void fprintMatrix(const char * fname,vector<vector<float> > & matrix);


int main(int argc, char *argv[])  
{
  // Make sure that the output filename argument has been provided
  if (argc < 3) {
    fprintf(stderr, "Usage: %s outfile infile [infile2 infile3 ...]", argv[0]);
    return 1;
  }

  string fprefix(argv[1]);
  stringstream pngfnamestr;
  pngfnamestr << fprefix << ".png";
  string pngfname(pngfnamestr.str());

  stringstream errfnamestr;
  errfnamestr << fprefix << ".err";
  string errfname(errfnamestr.str());

  vector<vector<float> > matrix;
  vector<vector<float> > errorHistory;

  //cout << "Attempting to read from " << argv[1] << endl;

  for (int idx=2;idx<argc;idx++)
    {
      vector<vector<float> > newMatrix;
      vector<float> errTrace;
      if (readMatrix(argv[idx],newMatrix)==1)
	{
	  cout << "Some error happened.\n";
	  return 1;
	}

      cout << "Got matrix with " << newMatrix.size() << " rows.\n";
      shiftMatrix(newMatrix,-1,-1);

      countErrorTrace(newMatrix,errTrace);
      errorHistory.push_back(errTrace);

      mergeMatrices(matrix,newMatrix);
    }
  int rows= matrix.size();
  int columns = matrix[0].size();
 
  cout << "Writing png file to " << pngfname << endl;
 
  int result = writeImage(pngfname.c_str(), columns, rows, matrix, "This is my test image");

  cout << "Writing error trace to " << errfname << endl;
  fprintMatrix(errfname.c_str(),errorHistory);
  return 0;
}

void mergeMatrices(vector<vector<float> > & matrix1, vector<vector<float> > & matrix2)
{
  int rows1 = matrix1.size();
  int rows2 = matrix2.size();
  int minsize=min(rows1,rows2);
  for (int idx=0; idx<minsize; idx++)
    for (int jdx=0; jdx<matrix1[idx].size(); jdx++)
      matrix1[idx][jdx] += matrix2[idx][jdx];
  if (rows1<rows2)
    for (int idx=rows1; idx<rows2; idx++)
      matrix1.push_back(matrix2[idx]);
}


int readMatrix(char * fname, vector<vector<float> > & matrix)
{
  ifstream infile(fname,ios::in);
  if (!infile)
    {
      cout << "Failed to open input file.\n";
      return 1;
    }

  //if (infile.is_open())
  //  cout << "It's open...\n";

  //infile.seekg(0);
  string line;
  while (getline(infile,line)) //!infile.eof())
    {
      vector<float> entries;

      /*
      string line;
      if (getline(infile,line))
	{
	  cout << "Acquired line: ";
	  cout << line;
	  cout << "\nRecovered values:\n";
      */
      stringstream ss(line);
      for(float num; ss >> num;)
	{
	  //cout << num << " ";
	  entries.push_back(num);
	}
      //cout << endl;
      /*}
      else
	{
	  cout << "Failed to read line from input file.\n";
	  return 1;
	  }*/
      if (entries.size() > 0)
	matrix.push_back(entries);
    }
  //cout << "Got " << matrix.size() << " rows.\n";
  infile.close();
  return 0;
}

void printMatrix(vector<vector<float> > matrix)
{
  for (int idx=0; idx<matrix.size(); idx++)
    {
      for (int jdx=0; jdx<matrix[idx].size(); jdx++)
	cout << matrix[idx][jdx] << " ";
      cout << endl;
    }
}


void shiftMatrix(vector<vector<float> > & matrix, float shiftval, float scaleVal)
{
  for (int idx=0; idx<matrix.size(); idx++)
    for (int jdx=0; jdx<matrix[idx].size(); jdx++)
      matrix[idx][jdx] = (matrix[idx][jdx] + shiftval)*scaleVal;
}


inline void setRGB(png_byte *ptr, float val)
{
  int v = (int)(val * 3);
  if (v < 0) v = 0;
  if (v > 767) v = 767;
  int offset = v % 256;

  if (v<256) {
    ptr[0] = 0; ptr[1] = 0; ptr[2] = offset;
  }
  else if (v<512) {
    ptr[0] = 0; ptr[1] = offset; ptr[2] = 255-offset;
  }
  else {
    ptr[0] = offset; ptr[1] = 255-offset; ptr[2] = 0;
  }
}

int writeImage(const char * filename, int width, int height, vector<vector<float> > & buffer, char* title)
{
  int code = 0;
  FILE *fp;
  png_structp png_ptr;
  png_infop info_ptr;
  png_bytep row;
	
  // Open file for writing (binary mode)
  fp = fopen(filename, "wb");
  if (fp == NULL) {
    fprintf(stderr, "Could not open file %s for writing\n", filename);
    code = 1;
    goto finalise;
  }

  // Initialize write structure
  png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  if (png_ptr == NULL) {
    fprintf(stderr, "Could not allocate write struct\n");
    code = 1;
    goto finalise;
  }

  // Initialize info structure
  info_ptr = png_create_info_struct(png_ptr);
  if (info_ptr == NULL) {
    fprintf(stderr, "Could not allocate info struct\n");
    code = 1;
    goto finalise;
  }

  // Setup Exception handling
  if (setjmp(png_jmpbuf(png_ptr))) {
    fprintf(stderr, "Error during png creation\n");
    code = 1;
    goto finalise;
  }

  png_init_io(png_ptr, fp);

  // Write header (8 bit colour depth)
  png_set_IHDR(png_ptr, info_ptr, width, height,
	       8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
	       PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

  // Set title
  if (title != NULL) {
    png_text title_text;
    title_text.compression = PNG_TEXT_COMPRESSION_NONE;
    title_text.key = "Title";
    title_text.text = title;
    png_set_text(png_ptr, info_ptr, &title_text, 1);
  }

  png_write_info(png_ptr, info_ptr);

  // Allocate memory for one row (3 bytes per pixel - RGB)
  row = (png_bytep) malloc(3 * width * sizeof(png_byte));

  // Write image data
  int x, y;
  for (y=0 ; y<height ; y++) {
    for (x=0 ; x<width ; x++) {
      setRGB(&(row[x*3]), buffer[y][x]);
    }
    png_write_row(png_ptr, row);
  }

  // End write
  png_write_end(png_ptr, NULL);

 finalise:
  if (fp != NULL) fclose(fp);
  if (info_ptr != NULL) png_free_data(png_ptr, info_ptr, PNG_FREE_ALL, -1);
  if (png_ptr != NULL) png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
  if (row != NULL) free(row);

  return code;
}


void countErrorTrace(vector<vector<float> > & matrix ,vector<float> & errTrace)
{
  for (int idx=0; idx<matrix.size(); idx++)
    {
      int err=0;
      for (int jdx=0; jdx<matrix[idx].size(); jdx++)
	err += matrix[idx][jdx];
      errTrace.push_back(err);
    }
}


void fprintMatrix(const char * fname,vector<vector<float> > & matrix)
{
  ofstream of(fname,ios::out);
  for (int idx=0; idx<matrix.size(); idx++)
    {
      for (int jdx=0; jdx<matrix[idx].size(); jdx++)
	of << matrix[idx][jdx] << "\t";
      of << endl;
    }  
  of.close();
}
