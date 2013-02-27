// CSCI 4810/6810 Program 4 (Assignment 3)
// Yibin Liao
// tigerlyb@uga.edu

#include <GL/glut.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <fstream>

#define PI 3.14159265
#define LEFT 1
#define RIGHT 2
#define BELOW 4
#define ABOVE 8

using namespace std;

string inputName;  	// define the datalines input file name
int lineNum;			// define the number of lines to be read from input.txt
int dashLength;			// define the dash length entered by user
double vx, vy, vz;		// define the viewpoint coordinate
double s;				// define the screen size
double d;				// define the distance of the screen
double screenPixel;		// define the coordinate system of the screen runs from 0 to the number of screenPixel
double angleRotate;		// define the angle of rotation
bool orginalLines;		// define a flag to determine if the line is drawing before or after transformation
double rx1, rx2, ry1, ry2, rz1, rz2;
// define and initial the translate matrix
double	matrixTranslate [4][4] = {
	{1, 0, 0, 0},
	{0, 1, 0, 0}, 
	{0, 0, 1, 0},
	{0, 0, 0, 1}
};

// define the initial scale matrix
double	matrixScale [4][4] = {
	{1, 0, 0, 0},
	{0, 1, 0, 0}, 
	{0, 0, 1, 0},
	{0, 0, 0, 1}
};

// define the initial rotate matrix
double matrixRotate [4][4] = {
	{1, 0, 0, 0},
	{0, 1, 0, 0}, 
	{0, 0, 1, 0},
	{0, 0, 0, 1}
};

// define the initial concatenate matrix
double matrixConcatenate [4][4] = {
	{1, 0, 0, 0},
	{0, 1, 0, 0}, 
	{0, 0, 1, 0},
	{0, 0, 0, 1}
};

// initial the matrix
double arbitraryRotation [4][4] = {
	{1, 0, 0, 0},
	{0, 1, 0, 0}, 
	{0, 0, 1, 0},
	{0, 0, 0, 1}
};

// initial the matrix
double V [4][4] = {
	{1, 0, 0, 0},
	{0, 1, 0, 0}, 
	{0, 0, 1, 0},
	{0, 0, 0, 1}
};

// initial the matrix
double N [4][4] = {
	{1, 0, 0, 0},
	{0, 1, 0, 0}, 
	{0, 0, 1, 0},
	{0, 0, 0, 1}
};

// initial the matrix
double T2 [4][4] = {
	{1, 0, 0, 0},
	{0, 0, -1, 0}, 
	{0, 1, 0, 0},
	{0, 0, 0, 1}
};

// initial the matrix
double T3 [4][4] = {
	{1, 0, 0, 0},
	{0, 1, 0, 0}, 
	{0, 0, 1, 0},
	{0, 0, 0, 1}
};

// initial the matrix
double T4 [4][4] = {
	{1, 0, 0, 0},
	{0, 1, 0, 0}, 
	{0, 0, 1, 0},
	{0, 0, 0, 1}
};

double dataW [100][6];					// define a dataline array to store the line points at the WCS
double dataC [100][6];					// define a dataline array to store the line points at the CCS
int dataS [100][4];					// define a dataline array to stroe the line points in the screen

void init (double x, double y)
{
	glClearColor (0.0, 0.0, 0.0, 0.0); // set display-window color to black
	glMatrixMode (GL_PROJECTION); // set projection parameters
	gluOrtho2D (0.0, x, 0.0, y); // set the x-coordinate from 0 to 399 and y-coordinate from 0 to 399
}

void setPixel (int x, int y)
{
	glColor3f (1.0, 1.0, 1.0); // activate the pixel by setting the point color to white, draw orginal lines
	glBegin (GL_POINTS);
		glVertex2i (x, y);	// set the point
	glEnd ();
	glFlush (); // process all openGL routines as quickly as possible
}

void setPixelRed (int x, int y)
{
	glColor3f (1.0, 0.0, 0.0); // activate the pixel by setting the point color to red, drwa lines after transformation
	glBegin (GL_POINTS);
		glVertex2i (x, y);	// set the point
	glEnd ();
	glFlush (); // process all openGL routines as quickly as possible
}

void offPixel (int x, int y)
{
	glColor3f (0.0, 0.0, 0.0); // deactivate the pixel by setting the point color to black
	glBegin (GL_POINTS);
		glVertex2i (x, y);	// set the point
	glEnd ();
	glFlush (); // process all openGL routines as quickly as possible
}

int round (double r)
{
	return int (r + 0.5); // round the value to the nearest integer
}

// bresenham line drawing algorithm
void bresenhamAlg (bool orginal, int dashLength, int x0, int y0, int x1, int y1)
{
	int dx = abs (x1 - x0);
	int dy = abs (y1 - y0);
	int x, y;
	
	// draw the solid line
	if (dashLength == 0) {
		if (dy == 0) { // draw the horizontal line
			if (x0 < x1) { // set the left point to the starting point
				x = x0;
				y = y0;
			}
			else { // exchange the right point to the starting point and left point to the end point
				x = x1;			
				y = y1;
				x1 = x0;
				y1 = y0;
			}

			if (orginal)
				setPixel (x, y); // activate the starting point
			else
				setPixelRed (x, y); // activate the starting point

			while (x < x1) {
				x++;
				if (orginal)
					setPixel (x, y); // activate the starting point
				else
					setPixelRed (x, y); // activate the starting point
			}
		}
		else if (dx == 0) {
			if (y0 < y1) { // set the left point to the starting point
				x = x0;
				y = y0;
			}
			else { // exchange the right point to the starting point and left point to the end point
				x = x1;			
				y = y1;
				y1 = y0;
				x1 = x0;
			}
			if (orginal)
				setPixel (x, y); // activate the starting point
			else
				setPixelRed (x, y); // activate the starting point
			while (y < y1) {
				y++;
				if (orginal)
					setPixel (x, y); // activate the starting point
				else
					setPixelRed (x, y); // activate the starting point
			}
		}
		else if (dx >= dy){ //slope |m| <= 1.0;
			int e = (dy << 1) - dx;
			int inc1 = dy << 1;
			int inc2 = (dy - dx) << 1;	
			if (x0 < x1) { // set the left point to the starting point
				x = x0;
				y = y0;
			}
			else { // exchange the right point to the starting point and left point to the end point
				x = x1;			
				y = y1;
				x1 = x0;
				y1 = y0;
			}

			if (orginal)
				setPixel (x, y); // activate the starting point
			else
				setPixelRed (x, y); // activate the starting point
			
			while (x < x1) {
				if (e < 0)
					e += inc1;
				else {
					if (y < y1) {
						y++;
						e += inc2;
					} else {
						y--;
						e += inc2;
					}
				}
				x++;
				if (orginal)
					setPixel (x, y); // activate the starting point
				else
					setPixelRed (x, y); // activate the starting point
			}
		} 
		else { //slope |m| > 1.0;
			int e = (dx << 1) - dy;
			int inc1 = dx << 1;
			int inc2 = (dx - dy) << 1;
			if (y0 < y1) { // set the left point to the starting point
				x = x0;
				y = y0;
			}
			else { // exchange the right point to the starting point and left point to the end point
				x = x1;			
				y = y1;
				y1 = y0;
				x1 = x0;
			}
			if (orginal)
				setPixel (x, y); // activate the starting point
			else
				setPixelRed (x, y); // activate the starting point
			while (y < y1) {
				if (e < 0)
					e += inc1;
				else {
					if (x > x1){
						x--;
						e += inc2;
					} else {
						x++;
						e += inc2;
					}
				}
				y++;
				if (orginal)
					setPixel (x, y); // activate the starting point
				else
					setPixelRed (x, y); // activate the starting point
			}
		}
	}
	// draw the dashed line
	else { 
		int dashed = 0; // define a dashed flag to determined if set or turn off the pixel
		int dot = 0; // used to count the activated pixels;

		if (dx >= dy){ //slope |m| <= 1.0;
			int e = (dy << 1) - dx;
			int inc1 = dy << 1;
			int inc2 = (dy - dx) << 1;
			// set the left point to the starting point
			if (x0 < x1) {
				x = x0;
				y = y0;
			}
			else { // exchange the right point to the starting point and left point to the end point
				x = x1;			
				y = y1;
				x1 = x0;
				y1 = y0;
			}
			if (orginal)
				setPixel (x, y); // activate the starting point
			else
				setPixelRed (x, y); // activate the starting point
			dot++;
			while (x < x1) {
				if (e < 0)
					e += inc1;
				else {
					if (y < y1) {
						y++;
						e += inc2;
					} else {
						y--;
						e += inc2;
					}
				}
				x++;
				// if the dashed flag is on, deactivate the number of pixels determined by the dash length
				if (dashed) { 
					offPixel (x, y);
					dot++;
					// if dot == dashLength, deactivated pixels are done, turn off the dashed flag
					if (dot%dashLength == 0) 
					{
						dashed = 0;
					}
				} else {
					if (orginal)
						setPixel (x, y); // activate the starting point
					else
						setPixelRed (x, y); // activate the starting point
					dot++;
					// if dot == dashLength, deactivated pixels are done, turn off the dashed flag
					if (dot%dashLength == 0) 
					{
						dashed = 1;
					}
				}		
			}
		} else { // slope |m| > 1.0
			int e = (dx << 1) - dy;
			int inc1 = dx << 1;
			int inc2 = (dx - dy) << 1;

			if (y0 < y1) { // set the left point to the starting point
				x = x0;
				y = y0;
			}
			else { // exchange the right point to the starting point and left point to the end point
				x = x1;			
				y = y1;
				y1 = y0;
				x1 = x0;
			}
			
			if (orginal)
				setPixel (x, y); // activate the starting point
			else
				setPixelRed (x, y); // activate the starting point
			dot++;
			while (y < y1) {
				if (e < 0)
					e += inc1;
				else {
					if (x > x1){
						x--;
						e += inc2;
					} else {
						x++;
						e += inc2;
					}
				}
				y++;
				if (dashed) { // if the dashed flag is on, deactivate the number of pixels determined by the dash length
					offPixel (x, y);
					dot++;
					// if dot == dashLength, deactivated pixels are done, turn off the dashed flag
					if (dot%dashLength == 0) 
					{
						dashed = 0;
					}
				} else {
					if (orginal)
						setPixel (x, y); // activate the starting point
					else
						setPixelRed (x, y); // activate the starting point
					dot++;
					// if dot == dashLength, deactivated pixels are done, turn off the dashed flag
					if (dot%dashLength == 0) 
					{
						dashed = 1;
					}
				}				
			} 
		}
	}
}

// basic translate function to move from one location to another
void basicTranslate (double Tx, double Ty, double Tz)
{
	matrixTranslate [3][0] = -Tx;	// assign the x axis translate displacement 
	matrixTranslate [3][1] = -Ty;	// assign the y axis translate displacement
	matrixTranslate [3][2] = -Tz;	// assign the z axis translate displacement

	// print out the whole translate matrix
	cout << "Translation Matrix" << endl;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) 
			cout << matrixTranslate [i][j] << ", ";
		cout << endl;
	}
	cout << endl;
}

// basic scale function to enlarger or  shrinker an object
void basicScale (double Sx, double Sy, double Sz)
{
	matrixScale [0][0] = Sx;	// assign the x scaling factor
	matrixScale [1][1] = Sy;	// assign the y scaling factor
	matrixScale [2][2] = Sz;	// assign the z scaling factor

	// print out the whole scaling matrix
	cout << "Scaling Matrix" << endl;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++)
			cout << matrixScale [i][j] << ", ";
		cout << endl;
	}
	cout << endl;
}

// basic rotation function about x axis
void basicRotateX (double angle)
{
	// assign the angle of rotation to the matrix
	matrixRotate [1][1] = cos (angle*PI/180);
	matrixRotate [2][1] = -sin (angle*PI/180);
	matrixRotate [1][2] = sin (angle*PI/180);
	matrixRotate [2][2] = cos (angle*PI/180);
	
	// print out the whole rotation matrix
	cout << "Rotation Matrix by X axis" << endl;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++)
			cout << matrixRotate [i][j] << ", ";
		cout << endl;
	}
	cout << endl;
}


// basic rotation function about y axis
void basicRotateY (double angle)
{
	// assign the angle of rotation to the matrix
	matrixRotate [0][0] = cos (angle*PI/180);
	matrixRotate [0][2] = -sin (angle*PI/180);
	matrixRotate [2][0] = sin (angle*PI/180);
	matrixRotate [2][2] = cos (angle*PI/180);
	
	// print out the whole rotation matrix
	cout << "Rotation Matrix by Y axis" << endl;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++)
			cout << matrixRotate [i][j] << ", ";
		cout << endl;
	}
	cout << endl;
}

// basic rotation function about z axis
void basicRotateZ (double angle)
{
	// assign the angle of rotation to the matrix
	matrixRotate [0][0] = cos (angle*PI/180);
	matrixRotate [1][0] = -sin (angle*PI/180);
	matrixRotate [0][1] = sin (angle*PI/180);
	matrixRotate [1][1] = cos (angle*PI/180);
	
	// print out the whole rotation matrix
	cout << "Rotation Matrix by Z axis" << endl;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++)
			cout << matrixRotate [i][j] << ", ";
		cout << endl;
	}
	cout << endl;
}

// combine two int transformation matrix into one, the result will be stored in matrix1
void concatenate (double matrix1 [4][4], double matrix2 [4][4])
{
	double matrix [4][4];

	// perform the calculation of matrix1 * matrix2
	for (int row = 0; row < 4; row++)
	{
		for (int column = 0; column < 4; column++)
		{
			matrix [row][column] = 0;
			for (int i = 0; i < 4; i++)
			{
				matrix [row][column] = matrix [row][column] + matrix1 [row][i] * matrix2 [i][column];				
			}			
		}
	}

	// print out the cancontenate matrix and assign to the matrixConcatenate applied by later transformation
	cout << "Concatenate Matrix for this step: " << endl;
	for (int row = 0; row < 4; row++)
	{
		for (int column = 0; column < 4; column++)
		{
			matrixConcatenate [row][column] = matrix [row][column];
			cout << matrixConcatenate [row][column] << ", ";
		}
		cout << endl;
	}
	cout << endl;
}

void matrixT3 (double x, double y)
{
	double cosC = y / sqrt((pow(x, 2) + pow(y, 2)));
	double sinS = x / sqrt((pow(x, 2) + pow(y, 2)));
	// assign the angle of rotation to the matrix
	T3 [0][0] = -cosC;
	T3 [0][2] = sinS;
	T3 [2][0] = -sinS;
	T3 [2][2] = -cosC;
	
	// print out the whole rotation matrix
	cout << "Matrix T3" << endl;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++)
			cout << T3 [i][j] << ", ";
		cout << endl;
	}
	cout << endl;
}

void matrixT4 (double x, double y)
{
	double cosC = y / sqrt((pow(x, 2) + pow(y, 2)));
	double sinS = x / sqrt((pow(x, 2) + pow(y, 2)));
	// assign the angle of rotation to the matrix
	T4 [1][1] = cosC;
	T4 [2][1] = -sinS;
	T4 [1][2] = sinS ;
	T4 [2][2] = cosC;
	
	// print out the whole rotation matrix
	cout << "Matrix T4" << endl;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++)
			cout << T4 [i][j] << ", ";
		cout << endl;
	}
	cout << endl;
}

void matrixV (double m[4][4])
{
	// print out the cancontenate matrix and assign to the matrixConcatenate applied by later transformation
	cout << "Matrix V: " << endl;
	for (int row = 0; row < 4; row++)
	{
		for (int column = 0; column < 4; column++)
		{
			V [row][column] = m [row][column];
			cout << V [row][column] << ", ";
		}
		cout << endl;
	}
	cout << endl;
}

void matrixN ()
{
	N [0][0] = 2*d / s;
	N [1][1] = 2*d / s;
	// print out the whole rotation matrix
	cout << "Matrix N" << endl;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++)
			cout << N [i][j] << ", ";
		cout << endl;
	}
	cout << endl;
}

// apply the int transformation matrix to the lines
void applyTransformation (double matrix [4][4], double datalines [100][6])
{
	double pointData [6];
	cout << "Transformation dataline result: " << endl;

	// perform the calculation of line point * matrix
	for (int i = 0; i < lineNum; i++) {
		pointData [0] = datalines [i][0] * matrix[0][0] + datalines[i][1] * matrix [1][0] + datalines [i][2] * matrix [2][0] + matrix [3][0] ; 
		pointData [1] = datalines [i][0] * matrix[0][1] + datalines[i][1] * matrix [1][1] + datalines [i][2] * matrix [2][1] + matrix [3][1] ; 
		pointData [2] = datalines [i][0] * matrix[0][2] + datalines[i][1] * matrix [1][2] + datalines [i][2] * matrix [2][2] + matrix [3][2] ; 
		
		pointData [3] = datalines [i][3] * matrix[0][0] + datalines[i][4] * matrix [1][0] + datalines [i][5] * matrix [2][0] + matrix [3][0] ; 
		pointData [4] = datalines [i][3] * matrix[0][1] + datalines[i][4] * matrix [1][1] + datalines [i][5] * matrix [2][1] + matrix [3][1] ; 
		pointData [5] = datalines [i][3] * matrix[0][2] + datalines[i][4] * matrix [1][2] + datalines [i][5] * matrix [2][2] + matrix [3][2] ; 
				
		// print out the result after transformation and store them in a new datalines array 
		for (int j = 0; j < 6; j++) {
			dataC [i][j] = pointData [j];
			cout << dataC [i][j] << ", ";
		}
		cout << endl;
	}
}


// read datalines from an external file
void inputLines (string datalines)
{
	ifstream inFile;
	inFile.open (datalines);	// open the input file 
	
	cout << endl;
	cout << "Original datalines: " << endl;
	for (int i = 0; i < lineNum; i++)
	{
		for (int j = 0; j < 6; j++)
		{
			inFile >> dataW [i][j]; // read the points in the file, and store them in a data array
			cout << dataW [i][j] << ", "; // print out the data
		}
		cout << endl;
	}
	cout << endl;
	inFile.close ();
}

void output ()
{
	double vcx, vcy, vsx, vsy;
	vcx = vcy = vsx = vsy = screenPixel / 2;
	cout << "Output: " << endl;
	for (int i = 0; i < lineNum; i++)
	{
		dataS [i][0] = round((dataC [i][0] / dataC [i][2]) * vsx + vcx);
		dataS [i][1] = round((dataC [i][1] / dataC [i][2]) * vsy + vcy);
		dataS [i][2] = round((dataC [i][3] / dataC [i][5]) * vsx + vcx);
		dataS [i][3] = round((dataC [i][4] / dataC [i][5]) * vsy + vcy);
	}
	for (int i = 0; i < lineNum; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cout << dataS [i][j] << ", "; // print out the data
		}
		cout << endl;
	}
	cout << endl;

}

// rotating a object about an arbitrary axis by theta degree
void Rotation (double angle, double x1, double y1, double z1, double x2, double y2, double z2)
{
	double x, y, z;
	x = x2 - x1;
	y = y2 - y1;
	z = z2 - z1;

	double D = sqrt(pow (x, 2) + pow (y, 2) + pow (z, 2));
	double a = x / D;
	double b = y / D;
	double c = z / D;

	double ad = sqrt(pow (b, 2) + pow (c, 2));
	double cosa = c / ad;
	double sina = b / ad;
	double cosb = ad;
	double sinb = a;
	double costheta = cos (angle*PI/180);
	double sintheta = sin (angle*PI/180);
	
	// initial the matrix
	double T [4][4] = {
		{1, 0, 0, 0},
		{0, 1, 0, 0}, 
		{0, 0, 1, 0},
		{-x1, -y1, -z1, 1}
	};

	// initial the matrix
	double NT [4][4] = {
		{1, 0, 0, 0},
		{0, 1, 0, 0}, 
		{0, 0, 1, 0},
		{x1, y1, z1, 1}
	};

	// initial the matrix
	double RA [4][4] = {
		{1, 0, 0, 0},
		{0, cosa, sina, 0}, 
		{0, -sina, cosa, 0},
		{0, 0, 0, 1}
	};

	// initial the matrix
	double NRA [4][4] = {
		{1, 0, 0, 0},
		{0, cosa, -sina, 0}, 
		{0, sina, cosa, 0},
		{0, 0, 0, 1}
	};

	// initial the matrix
	double RB [4][4] = {
		{cosb, 0, sinb, 0},
		{0, 1, 0, 0}, 
		{-sinb, 0, cosb, 0},
		{0, 0, 0, 1}
	};

	// initial the matrix
	double NRB [4][4] = {
		{cosb, 0, -sinb, 0},
		{0, 1, 0, 0}, 
		{sinb, 0, cosb, 0},
		{0, 0, 0, 1}
	};

	// initial the matrix
	double Rtheta [4][4] = {
		{costheta, sintheta, 0, 0},
		{-sintheta, costheta, 0, 0}, 
		{0, 0, 1, 0},
		{0, 0, 0, 1}
	};

	concatenate (T, RA);
	concatenate (matrixConcatenate, RB);
	concatenate (matrixConcatenate, Rtheta);
	concatenate (matrixConcatenate, NRB);
	concatenate (matrixConcatenate, NRA);
	concatenate (matrixConcatenate, NT);

	// print out the matrix 
	cout << "Arbitrary Rotation Matrix: " << endl;
	for (int row = 0; row < 4; row++)
	{
		for (int column = 0; column < 4; column++)
		{
			arbitraryRotation [row][column] = matrixConcatenate [row][column];
			cout << arbitraryRotation [row][column] << ", ";
		}
		cout << endl;
	}
	cout << endl;
}

// display the function implementation in the screen
void displayCube (void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	int turn;	
	inputLines (inputName);	// read datalines from an external file	
	orginalLines = true;
	basicTranslate (vx, vy, vz);
	matrixT3 (vx, vy);
	matrixT4 (vx, vy);
	basicScale (1, 1, -1);
	concatenate (matrixTranslate, T2);
	concatenate (matrixConcatenate, T3);
	concatenate (matrixConcatenate, T4);
	concatenate (matrixConcatenate, matrixScale);
	matrixV (matrixConcatenate);
	matrixN ();
	concatenate (V, N);
	applyTransformation (matrixConcatenate, dataW);
	output ();	
	for (int i = 0; i < lineNum; i++)
		bresenhamAlg (orginalLines, dashLength, dataS[i][0], dataS[i][1], dataS[i][2], dataS[i][3]);
}

// display the function implementation in the screen
void displayRotation ()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	int turn;	
	Rotation (angleRotate, rx1, ry1, rz1, rx2, ry2, rz2);
	inputLines (inputName);	// read datalines from an external file	
	orginalLines = true;
	basicTranslate (vx, vy, vz);
	matrixT3 (vx, vy);
	matrixT4 (vx, vy);
	basicScale (1, 1, -1);
	concatenate (matrixTranslate, T2);
	concatenate (matrixConcatenate, T3);
	concatenate (matrixConcatenate, T4);
	concatenate (matrixConcatenate, matrixScale);
	matrixV (matrixConcatenate);
	matrixN ();
	concatenate (V, N);
	applyTransformation (matrixConcatenate, dataW);
	applyTransformation (arbitraryRotation, dataC);
	output ();	
	for (int i = 0; i < lineNum; i++)
		bresenhamAlg (orginalLines, dashLength, dataS[i][0], dataS[i][1], dataS[i][2], dataS[i][3]);
}

int main (int argc, char** argv)
{
	
	// name of the input file name is provided by user
	cout << "Enter the datalines input file name: ";
	cin >> inputName;

	// number of lines is provided by user
	cout << "Enter the number of lines to transform from the datalines file: ";
	cin >> lineNum;

	// dash length is provided by user
	cout << "Enter the dash length: "; 
	cin >> dashLength;
	cout << endl;
	
	cout << "Enter the viewpoint x: ";
	cin >> vx;
	cout << "Enter the viewpoint y: ";
	cin >> vy;
	cout << "Enter the viewpoint z: ";
	cin >> vz;
	cout << "Enter the screen size: ";
	cin >> s;
	cout << "Enter the distance from the screen: ";
	cin >> d;
	cout << "Enter the screen resolution: ";
	cin >> screenPixel;
	/*
	cout << "Enter the arbitrary axis x1: ";
	cin >> rx1;
	cout << "Enter the arbitrary axis y1: ";
	cin >> ry1;
	cout << "Enter the arbitrary axis z1: ";
	cin >> rz1;
	cout << "Enter the arbitrary axis x2: ";
	cin >> rx2;
	cout << "Enter the arbitrary axis y2: ";
	cin >> ry2;
	cout << "Enter the arbitrary axis z2: ";
	cin >> rz2;
	cout << "Enter the rotation angle: ";
	cin >> angleRotate;
	*/
	
	glutInit (&argc, argv); // initialize GLUT
	glutInitDisplayMode (GLUT_SINGLE | GLUT_RGB); // set display mode
	glutInitWindowPosition (100, 100); // set top-left display-window position
	glutInitWindowSize (800, 600); // set display-window sidth and height
	glutCreateWindow ("Screen"); // create display window

	init (800, 600); // execute initialization procedure
	glutDisplayFunc (displayCube); // send graphics to display window
	//glutDisplayFunc (displayRotation); // send graphics to display window
	glutMainLoop (); // display everything and wait
	
}
