#include <Windows.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <stack>
#include <commctrl.h>
#include <shellapi.h>
#include <windowsx.h>
#include <string>
#include <tchar.h>
#include <ctime> // for time and strftime
#include <cstring> // for strcpy
#include <vector>
#include <algorithm>

using namespace std; 

// // Declare global variables
// Add this declaration to the top of the file

//HWND hwndText ;
//HWND hwndSubmit;
HWND hCombo = NULL;

// this declere what user select from drop down list
static int Selection = -1;


// colors
const COLORREF initColor = 100;
const COLORREF red = RGB(255, 0, 0);
const COLORREF blue = RGB(0, 0, 255);
const COLORREF green = RGB(0, 255, 0);
const COLORREF black = RGB(0, 0, 0);
static COLORREF newColor = initColor, lastColor = initColor;

// number of user click to draw 
static int number_of_click = 0;

// pos user click 
static int xx1, yy1, xx2, yy2,xx3,yy3;

// her we have an error **********************************************************
const int MAXSCANLINES = 1920;
//const int MAXSCANLINES = GetSystemMetrics(SM_CYSCREEN);



struct Vertex
{
	int x, y;
	Vertex(){}
	Vertex(int x, int y) :x(x), y(y)
	{
	}
};


// array of vertex of user click
const int pointSize = 5;
static  Vertex points[pointSize];



int Round(double x) {
    return (int)(x + 0.5);
}
 

// Implement line algorithms[DDA, Midpoint and parametric]

void LineDD(HDC hdc, int x1, int y1, int x2, int y2, COLORREF c) {
	int dx = x2 - x1, dy = y2 - y1;
	if (abs(dx) >= abs(dy)) {
		if (x1 > x2) {
			int temp = x1;
			x1 = x2;
			x2 = temp;
			temp = y1;
			y1 = y2;
			y2 = temp;
		}
		int x = x1;
		double y = y1;
		double m = (double)dy / dx;
		SetPixel(hdc, x1, y1, c);
		while (x < x2) {
			x++;
			y += m;
			SetPixel(hdc, x, Round(y), c);
		}
	}
	else {
		if (y1 > y2) {
			int temp = y1;
			y1 = y2;
			y2 = temp;
			temp = x1;
			x1 = x2;
			x2 = temp;
		}
		double x = x1;
		int y = y1;
		SetPixel(hdc, x1, y1, c);
		double m = (double)dx / dy;
		while (y < y2)
		{
			y++;
			x += m;
			SetPixel(hdc, Round(x), y, c);
		}

	}
}

void DrawLine1(HDC hdc, int x1, int y1, int x2, int y2, COLORREF c) {

	int dx = x2 - x1, dy = y2 - y1;
	if (abs(dx) >= abs(dy)) {
		if (x1 > x2) {
			int temp = x1;
			x1 = x2;
			x2 = temp;
			temp = y1;
			y1 = y2;
			y2 = temp;
		}

		int x = x1;
		double y;
		while (x <= x2) {
			y = (double)(((double)(y2 - y1) / (double)(x2 - x1)) * (double)(x - x1)) + y1;
			SetPixel(hdc, x, Round(y), c);
			x++;
		}
	}
	else if (abs(dx) < abs(dy)) {
		if (y1 > y2) {
			int temp = y1;
			y1 = y2;
			y2 = temp;
			temp = x1;
			x1 = x2;
			x2 = temp;

		}

		int y = y1;
		double x;
		while (y <= y2)
		{
			x = (double)(((double)(x2 - x1) / (double)(y2 - y1)) * (double)(y - y1)) + x1;
			SetPixel(hdc, Round(x), y, c);
			y++;
		}

	}
}

void DrawlineB(HDC hdc, int x1, int y1, int x2, int y2, COLORREF c) {

	int dx = x2 - x1, dy = y2 - y1;
	if (abs(dx) >= abs(dy)) {
		if (x1 > x2) {
			int temp = x1;
			x1 = x2;
			x2 = temp;
			int temp2 = y1;
			y1 = y2;
			y2 = temp2;
			dx = -dx; dy = -dy;
		}
		if (dy >= 0) {
			int x = x1, y = y1;
			SetPixel(hdc, x, y, c);
			int d_initial = dx - 2 * dy;
			int change_1 = -2 * dy;
			int change_2 = 2 * (dx - dy);
			SetPixel(hdc, x, y, c);
			while (x < x2) {

				if (d_initial > 0) {
					x++;
					d_initial += change_1;
				}
				else {
					x++, y++;
					d_initial += change_2;
				}
				SetPixel(hdc, x, y, c);
			}

		}
		else {
			int x = x1, y = y1;
			SetPixel(hdc, x, y, c);
			int d_initial = dx - 2 * abs(dy);
			int change_1 = -2 * abs(dy);
			int change_2 = 2 * (dx - abs(dy));
			SetPixel(hdc, x, y, c);
			while (x < x2) {
				if (d_initial > 0) {
					x++;
					d_initial += change_1;
				}
				else {
					x++, y--;
					d_initial += change_2;
				}
				SetPixel(hdc, x, y, c);
			}
		}



	}
	else {
		if (y1 > y2)
		{
			int temp = x1;
			x1 = x2;
			x2 = temp;
			int temp2 = y1;
			y1 = y2;
			y2 = temp2;
			dx = -dx; dy = -dy;
		}
		if (dx >= 0) {
			int x = x1, y = y1;
			SetPixel(hdc, x, y, c);
			int d_initial = dy - 2 * dx;
			int change_1 = -2 * dx;
			int change_2 = 2 * (dy - dx);
			SetPixel(hdc, x, y, c);
			while (y < y2) {

				if (d_initial >= 0) {
					y++;
					d_initial += change_1;
				}
				else {
					x++, y++;
					d_initial += change_2;
				}
				SetPixel(hdc, x, y, c);
			}
		}
		else {
			int x = x1, y = y1;
			SetPixel(hdc, x, y, c);
			int d_initial = dy - 2 * abs(dx);
			int change_1 = -2 * abs(dx);
			int change_2 = 2 * (dy - abs(dx));
			SetPixel(hdc, x, y, c);
			while (y < y2) {

				if (d_initial > 0) {
					y++;
					d_initial += change_1;
				}
				else {
					x--, y++;
					d_initial += change_2;
				}
				SetPixel(hdc, x, y, c);
			}
		}
	}
}

// ------------------------------------------------------------------------------------------------------------------------





// Implement Circle algorithms (Direct, Polar, iterative Polar, midpoint and modified Midpoint)

void Draw4Points(HDC hdc, int xc, int yc, int x, int y, COLORREF c)
{
	SetPixel(hdc, xc + x, yc + y, c);
	SetPixel(hdc, xc + x, yc - y, c);
	SetPixel(hdc, xc - x, yc - y, c);
	SetPixel(hdc, xc - x, yc + y, c);

}

void Draw8Points(HDC hdc, int xc, int yc, int a, int b, COLORREF color)
{
	SetPixel(hdc, xc + a, yc + b, color);
	SetPixel(hdc, xc - a, yc + b, color);
	SetPixel(hdc, xc - a, yc - b, color);
	SetPixel(hdc, xc + a, yc - b, color);
	SetPixel(hdc, xc + b, yc + a, color);
	SetPixel(hdc, xc - b, yc + a, color);
	SetPixel(hdc, xc - b, yc - a, color);
	SetPixel(hdc, xc + b, yc - a, color);
}

// Direct method to draw circle
void CircleDirect(HDC hdc,int xc , int yc, int radius, COLORREF color) {
	int radius_square = radius * radius;
	int x = 0, y = radius;
	
	Draw8Points(hdc, xc, yc, x, y, color);
	while (x < y) {
		x++;
		y = (int)sqrt(radius_square - (x * x)); // Compute new y
		
		Draw8Points(hdc, xc, yc, x, y, color);
	}
}

void CircleIterativePolar(HDC hdc, int xc, int yc, int R, COLORREF color)
{
	double x = R, y = 0;
	double dtheta = 1.0 / R;
	double cdtheta = cos(dtheta), sdtheta = sin(dtheta);
	Draw8Points(hdc, xc, yc, R, 0, color);

	while (x > y)
	{
		double x1 = x * cdtheta - y * sdtheta;
		y = x * sdtheta + y * cdtheta;
		x = x1;
		Draw8Points(hdc, xc, yc, round(x), round(y), color);
	}


}

void CircleFasterBresenham(HDC hdc, int xc, int yc, int R, COLORREF color)
{
	int x = 0, y = R;
	int d = 1 - R;
	int c1 = 3, c2 = 5 - 2 * R;
	Draw8Points(hdc, xc, yc, x, y, color);
	while (x < y)
	{
		if (d < 0)
		{
			d += c1;
			c2 += 2;
		}
		else
		{

			d += c2;
			c2 += 4;
			y--;
		}
		c1 += 2;
		x++;

		Draw8Points(hdc, xc, yc, x, y, color);
	}
}

// Direct Polar form to draw circle
void DirectPolar(HDC hdc, int xc , int yc, int radius, COLORREF color) {
	int x = radius, y = 0;
	double theta = 0.0, delta_theta = 1.0 / radius;
	Draw8Points(hdc, xc, yc, x, y, color);
	while (x > y) {
		x = (int)(radius * cos(theta));
		y = (int)(radius * sin(theta));
		theta += delta_theta;
		Draw8Points(hdc, xc, yc, x, y, color);
	}
}

void CircleBresenham(HDC hdc , int xc, int yc, int radius, COLORREF color) {
	int x = 0, y = radius;
	int d = 1 - radius;
	Draw8Points(hdc, xc, yc, x, y, color);
	while (x < y) {
		if (d < 0)
			d += 2 * x + 2;
		else {
			d += 2 * (x - y) + 5;
			y--;
		}
		x++;
		Draw8Points(hdc, xc, yc, x, y, color);
	}
}


// ------------------------------------------------------------------------------------------------------------------------




// Recursive and Non Recursive Flood Fill

void fFill(HDC hdc, int x, int y, COLORREF bc, COLORREF fc) {

	COLORREF c = GetPixel(hdc, x, y);
	if (c == bc || c == fc || c == red || c == blue || c == green || c == black || c == initColor) { return; }
	SetPixel(hdc, x, y, fc);
	fFill(hdc, x + 1, y, bc, fc);
	fFill(hdc, x - 1, y, bc, fc);
	fFill(hdc, x, y + 1, bc, fc);
	fFill(hdc, x, y - 1, bc, fc);
}

void NRFloodFill(HDC hdc, int x, int y, COLORREF Cb, COLORREF Cf)
{
	stack<Vertex> S;
	S.push(Vertex(x, y));
	while (!S.empty())
	{
		Vertex v = S.top();
		S.pop();
		COLORREF c = GetPixel(hdc, v.x, v.y);
		if (c == Cb || c == Cf || c == red || c == blue || c == green || c == black || c == initColor)continue;
		SetPixel(hdc, v.x, v.y, Cf);
		S.push(Vertex(v.x + 1, v.y));
		S.push(Vertex(v.x - 1, v.y));
		S.push(Vertex(v.x, v.y + 1));
		S.push(Vertex(v.x, v.y - 1));
	}
}

// ------------------------------------------------------------------------------------------------------------------------





// Filling Circle with other circles after taking filling quarter from user

void quarterDrawPoints(HDC hdc, int xc, int yc, int a, int b, int quarter_number, COLORREF color) {

	if (quarter_number==1)
	{
		SetPixel(hdc, xc + b, yc - a, color);
	}
	else if(quarter_number == 2)
	{
		
		SetPixel(hdc, xc + a, yc - b, color);
	}
	else if (quarter_number == 3)
	{
		SetPixel(hdc, xc - a, yc - b, color);
	}
	else if (quarter_number == 4)
	{
		
		SetPixel(hdc, xc - b, yc - a, color);
	}
	else if (quarter_number == 5)
	{
		
		SetPixel(hdc, xc - b, yc + a, color);
	}
	else if (quarter_number == 6)
	{
		
		SetPixel(hdc, xc - a, yc + b, color);
	}
	else if (quarter_number == 7)
	{
		
		SetPixel(hdc, xc + a, yc + b, color);
	}
	else if (quarter_number == 8)
	{
		
		SetPixel(hdc, xc + b, yc + a, color);
	}

}

void quarterCircleFaster(HDC hdc, int xc, int yc, int R, int quarter_number, COLORREF color) {
	int x = 0, y = R;
	int d = 1 - R;
	int c1 = 3, c2 = 5 - 2 * R;
	quarterDrawPoints(hdc, xc, yc, x, y, quarter_number,color);
	while (x < y)
	{
		if (d < 0)
		{
			d += c1;
			c2 += 2;
		}
		else
		{

			d += c2;
			c2 += 4;
			y--;
		}
		c1 += 2;
		x++;

		quarterDrawPoints(hdc, xc, yc, x, y,quarter_number ,color);
	}
}

void Filling_Circle_other_circles(HDC hdc, int xc, int yc, int R,int quarter_number,COLORREF color) {
	CircleFasterBresenham(hdc, xc, yc, R, color);
	while (R>0)
	{
		quarterCircleFaster(hdc, xc, yc, R, quarter_number, color);
		R = R - 2;
	}
	

}

// ------------------------------------------------------------------------------------------------------------------------





// Filling Circle with lines after taking filling quarter from user
void quarterDrawLines(HDC hdc, int xc, int yc, int x, int y, int quarter_number, COLORREF color) {
	if (quarter_number == 1)
	{
		
		DrawlineB(hdc, xc, yc, xc + y, yc - x, color);
	}
	else if (quarter_number == 2)
	{

		DrawlineB(hdc, xc, yc, xc + x, yc - y, color);
	}
	else if (quarter_number == 3)
	{
		DrawlineB(hdc, xc, yc, xc - x, yc - y, color);
	}
	else if (quarter_number == 4)
	{

		DrawlineB(hdc, xc, yc, xc - y, yc - x, color);
	}
	else if (quarter_number == 5)
	{

		DrawlineB(hdc, xc, yc, xc - y , yc + x, color);
	}
	else if (quarter_number == 6)
	{

		DrawlineB(hdc, xc, yc, xc - x, yc + y, color);
	}
	else if (quarter_number == 7)
	{


		DrawlineB(hdc, xc, yc, xc + x, yc + y, color);
	}
	else if (quarter_number == 8)
	{

		DrawlineB(hdc, xc, yc, xc + y, yc + x, color);
	}

}

void Filling_Circle_lines(HDC hdc, int xc, int yc, int R,int quarter_number, COLORREF color)
{
	int x = 0, y = R;
	int d = 1 - R;
	int c1 = 3, c2 = 5 - 2 * R;
	Draw8Points(hdc, xc, yc, x, y, color);
	while (x < y)
	{
		if (d < 0)
		{
			d += c1;
			c2 += 2;
		}
		else
		{

			d += c2;
			c2 += 4;
			y--;
		}
		c1 += 2;
		x++;
		Draw8Points(hdc, xc, yc, x, y, color);
		quarterDrawLines(hdc, xc, yc, x, y, quarter_number,color);

	}
}


// ------------------------------------------------------------------------------------------------------------------------






// Ellipse Algorithms [Direct, polar and midpoint] 

void DirectEllipse(HDC hdc, int xc, int yc, int a, int b, COLORREF c)
{

}

void DrawEllipseMidpoint(HDC hdc, int xc, int yc, int A, int B, COLORREF c)
{
	int x = 0;
	double y = B;
	Draw4Points(hdc, xc, yc, 0, B, c);
	while (x * B * B < y * A * A)
	{
		x++;
		y = B * sqrt(1.0 - (double)x * x / (A * A));
		Draw4Points(hdc, xc, yc, x, Round(y), c);
	}
	int y1 = 0;
	double x1 = A;
	Draw4Points(hdc, xc, yc, A, 0, c);
	while (x1 * B * B > y1 * A * A)
	{
		y1++;
		x1 = A * sqrt(1.0 - (double)y1 * y1 / (B * B));
		Draw4Points(hdc, xc, yc, Round(x1), y1, c);
	}

}


void DrawEllipsePoler(HDC hdc, int xc, int yc, int A, int B, COLORREF c)
{

	double theta = 0;
	const double PI = 3.14159265358979323846;
	while (theta < PI / 2)
	{
		int x = round(A * cos(theta));
		int y = round(B * sin(theta));
		Draw4Points(hdc, xc, yc, x, y, c);

		theta += 0.001;
	}
}




// ------------------------------------------------------------------------------------------------------------------------



// Convex and Non Convex Filling Algorithm




void SwapPoints(Vertex p1, Vertex p2) {
	int tmp;
	tmp = p1.x;
	p1.x = p2.x;
	p2.x = tmp;

	tmp = p1.y;
	p1.y = p2.y;
	p2.y = tmp;
}

// Each record in the table is Table object
class Table {
public:
	int xLeft, xRight;

	Table() {
		xLeft = MAXINT;
		xRight = -MAXINT;
	}
};

void DrawHorizontalLine(HDC hdc, Vertex start, Vertex end, COLORREF color) {
	if (start.x > end.x) SwapPoints(start, end);
	int x = start.x;
	int y = start.y;
	while (x <= end.x) {
		SetPixel(hdc, x++, y, color);
	}
}

// Define utility function that takes the end points of the edge,
// generates all points of intersection with scan lines and update the entry table
void Edge2Table(Vertex p1, Vertex p2, Table table[]) {
	// I have a point from the edge
	// and I want to know if it is on the left of the polygon
	// or on the right of the polygon
	if (p1.y == p2.y) return;

	if (p1.y > p2.y) {
		Vertex temp = p1;
		p1 = p2;
		p2 = temp;
	}

	// Inverse of the slope
	double mInverse = (double)(p2.x - p1.x) / (p2.y - p1.y);

	// Define start point (xMin, y) = (p1.xMin, p1.y)
	double x = p1.x;
	int y = p1.y;
	while (y < p2.y) {
		// ceil and floor with left and right to make the pixel inside polygon
		if (x < table[y].xLeft) table[y].xLeft = (int)ceil(x);
		if (x > table[y].xRight) table[y].xRight = (int)floor(x);
		y++;
		x += mInverse;
	}
}

void Polygon2Table(Vertex points[], int pointsTableSize, Table table[]) {
	// Start with last point
	Vertex p1 = points[pointsTableSize - 1];
	for (int i = 0; i < pointsTableSize; ++i) {
		Vertex p2 = points[i];
		Edge2Table(p1, p2, table);
		p1 = p2; // Last drawn point is the start point of the next iteration
	}
}

void Table2Screen(HDC hdc, Table table[], COLORREF color) {
	for (int y = 0; y < MAXSCANLINES; ++y) {
		if (table[y].xLeft < table[y].xRight) {
			Vertex p1(table[y].xLeft, y);
			Vertex p2(table[y].xRight, y);
			DrawHorizontalLine(hdc, p1, p2, color);
		}
	}
}

// Main function
void ConvexPolygonFilling(HDC hdc, Vertex points[], int pointsSize, COLORREF color) {
	Table table[MAXSCANLINES];
	Polygon2Table(points, pointsSize, table); // Fill the table
	Table2Screen(hdc, table, color);
}





// ------------------------------------------------------------------------------------------------------------------------





class EdgeRecord {
public:
	double xMin, mInverse;
	int yMax;

	// Default constructor
	EdgeRecord() {
		xMin = 0.0;
		mInverse = 0.0;
		yMax = 0;
	}

	void setEdgeRecord(Vertex v1, Vertex v2) {
		if (v1.y > v2.y)
			SwapPoints(v1, v2);
		this->xMin = v1.x;
		this->mInverse = (double)(v2.x - v1.x) / (v2.y - v1.y);
		this->yMax = v2.y;
	}

	void deleteEdgeRecord() {
		xMin = 0.0;
		mInverse = 0.0;
		yMax = 0;
	}
};

vector<vector<EdgeRecord>> table(MAXSCANLINES);

// Utility function used in sorting active list based on nodes' xMin value
bool sortNodes(EdgeRecord node1, EdgeRecord node2) {
	return node1.xMin < node2.xMin;
}

// This function records edges into table
void Edges2Table(Vertex polygon[], int polygonSize) {
	Vertex v1 = polygon[polygonSize - 1];
	for (int i = 0; i < polygonSize; ++i) {
		Vertex v2 = polygon[i];
		if (v1.y == v2.y) {
			v1 = v2;
			continue;
		}
		EdgeRecord edgeRecord;
		edgeRecord.setEdgeRecord(v1, v2);
		table[min(v1.y, v2.y)].push_back(edgeRecord);
		v1 = v2;
	}
}

// Append edge records to active list
void appendActiveList(vector<EdgeRecord>& active, int y) {
	for (auto it : table[y]) {
		active.push_back(it);
		it.deleteEdgeRecord();
	}
}

void GeneralPolygonFilling(HDC hdc, Vertex polygon[], int polygonSize, COLORREF color) {
	Edges2Table(polygon, polygonSize); // Fill table with edges records

	int y = 0; // The first index in table with non-empty edge vector entry

	// Get first non-empty vector entry
	// when this loop ends, y will be the index of first non-empty edge vector entry
	while (y < MAXSCANLINES && table[y].empty()) y++;

	if (y == MAXSCANLINES) return; // If there is any full entry

	int deleted = 0;
	vector<EdgeRecord> activeList = table[y]; // Define active list
	while (activeList.size() != deleted) {
		// 4.1. Sort nodes in active list based on their x values
		sort(activeList.begin(), activeList.end(), sortNodes);
		auto it = activeList.begin();

		// 4.2. Draw horizontal lines between points represented by successive pairs of nodes in activeList

		for (it = activeList.begin() + deleted; it != activeList.end(); it++) {
			int x1 = (int)ceil(it->xMin);
			it++;
			int x2 = (int)floor(it->xMin);
			for (int x = x1; x < x2; ++x) {
				SetPixel(hdc, x, y, color);
			}
		}

		// 4.3. Increment y by 1 (to go to a new scan line)
		y++;

		// 4.4. Delete from ActiveList those nodes with yMax = y
		it = activeList.begin();
		while (it != activeList.end()) {
			if (y == it->yMax) {
				it->deleteEdgeRecord();
				deleted++;
			}
			it++;
		}
		// 4.5. Increment by mInverse
		for (it = activeList.begin(); it != activeList.end(); it++) {
			it->xMin += it->mInverse;
		}
		// 4.6. Append to ActiveList the new nodes at table[y] (if any)
		appendActiveList(activeList, y);
	}
}



// ------------------------------------------------------------------------------------------------------------------------


//  Filling Square with Hermit Curve[Vertical]
//  Filling Rectangle with Bezier Curve[horizontal]
//  Cardinal Spline Curve


struct Vector2 {
	double x, y;

	Vector2(double a = 0, double b = 0) {
		x = a;
		y = b;
	}
};

class Vector4 {
	double v[4];
public:
	Vector4(double a = 0, double b = 0, double c = 0, double d = 0) {
		v[0] = a;
		v[1] = b;
		v[2] = c;
		v[3] = d;
	}

	Vector4(double a[]) {
		memcpy(v, a, 4 * sizeof(double));
	}

	double& operator[](int i) {
		return v[i];
	}
};

class Matrix4 {
	Vector4 M[4];
public:
	Matrix4(double A[]) {
		memcpy(M, A, 16 * sizeof(double));
	}

	Vector4& operator[](int i) {
		return M[i];
	}
};

Vector4 operator*(Matrix4 M, Vector4& b) // right multiplication of M by b
{
	Vector4 res;
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			res[i] += M[i][j] * b[j];
	return res;
}

double DotProduct(Vector4& a, Vector4& b) //multiplying a raw vector by a column vector
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3];
}

Vector4 GetHermiteCoeff(double x0, double s0, double x1, double s1) {
	static double H[16] = { 2, 1, -2, 1, -3, -2, 3, -1, 0, 1, 0, 0, 1, 0, 0, 0 };
	static Matrix4 basis(H);
	Vector4 v(x0, s0, x1, s1);
	return basis * v;
}

void DrawHermiteCurve(HDC hdc, Vector2& P0, Vector2& T0, Vector2& P1, Vector2& T1, int numpoints) {
	Vector4 xcoeff = GetHermiteCoeff(P0.x, T0.x, P1.x, T1.x);
	Vector4 ycoeff = GetHermiteCoeff(P0.y, T0.y, P1.y, T1.y);
	if (numpoints < 2)return;
	double dt = 1.0 / (numpoints - 1);
	for (double t = 0; t <= 1; t += dt) {
		Vector4 vt;
		vt[3] = 1;
		for (int i = 2; i >= 0; i--) {
			vt[i] = vt[i + 1] * t;
		}
		int x = (int)(DotProduct(xcoeff, vt));
		int y = (int)(DotProduct(ycoeff, vt));
		if (t == 0) MoveToEx(hdc, x, y, nullptr);
		else LineTo(hdc, x, y);
	}
}

void DrawBezierCurve(HDC hdc, Vector2& P0, Vector2& P1, Vector2& P2, Vector2& P3, int numpoints) {
	Vector2 T0(3 * (P1.x - P0.x), 3 * (P1.y - P0.y));
	Vector2 T1(3 * (P3.x - P2.x), 3 * (P3.y - P2.y));
	DrawHermiteCurve(hdc, P0, T0, P3, T1, numpoints);
}

void DrawCardinalSpline(HDC hdc, Vector2 P[], int n, double c, int numpix) {
	double c1 = 1 - c;
	Vector2 T0(c1 * (P[2].x - P[0].x), c1 * (P[2].y - P[0].y));
	int i = 1;
	for (; i < n - 1; i++) {
		Vector2 T1(c1 * (P[i + 1].x - P[i - 1].x), c1 * (P[i + 1].y - P[i - 1].y));
		DrawHermiteCurve(hdc, P[i - 1], T0, P[i], T1, numpix);
		T0 = T1;
	}
	Vector2 T1(c1 * (P[i + 1].x - P[i].x), c1 * (P[i + 1].y - P[i].y));
	DrawHermiteCurve(hdc, P[i - 1], T0, P[i], T1, numpix);
}

void FillSquareWithHermitVertical(HDC hdc, Vector2& start, Vector2& end, COLORREF color) {
	end.y = start.y + abs(start.x - end.x);
	int m = 3;
	while (end.y > start.y + m) {
		Vector2 P0(start.x + m - 3, start.y);
		Vector2 P1(start.x + m, start.y);
		Vector2 P2(start.x + m - 3, end.y);
		Vector2 P3(start.x + m, end.y);
		Vector2 T0(3 * (P1.x - P0.x), 3 * (P1.y - P0.y));
		Vector2 T1(3 * (P3.x - P2.x), 3 * (P3.y - P2.y));
		DrawHermiteCurve(hdc, P0, T0, P3, T1, 20);
		m += 3;
	}
}

void FillRectangleWithBezierHorizontal(HDC hdc, Vector2& start, Vector2& end, COLORREF color)
{
	int m = 3;
	while (end.y > start.y + m) {
		Vector2 P0(start.x, end.y - m + 3);
		Vector2 P1(start.x, end.y - m);
		Vector2 P2(end.x, end.y - m + 3);
		Vector2 P3(end.x, end.y - m);
		DrawBezierCurve(hdc, P0, P1, P2, P3, 20);
		m += 2;
	}
}



// ------------------------------------------------------------------------------------------------------------------------

// Clipping algorithms using Rectangle as Clipping Window[Point ,Line, Polygon]
// Clipping algorithms using Square as Clipping Window[Point, Line]


void ClipPoint(HDC hdc, int x, int y, int xLeft, int yTop, int xRight, int yBottom, COLORREF color) {
	if (x >= xLeft && x <= xRight && y >= yTop && y <= yBottom)
		SetPixel(hdc, x, y, color);
}

union OutCode {
	unsigned All : 4;
	struct {
		unsigned left : 1, top : 1, right : 1, bottom : 1;
	};
};

OutCode GetOutCode(double x, double y, int xLeft, int yTop, int xRight, int yBottom) {
	OutCode out;
	out.All = 0;
	if (x < xLeft)
		out.left = 1;
	else if (x > xRight)
		out.right = 1;
	if (y < yTop)
		out.top = 1;
	else if (y > yBottom)
		out.bottom = 1;
	return out;
}

void VIntersect(double xs, double ys, double xe, double ye, int x, double* xi, double* yi) {
	*xi = x;
	*yi = ys + (x - xs) * (ye - ys) / (xe - xs);
}

void HIntersect(double xs, double ys, double xe, double ye, int y, double* xi, double* yi) {
	*yi = y;
	*xi = xs + (y - ys) * (xe - xs) / (ye - ys);
}

void CohenSutherland(HDC hdc, int xs, int ys, int xe, int ye, int xLeft, int yTop, int xRight, int yBottom) {
	double x1 = xs, y1 = ys, x2 = xe, y2 = ye;
	OutCode out1 = GetOutCode(x1, y1, xLeft, yTop, xRight, yBottom);
	OutCode out2 = GetOutCode(x2, y2, xLeft, yTop, xRight, yBottom);
	while ((out1.All || out2.All) && !(out1.All & out2.All)) {
		double xi, yi;
		if (out1.All) {
			if (out1.left)
				VIntersect(x1, y1, x2, y2, xLeft, &xi, &yi);
			else if (out1.top)
				HIntersect(x1, y1, x2, y2, yTop, &xi, &yi);
			else if (out1.right)
				VIntersect(x1, y1, x2, y2, xRight, &xi, &yi);
			else if (out1.bottom)
				HIntersect(x1, y1, x2, y2, yBottom, &xi, &yi);
			x1 = xi;
			y1 = yi;
			out1 = GetOutCode(x1, y1, xLeft, yTop, xRight, yBottom);
		}
		else {
			if (out2.left)
				VIntersect(x1, y1, x2, y2, xLeft, &xi, &yi);
			else if (out2.top)
				HIntersect(x1, y1, x2, y2, yTop, &xi, &yi);
			else if (out2.right)
				VIntersect(x1, y1, x2, y2, xRight, &xi, &yi);
			else if (out2.bottom)
				HIntersect(x1, y1, x2, y2, yBottom, &xi, &yi);
			x2 = xi;
			y2 = yi;
			out2 = GetOutCode(x2, y2, xLeft, yTop, xRight, yBottom);
		}
	}
	if (!out1.All && !out2.All) {
		MoveToEx(hdc, (int)x1, (int)y1, nullptr);
		LineTo(hdc, (int)x2, (int)y2);
	}
}

void ClipPolygon(HDC hdc, vector<Vertex> polygon, int size, int xLeft, int yTop, int xRight, int yBottom) {
	Vertex p1 = polygon[size - 1];
	for (int i = 0; i < size; ++i) {
		Vertex p2 = polygon[i];
		CohenSutherland(hdc, p1.x, p1.y, p2.x, p2.y, xLeft, yTop, xRight, yBottom);
		p1 = p2; // Last drawn point is the start point of the next iteration
	}
}



// ------------------------------------------------------------------------------------------------------------------------

// Clipping algorithms using circle as a Clipping Window [Point, Line]

void Draw8Points(HDC hdc, Vertex center, Vertex p, COLORREF color) {
	SetPixel(hdc, center.x + p.x, center.y + p.y, color);

	SetPixel(hdc, center.x - p.x, center.y + p.y, color);
	SetPixel(hdc, center.x - p.x, center.y - p.y, color);
	SetPixel(hdc, center.x + p.x, center.y - p.y, color);

	SetPixel(hdc, center.x + p.y, center.y + p.x, color);
	SetPixel(hdc, center.x - p.y, center.y + p.x, color);
	SetPixel(hdc, center.x - p.y, center.y - p.x, color);
	SetPixel(hdc, center.x + p.y, center.y - p.x, color);
}

// This function computes the radius of specific circle
int ComputeRadius(Vertex center, Vertex p) {
	return (int)sqrt((p.y - center.y) * (p.y - center.y) + (p.x - center.x) * (p.x - center.x));
}

bool checkPoint(Vertex center, Vertex p, int r) {
	return (p.x - center.x) * (p.x - center.x) + (p.y - center.y) * (p.y - center.y) <= r * r;
}

// Modified Bresenham algorithm
void ModifiedMidpoint(HDC hdc, Vertex center, int radius, COLORREF color) {
	int x = 0, y = radius;
	int d = 1 - radius;
	int c1 = 3, c2 = 5 - 2 * radius;
	Vertex p(x, y);
	Draw8Points(hdc, center, p, color);
	while (x < y) {
		if (d < 0) {
			d += c1;
			c2 += 2;
		}
		else {
			d += c2;
			c2 += 4;
			y--;
		}
		c1 += 2;
		x++;
		p.x = x;
		p.y = y;
		Draw8Points(hdc, center, p, color);
	}
}

void LineDD(HDC hdc, Vertex center, int r, int x1, int y1, int x2, int y2, COLORREF c) {
	int dx = x2 - x1, dy = y2 - y1;
	Vertex p;
	if (abs(dx) >= abs(dy)) {
		if (x1 > x2) {
			int temp = x1;
			x1 = x2;
			x2 = temp;
			temp = y1;
			y1 = y2;
			y2 = temp;
		}
		int x = x1;
		double y = y1;
		double m = (double)dy / dx;
		p.x = x;
		p.y = (int)y;
		if (checkPoint(center, p, r)) SetPixel(hdc, x, (int)y, c);

		while (x < x2) {
			x++;
			y += m;
			p.x = x;
			p.y = (int)y;
			if (checkPoint(center, p, r)) SetPixel(hdc, x, (int)y, c);
		}
	}
	else {
		if (y1 > y2) {
			int temp = y1;
			y1 = y2;
			y2 = temp;
			temp = x1;
			x1 = x2;
			x2 = temp;
		}
		double x = x1;
		int y = y1;
		p.x = (int)x;
		p.y = y;
		if (checkPoint(center, p, r)) SetPixel(hdc, (int)x, y, c);
		double m = (double)dx / dy;
		while (y < y2) {
			y++;
			x += m;
			p.x = (int)x;
			p.y = y;
			if (checkPoint(center, p, r)) SetPixel(hdc, (int)x, y, c);
		}

	}
}


// ------------------------------------------------------------------------------------------------------------------------


// Implement save function for all data in screen
// Implement load function to load data from files


// Saving and Loading To-From File
bool HDCToFile(const char* FilePath, HWND& hWnd, HDC Context, RECT Area, uint16_t BitsPerPixel = 24)
{
	//MessageBox(hWnd,(LPCWSTR)FilePath, L"comfirmation", MB_OKCANCEL);
	uint32_t Width = Area.right - Area.left;
	uint32_t Height = Area.bottom - Area.top;
	BITMAPINFO Info;
	BITMAPFILEHEADER Header;
	memset(&Info, 0, sizeof(Info));
	memset(&Header, 0, sizeof(Header));
	Info.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
	Info.bmiHeader.biWidth = Width;
	Info.bmiHeader.biHeight = Height;
	Info.bmiHeader.biPlanes = 1;
	Info.bmiHeader.biBitCount = BitsPerPixel;
	Info.bmiHeader.biCompression = BI_RGB;
	Info.bmiHeader.biSizeImage = Width * Height * (BitsPerPixel > 24 ? 4 : 3);
	Header.bfType = 0x4D42;
	Header.bfOffBits = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER);
	char* Pixels = NULL;
	HDC MemDC = CreateCompatibleDC(Context);
	HBITMAP Section = CreateDIBSection(Context, &Info, DIB_RGB_COLORS, (void**)&Pixels, 0, 0);
	DeleteObject(SelectObject(MemDC, Section));
	BitBlt(MemDC, 0, 0, Width, Height, Context, Area.left, Area.top, SRCCOPY);
	DeleteDC(MemDC);
	std::fstream hFile(FilePath, std::ios::out | std::ios::binary);
	if (hFile.is_open())
	{
		hFile.write((char*)&Header, sizeof(Header));
		hFile.write((char*)&Info.bmiHeader, sizeof(Info.bmiHeader));
		hFile.write(Pixels, (((BitsPerPixel * Width + 31) & ~31) / 8) * Height);
		hFile.close();
		DeleteObject(Section);
		return true;
	}
	else {
		MessageBox(hWnd,L"Error, Can not open file", L"Error", MB_OKCANCEL);
	}
	DeleteObject(Section);
	return false;
}

void load(HWND hWnd, HDC& hdc)
{

	TCHAR szFile[MAX_PATH] = { 0 }; // buffer to hold file path
	OPENFILENAME ofn = { 0 }; // initialize structure
	ofn.lStructSize = sizeof(ofn);
	ofn.hwndOwner = hWnd;
	ofn.lpstrFilter = L"Bitmap Image (*.bmp)\0*.bmp\0All Files (*.*)\0*.*\0"; // filter for file types
	ofn.lpstrFile = szFile;
	ofn.nMaxFile = MAX_PATH;
	ofn.Flags = OFN_EXPLORER | OFN_OVERWRITEPROMPT; // options for dialog
	// open dialog to get file name
	if (GetOpenFileName(&ofn))
	{
		HDC hdc = GetDC(hWnd);
		int windowWidth;
		int windowHeight;
		RECT rect;
		if (GetWindowRect(hWnd, &rect))
		{
			windowWidth = rect.right - rect.left;
			windowHeight = rect.bottom - rect.top;
		}
		RECT rect1 = { 0, 0, windowWidth, windowHeight };
	}

	HBITMAP hBitmap;
	hBitmap = (HBITMAP)::LoadImage(NULL, szFile, IMAGE_BITMAP, 0, 0, LR_LOADFROMFILE);
	HDC hLocalDC;
	hLocalDC = CreateCompatibleDC(hdc);
	BITMAP qBitmap;
	int iReturn = GetObject(reinterpret_cast<HGDIOBJ>(hBitmap), sizeof(BITMAP), reinterpret_cast<LPVOID>(&qBitmap));
	HBITMAP hOldBmp = (HBITMAP)SelectObject(hLocalDC, hBitmap);
	BOOL qRetBlit = BitBlt(hdc, 0, 0, qBitmap.bmWidth, qBitmap.bmHeight, hLocalDC, 0, 0, SRCCOPY);
	SelectObject(hLocalDC, hOldBmp);
	DeleteDC(hLocalDC);
	DeleteObject(hBitmap);
}

// Function to generate a file name using date and time
wstring GenerateFileName()
{
	time_t currentTime = time(nullptr);
	tm localTime;
	localtime_s(&localTime, &currentTime);

	wchar_t buffer[400];
	wcsftime(buffer, sizeof(buffer), L"%Y-%m-%d_%H-%M-%S.bmp", &localTime);

	return wstring(buffer);
}

// Function to convert wstring to char*
char* ConvertWstringToChar(const wstring& wstr)
{
	int size = WideCharToMultiByte(CP_UTF8, 0, wstr.c_str(), -1, nullptr, 0, nullptr, nullptr);
	char* buffer = new char[size];
	WideCharToMultiByte(CP_UTF8, 0, wstr.c_str(), -1, buffer, size, nullptr, nullptr);
	return buffer;
}

void save(HWND& hWnd)
{
	HDC hdc = GetDC(hWnd);

	int windowWidth;
	int windowHeight;
	RECT rect;
	if (GetWindowRect(hWnd, &rect))
	{
		windowWidth = rect.right - rect.left;
		windowHeight = rect.bottom - rect.top;
	}
	RECT rect1 = { 0, 0, windowWidth, windowHeight };

	char* fileName = ConvertWstringToChar(GenerateFileName());

	HDCToFile((char*)fileName, hWnd, hdc, rect1);
}


// ------------------------------------------------------------------------------------------------------------------------


LRESULT WINAPI WndProc(HWND hWnd, UINT m, WPARAM wp, LPARAM lp)
{
    HDC hdc;
	static HINSTANCE hInstance;
    switch (m)
    {

    case WM_LBUTTONDOWN:
    {
        
        if (Selection == 5)
        {
            if (number_of_click==0)
            {
                hdc = GetDC(hWnd);
                xx1 = LOWORD(lp);
                yy1 = HIWORD(lp);
                number_of_click++;
            }
            else {
                hdc = GetDC(hWnd);
                xx2 = LOWORD(lp);
                yy2 = HIWORD(lp);
                LineDD(hdc, xx1, yy1, xx2, yy2, newColor);
                number_of_click = 0;
            }
        }
        else if (Selection == 6) {
            if (number_of_click == 0)
            {
                hdc = GetDC(hWnd);
                xx1 = LOWORD(lp);
                yy1 = HIWORD(lp);
                number_of_click++;
            }
            else {
                hdc = GetDC(hWnd);
                xx2 = LOWORD(lp);
                yy2 = HIWORD(lp);
                DrawlineB(hdc, xx1, yy1, xx2, yy2, newColor);
                number_of_click = 0;
            }
        }
		else if (Selection == 7) {
			if (number_of_click == 0)
			{
				hdc = GetDC(hWnd);
				xx1 = LOWORD(lp);
				yy1 = HIWORD(lp);
				number_of_click++;
			}
			else {
				hdc = GetDC(hWnd);
				xx2 = LOWORD(lp);
				yy2 = HIWORD(lp);
				DrawLine1(hdc, xx1, yy1, xx2, yy2, newColor);
				number_of_click = 0;
			}
		}
		else if (Selection == 8) {
			if (number_of_click == 0)
			{
				hdc = GetDC(hWnd);
				xx1 = LOWORD(lp);
				yy1 = HIWORD(lp);
				number_of_click++;
			}
			else {
				hdc = GetDC(hWnd);
				xx2 = LOWORD(lp);
				yy2 = HIWORD(lp);
				int R = sqrt(pow(xx2 - xx1, 2) + pow(yy2 - yy1, 2));
				CircleDirect(hdc, xx1, yy1, R, newColor);
				number_of_click = 0;
			}
		}
		else if (Selection == 9) {
			if (number_of_click == 0)
			{
				hdc = GetDC(hWnd);
				xx1 = LOWORD(lp);
				yy1 = HIWORD(lp);
				number_of_click++;
			}
			else {
				hdc = GetDC(hWnd);
				xx2 = LOWORD(lp);
				yy2 = HIWORD(lp);
				int R = sqrt(pow(xx2 - xx1, 2) + pow(yy2 - yy1, 2));
				DirectPolar(hdc, xx1, yy1, R, newColor);
				number_of_click = 0;
			}
		}
		else if (Selection == 10) {
			if (number_of_click == 0)
			{
				hdc = GetDC(hWnd);
				xx1 = LOWORD(lp);
				yy1 = HIWORD(lp);
				number_of_click++;
			}
			else {
				hdc = GetDC(hWnd);
				xx2 = LOWORD(lp);
				yy2 = HIWORD(lp);
				int R = sqrt(pow(xx2 - xx1, 2) + pow(yy2 - yy1, 2));
				CircleIterativePolar(hdc, xx1, yy1, R, newColor);
				number_of_click = 0;
			}
		}
		else if (Selection == 11) {
			if (number_of_click == 0)
			{
				hdc = GetDC(hWnd);
				xx1 = LOWORD(lp);
				yy1 = HIWORD(lp);
				number_of_click++;
			}
			else {
				hdc = GetDC(hWnd);
				xx2 = LOWORD(lp);
				yy2 = HIWORD(lp);
				int R = sqrt(pow(xx2 - xx1, 2) + pow(yy2 - yy1, 2));
				CircleBresenham(hdc, xx1, yy1, R, newColor);
				number_of_click = 0;
			}
		}
		else if (Selection == 12) {
			if (number_of_click == 0)
			{
				hdc = GetDC(hWnd);
				xx1 = LOWORD(lp);
				yy1 = HIWORD(lp);
				number_of_click++;
			}
			else {
				hdc = GetDC(hWnd);
				xx2 = LOWORD(lp);
				yy2 = HIWORD(lp);
				int R = sqrt(pow(xx2 - xx1, 2) + pow(yy2 - yy1, 2));
				CircleFasterBresenham(hdc, xx1, yy1, R, newColor);
				number_of_click = 0;
			}
		}
		else if (Selection == 13) {
			if (number_of_click == 0)
			{
				hdc = GetDC(hWnd);
				xx1 = LOWORD(lp);
				yy1 = HIWORD(lp);
				number_of_click++;
			}
			else {
				hdc = GetDC(hWnd);
				xx2 = LOWORD(lp);
				yy2 = HIWORD(lp);
				int R = sqrt(pow(xx2 - xx1, 2) + pow(yy2 - yy1, 2));
				Filling_Circle_lines(hdc, xx1, yy1, R, 3, newColor);
				number_of_click = 0;
			}

		}
		else if (Selection == 14) {
			if (number_of_click == 0)
			{
				
				xx1 = LOWORD(lp);
				yy1 = HIWORD(lp);
				number_of_click++;
			}
			else {
				hdc = GetDC(hWnd);
				xx2 = LOWORD(lp);
				yy2 = HIWORD(lp);
				int R = sqrt(pow(xx2 - xx1, 2) + pow(yy2 - yy1, 2));
				Filling_Circle_other_circles(hdc, xx1, yy1, R, 3, newColor);
				number_of_click = 0;
			}

		}
		else if (Selection == 15){
		// Filling Square with Hermit
		static vector<Vector2> v(2);
		if (number_of_click < 2)
		{
			v[number_of_click++] = Vector2(LOWORD(lp), HIWORD(lp));
		}else 
		{
			hdc = GetDC(hWnd);
			Rectangle(hdc, v[0].x, v[0].y, v[1].x, v[1].y);
			//FillSquareWithHermitVertical(hdc, Vector2(v[0].x, v[0].y),Vector2(v[1].x, v[1].y), newColor);
		}
		

        }
		else if (Selection == 16) { 
		// Filling Rectangle with Bezier
        }
		else if (Selection == 17) {
		if (number_of_click < pointSize -1)
		{
			points[number_of_click].x = LOWORD(lp);
			points[number_of_click].y = HIWORD(lp);
			number_of_click++;
		}
		else 
		{
			hdc = GetDC(hWnd);
			points[number_of_click].x = LOWORD(lp);
			points[number_of_click].y = HIWORD(lp);
			number_of_click = 0;
			ConvexPolygonFilling(hdc, points, pointSize, newColor);
		}

		}
		else if (Selection == 18) {
		if (number_of_click < pointSize - 1)
		{
			points[number_of_click].x = LOWORD(lp);
			points[number_of_click].y = HIWORD(lp);
			number_of_click++;
		}
		else
		{
			hdc = GetDC(hWnd);
			points[number_of_click].x = LOWORD(lp);
			points[number_of_click].y = HIWORD(lp);
			number_of_click = 0;
			GeneralPolygonFilling(hdc, points, pointSize, newColor);
		}

		}
		else if (Selection == 19)
		{
			hdc = GetDC(hWnd);
			xx1 = LOWORD(lp);
			yy1 = HIWORD(lp);
			fFill(hdc, xx1, yy1, lastColor, newColor);
		}
		else if (Selection == 20)
		{
			hdc = GetDC(hWnd);
			xx1 = LOWORD(lp);
			yy1 = HIWORD(lp);
			NRFloodFill(hdc, xx1, yy1,  lastColor , newColor);
		}
		else if (Selection == 21) {
		// cardinal spline
	    }
		else if (Selection == 22) {
		    // ellips direct
			if (number_of_click == 0)
			{
				hdc = GetDC(hWnd);
				xx1 = LOWORD(lp);
				yy1 = HIWORD(lp);
				number_of_click++;
			}
			else if (number_of_click == 1) {
				hdc = GetDC(hWnd);
				xx2 = LOWORD(lp);
				yy2 = HIWORD(lp);
				number_of_click++;
			}
			else if (number_of_click == 2) {
				hdc = GetDC(hWnd);
				xx3 = LOWORD(lp);
				yy3 = HIWORD(lp);
				number_of_click = 0;
				int A = sqrt(pow(xx2 - xx1, 2) + pow(yy2 - yy1, 2));
				int B = sqrt(pow(xx3 - xx1, 2) + pow(yy3 - yy1, 2));
				DirectEllipse(hdc, xx1, yy1, A, B, newColor);
			}
        }
		else if (Selection == 23) {
		// ellips poler
		if (number_of_click == 0)
		{
			hdc = GetDC(hWnd);
			xx1 = LOWORD(lp);
			yy1 = HIWORD(lp);
			number_of_click++;
		}
		else if (number_of_click == 1) {
			hdc = GetDC(hWnd);
			xx2 = LOWORD(lp);
			yy2 = HIWORD(lp);
			number_of_click++;
		}
		else if (number_of_click == 2) {
			hdc = GetDC(hWnd);
			xx3 = LOWORD(lp);
			yy3 = HIWORD(lp);
			number_of_click = 0;
			int A = sqrt(pow(xx2 - xx1, 2) + pow(yy2 - yy1, 2));
			int B = sqrt(pow(xx3 - xx1, 2) + pow(yy3 - yy1, 2));
			DrawEllipsePoler(hdc, xx1, yy1, A, B, newColor);
		}
		}
		else if (Selection == 24) {
		// ellips midpoint
		if (number_of_click == 0)
		{
			hdc = GetDC(hWnd);
			xx1 = LOWORD(lp);
			yy1 = HIWORD(lp);
			number_of_click++;
		}
		else if (number_of_click == 1) {
			hdc = GetDC(hWnd);
			xx2 = LOWORD(lp);
			yy2 = HIWORD(lp);
			number_of_click++;
		}
		else if (number_of_click == 2) {
			hdc = GetDC(hWnd);
			xx3 = LOWORD(lp);
			yy3 = HIWORD(lp);
			number_of_click = 0;
			int A = sqrt(pow(xx2 - xx1, 2) + pow(yy2 - yy1, 2));
			int B = sqrt(pow(xx3 - xx1, 2) + pow(yy3 - yy1, 2));
			DrawEllipseMidpoint(hdc, xx1, yy1, A, B, newColor);
		}
		}
		else if (Selection == 25) { 
        // "Clipping Rectangle as Clipping point
		static vector<Vector2> v(3);
		if (number_of_click < 2)
		{
			v[number_of_click++] = Vector2(LOWORD(lp), HIWORD(lp));
			if (number_of_click==2)
			{
				hdc = GetDC(hWnd);
				Rectangle(hdc, v[0].x, v[0].y, v[1].x, v[1].y);
			}
		}
		else
		{
			hdc = GetDC(hWnd);
			v[number_of_click++] = Vector2(LOWORD(lp), HIWORD(lp));
			ClipPoint(hdc, v[2].x, v[2].y, (int)v[0].x, (int)v[0].y, (int)v[1].x, (int)v[1].y, newColor);
			number_of_click = 2;

		}

        }
		else if (Selection == 26) {  
        // "Clipping Rectangle as Clipping line
        }
		else if (Selection == 27) {
		// "Clipping Rectangle as Clipping polygon
		}
		else if (Selection == 28) {
		// Clipping Square as Clipping point
		}
		else if (Selection == 29) {
		// Clipping Square as Clipping line
		}
		else if (Selection == 30) {
		//  Clipping circle as Clipping point
		}
		else if (Selection == 31) {
		// Clipping circle as Clipping line
		}
		else {
			MessageBox(hWnd, L"Please select function you want to do it from list.", L"Warning", MB_OK);
		}
		

        break;
    }
    case WM_COMMAND:
    {
        if (HIWORD(wp) == CBN_SELCHANGE && (HWND)lp== hCombo)
        {
            int index = SendMessage(hCombo, CB_GETCURSEL, 0, 0);
            if (index != CB_ERR)
            {
				// this code to get text of selection
                //wchar_t buffer[256];
                //SendMessage(hCombo, CB_GETLBTEXT, index, (LPARAM)buffer);
                // Do something with the selected text...
                // MessageBox(hWnd, buffer, L"Hello", MB_OK);

				if (index == 0) {
					
					// Change the background color to white
					SetClassLongPtr(hWnd, GCLP_HBRBACKGROUND, (LONG)CreateSolidBrush(RGB(255, 255, 255)));
					// Redraw the window to update the background color
					RedrawWindow(hWnd, NULL, NULL, RDW_ERASE | RDW_INVALIDATE);
				}
				else if(index == 1){
					// red color
					lastColor = newColor;
					newColor = red;
				}
				else if (index == 2) {
					// blue color 
					lastColor = newColor;
					newColor = blue;
				}
				else if (index == 3) {
					// green color 
					lastColor = newColor;
					newColor = green;
				}
				else if (index == 4) {
					// black color 
					lastColor = newColor;
					newColor = black;
				}
				else if (index == 32) {
					int result = MessageBox(hWnd, L"You want to clear window!?", L"comfirmation", MB_OKCANCEL);
					if (result == IDOK)
					{
						// User clicked OK

						hdc = GetDC(hWnd);
						COLORREF c = GetPixel(hdc, 201, 0);

						// Change the background color to blue
						SetClassLongPtr(hWnd, GCLP_HBRBACKGROUND, (LONG)CreateSolidBrush(c));
						// Redraw the window to update the background color
						RedrawWindow(hWnd, NULL, NULL, RDW_ERASE | RDW_INVALIDATE);

						// RESET SELECTION
						Selection = -1;

					}
					
				}
				else if (index == 33) {
					int result = MessageBox(hWnd, L"You want to save window!?", L"comfirmation", MB_OKCANCEL);
					if (result == IDOK)
					{
						// User clicked OK
						save(hWnd);

						// reset Selection
						Selection = -1;
			
					}

					
				}
				else if (index==34) {

					int result = MessageBox(hWnd, L"You want to load window from file!?", L"comfirmation", MB_OKCANCEL);
					if (result == IDOK)
					{
						// User clicked OK
					   hdc = GetDC(hWnd);
					   load(hWnd, hdc);

						// // reset Selection
						Selection = -1;

					}

					
				
				}
				else {
					Selection = index;
					number_of_click = 0;
				}


		 
            }
        }

        break;
       
    }
    case WM_CLOSE:
    {
        DestroyWindow(hWnd);
        break;
    }
    case WM_DESTROY:
    {
        PostQuitMessage(0);
        break;
    }
    default:
    {
        return DefWindowProc(hWnd, m, wp, lp);
    }
    }
    return 0;
    
}


    int APIENTRY WinMain(HINSTANCE h, HINSTANCE p, LPSTR cmd, int csh)
    {
		WNDCLASS wc;
		wc.cbClsExtra = wc.cbWndExtra = 0;
		wc.hbrBackground = (HBRUSH)GetStockObject(LTGRAY_BRUSH);
		wc.hCursor = LoadCursor(NULL, IDC_ARROW);
		wc.hIcon = LoadIcon(NULL, IDI_APPLICATION);
		wc.hInstance = h;
		wc.lpfnWndProc = WndProc;
		wc.lpszClassName = L"MyClass";
		wc.lpszMenuName = NULL;
		wc.style = CS_HREDRAW | CS_VREDRAW;
		RegisterClass(&wc);
        HWND hWnd = CreateWindow(L"MyClass", L"MMAAN", WS_OVERLAPPEDWINDOW, 250, 100, 1000, 700, NULL, NULL, h, 0);

        // Create a drop-down list with items
        hCombo = CreateWindow(L"COMBOBOX", L"", CBS_DROPDOWNLIST | WS_VISIBLE | WS_CHILD, 0, 0, 280,600, hWnd, NULL, h, NULL);
		// ---------------------------------------------------------------------------------------------/     index
		SendMessage(hCombo, CB_ADDSTRING, 0, (LPARAM)L"background white");                                   // 0   done
		SendMessage(hCombo, CB_ADDSTRING, 0, (LPARAM)L"color Red");                                          // 1   done
		SendMessage(hCombo, CB_ADDSTRING, 0, (LPARAM)L"color Blue");                                         // 2   done
		SendMessage(hCombo, CB_ADDSTRING, 0, (LPARAM)L"color Green");                                        // 3   done
		SendMessage(hCombo, CB_ADDSTRING, 0, (LPARAM)L"color Black");                                        // 4   done
        SendMessage(hCombo, CB_ADDSTRING, 0, (LPARAM)L"line DDA");                                           // 5   done
        SendMessage(hCombo, CB_ADDSTRING, 0, (LPARAM)L"line Midpoint");                                      // 6   done
		SendMessage(hCombo, CB_ADDSTRING, 0, (LPARAM)L"line parametric");                                    // 7   done
		SendMessage(hCombo, CB_ADDSTRING, 0, (LPARAM)L"Circle Direct");                                      // 8   done
		SendMessage(hCombo, CB_ADDSTRING, 0, (LPARAM)L"Circle Polar");                                       // 9   done
		SendMessage(hCombo, CB_ADDSTRING, 0, (LPARAM)L"Circle iterative Polar");                             // 10  done 
		SendMessage(hCombo, CB_ADDSTRING, 0, (LPARAM)L"Circle midpoint");                                    // 11  done
		SendMessage(hCombo, CB_ADDSTRING, 0, (LPARAM)L"Circle modified Midpoint");                           // 12  done
		SendMessage(hCombo, CB_ADDSTRING, 0, (LPARAM)L"Filling Circle with lines");                          // 13  done
		SendMessage(hCombo, CB_ADDSTRING, 0, (LPARAM)L"Filling Circle with other circles");                  // 14  done
		SendMessage(hCombo, CB_ADDSTRING, 0, (LPARAM)L"Filling Square with Hermit");                         // 15  done
		SendMessage(hCombo, CB_ADDSTRING, 0, (LPARAM)L"Filling Rectangle with Bezier");                      // 16  done
		SendMessage(hCombo, CB_ADDSTRING, 0, (LPARAM)L"Convex Filling");                                     // 17  done
		SendMessage(hCombo, CB_ADDSTRING, 0, (LPARAM)L"Non Convex Filling");                                 // 18  done
		SendMessage(hCombo, CB_ADDSTRING, 0, (LPARAM)L"FloodFill Recursive");                                // 19  done
		SendMessage(hCombo, CB_ADDSTRING, 0, (LPARAM)L"FloodFill Non Recursive");                            // 20  done
		SendMessage(hCombo, CB_ADDSTRING, 0, (LPARAM)L"Cardinal Spline Curve");                              // 21  done
		SendMessage(hCombo, CB_ADDSTRING, 0, (LPARAM)L"Ellipse Direct");                                     // 22
		SendMessage(hCombo, CB_ADDSTRING, 0, (LPARAM)L"Ellipse polar");                                      // 23  done
		SendMessage(hCombo, CB_ADDSTRING, 0, (LPARAM)L"Ellipse midpoint");                                   // 24  done
		SendMessage(hCombo, CB_ADDSTRING, 0, (LPARAM)L"Clipping Rectangle as Clipping point");               // 25
		SendMessage(hCombo, CB_ADDSTRING, 0, (LPARAM)L"Clipping Rectangle as Clipping line");                // 26
		SendMessage(hCombo, CB_ADDSTRING, 0, (LPARAM)L"Clipping Rectangle as Clipping Polygon");             // 27
		SendMessage(hCombo, CB_ADDSTRING, 0, (LPARAM)L"Clipping Square as Clipping point");                  // 28
		SendMessage(hCombo, CB_ADDSTRING, 0, (LPARAM)L"Clipping Square as Clipping line");                   // 29
		SendMessage(hCombo, CB_ADDSTRING, 0, (LPARAM)L"Clipping circle as Clipping point");                  // 30
		SendMessage(hCombo, CB_ADDSTRING, 0, (LPARAM)L"Clipping circle as Clipping line");                   // 31
        SendMessage(hCombo, CB_ADDSTRING, 0, (LPARAM)L"clear");                                              // 32  done
		SendMessage(hCombo, CB_ADDSTRING, 0, (LPARAM)L"save");                                               // 33  done
		SendMessage(hCombo, CB_ADDSTRING, 0, (LPARAM)L"load");                                               // 34  done

        ShowWindow(hWnd, csh);
        UpdateWindow(hWnd);
        MSG msg;
        while (GetMessage(&msg, NULL, 0, 0) > 0)
        {
            TranslateMessage(&msg);
            DispatchMessage(&msg);
        }

		return (int)msg.wParam;
    }












