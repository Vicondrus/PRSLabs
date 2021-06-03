// OpenCVApplication.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "common.h"
#include <vector>
#include <algorithm>
#include <random>

void choose(const int size, int& first, int& second)
{
	first = rand() % size;
	second = rand() % size;
	if (second == first)
	{
		++second;
	}
	if (second < first)
	{
		int aux = first;
		first = second;
		second = aux;
	}
}

std::vector<Point> findPoints(Mat img)
{
	std::vector<Point> pointsVector;
	for (int i = 0; i < img.rows; i++) 
	{
		for (int j = 0; j < img.cols; j++) 
		{
			if (img.at<uchar>(i,j) == 0)
			{
				pointsVector.push_back(Point(j, i));
			}
		}
	}
	return pointsVector;
}

std::vector<float> computeLineEquation(Point n, Point m)
{
	std::vector<float> paramVector;
	float a = n.y - m.y;
	float b = m.x - n.x;
	float c = n.x * m.y - m.x * n.y;
	paramVector.push_back(a);
	paramVector.push_back(b);
	paramVector.push_back(c);
	return paramVector;
}

std::vector<float> computeBetaAndRho(std::vector<Point> points)
{
	float beta, rho;

	float sumX = 0.f, sumY = 0.f, sumXY = 0.f, sumXSquared = 0.f, sumYSquared = 0.f;

	for (int i = 0; i < points.size(); i++) {
		sumX += points[i].x;
		sumY += points[i].y;
		sumXY += points[i].x * points[i].y;
		sumXSquared += points[i].x * points[i].x;
		sumYSquared += points[i].y * points[i].y;
	}

	beta = -atan2(2 * sumXY - (2 * sumX * sumY) / points.size(), sumYSquared - sumXSquared + (sumX * sumX) / points.size() - (sumY * sumY) / points.size()) / 2;
	rho = (cos(beta) * sumX + sin(beta) * sumY) / points.size();

	std::vector<float> params;
	params.push_back(beta);
	params.push_back(rho);
	return params;
}

float computeDistance(Point p, std::vector<float> paramVector)
{
	return fabs(p.x * paramVector[0] + p.y * paramVector[1] + paramVector[2]) / sqrt(paramVector[0] * paramVector[0] + paramVector[1] * paramVector[1]);
}

int main()
{
	int t = 10, s = 2;
	float p = 0.99, q = 0.8;
	printf("Input q: ");
	scanf("%f", &q);
	float T;
	int N;
	char fname[MAX_PATH];
	while (openFileDlg(fname)) 
	{
		Mat src = imread(fname, IMREAD_GRAYSCALE);
		std::vector<Point> pointsVector = findPoints(src);
		int n = pointsVector.size();
		T = q * n;
		N = (log(1 - p) / log(1 - pow(q, 2)) + 1);
		printf("N: %d, T: %f\n", N, T);
		Mat dst(src.rows, src.cols, CV_8UC1);
		int iter = 0;
		std::vector<float> lineVector;
		int consSetSize = 0;
		int maxConsSet = -1;
		std::vector<Point> consSetMax;
		do
		{
			consSetSize = 0;
			int pos1, pos2;
			choose(pointsVector.size(), pos1, pos2);
			Point p1 = pointsVector[pos1];
			Point p2 = pointsVector[pos2];
			std::vector<float> lv = computeLineEquation(p1, p2);
			std::vector<Point> consSet;
			for (Point p : pointsVector)
			{
				float d = computeDistance(p, lv);
				if (d < t)
				{
					consSetSize++;
					consSet.push_back(p);
				}
			}
			if (maxConsSet < consSetSize)
			{
				lineVector = lv;
				maxConsSet = consSetSize;
				consSetMax = consSet;
			}
			++iter;
		} while (iter < N && consSetSize < T);

		printf("Max Cons Size: %d\n", maxConsSet);

		std::vector<float> params = computeBetaAndRho(consSetMax);
		float beta = params[0];
		float rho = params[1];


		if (fabs(beta) < PI / 4 || fabs(beta) > 3 * PI / 4) {
			float x1 = 0, y1;
			y1 = (-cos(beta) * x1 + rho) / sin(beta);
			float x2 = src.cols - 1, y2;
			y2 = (-cos(beta) * x2 + rho) / sin(beta);
			line(src, Point(x1, y1), Point(x2, y2), Scalar(0, 255, 0));
		}
		else {
			float x1, y1 = 0;
			x1 = (-sin(beta) * y1 + rho) / cos(beta);
			float x2, y2 = src.cols - 1;
			x2 = (-sin(beta) * y2 + rho) / cos(beta);
			line(src, Point(x1, y1), Point(x2, y2), Scalar(0, 255, 0));
		}

		imshow("image", src);
		waitKey();
	}
	return 0;
}