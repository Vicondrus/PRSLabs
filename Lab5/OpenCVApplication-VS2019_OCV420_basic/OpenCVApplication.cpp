#include "stdafx.h"
#include "common.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <random>

void methodOne() {
	char fname[MAX_PATH];
	while (openFileDlg(fname))
	{
		FILE* f = fopen(fname, "r");
		int n;
		fscanf(f, "%d", &n);
		std::vector<Point> pointsVector;
		float xmax = -10000, ymin = -10000;
		for (int i = 0; i < n; i++) {
			float x, y;
			fscanf(f, "%f%f", &x, &y);
			pointsVector.push_back(Point(x, y));
		}
		fclose(f);
		Mat img(500, 500, CV_8UC3);
		for (int i = 0; i < n; i++) {
			if (pointsVector[i].x > 0 && pointsVector[i].y > 0) {
				circle(img, pointsVector[i], 1, Scalar(0, 0, 0), FILLED, 1, 0);
			}
		}

		float theta0, theta1;

		float sumX = 0.f, sumY = 0.f, sumXY = 0.f, sumXSquared = 0.f;

		for (int i = 0; i < n; i++) {
			sumX += pointsVector[i].x;
			sumY += pointsVector[i].y;
			sumXY += pointsVector[i].x * pointsVector[i].y;
			sumXSquared += pointsVector[i].x * pointsVector[i].x;
		}

		theta1 = (n * sumXY - sumX * sumY) / (n * sumXSquared - sumX * sumX);
		theta0 = (sumY - theta1 * sumX) / n;

		if (fabs(theta1) < 1) {
			int x1 = 0, y1;
			y1 = theta1 * x1 + theta0;
			int x2 = img.cols - 1, y2;
			y2 = theta1 * x2 + theta0;
			line(img, Point(x1, y1), Point(x2, y2), Scalar(0, 0, 255));
		}
		else {
			int x1, y1 = 0;
			x1 = (y1 - theta0) / theta1;
			int x2, y2 = img.cols - 1;
			x2 = (y2 - theta0) / theta1;
			line(img, Point(x1, y1), Point(x2, y2), Scalar(0, 0, 255));
		}

		printf("%f %f\n", theta0, theta1);

		imshow("image", img);
		waitKey();
	}
}

void methodTwo() {
	char fname[MAX_PATH];
	while (openFileDlg(fname))
	{
		FILE* f = fopen(fname, "r");
		int n;
		fscanf(f, "%d", &n);
		std::vector<Point> pointsVector;
		float xmax = -10000, ymin = -10000;
		for (int i = 0; i < n; i++) {
			float x, y;
			fscanf(f, "%f%f", &x, &y);
			pointsVector.push_back(Point(x, y));
		}
		fclose(f);
		Mat img(500, 500, CV_8UC3);
		for (int i = 0; i < n; i++) {
			if (pointsVector[i].x > 0 && pointsVector[i].y > 0) {
				circle(img, pointsVector[i], 1, Scalar(0, 0, 0), FILLED, 1, 0);
			}
		}

		float beta, rho;

		float sumX = 0.f, sumY = 0.f, sumXY = 0.f, sumXSquared = 0.f, sumYSquared = 0.f;

		for (int i = 0; i < n; i++) {
			sumX += pointsVector[i].x;
			sumY += pointsVector[i].y;
			sumXY += pointsVector[i].x * pointsVector[i].y;
			sumXSquared += pointsVector[i].x * pointsVector[i].x;
			sumYSquared += pointsVector[i].y * pointsVector[i].y;
		}

		beta = -atan2(2 * sumXY - (2 * sumX * sumY) / n, sumYSquared - sumXSquared + (sumX * sumX) / n - (sumY * sumY) / n) / 2;
		rho = (cos(beta) * sumX + sin(beta) * sumY) / n;

		if (fabs(beta) < PI / 4 || fabs(beta) > 3 * PI / 4) {
			float x1 = 0, y1;
			y1 = (-cos(beta) * x1 + rho) / sin(beta);
			float x2 = img.cols - 1, y2;
			y2 = (-cos(beta) * x2 + rho) / sin(beta);
			line(img, Point(x1, y1), Point(x2, y2), Scalar(0, 255, 0));
		}
		else {
			float x1, y1 = 0;
			x1 = (-sin(beta) * y1 + rho) / cos(beta);
			float x2, y2 = img.cols - 1;
			x2 = (-sin(beta) * y2 + rho) / cos(beta);
			line(img, Point(x1, y1), Point(x2, y2), Scalar(0, 255, 0));
		}

		imshow("image", img);
		waitKey();
	}
}

void methodThree() {
	char fname[MAX_PATH];
	while (openFileDlg(fname))
	{
		FILE* f = fopen(fname, "r");
		int n;
		fscanf(f, "%d", &n);
		std::vector<Point> pointsVector;
		float xmax = -10000, ymin = -10000;
		for (int i = 0; i < n; i++) {
			float x, y;
			fscanf(f, "%f%f", &x, &y);
			pointsVector.push_back(Point(x, y));
		}
		fclose(f);
		Mat img(500, 500, CV_8UC3);
		for (int i = 0; i < n; i++) {
			if (pointsVector[i].x > 0 && pointsVector[i].y > 0) {
				circle(img, pointsVector[i], 1, Scalar(0, 0, 0), FILLED, 1, 0);
			}
		}

		float theta0 = 300, theta1 = -3;
		float alpha = 0.0000001;
		float eps;

		do {
			Mat img(500, 500, CV_8UC3);
			for (int i = 0; i < n; i++) {
				if (pointsVector[i].x > 0 && pointsVector[i].y > 0) {
					circle(img, pointsVector[i], 1, Scalar(0, 0, 0), FILLED, 1, 0);
				}
			}

			float theta0Diff = 0, theta1Diff = 0;


			for (int i = 0; i < n; i++) {
				float diff = theta0 + theta1 * pointsVector[i].x - pointsVector[i].y;
				theta0Diff += diff;
				theta1Diff += diff * pointsVector[i].x;
			}

			theta0 -= alpha * 100 * theta0Diff;
			theta1 -= alpha * theta1Diff;

			eps = (abs(theta0Diff) + abs(theta1Diff));


			if (fabs(theta1) < 1) {
				float x1 = 0, y1;
				y1 = theta1 * x1 + theta0;
				float x2 = img.cols - 1, y2;
				y2 = theta1 * x2 + theta0;
				line(img, Point(x1, y1), Point(x2, y2), Scalar(0, 0, 255));
			}
			else {
				float x1, y1 = 0;
				x1 = (y1 - theta0) / theta1;
				float x2, y2 = img.cols - 1;
				x2 = (y2 - theta0) / theta1;
				line(img, Point(x1, y1), Point(x2, y2), Scalar(0, 0, 255));
			}

			imshow("image", img);
			waitKey(1);

		} while (eps > 0.0001);

		imshow("image", img);
		waitKey();
	}
}

void methodFour() {
	char fname[MAX_PATH];
	while (openFileDlg(fname))
	{
		FILE* f = fopen(fname, "r");
		int n;
		fscanf(f, "%d", &n);
		std::vector<Point> pointsVector;
		float xmax = -10000, ymin = -10000;
		for (int i = 0; i < n; i++) {
			float x, y;
			fscanf(f, "%f%f", &x, &y);
			pointsVector.push_back(Point(x, y));
		}
		fclose(f);
		Mat img(500, 500, CV_8UC3);
		for (int i = 0; i < n; i++) {
			if (pointsVector[i].x > 0 && pointsVector[i].y > 0) {
				circle(img, pointsVector[i], 1, Scalar(0, 0, 0), FILLED, 1, 0);
			}
		}

		float beta = 1, rho = 1;
		float alpha = 0.00000001;
		float eps;

		do {
			Mat img(500, 500, CV_8UC3);
			for (int i = 0; i < n; i++) {
				if (pointsVector[i].x > 0 && pointsVector[i].y > 0) {
					circle(img, pointsVector[i], 1, Scalar(0, 0, 0), FILLED, 1, 0);
				}
			}

			float betaDiff = 0, rhoDiff = 0;


			for (int i = 0; i < n; i++) {
				float diff = pointsVector[i].x * cos(beta) + pointsVector[i].y * sin(beta) - rho;
				betaDiff += diff * (-pointsVector[i].x * sin(beta) + pointsVector[i].y * cos(beta));
				rhoDiff += -diff;
			}

			beta -= alpha * betaDiff;
			rho -= alpha * 100 * rhoDiff;

			eps = (abs(beta) + abs(rho));


			if (fabs(beta) < PI / 4 || fabs(beta) > 3 * PI / 4) {
				float x1 = 0, y1;
				y1 = (-cos(beta) * x1 + rho) / sin(beta);
				float x2 = img.cols - 1, y2;
				y2 = (-cos(beta) * x2 + rho) / sin(beta);
				line(img, Point(x1, y1), Point(x2, y2), Scalar(0, 255, 0));
			}
			else {
				float x1, y1 = 0;
				x1 = (-sin(beta) * y1 + rho) / cos(beta);
				float x2, y2 = img.cols - 1;
				x2 = (-sin(beta) * y2 + rho) / cos(beta);
				line(img, Point(x1, y1), Point(x2, y2), Scalar(0, 255, 0));
			}

			imshow("image", img);
			waitKey(1);

		} while (eps > 0.0001);

		imshow("image", img);
		waitKey();
	}
}

void lab1() {
	printf("Choose method:\n");
	int op;
	do
	{
		system("cls");
		destroyAllWindows();
		printf("Menu:\n");
		printf(" 1 - Method One\n");
		printf(" 2 - Method Two\n");
		printf(" 3 - Method One - Gradient Descent\n");
		printf(" 4 - Method Two - Gradient Descent\n");
		printf(" 0 - Exit\n\n");
		printf("Option: ");
		scanf("%d", &op);
		switch (op)
		{
		case 1:
			methodOne();
			break;
		case 2:
			methodTwo();
			break;
		case 3:
			methodThree();
			break;
		case 4:
			methodFour();
			break;
		}
	} while (op != 0);
}

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
			if (img.at<uchar>(i, j) == 0)
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

void ransac() {
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
}

void hough(int n, int k) {
	Mat src;
	Mat res;
	char fname[MAX_PATH];
	while (openFileDlg(fname))
	{
		src = imread(fname, IMREAD_GRAYSCALE);
		char fname2[MAX_PATH];
		openFileDlg(fname2);
		res = imread(fname2, IMREAD_COLOR);
		int d = sqrt(src.cols * src.cols + src.rows * src.rows);
		Mat hough(d + 1, 360, CV_32SC1);
		hough.setTo(0);
		int theta = 0, rho = 0;
		float beta = 0;
		int maxx = 0;
		int thetaStep = 1;

		for (int i = 0; i < src.rows; i++) {
			for (int j = 0; j < src.cols; j++) {
				if (src.at<uchar>(i, j) == 255) {
					for (theta = 0; theta < 360; theta += thetaStep) {

						beta = (float)(theta * PI) / 180.0f;
						float x = j * cos(beta) + i * sin(beta);

						if (x >= 0 && x <= d) {
							rho = (int)x;
							hough.at<int>(rho, theta)++;
							if (hough.at<int>(rho, theta) > maxx) {
								maxx = hough.at<int>(rho, theta);
							}
						}
					}
				}
			}
		}

		Mat houghImg;
		hough.convertTo(houghImg, CV_8UC1, 255.f / maxx);
		imshow("houghspace", houghImg);

		struct peak {
			int theta, rho, hval;
			bool operator < (const peak& o) const {
				return hval > o.hval;
			}
		};

		std::vector<peak> peaks;

		for (int i = 0; i < d + 1; i++) {
			for (int j = 0; j < 360; j++) {
				peak p;
				bool maxim = true;
				for (int k = -n / 2; k <= n / 2; k++) {
					for (int l = -n / 2; l <= n / 2; l++) {
						if (hough.at<int>(i, j) < hough.at<int>(min(max(i + k, 0), d), min(max(j + l, 0), 360 - 1))) {
							maxim = false;
						}
					}
				}
				if (maxim == true)
				{
					p.theta = j;
					p.rho = i;
					p.hval = hough.at<int>(i, j);

					peaks.push_back(p);
				}
			}
		}

		std::sort(peaks.begin(), peaks.end());
		Point p1, p2;
		for (int i = 0; i < k; i++) {
			rho = peaks.at(i).rho;
			theta = peaks.at(i).theta;
			beta = (float)(theta * PI) / 180.0f;
			if (fabs(beta) < PI / 4 || fabs(beta) > 3 * PI / 4) {
				float x1 = 0, y1;
				y1 = (-cos(beta) * x1 + rho) / sin(beta);
				float x2 = res.cols - 1, y2;
				y2 = (-cos(beta) * x2 + rho) / sin(beta);
				line(res, Point(x1, y1), Point(x2, y2), Scalar(0, 255, 0));
			}
			else {
				float x1, y1 = 0;
				x1 = (-sin(beta) * y1 + rho) / cos(beta);
				float x2, y2 = res.cols - 1;
				x2 = (-sin(beta) * y2 + rho) / cos(beta);
				line(res, Point(x1, y1), Point(x2, y2), Scalar(0, 255, 0));
			}
		}

		imshow("result", res);
		waitKey(0);
	}

}

void houghTransform() {
	int k, n;
	printf("Input window size: ");
	scanf("%d", &n);
	printf("Input number of local maxima: ");
	scanf("%d", &k);
	hough(n, k);
}

Mat chamferDT(Mat img) {
	int di[9] = { -1, -1, -1, 0, 0, 0, 1, 1, 1 };
	int dj[9] = { -1, 0, 1, -1, 0, 1, -1, 0, 1 };
	int w[9] = { 7, 5, 7, 5, 0, 5, 7, 5, 7 };

	Mat src = img.clone();

	for (int i = 1; i < src.rows; i++) {
		for (int j = 1; j < src.cols - 1; j++) {
			uchar min = src.at<uchar>(i, j);
			for (int k = 0; k < 5; k++) {
				if (min > (src.at<uchar>(i + di[k], j + dj[k]) + w[k])) {
					min = src.at<uchar>(i + di[k], j + dj[k]) + w[k];
				}
			}
			src.at<uchar>(i, j) = min;
		}
	}

	for (int i = src.rows - 2; i >= 0; i--) {
		for (int j = src.cols - 2; j > 0; j--) {
			uchar min = src.at<uchar>(i, j);
			for (int k = 4; k < 9; k++) {
				if (min > (src.at<uchar>(i + di[k], j + dj[k]) + w[k])) {
					min = src.at<uchar>(i + di[k], j + dj[k]) + w[k];
				}
			}
			src.at<uchar>(i, j) = min;
		}
	}

	return src;
}

int computeArea(Mat src) {
	int area = 0;
	for (int i = 0; i < src.rows; i++) {
		for (int j = 0; j < src.cols; j++) {
			if (src.at<uchar>(i, j) == 0)
				area++;
		}
	}
	return area;
}

Point computeCenterOfMass(Mat src) {
	int area = computeArea(src);
	int sumr = 0;
	int sumc = 0;
	for (int i = 0; i < src.rows; i++) {
		for (int j = 0; j < src.cols; j++) {
			if (src.at<uchar>(i, j) == 0) {
				sumr += i;
				sumc += j;
			}
		}
	}
	Point result(sumr / area, sumc / area);
	return result;
}

float patternMatching(Mat pattern, Mat unknown) {
	Point center1 = computeCenterOfMass(pattern);
	Point center2 = computeCenterOfMass(unknown);
	Mat trans = chamferDT(pattern);
	int sum = 0;
	int count = 0;
	for (int i = 0; i < unknown.rows; i++) {
		for (int j = 0; j < unknown.cols; j++) {
			if (unknown.at<uchar>(i, j) == 0) {
				if (i - center1.x + center2.x >= trans.rows || j - center1.y + center2.y >= trans.cols ||
					i - center1.x + center2.x < 0 || j - center1.y + center2.y < 0)
					continue;
				sum += trans.at<uchar>(i - center1.x + center2.x, j - center1.y + center2.y);
				count++;
			}
		}
	}
	float avg = (float)sum / (float)count;
	return avg;
}

void lab4()
{
	printf("Chamfer transform (1) or pattern matching (2)?\n");
	int choice;
	scanf("%d", &choice);
	char fname[MAX_PATH];
	switch (choice) {
	case 1:

		while (openFileDlg(fname))
		{
			Mat src = imread(fname, IMREAD_GRAYSCALE);
			Mat trans = chamferDT(src);
			imshow("src", src);
			imshow("trans", trans);
			waitKey(0);
		}
		break;
	case 2:
		char fname2[MAX_PATH];
		while (openFileDlg(fname))
		{
			Mat src = imread(fname, IMREAD_GRAYSCALE);
			openFileDlg(fname2);
			Mat ukn = imread(fname2, IMREAD_GRAYSCALE);
			float score = patternMatching(src, ukn);
			printf("%f\n", score);
			waitKey(0);
		}
		break;
	default:
		break;
	}
}

Mat computeFeatureMatrix(int noImages, int imageW, int imageH) {
	char folder[256] = "faces";
	char fname[256];
	Mat I = Mat(noImages, imageW * imageH, CV_8UC1);
	for (int i = 1; i <= noImages; i++) {
		sprintf(fname, "%s/face%05d.bmp", folder, i);
		Mat img = imread(fname, 0);
		for (int j = 0; j < img.rows; j++) {
			for (int k = 0; k < img.cols; k++) {
				I.at<uchar>(i - 1, j * img.cols + k) = img.at<uchar>(j, k);
			}
		}
	}
	return I;
}

Mat computeMeans(Mat I, int noImages, int imageW, int imageH) {
	Mat means = Mat(1, imageW * imageH, CV_32FC1);
	for (int i = 0; i < I.cols; i++) {
		int sum = 0;
		for (int k = 0; k < I.rows; k++) {
			sum += I.at<uchar>(k, i);
		}
		means.at<float>(0, i) = (float)sum / noImages;
	}
	return means;
}

void saveMatToCsv(Mat mat, char csvName[]) {
	std::ofstream csvFile;
	csvFile.open(csvName);
	for (int i = 0; i < mat.rows; i++) {
		for (int j = 0; j < mat.cols; j++) {
			csvFile << mat.at<float>(i, j);
			csvFile << ", ";
		}
		csvFile << "\n";
	}
	csvFile.close();
}

Mat computeCovariance(Mat means, Mat I, int noImages, int imageW, int imageH) {
	Mat covariance = Mat(imageW * imageH, imageW * imageH, CV_32FC1);
	for (int i = 0; i < covariance.rows; i++) {

		for (int j = 0; j < covariance.cols; j++) {
			int sum = 0;
			for (int k = 0; k < noImages; k++) {
				sum += (I.at<uchar>(k, i) - means.at<float>(0,i)) * (I.at<uchar>(k, j) - means.at<float>(0, j));
			}
			covariance.at<float>(i, j) = (float)sum / noImages;
		}
	}
	return covariance;
}

Mat computeStdDev(Mat means, Mat I, int noImages, int imageW, int imageH) {
	Mat stdDev = Mat(1, imageW * imageH, CV_32FC1);
	for (int i = 0; i < I.cols; i++) {
		int sum = 0;
		for (int j = 0; j < I.rows; j++) {
			sum += (I.at<uchar>(j, i) - means.at<float>(0, i)) * (I.at<uchar>(j, i) - means.at<float>(0, i));
		}
		stdDev.at<float>(0, i) = sqrt((float)sum / noImages);
	}
	return stdDev;
}

Mat computeCorrelation(Mat stdDev, Mat covariance, int noImages, int imageW, int imageH) {
	Mat correlation = Mat(361, 361, CV_32FC1);
	for (int i = 0; i < covariance.rows; i++) {

		for (int j = 0; j < covariance.cols; j++) {
			correlation.at<float>(i, j) = covariance.at<float>(i, j) / (stdDev.at<float>(0, i) * stdDev.at<float>(0, j));
		}
	}
	return correlation;
}

Mat showCorrelationChart(Mat I, Mat correlation, int noImage, int imageW, int imageH, int posX1, int posY1, int posX2, int posY2) {
	Mat correlationChart = Mat(256, 256, CV_8UC1);
	for (int i = 0; i < 256; i++)
		for (int j = 0; j < 256; j++)
			correlationChart.at<uchar>(i, j) = 255;
	for (int i = 0; i < noImage; i++) {
		correlationChart.at<uchar>(I.at<uchar>(i, posY1 * imageW + posX1), I.at<uchar>(i, posY2 * imageW + posX2)) = 0;
	}
	printf("The correlation between %d and %d is of %f\n", posY1 * imageW + posX1, posY2 * imageW + posX2, correlation.at<float>(posY1 * imageW + posX1, posY2 * imageW + posX2));
	imshow("chart", correlationChart);
	waitKey(0);
	return correlationChart;
}

Mat drawGaussian(int feature, Mat stdDev, Mat means) {
	Mat gaussianChart = Mat(256, 256, CV_8UC1);
	for (int i = 0; i < 256; i++)
		for (int j = 0; j < 256; j++)
			gaussianChart.at<uchar>(i, j) = 255;
	for (int i = 0; i < 255; i++) {
		int j = 255 * exp(-(i - means.at<float>(0, feature)) * (i - means.at<float>(0, feature)) / (2 * stdDev.at<float>(0, feature) * stdDev.at<float>(0, feature)));
		gaussianChart.at<uchar>(255 - j, i) = 0;
	}
	imshow("gaussian", gaussianChart);
	waitKey(0);
	return gaussianChart;
}

int main()
{
	printf("Pick the laboratory number (1-5):\n");
	int choice;
	scanf("%d", &choice);
	switch (choice) {
	case 1:
		lab1();
		break;
	case 2:
		ransac();
		break;
	case 3:
		houghTransform();
		break;
	case 4:
		lab4();
		break;
	case 5:
		Mat I = computeFeatureMatrix(400, 19, 19);
		Mat means = computeMeans(I, 400, 19, 19);
		Mat covariance = computeCovariance(means, I, 400, 19, 19);
		Mat stdDev = computeStdDev(means, I, 400, 19, 19);
		Mat correlation = computeCorrelation(stdDev, covariance, 400, 19, 19);
		saveMatToCsv(means, "means.csv");
		saveMatToCsv(covariance, "covariance.csv");
		saveMatToCsv(correlation, "correlation.csv");
		showCorrelationChart(I, correlation, 400, 19, 19, 4, 5, 14, 5);
		showCorrelationChart(I, correlation, 400, 19, 19, 3, 10, 15, 9);
		showCorrelationChart(I, correlation, 400, 19, 19, 4, 5, 0, 18);
		int feature;
		do {
			scanf("%d", &feature);
			drawGaussian(feature, stdDev, means);
		} while (feature != -1);
		break;
	}
	return 0;
}