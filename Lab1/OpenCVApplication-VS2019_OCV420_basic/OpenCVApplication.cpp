// OpenCVApplication.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "common.h"
#include <vector>

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

		float theta0 = 250, theta1 = 1;
		float alpha = 0.000001;
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

			theta0 -= alpha * theta0Diff;
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
		float alpha = 0.0000001;
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
			rho -= alpha * rhoDiff;

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

int main()
{
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
	return 0;
}