#include "stdafx.h"
#include "common.h"
#include <vector>
#include <algorithm>
#include <random>

int getClass(Vec3b x) {
	if (x == Vec3b(0, 0, 255)) {
		return 1;
	}
	if (x == Vec3b(255, 0, 0)) {
		return -1;
	}
	return 0;
}

void draw(Mat src, std::vector<float> w, int wait = 0) {
	float theta0 = - w[0] / w[2], theta1 = - w[1] / w[2];
	Mat clone = src.clone();
	if (fabs(theta1) < 1) {
		int x1 = 0, y1;
		y1 = theta1 * x1 + theta0;
		int x2 = clone.cols - 1, y2;
		y2 = theta1 * x2 + theta0;
		line(clone, Point(x1, y1), Point(x2, y2), Scalar(0, 255, 0));
	}
	else {
		int x1, y1 = 0;
		x1 = (y1 - theta0) / theta1;
		int x2, y2 = clone.cols - 1;
		x2 = (y2 - theta0) / theta1;
		line(clone, Point(x1, y1), Point(x2, y2), Scalar(0, 255, 0));
	}
	imshow("line", clone);
	waitKey(wait);
}

std::pair<Mat, Mat> getTrainingSet(Mat src) {
	Mat X(0, 3, CV_32SC1);
	Mat y(0, 1, CV_32SC1);
	for (int i = 0; i < src.rows; i++) {
		for (int j = 0; j < src.cols; j++) {
			int currClass = getClass(src.at<Vec3b>(i, j));
			if (currClass != 0) {
				Mat aux1 = (Mat_<int>(1, 3) << 1, j, i);
				X.push_back(aux1);
				Mat aux2(1, 1, CV_32SC1, currClass);
				y.push_back(aux2);
			}
		}
	}
	for (int i = 0; i < X.rows; i++) {
		printf("%d %d %d		%d\n", X.at<int>(i, 0), X.at<int>(i, 1), X.at<int>(i, 2), y.at<int>(i, 0));
	}
	printf("\n");
	return std::pair<Mat, Mat>(X, y);
}

std::vector<float> onlineLearning(Mat X, Mat y, int max_it = 100000, double e_limit = 0.00001, double eta = 0.0001, std::vector<float> w = {0.1, 0.1, 0.1}) {
	for (int i = 0; i < max_it; i++) {
		float e = 0;
		for (int j = 0; j < X.rows; j++) {
			float z = 0.0;
			for (int k = 0; k < X.cols; k++) {
				z += w[k] * X.at<int>(j, k);
			}
			if (z * y.at<int>(j, 0) <= 0) {
				for (int k = 0; k < X.cols; k++)
				{
					w[k] = w[k] + eta * X.at<int>(j, k) * y.at<int>(j, 0);
				}
				e++;
			}
		}
		e = e / X.rows;
		if (e < e_limit) {
			break;
		}
	}
	return w;
}

std::vector<float> batchLearning(Mat X, Mat y, Mat src, int max_it = 100000, double e_limit = 0.00001, double eta = 0.00001, std::vector<float> w = { 0.1, 0.1, 0.1 }) {
	for (int iter = 0; iter < max_it; iter++) {
		float e = 0, l = 0;
		std::vector<float> grad = { 0, 0, 0 };
		for (int i = 0; i < X.rows; i++) {
			float z = 0.0;
			for (int j = 0; j < X.cols; j++) {
				z += w[j] * X.at<int>(i, j);
			}
			if (z * y.at<int>(i, 0) <= 0) {
				for (int j = 0; j < X.cols; j++) {
					grad[j] -= X.at<int>(i, j) * y.at<int>(i, 0);
				}
				e++;
				l -= z * y.at<int>(i, 0);
			}
		}
		e = e / X.rows;
		l = l / X.rows;
		printf("\rLoss: %f", l);
		for (int j = 0; j < X.cols; j++) {
			grad[j] = grad[j] / X.rows;
		}
		if (e < e_limit) {
			break;
		}
		for (int j = 0; j < X.cols; j++) {
			w[j] -= eta * grad[j];
		}
		draw(src, w, 1);
	}
	return w;
}

int main()
{
	Mat src;
	char fname[MAX_PATH];
	while (openFileDlg(fname)) {
		src = imread(fname, IMREAD_COLOR);
		std::pair<Mat, Mat> Xy = getTrainingSet(src);
		printf("Choose online (1) or batch (2) learning\nChoice: ");
		int choice = 0;
		do {
			scanf("%d", &choice);
		} while (choice != 1 && choice != 2);
		if (choice == 1) {
			std::vector<float> w = onlineLearning(Xy.first, Xy.second);
			draw(src, w);
		}
		if (choice == 2) {
			std::vector<float> w = batchLearning(Xy.first, Xy.second, src);
			draw(src, w);
		}
	}
	return 0;
}