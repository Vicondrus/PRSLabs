#include "stdafx.h"
#include "common.h"
#include <vector>
#include <algorithm>
#include <random>

Mat readFile(char fname[MAX_PATH]) {
	int n, d;
	FILE* fp;
	fp = fopen(fname, "r");
	fscanf(fp, "%d %d", &n, &d);
	Mat X(n, d, CV_64FC1);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < d; j++) {
			double x;
			fscanf(fp, "%lf", &x);
			X.at<double>(i, j) = x;
		}
	}
	return X;
}

void printMat(Mat X) {
	for (int i = 0; i < X.rows; i++) {
		for (int j = 0; j < X.cols; j++) {
			printf("%lf ", X.at<double>(i, j));
		}
		printf("\n");
	}
}

std::vector<Mat> getPcaAndApproximate(Mat Q, Mat X, int k) {
	Mat Qk(X.cols, k, CV_64FC1);
	for (int i = 0; i < X.cols; i++) {
		for (int j = 0; j < k; j++) {
			Qk.at<double>(i, j) = Q.at<double>(i, j);
		}
	}
	
	Mat Xcoef(X.cols, k, CV_64FC1);
	Xcoef = X * Q;

	Mat Xk(X.rows, X.cols, CV_64FC1);
	Xk = X * Qk * Qk.t();

	std::vector<Mat> result;
	result.push_back(Xcoef);
	result.push_back(Xk);
	return result;
}

Mat subtractMeanVector(Mat X) {
	int n = X.rows;
	int d = X.cols;
	std::vector<double> means(d);
	for (int i = 0; i < d; i++) {
		double sum = 0;
		for (int j = 0; j < n; j++) {
			sum += X.at<double>(j, i);
		}
		means.at(i) = sum / n;
	}

	Mat X2(n, d, CV_64FC1);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < d; j++) {
			X2.at<double>(i, j) = X.at<double>(i, j) - means.at(j);
		}
	}

	return X2;
}

double getMeanAbsoluteDifference(Mat X, Mat Xk) {
	double sum = 0;
	for (int i = 0; i < X.rows; i++) {
		for (int j = 0; j < X.cols; j++) {
			sum += abs(X.at<double>(i, j) - Xk.at<double>(i, j));
		}
	}

	return sum / (X.cols * X.rows);
}

void printEigenvaluesAndMad(Mat Lambda, Mat X, Mat Q) {
	std::vector<Mat> aux = getPcaAndApproximate(Q, X, 1);
	printf("Eigenvalues:\n");
	for (int i = 0; i < Lambda.rows; i++) {
		printf("%f\n", Lambda.at<double>(i, 0));
	}
	double mad = getMeanAbsoluteDifference(X, aux[1]);
	printf("Mean absolute difference: %lf\n", mad);
}

std::vector<std::vector<double>> getMinAndMax(Mat Xcoef) {
	std::vector<double> mins;
	std::vector<double> maxs;
	for (int i = 0; i < Xcoef.cols; i++) {
		double min = FLT_MAX;
		double max = -FLT_MAX;
		for (int j = 0; j < Xcoef.rows; j++) {
			double value = Xcoef.at<double>(j, i);
			if (min > value) {
				min = value;
			}
			else if (max < value) {
				max = value;
			}
		}
		mins.push_back(min);
		maxs.push_back(max);
	}
	std::vector<std::vector<double>> result;
	result.push_back(mins);
	result.push_back(maxs);
	return result;
}

void showPca2D(Mat Xcoef) {
	std::vector<std::vector<double>> minsAndMaxs = getMinAndMax(Xcoef);
	std::vector<double> mins = minsAndMaxs[0];
	std::vector<double> maxs = minsAndMaxs[1];

	Mat res(maxs[0] - mins[0] + 1, maxs[1] - mins[1] + 1, CV_8UC1);

	for (int i = 0; i < res.rows; i++) {
		for (int j = 0; j < res.cols; j++) {
			res.at<uchar>(i, j) = 255;
		}
	}

	for (int i = 0; i < Xcoef.rows; i++) {
		res.at<uchar>(Xcoef.at<double>(i, 1) - mins[1], Xcoef.at<double>(i, 0) - mins[0]) = 0;
	}

	imshow("2D", res);
	waitKey();
}

void showPca3D(Mat Xcoef) {
	std::vector<std::vector<double>> minsAndMaxs = getMinAndMax(Xcoef);
	std::vector<double> mins = minsAndMaxs[0];
	std::vector<double> maxs = minsAndMaxs[1];

	Mat res(maxs[0] - mins[0] + 1, maxs[1] - mins[1] + 1, CV_8UC1);

	for (int i = 0; i < res.rows; i++) {
		for (int j = 0; j < res.cols; j++) {
			res.at<uchar>(i, j) = 255;
		}
	}

	for (int i = 0; i < Xcoef.rows; i++) {
		uchar color = (Xcoef.at<double>(i, 2) - mins[2]) * 256 / (maxs[2] - mins[2]);
		res.at<uchar>(Xcoef.at<double>(i, 1) - mins[1], Xcoef.at<double>(i, 0) - mins[0]) = color;
	}

	imshow("3D", res);
	waitKey();
}

int findK(int percentage, Mat Lambda) {
	double sumK = 0, sumD = 0;
	for (int i = 0; i < Lambda.rows; i++) {
		sumD += Lambda.at<double>(i, 0);
	}

	for (int k = 0; k < Lambda.rows; k++) {
		sumK += Lambda.at<double>(k, 0);
		if ((sumK * 100 / sumD) > percentage) {
			return k + 1;
		}
	}
}

int main()
{
	char fname[MAX_PATH];
	while (openFileDlg(fname)) {
		Mat X = readFile(fname);
		X = subtractMeanVector(X);
		Mat C = X.t() * X / (X.rows - 1);
		Mat Lambda, Q;
		eigen(C, Lambda, Q);
		Q = Q.t();
		printEigenvaluesAndMad(Lambda, X, Q);
		std::vector<Mat> aux = getPcaAndApproximate(Q, X, 3);
		showPca2D(aux[0]);
		showPca3D(aux[0]);
		printf("%d\n", findK(99, Lambda));
	}
	return 0;
}