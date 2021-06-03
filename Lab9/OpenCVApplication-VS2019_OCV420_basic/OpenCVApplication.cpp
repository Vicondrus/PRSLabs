#include "stdafx.h"
#include "common.h"
#include <vector>
#include <algorithm>
#include <random>

Mat binarize(Mat img, uchar threshold) {
	for (int i = 0; i < img.rows; i++) {
		for (int j = 0; j < img.cols; j++) {
			if (img.at<uchar>(i, j) < threshold) {
				img.at<uchar>(i, j) = 0;
			}
			else {
				img.at<uchar>(i, j) = 255;
			}

		}
	}
	return img;
}

std::vector<Mat> getXYP(int limit, int classes) {
	char fname[256]; 
	int index = 0;
	Mat X(limit*classes, 28*28, CV_8UC1);
	Mat y(limit*classes, 1, CV_8UC1);
	Mat p(classes, 1, CV_64FC1);

	for (int c = 0; c < classes; c++) {
		for (int ind=0; ind < limit; ind++) {
			sprintf(fname, "Data/train/%d/%06d.png", c, ind);
			Mat img = imread(fname, 0);
			if (img.cols == 0) break;
			
			img = binarize(img, 128);

			int x = 0;
			for (int i = 0; i < img.rows; i++) {
				for (int j = 0; j < img.cols; j++) {
					X.at<uchar>(index, x) = img.at<uchar>(i, j);
					x++;
				}
			}
			y.at<uchar>(index, 0) = c;
			index++;
			printf("\rRead %d", index);
		}
		p.at<double>(c, 0) = (double)limit / (limit * classes);
	}

	std::vector<Mat> res;
	res.push_back(X);
	res.push_back(y);
	res.push_back(p);
	return res;
}

Mat computeLikelihood(int classes, int limit, Mat X, Mat y) {
	const int d = 28 * 28;
	Mat likelihood(classes, d, CV_64FC1);
	for (int i = 0; i < d; i++) {
		for (int j = 0; j < classes; j++) {
			likelihood.at<double>(j, i) = 0;
		}
	}

	for (int i = 0; i < d; i++) {
		for (int j = 0; j < limit * classes; j++) {
			if (X.at<uchar>(j, i) == 255) {
				likelihood.at<double>(y.at<uchar>(j), i)++;
			}
		}
	}

	for (int i = 0; i < likelihood.rows; i++) {
		for (int j = 0; j < likelihood.cols; j++) {
			double l = likelihood.at<double>(i, j);
			likelihood.at<double>(i, j) = (l + 1) / (limit + classes);
		}
	}

	//for (int i = 0; i < d; i++) {
	//	for (int j = 0; j < classes; j++) {
	//		printf("%lf ", likelihood.at<double>(j, i));
	//	}
	//	printf("\n");
	//}

	return likelihood;
}

std::pair<int, std::vector<double>> classifyImage(Mat img, Mat likelihood, Mat p, int classes) {
	img = binarize(img, 128);

	Mat features(1, 28 * 28, CV_8UC1);
	int x = 0;
	for (int i = 0; i < img.rows; i++) {
		for (int j = 0; j < img.cols; j++) {
			features.at<uchar>(0, x) = img.at<uchar>(i, j);
			x++;
		}
	}

	std::vector<double> probabilities(classes);
	for (int c = 0; c < classes; c++) {
		double probability = 0;
		for (int i = 0; i < 28 * 28; i++) {
			if (features.at<uchar>(0, i) == 255) {
				double q = likelihood.at<double>(c, i);
				probability += log(q);
			}
			else {
				double q = likelihood.at<double>(c, i);
				probability += log(1 - q);
			}
		}
		probability += log(p.at<double>(c, 0));
		probabilities[c] = probability;
	}
	int maxc = -1;
	double maxx = -10000;
	for (int i = 0; i < classes; i++) {
		if (probabilities[i] > maxx) {
			maxc = i;
			maxx = probabilities[i];
		}
	}

	return std::pair<int, std::vector<double>>(maxc, probabilities);
}

void printErrorAndConfusion() {
	const int classes = 10;
	const int limit = 5800;
	const int testLimit = 840;
	char fname[256];
	std::vector<Mat> xyp = getXYP(limit, classes);
	Mat likelihood = computeLikelihood(classes, limit, xyp[0], xyp[1]);
	Mat confusion(classes, classes, CV_32FC1);
	int total = 0;
	double wrong = 0;
	for (int i = 0; i < classes; i++) {
		for (int j = 0; j < classes; j++) {
			confusion.at<int>(i, j) = 0;
		}
	}
	for (int c = 0; c < classes; c++) {
		for (int ind = 0; ind < testLimit; ind++) {
			sprintf(fname, "Data/test/%d/%06d.png", c, ind);
			Mat img = imread(fname, 0);
			if (img.cols == 0) break;

			img = binarize(img, 128);
			std::pair<int, std::vector<double>> pair = classifyImage(img, likelihood, xyp[2], classes);

			confusion.at<int>(c, pair.first)++;
			if (c != pair.first) {
				wrong++;
			}
			total++;
			printf("\rPredicted %d", total);
		}
	}

	printf("\nError: %lf\n", (double)(wrong / total));
	for (int i = 0; i < classes; i++) {
		for (int j = 0; j < classes; j++) {
			printf("%3d ", confusion.at<int>(i, j));
		}
		printf("\n");
	}
}

void classifySingleImage() {
	const int classes = 10;
	const int limit = 5800;
	std::vector<Mat> xyp = getXYP(limit, classes);
	Mat likelihood = computeLikelihood(classes, limit, xyp[0], xyp[1]);
	char fname[MAX_PATH];
	while (openFileDlg(fname))
	{
		Mat src;
		src = imread(fname, IMREAD_GRAYSCALE);
		std::pair<int, std::vector<double>> pair = classifyImage(src, likelihood, xyp[2], classes);
		for (int i = 0; i < classes; i++) {
			printf("Probability for class %d is %f \n", i, pair.second[i]);
		}
		printf("Class is %d\n", pair.first);
	}
}

int main()
{
	printf("Options:\n");
	printf("Classify Single Image - 1\n");
	printf("Error and Confusion for Test Dataset - 2\n");
	printf("Choice: ");
	int choice;
	scanf("%d", &choice);
	if (choice == 1) {
		classifySingleImage();
	}
	else if (choice == 2) {
		printErrorAndConfusion();
	}
	return 0;
}