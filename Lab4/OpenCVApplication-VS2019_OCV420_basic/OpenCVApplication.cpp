#include "stdafx.h"
#include "common.h"
#include <vector>
#include <algorithm>
#include <random>

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
			if (src.at<uchar>(i,j) == 0) {
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

int main()
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
		return 0;
		break;
	default:
		break;
	}
	return 0;
}