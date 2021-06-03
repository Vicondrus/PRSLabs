#include "stdafx.h"
#include "common.h"
#include <vector>
#include <algorithm>
#include <random>

struct weaklearner {
	int feature_i;
	int threshold;
	int class_label;
	float error;
	int classify(Mat X) {
		if (X.at<int>(feature_i) < threshold)
			return class_label;
		else
			return -class_label;
	}
};

int getClass(Vec3b x) {
	if (x == Vec3b(0, 0, 255)) {
		return 1;
	}
	if (x == Vec3b(255, 0, 0)) {
		return -1;
	}
	return 0;
}

std::pair<Mat, Mat> getTrainingSet(Mat src) {
	Mat X(0, 2, CV_32SC1);
	Mat y(0, 1, CV_32SC1);
	for (int i = 0; i < src.rows; i++) {
		for (int j = 0; j < src.cols; j++) {
			int currClass = getClass(src.at<Vec3b>(i, j));
			if (currClass != 0) {
				Mat aux1 = (Mat_<int>(1, 2) << i, j);
				X.push_back(aux1);
				Mat aux2(1, 1, CV_32SC1, currClass);
				y.push_back(aux2);
			}
		}
	}
	return std::pair<Mat, Mat>(X, y);
}

weaklearner findWeakLearner(Mat X, Mat y, Mat w, MatSize imgSize) {
	weaklearner bestWeak;
	float minError = FLT_MAX;
	for (int j = 0; j < X.cols; j++) {
		for (int t = 0; t < imgSize[j]; t++) {
			int cls = -1;
			while (cls <= 1) {
				float err = 0;
				for (int i = 0; i < X.rows; i++) {
					int aux_cls;
					if (X.at<int>(i, j) < t) {
						aux_cls = cls;
					}
					else {
						aux_cls = -cls;
					}
					if (aux_cls * y.at<int>(i, 0) < 0) {
						err += w.at<float>(i, 0);
					}
				}
				if (err < minError) {
					minError = err;
					bestWeak = { j, t, cls, err };
				}
				cls += 2;
			}
		}
	}

	return bestWeak;
}

struct classifier {
	int T;
	std::vector<weaklearner> hs;
	std::vector<float> alphas;
	int classify(Mat X) {
		float sum = 0;
		for (int i = 0; i < T; i++) {
			sum += alphas[i] * hs[i].classify(X);
		}

		if (sum < 0) {
			return -1;
		}
		else {
			return 1;
		}
	}
};

std::pair<Mat, int> drawBoundary(Mat img, classifier clf) {
	Mat src(img.rows, img.cols, CV_8UC3);
	int error = 0;

	for (int i = 0; i < img.rows; i++) {
		for (int j = 0; j < img.cols; j++) {
			int v[2] = { i,j };
			Mat a(1, 2, CV_32SC1, v);
			int decision = clf.classify(a);
			if (decision * getClass(img.at<Vec3b>(i, j)) < 0) {
				error++;
			}
			if (img.at<Vec3b>(i, j) == Vec3b(0, 0, 255)) {
				src.at<Vec3b>(i, j) = img.at<Vec3b>(i, j);
			}
			else if (img.at<Vec3b>(i, j) == Vec3b(255, 0, 0)) {
				src.at<Vec3b>(i, j) = img.at<Vec3b>(i, j);
			}
			else {
				if (decision < 0) {
					src.at<Vec3b>(i, j) = Vec3b(255, 255, 0);
				}
				else {
					src.at<Vec3b>(i, j) = Vec3b(0, 255, 255);
				}
			}
		}
	}

	return std::pair<Mat, int>(src, error);
}

classifier adaboost(Mat src, int T) {
	std::pair<Mat, Mat> Xy = getTrainingSet(src);
	Mat X = Xy.first;
	Mat y = Xy.second;
	Mat W(X.rows, 2, CV_32FC1);
	for (int i = 0; i < X.rows; i++) {
		W.at<float>(i, 0) = (float)1 / X.rows;
	}

	classifier adaboost;
	adaboost.T = T;
	for (int t = 0; t < T; t++) {
		weaklearner wl = findWeakLearner(X, y, W, src.size);
		adaboost.hs.push_back(wl);
		adaboost.alphas.push_back(0.5f * log((1 - wl.error) / wl.error));
		float s = 0;
		for (int i = 0; i < X.rows; i++) {
			W.at<float>(i, 0) = W.at<float>(i, 0) * exp(-adaboost.alphas[t] * y.at<int>(i, 0) * wl.classify(X.row(i)));
			s += W.at<float>(i, 0);
		}
		for (int i = 0; i < X.rows; i++) {
			W.at<float>(i, 0) = W.at<float>(i, 0) / s;
		}
	}
	return adaboost;
}

int main()
{
	char fname[MAX_PATH];
	while (openFileDlg(fname))
	{
		Mat src = imread(fname, IMREAD_COLOR);
		classifier classif = adaboost(src, 1);
		std::pair<Mat, int> imageAndError = drawBoundary(src, classif);
		imshow("one learner", imageAndError.first);
		waitKey(0);
		int min = MAXINT;
		int bestT = -1;
		for (int i = 0; i < 100; i++) {
			classif = adaboost(src, i);
			imageAndError = drawBoundary(src, classif);
			if (min > imageAndError.second) {
				min = imageAndError.second;
				bestT = i;
			}
			imshow("lowest error", imageAndError.first);
			waitKey(1);
			printf("\r%d", i);
		}
		printf("\nResult: %d", bestT);
		classif = adaboost(src, 2000);
		imageAndError = drawBoundary(src, classif);
		imshow("lowest error", imageAndError.first);
		waitKey(0);
	}
	return 0;
}