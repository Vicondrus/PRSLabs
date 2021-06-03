#include "stdafx.h"
#include "common.h"
#include <vector>
#include <random>

using namespace std;

const int nrclasses = 6;
char classes[nrclasses][10] =
{ "beach", "city", "desert", "forest", "landscape", "snow" };

vector<int> calcHist(Mat img, int nr_bins, vector<int> hist)
{
	vector<int> histogram;

	int height = img.rows;
	int width = img.cols;

	float D = 256 / nr_bins;

	for (int color = 0; color < 3; color++)
	{
		vector<int> hist(nr_bins);
		for (int i = 0; i < height; i++)
		{
			for (int j = 0; j < width; j++)
			{
				hist.at(img.at<Vec3b>(i, j)[color] / D)++;
			}
		}
		histogram.insert(histogram.end(), hist.begin(), hist.end());
	}

	return histogram;
}

vector<int> computeHistogram(Mat src)
{
	vector<int> histogram(257);

	int height = src.rows;
	int width = src.cols;

	Mat_<uchar> binary(height, width);

	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			uchar value = src.at<uchar>(i, j);
			histogram[value]++;
		}
	}

	return histogram;
}

pair<Mat, Mat> computeFeatureMatrixAndClassLabel(int m)
{
	int nrinst = 672;
	int feature_dim = m * 3;

	Mat X(nrinst, feature_dim, CV_32FC1);
	Mat y(nrinst, 1, CV_8UC1);

	char fname[MAX_PATH];

	int fileNr = 0, rowX = 0;

	for (int c = 0; c < nrclasses; c++)
	{
		while (1) {
			sprintf(fname, "prs_res_KNN/train/%s/%06d.jpeg", classes[c], fileNr++);
			Mat img = imread(fname, IMREAD_COLOR);
			if (img.cols == 0) break;

			//calculate the histogram in hist
			vector<int> histogram = computeHistogram(img);
			vector<int> hist = calcHist(img, m, histogram);

			for (int d = 0; d < hist.size(); d++)
				X.at<float>(rowX, d) = hist[d];
			y.at<uchar>(rowX) = c;
			rowX++;
		}
	}

	return pair<Mat, Mat>(X, y);
}

bool compareByDistance(const pair<float, int>& a, const pair<float, int>& b)
{
	return (a.second < b.second);
}

int knn(Mat img, int K, int m, Mat X, Mat y)
{
	vector<int> histogram = computeHistogram(img);
	vector<int> hist = calcHist(img, m, histogram);

	int height = X.rows;
	int width = X.cols;

	vector<pair<float, int>> distances;

	for (int i = 0; i < height; i++)
	{
		float distance = 0;
		for (int j = 0; j < width; j++)
		{
			distance += (abs(hist.at(j) - X.at<float>(i, j)));
		}

		distances.push_back(pair<float, int>(distance, y.at<uchar>(i)));
	}

	sort(distances.begin(), distances.end(), compareByDistance);

	vector<int> votes(K);
	for (int i = 0; i < K; i++)
	{
		votes[(int)(distances.at(i).second)]++;
	}

	int nearestClass = -1;
	int maxVote = -1;

	for (int i = 0; i < K; i++)
	{
		if (votes.at(i) > maxVote)
		{
			maxVote = votes.at(i);
			nearestClass = i;
		}
	}

	return nearestClass;
}

void knnTestImages(int K, int m)
{
	char fname[MAX_PATH];

	int fileNr = 0;

	pair<Mat, Mat> res = computeFeatureMatrixAndClassLabel(m);

	//create confusion matrix
	Mat_<float> C(nrclasses, nrclasses);

	for (int i = 0; i < nrclasses; i++)
	{
		for (int j = 0; j < nrclasses; j++)
		{
			C(i, j) = 0;
		}
	}

	int nrClassifications = 0;
	float correctlyClassified = 0;

	for (int c = 0; c < nrclasses; c++)
	{
		while (1) {
			sprintf(fname, "prs_res_KNN/test/%s/%06d.jpeg", classes[c], fileNr++);
			Mat img = imread(fname, IMREAD_COLOR);
			if (img.cols == 0) break;

			int nearestC = knn(img, K, m, res.first, res.second);

			C.at<float>(c, nearestC) += 1;

			if (c == nearestC)
			{
				correctlyClassified++;
			}
			nrClassifications++;
		}
	}

	//print confusion matrix
	for (int i = 0; i < nrclasses; i++)
	{
		for (int j = 0; j < nrclasses; j++)
		{
			cout << C(i, j) << " ";
		}
		cout << endl;
	}

	float accuracy = (float)correctlyClassified / nrClassifications;
	cout << "Acc: " << accuracy << endl;

	cin >> accuracy;
}

int main()
{
	printf("Choose method:\n");
	Mat_<uchar> img_grayscale;
	Mat_<int> points;
	Mat_<Vec3b> img;
	char fname[MAX_PATH];
	int op;
	do
	{
		system("cls");
		destroyAllWindows();
		printf("Menu:\n");
		printf(" 1 - KNN\n");
		printf(" 0 - Exit\n\n");
		printf("Option: ");
		scanf("%d", &op);
		switch (op)
		{
		case 1:
			knnTestImages(8, 16);
			break;
		default:

			break;
		}
	} while (op != 0);
	return 0;
}