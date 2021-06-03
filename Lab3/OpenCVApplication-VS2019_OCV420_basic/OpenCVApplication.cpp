#include "stdafx.h"
#include "common.h"
#include <vector>
#include <algorithm>
#include <random>

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
		int theta = 0, rho=0;
		float beta = 0;
		int maxx = 0;
		int thetaStep = 1;

		for (int i = 0; i < src.rows; i++) {
			for (int j = 0; j < src.cols; j++) {
				if (src.at<uchar>(i, j) == 255) {
					for (theta = 0; theta < 360; theta+=thetaStep) {

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

int main()
{
	int k, n;
	printf("Input window size: ");
	scanf("%d", &n);
	printf("Input number of local maxima: ");
	scanf("%d", &k);
	hough(n, k);
	return 0;
}