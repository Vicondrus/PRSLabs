#include "stdafx.h"
#include "common.h"
#include <vector>
#include <algorithm>
#include <random>

#define MAX_ITERATIONS 50

std::default_random_engine gen;

struct data_point {
	std::vector<int> dimensions;
};

struct clustered_point {
	data_point point;
	data_point cluster;
};

float computeDistance(data_point a, data_point b) {
	float dist = 0;
	for (int i = 0; i < a.dimensions.size();i++) {
		dist += (a.dimensions[i] - b.dimensions[i]) * (a.dimensions[i] - b.dimensions[i]);
	}
	dist = sqrtf(dist);
	return dist;
}

data_point cloneData_point(data_point x) {
	data_point y;
	for (int dimension : x.dimensions) {
		y.dimensions.push_back(dimension);
	}
	return y;
}

bool compareData_points(data_point a, data_point b) {
	for (int i = 0; i < a.dimensions.size(); i++) {
		if (a.dimensions[i] != b.dimensions[i]) {
			return false;
		}
	}
	return true;
}

bool checkIfIn(data_point point, std::vector<data_point> vector) {
	for (int i = 0; i < vector.size(); i++) {
		if (compareData_points(point, vector[i])) {
			return true;
		}
	}
	return false;
}

std::vector<clustered_point> kMeans(std::vector<data_point> points, int K, int limit) {
	std::vector<data_point> clusters;

	int d = points[0].dimensions.size();

	int n = points.size();

	std::uniform_int_distribution<int> dist_gen(0, n - 1);

	for (int k = 0; k < K; k++) {
		int index = dist_gen(gen);
		clusters.push_back(cloneData_point(points.at(index)));
	}

	std::vector<clustered_point> clustered_points;
	float dist = 0;

	bool modified = true;
	int counter = 0;
	while (modified) {
		clustered_points.clear();
		for (int i = 0; i < points.size(); i++) {
			float min_dist = 100000;
			data_point closest_cluster;
			for (int k = 0; k < K; k++) {
				data_point cluster = clusters.at(k);
				data_point current_point = points.at(i);
				float dist = computeDistance(current_point, cluster);
				if (dist < min_dist) {
					min_dist = dist;
					closest_cluster = cluster;
				}
			}
			clustered_point p;
			p.point = cloneData_point(points.at(i));
			p.cluster = closest_cluster;
			clustered_points.push_back(p);
		}

		std::vector<data_point> new_clusters(K);
		for (int k = 0; k < K; k++) {
			int count = 0;
			std::vector<int> sums(d);
			for (int i = 0; i < clustered_points.size(); i++) {
				if (compareData_points(clustered_points[i].cluster, clusters[k])) {
					for (int j = 0; j < d;j++) {
						sums[j] += clustered_points[i].point.dimensions[j];
					}
					count++;
				}
			}

			for (int j = 0; j < d; j++) {
				new_clusters[k].dimensions.push_back(sums[j] / count);
			}
		}

		for (int k = 0; k < K; k++) {
			if (compareData_points(clusters[k], new_clusters[k])) {
				modified = false;
			}
			else {
				clusters[k] = new_clusters[k];
				modified = true;
			}
		}

		if (counter >= limit) {
			modified = false;
		}
		counter++;
		printf("%d ", counter);
	}

	return clustered_points;
}

std::vector<data_point> getClusters(std::vector<clustered_point> points) {
	std::vector<data_point> clusters;
	for (int i = 0; i < points.size(); i++) {
		data_point cluster = points[i].cluster;
		if (!checkIfIn(cluster, clusters)) {
			clusters.push_back(cloneData_point(cluster));
		}
	}
	return clusters;
}

data_point getClosestCluster(data_point point, std::vector<data_point> clusters) {
	float min_dist = 100000;
	data_point closest_cluster;
	for (int k = 0; k < clusters.size(); k++) {
		data_point cluster = clusters.at(k);
		float dist = computeDistance(point, cluster);
		if (dist < min_dist) {
			min_dist = dist;
			closest_cluster = cluster;
		}
	}
	return closest_cluster;
}

void kMeans2d() {
	Mat src;

	char fname[MAX_PATH];
	while (openFileDlg(fname)) {
		src = imread(fname, IMREAD_GRAYSCALE);
		Mat colored(src.rows, src.cols, CV_8UC3);

		int K;
		printf("Input K: ");
		scanf("%d", &K);
		std::vector<data_point> points;
		for (int i = 0; i < src.rows; i++) {
			for (int j = 0; j < src.cols; j++) {
				if (src.at<uchar>(i, j) == 0) {
					data_point p;
					p.dimensions.push_back(i);
					p.dimensions.push_back(j);
					points.push_back(p);
				}
			}
		}

		std::vector<clustered_point> clustered_points = kMeans(points, K, MAX_ITERATIONS);

		for (int i = 0; i < clustered_points.size(); i++) {
			data_point cluster = clustered_points[i].cluster;
			int color = (cluster.dimensions[0] + cluster.dimensions[1]) * 13 % 255;
			colored.at<Vec3b>(clustered_points[i].point.dimensions[0], clustered_points[i].point.dimensions[1]) = Vec3b(50 * color % 255, 20 * color % 255, 10 * color % 255);
			circle(colored, Point(clustered_points[i].point.dimensions[1], clustered_points[i].point.dimensions[0]), 1, Scalar(50 * color % 255, 20 * color % 255, 10 * color % 255));
		}

		imshow("colored", colored);
		waitKey(0);

		Mat voronoi(src.rows, src.cols, CV_8UC3);
		std::vector<data_point> clusters = getClusters(clustered_points);
		for (int i = 0; i < src.rows; i++) {
			for (int j = 0; j < src.cols; j++) {
				data_point p;
				p.dimensions.push_back(i);
				p.dimensions.push_back(j);
				data_point cluster = getClosestCluster(p, clusters);
				int color = (cluster.dimensions[0] + cluster.dimensions[1]) * 13 % 255;
				voronoi.at<Vec3b>(i, j) = Vec3b(50 * color % 255, 20 * color % 255, 10 * color % 255);
			}
		}

		imshow("voronoi", voronoi);
		waitKey(0);
	}
}

void grayScaleClustering() {
	Mat src;

	char fname[MAX_PATH];
	while (openFileDlg(fname)) {
		src = imread(fname, IMREAD_GRAYSCALE);

		Mat dst(src.rows, src.cols, CV_8UC1);

		std::vector<data_point> points;

		for (int i = 0; i < src.rows; i++) {
			for (int j = 0; j < src.cols; j++) {
				data_point p;
				p.dimensions.push_back(src.at<uchar>(i, j));
				points.push_back(p);
			}
		}

		int K;
		printf("Input K: ");
		scanf("%d", &K);
		std::vector<clustered_point> clustered = kMeans(points, K, MAX_ITERATIONS);

		for (int i = 0; i < src.rows; i++) {
			for (int j = 0; j < src.cols; j++) {
				dst.at<uchar>(i,j) = clustered[i * ((int)src.cols) + j].cluster.dimensions[0];
			}
		}

		imshow("clustered", dst);
		waitKey(0);
	}
}

void colorClustering() {
	Mat src;

	char fname[MAX_PATH];
	while (openFileDlg(fname)) {
		src = imread(fname, IMREAD_COLOR);

		Mat dst(src.rows, src.cols, CV_8UC3);

		std::vector<data_point> points;

		for (int i = 0; i < src.rows; i++) {
			for (int j = 0; j < src.cols; j++) {
				data_point p;
				p.dimensions.push_back(src.at<Vec3b>(i, j)[0]);
				p.dimensions.push_back(src.at<Vec3b>(i, j)[1]);
				p.dimensions.push_back(src.at<Vec3b>(i, j)[2]);
				points.push_back(p);
			}
		}

		int K;
		printf("Input K: ");
		scanf("%d", &K);
		std::vector<clustered_point> clustered = kMeans(points, K, MAX_ITERATIONS);

		for (int i = 0; i < src.rows; i++) {
			for (int j = 0; j < src.cols; j++) {
				dst.at<Vec3b>(i, j)[0] = clustered[i * ((int)src.cols) + j].cluster.dimensions[0];
				dst.at<Vec3b>(i, j)[1] = clustered[i * ((int)src.cols) + j].cluster.dimensions[1];
				dst.at<Vec3b>(i, j)[2] = clustered[i * ((int)src.cols) + j].cluster.dimensions[2];
			}
		}

		imshow("clustered", dst);
		waitKey(0);
	}
}

int main()
{
	int choice = -1;
	printf("Choose (1) for 2dKmeans, (2) for grayscaleClustering, and (3) for colorClustering\n");
	scanf("%d", &choice);
	switch (choice) {
	case 1:
		kMeans2d();
		break;
	case 2:
		grayScaleClustering();
		break;
	case 3:
		colorClustering();
		break;
	default:
		break;
	}
	return 0;
}