#include "TableEdgeSkew.h"
#include "TaskStatus.h"
#include "imageproc/Constants.h"
#include <QImage>
#include <QPoint>
#include <QTransform>
#include <algorithm>
#include <cmath>
#include <vector>

namespace deskew
{
namespace
{

int const strip_count = 25;
int const minimum_rows = 8;

struct Edge
{
	double y;
	double angle;
	int support;
};

double median(std::vector<double>& values)
{
	std::sort(values.begin(), values.end());
	size_t const middle = values.size() / 2;
	return values.size() % 2 ? values[middle] : (values[middle-1] + values[middle]) / 2;
}

bool moreSupport(Edge const& first, Edge const& second)
{
	return first.support > second.support;
}

bool separated(std::vector<Edge> const& edges, double y, int distance)
{
	for (size_t i = 0; i < edges.size(); ++i) {
		if (std::fabs(edges[i].y - y) <= distance) {
			return false;
		}
	}
	return true;
}

void fitEdge(std::vector<QPoint> const& points, int width, int y,
	double resolution_ratio, std::vector<Edge>& edges)
{
	if (points.size() < 9 || points.back().x() - points.front().x() < .34 * width) {
		return;
	}
	std::vector<double> slopes;
	for (size_t i = 0; i < points.size(); ++i) {
		for (size_t j = i + 1; j < points.size(); ++j) {
			int const dx = points[j].x() - points[i].x();
			if (dx > .16 * width) {
				slopes.push_back(double(points[j].y() - points[i].y()) / dx);
			}
		}
	}
	double const initial_slope = median(slopes);
	std::vector<double> offsets;
	for (size_t i = 0; i < points.size(); ++i) {
		offsets.push_back(points[i].y() - initial_slope * points[i].x());
	}
	double const initial_offset = median(offsets);

	int support = 0;
	int first_x = width;
	int last_x = 0;
	double sum_x = 0, sum_y = 0, sum_xx = 0, sum_xy = 0;
	for (size_t i = 0; i < points.size(); ++i) {
		double const x = points[i].x();
		double const dy = points[i].y();
		if (std::fabs(dy - initial_slope * x - initial_offset) <= 2.0) {
			first_x = std::min(first_x, points[i].x());
			last_x = points[i].x();
			++support;
			sum_x += x;
			sum_y += dy;
			sum_xx += x * x;
			sum_xy += x * dy;
		}
	}
	if (support < 8 || support < .70 * points.size() || last_x - first_x < .34 * width) {
		return;
	}
	double const slope = (support * sum_xy - sum_x * sum_y) /
		(support * sum_xx - sum_x * sum_x);
	double const angle = std::atan(slope * resolution_ratio) / imageproc::constants::DEG2RAD;
	if (std::fabs(angle) <= 1.25) {
		Edge const edge = { y + (sum_y - slope * sum_x) / support + slope * width / 2, angle, support };
		edges.push_back(edge);
	}
}

void collectEdges(std::vector<double> const (&gradient)[strip_count],
	int const (&centers)[strip_count], int reference, int width, int height,
	int radius, double resolution_ratio, TaskStatus const& status, std::vector<Edge>& edges)
{
	int const skip = std::max(8, height / 30);
	std::vector<double> const& signal = gradient[reference];
	std::vector<int> peaks;
	for (int y = skip + radius + 1; y < int(signal.size()) - skip - radius; ++y) {
		if (signal[y] >= 15) {
			peaks.push_back(y);
		}
	}
	// Rank by strength before suppressing repeated detections of the same edge.
	std::vector<std::pair<double, int> > ranked;
	for (size_t i = 0; i < peaks.size(); ++i) {
		ranked.push_back(std::make_pair(-signal[peaks[i]], peaks[i]));
	}
	std::sort(ranked.begin(), ranked.end());
	std::vector<Edge> processed;
	for (size_t i = 0; i < ranked.size(); ++i) {
		status.throwIfCancelled();
		int const y = ranked[i].second;
		if (!separated(processed, y, 2 * radius)) {
			continue;
		}
		Edge const peak = { double(y), 0, 0 };
		processed.push_back(peak);
		std::vector<QPoint> points;
		for (int j = 0; j < strip_count; ++j) {
			std::vector<double>::const_iterator const first = gradient[j].begin() + y - radius;
			std::vector<double>::const_iterator const best = std::max_element(first, first + 2 * radius + 1);
			if (*best >= std::max(12.0, .25 * signal[y])) {
				points.push_back(QPoint(centers[j], int(best - first) - radius));
			}
		}
		fitEdge(points, width, y, resolution_ratio, edges);
	}
}

}

bool
findTableEdgeCorrection(QImage const& image, QRect const& area, int rotation,
	double resolution_ratio, TaskStatus const& status, double& correction)
{
	if (area.width() < 300 || area.height() < 300 || image.isGrayscale()) {
		return false;
	}
	QImage const color(image.copy(area).transformed(QTransform().rotate(rotation)).convertToFormat(QImage::Format_RGB32));
	int const width = color.width();
	int const height = color.height();
	int const half_width = std::max(5, width / 125);
	int centers[strip_count];
	std::vector<double> gradient[strip_count];
	for (int i = 0; i < strip_count; ++i) {
		status.throwIfCancelled();
		centers[i] = int(width * (.08 + .84 * i / (strip_count - 1)));
		std::vector<int> chroma(height);
		for (int y = 0; y < height; ++y) {
			QRgb const* row = reinterpret_cast<QRgb const*>(color.constScanLine(y));
			for (int x = centers[i] - half_width; x < centers[i] + half_width; ++x) {
				chroma[y] += std::max(qRed(row[x]), qGreen(row[x])) - qBlue(row[x]);
			}
		}
		gradient[i].resize(height - 4);
		for (int y = 0; y < height - 4; ++y) {
			gradient[i][y] = (chroma[y+4] - chroma[y]) / (2.0 * half_width);
		}
	}

	int const radius = std::max(10, width / 110);
	std::vector<Edge> edges;
	int const references[] = {7, 11, 15};
	for (size_t i = 0; i < sizeof(references) / sizeof(references[0]); ++i) {
		collectEdges(gradient, centers, references[i], width, height,
			radius, resolution_ratio, status, edges);
	}
	std::stable_sort(edges.begin(), edges.end(), moreSupport);
	std::vector<Edge> distinct;
	std::vector<double> angles;
	for (size_t i = 0; i < edges.size(); ++i) {
		if (separated(distinct, edges[i].y, radius)) {
			distinct.push_back(edges[i]);
			angles.push_back(edges[i].angle);
		}
	}
	if (angles.size() < minimum_rows) {
		return false;
	}

	double const consensus = median(angles);
	std::vector<double> agreeing;
	double top = height;
	double bottom = 0;
	for (size_t i = 0; i < distinct.size(); ++i) {
		if (std::fabs(distinct[i].angle - consensus) <= .15) {
			agreeing.push_back(distinct[i].angle);
			top = std::min(top, distinct[i].y);
			bottom = std::max(bottom, distinct[i].y);
		}
	}
	if (agreeing.size() < minimum_rows || agreeing.size() < .75 * distinct.size() ||
		bottom - top < .40 * height) {
		return false;
	}
	correction = median(agreeing);
	return true;
}

}
