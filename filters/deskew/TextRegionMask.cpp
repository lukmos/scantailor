#include "TextRegionMask.h"
#include "Dpi.h"
#include "TaskStatus.h"
#include "imageproc/RasterOp.h"
#include "imageproc/BWColor.h"
#include <QImage>
#include <algorithm>
#include <vector>

namespace deskew
{
namespace
{

#include "TextRegionWeights.inc"

using namespace imageproc;

int const patch_size = 96;
int const input_size = patch_size / 2;

template<int InputSize, int InputChannels, int OutputChannels, int KernelSize>
void convolve(float const* input, float const* weights, float const* biases,
	float* output)
{
	int const output_size = InputSize / 2;
	int const padding = KernelSize / 2;
	for (int channel = 0; channel < OutputChannels; ++channel) {
		float* plane = output + channel * output_size * output_size;
		std::fill(plane, plane + output_size * output_size, biases[channel]);
		for (int in = 0; in < InputChannels; ++in) {
			float const* source = input + in * InputSize * InputSize;
			for (int ky = 0; ky < KernelSize; ++ky) {
				int const offset_y = ky - padding;
				for (int kx = 0; kx < KernelSize; ++kx) {
					int const offset_x = kx - padding;
					float const weight = *weights++;
					int const first_x = offset_x < 0 ? (1 - offset_x) / 2 : 0;
					int const last_x = std::min(output_size, (InputSize - offset_x + 1) / 2);
					for (int y = 0; y < output_size; ++y) {
						int const sy = 2*y + offset_y;
						if (sy < 0 || sy >= InputSize) {
							continue;
						}
						float const* row = source + sy * InputSize;
						for (int x = first_x; x < last_x; ++x) {
							plane[y * output_size + x] += row[2*x + offset_x] * weight;
						}
					}
				}
			}
		}
		for (int i = 0; i < output_size * output_size; ++i) {
			plane[i] = std::max(0.0f, plane[i]);
		}
	}
}

struct PatchClassifier
{
	float input[input_size * input_size];
	float first[16 * 24 * 24];
	float second[32 * 12 * 12];
	float third[48 * 6 * 6];

	bool selectsText()
	{
		convolve<48, 1, 16, 5>(input, conv1_weights, conv1_bias, first);
		convolve<24, 16, 32, 3>(first, conv2_weights, conv2_bias, second);
		convolve<12, 32, 48, 3>(second, conv3_weights, conv3_bias, third);
		float logit = classifier_bias[0];
		for (int channel = 0; channel < 48; ++channel) {
			float sum = 0;
			for (int i = 0; i < 36; ++i) {
				sum += third[channel * 36 + i];
			}
			logit += sum / 36.0f * classifier_weights[channel];
		}
		// log(0.90 / 0.10), matching the probability cutoff in training evaluation.
		return logit >= 2.1972246f;
	}
};

void sampleOffsets(int extent, int (&offsets)[patch_size])
{
	for (int i = 0; i < patch_size; ++i) {
		offsets[i] = static_cast<int>((i + 0.5) * extent / patch_size);
	}
}

bool readPatch(QImage const& gray, int left, int top,
	int const (&xs)[patch_size], int const (&ys)[patch_size], float* input)
{
	bool has_ink = false;
	for (int y = 0; y < patch_size; y += 2) {
		uchar const* first = gray.constScanLine(top + ys[y]);
		uchar const* second = gray.constScanLine(top + ys[y+1]);
		for (int x = 0; x < patch_size; x += 2) {
			int const a = left + xs[x];
			int const b = left + xs[x+1];
			int const ink = 1020 - first[a] - first[b] - second[a] - second[b];
			*input++ = ink / 1020.0f;
			has_ink = has_ink || ink != 0;
		}
	}
	return has_ink;
}

int tileExtent(int dpi, int limit)
{
	if (dpi <= 0) {
		return patch_size;
	}
	double const extent = patch_size * (dpi / 300.0);
	return extent > limit ? limit + 1 : std::max(2, static_cast<int>(extent + 0.5));
}

std::vector<int> tilePositions(int limit, int extent)
{
	std::vector<int> positions;
	int const last = limit - extent;
	for (int start = 0; start < last; start += extent / 2) {
		positions.push_back(start);
	}
	positions.push_back(last);
	return positions;
}

}

imageproc::BinaryImage
retainTextRegions(imageproc::BinaryImage const& image, QImage const& gray,
	Dpi const& dpi, TaskStatus const& status)
{
	int const width = tileExtent(dpi.horizontal(), image.width());
	int const height = tileExtent(dpi.vertical(), image.height());
	if (width > image.width() || height > image.height()) {
		return BinaryImage();
	}
	std::vector<int> const columns = tilePositions(image.width(), width);
	std::vector<int> const rows = tilePositions(image.height(), height);
	int xs[patch_size], ys[patch_size];
	sampleOffsets(width, xs);
	sampleOffsets(height, ys);

	BinaryImage selected(image.size(), WHITE);
	PatchClassifier classifier;
	for (size_t y = 0; y < rows.size(); ++y) {
		for (size_t x = 0; x < columns.size(); ++x) {
			status.throwIfCancelled();
			int const left = columns[x];
			int const top = rows[y];
			if (!readPatch(gray, left, top, xs, ys, classifier.input) ||
				!classifier.selectsText()) {
				continue;
			}
			QRect region(left, top, width, height);
			region.adjust(x == 0 ? 0 : width/6, y == 0 ? 0 : height/6,
				x+1 == columns.size() ? 0 : -width/6, y+1 == rows.size() ? 0 : -height/6);
			rasterOp<RopSrc>(selected, region, image, region.topLeft());
		}
	}
	return selected;
}

}
