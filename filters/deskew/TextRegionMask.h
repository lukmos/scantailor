#ifndef SCANTAILOR_DESKEW_TEXT_REGION_MASK_H_
#define SCANTAILOR_DESKEW_TEXT_REGION_MASK_H_

#include "imageproc/BinaryImage.h"

class Dpi;
class QImage;
class TaskStatus;

namespace deskew
{

// Returns a null image if the page is too small for a model tile.
imageproc::BinaryImage retainTextRegions(imageproc::BinaryImage const& image,
	QImage const& gray, Dpi const& dpi, TaskStatus const& status);

}

#endif
