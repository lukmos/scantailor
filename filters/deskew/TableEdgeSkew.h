#ifndef SCANTAILOR_DESKEW_TABLE_EDGE_SKEW_H_
#define SCANTAILOR_DESKEW_TABLE_EDGE_SKEW_H_

class QImage;
class QRect;
class TaskStatus;

namespace deskew
{

// Returns a counter-clockwise correction when repeated color boundaries agree.
bool findTableEdgeCorrection(QImage const& image, QRect const& area, int rotation,
	double resolution_ratio, TaskStatus const& status, double& correction);

}

#endif
