//
// Created by 徐溶延 on 2020/10/12.
//

#include "MathUtils.h"
#include <Eigen/Dense>

namespace Surface_Mesh {
	double angle(const Point &p1, const Point& p2) {
		return atan2(p1.cross(p2).norm(), p1.dot(p2));
	}
}