#include <vector>

#include "UsingDeclaration.h"

#ifndef CGALSTUDY_MYBBOX_H
#define CGALSTUDY_MYBBOX_H

struct MyBbox {
    MyBbox(const std::vector<Point_3> & points) {
        for(auto point : points) {
            minx = std::min(minx, point.x());
            miny = std::min(miny, point.y());
            minz = std::min(minz, point.z());
            maxx = std::max(maxx, point.x());
            maxy = std::max(maxy, point.y());
            maxz = std::max(maxz, point.z());
        }
    }

    double minx {0.}, miny {0.}, minz {0.}, maxx {0.}, maxy {0.}, maxz {0.};
};

#endif //CGALSTUDY_MYBBOX_H
