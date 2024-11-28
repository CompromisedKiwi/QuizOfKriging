#include <CGAL/Surface_mesh.h>

#include <Eigen/Dense>

#include <vector>

#include "UsingDeclaration.h"

using Mesh = CGAL::Surface_mesh<Point_3>;

#ifndef CGALSTUDY_KRIGING_H
#define CGALSTUDY_KRIGING_H

struct Semivariance {
    Semivariance(double _c0, double _c1, double _maxDistance) :
            c0(_c0), c1(_c1), maxDistance(_maxDistance) {}

    double c0;
    double c1;
    double maxDistance;
    double Evaluate(const double & distance) const {
        double tmp = distance / maxDistance;
        return c0 + c1 * (1.5 * tmp - 0.5 * pow(tmp, 3));
    }
};

void Kriging(Mesh & dsm_mesh, const std::vector<Point_3> & points, int sampleLevel);

#endif //CGALSTUDY_KRIGING_H
