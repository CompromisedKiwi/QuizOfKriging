#include <iostream>
#include <array>
#include <cmath>

#include "Kriging.h"
#include "MyBbox.h"

void Kriging(Mesh & remeshedMesh, const std::vector<Point_3> & points, int sampleLevel) {
    // Ordinary Kriging Interpolation
    // https://pro.arcgis.com/en/pro-app/latest/tool-reference/3d-analyst/how-kriging-works.htm

    // Distance and semi-variances as sample points.
    size_t pointAmount = points.size();
    std::vector<std::array<double, 2>> samples;
    double maxDistance = 0;
    for (size_t i = 0; i < pointAmount - 1; ++i) {
        Point_2 pi(points[i].x(), points[i].y());

        for (size_t j = i + 1; j < pointAmount; ++j) {
            Point_2 pj(points[j].x(), points[j].y());
            double distance = CGAL::sqrt(CGAL::squared_distance(pi, pj));
            maxDistance = std::max(distance, maxDistance);
            double semiVariance = 0.5 * pow(points[i].z() - points[j].z(), 2);

            samples.push_back(std::array<double, 2>({distance, semiVariance}));
        }
    }

    // Fitting of Semi-Variance function.
    std::cout << "Start fitting." << std::endl;
    Eigen::MatrixXd A(samples.size(), 2);
    Eigen::VectorXd b(samples.size());
    for (int i = 0; i < samples.size(); ++i) {
        double tmp = samples[i][0] / maxDistance;
        A(i, 0) = 1.;
        A(i, 1) = 1.5 * tmp - 0.5 * pow(tmp, 3);
        b(i) = samples[i][1];
    }
    Eigen::Vector2d solution = (A.transpose() * A).ldlt().solve(A.transpose() * b);
    Semivariance sv(solution[0], solution[1], maxDistance);
    std::cout << "Fitting solution (C0, C1):\n" << solution << std::endl;

    // Construct the semi-variance matrix.
    std::cout << "Construct the semi-variance matrix, point amount: " << pointAmount << std::endl;
    std::vector<double> heightsVector;
    Eigen::MatrixXd SVMatrix = Eigen::MatrixXd::Ones(pointAmount + 1, pointAmount + 1);
    SVMatrix(pointAmount, pointAmount) = 0.;
    Eigen::VectorXd heights = Eigen::VectorXd::Zero(pointAmount + 1);

    for (int indexI = 0; indexI < pointAmount; ++indexI) {
        heights(indexI) = points[indexI].z();

        for (int indexJ = indexI; indexJ < pointAmount; ++indexJ) {
            double distance = CGAL::sqrt(CGAL::squared_distance(
                    Point_2(points[indexI].x(), points[indexI].y()),
                    Point_2(points[indexJ].x(), points[indexJ].y())));
            double svij = sv.Evaluate(distance);
            SVMatrix(indexI, indexJ) = svij;
            SVMatrix(indexJ, indexI) = svij;
        }
    }
    auto decomp = SVMatrix.ldlt();

    // Interpolation
    std::cout << "Interpolation" << std::endl;
    int count = 0;
    Eigen::VectorXd tempB(pointAmount + 1);  // b in Ax=b
    tempB(pointAmount) = 1.;

    int remeshedPointAmount = remeshedMesh.num_vertices();
    for (int i = 0; i < remeshedPointAmount; ++i) {
        Point_3 &p = remeshedMesh.point(static_cast<CGAL::SM_Vertex_index>(i));

        // Construct b in Ax=b
        for (int j = 0; j < pointAmount; ++j) {
            double d =  CGAL::sqrt(CGAL::squared_distance(
                    Point_2(points[j].x(), points[j].y()), Point_2(p.x(), p.y())));
            tempB(j) = sv.Evaluate(d);
        }
        Eigen::VectorXd lbd = decomp.solve(tempB);
        remeshedMesh.point(static_cast<CGAL::SM_Vertex_index>(i)) = Point_3(
                p.x(), p.y(), lbd.dot(heights));

        count++;
        if (count % 100 == 0) {
            std::cout << "Lambda sum: " << lbd.sum()
                      << "\t\tInterpolation progress: " << count << " / " << remeshedPointAmount << std::endl;
        }
    }
}