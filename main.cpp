#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/boost/graph/graph_traits_Delaunay_triangulation_2.h>
#include <CGAL/Polygon_mesh_processing/clip.h>

#include <CGAL/Qt/Basic_viewer.h>
#include <CGAL/draw_triangulation_2.h>
#include <CGAL/draw_surface_mesh.h>
#include <QMainWindow>

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <cmath>

#include "UsingDeclaration.h"
#include "Kriging.h"
#include "MyBbox.h"


using Projection_traits = CGAL::Projection_traits_xy_3<Kernel>;
using TIN = CGAL::Delaunay_triangulation_2<Projection_traits>;
using Plane_3 = Kernel::Plane_3;

using Primitive = CGAL::AABB_face_graph_triangle_primitive<Mesh>;
using AABBTree = CGAL::AABB_tree<CGAL::AABB_traits_3<Kernel, Primitive>>;
using Ray = Kernel::Ray_3;

std::vector<Point_3> LoadPoints(const char * filename);
std::vector<double> StringSplit(std::string& str, const char delimiter);


int main(int argc, char** argv)
{
    // CGAL GIS TUTORIAL: https://doc.cgal.org/latest/Manual/tuto_gis.html

    CGAL::Graphics_scene scene1, scene2, scene3, scene4;

    // Read point data.
    std::cout << "Reading points..." << std::endl;
    std::vector<Point_3> points = LoadPoints("../resources/points.txt");
    MyBbox bbox(points);

    // Delaunay triangulation.
    std::cout << "Delaunay triangulation..." << std::endl;
    TIN dsm(points.begin(), points.end());
    CGAL::add_to_graphics_scene(dsm, scene1);

    // Convert triangulation to a mesh.
    std::cout << "Converting triangulation to mesh..." << std::endl;
    Mesh dsm_mesh;
    CGAL::copy_face_graph(dsm, dsm_mesh);

    // Remesh (to square mesh)
    std::cout << "Remeshing..." << std::endl;
    constexpr double interval { 2. };
    std::vector<Point_3> newMeshPoints;

    Kernel::Vector_3 direction(0, 0, 1);
    AABBTree tree(faces(dsm_mesh).first, faces(dsm_mesh).second, dsm_mesh);
    for(double x = bbox.minx; x <= bbox.maxx; x += interval) {
        for(double y = bbox.miny; y <= bbox.maxy; y += interval) {
            // Calculate intersection.
            auto intersection = tree.first_intersection(Ray(Point_3(x, y, -10), direction));
            if (intersection) {
                if(std::get_if<Point_3>(&(intersection->first))){
                    const Point_3 * hitPoint =  std::get_if<Point_3>(&(intersection->first) );
                    newMeshPoints.emplace_back(x, y, hitPoint->z());
                }
            }
        }
    }

    TIN remeshedDSM(newMeshPoints.begin(), newMeshPoints.end());
    Mesh remeshedMesh;
    CGAL::copy_face_graph(remeshedDSM, remeshedMesh);
    CGAL::add_to_graphics_scene(remeshedMesh, scene2);

    std::cout << "Kriging..." << std::endl;
    Kriging(remeshedMesh, points, 3000);
    CGAL::add_to_graphics_scene(remeshedMesh, scene3);


    // Contouring
    std::cout << "Contouring..." << std::endl;
    double minHeight = bbox.minz;
    double maxHeight = bbox.maxz;

    Mesh meshToBeSplited;
    CGAL::copy_face_graph(remeshedMesh, meshToBeSplited);
    constexpr int contourAmount { 16 };
    for (int i = 1; i <= contourAmount; ++i) {
        double isovalues = minHeight + (maxHeight - minHeight) / contourAmount * i;
        CGAL::Polygon_mesh_processing::split(meshToBeSplited, Plane_3(0, 0, 1, -isovalues));
    }
    CGAL::add_to_graphics_scene(meshToBeSplited, scene4);


    QApplication app(argc, argv);
    QMainWindow* mainWindow=new QMainWindow;
    QWidget *centralWidget = new QWidget(mainWindow);
    QGridLayout *layout = new QGridLayout(mainWindow);

    CGAL::Qt::Basic_viewer bv1(mainWindow, scene1);
    layout->addWidget(&bv1, 0, 0);
    CGAL::Qt::Basic_viewer bv2(mainWindow, scene2);
    layout->addWidget(&bv2, 0, 1);
    CGAL::Qt::Basic_viewer bv3(mainWindow, scene3);
    layout->addWidget(&bv3, 1, 0);
    CGAL::Qt::Basic_viewer bv4(mainWindow, scene4);
    layout->addWidget(&bv4, 1, 1);

    centralWidget->setLayout(layout);
    mainWindow->setCentralWidget(centralWidget);
    mainWindow->show();
    app.exec();

    return 0;
}

std::vector<Point_3> LoadPoints(const char * filename)
{
    std::ifstream file(filename);
    std::vector<Point_3> points;
    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::vector<double> row = StringSplit(line, ',');
        points.emplace_back(row[1], row[2], row[3]);
    }
    file.close();

    return points;
}

std::vector<double> StringSplit(std::string& str, const char delimiter) {
    str.erase(std::remove_if(str.begin(), str.end(), isspace), str.end());
    str = str + delimiter;
    size_t pos = 0;
    std::vector<double> res;

    std::string token;
    while ((pos = str.find(delimiter)) != std::string::npos) {
        res.emplace_back(std::stod(str.substr(0, pos)));
        str.erase(0, pos + 1);
    }

    return res;
}
