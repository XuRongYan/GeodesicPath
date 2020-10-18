#include <iostream>
#include <fstream>

#include <Eigen/Dense>
#include <SurfaceMesh/SurfaceMesh.h>

#include "vtk.h"
#include "Geodesic.h"

#include <gtest/gtest.h>

using namespace Surface_Mesh;

void path2vtk(const std::string &file_name, const SurfaceMesh& mesh, const std::vector<SurfaceMesh::Vertex>& path) {
	std::vector<double> points;
	std::vector<int> lines;
	for (const auto &v : path) {
		assert(mesh.is_valid(v));
		const auto &p = mesh.position(v);
		points.emplace_back(p[0]);
		points.emplace_back(p[1]);
		points.emplace_back(p[2]);
	}
	for (size_t i = 0; i < path.size() - 1; i++) {
		lines.emplace_back(i);
		lines.emplace_back(i + 1);
	}
	std::ofstream ofs(file_name);
	line2vtk(ofs, points.data(), points.size() / 3, lines.data(), lines.size() / 2);
}

int main(int argc, char* argv[]) {
//    ::testing::InitGoogleTest(&argc, argv);
//    return RUN_ALL_TESTS();
	// load mesh
	SurfaceMesh mesh;
	mesh.read("opt_pmesh000.obj");

	// read path
	std::ifstream ifs("opt_pmesh_path2.idx");
	std::vector<SurfaceMesh::Vertex> original_path;
	while (!ifs.eof()) {
		int id;
		ifs >> id;
		original_path.emplace_back(SurfaceMesh::Vertex(id));
	}
	path2vtk("originalPath.vtk", mesh, original_path);

	std::chrono::system_clock::time_point start = std::chrono::system_clock::now();
	Geodesic geodesic(mesh, original_path);
	double dist = geodesic.makePathGeodesic();
	std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
	std::cout << "duration: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;
	auto new_path = geodesic.getPath();
	std::cout << dist << std::endl;
	geodesic.getMesh().write("flippedMesh.obj");
	path2vtk("newPath.vtk", mesh, new_path);
	return 0;
}
