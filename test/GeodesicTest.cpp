//
// Created by 徐溶延 on 2020/10/14.
//
#include <gtest/gtest.h>
#include <fstream>
#include <vector>

#include <SurfaceMesh.h>

#include "../Geodesic.h"
#include "../MathUtils.h"

using namespace Surface_Mesh;

class GeodesicTest : public ::testing::Test {

protected:
	void SetUp() override {
		Test::SetUp();
		mesh.read("opt_pmesh000.obj");
		std::ifstream ifs("path.idx");
		while (!ifs.eof()) {
			int id;
			ifs >> id;
			path.emplace_back(SurfaceMesh::Vertex(id));
		}
	}

protected:
	SurfaceMesh mesh;
	std::vector<SurfaceMesh::Vertex> path;
};

TEST_F(GeodesicTest, JointInitTest) {
	Joint joint(&mesh, 1, path[0], path[1], path[2]);
	std::vector<SurfaceMesh::Vertex> ideal_outer_arc = {SurfaceMesh::Vertex(427), SurfaceMesh::Vertex(520),
														SurfaceMesh::Vertex(139), SurfaceMesh::Vertex(355)};
	const Point &p1 = mesh.position(path[0]);
	const Point &p2 = mesh.position(path[1]);
	const Point &p3 = mesh.position(path[2]);
	double ideal_angle = angle(p1 - p2, p3 - p2);
	const auto &outer_arc = joint.getOuterArc();
	EXPECT_LE(fabs(ideal_angle - joint.getAlpha()), 1e-3);
	EXPECT_EQ(outer_arc.size(), 4);
	EXPECT_EQ(ideal_outer_arc, outer_arc);
	EXPECT_TRUE(joint.isFlexible());
	EXPECT_EQ(2, joint.getFlippableEdges().size());
}

TEST_F(GeodesicTest, GeodesicInitTest) {
	Geodesic geodesic(mesh, path);
	auto is_path = geodesic.getMesh().get_edge_property<bool>("e:is_path");
	EXPECT_TRUE(is_path);
	for (size_t i = 0; i < path.size() - 1; i++) {
		EXPECT_TRUE(is_path[mesh.find_edge(path[i], path[i + 1])]);
	}
	Joint joint = geodesic.getVecJoints()[path[1].idx()];
	Joint* pJoint = &joint;
	EXPECT_TRUE(pJoint->isFlexible());
	int count = 1;
	while (pJoint->next_) {
		count++;
		pJoint = pJoint->next_;
	}
	EXPECT_EQ(count, path.size() - 2);
	while (pJoint->prev_) {
		count--;
		pJoint = pJoint->prev_;
	}
	EXPECT_EQ(count, 1);
	auto q = geodesic.getJoints();
	double last = std::numeric_limits<double>::min();
	while (!q.empty()) {
		auto j = q.top();
		q.pop();
		EXPECT_LT(last, j->getAlpha());
	}
}

TEST_F(GeodesicTest, GeodesicPathTest) {
	Geodesic geodesic(mesh, path);
	Joint joint = geodesic.getVecJoints()[path[1].idx()];
	auto outer_arc = joint.getOuterArc();
	Joint shorter = geodesic.flipOut(joint);
	EXPECT_EQ(shorter.getOuterArc().size(), 2);
	EXPECT_EQ(shorter.getOuterArc().front(), shorter.getA());
	EXPECT_EQ(shorter.getOuterArc().back(), shorter.getC());
	EXPECT_DOUBLE_EQ(shorter.getAlpha(), M_PI);
	EXPECT_FALSE(shorter.isFlexible());
}

