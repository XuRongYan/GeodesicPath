//
// Created by 徐溶延 on 2020/10/12.
//

#ifndef GEODESICPATH_GEODESIC_H
#define GEODESICPATH_GEODESIC_H

#include <queue>

#include <SurfaceMesh.h>
#include <Eigen.h>

using namespace Surface_Mesh;

class Joint {
public:
	Joint(SurfaceMesh *mesh, const SurfaceMesh::Vertex &a, const SurfaceMesh::Vertex &b, const SurfaceMesh::Vertex &c);

	SurfaceMesh *getMesh() const;

	void deleteArcPoint(size_t i);

	void updateOuterArc();

	double updateAlpha();

	void updateFlexibleState();

	const SurfaceMesh::Vertex &getA() const;

	const SurfaceMesh::Vertex &getB() const;

	const SurfaceMesh::Vertex &getC() const;

	double getAlpha() const;

	bool isFlexible() const;

	const std::vector<bool> &getDeletedArc() const;

	const std::vector<double> &getBetas() const;

	const std::vector<SurfaceMesh::Vertex> &getOuterArc() const;

	const std::vector<SurfaceMesh::Edge> &getFlippableEdges() const;

	bool operator<(const Joint &rhs) const;

	bool operator>(const Joint &rhs) const;

	bool operator<=(const Joint &rhs) const;

	bool operator>=(const Joint &rhs) const;


private:
	void init();

	std::pair<double, std::vector<SurfaceMesh::Vertex>> computeAlphaAndOuterArc(const SurfaceMesh *mesh,
																				const SurfaceMesh::Vertex &a,
																				const SurfaceMesh::Vertex &b,
																				const SurfaceMesh::Vertex &c);

	void computeFlippableEdges();

private:
	SurfaceMesh *mesh_;
	SurfaceMesh::Vertex a_, b_, c_;
	bool flexible_{true};
	double alpha_{0.0};
	std::vector<double> betas_;
	std::vector<SurfaceMesh::Vertex> outer_arc_;
	std::vector<SurfaceMesh::Edge> flippable_edges_;
	std::vector<bool> deleted_arc_;
};

class Geodesic {
public:
	Geodesic(const SurfaceMesh &mesh, const std::vector<SurfaceMesh::Vertex> &path);

	double makePathGeodesic();

private:
	void init();

	void computeJoints();

	Joint flipOut(const Joint& joint);

private:
	SurfaceMesh mesh_;
	std::priority_queue<Joint> joints_;
	std::vector<SurfaceMesh::Vertex> path_;
};


#endif //GEODESICPATH_GEODESIC_H
