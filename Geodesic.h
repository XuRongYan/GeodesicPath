//
// Created by 徐溶延 on 2020/10/12.
//

#ifndef GEODESICPATH_GEODESIC_H
#define GEODESICPATH_GEODESIC_H

#include <queue>
#include <unordered_map>

#include <SurfaceMesh.h>
#include <Eigen.h>

using namespace Surface_Mesh;

class Joint {
public:
	Joint();

	Joint(SurfaceMesh *mesh, size_t path_idx, const SurfaceMesh::Vertex &a, const SurfaceMesh::Vertex &b, const SurfaceMesh::Vertex &c);

	SurfaceMesh *getMesh() const;

	void deleteArcPoint(size_t i);

	void updateOuterArc();

	double updateAlpha();

	void updateFlexibleState();

	size_t getPathIdx() const;

	const SurfaceMesh::Vertex &getA() const;

	const SurfaceMesh::Vertex &getB() const;

	const SurfaceMesh::Vertex &getC() const;

	double getAlpha() const;

	bool isFlexible() const;

	void setFlexible(bool flexible);

	const std::vector<bool> &getDeletedArc() const;

	const std::vector<double> &getBetas() const;

	const std::vector<SurfaceMesh::Vertex> &getOuterArc() const;

	const std::vector<std::pair<size_t, SurfaceMesh::Edge>> &getFlippableEdges() const;

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

public:
	Joint *prev_{nullptr}, *next_{nullptr};		// 将Joint维护为一个双向链表

private:
	SurfaceMesh *mesh_;
	SurfaceMesh::Vertex a_, b_, c_;
	size_t path_idx_{0};
	bool flexible_{false};
	double alpha_{0.0};
	std::vector<double> betas_;
	std::vector<SurfaceMesh::Vertex> outer_arc_;
	std::vector<std::pair<size_t, SurfaceMesh::Edge>> flippable_edges_;
	std::vector<bool> deleted_arc_;
};

class Geodesic {
public:
	Geodesic(const SurfaceMesh &mesh, const std::vector<SurfaceMesh::Vertex> &path);

	double makePathGeodesic();

private:
	void init();

	void indexEdges();

	void computeJoints(size_t start, size_t end);

	Joint flipOut(Joint& joint);

	void updatePath(const Joint &joint);

private:
	SurfaceMesh mesh_;
	std::priority_queue<Joint*> joints_;
	std::vector<Joint> vec_joints_;
	std::vector<SurfaceMesh::Vertex> path_;
};


#endif //GEODESICPATH_GEODESIC_H
