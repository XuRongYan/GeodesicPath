//
// Created by 徐溶延 on 2020/10/12.
//

#ifndef GEODESICPATH_GEODESIC_H
#define GEODESICPATH_GEODESIC_H

#include <queue>
#include <unordered_map>

#include <SurfaceMesh.h>
#include <Eigen.h>
#include <ostream>

using namespace Surface_Mesh;



class Joint {
public:
	Joint();

	Joint(SurfaceMesh *mesh, size_t path_idx, const SurfaceMesh::Vertex &a, const SurfaceMesh::Vertex &b, const SurfaceMesh::Vertex &c);

	friend std::ostream &operator<<(std::ostream &os, const Joint &joint);

	SurfaceMesh *getMesh() const;

	void deleteArcPoint(size_t i);

	void updateOuterArc();

	double updateAlpha();

	void updateFlexibleState();

	size_t getPathIdx() const;

	void setPathIdx(size_t pathIdx);

	const SurfaceMesh::Vertex &getA() const;

	const SurfaceMesh::Vertex &getB() const;

	const SurfaceMesh::Vertex &getC() const;

	double getAlpha() const;

	bool empty();

	bool isFlexible() const;

	void setFlexible(bool flexible);

	bool isFirst() const;

	void setFirst(bool first);

	const std::vector<bool> &getDeletedArc() const;

	const std::vector<double> &getBetas() const;

	const std::vector<SurfaceMesh::Vertex> &getOuterArc() const;

	const std::vector<std::pair<size_t, SurfaceMesh::Edge>> &getFlippableEdges() const;

	bool operator<(const Joint &rhs) const;

	bool operator>(const Joint &rhs) const;

	bool operator<=(const Joint &rhs) const;

	bool operator>=(const Joint &rhs) const;

	std::pair<double, std::vector<SurfaceMesh::Vertex>> computeAlphaAndOuterArcOnBound(const SurfaceMesh *mesh,
																					   const SurfaceMesh::Vertex &a,
																					   const SurfaceMesh::Vertex &b,
																					   const SurfaceMesh::Vertex &c);

	std::pair<double, std::vector<SurfaceMesh::Vertex>> computeAlphaAndOuterArc(const SurfaceMesh *mesh,
																				const SurfaceMesh::Vertex &a,
																				const SurfaceMesh::Vertex &b,
																				const SurfaceMesh::Vertex &c);


	void computeFlippableEdges();

private:
	void init();

private:
	SurfaceMesh *mesh_;
	SurfaceMesh::Vertex a_, b_, c_;
	size_t path_idx_{0};
	bool flexible_{false};
	bool first_{true};
	double alpha_{0.0};
	std::vector<double> betas_;
	std::vector<SurfaceMesh::Vertex> outer_arc_;
	std::vector<std::pair<size_t, SurfaceMesh::Edge>> flippable_edges_;
	std::vector<bool> deleted_arc_;
};

struct cmp {
	bool operator()(Joint *&a, Joint *&b) const
	{
		return a->getAlpha() > b->getAlpha();
	}
};

class Geodesic {
public:
	Geodesic(const SurfaceMesh &mesh, const std::vector<SurfaceMesh::Vertex> &path);

	double makePathGeodesic();

	const SurfaceMesh &getMesh() const;

	const std::priority_queue<Joint*, std::vector<Joint*>, cmp> &getJoints() const;

	const std::vector<Joint> &getVecJoints() const;

	const std::vector<SurfaceMesh::Vertex> &getPath() const;

	Joint flipOut(Joint& joint);

	void updatePath(const Joint &joint);

	void updatePathState(const Joint &joint);

private:
	void init();

	void indexEdges();

	void computeJoints(int start, int end);



private:
	SurfaceMesh mesh_;
	std::priority_queue<Joint*, std::vector<Joint*>, cmp> joints_;
	std::vector<Joint> vec_joints_;
	std::vector<SurfaceMesh::Vertex> path_;
};


#endif //GEODESICPATH_GEODESIC_H
