//
// Created by 徐溶延 on 2020/10/12.
//

#include "Geodesic.h"
#include "MathUtils.h"


Joint::Joint(SurfaceMesh *mesh, const SurfaceMesh::Vertex &a, const SurfaceMesh::Vertex &b,
			 const SurfaceMesh::Vertex &c) : mesh_(mesh), a_(a), b_(b), c_(c) {
	assert(mesh_);
	init();
}

void Joint::init() {
	auto[alpha, outer_arc] = computeAlphaAndOuterArc(mesh_, a_, b_, c_);
	alpha_ = alpha;
	outer_arc_ = outer_arc;
	deleted_arc_.resize(outer_arc_.size());
	updateFlexibleState();
	computeFlippableEdges();
}

std::pair<double, std::vector<SurfaceMesh::Vertex>> Joint::computeAlphaAndOuterArc(const SurfaceMesh *mesh,
																				   const SurfaceMesh::Vertex &a,
																				   const SurfaceMesh::Vertex &b,
																				   const SurfaceMesh::Vertex &c) {
	bool flag = false;
	double alpha1 = 0, alpha2 = 0;
	std::vector<SurfaceMesh::Vertex> side1, side2;
	const auto &pb = mesh->position(b);
	side1.emplace_back(a);
	side2.emplace_back(c);
	// use halfedge b->a as start edge.
	SurfaceMesh::Halfedge start = mesh->find_halfedge(b, a);
	SurfaceMesh::Halfedge next{};
	while (next != start) {
		// traverse halfedge counterclockwise.
		next = mesh->ccw_rotated_halfedge(start);
		SurfaceMesh::Vertex v1 = mesh->to_vertex(start);
		SurfaceMesh::Vertex v2 = mesh->to_vertex(next);
		assert(v1 != v2);
		const auto &p1 = mesh->position(v1);
		const auto &p2 = mesh->position(v2);
		double tmp = angle(p1 - pb, p2 - pb);
		if (!flag) {
			alpha1 += tmp;
			side1.emplace_back(v2);
		} else {
			alpha2 += tmp;
			side2.emplace_back(v2);
		}
		// if come across c, switch the angle counter and outer arc path.
		if (v2 == c) flag = true;
		start = next;
	}
	return alpha1 < alpha2 ? std::make_pair(alpha1, side1) : std::make_pair(alpha2, side2);
}

void Joint::computeFlippableEdges() {
	flippable_edges_.clear();
	for (int i = 1; i < outer_arc_.size() - 1; i++) {
		Point p1, p2, p3, pb;
		pb = mesh_->position(b_);
		p1 = mesh_->position(outer_arc_[i - 1]);
		p2 = mesh_->position(outer_arc_[i]);
		p3 = mesh_->position(outer_arc_[i + 1]);
		double beta = 0;
		beta += angle(p1 - p2, pb - p2);
		beta += angle(p3 - p2, pb - p2);
		if (beta < M_PI) flippable_edges_.emplace_back(mesh_->find_edge(outer_arc_[i], b_));
	}
}

void Joint::deleteArcPoint(size_t i) {
	assert(deleted_arc_.size() > i);
	deleted_arc_[i] = true;
}

void Joint::updateOuterArc() {
	std::vector<SurfaceMesh::Vertex> outer_arc;
	outer_arc.reserve(outer_arc_.size());
	for (int i = 0; i < outer_arc_.size(); i++) {
		if (!deleted_arc_[i]) {
			outer_arc.emplace_back(outer_arc_[i]);
		}
	}
	outer_arc_ = outer_arc;
	deleted_arc_.clear();
	deleted_arc_.resize(outer_arc_.size());
	alpha_ = updateAlpha();
	updateFlexibleState();
	computeFlippableEdges();
}

void Joint::updateFlexibleState() {
	flexible_ = alpha_ < M_PI || outer_arc_.size() > 3;
}

SurfaceMesh *Joint::getMesh() const {
	return mesh_;
}

const SurfaceMesh::Vertex &Joint::getA() const {
	return a_;
}

const SurfaceMesh::Vertex &Joint::getB() const {
	return b_;
}

const SurfaceMesh::Vertex &Joint::getC() const {
	return c_;
}

double Joint::getAlpha() const {
	return alpha_;
}

bool Joint::isFlexible() const {
	return flexible_;
}

const std::vector<bool> &Joint::getDeletedArc() const {
	return deleted_arc_;
}

const std::vector<double> &Joint::getBetas() const {
	return betas_;
}

const std::vector<SurfaceMesh::Vertex> &Joint::getOuterArc() const {
	return outer_arc_;
}

const std::vector<SurfaceMesh::Edge> &Joint::getFlippableEdges() const {
	return flippable_edges_;
}

bool Joint::operator<(const Joint &rhs) const {
	return alpha_ < rhs.alpha_;
}

bool Joint::operator>(const Joint &rhs) const {
	return rhs < *this;
}

bool Joint::operator<=(const Joint &rhs) const {
	return !(rhs < *this);
}

bool Joint::operator>=(const Joint &rhs) const {
	return !(*this < rhs);
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

Geodesic::Geodesic(const SurfaceMesh &mesh, const std::vector<SurfaceMesh::Vertex> &path) : mesh_(mesh), path_(path) {
	init();
}

void Geodesic::init() {
	computeJoints();
}

void Geodesic::computeJoints() {
	for (size_t i = 1; i < path_.size() - 1; i++) {
		Joint joint(Joint(&mesh_, path_[i - 1], path_[i], path_[i + 1]));
		if (joint.getAlpha() < M_PI)
			joints_.push(joint);
	}
}


double Geodesic::makePathGeodesic() {
	return 0;
}

Joint Geodesic::flipOut(const Joint &joint) {

	return Joint(nullptr, SurfaceMesh::Vertex(), SurfaceMesh::Vertex(), SurfaceMesh::Vertex());
}
