//
// Created by 徐溶延 on 2020/10/12.
//

#include "Geodesic.h"
#include "MathUtils.h"


Joint::Joint() {}

Joint::Joint(SurfaceMesh *mesh, size_t path_idx, const SurfaceMesh::Vertex &a, const SurfaceMesh::Vertex &b,
			 const SurfaceMesh::Vertex &c) : mesh_(mesh), path_idx_(path_idx), a_(a), b_(b), c_(c) {
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
		if (beta < M_PI) flippable_edges_.emplace_back(i, mesh_->find_edge(outer_arc_[i], b_));
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

/**
 * 条件1：alpha < pi
 * 条件2：outerArc中的顶点数量小于3
 * 条件3：入射点之中是否有别的路径存在
 * 条件4：当前joint是否为最左的joint
 */
void Joint::updateFlexibleState() {
	flexible_ = alpha_ < M_PI || outer_arc_.size() >= 3;
	if (!flexible_) return;
	auto is_path = mesh_->get_edge_property<bool>("e:is_path");
	for (auto v : outer_arc_) {
		if (v == a_ || v == c_) continue;
		if (is_path[mesh_->find_edge(v, b_)]) flexible_ = false;
	}
}

double Joint::updateAlpha() {
	double sum = 0;
	for (size_t i = 0; i < outer_arc_.size() - 1; i++) {
		const Point &p1 = mesh_->position(outer_arc_[i]);
		const Point &p2 = mesh_->position(outer_arc_[i + 1]);
		const Point &pb = mesh_->position(b_);
		sum += angle(p1 - pb, p2 - pb);
	}
	return sum;
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

void Joint::setFlexible(bool flexible) {
	flexible_ = flexible;
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

const std::vector<std::pair<size_t, SurfaceMesh::Edge>> &Joint::getFlippableEdges() const {
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

size_t Joint::getPathIdx() const {
	return path_idx_;
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

Geodesic::Geodesic(const SurfaceMesh &mesh, const std::vector<SurfaceMesh::Vertex> &path) : mesh_(mesh), path_(path) {
	init();
}

double Geodesic::makePathGeodesic() {
	while (!joints_.empty()) {
		auto joint = joints_.top();
		joints_.pop();
		if (!joint->isFlexible()) continue;
		Joint joint_shorter = flipOut(*joint);
		updatePath(*joint);
	}
	double dist = 0;
	for (size_t i = 0; i < path_.size() - 1; i++) {
		const Point &p1 = mesh_.position(path_[i]);
		const Point &p2 = mesh_.position(path_[i + 1]);
		dist += (p1 - p2).norm();
	}
	return dist;
}

void Geodesic::init() {
	vec_joints_.resize(mesh_.n_vertices());
	indexEdges();
	computeJoints(1, path_.size() - 1);
}

void Geodesic::indexEdges() {
	auto is_path = mesh_.add_edge_property<bool>("e:is_path");
	for (size_t i = 0; i < path_.size() - 1; i++) {
		is_path[mesh_.find_edge(path_[i], path_[i + 1])] = true;
	}
}

void Geodesic::computeJoints(size_t start, size_t end) {
	for (size_t i = start; i < end; i++) {
		Joint joint(&mesh_, i, path_[i - 1], path_[i], path_[i + 1]);
		vec_joints_[path_[i].idx()] = joint;
		// 构建双向链表
		if (i > 1) {
			vec_joints_[path_[i].idx()].prev_ = &vec_joints_[path_[i - 1].idx()];
			vec_joints_[path_[i - 1].idx()].next_ = &vec_joints_[path_[i].idx()];
		}

		if (joint.isFlexible()) {
			joints_.push(&vec_joints_[path_[i].idx()]);
		}
	}
	// 有尾指针
	if (end < path_.size() && vec_joints_[end].isFlexible()) {
		vec_joints_[path_[end - 1].idx()].next_ = &vec_joints_[path_[end].idx()];
		vec_joints_[path_[end].idx()].prev_ = &vec_joints_[path_[end - 1].idx()];
	}
}

Joint Geodesic::flipOut(Joint &joint) {
	assert(joint.isFlexible());
	while (joint.isFlexible()) {
		auto edges = joint.getFlippableEdges();
		for (const auto &p : edges) {
			assert(mesh_.is_flip_ok(p.second));
			mesh_.flip(p.second);
			joint.deleteArcPoint(p.first);
		}
		joint.updateOuterArc();
	}
	return joint;
}

void Geodesic::updatePath(const Joint &joint) {
	// 最佳的outer arc作为新的路径
	auto new_joint = joint.getOuterArc();
	size_t start = joint.getPathIdx(), end = start;
	if (joint.prev_) {
		joint.prev_->setFlexible(false);
		start = joint.prev_->getPathIdx();
	}
	if (joint.next_) {
		joint.next_->setFlexible(false);
		end = joint.next_->getPathIdx();
	}
	std::vector<SurfaceMesh::Vertex> new_path;
	new_path.reserve(path_.size() + new_joint.size() - 3);
	for (size_t i = 0; i < start; i++)
		new_path.emplace_back(path_[i]);
	size_t new_end = new_path.size();
	for (const auto &v : new_joint)
		new_path.emplace_back(v);
	for (size_t i = end + 1; i < path_.size(); i++)
		new_path.emplace_back(path_[i]);
	path_ = new_path;
	indexEdges();
	computeJoints(start, new_end);
}
