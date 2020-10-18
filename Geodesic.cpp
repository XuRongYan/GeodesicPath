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
	if (mesh->is_boundary(b)) return computeAlphaAndOuterArcOnBound(mesh, a, b, c);
	bool flag = false;
	double alpha1 = 0, alpha2 = 0;
	std::vector<SurfaceMesh::Vertex> side1, side2;
	const auto &pb = mesh->position(b);
	side1.emplace_back(a);
	side2.emplace_back(c);
	// use halfedge b->a as start edge.
	SurfaceMesh::Halfedge start = mesh->find_halfedge(b, a);
	SurfaceMesh::Halfedge next = mesh->ccw_rotated_halfedge(start);
	SurfaceMesh::Halfedge tower = start;
	while (next != tower) {
		// traverse halfedge counterclockwise.
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
		next = mesh->ccw_rotated_halfedge(start);
	}

	const auto &p1 = mesh->position(a_);
	const auto &p2 = mesh->position(side2.back());
	side2.emplace_back(a_);
	alpha2 += angle(p1 - pb, p2 - pb);
	std::reverse(side2.begin(), side2.end());

	return alpha1 < alpha2 ? std::make_pair(alpha1, side1) : std::make_pair(alpha2, side2);
}

std::pair<double, std::vector<SurfaceMesh::Vertex>> Joint::computeAlphaAndOuterArcOnBound(const SurfaceMesh *mesh,
																						  const SurfaceMesh::Vertex &a,
																						  const SurfaceMesh::Vertex &b,
																						  const SurfaceMesh::Vertex &c) {
	assert(mesh->is_boundary(b));
	auto is_path = mesh->get_edge_property<bool>("e:is_path");
	const Point &p1 = mesh->position(a);
	const Point &p2 = mesh->position(b);
	const Point &p3 = mesh->position(c);
	double alpha = angle(p1 - p2, p3 - p2);
	std::vector<SurfaceMesh::Vertex> outer_arc;
	auto start = mesh->find_halfedge(b, a);
	auto next = mesh->ccw_rotated_halfedge(start);
	// 1. 若逆时针旋转起始就遇到两个boundary
	if (mesh->is_boundary(mesh->edge(start)) && mesh->is_boundary(mesh->edge(next))) {
		start = mesh->find_halfedge(b, c);
		next = mesh->ccw_rotated_halfedge(start);
	}
	auto tower = start;
	auto end = mesh->to_vertex(tower) == a ? c : a;
	outer_arc.emplace_back(mesh->to_vertex(start));
	while (next != tower) {
		SurfaceMesh::Vertex v1 = mesh->to_vertex(start);
		SurfaceMesh::Vertex v2 = mesh->to_vertex(next);
		outer_arc.emplace_back(v2);
		if (v2 == end) break;
		start = next;
		next = mesh->ccw_rotated_halfedge(start);
	}

	if (outer_arc[0] == c) std::reverse(outer_arc.begin(), outer_arc.end());
	return {alpha, outer_arc};
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
	computeFlippableEdges();
	updateFlexibleState();

}

/**
 * 条件1：alpha < pi
 * 条件2：outerArc中的顶点数量小于3
 * 条件3：入射点之中是否有别的路径存在
 * 条件4：当前joint是否为最左的joint
 */
void Joint::updateFlexibleState() {
	flexible_ = alpha_ < M_PI - 1e-6 && outer_arc_.size() >= 3 && !flippable_edges_.empty();
	if (!flexible_) return;
	auto is_path = mesh_->get_edge_property<bool>("e:is_path");
	if (is_path) {
		for (auto v : outer_arc_) {
			if (v == a_ || v == c_) continue;
			if (is_path[mesh_->find_edge(v, b_)]) flexible_ = false;
		}
	}
}

double Joint::updateAlpha() {
	double sum = 0;
	if (outer_arc_.size() < 3) return M_PI;
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

bool Joint::empty() {
	return mesh_ == nullptr;
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

void Joint::setPathIdx(size_t pathIdx) {
	path_idx_ = pathIdx;
}

std::ostream &operator<<(std::ostream &os, const Joint &joint) {
	os << "path idx:" << joint.path_idx_ << " a_: " << joint.a_ << " b_: " << joint.b_ << " c_: " << joint.c_ << " flexible_: " << joint.flexible_
	   << " alpha_: " << joint.alpha_;
	return os;
}



//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

Geodesic::Geodesic(const SurfaceMesh &mesh, const std::vector<SurfaceMesh::Vertex> &path) : mesh_(mesh), path_(path) {
	init();
}

double Geodesic::makePathGeodesic() {
	// std::cout << "compute path..." << std::endl;
	size_t iter = 0;
	// std::cout << "iter = " << iter << std::endl;
	while (!joints_.empty()) {
		auto joint = joints_.top();
		joints_.pop();
		if (!joint->isFlexible()) continue;
		Joint joint_shorter = flipOut(*joint);
		mesh_.write("flippedMesh.obj");
		updatePath(joint_shorter);
		iter++;
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
	SurfaceMesh::Edge_property<bool> is_path = mesh_.get_edge_property<bool>("e:is_path");
	if (!is_path) {
		is_path = mesh_.add_edge_property<bool>("e:is_path");
	}
	for (size_t i = 0; i < path_.size() - 1; i++) {
		auto e = mesh_.find_edge(path_[i], path_[i + 1]);
		assert(mesh_.is_valid(e));
		is_path[e] = true;
	}
}

void Geodesic::computeJoints(int start, int end) {
	if (start < 1) start = 1;
	if (end > path_.size() - 1) end = path_.size() - 1;
	for (size_t i = start; i < end; i++) {
		Joint joint(&mesh_, i, path_[i - 1], path_[i], path_[i + 1]);
		joint.updateFlexibleState();
		// std::cout << joint << std::endl;
		vec_joints_[path_[i].idx()] = joint;
		// 构建双向链表
		if (joint.isFlexible()) {
			joints_.push(&vec_joints_[path_[i].idx()]);
		}
	}
	for (size_t i = end; i < path_.size(); i++) {
		vec_joints_[path_[i].idx()].setPathIdx(i);
	}
	// 有尾指针
	if (end < path_.size() && vec_joints_[end].isFlexible()) {
		vec_joints_[path_[end - 1].idx()].next_ = &vec_joints_[path_[end].idx()];
		vec_joints_[path_[end].idx()].prev_ = &vec_joints_[path_[end - 1].idx()];
	}
}

Joint Geodesic::flipOut(Joint &joint) {
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
	// std::cout << "update path: start = " << joint.getA() << " end = " << joint.getC() << std::endl;
	// 最佳的outer arc作为新的路径
	auto new_joint = joint.getOuterArc();
	int start = joint.getPathIdx(), end = start + 2;
	if (start - 1 >= 0)
		vec_joints_[start - 1].setFlexible(false);
	if (start + 1 < vec_joints_.size())
		vec_joints_[start + 1].setFlexible(false);
	std::vector<SurfaceMesh::Vertex> new_path;
	new_path.reserve(path_.size() + new_joint.size() - 3);
	for (size_t i = 0; i < start - 1; i++)
		new_path.emplace_back(path_[i]);
	size_t new_end = new_path.size();
	for (const auto &v : new_joint)
		new_path.emplace_back(v);
	for (size_t i = end; i < path_.size(); i++)
		new_path.emplace_back(path_[i]);
	path_ = new_path;
//	for (auto v : path_) {
//		std::cout << v.idx() << " ";
//	}
//	std::cout << std::endl;
	updatePathState(joint);
	start -= 2;
	computeJoints(start, start + new_joint.size() + 2);
}

void Geodesic::updatePathState(const Joint &joint) {
	auto is_path = mesh_.get_edge_property<bool>("e:is_path");
	SurfaceMesh::Edge e1, e2;
	e1 = mesh_.find_edge(joint.getA(), joint.getB());
	e2 = mesh_.find_edge(joint.getB(), joint.getC());
	assert(mesh_.is_valid(e1));
	is_path[mesh_.find_edge(joint.getA(), joint.getB())] = false;
	assert(mesh_.is_valid(e2));
	is_path[mesh_.find_edge(joint.getB(), joint.getC())] = false;
	indexEdges();
}

const SurfaceMesh &Geodesic::getMesh() const {
	return mesh_;
}

const std::priority_queue<Joint*, std::vector<Joint*>, cmp> &Geodesic::getJoints() const {
	return joints_;
}

const std::vector<Joint> &Geodesic::getVecJoints() const {
	return vec_joints_;
}

const std::vector<SurfaceMesh::Vertex> &Geodesic::getPath() const {
	return path_;
}
