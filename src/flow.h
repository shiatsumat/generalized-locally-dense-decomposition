#pragma once

#include <cassert>
#include <limits>
#include <vector>
using namespace std;

#ifdef DEBUG
#  define AT(index) at(index)
#else
#  define AT(index) operator[](index)
#endif

using vid_t = int; using eid_t = long long;
using flow_t = long long;
constexpr flow_t INF_FLOW = numeric_limits<flow_t>::max();
struct flowedge_t {
  vid_t from, to;
  flow_t fw_cap, bw_cap, flow;
  inline flowedge_t(vid_t from, vid_t to,
    flow_t fw_cap, flow_t bw_cap = 0, flow_t flow = 0) :
    from(from), to(to),
    fw_cap(fw_cap), bw_cap(bw_cap), flow(flow) {
      #ifdef DEBUG
        assert(-bw_cap <= flow && flow <= fw_cap);
      #endif
    }
  inline bool operator<(const flowedge_t& that) const {
    return this->from < that.from
      || (this->from == that.from && this->to < that.to);
  }
};
struct flowrefedge_t {
  flowedge_t& flowedge;
  bool forward;
  inline flowrefedge_t(flowedge_t& flowedge, bool forward) :
    flowedge(flowedge), forward(forward) {}
  inline vid_t to() const {
    return forward ?
      flowedge.to : flowedge.from;
  }
  inline flow_t residue() const {
    return forward ?
      flowedge.fw_cap - flowedge.flow :
      flowedge.bw_cap + flowedge.flow;
  }
  inline flow_t back_residue() const {
    return forward ?
      flowedge.bw_cap + flowedge.flow :
      flowedge.fw_cap - flowedge.flow;
  }
  inline void augment(flow_t flow) {
    if(forward) flowedge.flow += flow;
    else flowedge.flow -= flow;
  }
};

using flowadj_t = vector<vector<flowrefedge_t>>;

inline bool dinic_levelize(const flowadj_t& flowadj,
  vector<vid_t>& level, vector<vid_t>& que) {
  vid_t n_vertices = (vid_t)flowadj.size() - 2;
  vid_t s = n_vertices, t = n_vertices + 1;
  fill(level.begin(), level.end(), -1);
  vid_t l = 0, r = 0; level.AT(s) = 0; que.AT(r++) = s;
  while(l != r) {
    vid_t v = que.AT(l++); if(v == t) break;
    for(const flowrefedge_t& flowrefedge : flowadj.AT(v)) {
      vid_t w = flowrefedge.to();
      flow_t residue = flowrefedge.residue();
      if(level.AT(w) != -1 || residue == 0) continue;
      level.AT(w) = level.AT(v) + 1; que.AT(r++) = w;
    }
  }
  return level.AT(t) != -1;
}
flow_t dinic_augment(vid_t v, flow_t limit,
  flowadj_t& flowadj,
  const vector<vid_t>& level, vector<vid_t>& progress) {
  vid_t n_vertices = (vid_t)flowadj.size() - 2;
  vid_t t = n_vertices + 1; if(v == t) return limit;
  flow_t res = 0;
  for(vid_t& prog = progress.AT(v);
    limit > 0 && prog < (vid_t)flowadj.AT(v).size(); ++prog) {
    flowrefedge_t& flowrefedge = flowadj.AT(v).AT(prog);
    vid_t w = flowrefedge.to();
    flow_t residue = flowrefedge.residue();
    if(residue == 0 || level.AT(v) >= level.AT(w)) continue;
    flow_t flow = dinic_augment(w, min(limit, residue),
      flowadj, level, progress);
    flowrefedge.augment(flow); res += flow; limit -= flow;
  }
  return res;
}
flow_t dinic(flowadj_t& flowadj) {
  vid_t n_vertices = (vid_t)flowadj.size() - 2;
  vector<vid_t> level(n_vertices + 2), progress(n_vertices + 2);
  vector<vid_t> que(n_vertices + 2);
  flow_t res = 0; vid_t s = n_vertices;
  while(dinic_levelize(flowadj, level, que)) {
    fill(progress.begin(), progress.end(), 0);
    flow_t aug = dinic_augment(s, INF_FLOW,
      flowadj, level, progress);
    #ifdef DEBUG
      assert(aug > 0);
    #endif
    res += aug;
  }
  return res;
}

void hipr(flowadj_t& flowadj) {
}

inline void make_flowadj(
  vid_t n_vertices, vector<flowedge_t>& flowedges,
  flowadj_t& flowadj) {
  flowadj.resize(n_vertices + 2);
  for(flowedge_t& flowedge : flowedges) {
    flowadj.AT(flowedge.from).emplace_back(flowedge, true);
    flowadj.AT(flowedge.to).emplace_back(flowedge, false);
  }
}
