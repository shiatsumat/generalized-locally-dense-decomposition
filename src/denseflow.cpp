#include <algorithm>
#include <array>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <unordered_map>
#include "flow.h"
using namespace std;

bool SORT = true;

void solve_body(vector<flowedge_t>& flowedges,
  vector<pair<vid_t, vid_t>>& deg, vector<vid_t>& dict,
  vector<vector<vid_t>>& layers,
  int call_depth, int& max_call_depth) {
  max_call_depth = max(max_call_depth, call_depth);
  vid_t n_vertices = (vid_t)dict.size();
  cout << string(call_depth * 2, ' ')
    << "solve: #vertices " << n_vertices << endl;
  flowadj_t flowadj;
  make_flowadj(n_vertices, flowedges, flowadj);

  #ifdef DEBUG
  flow_t flow =
  #endif
  dinic(flowadj);

  vector<int> side(n_vertices + 2, 1);
  vector<vid_t> que(n_vertices + 2);
  vid_t s = n_vertices;
  vid_t l = 0, r = 0; que.AT(r++) = s; side.AT(s) = 0;
  array<vid_t, 2> n_subvertices = {{0, 0}};
  while(l != r) {
    vid_t v = que.AT(l++);
    for(const flowrefedge_t& flowrefedge : flowadj.AT(v)) {
      vid_t w = flowrefedge.to();
      if(side.AT(w) == 0 || flowrefedge.residue() == 0) continue;
      que.AT(r++) = w; side.AT(w) = 0; ++n_subvertices.AT(0);
    }
  }
  que.clear();
  n_subvertices.AT(1) = n_vertices - n_subvertices.AT(0);
  #ifdef DEBUG
    assert(n_subvertices.AT(1) > 0);
    flow_t cut = 0;
    for(const flowedge_t& flowedge : flowedges) {
      int sd = side.AT(flowedge.from), sd2 = side.AT(flowedge.to);
      if(sd < sd2) cut += flowedge.flow;
      else if(sd > sd2) cut += -flowedge.flow;
    }
    assert(cut == flow);
  #endif

  if(n_subvertices.AT(0) == 0) {
    #ifdef DEBUG
      vid_t t = n_vertices + 1;
      for(const flowrefedge_t& flowrefedge : flowadj.AT(t)) {
        assert(flowrefedge.back_residue() == 0);
      }
    #endif
    layers.emplace_back(); dict.swap(layers.back());
    flowedges.clear(); deg.clear();
    return;
  }
  flowadj.clear();
  cout << string(call_depth * 2, ' ')
    << " divided into size " << n_subvertices.AT(0)
    << " and size " << n_subvertices.AT(1) << endl;

  array<vid_t, 2> count = {{0, 0}};
  vector<vid_t> idict(n_vertices);
  array<vector<vid_t>, 2> subdict;
  if(SORT) {
    for(int vx = 0; vx < n_vertices; ++vx) {
      idict.AT(vx) = count.AT(side.AT(vx))++;
    }
    array<vector<pair<vid_t, vid_t>>, 2> minideg;
    for(int vx = 0; vx < n_vertices; ++vx) {
      minideg.AT(side.AT(vx)).emplace_back(0, -vx);
    }
    for(flowedge_t& flowedge : flowedges) {
      vid_t fromx = flowedge.from; if(fromx == s) break;
      vid_t tox = flowedge.to;
      int sd = side.AT(fromx); if(sd != side.AT(tox)) continue;
      ++minideg.AT(sd).AT(idict.AT(fromx)).first;
      ++minideg.AT(sd).AT(idict.AT(tox)).first;
    }
    for(int sd = 0; sd < 2; ++sd) {
      sort(minideg.AT(sd).begin(), minideg.AT(sd).end(),
        greater<pair<vid_t, vid_t>>());
      subdict.AT(sd).resize(n_subvertices.AT(sd));
      for(vid_t v = 0; v < (int)minideg.AT(sd).size(); ++v) {
        vid_t vx = -minideg.AT(sd).AT(v).second;
        idict.AT(vx) = v; subdict.AT(sd).AT(v) = vx;
      }
      minideg.AT(sd).clear();
    }
  } else {
    for(int vx = 0; vx < n_vertices; ++vx) {
      int sd = side.AT(vx);
      idict.AT(vx) = (vid_t)subdict.AT(sd).size();
      subdict.AT(sd).push_back(vx);
    }
  }

  array<eid_t, 2> n_suboutedges = {{0, 0}};
  array<vector<flowedge_t>, 2> flowsubedges;
  for(flowedge_t& flowedge : flowedges) {
    vid_t fromx = flowedge.from; if(fromx == s) break;
    vid_t tox = flowedge.to;
    vid_t from = idict.AT(fromx), to = idict.AT(tox);
    int sd = side.AT(fromx), sd2 = side.AT(tox);
    if(sd == sd2) {
      ++n_suboutedges.AT(sd);
      if(from > to) swap(from, to);
      flowsubedges.AT(sd).emplace_back(
        from, to, n_subvertices.AT(sd), n_subvertices.AT(sd));
    } else {
      if(sd > sd2) swap(fromx, tox);
      --deg.AT(fromx).second; ++deg.AT(tox).first;
    }
  }
  idict.clear(); side.clear(); flowedges.clear();

  array<vector<pair<vid_t, vid_t>>, 2> subdeg;
  for(int sd = 0; sd < 2; ++sd) {
    if(SORT) {
      sort(flowsubedges.AT(sd).begin(), flowsubedges.AT(sd).end());
    }
    subdeg.AT(sd).resize(n_subvertices.AT(sd));
    for(vid_t v = 0; v < n_subvertices.AT(sd); ++v) {
      subdeg.AT(sd).AT(v) = deg.AT(subdict.AT(sd).AT(v));
      n_suboutedges.AT(sd) += subdeg.AT(sd).AT(v).first;
      subdict.AT(sd).AT(v) = dict.AT(subdict.AT(sd).AT(v));
    }
    for(vid_t v = 0; v < n_subvertices.AT(sd); ++v) {
      flowsubedges.AT(sd).emplace_back(
        n_subvertices.AT(sd), v,
        (flow_t) n_subvertices.AT(sd) *
          (subdeg.AT(sd).AT(v).first
          + subdeg.AT(sd).AT(v).second));
      flowsubedges.AT(sd).emplace_back(
        v, n_subvertices.AT(sd) + 1,
        2 * n_suboutedges.AT(sd));
    }
  }
  deg.clear(); dict.clear();
  for(int sd = 1; sd >= 0; --sd) {
    solve_body(
      flowsubedges.AT(sd), subdeg.AT(sd), subdict.AT(sd), layers,
      call_depth + 1, max_call_depth);
  }
}

struct edge_t {
  vid_t from, to;
  inline edge_t(vid_t from, vid_t to) : from(from), to(to) {}
  inline bool operator<(const edge_t& that) const {
    return this->from < that.from
      || (this->from == that.from && this->to < that.to);
  }
};

void solve(vid_t n_vertices, vector<edge_t>& edges,
  vector<vector<vid_t>>& layers, int& max_call_depth) {
  vid_t s = n_vertices, t = n_vertices + 1;
  eid_t n_edges = edges.size();
  vector<pair<vid_t, vid_t>> deg(n_vertices);
  vector<vid_t> dict(n_vertices);
  vector<flowedge_t> flowedges;
  for(const edge_t& edge : edges) {
    flowedges.emplace_back(
      edge.from, edge.to, n_vertices, n_vertices);
    ++deg.AT(edge.from).second; ++deg.AT(edge.to).second;
  }
  for(int v = 0; v < n_vertices; ++v) {
    flowedges.emplace_back(s, v,
      (eid_t) n_vertices * deg.AT(v).second);
    flowedges.emplace_back(v, t,
      2 * n_edges);
  }
  for(int v = 0; v < n_vertices; ++v) dict.AT(v) = v;
  max_call_depth = 0;
  solve_body(flowedges, deg, dict, layers, 0, max_call_depth);
  int n_layers = (int)layers.size();
  for(int i = 0; i < n_layers - 1 - i; ++i) {
    layers.AT(i).swap(layers.AT(n_layers - 1 - i));
  }
}

void sort_graph(vid_t n_vertices, vector<edge_t>& edges) {
  vector<pair<vid_t, vid_t>> predeg(n_vertices);
  vector<vid_t> idict(n_vertices);
  for(vid_t v = 0; v < n_vertices; ++v) predeg.AT(v) = {0, -v};
  for(const edge_t& edge : edges) {
    ++predeg.AT(edge.from).first; ++predeg.AT(edge.to).first;
  }
  sort(predeg.begin(), predeg.end(),
    greater<pair<vid_t, vid_t>>());
  for(vid_t v = 0; v < n_vertices; ++v) {
    idict.AT(-predeg.AT(v).second) = v;
  }
  for(edge_t& edge : edges) {
    edge.from = idict.AT(edge.from); edge.to = idict.AT(edge.to);
    if(edge.from > edge.to) swap(edge.from, edge.to);
  }
  sort(edges.begin(), edges.end());
}

vid_t get_vid(vid_t vidx,
  vid_t& n_vertices, unordered_map<vid_t, vid_t>& idict) {
  if(idict.find(vidx) != idict.end()) {
    return idict.AT(vidx);
  }
  vid_t res = idict[vidx] = n_vertices; ++n_vertices;
  return res;
}
void load_graph(basic_istream<char>& istrm,
  vid_t& n_vertices, vector<edge_t>& edges) {
  unordered_map<vid_t, vid_t> idict;
  n_vertices = 0; edges.clear(); string line;
  while(getline(istrm, line)) {
    if(line.AT(0) == '#') continue;
    vid_t fromx, tox;
    sscanf(line.c_str(), "%d%d", &fromx, &tox);
    edges.emplace_back(
      get_vid(fromx, n_vertices, idict),
      get_vid(tox, n_vertices, idict));
  }
}

int count_digits(size_t a) {
  int d = 0;
  for(size_t b = 1; b <= a; ++d, b *= 10) {}
  return d;
}
int main(int argc, char* argv[]) {
  ifstream in;
  for(int i = 1; i < argc; ++i) {
    string arg = argv[i];
    if(arg == "--no-sort") {
      SORT = false;
      cout << "no sort" << endl;
    } else {
      in.open(arg);
      if(in.fail()) {
        cout << "Non-existent file: " << arg << endl;
        exit(EXIT_FAILURE);
      }
      cout << "opened " << arg << endl;
      cin.rdbuf(in.rdbuf());
    }
  }

  vid_t n_vertices; vector<edge_t> edges;
  cout << "Loading graph ..." << flush;
  load_graph(cin, n_vertices, edges);
  if(n_vertices == 0) {
    cout << endl << "Empty graph!" << endl;
    exit(EXIT_FAILURE);
  }
  if(SORT) sort_graph(n_vertices, edges);
  cout << "\rLoaded graph: "
    << n_vertices << " vertices, "
    << edges.size() << " edges" << endl;

  cout << "Computing ..." << endl;
  auto start = chrono::system_clock::now();
  vector<vector<vid_t>> layers;
  int max_call_depth;
  solve(n_vertices, edges, layers, max_call_depth);
  auto finish = chrono::system_clock::now();
  cout << "Computation finished!" << endl;
  int n_layers = (int)layers.size();
  cout << "#layers: " << n_layers << endl;

  unordered_map<vid_t, int> layer_dict;
  vector<vid_t> layer_size; vector<eid_t> layer_mass;
  for(int l = 0; l < n_layers; ++l) {
    layer_size.push_back((int)layers.AT(l).size());
    for(int vx : layers.AT(l)) {
      #ifdef DEBUG
        assert(layer_dict.find(vx) == layer_dict.end());
      #endif
      layer_dict[vx] = l;
    }
  }
  layer_mass.resize(layers.size());
  for(const edge_t& edge : edges) {
    ++layer_mass[max(
      layer_dict.AT(edge.from), layer_dict.AT(edge.to))];
  }
  int digit_lid = count_digits(n_layers);
  int digit_mass = 0, digit_size = 0;
  for(int l = 0; l < n_layers; ++l) {
    digit_mass = max(digit_mass, count_digits(layer_mass.AT(l)));
    digit_size = max(digit_size, count_digits(layer_size.AT(l)));
  }
  for(int l = 0; l < n_layers; ++l) {
    cout << "layer " << setw(digit_lid) << l+1
      << ": outer mass " << setw(digit_mass) << layer_mass.AT(l)
      << " outer size " << setw(digit_size) << layer_size.AT(l)
      << " outer density "
      << layer_mass.AT(l) * 1.0 / layer_size.AT(l) << endl;
  }
  cout << "Max call depth: " << max_call_depth << endl;
  cout << "Elapsed time: " <<
    chrono::duration_cast<chrono::milliseconds>
      (finish - start).count() / 1e3 << " s" << endl;
  return 0;
}
