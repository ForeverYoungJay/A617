#!/usr/bin/env python3
"""Export topology features from validated grain-boundary tables.

Stdlib-only; intentionally works on CSV contracts so it can run after the
server-side ingestion/reconstruction jobs without importing a ML framework.
"""
from __future__ import annotations
import argparse, csv, math
from collections import defaultdict, deque
from pathlib import Path

def read(path):
    with path.open(newline="") as f: return list(csv.DictReader(f))

def graph(edges):
    g = defaultdict(set)
    for e in edges:
        a, b = int(e["grain_a"]), int(e["grain_b"])
        if a != b: g[a].add(b); g[b].add(a)
    return g

def components(g):
    seen, out = set(), []
    for root in g:
        if root in seen: continue
        q, seen_now = [root], set()
        while q:
            x = q.pop()
            if x in seen: continue
            seen.add(x); seen_now.add(x); q.extend(g[x] - seen)
        out.append(sorted(seen_now))
    return out

def betweenness(g):
    # Brandes for an unweighted undirected graph; ponytail: O(VE), replace with
    # a graph library only when volumes make this measurable.
    score = defaultdict(float)
    for s in g:
        stack, pred, sigma, dist = [], defaultdict(list), {s: 1.0}, {s: 0}
        q = deque([s])
        while q:
            v = q.popleft(); stack.append(v)
            for w in g[v]:
                if w not in dist: dist[w] = dist[v] + 1; q.append(w)
                if dist[w] == dist[v] + 1:
                    sigma[w] = sigma.get(w, 0) + sigma[v]; pred[w].append(v)
        delta = defaultdict(float)
        while stack:
            w = stack.pop()
            for v in pred[w]: delta[v] += sigma[v] / sigma[w] * (1 + delta[w])
            if w != s: score[w] += delta[w]
    return {v: score[v] / 2 for v in g}

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--nodes", type=Path, required=True)
    p.add_argument("--edges", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args(); args.out.mkdir(parents=True, exist_ok=True)
    nodes, edges = read(args.nodes), read(args.edges); g = graph(edges); bc = betweenness(g)
    comps = components(g); comp_id = {v: i for i, c in enumerate(comps) for v in c}
    with (args.out / "graph_nodes.csv").open("w", newline="") as f:
        fields = ["grain_id", "degree", "betweenness", "component_id"]
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader()
        for n in nodes:
            v = int(n["grain_id"]); w.writerow({"grain_id": v, "degree": len(g[v]), "betweenness": bc.get(v, 0.0), "component_id": comp_id.get(v, -1)})
    with (args.out / "graph_edges.csv").open("w", newline="") as f:
        fields = list(edges[0]) + ["component_id"] if edges else ["grain_a", "grain_b", "component_id"]
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader()
        for e in edges:
            e = dict(e); e["component_id"] = comp_id.get(int(e["grain_a"]), -1); w.writerow(e)
    (args.out / "graph_summary.csv").write_text("nodes,edges,components,mean_degree\n" + f"{len(g)},{len(edges)},{len(comps)},{2*len(edges)/len(g) if g else math.nan}\n")
    print(f"nodes={len(g)} edges={len(edges)} components={len(comps)}")

if __name__ == "__main__": main()
