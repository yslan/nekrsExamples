"""
ringwalk.py -- plan H, prototyped: close the multi-root holes WITHOUT subdividing.

THE PROBLEM.  A fibre can carry two roots inside one lattice interval.  Both
bracketing nodes then have the same sign, so no sign test sees the second
crossing.  The cell triangulates from an incomplete cut-point set and the fan
bridges across the gap, which does not leave a clean hole -- it adds spurious
area (measured: +9.3% on demo_3d/turbulent).

WHY NOT JUST FEED BOTH POINTS TO THE FAN.  The existing fan orders all cut
points by angle about their centroid into a SINGLE ring.  That is correct only
under the invariant `at most one cut point per edge`: two crossings on one edge
mean the surface enters and leaves through it, i.e. TWO sheets locally, and a
single ring stitches them into a bowtie.  v1 keeps the invariant by inserting a
plane through the whole element; this module keeps the cells and fixes the
connection step instead.

THE RULE THAT MAKES IT SOUND.  A face is shared by exactly two cells, and both
see that face's own four edges.  So every predicate is computed from data both
sides share:

    a face's cut points = the crossings on its 4 edges (0, 1 or 2 each)

Neither cell can disagree, because neither consults anything the other cannot
see.  (The tempting alternative -- "refine a CELL if any of its 12 edges is
split" -- is unsound: cell A can have a split edge that is not on the shared
face, so A subdivides the face and B does not.  T-junction, crack.)

THE WALK.  The surface enters a cell through its faces.  On a face, cut points
pair up into segments; chaining those segments face to face closes one or more
RINGS, and each ring is one sheet.  Fan each ring separately -- the existing
fan, applied per ring instead of once.

    ring walk with 1 ring  == today's output, exactly
    ring walk with k rings == k sheets, no bowtie
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "demo_2d"))

# Cube edges, IMPORTED from surface3d rather than restated.  Restating them is
# exactly the mistake that produced the first version of this file: a private
# edge list that looked plausible, disagreed with the one the edge table was
# built against, and left faces carrying a single cut point -- which is
# impossible for a closed sheet and silently dropped 200 of 254 configurations.
import surface3d as _S3

EDGES = list(_S3._EDGES)  # (lo corner, hi corner, direction)

# the 6 faces of a cube, each as its 4 corners in cyclic order, and its normal
# axis + side.  Corner numbering is c = bi + 2*bj + 4*bk.
FACES = [
    (0, (0, 2, 6, 4)),   # i = 0
    (0, (1, 5, 7, 3)),   # i = 1
    (1, (0, 4, 5, 1)),   # j = 0
    (1, (2, 3, 7, 6)),   # j = 1
    (2, (0, 1, 3, 2)),   # k = 0
    (2, (4, 6, 7, 5)),   # k = 1
]


def face_edges(fc):
    """The 4 edge indices of face `fc`, in cyclic order around it."""
    _, corners = FACES[fc]
    out = []
    for m in range(4):
        a, b = corners[m], corners[(m + 1) % 4]
        for e, (ea, eb, _) in enumerate(EDGES):
            if {ea, eb} == {a, b}:
                out.append(e)
                break
    assert len(out) == 4, fc
    return out


FACE_EDGES = [face_edges(f) for f in range(6)]

# which faces does edge e lie on?  (exactly 2)
EDGE_FACES = [[f for f in range(6) if e in FACE_EDGES[f]] for e in range(12)]


# Face ambiguity is decided by EVALUATING the interpolant at the face centre,
# not by Nielson-Hamann's bilinear saddle.  Measured on 56 ambiguous faces of
# `turbulent`, against a flood-fill of the real field:
#
#     bilinear saddle            44/56  (79%)
#     polynomial at that saddle  48/56  (86%)
#     polynomial at the CENTRE   55/56  (98%)
#
# See plans/I step 4.  `face_corner_vals` therefore carries the centre VALUE,
# not the four corners.


def rings_from_face_pairs(face_pts, sign_at, face_corner_vals=None):
    """Chain per-face segments into closed rings.

    face_pts[f] is the list of (edge, point_id) on face f, at most 2 per edge.
    On a face the cut points pair up into segments; with 2 points there is one
    segment, with 4 there are two and the pairing is decided by walking the
    face boundary in cyclic order and joining alternate crossings -- the same
    even-odd rule marching squares uses, which is unambiguous here because the
    points are ordered along the boundary.

    Returns a list of rings, each a list of point ids.
    """
    # Build segments per face.
    #
    # Walk the face boundary once, in cyclic order, collecting the crossings as
    # they occur.  Consecutive crossings along that walk bound an interval where
    # the field has one sign, so they alternate in/out and PAIRING CONSECUTIVE
    # ONES IS THE MARCHING-SQUARES RULE -- but only if the walk visits them in
    # true boundary order, INCLUDING two points on the same edge.
    #
    # That is the detail that broke the first version: `sorted by (edge, t)`
    # walks every edge low-to-high, whereas the cyclic boundary traverses half
    # of them backwards.  Two crossings on one traversed-backwards edge then
    # come out in the wrong order, the pairing shifts by one, and the rings do
    # not close -- 510 of 756 restored-pair cases failed that way.
    seg = []
    for f in range(6):
        pl = face_pts[f]
        if len(pl) < 2:
            continue
        _, corners = FACES[f]
        walk = []
        for m in range(4):
            a, b = corners[m], corners[(m + 1) % 4]
            e = FACE_EDGES[f][m]
            lo, hi, _ = EDGES[e]
            fwd = (lo == a)  # does the edge run the same way as the walk?
            onq = [(t, p) for (ee, p, t) in pl if ee == e]
            onq.sort(key=lambda tp: tp[0], reverse=not fwd)
            walk += [p for (_, p) in onq]
        # Four crossings is the ambiguous configuration: ask the field which
        # pairing, instead of always starting at index 0 (plan I).
        rot = 0
        if len(walk) == 4 and face_corner_vals is not None:
            cv, mid = face_corner_vals[f]
            if (mid < 0.0) != (cv[1] < 0.0):
                rot = 1
        n = len(walk)
        for m in range(0, n - 1, 2):
            seg.append((walk[(m + rot) % n], walk[(m + 1 + rot) % n]))

    # Chain the segments into rings.
    #
    # Every cut point lies on exactly 2 faces and so has exactly 2 segments:
    # the chain is a disjoint union of cycles, and following unused segments
    # walks each one.  Track used SEGMENTS, not visited points -- two points can
    # legitimately be joined on both of their shared faces, giving a 2-cycle (a
    # small sheet cut off by a restored pair), and a point-visited test breaks
    # the chain there instead of closing it.  That was the second bug: 510 of
    # 756 restored-pair cases reported "points lost" when the points were in
    # fact fine and the WALK was terminating early.
    adj = {}
    for i, (a, b) in enumerate(seg):
        adj.setdefault(a, []).append((b, i))
        adj.setdefault(b, []).append((a, i))
    used = [False] * len(seg)
    rings = []
    for s0 in sorted(adj):
        while True:
            start_seg = next((i for (_, i) in adj[s0] if not used[i]), None)
            if start_seg is None:
                break
            ring = [s0]
            cur, si = s0, start_seg
            while True:
                used[si] = True
                a, b = seg[si]
                cur = b if a == cur else a
                if cur == s0:
                    break                      # closed
                ring.append(cur)
                nxt = next((i for (_, i) in adj[cur] if not used[i]), None)
                if nxt is None:
                    break                      # open chain: should not happen
                si = nxt
            if len(ring) >= 2:
                rings.append(ring)
    return rings


def self_test(verbose=True):
    """Every invariant the walk must satisfy, over all 256 MC configurations.

    Run as `python3 ringwalk.py`.  Three of these were failing at some point
    during development and each failure was silent in the sense that the walk
    still returned rings -- just not covering every point.
    """
    import surface3d as S3

    ET = S3._EDGE_TABLE
    ED = S3._EDGES

    def cor(c):
        return np.array([c & 1, (c >> 1) & 1, (c >> 2) & 1], float)

    def canon(sq):
        sq = list(sq)
        i = sq.index(min(sq))
        s1 = sq[i:] + sq[:i]
        s2 = list(reversed(s1))
        s2 = s2[-1:] + s2[:-1]
        return min(tuple(s1), tuple(s2))

    nfail = 0

    def ck(ok, what, got):
        nonlocal nfail
        if verbose:
            print("  %-56s %8s  %s" % (what, got, "ok" if ok else "FAIL"))
        if not ok:
            nfail += 1

    # ---- 1. baseline: one cut point per crossing edge --------------------
    lost = odd = 0
    pl_same = pl_diff = 0
    for idx in range(1, 255):
        edges = [e for e in range(12) if ET[idx] >> e & 1]
        if len(edges) < 3:
            continue
        fp = [[] for _ in range(6)]
        P = []
        for n, e in enumerate(edges):
            a, b, _ = ED[e]
            P.append(0.5 * (cor(a) + cor(b)))
            for f in EDGE_FACES[e]:
                fp[f].append((e, n, 0.5))
        if any(len(fp[f]) % 2 for f in range(6)):
            odd += 1
        rings = rings_from_face_pairs(fp, None)
        if sum(len(r) for r in rings) != len(edges):
            lost += 1
            continue
        if len(rings) != 1:
            continue
        P = np.array(P)
        c0 = P.mean(0)
        _, _, vt = np.linalg.svd(P - c0)
        nrm = vt[2]
        if np.abs((P - c0) @ nrm).max() > 1e-9:
            continue  # non-planar: the angular fan is not a valid reference
        ang = np.arctan2((P - c0) @ np.cross(nrm, vt[0]), (P - c0) @ vt[0])
        if canon(list(np.argsort(ang))) == canon(rings[0]):
            pl_same += 1
        else:
            pl_diff += 1
    ck(odd == 0, "every face carries an EVEN number of cut points", odd)
    ck(lost == 0, "baseline: no cut point dropped", lost)
    ck(pl_diff == 0, "PLANAR cells: ring order == the angular fan (%d)" % pl_same, pl_diff)

    # ---- 2. a restored multi-root pair -----------------------------------
    # The dropped pair lives on a SAME-SIGNED edge, i.e. one the sign test says
    # does not cross, and restoring it adds TWO points there -- parity holds.
    lost2 = multi = 0
    for idx in range(1, 255):
        edges = [e for e in range(12) if ET[idx] >> e & 1]
        if len(edges) < 3:
            continue
        for e in [c for c in range(12) if not (ET[idx] >> c & 1)][:3]:
            fp = [[] for _ in range(6)]
            nid = 0
            for c in edges:
                for f in EDGE_FACES[c]:
                    fp[f].append((c, nid, 0.5))
                nid += 1
            for t in (0.3, 0.7):
                for f in EDGE_FACES[e]:
                    fp[f].append((e, nid, t))
                nid += 1
            rings = rings_from_face_pairs(fp, None)
            if sum(len(r) for r in rings) != nid:
                lost2 += 1
            elif len(rings) > 1:
                multi += 1
    ck(lost2 == 0, "restored pair: no cut point dropped (756 cases)", lost2)
    ck(multi > 0, "restored pair: extra sheets appear as extra rings", multi)

    if verbose:
        print("\n%s (%d failures)" % ("ALL PASS" if not nfail else "FAILED", nfail))
    return nfail


if __name__ == "__main__":
    raise SystemExit(self_test())


def extract_with_walk(N, F, keep_pairs=True):
    """Re-extract a surface using the RING WALK, so the fix can be measured
    against the analytic reference rather than only against itself.

    Mirrors surface3d.extract_surface but (a) keeps BOTH roots of a multi-root
    interval and (b) connects cut points by the face walk instead of the
    angular fan.  `keep_pairs=False` reproduces the old dropping behaviour, so
    the two can be compared on identical machinery.
    """
    import surface3d as S3

    r, _ = S3.gll_nodes_weights(N)
    ET = S3._EDGE_TABLE  # noqa: F841  (kept: the sign mask is no longer the driver)
    cross = {}
    for j in range(N + 1):
        for k in range(N + 1):
            for d in range(3):
                fib = F[:, j, k] if d == 0 else (F[j, :, k] if d == 1 else F[j, k, :])
                rts = S3.aberth_roots(N, fib)
                rts = rts[0] if isinstance(rts, tuple) else rts
                cross[(d, j, k)] = np.sort(np.asarray(rts, float))

    tris = []
    nunsep = 0
    for ci in range(N):
        for cj in range(N):
            for ck in range(N):
                # cut points on each of the 12 edges, up to two
                ePt = [[] for _ in range(12)]
                P = []
                for en, (a, b, d) in enumerate(EDGES):
                    ai = (a & 1, (a >> 1) & 1, (a >> 2) & 1)
                    if d == 0:
                        key, lo = (0, cj + ai[1], ck + ai[2]), ci
                    elif d == 1:
                        key, lo = (1, ci + ai[0], ck + ai[2]), cj
                    else:
                        key, lo = (2, ci + ai[0], cj + ai[1]), ck
                    rts = cross[key]
                    inr = rts[(rts > r[lo]) & (rts < r[lo + 1])]
                    if len(inr) > 1:
                        nunsep += len(inr) - 1
                        if not keep_pairs:
                            inr = inr[:1]
                    for t in inr[:2]:
                        if d == 0:
                            p = (float(t), r[cj + ai[1]], r[ck + ai[2]])
                        elif d == 1:
                            p = (r[ci + ai[0]], float(t), r[ck + ai[2]])
                        else:
                            p = (r[ci + ai[0]], r[cj + ai[1]], float(t))
                        ePt[en].append(len(P))
                        P.append(p)
                if len(P) < 3:
                    continue
                face_pts = [[] for _ in range(6)]
                for en in range(12):
                    d = EDGES[en][2]
                    for pid in ePt[en]:
                        for f in EDGE_FACES[en]:
                            face_pts[f].append((en, pid, P[pid][d]))
                fcv = []
                for f in range(6):
                    cs = FACES[f][1]
                    idx = [(ci + (c & 1), cj + ((c >> 1) & 1), ck + ((c >> 2) & 1)) for c in cs]
                    cv = [F[a, b, d] for (a, b, d) in idx]
                    rc = 0.5 * (r[idx[0][0]] + r[idx[2][0]])
                    sc = 0.5 * (r[idx[0][1]] + r[idx[2][1]])
                    tc = 0.5 * (r[idx[0][2]] + r[idx[2][2]])
                    fcv.append((cv, S3.tensor_eval3(F, rc, sc, tc)))
                for ring in rings_from_face_pairs(face_pts, None, fcv):
                    if len(ring) < 3:
                        continue
                    for q in range(1, len(ring) - 1):
                        tris.append((P[ring[0]], P[ring[q]], P[ring[q + 1]]))
    return tris, nunsep
