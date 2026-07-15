"""
generate_database.py

Run this locally to generate graphs_slim.json, then paste the data
into digraph_database.html to update the website.

Usage:
    python generate_database.py

Requirements: dichromatic.py in the same folder, plus sympy, matplotlib, numpy.
Edit D6_FILES to include whichever .d6 files you have downloaded.
"""

import json
import math
import sys
import functools
import sympy as sp
from sympy import nroots, Symbol, factor, roots as sympy_roots
from itertools import permutations
from multiprocessing import Pool, cpu_count

from sympy import symbols, Poly
from itertools import combinations
import matplotlib.pyplot as plt
import numpy as np

k = symbols('k')

D6_FILES = []

TOURNAMENT_FILES = [
    "tour9.d6",
]

OUTPUT_JSON = "graphs_tour_n9_169100-.json"
OUTPUT_HTML = "index.html"

CHECKPOINT_EVERY = 100   # save progress every N graphs

# To skip the first N lines in a specific file, set its value here.
# Useful to jump ahead when you know a portion is already processed.
# 0 = start from the beginning (default).
# The checkpoint system handles resumption automatically in most cases —
# use this as a manual override when you know exactly where to start.
START_FROM_LINE = {
    # "tour2.d6": 0,
    # "tour3.d6": 0,
    # "tour4.d6": 0,
    # "tour5.d6": 0,
    # "tour6.d6": 0,
    # "tour7.d6": 0,
    # "tour8.d6": 0,
    "tour9.d6": 169100,
}

# ============================================================
# digraph6 parsing / encoding
# ============================================================

def _n_from_bytes(data, pos):
    if data[pos] != 126:
        return data[pos] - 63, pos + 1
    pos += 1
    if data[pos] != 126:
        x = 0
        for _ in range(3):
            x = (x << 6) | (data[pos] - 63); pos += 1
        return x, pos
    pos += 1
    x = 0
    for _ in range(6):
        x = (x << 6) | (data[pos] - 63); pos += 1
    return x, pos

def parse_d6(s):
    """digraph6 string → n×n adjacency matrix."""
    s = s.strip()
    if s.startswith(">>digraph6<<"): s = s[len(">>digraph6<<"):]
    if s.startswith("&"): s = s[1:]
    data = s.encode("ascii")
    n, pos = _n_from_bytes(data, 0)
    bits = []
    for byte in data[pos:]:
        v = byte - 63
        for shift in range(5, -1, -1):
            bits.append((v >> shift) & 1)
    adj = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(n):
            if idx < len(bits):
                adj[i][j] = bits[idx]; idx += 1
    return adj

def to_d6(adj):
    """n×n adjacency matrix → digraph6 string."""
    n = len(adj)
    if n <= 62:
        n_bytes = bytes([n + 63])
    elif n <= 258047:
        x = n; b = []
        for _ in range(3): b.append((x & 0x3F) + 63); x >>= 6
        n_bytes = bytes([126]) + bytes(reversed(b))
    else:
        x = n; b = []
        for _ in range(6): b.append((x & 0x3F) + 63); x >>= 6
        n_bytes = bytes([126, 126]) + bytes(reversed(b))
    bits = [adj[i][j] for i in range(n) for j in range(n)]
    while len(bits) % 6: bits.append(0)
    r_bytes = []
    for i in range(0, len(bits), 6):
        v = 0
        for b in bits[i:i+6]: v = (v << 1) | b
        r_bytes.append(v + 63)
    return "&" + (n_bytes + bytes(r_bytes)).decode("ascii")

# ============================================================
# Graph utilities
# ============================================================

def _n(adj): return len(adj)
def _adj_copy(adj): return [row[:] for row in adj]

def _is_acyclic(adj):
    n = _n(adj)
    WHITE, GRAY, BLACK = 0, 1, 2
    color = [WHITE] * n
    def dfs(u):
        color[u] = GRAY
        for v in range(n):
            if adj[u][v]:
                if color[v] == GRAY: return True
                if color[v] == WHITE and dfs(v): return True
        color[u] = BLACK
        return False
    return not any(dfs(u) for u in range(n) if color[u] == 0)

def _strong_components(adj):
    n = _n(adj)
    visited = [False]*n; order = []
    def dfs1(u):
        stack = [(u, 0)]
        while stack:
            v, idx = stack[-1]
            if idx == 0: visited[v] = True
            found = False
            for w in range(idx, n):
                stack[-1] = (v, w+1)
                if adj[v][w] and not visited[w]:
                    visited[w] = True; stack.append((w, 0)); found = True; break
            if not found: stack.pop(); order.append(v)
    for u in range(n):
        if not visited[u]: dfs1(u)
    T = [[adj[j][i] for j in range(n)] for i in range(n)]
    visited2 = [False]*n; comps = []
    def dfs2(u):
        comp = []; stack = [u]; visited2[u] = True
        while stack:
            v = stack.pop(); comp.append(v)
            for w in range(n):
                if T[v][w] and not visited2[w]:
                    visited2[w] = True; stack.append(w)
        return comp
    for u in reversed(order):
        if not visited2[u]: comps.append(frozenset(dfs2(u)))
    return comps

def _add_arc(adj, u, v):
    a = _adj_copy(adj); a[u][v] = 1; return a

def _contract(adj, vertex_set):
    n = _n(adj); S = set(vertex_set)
    others = sorted(set(range(n)) - S)
    new_order = [min(S)] + others
    m = 1 + len(others)
    new_adj = [[0]*m for _ in range(m)]
    for i, u in enumerate(new_order):
        for j, v in enumerate(new_order):
            if i == j: continue
            u_group = S if i == 0 else {u}
            v_group = S if j == 0 else {v}
            new_adj[i][j] = 1 if any(adj[a][b] for a in u_group for b in v_group) else 0
    return new_adj

def _find_non_symmetric_pair(adj):
    n = _n(adj)
    for u in range(n):
        for v in range(u+1, n):
            if not (adj[u][v] and adj[v][u]): return u, v
    return None

def _paths_of_length_ge2(adj, u, v):
    n = _n(adj); results = []
    def dfs(current, path, visited):
        for nxt in range(n):
            if not adj[current][nxt]: continue
            if nxt == v and len(path) >= 2:
                full = path + [v]
                induced = [[adj[a][b] for b in full] for a in full]
                if _is_acyclic(induced): results.append(tuple(full))
            elif nxt != v and nxt not in visited:
                visited.add(nxt); dfs(nxt, path+[nxt], visited); visited.remove(nxt)
    dfs(u, [u], {u})
    return results

def _build_equivalence_classes(paths, u, v):
    if not paths: return []
    m = len(paths)
    seen = {}
    for r in range(1, m+1):
        for combo in combinations(range(m), r):
            vset = frozenset(vtx for idx in combo for vtx in paths[idx])
            if vset not in seen: seen[vset] = combo
    return list(seen.keys())

def _poly_complete_symmetric(n):
    result = Poly(1, k, domain='ZZ')
    for i in range(n): result = result * Poly(k - i, k, domain='ZZ')
    return result

def _rec(adj):
    n = _n(adj)
    if n <= 1:                    return Poly(k**n, k, domain='ZZ')
    if _is_acyclic(adj):          return Poly(k**n, k, domain='ZZ')
    comps = _strong_components(adj)
    if len(comps) > 1:
        result = Poly(1, k, domain='ZZ')
        for comp in comps:
            comp = sorted(comp)
            sub = [[adj[u][v] for v in comp] for u in comp]
            result = result * _rec(sub)
        return result
    pair = _find_non_symmetric_pair(adj)
    if pair is None: return _poly_complete_symmetric(n)
    u, v = pair
    p1 = _rec(_add_arc(_add_arc(adj, u, v), v, u))
    p2 = _rec(_contract(adj, {u, v}))
    all_paths = _paths_of_length_ge2(adj, u, v) + _paths_of_length_ge2(adj, v, u)
    R_stars   = _build_equivalence_classes(all_paths, u, v)
    p3 = Poly(0, k, domain='ZZ')
    for R_star in R_stars:
        R_list  = sorted(R_star)
        induced = [[adj[a][b] for b in R_list] for a in R_list]
        if not _is_acyclic(induced):
            continue
        p3 = p3 + _rec(_contract(adj, R_star))
    return p1 + p2 + p3

@functools.lru_cache(maxsize=None)
def _rec_cached(adj_tuple):
    adj = [list(row) for row in adj_tuple]
    return _rec(adj)

def dichromatic_poly(adj):
    adj_tuple = tuple(tuple(row) for row in adj)
    return _rec_cached(adj_tuple)

def dichromatic_number(adj):
    if _is_acyclic(adj): return 1
    poly = dichromatic_poly(adj)
    for kv in range(1, _n(adj) + 2):
        if poly.eval(kv) > 0: return kv
    return _n(adj)

def process_d6_string(s):
    s = s.strip()
    if not s or s.startswith("#"): return
    adj  = parse_d6(s)
    poly = dichromatic_poly(adj)
    dc   = dichromatic_number(adj)
    print(f"d6     : {s}")
    print(f"n      : {_n(adj)}")
    print(f"P(D,k) : {poly.as_expr()}")
    print(f"dc(D)  : {dc}")
    print()

def process_d6_file(path):
    with open(path) as fh:
        for line in fh:
            process_d6_string(line)

def _n_from_tri_len(L):
    n = int((1 + math.sqrt(1 + 8*L)) / 2)
    assert n*(n-1)//2 == L, f"String length {L} is not a valid upper-triangle length"
    return n

def parse_tournament(s):
    s = s.strip()
    n = _n_from_tri_len(len(s))
    adj = [[0]*n for _ in range(n)]
    idx = 0

    for i in range(n):
        for j in range(i+1, n):
            if s[idx] == '1':
                adj[i][j] = 1
            else:
                adj[j][i] = 1
            idx += 1

    return adj

def process_tournament_string(s):
    s = s.strip()
    if not s or s.startswith("#"): return
    adj  = parse_tournament(s)
    poly = dichromatic_poly(adj)
    dc   = dichromatic_number(adj)
    print(f"tournament : {s}")
    print(f"  n      = {_n(adj)}")
    print(f"  P(D,k) = {poly.as_expr()}")
    print(f"  dc(D)  = {dc}")
    print()

def process_tournament_file(path):
    with open(path) as fh:
        for line in fh:
            process_tournament_string(line)

def print_graph(adj, labels=None):
    n = len(adj)
    if labels is None:
        labels = [str(i) for i in range(n)]
    print("Digraph:")
    for i in range(n):
        arcs = [labels[j] for j in range(n) if adj[i][j]]
        if arcs:
            print(f"  {labels[i]} → {', '.join(arcs)}")
        else:
            print(f"  {labels[i]} → (none)")

def draw_graph(adj, labels=None, title="Digraph"):
    n = len(adj)
    if labels is None:
        labels = [str(i) for i in range(n)]
    angles = [2 * np.pi * i / n for i in range(n)]
    pos = {i: (np.cos(a), np.sin(a)) for i, a in enumerate(angles)}
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title(title, fontsize=13)
    for i in range(n):
        for j in range(n):
            if not adj[i][j]: continue
            xi, yi = pos[i]
            xj, yj = pos[j]
            if adj[j][i]:
                dx, dy = yj - yi, -(xj - xi)
                norm = np.sqrt(dx**2 + dy**2) + 1e-9
                off = 0.08
                xi2 = xi + off*dx/norm; yi2 = yi + off*dy/norm
                xj2 = xj + off*dx/norm; yj2 = yj + off*dy/norm
            else:
                xi2, yi2, xj2, yj2 = xi, yi, xj, yj
            ax.annotate("",
                xy=(xj2, yj2), xytext=(xi2, yi2),
                arrowprops=dict(arrowstyle="-|>", color="steelblue",
                                lw=1.5, mutation_scale=15,
                                shrinkA=14, shrinkB=14))
    for i in range(n):
        x, y = pos[i]
        ax.add_patch(plt.Circle((x, y), 0.1, color='steelblue', zorder=3))
        ax.text(x, y, labels[i], ha='center', va='center',
                color='white', fontsize=11, fontweight='bold', zorder=4)
    ax.set_xlim(-1.4, 1.4); ax.set_ylim(-1.4, 1.4)
    plt.tight_layout()
    plt.show()

# ============================================================
# Tournament-specific matrices (skew, walk, hermitian)
# ============================================================

def skew_matrix(adj):
    """
    Skew-adjacency matrix S of a tournament.
      S[i][j] =  1  if i -> j
                -1  if j -> i
                 0  if i == j
    Skew-symmetric by construction (S^T = -S).
    """
    n = len(adj)
    S = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            if adj[i][j]:
                S[i][j] = 1
            elif adj[j][i]:
                S[i][j] = -1
    return S

def adjacency_walk_matrix(adj):
    """
    Adjacency walk matrix W_A = [e | Ae | A^2 e | ... | A^{n-1} e]
    where e is the all-ones vector and A is the adjacency matrix.
    det(W_A) encodes the number of edges/arcs via its spectrum.
    """
    n = len(adj)
    A = np.array(adj, dtype=float)
    e = np.ones(n, dtype=float)
    W = np.zeros((n, n), dtype=float)
    Ak_e = e.copy()
    for k in range(n):
        W[:, k] = Ak_e
        Ak_e = A @ Ak_e
    return W

def skew_walk_matrix(adj):
    """
    Skew walk matrix W_S = [e | Se | S^2 e | ... | S^{n-1} e]
    where S is the skew-adjacency matrix.
    Used in DGSS (determined by generalised skew-adjacency spectrum) theory.
    """
    n = len(adj)
    S = skew_matrix(adj).astype(float)
    e = np.ones(n, dtype=float)
    W = np.zeros((n, n), dtype=float)
    Sk_e = e.copy()
    for k in range(n):
        W[:, k] = Sk_e
        Sk_e = S @ Sk_e
    return W

def hermitian_adj(adj):
    """
    Hermitian adjacency matrix H.
      H[i][j] =  1   if i ~ j  (undirected edge — absent in tournaments)
                 i   if i -> j
                -i   if j -> i
                 0   otherwise
    H is Hermitian (H^* = H) so its eigenvalues are real.
    """
    n = len(adj)
    H = np.zeros((n, n), dtype=complex)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            if adj[i][j] and adj[j][i]:
                H[i][j] = 1.0
            elif adj[i][j]:
                H[i][j] = 1j
            elif adj[j][i]:
                H[i][j] = -1j
    return H

def hermitian_walk_matrix(adj):
    """
    Hermitian walk matrix W_H = [e | He | H^2 e | ... | H^{n-1} e]
    where H is the Hermitian adjacency matrix.
    """
    n = len(adj)
    H = hermitian_adj(adj).astype(complex)
    e = np.ones(n, dtype=complex)
    W = np.zeros((n, n), dtype=complex)
    Hk_e = e.copy()
    for k in range(n):
        W[:, k] = Hk_e
        Hk_e = H @ Hk_e
    return W


# ============================================================
# Spectrum / determinant helpers
# ============================================================

def _clean_val(val, tol=1e-9):
    """Snap near-integer / near-real complex values for clean display."""
    val = complex(val)
    re = int(round(val.real)) if abs(val.real - round(val.real)) < tol else round(val.real, 10)
    im = int(round(val.imag)) if abs(val.imag - round(val.imag)) < tol else round(val.imag, 10)
    if abs(im) < tol:
        return re
    return complex(re, im)

def spectrum_of(M, tol=1e-9):
    """
    Return sorted eigenvalues of M.
    Uses eigvalsh for Hermitian/symmetric matrices (guaranteed real output),
    eigvals otherwise.
    """
    if np.allclose(M, M.conj().T, atol=tol):
        eigs = np.linalg.eigvalsh(M)
        return sorted([_clean_val(e, tol) for e in eigs])
    else:
        eigs = np.linalg.eigvals(M)
        return sorted(
            [_clean_val(e, tol) for e in eigs],
            key=lambda x: (x.real if isinstance(x, complex) else x,
                           x.imag if isinstance(x, complex) else 0)
        )

def det_of(M, tol=1e-9):
    """
    Determinant of M, snapped to integer when close.
    Singular matrices (det=0) are handled cleanly without NumPy warnings.
    """
    M = M.astype(complex)
    # Check for singularity first via rank to avoid NumPy divide-by-zero warnings
    rank = np.linalg.matrix_rank(M)
    if rank < M.shape[0]:
        return 0
    with np.errstate(divide='ignore', invalid='ignore'):
        d = np.linalg.det(M)
    return _clean_val(d, tol)

def tournament_matrices(adj):
    """
    Compute all four tournament matrices, their determinants and spectra.

    Returns a dict with keys:
        S, WA, WS, WH          — the matrices (np.ndarray)
        det_S, det_WA, det_WS, det_WH  — determinants
        spec_S, spec_WA, spec_WS, spec_WH  — eigenvalue lists
    """
    S  = skew_matrix(adj)
    WA = adjacency_walk_matrix(adj)
    WS = skew_walk_matrix(adj)
    WH = hermitian_walk_matrix(adj)

    return {
        "S":       S,
        "WA":      WA,
        "WS":      WS,
        "WH":      WH,
        "det_S":   det_of(S.astype(float)),
        "det_WA":  det_of(WA),
        "det_WS":  det_of(WS),
        "det_WH":  det_of(WH),
        "spec_S":  spectrum_of(S.astype(float)),
        "spec_WA": spectrum_of(WA),
        "spec_WS": spectrum_of(WS),
        "spec_WH": spectrum_of(WH),
    }


# ============================================================
# Existing helpers (unchanged)
# ============================================================

def is_tournament(adj):
    n = len(adj)
    for i in range(n):
        for j in range(i+1, n):
            if adj[i][j] + adj[j][i] != 1:
                return False
    return True

def to_tournament_str(adj):
    n = len(adj)
    return "".join(str(adj[i][j]) for i in range(n) for j in range(i+1, n))

def count_edges(adj):
    n = len(adj)
    return sum(adj[i][j] for i in range(n) for j in range(n))

def clean_number(x, tol=1e-10):
    x = complex(x)
    if abs(x.imag) < tol:
        real = x.real
        if abs(real - round(real)) < tol:
            return int(round(real))
        return round(real, 12)
    real = x.real
    imag = x.imag
    if abs(real - round(real)) < tol:
        real = int(round(real))
    else:
        real = round(real, 12)
    if abs(imag - round(imag)) < tol:
        imag = int(round(imag))
    else:
        imag = round(imag, 12)
    return complex(real, imag)

def format_roots(poly_expr, var=None, prec=50, tol=1e-10):
    if var is None:
        var = list(poly_expr.free_symbols)[0]
    poly = sp.Poly(poly_expr, var)
    roots_out = []
    try:
        factors = sp.factor_list(poly_expr)[1]
        for factor, mult in factors:
            rdict = sp.roots(factor, var)
            for root, m in rdict.items():
                val = sp.N(root)
                roots_out.extend([val] * (m * mult))
        if not roots_out:
            return ["N/A"]
        return [clean_number(r, tol) for r in roots_out]
    except Exception:
        pass
    roots = sp.nroots(poly_expr, n=prec)
    if not roots:
        return ["N/A"]
    return [clean_number(r, tol) for r in roots]

def hermitian_eigenvalues(adj, tol=1e-8):
    H = hermitian_adj(adj)
    eigs = np.linalg.eigvalsh(H)
    result = []
    for e in sorted(eigs):
        e_real = float(e.real)
        if abs(e_real - round(e_real)) < tol:
            result.append(int(round(e_real)))
        else:
            result.append(round(e_real, 8))
    return result

def _fmt_eig(e, tol=1e-8):
    re, im = float(e.real), float(e.imag)
    if abs(im) < tol:
        if abs(re - round(re)) < tol:
            return str(int(round(re)))
        return f"{re:.8f}"
    mod = np.sqrt(re**2 + im**2)
    theta = np.arctan2(im, re)
    return f"{round(mod,8)}e^(i*{round(theta,8)})"

def adjacency_eigenvalues(adj, tol=1e-8):
    A = np.array(adj, dtype=float)
    eigs = np.linalg.eigvals(A)
    result = []
    for e in sorted(eigs, key=lambda x: (x.real, x.imag)):
        result.append(_fmt_eig(e, tol))
    return result

def format_eigenvalues(eig_list):
    return ", ".join(str(e) for e in eig_list)

def _transpose(adj):
    n = len(adj)
    return [[adj[j][i] for j in range(n)] for i in range(n)]

def is_self_converse(adj):
    n = len(adj)
    conv = _transpose(adj)
    if adj == conv:
        return True
    if n > 8:
        return None
    for perm in permutations(range(n)):
        if all(conv[perm[i]][perm[j]] == adj[i][j]
               for i in range(n) for j in range(n)):
            return True
    return False

# ============================================================
# Graph processing (extended with tournament matrices)
# ============================================================

def format_root_val(r):
    if isinstance(r, complex):
        re_s = f"{r.real:.8f}" if abs(r.real - round(r.real)) > 1e-8 else str(int(round(r.real)))
        im_s = f"{abs(r.imag):.8f}" if abs(r.imag - round(r.imag)) > 1e-8 else str(int(round(abs(r.imag))))
        sign = "+" if r.imag >= 0 else "-"
        return f"{re_s}{sign}{im_s}i"
    elif isinstance(r, float):
        return f"{r:.8f}" if abs(r - round(r)) > 1e-8 else str(int(round(r)))
    else:
        return str(r)

def _fmt_complex_val(v):
    """Convert a single complex/float/int value to a JSON-safe string."""
    if isinstance(v, complex) or (hasattr(v, 'imag') and v.imag != 0):
        re = v.real; im = v.imag
        re_s = str(int(round(re))) if abs(re - round(re)) < 1e-8 else f"{re:.8f}"
        im_s = str(int(round(abs(im)))) if abs(abs(im) - round(abs(im))) < 1e-8 else f"{abs(im):.8f}"
        sign = "+" if im >= 0 else "-"
        return f"{re_s}{sign}{im_s}i"
    else:
        v = v.real if hasattr(v, 'real') else v
        return str(int(round(v))) if abs(v - round(v)) < 1e-8 else f"{v:.8f}"

def _matrix_to_json(M):
    """
    Convert a matrix (possibly complex) to a JSON-serializable list of lists.
    Real matrices stay as numbers; complex entries are serialized as strings.
    """
    if np.iscomplexobj(M) and not np.allclose(M.imag, 0):
        return [[_fmt_complex_val(v) for v in row] for row in M]
    else:
        M = M.real
        return [[int(round(v)) if abs(v - round(v)) < 1e-8 else round(float(v), 10)
                 for v in row] for row in M]

def _fmt_spec(spec):
    """Format a spectrum list (may contain int, float, or complex) as a string."""
    parts = []
    for v in spec:
        if isinstance(v, complex):
            re_s = str(int(round(v.real))) if abs(v.real - round(v.real)) < 1e-8 else f"{v.real:.8f}"
            im_s = str(int(round(abs(v.imag)))) if abs(v.imag - round(v.imag)) < 1e-8 else f"{abs(v.imag):.8f}"
            sign = "+" if v.imag >= 0 else "-"
            parts.append(f"{re_s}{sign}{im_s}i")
        elif isinstance(v, float):
            parts.append(str(int(round(v))) if abs(v - round(v)) < 1e-8 else f"{v:.8f}")
        else:
            parts.append(str(v))
    return ", ".join(parts)

def process_graph(args):
    """Process a single graph. args = (adj, upper_tri_string_or_None)"""
    adj, upper_tri = args
    n = _n(adj)
    e = count_edges(adj)
    poly = dichromatic_poly(adj)
    poly_expr = poly.as_expr()
    dc = dichromatic_number(adj)
    tourn = is_tournament(adj)
    acyc = _is_acyclic(adj)

    roots_list = format_roots(poly_expr)
    roots_str = ", ".join(format_root_val(r) for r in roots_list)

    h_eigs = format_eigenvalues(hermitian_eigenvalues(adj))
    a_eigs = format_eigenvalues(adjacency_eigenvalues(adj))

    h_eigs_list = hermitian_eigenvalues(adj)
    a_eigs_raw = np.linalg.eigvals(np.array(adj, dtype=float))

    h_spectral_radius = max(abs(float(v)) for v in h_eigs_list) if h_eigs_list else 0
    a_spectral_radius = float(max(abs(v) for v in a_eigs_raw)) if len(a_eigs_raw) > 0 else 0

    h_spectral_radius = int(round(h_spectral_radius)) if abs(h_spectral_radius - round(h_spectral_radius)) < 1e-8 else round(h_spectral_radius, 8)
    a_spectral_radius = int(round(a_spectral_radius)) if abs(a_spectral_radius - round(a_spectral_radius)) < 1e-8 else round(a_spectral_radius, 8)

    comps = _strong_components(adj)
    n_components = len(comps)
    is_strongly_connected = n_components == 1

    self_conv = is_self_converse(adj)

    # --- tournament matrices (only computed for tournaments) ---
    tm = {}
    if tourn:
        tm = tournament_matrices(adj)

    return {
        "d6":               str(to_d6(adj)),
        "n":                int(n),
        "e":                int(e),
        "poly":             str(poly_expr),
        "dc":               int(dc),
        "roots":            roots_str,
        "tournament":       bool(tourn),
        "acyclic":          bool(acyc),
        "adj":              [[int(x) for x in row] for row in adj],
        "upper_tri":        upper_tri,
        "h_eigs":           h_eigs,
        "a_eigs":           a_eigs,
        "h_spectral_radius": h_spectral_radius,
        "a_spectral_radius": a_spectral_radius,
        "n_scc":            int(n_components),
        "strongly_connected": bool(is_strongly_connected),
        "self_converse":    self_conv,
        # tournament-specific fields (empty strings if not a tournament)
        "skew_matrix":      _matrix_to_json(tm["S"])  if tm else [],
        "det_S":            str(tm["det_S"])           if tm else "",
        "spec_S":           _fmt_spec(tm["spec_S"])    if tm else "",
        "adj_walk_matrix":  _matrix_to_json(tm["WA"]) if tm else [],
        "det_WA":           str(tm["det_WA"])          if tm else "",
        "spec_WA":          _fmt_spec(tm["spec_WA"])   if tm else "",
        "skew_walk_matrix": _matrix_to_json(tm["WS"]) if tm else [],
        "det_WS":           str(tm["det_WS"])          if tm else "",
        "spec_WS":          _fmt_spec(tm["spec_WS"])   if tm else "",
        "herm_walk_matrix": _matrix_to_json(tm["WH"]) if tm else [],
        "det_WH":           str(tm["det_WH"])          if tm else "",
        "spec_WH":          _fmt_spec(tm["spec_WH"])   if tm else "",
    }

# ============================================================
# File / generation config
# ============================================================

def _save_checkpoint(graphs, path):
    """Write current results to disk atomically (write to .tmp then rename)."""
    tmp = path + ".tmp"
    with open(tmp, 'w') as f:
        json.dump(graphs, f, separators=(',', ':'))
    import os; os.replace(tmp, path)
    print(f"    [checkpoint] {len(graphs)} graphs → {path}", flush=True)

def _load_checkpoint(path):
    """Load existing checkpoint if present, return (graphs, processed_d6_set)."""
    import os
    if not os.path.exists(path):
        return [], set()
    try:
        with open(path) as f:
            graphs = json.load(f)
        done = {g['d6'] for g in graphs}
        print(f"  Resuming from checkpoint: {len(graphs)} graphs already done.")
        return graphs, done
    except Exception as e:
        print(f"  WARNING: could not load checkpoint ({e}), starting fresh.")
        return [], set()

def generate():
    total_files = 0

    # resume from checkpoint if one exists
    graphs, done_d6s = _load_checkpoint(OUTPUT_JSON)

    for fname in D6_FILES:
        try:
            with open(fname) as f:
                lines = [l.strip() for l in f if l.strip()]
            print(f"  {fname}: {len(lines)} graphs", flush=True)
            total_files += 1
            for i, line in enumerate(lines):
                if i % 100 == 0:
                    print(f"    {i}/{len(lines)}...", end='\r', flush=True)
                adj = parse_d6(line)
                d6  = to_d6(adj)
                if d6 in done_d6s:
                    continue
                upper_tri = str(to_tournament_str(adj)) if is_tournament(adj) else None
                graphs.append(process_graph((adj, upper_tri)))
                done_d6s.add(d6)
                if len(graphs) % CHECKPOINT_EVERY == 0:
                    _save_checkpoint(graphs, OUTPUT_JSON)
            print(f"    done.           ", flush=True)
        except FileNotFoundError:
            print(f"  WARNING: {fname} not found, skipping.")

    for fname in TOURNAMENT_FILES:
        try:
            with open(fname) as f:
                lines = [l.strip() for l in f if l.strip()]
            print(f"  {fname}: {len(lines)} graphs", flush=True)
            total_files += 1

            # skip first N lines if configured
            skip = START_FROM_LINE.get(fname, 0)
            if skip > 0:
                print(f"    Skipping first {skip} graphs as configured.", flush=True)
                lines = lines[skip:]

            # filter out already-processed lines
            pending = [(parse_tournament(line), line) for line in lines
                       if str(to_d6(parse_tournament(line))) not in done_d6s]
            print(f"    {len(pending)} remaining (skipping {len(lines)-len(pending)} done)",
                  flush=True)

            if not pending:
                print(f"    all done, skipping.", flush=True)
                continue

            ncpus = cpu_count()
            print(f"    Processing with {ncpus} CPUs...", flush=True)
            with Pool(processes=ncpus) as pool:
                batch = []
                for i, result in enumerate(pool.imap(process_graph, pending,
                                                      chunksize=50)):
                    batch.append(result)
                    done_d6s.add(result['d6'])
                    if i % 100 == 0:
                        print(f"    {i}/{len(pending)}...", end='\r', flush=True)
                    if len(batch) >= CHECKPOINT_EVERY:
                        graphs.extend(batch)
                        batch = []
                        _save_checkpoint(graphs, OUTPUT_JSON)
                graphs.extend(batch)
            print(f"    done.           ", flush=True)
            _save_checkpoint(graphs, OUTPUT_JSON)
        except FileNotFoundError:
            print(f"  WARNING: {fname} not found, skipping.")

    print(f"\nTotal: {len(graphs)} graphs from {total_files} files.")

    with open(OUTPUT_JSON, 'w') as f:
        json.dump(graphs, f, separators=(',', ':'))
    print(f"Wrote {OUTPUT_JSON} ({len(json.dumps(graphs))//1024}KB)")

    if OUTPUT_HTML:
        update_html(graphs)

def update_html(graphs):
    try:
        with open(OUTPUT_HTML) as f:
            html = f.read()
    except FileNotFoundError:
        print(f"WARNING: {OUTPUT_HTML} not found, skipping HTML update.")
        return

    DATA = json.dumps(graphs, separators=(',', ':'))

    import re
    start = html.find('const RAW=[')
    if start == -1:
        start = html.find('const RAW = [')
    if start == -1:
        print("WARNING: Could not find 'const RAW = ...' in HTML to replace.")
        print(f"  Manually replace it with the contents of {OUTPUT_JSON}")
        return
    end = html.rfind('];', start) + 2
    new_html = html[:start] + 'const RAW=' + DATA + ';' + html[end:]

    with open(OUTPUT_HTML, 'w') as f:
        f.write(new_html)
    print(f"Updated {OUTPUT_HTML} ({len(new_html)//1024}KB)")

if __name__ == "__main__":
    generate()
