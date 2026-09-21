"""Per-gene anglemania scores (``signal`` / ``R``) and expression-binned selection.

Implements the proposal in ``docs/making_anglemania_per_gene_level_score.md``.
anglemania scores *gene pairs* and turns that into a gene list via
``extract_unique_genes``, so a gene's effective score is the rank of the single
best pair it participates in -- a max-statistic over ~p pairs. This module
replaces that with row-sum statistics over the same accumulators.

From the cross-batch weighted mean ``M`` and sd ``S`` of the per-batch z-score
matrices (``StreamingZscoreStats.finalize()``), and using
``sum_b w_b (z_ij^b)^2 == M_ij^2 + S_ij^2``:

    signal_i = sum_j M_ij^2          # reproducible relational energy
    noise_i  = sum_j S_ij^2          # batch-specific relational energy
    R_i      = signal_i / (signal_i + noise_i)

``R_i`` is a repeatability coefficient (an ICC): the fraction of gene i's total
relational energy that is shared across batches rather than batch-specific.
Aggregating before dividing (rather than averaging the per-pair SNR ``|M|/S``)
matters because ``S`` has only ``B - 1`` degrees of freedom, so a per-pair ``S``
that happens to come out small inflates that pair's SNR arbitrarily.

The two axes are deliberately kept separate: high signal / high R are the
integration features, high signal / low R are the batch-sensitive genes, low
signal / high R are stable but uninformative.

**Expression binning.** Highly expressed genes have less sampling noise in their
correlations, hence smaller ``S`` and larger raw correlations -- so both scores
are partly expression scores. DUBStepR bins genes by mean expression before
z-scoring its correlation-range statistic for this reason, and the same fix is
applied here: genes go into equal-frequency bins of mean log-normalized
expression, and ``log10(signal)`` and ``R`` are z-scored *within* each bin.

**Asymmetry.** ``factorise``'s z-score matrix is deliberately asymmetric (entry
(i, j) is standardized against gene j's null), so row sums and column sums of
``M`` are different objects. The headline scores use the symmetrized
``(M + M^T) / 2``; the raw row/column margins are computed too, for diagnostics.

Requires the ``pyanglemania`` conda env.
"""

from __future__ import annotations

import sys
import time
from collections import Counter

import numpy as np
import pandas as pd

from pyanglemania.preprocessing._angles import factorise, factorise_chunked
from pyanglemania.preprocessing._batches import (
    add_unique_batch_key,
    align_to_common_genes,
    compute_dataset_weights,
    genes_passing_min_cells,
    intersect_genes,
    split_obs_indices_by_batch,
)
from pyanglemania.preprocessing._stats import StreamingZscoreStats

NORMALIZATION_METHODS = {
    "classical": "divide_by_total_counts",
    "pflog1ppf": "pflog1ppf",
}


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True, file=sys.stderr)


# ---------------------------------------------------------------------------
# cross-batch M / S
# ---------------------------------------------------------------------------
def compute_MS(
    adata,
    *,
    batch_key,
    dataset_key=None,
    layer=None,
    method="cosine",
    normalization_mode="classical",
    permutation_function="sample",
    min_cells_per_gene=1,
    min_samples_per_gene=2,
    allow_missing_features=False,
    cell_chunk_size=0,
    seed=1,
    gpu=False,
):
    """``anglemania()``'s own pipeline, stopped before ``prefilter_gene_pairs``.

    Returns ``(mean_zscore, sds_zscore, common_genes)`` -- the two full
    ``(genes x genes)`` matrices that ``anglemania()`` consumes internally and
    never exposes, plus the gene names they are indexed by. Both are host
    (numpy) arrays: ``StreamingZscoreStats`` keeps its accumulators on the host
    regardless of backend, so ``gpu=True`` only moves the per-batch angle
    computation onto the device.

    Unlike ``anglemania()``, this does *not* free ``adata.X`` on the GPU path
    (it copies to device instead of replacing the layer), because callers here
    still need the counts afterwards for the expression covariate that the
    binning is computed on.
    """
    if gpu:
        import cupy as cp
        from cupyx.scipy import sparse as csp
        from scipy import sparse as sp

        xp, sp_mod = cp, csp
        X_src = adata.X if layer is None else adata.layers[layer]
        X_full = csp.csr_matrix(X_src) if sp.issparse(X_src) else cp.asarray(X_src)
        del X_src
    else:
        from scipy import sparse as sp

        xp, sp_mod = np, sp
        X_full = adata.X if layer is None else adata.layers[layer]

    add_unique_batch_key(adata, batch_key, dataset_key)
    weights = compute_dataset_weights(adata.obs, batch_key, dataset_key)
    batch_indices = split_obs_indices_by_batch(adata)
    var_names = np.asarray(adata.var_names)
    log(f"{len(batch_indices)} batches: "
        + ", ".join(f"{b}(n={len(i)})" for b, i in batch_indices.items()))

    batch_X, batch_genes = {}, {}
    for label, idx in batch_indices.items():
        X_b = X_full[idx]
        if sp_mod.issparse(X_b):
            X_b = X_b.tocsc()
        mask = genes_passing_min_cells(X_b, min_cells_per_gene, xp, sp_mod)
        batch_genes[label] = list(var_names[mask])
        batch_X[label] = X_b[:, mask]
    del X_full

    common_genes = intersect_genes(
        list(batch_genes.values()), allow_missing_features, min_samples_per_gene,
        verbose=True,
    )
    mode = (f">={min_samples_per_gene}-of-{len(batch_genes)}-batches"
            if allow_missing_features else "all-batch intersection")
    log(f"{len(common_genes)} genes pass the {mode} filter (of {adata.n_vars})")
    # How many batches each common gene is actually present in. With
    # allow_missing_features=True a gene can be absent from up to
    # len(batches) - min_samples_per_gene of them, and align_to_common_genes
    # ZERO-FILLS those batches rather than masking them -- the accumulator then
    # averages over every batch regardless. So a gene's S (and hence `noise`)
    # picks up spurious cross-batch spread purely from how often it was missing.
    # Returned so callers can test for that confound instead of assuming it away.
    present = Counter(g for genes in batch_genes.values() for g in genes)
    n_batches_present = np.array([present[g] for g in common_genes], dtype=np.int32)

    stats = StreamingZscoreStats(len(common_genes))
    for label in batch_indices:
        t0 = time.time()
        if cell_chunk_size > 0:
            # Memory-bounded path, required at the gene counts
            # allow_missing_features produces (~14k genes): the unchunked route
            # materializes several (cells x genes) dense arrays at once and OOMs a
            # 24 GB card. factorise_chunked aligns per chunk, so it takes X still
            # restricted to this batch's own genes.
            z = factorise_chunked(
                batch_X[label], batch_genes[label], common_genes, xp,
                cell_chunk_size=cell_chunk_size,
                method=method,
                seed=seed,
                permute_row_or_column="column",
                permutation_function=permutation_function,
                normalization_method=NORMALIZATION_METHODS[normalization_mode],
                do_normalize=True,
            )
        else:
            X_dense = align_to_common_genes(batch_X[label], batch_genes[label],
                                            common_genes, xp)
            z = factorise(
                X_dense, xp,
                method=method,
                seed=seed,
                permute_row_or_column="column",
                permutation_function=permutation_function,
                normalization_method=NORMALIZATION_METHODS[normalization_mode],
                do_normalize=True,
            )
            del X_dense
        stats.update(z, float(weights[label]))
        del z
        if gpu:
            xp.get_default_memory_pool().free_all_blocks()
        log(f"  batch {label} folded in ({time.time() - t0:.1f}s)")
    del batch_X

    mean_z, sd_z, _sn = stats.finalize()
    return mean_z, sd_z, np.asarray(common_genes), n_batches_present


# ---------------------------------------------------------------------------
# gene-level reduction
# ---------------------------------------------------------------------------
def _row_sumsq(A, chunk=2000, symmetrize=False):
    """``sum_j A_ij^2`` per row, without materializing a second (genes x genes) array.

    With ``symmetrize=True`` the sum is over ``((A_ij + A_ji) / 2)^2`` instead;
    the symmetrized block is built one row-chunk at a time, so this still costs
    ``O(chunk x genes)`` rather than another full matrix (~3 GB at 19k genes).
    """
    out = np.empty(A.shape[0], dtype=np.float64)
    for s in range(0, A.shape[0], chunk):
        e = s + chunk
        block = A[s:e]
        if symmetrize:
            block = 0.5 * (block + A[:, s:e].T)
        out[s:e] = np.einsum("ij,ij->i", block, block)
    return out


def gene_level_scores(mean_z, sd_z):
    """``signal``/``noise``/``R`` per gene, symmetrized plus each raw margin.

    The diagonal is excluded throughout: ``M``'s diagonal is a gene's angle
    with itself (cosine == 1, a large z-score carrying no relational
    information) and ``S``'s diagonal is NaN by construction in ``finalize()``.
    Both inputs are modified in place (diagonal zeroed).
    """
    n = mean_z.shape[0]
    idx = np.arange(n)
    mean_z[idx, idx] = 0.0
    sd_z[idx, idx] = 0.0  # was NaN

    out = {}
    out["signal_row"] = _row_sumsq(mean_z)
    out["noise_row"] = _row_sumsq(sd_z)
    out["signal_col"] = _row_sumsq(mean_z.T)
    out["noise_col"] = _row_sumsq(sd_z.T)
    out["signal"] = _row_sumsq(mean_z, symmetrize=True)
    out["noise"] = _row_sumsq(sd_z, symmetrize=True)

    for tag in ("", "_row", "_col"):
        s, nz = out["signal" + tag], out["noise" + tag]
        out["R" + tag] = s / (s + nz)
    return out


def expression_covariates(adata, gene_names, layer=None, chunk=2000):
    """Mean CP10K+log1p expression and detection rate per gene, over all cells.

    This is the covariate the binning conditions on, so it is deliberately the
    plain CP10K+log1p mean regardless of which ``normalization_mode``
    anglemania itself ran with -- it stands in for "how well measured is this
    gene", not for anything in the angle computation.
    """
    from scipy import sparse as sp

    X = adata.X if layer is None else adata.layers[layer]
    gene_pos = pd.Index(adata.var_names).get_indexer(gene_names)
    if (gene_pos < 0).any():
        raise ValueError("some scored genes are missing from adata.var_names")

    if sp.issparse(X):
        X = X.tocsr()
        total = np.asarray(X.sum(axis=1)).ravel().astype(np.float64)
    else:
        X = np.asarray(X)
        total = X.sum(axis=1, dtype=np.float64)
    total[total == 0] = 1.0

    mean_log = np.zeros(len(gene_pos))
    det = np.zeros(len(gene_pos))
    for s in range(0, len(gene_pos), chunk):
        cols = gene_pos[s:s + chunk]
        sub = X[:, cols]
        sub = sub.toarray() if sp.issparse(sub) else np.asarray(sub)
        det[s:s + chunk] = (sub > 0).mean(axis=0)
        mean_log[s:s + chunk] = np.log1p(sub / total[:, None] * 1e4).mean(axis=0)
    return mean_log, det


def add_binned_zscores(df, n_bins=20, weight_signal=0.5, weight_R=0.5):
    """Add within-expression-bin z-scores of ``log10(signal)`` and ``R``, and scores.

    Equal-frequency bins on ``mean_lognorm`` (ranked first, so that datasets
    with many tied near-zero expression values still split into ``n_bins``
    equal-sized bins rather than collapsing). Both the binned score and its
    unbinned counterpart are added, so the effect of binning stays measurable
    rather than assumed.
    """
    df = df.copy()
    df["expr_bin"] = pd.qcut(
        df["mean_lognorm"].rank(method="first"), n_bins, labels=False, duplicates="drop"
    )
    df["log_signal"] = np.log10(df["signal"] + 1e-300)

    def _z(s):
        sd = s.std(ddof=1)
        return (s - s.mean()) / sd if sd > 0 else s * 0.0

    for col, new in (("log_signal", "signal_z_binned"), ("R", "R_z_binned")):
        df[new] = df.groupby("expr_bin")[col].transform(_z)
        df[new.replace("_binned", "_global")] = _z(df[col])

    df["score_binned"] = weight_signal * df["signal_z_binned"] + weight_R * df["R_z_binned"]
    df["score_unbinned"] = weight_signal * df["signal_z_global"] + weight_R * df["R_z_global"]
    return df


# ---------------------------------------------------------------------------
# residual axes: the 2-D (signal_resid, noise_resid) plane
# ---------------------------------------------------------------------------
def _loess_fit(x, y, span=0.5, degree=2):
    """LOESS fitted values, returned in the caller's original point order."""
    from skmisc.loess import loess

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    order = np.argsort(x)
    m = loess(x[order], y[order], span=span, degree=degree)
    m.fit()
    out = np.empty_like(y)
    out[order] = m.outputs.fitted_values
    return out


def _z(v):
    v = np.asarray(v, dtype=float)
    sd = v.std(ddof=1)
    return (v - v.mean()) / sd if sd > 0 else v * 0.0


def _z_series(s):
    """``_z`` for a pandas groupby-transform (must return a Series, keeps index)."""
    sd = s.std(ddof=1)
    return (s - s.mean()) / sd if sd > 0 else s * 0.0


def add_residual_axes(df, span=0.5, suffix="", signal_col="signal", noise_col="noise"):
    """The two decoupled axes, plus the ``anchor``/``biology`` selectors built on them.

    ``R = signal / (signal + noise)`` is *not* a second axis: ``log noise`` rises
    as ``log signal`` with a slope ``b`` of only 0.17-0.32, so ``noise/signal ∝
    signal^(b-1)`` and 96-99% of the variance of ``log(noise/signal)`` is
    explained by ``signal`` alone (measured on three NBAtlas compartments; see
    ``anglemania_analysis/tasks/nmf/docs/tasks/pergene_signal_noise_R_axes.md``).
    ``rho(signal, R)`` is 0.98-0.997 and survives expression-trend correction.

    What does work is residualising ``noise`` *on* ``signal``, and doing it with
    a LOESS rather than a line -- the plane is curved, so OLS is
    Pearson-orthogonal by construction but still leaves rank structure behind
    (``spearman(log signal, resid)`` -0.31 to -0.42 for OLS vs -0.08 to -0.11
    for LOESS). Residualising the other way round (``signal ~ noise``) is not
    orthogonal to ``signal`` at all, and a symmetric PC2 is dominated by the ~3x
    variance imbalance between the two axes.

        signal_resid = log10 signal - LOESS(log10 signal | mean_lognorm)
        noise_resid  = log10 noise  - LOESS(log10 noise  | log10 signal)

    These two are mutually orthogonal (``rho`` -0.12 / +0.002 / -0.055 on the
    strict universe), so the 2-D framing is well posed, and

        anchor_score  = z(signal_resid) - z(noise_resid)   # structured AND batch-invariant
        biology_score = z(signal_resid) + z(noise_resid)   # structured AND patient-variable

    are its two diagonals. ``suffix`` names the column set (``"_mod"`` for the
    moderated-``S`` variant), so both can live in one table.
    """
    df = df.copy()
    s = f"log_signal{suffix}"
    n = f"log_noise{suffix}"
    df[s] = np.log10(np.maximum(df[signal_col].to_numpy(dtype=float), 1e-300))
    df[n] = np.log10(np.maximum(df[noise_col].to_numpy(dtype=float), 1e-300))

    df[f"signal_resid{suffix}"] = df[s] - _loess_fit(
        np.log10(df["mean_lognorm"].to_numpy(dtype=float) + 1e-6), df[s], span)
    df[f"noise_resid{suffix}"] = df[n] - _loess_fit(df[s], df[n], span)

    zs = _z(df[f"signal_resid{suffix}"])
    zn = _z(df[f"noise_resid{suffix}"])
    df[f"anchor_score{suffix}"] = zs - zn
    df[f"biology_score{suffix}"] = zs + zn
    return df


# ---------------------------------------------------------------------------
# limma-style empirical-Bayes moderation of the per-pair sd
# ---------------------------------------------------------------------------
def _trigamma_inverse(x):
    """Solve ``trigamma(y) = x`` for ``y > 0`` (limma's ``trigammaInverse``)."""
    from scipy.special import polygamma

    x = float(x)
    if not np.isfinite(x) or x <= 0:
        return np.inf
    if x > 1e7:
        return 1.0 / np.sqrt(x)
    if x < 1e-6:
        return 1.0 / x
    y = 0.5 + 1.0 / x
    for _ in range(50):
        tri = polygamma(1, y)
        dif = tri * (1.0 - tri / x) / polygamma(2, y)
        y += dif
        if abs(dif) / max(y, 1e-12) < 1e-8:
            break
    return y


def moderate_sd(sd_z, n_batches, mean_lognorm, n_bins=10, chunk=2000, inplace=True):
    """Shrink each pair's empirical sd toward a trend-predicted prior.

    ``S_ij`` is the sd of ``z_ij`` across batches, so it carries only ``B - 1``
    degrees of freedom -- 3 on a 4-batch simulation, 15 on pancreas. A pair whose
    ``S`` happens to come out small therefore looks far more reproducible than it
    is, and ``noise_i = sum_j S_ij^2`` sums ~p such estimates per gene. This is
    the same problem limma's ``eBayes`` solves for per-gene residual variances,
    and the same fix applies: assume ``s^2_ij ~ s0^2(x_ij) * chisq_d / d`` with
    ``s0^2`` drawn from a scaled inverse-chi-square with ``d0`` df, then use the
    posterior

        s^2_post = (d0 * s0^2(x) + d * s^2) / (d0 + d)

    ``d0`` is estimated globally by Smyth (2004) method-of-moments on
    ``e = log(s^2) - digamma(d/2) + log(d/2)``: ``trigamma(d0/2) = var(e) -
    trigamma(d/2)``, with ``d0 = inf`` (complete shrinkage to the trend) when
    that is non-positive.

    The prior *location* ``s0^2(x)`` is trended, because reliability of a
    correlation depends on how well both genes are measured. The trend covariate
    is the unordered pair of the two genes' mean-expression bins, which is
    non-parametric and needs no choice of how to combine the two means -- unlike
    ``min`` or ``mean``, which would each be an assumption. ``n_bins=10`` gives
    55 pair cells, each holding ~1e5-1e6 pairs at realistic gene counts.

    Modifies ``sd_z`` in place by default (it is a full genes x genes float64
    array). The diagonal is expected to be zeroed already by
    ``gene_level_scores``; it is excluded from the moment estimates and left
    untouched.
    """
    from scipy.special import digamma, polygamma

    d = float(n_batches - 1)
    if d < 2:
        raise ValueError(f"need >=3 batches to moderate S, got {n_batches}")
    p = sd_z.shape[0]

    gene_bin = pd.qcut(pd.Series(np.asarray(mean_lognorm, dtype=float)).rank(method="first"),
                       n_bins, labels=False, duplicates="drop").to_numpy()
    n_bins = int(gene_bin.max()) + 1
    # unordered pair of bins -> a single cell id, so the trend is symmetric in i, j
    cell_id = np.zeros((n_bins, n_bins), dtype=np.int64)
    k = 0
    for a in range(n_bins):
        for b in range(a, n_bins):
            cell_id[a, b] = cell_id[b, a] = k
            k += 1
    n_cells = k
    off = digamma(d / 2.0) - np.log(d / 2.0)

    # --- pass 1: moments of e, globally and per trend cell ---
    tot_n = 0
    tot_s = 0.0
    tot_s2 = 0.0
    cell_n = np.zeros(n_cells, dtype=np.int64)
    cell_s = np.zeros(n_cells, dtype=np.float64)
    for s0 in range(0, p, chunk):
        e0 = min(s0 + chunk, p)
        blk = sd_z[s0:e0]
        s2 = blk * blk
        e = np.log(s2, out=np.full_like(s2, -np.inf), where=s2 > 0) - off
        ok = np.isfinite(e)
        ok[np.arange(e0 - s0), np.arange(s0, e0)] = False  # diagonal
        cells = cell_id[np.ix_(gene_bin[s0:e0], gene_bin)]
        ev = e[ok]
        tot_n += ev.size
        tot_s += float(ev.sum())
        tot_s2 += float(np.dot(ev, ev))
        cell_n += np.bincount(cells[ok], minlength=n_cells)
        cell_s += np.bincount(cells[ok], weights=ev, minlength=n_cells)
        del blk, s2, e, ok, cells, ev

    if tot_n < 100:
        raise ValueError("too few finite pair variances to moderate")
    var_e = tot_s2 / tot_n - (tot_s / tot_n) ** 2
    target = var_e - float(polygamma(1, d / 2.0))
    d0 = 2.0 * _trigamma_inverse(target) if target > 0 else np.inf

    # empty/singleton cells fall back to the global prior location
    mean_e_global = tot_s / tot_n
    mean_e = np.where(cell_n > 0, cell_s / np.maximum(cell_n, 1), mean_e_global)
    # E[e] = log s0^2 + log(d0/2) - digamma(d0/2), so invert that. As d0 -> inf
    # the correction vanishes (digamma(x) -> log x) and the prior is exp(mean_e).
    if np.isfinite(d0):
        log_s0_2 = mean_e + digamma(d0 / 2.0) - np.log(d0 / 2.0)
    else:
        log_s0_2 = mean_e
    s0_2 = np.exp(log_s0_2)

    info = {"d": d, "d0": float(d0), "var_e": float(var_e),
            "trigamma_d_over_2": float(polygamma(1, d / 2.0)),
            "n_pairs": int(tot_n), "n_trend_cells": int(n_cells),
            "s0_2_min": float(s0_2.min()), "s0_2_max": float(s0_2.max()),
            "s0_2_median": float(np.median(s0_2)),
            "shrinkage_weight_prior": (1.0 if not np.isfinite(d0) else d0 / (d0 + d))}
    log(f"moderate_sd: d={d:.0f}, d0={d0:.3g}, prior weight "
        f"{info['shrinkage_weight_prior']:.3f}, s0^2 trend range "
        f"{s0_2.min():.3g}-{s0_2.max():.3g} over {n_cells} cells")

    # --- pass 2: posterior sd ---
    out = sd_z if inplace else np.empty_like(sd_z)
    for s0 in range(0, p, chunk):
        e0 = min(s0 + chunk, p)
        blk = sd_z[s0:e0]
        s2 = blk * blk
        prior = s0_2[cell_id[np.ix_(gene_bin[s0:e0], gene_bin)]]
        if np.isfinite(d0):
            post = (d0 * prior + d * s2) / (d0 + d)
        else:
            post = prior
        np.sqrt(post, out=post)
        post[np.arange(e0 - s0), np.arange(s0, e0)] = 0.0  # keep diagonal zeroed
        out[s0:e0] = post
        del blk, s2, prior, post
    return out, info


def score_genes(
    adata,
    *,
    batch_key,
    dataset_key=None,
    layer=None,
    method="cosine",
    normalization_mode="classical",
    permutation_function="sample",
    min_cells_per_gene=1,
    min_samples_per_gene=2,
    allow_missing_features=False,
    cell_chunk_size=0,
    seed=1,
    gpu=False,
    n_bins=20,
    weight_signal=0.5,
    weight_R=0.5,
    loess_span=0.5,
    moderate=False,
    moderate_n_bins=10,
):
    """End-to-end: cross-batch M/S -> per-gene signal/R -> ranking scores.

    Returns a per-gene DataFrame with ``score_binned`` / ``score_unbinned``
    (higher is better) alongside every intermediate, one row per gene that
    passed anglemania's all-batch intersection.

    ``add_residual_axes`` is always applied, so ``signal_resid``,
    ``noise_resid``, ``anchor_score`` and ``biology_score`` are present too.
    With ``moderate=True`` the per-pair sd is additionally shrunk toward a
    trended prior (``moderate_sd``) and the whole ``noise``-derived half of the
    table is recomputed from it into ``*_mod`` columns -- ``signal`` is
    untouched, since moderation is a statement about ``S`` only.
    """
    mean_z, sd_z, common_genes, n_batches_present = compute_MS(
        adata,
        batch_key=batch_key, dataset_key=dataset_key, layer=layer, method=method,
        normalization_mode=normalization_mode, permutation_function=permutation_function,
        min_cells_per_gene=min_cells_per_gene, min_samples_per_gene=min_samples_per_gene,
        allow_missing_features=allow_missing_features,
        cell_chunk_size=cell_chunk_size, seed=seed, gpu=gpu,
    )
    log("cross-batch M/S done; reducing to gene level")
    scores = gene_level_scores(mean_z, sd_z)

    mean_log, det = expression_covariates(adata, common_genes, layer)
    df = pd.DataFrame({"gene": common_genes, "mean_lognorm": mean_log,
                       "detection_rate": det,
                       "n_batches_present": n_batches_present, **scores})
    df = add_binned_zscores(df, n_bins, weight_signal, weight_R)
    df = add_residual_axes(df, span=loess_span)

    if moderate:
        # compute_MS already ran add_unique_batch_key on adata, so this is the
        # same batch partition the accumulator saw.
        n_batches = len(split_obs_indices_by_batch(adata))
        log(f"moderating S ({n_batches} batches)")
        sd_mod, info = moderate_sd(sd_z, n_batches, df["mean_lognorm"].to_numpy(),
                                   n_bins=moderate_n_bins, inplace=True)
        df["noise_mod"] = _row_sumsq(sd_mod, symmetrize=True)
        df["noise_row_mod"] = _row_sumsq(sd_mod)
        df["noise_col_mod"] = _row_sumsq(sd_mod.T)
        df["R_mod"] = df["signal"] / (df["signal"] + df["noise_mod"])
        # same binned score as `anglgene`, but with the moderated R
        df["R_z_binned_mod"] = df.groupby("expr_bin")["R_mod"].transform(_z_series)
        df["score_binned_mod"] = (weight_signal * df["signal_z_binned"]
                                  + weight_R * df["R_z_binned_mod"])
        df = add_residual_axes(df, span=loess_span, suffix="_mod",
                               signal_col="signal", noise_col="noise_mod")
        df.attrs["moderation"] = info
        del sd_mod

    del mean_z, sd_z
    return df


# ---------------------------------------------------------------------------
# the selection variants the benchmark runs as separate arms
# ---------------------------------------------------------------------------
# tag -> (score column, ascending, needs the moderated-S pass). Defined here so
# prepare_inputs.py (pipeline) and make_genelevel_variant_lists.py (one-shot,
# shares a single M/S computation across all arms) cannot drift apart.
VARIANT_SPECS = {
    "anglgene":             ("score_binned",       False, False),
    "anglgene_nobin":       ("score_unbinned",     False, False),
    # the 2-D plane's two diagonals
    "anglgene_anchor":      ("anchor_score",       False, False),
    "anglgene_biology":     ("biology_score",      False, False),
    # the second axis on its own, both directions: low = patient-invariant
    # relational context, high = patient-heterogeneous / HVG-like
    "anglgene_noiselo":     ("noise_resid",        True,  False),
    "anglgene_noisehi":     ("noise_resid",        False, False),
    # same three scores, but with S empirically-Bayes moderated first
    "anglgene_mod":         ("score_binned_mod",   False, True),
    "anglgene_anchor_mod":  ("anchor_score_mod",   False, True),
    "anglgene_biology_mod": ("biology_score_mod",  False, True),
}

def top_genes(df, n_genes, binned=True, column=None, ascending=False):
    """The ``n_genes`` best genes by ``column`` (default: the binned score).

    ``ascending=True`` takes the *bottom* of the column instead -- used by the
    ``noiselo`` arm, which ranks on the low-``noise_resid`` direction.
    """
    col = column or ("score_binned" if binned else "score_unbinned")
    if col not in df.columns:
        raise KeyError(f"no score column {col!r}; have {sorted(df.columns)}")
    return (df.sort_values(col, ascending=ascending)["gene"]
              .head(min(n_genes, len(df))).tolist())


# ---------------------------------------------------------------------------
# combining R with HVG: the anglemania intersection, and HVG's own per-gene score
# ---------------------------------------------------------------------------
def intersect_gene_names(
    adata,
    *,
    batch_key,
    dataset_key=None,
    layer=None,
    min_cells_per_gene=1,
    min_samples_per_gene=2,
):
    """anglemania's all-batch gene intersection, without computing any angles.

    ``compute_MS`` only ever scores genes that clear this, so any selection
    built on ``R`` is implicitly restricted to it. Exposing it on its own makes
    the "HVG on the intersection, no R" control possible -- without that
    control, a selection that filters on ``R`` is confounded with simply having
    dropped every gene that isn't detected in all batches (30% of genes on
    pancreas, 0.3% on the simulation).
    """
    from scipy import sparse as sp

    add_unique_batch_key(adata, batch_key, dataset_key)
    batch_indices = split_obs_indices_by_batch(adata)
    X_full = adata.X if layer is None else adata.layers[layer]
    var_names = np.asarray(adata.var_names)

    batch_genes = []
    for idx in batch_indices.values():
        X_b = X_full[idx]
        if sp.issparse(X_b):
            X_b = X_b.tocsc()
        mask = genes_passing_min_cells(X_b, min_cells_per_gene, np, sp)
        batch_genes.append(list(var_names[mask]))
    return np.asarray(
        intersect_genes(batch_genes, False, min_samples_per_gene, verbose=True)
    )


def hvg_dispersion(adata, batch_key, n_top_genes=2000, genes=None):
    """Per-gene normalized dispersion from scanpy's ``seurat`` HVG flavor.

    Returns a DataFrame with ``gene`` and ``dispersions_norm``, computed on a
    copy so the caller's ``adata`` keeps its raw counts.

    ``dispersions_norm`` is the continuous score the ``hvg`` selection ranks on,
    and it is constructed the same way as this module's ``signal_z_binned`` /
    ``R_z_binned``: scanpy's ``seurat`` flavor bins genes by mean expression and
    z-scores dispersion *within* each bin. So adding it to ``R_z_binned`` adds
    two quantities on the same scale, both already stripped of the expression
    confound -- which is what makes the R-weighted variant well posed.

    With a ``batch_key`` scanpy computes this per batch and averages the
    normalized dispersions, so the score is batch-aware to begin with.
    ``n_top_genes`` only sets the ``highly_variable`` flag, which is ignored
    here; the returned column does not depend on it.
    """
    import scanpy as sc

    sub = adata[:, genes].copy() if genes is not None else adata.copy()
    sc.pp.normalize_total(sub, target_sum=1e4)
    sc.pp.log1p(sub)
    sc.pp.highly_variable_genes(
        sub, flavor="seurat", n_top_genes=min(n_top_genes, sub.n_vars),
        batch_key=batch_key,
    )
    disp = sub.var["dispersions_norm"].to_numpy(dtype=float)
    # a gene with zero variance in some batch can come back NaN; rank it last
    disp = np.where(np.isfinite(disp), disp, -np.inf)
    return pd.DataFrame({"gene": np.asarray(sub.var_names), "dispersions_norm": disp})
