"""
Thin Python wrappers around kmer_seed.so and wfa_align.so.

Usage
-----
    from algorithms.c_libs import KmerLib, WFALib
    from pathlib import Path

    kmer = KmerLib(Path('algorithms/kmer_seed.so'))
    wfa  = WFALib(Path('algorithms/wfa_align.so'))
"""

from __future__ import annotations

import ctypes
import os
import subprocess
from pathlib import Path
from typing import List, Tuple

import numpy as np

from utilities import log


_ALGO_DIR     = Path(__file__).parent
_BUILD_SCRIPT = _ALGO_DIR / 'build_feature_align.sh'


def _ensure_built(so_path: Path) -> None:
    """Compile the .so files via the build script if not already present.

    """
    if so_path.exists():
        return
    if os.environ.get('HIREPEAT_IN_CONTAINER'):
        raise RuntimeError(
            f"{so_path} is missing from the container image, which should be "
            "impossible — `make` runs at image build time."
        )
    log.step(f"  Building {so_path.name} ...")
    r = subprocess.run(
        ['bash', str(_BUILD_SCRIPT)],
        capture_output=True, text=True, cwd=str(_ALGO_DIR),
    )
    if r.returncode != 0:
        raise RuntimeError(f"Build failed:\n{r.stderr}")


# ---------------------------------------------------------------------------
# KmerLib — kmer_seed.so
# ---------------------------------------------------------------------------

class KmerLib:
    """
    Wrapper for kmer_seed.so.

    Methods
    -------
    emit_pairs_pos(seqs, lens, kmer_size, max_bucket) -> (N, 4) int32
    seed_and_anchor_v2(...) -> chained anchors for one pair
    """

    def __init__(self, so_path: Path) -> None:
        _ensure_built(so_path)
        lib = ctypes.CDLL(str(so_path))

        # void kmer_free(void *ptr)
        lib.kmer_free.restype  = None
        lib.kmer_free.argtypes = [ctypes.c_void_p]

        # int kmer_emit_pairs_pos_alloc(seqs, lens, n, k, max_bucket,
        #                               **out_i, **out_j, **out_pi, **out_pj)
        lib.kmer_emit_pairs_pos_alloc.restype  = ctypes.c_int
        lib.kmer_emit_pairs_pos_alloc.argtypes = [
            ctypes.POINTER(ctypes.c_char_p),                # seqs
            ctypes.POINTER(ctypes.c_int),                   # lens
            ctypes.c_int,                                   # n_seqs
            ctypes.c_int,                                   # kmer_size
            ctypes.c_int,                                   # max_bucket
            ctypes.POINTER(ctypes.POINTER(ctypes.c_int)),   # out_i
            ctypes.POINTER(ctypes.POINTER(ctypes.c_int)),   # out_j
            ctypes.POINTER(ctypes.POINTER(ctypes.c_int)),   # out_pi
            ctypes.POINTER(ctypes.POINTER(ctypes.c_int)),   # out_pj
        ]

        # int kmer_seed_and_anchor_v2_alloc(...)
        lib.kmer_seed_and_anchor_v2_alloc.restype  = ctypes.c_int
        lib.kmer_seed_and_anchor_v2_alloc.argtypes = [
            ctypes.POINTER(ctypes.c_char_p),                # seqs
            ctypes.POINTER(ctypes.c_int),                   # lens
            ctypes.c_int,                                   # n_seqs
            ctypes.c_int,                                   # kmer_size
            ctypes.c_int,                                   # max_bucket
            ctypes.c_int,                                   # min_shared
            ctypes.c_float,                                 # max_len_ratio
            ctypes.c_int,                                   # filter_n
            ctypes.c_int,                                   # filter_mode
            ctypes.POINTER(ctypes.POINTER(ctypes.c_int)),   # out_i
            ctypes.POINTER(ctypes.POINTER(ctypes.c_int)),   # out_j
            ctypes.POINTER(ctypes.POINTER(ctypes.c_int)),   # out_pi
            ctypes.POINTER(ctypes.POINTER(ctypes.c_int)),   # out_pj
            ctypes.POINTER(ctypes.POINTER(ctypes.c_int)),   # out_rc
            ctypes.POINTER(ctypes.c_long),                  # out_n_raw_events
            ctypes.POINTER(ctypes.c_long),                  # out_n_unique_pairs
            ctypes.POINTER(ctypes.c_int),                   # out_bucket_hist
            ctypes.POINTER(ctypes.c_int),                   # out_hist_overflow
        ]

        self._lib = lib

    # ------------------------------------------------------------------

    def seed_and_anchor_v2(
        self,
        seqs:           List[bytes],
        lens:           List[int],
        kmer_size:      int,
        max_bucket:     int,
        min_shared:     int,
        max_len_ratio:  float = 0.0,
        filter_n:       int   = 0,
        filter_mode:    int   = 0,
    ) -> Tuple[np.ndarray, int, int, np.ndarray, int]:
        """
        fused k-mer seeding + anchor selection).

        Returns
        -------
        anchors : np.ndarray, shape (N, 5), dtype int32
            Columns: [i, j, anchor_pos_i, anchor_pos_j, is_rc].
        n_raw_events : int
        n_unique_pairs : int
        bucket_hist : np.ndarray, shape (max_bucket+2,), dtype int32
        hist_overflow : int
        """
        n = len(seqs)
        seq_arr = (ctypes.c_char_p * n)(*seqs)
        len_arr = (ctypes.c_int    * n)(*lens)
        hist_cap = max_bucket + 2

        ptr_i  = ctypes.POINTER(ctypes.c_int)()
        ptr_j  = ctypes.POINTER(ctypes.c_int)()
        ptr_pi = ctypes.POINTER(ctypes.c_int)()
        ptr_pj = ctypes.POINTER(ctypes.c_int)()
        ptr_rc = ctypes.POINTER(ctypes.c_int)()
        n_raw    = ctypes.c_long(0)
        n_unique = ctypes.c_long(0)
        hist_buf = (ctypes.c_int * hist_cap)()
        hist_overflow = ctypes.c_int(0)

        n_anchors = self._lib.kmer_seed_and_anchor_v2_alloc(
            seq_arr, len_arr, n, kmer_size, max_bucket, min_shared,
            ctypes.c_float(max_len_ratio),
            filter_n, filter_mode,
            ctypes.byref(ptr_i), ctypes.byref(ptr_j),
            ctypes.byref(ptr_pi), ctypes.byref(ptr_pj), ctypes.byref(ptr_rc),
            ctypes.byref(n_raw),
            ctypes.byref(n_unique),
            hist_buf,
            ctypes.byref(hist_overflow),
        )

        bucket_hist = np.frombuffer(hist_buf, dtype=np.int32, count=hist_cap).copy()

        if n_anchors == 0:
            return (np.empty((0, 5), dtype=np.int32), int(n_raw.value),
                    int(n_unique.value), bucket_hist, int(hist_overflow.value))

        result = np.empty((n_anchors, 5), dtype=np.int32)
        result[:, 0] = np.ctypeslib.as_array(ptr_i,  shape=(n_anchors,))
        result[:, 1] = np.ctypeslib.as_array(ptr_j,  shape=(n_anchors,))
        result[:, 2] = np.ctypeslib.as_array(ptr_pi, shape=(n_anchors,))
        result[:, 3] = np.ctypeslib.as_array(ptr_pj, shape=(n_anchors,))
        result[:, 4] = np.ctypeslib.as_array(ptr_rc, shape=(n_anchors,))

        self._lib.kmer_free(ptr_i)
        self._lib.kmer_free(ptr_j)
        self._lib.kmer_free(ptr_pi)
        self._lib.kmer_free(ptr_pj)
        self._lib.kmer_free(ptr_rc)

        return (result, int(n_raw.value), int(n_unique.value),
                bucket_hist, int(hist_overflow.value))

    def emit_pairs_pos(
        self,
        seqs:       List[bytes],
        lens:       List[int],
        kmer_size:  int,
        max_bucket: int,
    ) -> np.ndarray:
        """
        Emit canonical k-mer seed pairs with hit positions.

        Returns
        -------
        np.ndarray, shape (N, 4), dtype int32
            Columns: local_i, local_j, pos_in_i, pos_in_j.
        """
        n = len(seqs)
        seq_arr = (ctypes.c_char_p * n)(*seqs)
        len_arr = (ctypes.c_int    * n)(*lens)

        ptr_i  = ctypes.POINTER(ctypes.c_int)()
        ptr_j  = ctypes.POINTER(ctypes.c_int)()
        ptr_pi = ctypes.POINTER(ctypes.c_int)()
        ptr_pj = ctypes.POINTER(ctypes.c_int)()

        n_ev = self._lib.kmer_emit_pairs_pos_alloc(
            seq_arr, len_arr, n, kmer_size, max_bucket,
            ctypes.byref(ptr_i), ctypes.byref(ptr_j),
            ctypes.byref(ptr_pi), ctypes.byref(ptr_pj),
        )

        if n_ev == 0:
            return np.empty((0, 4), dtype=np.int32)

        result = np.empty((n_ev, 4), dtype=np.int32)
        result[:, 0] = np.ctypeslib.as_array(ptr_i,  shape=(n_ev,))
        result[:, 1] = np.ctypeslib.as_array(ptr_j,  shape=(n_ev,))
        result[:, 2] = np.ctypeslib.as_array(ptr_pi, shape=(n_ev,))
        result[:, 3] = np.ctypeslib.as_array(ptr_pj, shape=(n_ev,))
        self._lib.kmer_free(ptr_i)
        self._lib.kmer_free(ptr_j)
        self._lib.kmer_free(ptr_pi)
        self._lib.kmer_free(ptr_pj)
        return result


# ---------------------------------------------------------------------------
# WFALib — wfa_align.so
# ---------------------------------------------------------------------------

class WFALib:
    """
    Wrapper for wfa_align.so.

    Methods
    -------
    align_batch(seqs, lens, pair_i, pair_j, max_edit_frac, min_coverage)
        -> (out_id float32, out_al int32) numpy arrays

    align_batch_anchored(seqs, lens, pair_i, pair_j, api, apj,
                         max_edit_frac, min_coverage)
        -> (out_id float32, out_al int32) numpy arrays

    align_cigar(a, b, max_edit=0, min_coverage=0.0)
        -> (rc: int, identity: float, cigar: str)

    align_cigar_anchored(a, b, anchor_a, anchor_b, max_edit=0, min_coverage=0.0)
        -> (rc: int, identity: float, cigar: str)
    """

    def __init__(self, so_path: Path) -> None:
        _ensure_built(so_path)
        lib = ctypes.CDLL(str(so_path))

        # void wfa_align_batch(seqs, lens, pi, pj, n, max_edit_frac, min_cov,
        #                      *out_id, *out_al)
        lib.wfa_align_batch.restype  = None
        lib.wfa_align_batch.argtypes = [
            ctypes.POINTER(ctypes.c_char_p),  # seqs
            ctypes.POINTER(ctypes.c_int),     # lens
            ctypes.POINTER(ctypes.c_int),     # pair_i
            ctypes.POINTER(ctypes.c_int),     # pair_j
            ctypes.c_int,                     # n_pairs
            ctypes.c_float,                   # max_edit_frac
            ctypes.c_float,                   # min_coverage
            ctypes.POINTER(ctypes.c_float),   # out_identity
            ctypes.POINTER(ctypes.c_int),     # out_aln_len
        ]

        # int wfa_align_cigar(a, la, b, lb, max_edit, min_cov,
        #                     *out_id, *out_al, *cigar, cap)
        lib.wfa_align_cigar.restype  = ctypes.c_int
        lib.wfa_align_cigar.argtypes = [
            ctypes.c_char_p,                  # a
            ctypes.c_int,                     # la
            ctypes.c_char_p,                  # b
            ctypes.c_int,                     # lb
            ctypes.c_int,                     # max_edit
            ctypes.c_float,                   # min_coverage
            ctypes.POINTER(ctypes.c_float),   # out_identity
            ctypes.POINTER(ctypes.c_int),     # out_aln_len
            ctypes.c_char_p,                  # out_cigar (pre-allocated)
            ctypes.c_int,                     # cigar_cap
        ]

        # void wfa_align_batch_anchored(seqs, lens, pi, pj, api, apj, n,
        #                               frac, cov, *out_id, *out_al)
        lib.wfa_align_batch_anchored.restype  = None
        lib.wfa_align_batch_anchored.argtypes = [
            ctypes.POINTER(ctypes.c_char_p),  # seqs
            ctypes.POINTER(ctypes.c_int),     # lens
            ctypes.POINTER(ctypes.c_int),     # pair_i
            ctypes.POINTER(ctypes.c_int),     # pair_j
            ctypes.POINTER(ctypes.c_int),     # anchor_pi
            ctypes.POINTER(ctypes.c_int),     # anchor_pj
            ctypes.c_int,                     # n_pairs
            ctypes.c_float,                   # max_edit_frac
            ctypes.c_float,                   # min_coverage
            ctypes.POINTER(ctypes.c_float),   # out_identity
            ctypes.POINTER(ctypes.c_int),     # out_aln_len
        ]

        # int wfa_align_cigar_anchored(a, la, b, lb, anch_a, anch_b,
        #                              max_edit, min_cov,
        #                              *out_id, *out_al, *cigar, cap)
        lib.wfa_align_cigar_anchored.restype  = ctypes.c_int
        lib.wfa_align_cigar_anchored.argtypes = [
            ctypes.c_char_p,                  # a
            ctypes.c_int,                     # la
            ctypes.c_char_p,                  # b
            ctypes.c_int,                     # lb
            ctypes.c_int,                     # anchor_a
            ctypes.c_int,                     # anchor_b
            ctypes.c_int,                     # max_edit
            ctypes.c_float,                   # min_coverage
            ctypes.POINTER(ctypes.c_float),   # out_identity
            ctypes.POINTER(ctypes.c_int),     # out_aln_len
            ctypes.c_char_p,                  # out_cigar
            ctypes.c_int,                     # cigar_cap
        ]

        self._lib = lib

    # ------------------------------------------------------------------

    @staticmethod
    def _to_int32(arr, m: int):
        """Convert list or numpy array to contiguous int32 numpy array."""
        if isinstance(arr, np.ndarray):
            return np.ascontiguousarray(arr, dtype=np.int32)
        return np.array(arr, dtype=np.int32)

    def align_batch(
        self,
        seqs:          List[bytes],
        lens:          List[int],
        pair_i,        # List[int] or np.ndarray int32
        pair_j,        # List[int] or np.ndarray int32
        max_edit_frac: float,
        min_coverage:  float,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Unanchored WFA on *m* pairs.

        Returns
        -------
        out_id : np.ndarray float32, shape (m,)  — identity; -1 on failure
        out_al : np.ndarray int32,   shape (m,)  — alignment length; 0 on failure
        """
        n = len(seqs)
        pi_np = self._to_int32(pair_i, 0)
        pj_np = self._to_int32(pair_j, 0)
        m = len(pi_np)
        seq_arr = (ctypes.c_char_p * n)(*seqs)
        len_arr = (ctypes.c_int    * n)(*lens)
        out_id  = np.zeros(m, dtype=np.float32)
        out_al  = np.zeros(m, dtype=np.int32)

        self._lib.wfa_align_batch(
            seq_arr, len_arr,
            pi_np.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            pj_np.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            m,
            ctypes.c_float(max_edit_frac),
            ctypes.c_float(min_coverage),
            out_id.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            out_al.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        )
        return out_id, out_al

    def align_batch_anchored(
        self,
        seqs:          List[bytes],
        lens:          List[int],
        pair_i,        # List[int] or np.ndarray int32
        pair_j,        # List[int] or np.ndarray int32
        api,           # List[int] or np.ndarray int32
        apj,           # List[int] or np.ndarray int32
        max_edit_frac: float,
        min_coverage:  float,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Anchored WFA on *m* pairs.
        
        Returns
        -------
        out_id : np.ndarray float32, shape (m,)
        out_al : np.ndarray int32,   shape (m,)
        """
        n = len(seqs)
        pi_np  = self._to_int32(pair_i, 0)
        pj_np  = self._to_int32(pair_j, 0)
        api_np = self._to_int32(api, 0)
        apj_np = self._to_int32(apj, 0)
        m = len(pi_np)
        seq_arr = (ctypes.c_char_p * n)(*seqs)
        len_arr = (ctypes.c_int    * n)(*lens)
        out_id  = np.zeros(m, dtype=np.float32)
        out_al  = np.zeros(m, dtype=np.int32)

        self._lib.wfa_align_batch_anchored(
            seq_arr, len_arr,
            pi_np.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            pj_np.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            api_np.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            apj_np.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            m,
            ctypes.c_float(max_edit_frac),
            ctypes.c_float(min_coverage),
            out_id.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            out_al.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        )
        return out_id, out_al

    def align_cigar(
        self,
        a:            bytes,
        b:            bytes,
        max_edit:     int   = 0,
        min_coverage: float = 0.0,
    ) -> Tuple[int, float, str]:
        """
        Single-pair WFA with CIGAR output.

        Returns
        -------
        rc       : int    — edit distance (≥ 0 on success, -1/-2 on failure)
        identity : float  — 1 - d/max(la,lb); -1.0 on failure
        cigar    : str    — RLE CIGAR string ('nM', 'nD', 'nI'); '' on failure
        """
        la, lb    = len(a), len(b)
        cap       = 12 * (la + lb) + 16
        cigar_buf = ctypes.create_string_buffer(cap)
        identity  = ctypes.c_float(-1.0)
        aln_len   = ctypes.c_int(0)

        rc = self._lib.wfa_align_cigar(
            a, la, b, lb, max_edit,
            ctypes.c_float(min_coverage),
            ctypes.byref(identity),
            ctypes.byref(aln_len),
            cigar_buf, cap,
        )
        if rc < 0:
            return rc, -1.0, ''
        return rc, identity.value, cigar_buf.value.decode('ascii')

    def align_cigar_anchored(
        self,
        a:            bytes,
        b:            bytes,
        anchor_a:     int,
        anchor_b:     int,
        max_edit:     int   = 0,
        min_coverage: float = 0.0,
    ) -> Tuple[int, float, str]:
        """
        Anchored single-pair WFA with CIGAR output.

        Returns
        -------
        rc, identity, cigar  (same semantics as align_cigar)
        """
        la, lb    = len(a), len(b)
        cap       = 12 * (la + lb) + 16
        cigar_buf = ctypes.create_string_buffer(cap)
        identity  = ctypes.c_float(-1.0)
        aln_len   = ctypes.c_int(0)

        rc = self._lib.wfa_align_cigar_anchored(
            a, la, b, lb, anchor_a, anchor_b, max_edit,
            ctypes.c_float(min_coverage),
            ctypes.byref(identity),
            ctypes.byref(aln_len),
            cigar_buf, cap,
        )
        if rc < 0:
            return rc, -1.0, ''
        return rc, identity.value, cigar_buf.value.decode('ascii')

