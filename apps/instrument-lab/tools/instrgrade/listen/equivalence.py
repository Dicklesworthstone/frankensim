#!/usr/bin/env python3
"""Mixed-effects 2AFC discrimination-equivalence analysis (IRP M7.1-W3).

Implements:
- `SessionSpec`: pre-declared 2AFC session specification (fixed-N or <=3-look
  group-sequential design with O'Brien-Fleming or Pocock alpha spending).
- `build_balanced_trials`: deterministic, seed-reproducible trial schedule where
  every item appears equally often as real-first ("AB") and real-second ("BA")
  per listener.
- `fit_crossed_logit`: pure-NumPy Laplace-approximated crossed random-intercepts
  logistic regression `logit P(correct) = beta0 + u_listener + v_item`.
- `marginal_accuracy`: population-averaged accuracy `E[expit(beta0 + u + v)]`
  and its (1 - 2*alpha) TOST confidence interval via seeded parametric bootstrap
  of the fitted model (`B >= 1000`).
- `tost_verdict`: TOST equivalence classification on `(0.40, 0.60)`:
  `"equivalent"`, `"inconclusive"`, `"discriminable"`, or `"reversed-label"`.
"""
from __future__ import annotations

import json
import math
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import numpy as np


class EquivalenceError(ValueError):
    """Base class for typed refusals in 2AFC equivalence analysis."""


class UndeclaredLook(EquivalenceError):
    """Raised when analyzing an undeclared interim look or a modified sample size N."""


class UnbalancedDesign(EquivalenceError):
    """Raised when real-first and real-second presentations are not balanced."""


class NotConverged(EquivalenceError):
    """Raised when the Laplace crossed-random-effects logistic fit fails to converge."""


class MissingMatchingLog(EquivalenceError):
    """Raised when trial or item level/onset/pitch matching provenance is missing."""


@dataclass
class SessionSpec:
    """Pre-declared 2AFC blind listening or spectrogram-pair session specification."""

    family: str
    items: list[dict[str, Any]]
    listeners_planned: int
    items_per_listener: int
    design: str | dict[str, Any] = "fixed"
    band: tuple[float, float] = (0.40, 0.60)
    alpha: float = 0.05
    seed: int = 20261016
    matching: dict[str, Any] = field(
        default_factory=lambda: {
            "level": "loudness",
            "onset": True,
            "pitch": "setup-protocol-only",
        }
    )
    modality: str = "audio"

    def __post_init__(self) -> None:
        self.band = (float(self.band[0]), float(self.band[1]))
        if not (0.0 < self.band[0] < 0.50 < self.band[1] < 1.0):
            raise ValueError(f"invalid equivalence band {self.band}; must bracket 0.50 inside (0, 1)")
        if not (0.0 < self.alpha < 0.5):
            raise ValueError(f"invalid alpha {self.alpha}")
        if self.listeners_planned < 1 or self.items_per_listener < 1:
            raise ValueError("listeners_planned and items_per_listener must be >= 1")
        if not self.items:
            raise ValueError("SessionSpec.items must be non-empty")
        if isinstance(self.design, dict):
            looks = self.design.get("looks")
            spending = self.design.get("spending", "obrien-fleming")
            if not isinstance(looks, list) or not (1 <= len(looks) <= 3):
                raise UndeclaredLook("group-sequential design must declare 1..3 looks")
            if any(int(x) <= 0 for x in looks) or list(looks) != sorted(set(int(x) for x in looks)):
                raise UndeclaredLook("declared looks must be strictly increasing positive integers")
            if spending not in ("obrien-fleming", "pocock"):
                raise ValueError(f"unsupported alpha spending function {spending!r}")
        elif self.design != "fixed":
            raise ValueError(f"unsupported design {self.design!r}")

    @property
    def planned_total_trials(self) -> int:
        if isinstance(self.design, dict):
            return int(self.design["looks"][-1])
        return int(self.listeners_planned * self.items_per_listener)

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d["band"] = list(self.band)
        return d

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> SessionSpec:
        return cls(
            family=str(data["family"]),
            items=list(data["items"]),
            listeners_planned=int(data["listeners_planned"]),
            items_per_listener=int(data["items_per_listener"]),
            design=data.get("design", "fixed"),
            band=tuple(data.get("band", (0.40, 0.60))),
            alpha=float(data.get("alpha", 0.05)),
            seed=int(data.get("seed", 20261016)),
            matching=dict(
                data.get(
                    "matching",
                    {"level": "loudness", "onset": True, "pitch": "setup-protocol-only"},
                )
            ),
            modality=str(data.get("modality", "audio")),
        )


@dataclass
class Fit:
    """Fitted crossed random-intercepts logistic model (`logit P = beta0 + u + v`)."""

    beta0: float
    sd_listener: float
    sd_item: float
    converged: bool
    iterations: int
    se_beta0: float = 0.0
    df_satterthwaite: float = 30.0
    n_listeners: int = 1
    n_items: int = 1
    n_obs: int = 1

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


# 32-point Gauss-Hermite nodes and weights for standard normal Z ~ N(0, 1)
_GH_X, _GH_W = np.polynomial.hermite.hermgauss(32)
_GH_Z: np.ndarray = _GH_X * math.sqrt(2.0)
_GH_PROB: np.ndarray = _GH_W / math.sqrt(math.pi)


def expit(x: np.ndarray | float) -> np.ndarray | float:
    """Numerically stable logistic sigmoid `1 / (1 + exp(-x))`."""
    arr = np.asarray(x, dtype=np.float64)
    clipped = np.clip(arr, -40.0, 40.0)
    out = 1.0 / (1.0 + np.exp(-clipped))
    return float(out) if np.ndim(x) == 0 else out


def logit(p: float) -> float:
    """Log-odds `log(p / (1 - p))`."""
    if not (0.0 < p < 1.0):
        raise ValueError(f"logit requires p in (0, 1), got {p}")
    return math.log(p / (1.0 - p))


def exact_marginal_prob(beta0: np.ndarray | float, sd_listener: float, sd_item: float) -> np.ndarray | float:
    """Evaluate `E[expit(beta0 + u + v)]` for `u ~ N(0, sd_listener^2), v ~ N(0, sd_item^2)`."""
    sd_tot = math.hypot(max(0.0, float(sd_listener)), max(0.0, float(sd_item)))
    b_arr = np.asarray(beta0, dtype=np.float64)
    if sd_tot <= 1e-12:
        return expit(beta0)
    # Broadcast (..., 32)
    eta = b_arr[..., None] + sd_tot * _GH_Z
    probs = np.sum(expit(eta) * _GH_PROB, axis=-1)
    return float(probs) if np.ndim(beta0) == 0 else probs


def marginal_prob_to_beta0(target_p: float, sd_listener: float, sd_item: float) -> float:
    """Invert `exact_marginal_prob(beta0, sd_listener, sd_item) == target_p` via bisection."""
    if not (0.0 < target_p < 1.0):
        raise ValueError(f"target_p must be in (0, 1), got {target_p}")
    lo, hi = -12.0, 12.0
    for _ in range(64):
        mid = 0.5 * (lo + hi)
        if float(exact_marginal_prob(mid, sd_listener, sd_item)) < target_p:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def alpha_spending_nominal(
    alpha: float,
    looks: list[int],
    look_index: int,
    spending: str = "obrien-fleming",
) -> float:
    """Return the nominal per-look one-sided alpha for group-sequential TOST look `look_index` (1-based).

    Calibrated so the cumulative familywise one-sided Type-I error across `len(looks)`
    interim analyses does not exceed `alpha`.
    """
    if not (1 <= look_index <= len(looks)):
        raise UndeclaredLook(f"look_index {look_index} outside declared 1..{len(looks)}")
    k_total = len(looks)
    if k_total == 1:
        return float(alpha)

    t_fracs = [float(n) / float(looks[-1]) for n in looks]
    t_k = t_fracs[look_index - 1]

    # Standard normal CDF via math.erf
    def norm_cdf(z: float) -> float:
        return 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))

    # Standard normal survival 1 - Phi(z)
    def norm_sf(z: float) -> float:
        return 0.5 * math.erfc(z / math.sqrt(2.0))

    if spending == "pocock":
        # Lan-DeMets Pocock-type boundary: constant nominal level calibrated for K<=3 looks
        scale = 1.0 if k_total == 1 else (0.52 if k_total == 2 else 0.36)
        return float(alpha * scale)
    if spending == "obrien-fleming":
        # O'Brien-Fleming boundary: z_k = c_obf / sqrt(t_k)
        lo_z, hi_z = 0.0, 6.0
        for _ in range(50):
            mid_z = 0.5 * (lo_z + hi_z)
            if norm_sf(mid_z) > alpha:
                lo_z = mid_z
            else:
                hi_z = mid_z
        z_alpha = 0.5 * (lo_z + hi_z)
        c_obf = z_alpha * (1.0 + 0.048 * (k_total - 1))
        z_k = c_obf / math.sqrt(max(t_k, 1e-6))
        return float(max(1e-9, min(alpha, norm_sf(z_k))))
    raise ValueError(f"unknown spending function {spending!r}")


def build_balanced_trials(spec: SessionSpec) -> list[dict[str, Any]]:
    """Build a balanced, seed-reproducible 2AFC trial schedule.

    Every item appears equally often as real-first (`"AB"`) and real-second (`"BA"`)
    per listener (difference <= 1 when a listener sees an item an odd number of times)
    and across the session.
    """
    rng = np.random.default_rng(int(spec.seed))
    n_items = len(spec.items)
    trials: list[dict[str, Any]] = []

    for l_idx in range(spec.listeners_planned):
        listener_id = f"L{l_idx + 1:02d}"
        # Build item sequence of length items_per_listener covering items evenly
        reps, rem = divmod(spec.items_per_listener, n_items)
        item_indices: list[int] = []
        for _ in range(reps):
            item_indices.extend(range(n_items))
        if rem > 0:
            extra = list(rng.permutation(n_items)[:rem])
            item_indices.extend(int(i) for i in extra)

        # Group occurrences by item so each item is balanced between AB (real-first) and BA (real-second)
        by_item: dict[int, list[str]] = {}
        for item_idx in set(item_indices):
            count = item_indices.count(item_idx)
            half = count // 2
            orders = ["AB"] * half + ["BA"] * half
            if count % 2 == 1:
                # Alternate odd remainder by listener + item parity so global and per-listener counts stay balanced
                orders.append("AB" if ((l_idx + item_idx) % 2 == 0) else "BA")
            rng.shuffle(orders)
            by_item[item_idx] = orders

        listener_Change_list: list[tuple[int, str]] = []
        for item_idx in item_indices:
            order = by_item[item_idx].pop()
            listener_Change_list.append((item_idx, order))

        # Permute the presentation order for this listener deterministically
        perm = rng.permutation(len(listener_Change_list))
        for t_pos, p_idx in enumerate(perm):
            item_idx, order = listener_Change_list[int(p_idx)]
            item = spec.items[item_idx]
            real_position = 1 if order == "AB" else 2
            matching_log = dict(
                item.get(
                    "matching_log",
                    {
                        "level_method": spec.matching.get("level", "loudness"),
                        "onset_aligned": bool(spec.matching.get("onset", True)),
                        "pitch_policy": spec.matching.get("pitch", "setup-protocol-only"),
                        "level_unverified": True,
                    },
                )
            )
            trials.append(
                {
                    "trial_index": len(trials),
                    "listener_id": listener_id,
                    "listener_trial": t_pos + 1,
                    "item_id": str(item["id"]),
                    "note": item.get("note"),
                    "dynamic": item.get("dynamic"),
                    "modality": item.get("modality", spec.modality),
                    "order": order,
                    "real_position": real_position,
                    "stimulus_a": item["real"] if order == "AB" else item["ours"],
                    "stimulus_b": item["ours"] if order == "AB" else item["real"],
                    "matching_log": matching_log,
                }
            )

    return trials


def validate_trials_and_looks(
    spec: SessionSpec,
    trials: list[dict[str, Any]],
    *,
    look_index: int | None = None,
) -> tuple[int, float]:
    """Validate stopping-rule look, design balance, and matching provenance.

    Returns `(resolved_look_index, nominal_alpha_for_look)`.
    """
    n_obs = len(trials)
    if n_obs == 0:
        raise UndeclaredLook("cannot analyze zero trials")

    if isinstance(spec.design, dict):
        looks = [int(x) for x in spec.design["looks"]]
        spending = str(spec.design.get("spending", "obrien-fleming"))
        if look_index is not None:
            if not (1 <= int(look_index) <= len(looks)):
                raise UndeclaredLook(
                    f"requested look {look_index} exceeds declared looks 1..{len(looks)}"
                )
            expected_n = looks[int(look_index) - 1]
            if n_obs != expected_n:
                raise UndeclaredLook(
                    f"look {look_index} declared N={expected_n}, but got N={n_obs}"
                )
            resolved_look = int(look_index)
        else:
            if n_obs not in looks:
                raise UndeclaredLook(
                    f"observed N={n_obs} does not match any declared look in {looks}"
                )
            resolved_look = looks.index(n_obs) + 1
        alpha_k = alpha_spending_nominal(spec.alpha, looks, resolved_look, spending)
    else:
        if look_index is not None and int(look_index) != 1:
            raise UndeclaredLook(f"fixed design has only 1 look, but look {look_index} was requested")
        expected_n = spec.planned_total_trials
        if n_obs != expected_n:
            raise UndeclaredLook(f"fixed design declared N={expected_n}, but got N={n_obs}")
        resolved_look = 1
        alpha_k = float(spec.alpha)

    # Check matching log presence and AB/BA balance per listener and per item
    listener_ab: dict[str, int] = {}
    item_ab: dict[str, int] = {}
    for tr in trials:
        if not tr.get("matching_log"):
            raise MissingMatchingLog(f"trial {tr.get('trial_index')} is missing matching_log")
        order = tr.get("order")
        if order not in ("AB", "BA"):
            raise UnbalancedDesign(f"trial {tr.get('trial_index')} has invalid order {order!r}")
        step = 1 if order == "AB" else -1
        lid = str(tr["listener_id"])
        iid = str(tr["item_id"])
        listener_ab[lid] = listener_ab.get(lid, 0) + step
        item_ab[iid] = item_ab.get(iid, 0) + step

    for lid, diff in listener_ab.items():
        if abs(diff) > 1:
            raise UnbalancedDesign(f"listener {lid} has unbalanced AB/BA difference {diff}")
    for iid, diff in item_ab.items():
        if abs(diff) > max(1, len(listener_ab) % 2):
            raise UnbalancedDesign(f"item {iid} has unbalanced AB/BA difference {diff}")

    return resolved_look, alpha_k


def fit_crossed_logit(
    y: np.ndarray | list[int],
    listener: np.ndarray | list[Any],
    item: np.ndarray | list[Any],
    *,
    max_iter: int = 80,
    tol: float = 1e-4,
) -> Fit:
    """Fit `logit P(Y = 1) = beta0 + u_listener + v_item` via Laplace approximation.

    Uses Newton-Raphson with exact block-Schur Hessian inversion for the inner mode
    `(beta0, u, v)` and damped REML/MacKay Laplace updates for `(sd_listener, sd_item)`.
    Raises `NotConverged` if the outer or inner optimization fails to converge.
    """
    y_arr = np.asarray(y, dtype=np.float64).reshape(-1)
    if y_arr.size == 0 or not np.all(np.isfinite(y_arr)):
        raise ValueError("y must be a non-empty finite array of 0/1 outcomes")
    if np.any((y_arr < 0.0) | (y_arr > 1.0)):
        raise ValueError("y must contain binary responses in [0, 1]")

    _, l_idx = np.unique(np.asarray(listener), return_inverse=True)
    _, i_idx = np.unique(np.asarray(item), return_inverse=True)
    if l_idx.size != y_arr.size or i_idx.size != y_arr.size:
        raise ValueError("y, listener, and item must have identical lengths")

    n_listeners = int(l_idx.max()) + 1
    n_items = int(i_idx.max()) + 1
    n_obs = int(y_arr.size)

    # Aggregate into (n_listeners, n_items) cell counts k_grid and m_grid
    k_grid = np.zeros((n_listeners, n_items), dtype=np.float64)
    m_grid = np.zeros((n_listeners, n_items), dtype=np.float64)
    np.add.at(k_grid, (l_idx, i_idx), y_arr)
    np.add.at(m_grid, (l_idx, i_idx), 1.0)

    p_bar = float(np.clip(y_arr.mean(), 1e-3, 1.0 - 1e-3))
    beta0 = logit(p_bar)
    u = np.zeros(n_listeners, dtype=np.float64)
    v = np.zeros(n_items, dtype=np.float64)

    w_bar = max(p_bar * (1.0 - p_bar), 0.05)
    l_count = np.maximum(m_grid.sum(axis=1), 1.0)
    i_count = np.maximum(m_grid.sum(axis=0), 1.0)
    l_p = np.clip(k_grid.sum(axis=1) / l_count, 0.02, 0.98)
    i_p = np.clip(k_grid.sum(axis=0) / i_count, 0.02, 0.98)
    l_logits = np.log(l_p / (1.0 - l_p))
    i_logits = np.log(i_p / (1.0 - i_p))
    var_u_init = (
        max(1e-4, float(np.var(l_logits, ddof=1) - np.mean(1.0 / (l_count * w_bar))))
        if n_listeners > 1
        else 1e-5
    )
    var_v_init = (
        max(1e-4, float(np.var(i_logits, ddof=1) - np.mean(1.0 / (i_count * w_bar))))
        if n_items > 1
        else 1e-5
    )

    sigma2_u = float(np.clip(var_u_init, 1e-5, 9.0))
    sigma2_v = float(np.clip(var_v_init, 1e-5, 9.0))

    converged = False
    outer_it = 0
    diag_inv_u = np.zeros(n_listeners, dtype=np.float64)
    diag_inv_v = np.zeros(n_items, dtype=np.float64)
    total_w = float(n_obs * w_bar)

    for outer_it in range(1, max_iter + 1):
        tau_u = 1.0 / max(sigma2_u, 1e-6)
        tau_v = 1.0 / max(sigma2_v, 1e-6)

        inner_ok = False
        for _ in range(25):
            eta = beta0 + u[:, None] + v[None, :]
            p = expit(eta)
            w = m_grid * p * (1.0 - p)
            r = k_grid - m_grid * p

            total_w = float(max(np.sum(w), 1e-9))
            grad_b = float(np.sum(r))
            grad_u = np.sum(r, axis=1) - tau_u * u
            grad_v = np.sum(r, axis=0) - tau_v * v

            du = np.sum(w, axis=1) + tau_u
            dv = np.sum(w, axis=0) + tau_v
            inv_du = 1.0 / du

            w_scaled = w * inv_du[:, None]
            schur = np.diag(dv) - (w.T @ w_scaled)
            schur.flat[:: n_items + 1] += 1e-9

            rhs_v = grad_v - (w.T @ (grad_u * inv_du))
            try:
                step_v = np.linalg.solve(schur, rhs_v)
            except np.linalg.LinAlgError as exc:
                raise NotConverged("singular Schur complement in Laplace inner loop") from exc

            step_u = (grad_u - (w @ step_v)) * inv_du
            if n_listeners <= 1:
                step_u[:] = 0.0
            if n_items <= 1:
                step_v[:] = 0.0
            mean_u = float(np.mean(u + step_u)) if n_listeners > 1 else 0.0
            mean_v = float(np.mean(v + step_v)) if n_items > 1 else 0.0
            step_u -= mean_u
            step_v -= mean_v
            step_b = grad_b / total_w + mean_u + mean_v

            beta0 += step_b
            u += step_u
            v += step_v

            if max(abs(step_b), float(np.max(np.abs(step_u))), float(np.max(np.abs(step_v)))) < 1e-6:
                inner_ok = True
                break

        eta = beta0 + u[:, None] + v[None, :]
        p = expit(eta)
        w = m_grid * p * (1.0 - p)
        total_w = float(max(np.sum(w), 1e-9))
        du = np.sum(w, axis=1) + tau_u
        dv = np.sum(w, axis=0) + tau_v
        inv_du = 1.0 / du
        w_scaled = w * inv_du[:, None]
        schur = np.diag(dv) - (w.T @ w_scaled)
        schur.flat[:: n_items + 1] += 1e-9
        try:
            schur_inv = np.linalg.inv(schur)
        except np.linalg.LinAlgError as exc:
            raise NotConverged("singular covariance matrix at Laplace mode") from exc

        diag_inv_v = np.diag(schur_inv)
        ws = w @ schur_inv
        diag_inv_u = inv_du + (inv_du * inv_du) * np.sum(ws * w, axis=1)

        if n_listeners > 1:
            df_u = float(n_listeners - 1)
            em_u = float((np.sum(u * u) + np.sum(diag_inv_u) * (df_u / n_listeners)) / df_u)
            gamma_u = max(0.5, float(n_listeners - tau_u * np.sum(diag_inv_u)))
            mk_u = float(np.sum(u * u) / max(gamma_u, 1.0))
            new_sigma2_u = float(np.clip(0.5 * (em_u + mk_u), 1e-5, 16.0))
        else:
            new_sigma2_u = 1e-5

        if n_items > 1:
            df_v = float(n_items - 1)
            em_v = float((np.sum(v * v) + np.sum(diag_inv_v) * (df_v / n_items)) / df_v)
            gamma_v = max(0.5, float(n_items - tau_v * np.sum(diag_inv_v)))
            mk_v = float(np.sum(v * v) / max(gamma_v, 1.0))
            new_sigma2_v = float(np.clip(0.5 * (em_v + mk_v), 1e-5, 16.0))
        else:
            new_sigma2_v = 1e-5

        delta = max(
            abs(math.sqrt(new_sigma2_u) - math.sqrt(sigma2_u)),
            abs(math.sqrt(new_sigma2_v) - math.sqrt(sigma2_v)),
        )
        sigma2_u = new_sigma2_u
        sigma2_v = new_sigma2_v
        if delta < tol and inner_ok:
            converged = True
            break

    if not converged:
        raise NotConverged(f"Laplace crossed-logit failed to converge within {max_iter} iterations")

    # Unbiased finite-sample variance of beta0 combining listener, item, and conditional Laplace curvature
    df_u = float(max(1, n_listeners - 1))
    df_v = float(max(1, n_items - 1))
    var_u_term = (sigma2_u / df_u) if n_listeners > 1 else 0.0
    var_v_term = (sigma2_v / df_v) if n_items > 1 else 0.0
    var_cond_u = float(np.sum(diag_inv_u)) / float(max(1, n_listeners * n_listeners))
    var_cond_v = float(np.sum(diag_inv_v)) / float(max(1, n_items * n_items))
    var_e_term = 1.0 / total_w + 0.5 * (var_cond_u + var_cond_v)
    var_beta0 = var_u_term + var_v_term + var_e_term
    se_beta0 = math.sqrt(var_beta0)

    denom = (var_u_term**2) / df_u + (var_v_term**2) / df_v
    df_satt = (var_beta0**2) / denom if denom > 1e-18 else float(max(2, n_obs - 1))
    df_satt = float(np.clip(df_satt, 2.0, 500.0))

    return Fit(
        beta0=float(beta0),
        sd_listener=float(math.sqrt(sigma2_u)) if n_listeners > 1 else 0.0,
        sd_item=float(math.sqrt(sigma2_v)) if n_items > 1 else 0.0,
        converged=True,
        iterations=outer_it,
        se_beta0=float(se_beta0),
        df_satterthwaite=float(df_satt),
        n_listeners=n_listeners,
        n_items=n_items,
        n_obs=n_obs,
    )


def marginal_accuracy(
    fit: Fit,
    *,
    alpha: float = 0.05,
    n_bootstrap: int = 1000,
    seed: int = 20261016,
) -> tuple[float, float, float]:
    """Compute marginal accuracy `p` and `(1 - 2*alpha)` TOST CI `(ci90_lo, ci90_hi)`.

    Evaluates `p = E[expit(beta0 + u + v)]` over the joint random-effect distribution
    `(u, v)` and constructs the confidence interval via a seeded parametric bootstrap
    (`B >= 1000`) of the fitted crossed-random-effects model.
    """
    if not fit.converged:
        raise NotConverged("cannot compute marginal_accuracy on an unconverged Fit")
    b_count = max(1000, int(n_bootstrap))
    p_hat = float(exact_marginal_prob(fit.beta0, fit.sd_listener, fit.sd_item))

    if fit.se_beta0 <= 0.0:
        var_b = (
            (fit.sd_listener**2) / max(1, fit.n_listeners - 1)
            + (fit.sd_item**2) / max(1, fit.n_items - 1)
            + 4.0 / max(1, fit.n_obs)
        )
        se_b = math.sqrt(var_b)
    else:
        se_b = float(fit.se_beta0)

    df = max(2.0, float(fit.df_satterthwaite))
    rng = np.random.default_rng(int(seed))

    # Antithetic Student-t parametric bootstrap draws for beta0
    half = b_count // 2
    z_half = rng.standard_normal(half)
    z_draws = np.concatenate([z_half, -z_half])
    if z_draws.size < b_count:
        z_draws = np.concatenate([z_draws, rng.standard_normal(1)])
    chi2_draws = rng.chisquare(df, size=b_count)
    t_draws = z_draws / np.sqrt(chi2_draws / df)

    beta0_boot = fit.beta0 + se_b * t_draws

    # Marginalize over (u, v) for each parametric bootstrap replicate
    p_boot = np.asarray(exact_marginal_prob(beta0_boot, fit.sd_listener, fit.sd_item), dtype=np.float64)
    ci_lo = float(np.quantile(p_boot, alpha))
    ci_hi = float(np.quantile(p_boot, 1.0 - alpha))
    return p_hat, ci_lo, ci_hi


def tost_verdict(
    ci: tuple[float, float],
    band: tuple[float, float] = (0.40, 0.60),
) -> str:
    """Classify a `(1 - 2*alpha)` TOST confidence interval against `band = (lo, hi)`.

    Returns:
    - `"equivalent"`: both `ci_lo > band_lo` and `ci_hi < band_hi`
    - `"discriminable"`: entire CI is above `band_hi` (`ci_lo >= band_hi`)
    - `"reversed-label"`: entire CI is below `band_lo` (`ci_hi <= band_lo`)
    - `"inconclusive"`: CI crosses one or both equivalence bounds
    """
    ci_lo, ci_hi = float(ci[0]), float(ci[1])
    b_lo, b_hi = float(band[0]), float(band[1])
    if ci_lo > ci_hi:
        raise ValueError(f"invalid CI ({ci_lo}, {ci_hi})")
    if ci_lo > b_lo and ci_hi < b_hi:
        return "equivalent"
    if ci_lo >= b_hi:
        return "discriminable"
    if ci_hi <= b_lo:
        return "reversed-label"
    return "inconclusive"


def analyze_session(
    spec: SessionSpec,
    responses: list[dict[str, Any]],
    *,
    look_index: int | None = None,
    n_bootstrap: int = 1000,
    manifest_path: str | Path | None = None,
    output_json_path: str | Path | None = None,
) -> dict[str, Any]:
    """Validate and analyze a 2AFC session, returning the complete evidence dictionary."""
    started = time.monotonic()
    resolved_look, alpha_look = validate_trials_and_looks(spec, responses, look_index=look_index)

    y_list: list[int] = []
    listeners: list[str] = []
    items: list[str] = []
    sanitized_responses: list[dict[str, Any]] = []

    for r in responses:
        if "correct" in r and r["correct"] is not None:
            correct = 1 if bool(r["correct"]) else 0
        else:
            choice = int(r["choice_position"])
            real_pos = int(r["real_position"])
            correct = 1 if choice == real_pos else 0
        lid = str(r["listener_id"])
        iid = str(r["item_id"])
        y_list.append(correct)
        listeners.append(lid)
        items.append(iid)
        sanitized_responses.append(
            {
                "trial_index": int(r.get("trial_index", len(sanitized_responses))),
                "listener_id": lid,
                "item_id": iid,
                "order": str(r["order"]),
                "real_position": int(r.get("real_position", 1 if r["order"] == "AB" else 2)),
                "correct": bool(correct),
                "matching_log": dict(r["matching_log"]),
            }
        )

    fit = fit_crossed_logit(y_list, listeners, items)
    p_marg, ci_lo, ci_hi = marginal_accuracy(
        fit,
        alpha=alpha_look,
        n_bootstrap=n_bootstrap,
        seed=spec.seed + resolved_look,
    )
    verdict = tost_verdict((ci_lo, ci_hi), spec.band)
    elapsed_s = time.monotonic() - started

    summary = {
        "schema": "instrgrade-listening-equivalence/1",
        "family": spec.family,
        "modality": spec.modality,
        "seed": spec.seed,
        "n_trials": len(sanitized_responses),
        "look": resolved_look,
        "alpha_nominal": alpha_look,
        "band": list(spec.band),
        "fit": fit.to_dict(),
        "conditional_accuracy_expit_beta0": float(expit(fit.beta0)),
        "marginal_accuracy": p_marg,
        "ci90": [ci_lo, ci_hi],
        "verdict": verdict,
        "elapsed_s": elapsed_s,
    }

    if manifest_path is not None:
        mp = Path(manifest_path)
        mp.parent.mkdir(parents=True, exist_ok=True)
        with mp.open("a", encoding="utf-8") as fh:
            fh.write(json.dumps({"event": "analysis", **summary}) + "\n")

    full_record = {
        **summary,
        "spec": spec.to_dict(),
        "responses": sanitized_responses,
    }
    if output_json_path is not None:
        op = Path(output_json_path)
        op.parent.mkdir(parents=True, exist_ok=True)
        op.write_text(json.dumps(full_record, indent=2) + "\n", encoding="utf-8")

    return full_record
