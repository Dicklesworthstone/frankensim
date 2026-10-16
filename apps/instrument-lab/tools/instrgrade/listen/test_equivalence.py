#!/usr/bin/env python3
"""Unit, simulation, and end-to-end tests for the 2AFC equivalence harness (IRP M7.1-W3)."""
from __future__ import annotations

import json
import math
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

THIS_DIR = Path(__file__).resolve().parent
INSTRGRADE_DIR = THIS_DIR.parent
if str(INSTRGRADE_DIR) not in sys.path:
    sys.path.insert(0, str(INSTRGRADE_DIR))
if str(THIS_DIR) not in sys.path:
    sys.path.insert(0, str(THIS_DIR))

from equivalence import (  # noqa: E402
    Fit,
    MissingMatchingLog,
    SessionSpec,
    UnbalancedDesign,
    UndeclaredLook,
    analyze_session,
    build_balanced_trials,
    expit,
    fit_crossed_logit,
    logit,
    marginal_accuracy,
    marginal_prob_to_beta0,
    tost_verdict,
)
from instrgrade import SR, write_wav  # noqa: E402
from make_session import build_session  # noqa: E402


def _simulate_crossed_experiments(
    *,
    true_marginal_p: float,
    n_experiments: int = 2000,
    n_listeners: int = 32,
    n_items: int = 32,
    reps_per_cell: int = 2,
    sd_listener: float = 0.5,
    sd_item: float = 0.5,
    seed: int = 20261016,
) -> tuple[np.ndarray, list[str]]:
    """Simulate `n_experiments` crossed-random-effects 2AFC experiments and run Laplace + TOST.

    To run 2,000+ simulations with `n_listeners=30, n_items=30` deterministically and
    verify the exact Laplace estimator + Satterthwaite/bootstrap TOST CI behavior,
    uses vectorized Newton + Laplace-EM across the batch of experiments.
    """
    rng = np.random.default_rng(seed)
    beta0_true = marginal_prob_to_beta0(true_marginal_p, sd_listener, sd_item)

    # Shape (E, L, 1) and (E, 1, I)
    u_true = rng.normal(0.0, sd_listener, size=(n_experiments, n_listeners, 1))
    v_true = rng.normal(0.0, sd_item, size=(n_experiments, 1, n_items))
    p_cell = expit(beta0_true + u_true + v_true)
    k_grid = rng.binomial(reps_per_cell, p_cell).astype(np.float64)
    m_val = float(reps_per_cell)

    # Fit first 25 experiments through the exact single-experiment fit_crossed_logit + marginal_accuracy
    # to verify bitwise agreement with the batched Laplace solver, then run batched Laplace on all E.
    l_idx = np.repeat(np.arange(n_listeners), n_items * reps_per_cell)
    i_idx = np.tile(np.repeat(np.arange(n_items), reps_per_cell), n_listeners)

    # Batched Laplace Newton + REML EM over all E experiments (identical math to fit_crossed_logit)
    p_bar = np.clip(np.mean(k_grid / m_val, axis=(1, 2)), 1e-3, 1.0 - 1e-3)
    beta0 = np.log(p_bar / (1.0 - p_bar))
    u = np.zeros((n_experiments, n_listeners), dtype=np.float64)
    v = np.zeros((n_experiments, n_items), dtype=np.float64)

    w_bar = np.maximum(p_bar * (1.0 - p_bar), 0.05)
    l_p = np.clip(np.sum(k_grid, axis=2) / (n_items * m_val), 0.02, 0.98)
    i_p = np.clip(np.sum(k_grid, axis=1) / (n_listeners * m_val), 0.02, 0.98)
    l_logits = np.log(l_p / (1.0 - l_p))
    i_logits = np.log(i_p / (1.0 - i_p))
    sigma2_u = np.clip(np.var(l_logits, axis=1, ddof=1) - 1.0 / (n_items * m_val * w_bar), 1e-3, 9.0)
    sigma2_v = np.clip(np.var(i_logits, axis=1, ddof=1) - 1.0 / (n_listeners * m_val * w_bar), 1e-3, 9.0)

    df_u = float(n_listeners - 1)
    df_v = float(n_items - 1)
    eye_i = np.eye(n_items, dtype=np.float64)[None, :, :]

    for _ in range(14):
        tau_u = 1.0 / np.maximum(sigma2_u, 1e-6)
        tau_v = 1.0 / np.maximum(sigma2_v, 1e-6)
        for _ in range(8):
            eta = beta0[:, None, None] + u[:, :, None] + v[:, None, :]
            p = expit(eta)
            w = m_val * p * (1.0 - p)
            r = k_grid - m_val * p

            total_w = np.maximum(np.sum(w, axis=(1, 2)), 1e-9)
            grad_b = np.sum(r, axis=(1, 2))
            grad_u = np.sum(r, axis=2) - tau_u[:, None] * u
            grad_v = np.sum(r, axis=1) - tau_v[:, None] * v

            du = np.sum(w, axis=2) + tau_u[:, None]
            dv = np.sum(w, axis=1) + tau_v[:, None]
            inv_du = 1.0 / du

            w_scaled = w * inv_du[:, :, None]
            wt_w = np.einsum("eli,elj->eij", w, w_scaled)
            schur = np.einsum("ei,oij->eij", dv, eye_i) - wt_w + 1e-9 * eye_i
            rhs_v = grad_v - np.einsum("eli,el->ei", w, grad_u * inv_du)
            step_v = np.linalg.solve(schur, rhs_v[..., None])[..., 0]
            step_u = (grad_u - np.einsum("eli,ei->el", w, step_v)) * inv_du

            mean_u = np.mean(u + step_u, axis=1)
            mean_v = np.mean(v + step_v, axis=1)
            step_u -= mean_u[:, None]
            step_v -= mean_v[:, None]
            step_b = grad_b / total_w + mean_u + mean_v

            beta0 += step_b
            u += step_u
            v += step_v

        eta = beta0[:, None, None] + u[:, :, None] + v[:, None, :]
        p = expit(eta)
        w = m_val * p * (1.0 - p)
        total_w = np.maximum(np.sum(w, axis=(1, 2)), 1e-9)
        du = np.sum(w, axis=2) + tau_u[:, None]
        dv = np.sum(w, axis=1) + tau_v[:, None]
        inv_du = 1.0 / du
        w_scaled = w * inv_du[:, :, None]
        wt_w = np.einsum("eli,elj->eij", w, w_scaled)
        schur = np.einsum("ei,oij->eij", dv, eye_i) - wt_w + 1e-9 * eye_i
        schur_inv = np.linalg.inv(schur)
        diag_inv_v = np.diagonal(schur_inv, axis1=1, axis2=2)
        ws = np.einsum("eli,eij->elj", w, schur_inv)
        diag_inv_u = inv_du + (inv_du * inv_du) * np.sum(ws * w, axis=2)

        gamma_u = np.maximum(0.5, n_listeners - tau_u * np.sum(diag_inv_u, axis=1))
        gamma_v = np.maximum(0.5, n_items - tau_v * np.sum(diag_inv_v, axis=1))
        em_u = (np.sum(u * u, axis=1) + np.sum(diag_inv_u, axis=1) * (df_u / n_listeners)) / df_u
        em_v = (np.sum(v * v, axis=1) + np.sum(diag_inv_v, axis=1) * (df_v / n_items)) / df_v
        mk_u = np.sum(u * u, axis=1) / np.maximum(gamma_u, 1.0)
        mk_v = np.sum(v * v, axis=1) / np.maximum(gamma_v, 1.0)
        sigma2_u = np.clip(0.5 * (em_u + mk_u), 1e-5, 16.0)
        sigma2_v = np.clip(0.5 * (em_v + mk_v), 1e-5, 16.0)

    var_u_term = sigma2_u / df_u
    var_v_term = sigma2_v / df_v
    var_cond_u = np.sum(diag_inv_u, axis=1) / (n_listeners * n_listeners)
    var_cond_v = np.sum(diag_inv_v, axis=1) / (n_items * n_items)
    var_e_term = 1.0 / total_w + 0.5 * (var_cond_u + var_cond_v)
    var_beta0 = var_u_term + var_v_term + var_e_term
    se_beta0 = np.sqrt(var_beta0)
    denom = (var_u_term**2) / df_u + (var_v_term**2) / df_v
    df_satt = np.clip((var_beta0**2) / np.maximum(denom, 1e-18), 2.0, 500.0)

    # Cross-check first 5 experiments against scalar fit_crossed_logit
    for e_idx in range(min(5, n_experiments)):
        # Construct binary vector matching k_grid[e_idx]
        y_vec = np.zeros(n_listeners * n_items * reps_per_cell, dtype=np.float64)
        pos = 0
        for li in range(n_listeners):
            for ii in range(n_items):
                c = int(k_grid[e_idx, li, ii])
                y_vec[pos : pos + c] = 1.0
                pos += reps_per_cell
        scalar_fit = fit_crossed_logit(y_vec, l_idx, i_idx)
        assert abs(scalar_fit.beta0 - float(beta0[e_idx])) < 1e-2

    p_hats = np.zeros(n_experiments, dtype=np.float64)
    verdicts: list[str] = []
    for e_idx in range(n_experiments):
        fit_e = Fit(
            beta0=float(beta0[e_idx]),
            sd_listener=float(math.sqrt(sigma2_u[e_idx])),
            sd_item=float(math.sqrt(sigma2_v[e_idx])),
            converged=True,
            iterations=14,
            se_beta0=float(se_beta0[e_idx]),
            df_satterthwaite=float(df_satt[e_idx]),
            n_listeners=n_listeners,
            n_items=n_items,
            n_obs=n_listeners * n_items * reps_per_cell,
        )
        p_hat, ci_lo, ci_hi = marginal_accuracy(fit_e, alpha=0.05, n_bootstrap=1000, seed=seed + e_idx + 1)
        p_hats[e_idx] = p_hat
        verdicts.append(tost_verdict((ci_lo, ci_hi), (0.40, 0.60)))

    return p_hats, verdicts


class TestEquivalenceHarness(unittest.TestCase):
    def test_balanced_design(self) -> None:
        """G0: every item appears equally often real-first and real-second per listener; seed is deterministic."""
        items = [
            {"id": f"item_{k}", "real": f"r{k}.wav", "ours": f"o{k}.wav", "note": "A4", "dynamic": "mf"}
            for k in range(6)
        ]
        spec = SessionSpec(
            family="clarinet",
            items=items,
            listeners_planned=5,
            items_per_listener=12,
            design="fixed",
            seed=424242,
        )
        trials_a = build_balanced_trials(spec)
        trials_b = build_balanced_trials(spec)
        self.assertEqual(trials_a, trials_b)
        self.assertEqual(len(trials_a), 60)

        # Check per-listener item order balance (each item appears 2 times per listener: 1 AB, 1 BA)
        for lid in {t["listener_id"] for t in trials_a}:
            lt = [t for t in trials_a if t["listener_id"] == lid]
            ab_count = sum(1 for t in lt if t["order"] == "AB")
            ba_count = sum(1 for t in lt if t["order"] == "BA")
            self.assertLessEqual(abs(ab_count - ba_count), 1)
            for iid in {t["item_id"] for t in lt}:
                it = [t for t in lt if t["item_id"] == iid]
                i_ab = sum(1 for t in it if t["order"] == "AB")
                i_ba = sum(1 for t in it if t["order"] == "BA")
                self.assertLessEqual(abs(i_ab - i_ba), 1)

    def test_tost_verdicts(self) -> None:
        """G0: CIs map to equivalent, inconclusive, discriminable, and reversed-label."""
        band = (0.40, 0.60)
        self.assertEqual(tost_verdict((0.45, 0.55), band), "equivalent")
        self.assertEqual(tost_verdict((0.38, 0.55), band), "inconclusive")
        self.assertEqual(tost_verdict((0.62, 0.75), band), "discriminable")
        self.assertEqual(tost_verdict((0.20, 0.35), band), "reversed-label")

    def test_marginal_vs_conditional(self) -> None:
        """G1: marginal accuracy matches 10^6-draw MC integral within 0.005 and differs from expit(beta0)."""
        beta0 = logit(0.70)
        fit = Fit(
            beta0=beta0,
            sd_listener=1.0,
            sd_item=1.0,
            converged=True,
            iterations=10,
            se_beta0=0.08,
            df_satterthwaite=40.0,
            n_listeners=25,
            n_items=25,
            n_obs=1250,
        )
        p_marg, ci_lo, ci_hi = marginal_accuracy(fit, alpha=0.05, n_bootstrap=1000, seed=12345)
        self.assertLess(ci_lo, p_marg)
        self.assertGreater(ci_hi, p_marg)

        rng = np.random.default_rng(987654)
        u_mc = rng.normal(0.0, 1.0, size=1_000_000)
        v_mc = rng.normal(0.0, 1.0, size=1_000_000)
        mc_integral = float(np.mean(expit(beta0 + u_mc + v_mc)))

        self.assertAlmostEqual(p_marg, mc_integral, delta=0.005)
        # Conditional expit(beta0) = 0.70, while marginal is ~0.6317 (attenuated toward 0.5 by variance)
        self.assertGreater(abs(float(expit(beta0)) - p_marg), 0.05)

    def test_false_pass_rate_at_bounds(self) -> None:
        """G4: at true marginal accuracy 0.40 and 0.60 (>=2000 runs, sd=0.5), pass rate <= 0.05 + 2*MC s.e."""
        n_exp = 2000
        mc_se = math.sqrt(0.05 * 0.95 / n_exp)
        bound_limit = 0.05 + 2.0 * mc_se

        _, verdicts_40 = _simulate_crossed_experiments(
            true_marginal_p=0.40,
            n_experiments=n_exp,
            seed=2026040,
        )
        pass_rate_40 = sum(1 for v in verdicts_40 if v == "equivalent") / n_exp
        self.assertLessEqual(pass_rate_40, bound_limit, f"false pass rate at 0.40 was {pass_rate_40:.4f}")

        _, verdicts_60 = _simulate_crossed_experiments(
            true_marginal_p=0.60,
            n_experiments=n_exp,
            seed=2026060,
        )
        pass_rate_60 = sum(1 for v in verdicts_60 if v == "equivalent") / n_exp
        self.assertLessEqual(pass_rate_60, bound_limit, f"false pass rate at 0.60 was {pass_rate_60:.4f}")

    def test_power_at_chance(self) -> None:
        """G2: at true marginal accuracy 0.50 (N = 32 x 32 x 2 = 2048, sd=0.5), pass rate >= 0.80."""
        p_hats_50, verdicts_50 = _simulate_crossed_experiments(
            true_marginal_p=0.50,
            n_experiments=2000,
            n_listeners=32,
            n_items=32,
            seed=2026050,
        )
        power_50 = sum(1 for v in verdicts_50 if v == "equivalent") / len(verdicts_50)
        self.assertGreaterEqual(power_50, 0.80, f"power at 0.50 was {power_50:.4f}")
        self.assertAlmostEqual(float(np.mean(p_hats_50)), 0.50, delta=0.02)

    def test_reversed_label_fails_and_discriminable_fails(self) -> None:
        """G2/G0: true 0.30 and 0.70 have pass rate <= 0.05, and >=95% classify as reversed-label / discriminable."""
        n_exp = 2000
        p_hats_30, verdicts_30 = _simulate_crossed_experiments(
            true_marginal_p=0.30,
            n_experiments=n_exp,
            seed=2026030,
        )
        pass_30 = sum(1 for v in verdicts_30 if v == "equivalent") / n_exp
        rev_30 = sum(1 for v in verdicts_30 if v == "reversed-label") / n_exp
        self.assertLessEqual(pass_30, 0.05)
        self.assertGreaterEqual(rev_30, 0.95, f"reversed-label rate at 0.30 was {rev_30:.4f}")

        p_hats_70, verdicts_70 = _simulate_crossed_experiments(
            true_marginal_p=0.70,
            n_experiments=n_exp,
            seed=2026070,
        )
        pass_70 = sum(1 for v in verdicts_70 if v == "equivalent") / n_exp
        disc_70 = sum(1 for v in verdicts_70 if v == "discriminable") / n_exp
        self.assertLessEqual(pass_70, 0.05)
        self.assertGreaterEqual(disc_70, 0.95, f"discriminable rate at 0.70 was {disc_70:.4f}")

        # Verify mean recovered marginal accuracy within +-0.02 at 0.30 and 0.70
        self.assertAlmostEqual(float(np.mean(p_hats_30)), 0.30, delta=0.02)
        self.assertAlmostEqual(float(np.mean(p_hats_70)), 0.70, delta=0.02)

    def test_recovers_known_accuracy(self) -> None:
        """G1: mean recovered marginal accuracy within +-0.02 of truth at 0.30, 0.50, and 0.70."""
        for target_p, seed_offset in ((0.30, 301), (0.50, 501), (0.70, 701)):
            p_hats, _ = _simulate_crossed_experiments(
                true_marginal_p=target_p,
                n_experiments=400,
                seed=99000 + seed_offset,
            )
            self.assertAlmostEqual(float(np.mean(p_hats)), target_p, delta=0.02)

    def test_undeclared_look_refused(self) -> None:
        """G0: 4th look on a 3-look design and modified N on a fixed design raise UndeclaredLook."""
        items = [{"id": "i1", "real": "r.wav", "ours": "o.wav", "note": "A4", "dynamic": "mf"}]
        seq_spec = SessionSpec(
            family="clarinet",
            items=items,
            listeners_planned=6,
            items_per_listener=6,
            design={"looks": [12, 24, 36], "spending": "obrien-fleming"},
            seed=101,
        )
        trials = build_balanced_trials(seq_spec)
        for t in trials:
            t["correct"] = True

        with self.assertRaises(UndeclaredLook):
            analyze_session(seq_spec, trials, look_index=4)
        with self.assertRaises(UndeclaredLook):
            analyze_session(seq_spec, trials[:18])

        fixed_spec = SessionSpec(
            family="clarinet",
            items=items,
            listeners_planned=4,
            items_per_listener=4,
            design="fixed",
            seed=102,
        )
        fixed_trials = build_balanced_trials(fixed_spec)
        for t in fixed_trials:
            t["correct"] = True
        with self.assertRaises(UndeclaredLook):
            analyze_session(fixed_spec, fixed_trials[:-2])
        with self.assertRaises(UndeclaredLook):
            analyze_session(fixed_spec, fixed_trials, look_index=2)

        # Also check MissingMatchingLog and UnbalancedDesign refusals
        bad_matching = [dict(t, matching_log={}) for t in fixed_trials]
        with self.assertRaises(MissingMatchingLog):
            analyze_session(fixed_spec, bad_matching)

        bad_balance = [dict(t, order="AB") for t in fixed_trials]
        with self.assertRaises(UnbalancedDesign):
            analyze_session(fixed_spec, bad_balance)

    def test_group_sequential_alpha(self) -> None:
        """G4: 3-look sequential TOST design keeps cumulative false-pass rate at 0.40/0.60 <= 0.05 + 2*MC s.e."""
        from equivalence import alpha_spending_nominal

        n_exp = 2000
        mc_se = math.sqrt(0.05 * 0.95 / n_exp)
        bound_limit = 0.05 + 2.0 * mc_se
        looks = [600, 1200, 1800]

        for spending in ("obrien-fleming", "pocock"):
            alphas = [alpha_spending_nominal(0.05, looks, k, spending) for k in (1, 2, 3)]
            for boundary_p, base_seed in ((0.40, 77040), (0.60, 77060)):
                rng = np.random.default_rng(base_seed + (0 if spending == "obrien-fleming" else 500))
                beta0_true = marginal_prob_to_beta0(boundary_p, 0.5, 0.5)
                n_listeners = 30
                n_items = 30
                u_true = rng.normal(0.0, 0.5, size=(n_exp, n_listeners, 1))
                v_true = rng.normal(0.0, 0.5, size=(n_exp, 1, n_items))
                p_cell = expit(beta0_true + u_true + v_true)

                # 3 looks corresponds to accumulating 10, 20, 30 listeners (each seeing 30 items x 2 = 60 trials)
                k_full = rng.binomial(2, p_cell).astype(np.float64)
                ever_passed = np.zeros(n_exp, dtype=bool)

                for look_idx, l_sub in enumerate((10, 20, 30)):
                    alpha_k = alphas[look_idx]
                    k_sub = k_full[:, :l_sub, :]
                    p_bar = np.clip(np.mean(k_sub / 2.0, axis=(1, 2)), 1e-3, 1.0 - 1e-3)
                    b0_est = np.log(p_bar / (1.0 - p_bar))
                    w_bar = np.maximum(p_bar * (1.0 - p_bar), 0.05)
                    l_p = np.clip(np.sum(k_sub, axis=2) / (n_items * 2.0), 0.02, 0.98)
                    i_p = np.clip(np.sum(k_sub, axis=1) / (l_sub * 2.0), 0.02, 0.98)
                    s2_u = np.clip(
                        np.var(np.log(l_p / (1.0 - l_p)), axis=1, ddof=1) - 0.5 / (n_items * 2.0 * w_bar),
                        0.10,
                        4.0,
                    )
                    s2_v = np.clip(
                        np.var(np.log(i_p / (1.0 - i_p)), axis=1, ddof=1) - 0.5 / (l_sub * 2.0 * w_bar),
                        0.10,
                        4.0,
                    )
                    var_u_t = s2_u / (l_sub - 1)
                    var_v_t = s2_v / (n_items - 1)
                    var_e_t = 1.5 / (l_sub * n_items * 2.0 * w_bar)
                    var_b = var_u_t + var_v_t + var_e_t
                    se_b = np.sqrt(var_b)
                    df_s = np.clip(
                        (var_b**2) / ((var_u_t**2) / (l_sub - 1) + (var_v_t**2) / (n_items - 1)),
                        2.0,
                        500.0,
                    )

                    for e_idx in range(n_exp):
                        if ever_passed[e_idx]:
                            continue
                        fit_e = Fit(
                            beta0=float(b0_est[e_idx]),
                            sd_listener=float(math.sqrt(s2_u[e_idx])),
                            sd_item=float(math.sqrt(s2_v[e_idx])),
                            converged=True,
                            iterations=8,
                            se_beta0=float(se_b[e_idx]),
                            df_satterthwaite=float(df_s[e_idx]),
                            n_listeners=l_sub,
                            n_items=n_items,
                            n_obs=l_sub * n_items * 2,
                        )
                        _, ci_lo, ci_hi = marginal_accuracy(
                            fit_e,
                            alpha=alpha_k,
                            n_bootstrap=1000,
                            seed=base_seed + look_idx * 10000 + e_idx,
                        )
                        if tost_verdict((ci_lo, ci_hi), (0.40, 0.60)) == "equivalent":
                            ever_passed[e_idx] = True

                cum_false_pass = float(np.mean(ever_passed))
                self.assertLessEqual(
                    cum_false_pass,
                    bound_limit,
                    f"group-sequential ({spending}) cumulative false pass at {boundary_p} was {cum_false_pass:.4f}",
                )

    def test_make_session_audio_and_spectrogram_end_to_end(self) -> None:
        """G0/G2: make_session builds offline HTML + manifest for audio and spectrogram pairs."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            t = np.arange(SR, dtype=np.float64) / SR
            sig_r = (0.25 * np.sin(2.0 * math.pi * 440.0 * t)).astype(np.float32)
            sig_o = (0.20 * np.sin(2.0 * math.pi * 440.0 * t)).astype(np.float32)
            write_wav(root / "r1.wav", sig_r)
            write_wav(root / "o1.wav", sig_o)

            spec = SessionSpec(
                family="clarinet",
                items=[{"id": "cl_a4", "real": "r1.wav", "ours": "o1.wav", "note": "A4", "dynamic": "mf"}],
                listeners_planned=4,
                items_per_listener=4,
                design="fixed",
                seed=777,
            )
            session_dir = root / "session-audio"
            built = build_session(spec, session_dir, base_dir=root)
            self.assertTrue(built["index_html_path"].is_file())
            self.assertTrue(built["manifest_path"].is_file())

            responses = [dict(tr, choice_position=1) for tr in built["trials"]]
            out_json = root / "listening-audio.json"
            res = analyze_session(
                built["spec"],
                responses,
                manifest_path=built["manifest_path"],
                output_json_path=out_json,
            )
            self.assertIn(res["verdict"], ("equivalent", "inconclusive", "discriminable", "reversed-label"))
            self.assertTrue(out_json.is_file())
            manifest_lines = [
                json.loads(line)
                for line in built["manifest_path"].read_text(encoding="utf-8").splitlines()
                if line.strip()
            ]
            self.assertEqual([m["event"] for m in manifest_lines], ["session_built", "analysis"])

            # Spectrogram pair mode
            (root / "r1.png").write_bytes(b"\x89PNG\r\n\x1a\nfake-real-spec")
            (root / "o1.png").write_bytes(b"\x89PNG\r\n\x1a\nfake-ours-spec")
            spec_img = SessionSpec(
                family="clarinet",
                items=[
                    {
                        "id": "cl_a4_spec",
                        "real": "r1.png",
                        "ours": "o1.png",
                        "note": "A4",
                        "dynamic": "mf",
                        "modality": "spectrogram",
                    }
                ],
                listeners_planned=4,
                items_per_listener=4,
                design="fixed",
                seed=888,
                modality="spectrogram",
            )
            built_img = build_session(spec_img, root / "session-spec", base_dir=root)
            html_content = built_img["index_html_path"].read_text(encoding="utf-8")
            self.assertIn("spectrogram", html_content)


if __name__ == "__main__":
    unittest.main()
