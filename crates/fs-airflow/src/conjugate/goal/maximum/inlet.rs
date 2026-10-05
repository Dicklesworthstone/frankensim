//! Pull back full-system reference-offset derivatives to prescribed air inlets.
//! The caller supplies the COMPLETE coupled adjoint, never frozen-solid bars.

use super::{AirPath, BTreeSet, Cx, Result, admitted_exchange_terms, bad, finite, poll, zeros};

/// Derivatives with respect to independent branch inlet temperatures, in path order.
///
/// `regions` and `reference_bars` must be the exact port order returned by the
/// full solid/air adjoint (`pullback_affine_controls`). For segment j the
/// reference offset changes by g_j times the product of upstream (1-eps).
/// The wall feedback is already in the supplied adjoint and MUST NOT be
/// applied a second time. Every branch starts with its own unit inlet change;
/// branches exchange heat through the solid, but their inlets are not mixed.
///
/// Geometry, h, mass flow and heat capacity remain fixed. This is an Estimated
/// derivative of the stored exponential exchange model, not a flow derivative,
/// uncertainty bound, or derivative of the location of a temperature maximum.
/// Uses O(ports) name storage and O(branches) output; no dense transfer matrix,
/// temperature division, perturbed primal, or iteration-trace derivative.
///
/// # Errors
/// Empty paths, exceeded port budget, mismatched/duplicate port ownership,
/// nonfinite bars or contractions, and cancellation return no partial result.
pub fn pullback_inlet_temperatures(
    cx: &Cx<'_>, paths: &[AirPath], regions: &[&str], reference_bars: &[f64],
    max_ports: usize,
) -> Result<Vec<f64>> {
    poll(cx)?;
    if paths.is_empty() || paths.len() > max_ports {
        return Err(bad("air-inlet pullback needs nonempty paths within the port budget"));
    }
    let mut n = 0_usize;
    for path in paths {
        poll(cx)?;
        n = n.checked_add(path.segments().len()).ok_or_else(|| bad("air-inlet port count overflow"))?;
        if n > max_ports { return Err(bad("air-inlet pullback exceeds the port budget")); }
    }
    if regions.len() != n || reference_bars.len() != n {
        return Err(bad("air-inlet pullback needs one named full-system derivative per segment"));
    }
    // Validate every coordinate before arithmetic; zip must never hide a
    // missing segment or reinterpret a derivative for another branch.
    let mut seen = BTreeSet::new();
    for ((segment, &name), &bar) in paths.iter().flat_map(|p| p.segments())
        .zip(regions).zip(reference_bars)
    {
        poll(cx)?;
        if segment.region() != name || !seen.insert(name) {
            return Err(bad("air-inlet derivative order or unique solid-trace ownership differs"));
        }
        finite(bar)?;
    }
    let mut result = zeros(cx, paths.len())?;
    let mut index = 0;
    for (branch, path) in paths.iter().enumerate() {
        let mut inlet_change = 1.0;
        for segment in path.segments() {
            poll(cx)?;
            let (_, eps, g) = admitted_exchange_terms(segment, path.capacity_rate_w_per_k())?;
            let weight = finite(g * inlet_change)?;
            result[branch] = finite(reference_bars[index].mul_add(weight, result[branch]))?;
            inlet_change = finite((1.0 - eps) * inlet_change)?;
            index += 1;
        }
    }
    poll(cx)?;
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::conjugate::AirSegment;
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};

    fn with_cx<T>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> T) -> T {
        ArenaPool::new(ArenaConfig::default()).scope(|arena| f(&Cx::new(gate, arena,
            StreamKey { seed: 7, kernel_id: 831, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic)))
    }
    fn path(inlet: f64, names: &[&str], ntus: &[f64]) -> AirPath {
        AirPath::new(inlet, 0.01, 1000.0, names.iter().zip(ntus)
            .map(|(name, ntu)| AirSegment::new(name, 0.1, ntu * 100.0).unwrap())
            .collect()).unwrap()
    }

    #[test]
    fn inlet_pullback_matches_multisegment_independent_branch_marches() {
        with_cx(&CancelGate::new_clock_free(), |cx| {
            let paths = [path(290.0, &["a0", "a1", "a2"], &[0.1, 0.8, 2.0]),
                path(315.0, &["b0", "b1"], &[1.5, 0.3])];
            let names = ["a0", "a1", "a2", "b0", "b1"];
            let bars = [0.7, -1.2, 2.3, -0.4, 1.1];
            let actual = pullback_inlet_temperatures(cx, &paths, &names, &bars, 5).unwrap();
            let evaluate = |a, b| {
                let a = path(a, &names[..3], &[0.1, 0.8, 2.0]).march(&[330.0, 320.0, 340.0]).unwrap();
                let b = path(b, &names[3..], &[1.5, 0.3]).march(&[305.0, 350.0]).unwrap();
                a.segments.iter().chain(&b.segments).zip(bars)
                    .map(|(s, bar)| s.reference_temperature_k * bar).sum::<f64>()
            };
            let h = 0.01;
            let expected = [(evaluate(290.0+h,315.0)-evaluate(290.0-h,315.0))/(2.0*h),
                (evaluate(290.0,315.0+h)-evaluate(290.0,315.0-h))/(2.0*h)];
            for (a, e) in actual.iter().zip(expected) { assert!((a-e).abs() < 1e-9, "{a} vs {e}"); }
            assert!((actual[0]-bars[0]).abs() > 0.01, "not just the first reference bar");
            let other = [paths[0].clone(), path(400.0, &names[3..], &[1.5, 0.3])];
            assert_eq!(actual, pullback_inlet_temperatures(cx, &other, &names, &bars, 5).unwrap());
            assert_eq!(actual, pullback_inlet_temperatures(cx, &paths, &names, &bars, 5).unwrap());
        });
    }

    #[test]
    fn inlet_pullback_keeps_tiny_ntu_and_zero_downstream_carry_defined() {
        with_cx(&CancelGate::new_clock_free(), |cx| {
            let paths = [path(300.0, &["tiny"], &[1e-12]),
                path(300.0, &["saturated", "downstream"], &[1000.0, 0.2])];
            let bars = [1.0, 1.0, 1e20];
            let actual = pullback_inlet_temperatures(cx, &paths,
                &["tiny", "saturated", "downstream"], &bars, 3).unwrap();
            assert!((actual[0]-1.0).abs() < 1e-11);
            assert!((actual[1]-0.001).abs() < 1e-15);
        });
    }

    #[test]
    fn inlet_pullback_refuses_wrong_ownership_shapes_work_and_cancellation() {
        let paths = [path(300.0, &["a", "b"], &[0.2, 0.3])];
        with_cx(&CancelGate::new_clock_free(), |cx| {
            for (names, bars, budget) in [(&["b", "a"][..], &[1.0, 1.0][..], 2),
                (&["a"][..], &[1.0, 1.0][..], 2),
                (&["a", "b"][..], &[1.0][..], 2),
                (&["a", "b"][..], &[1.0, f64::NAN][..], 2),
                (&["a", "b"][..], &[1.0, 1.0][..], 1)] {
                assert!(pullback_inlet_temperatures(cx, &paths, names, bars, budget).is_err());
            }
            assert!(pullback_inlet_temperatures(cx, &[], &[], &[], 2).is_err());
            let duplicate = [paths[0].clone(), paths[0].clone()];
            assert!(pullback_inlet_temperatures(cx, &duplicate, &["a", "b", "a", "b"], &[1.0;4], 4).is_err());
            assert_eq!(pullback_inlet_temperatures(cx, &paths, &["a", "b"], &[0.0;2], 2).unwrap(), [0.0]);
        });
        let gate = CancelGate::new_clock_free(); gate.request();
        with_cx(&gate, |cx| assert!(matches!(pullback_inlet_temperatures(cx, &paths,
            &["a", "b"], &[1.0;2], 2), Err(super::super::super::CoupledGoalError::Interrupted))));
    }
}
