//! Optimizer load declarations must not shrink with the active PDE domain.
use super::*;
use fs_cutfem::CutSdf;

// CutFEM deliberately assembles traction only on active material boundaries.
// A fixed-load optimization problem has a stronger contract: the COMPLETE
// declared band must stay inside material, even when its objective weight is
// zero. Otherwise deleting a loaded component can look like perfect compliance.
// Kernel admission has already checked that every level-set node is finite.
pub(super) fn require(phi: &GridSdf, support: EdgeBand, case: usize) -> Result<(), CutFemError> {
    let (a, b) = (support.start(), support.end());
    let (lo, hi) = match support.edge() {
        DesignBoxEdge::Bottom => ([a, 0.0], [b, 0.0]),
        DesignBoxEdge::Top => ([a, 1.0], [b, 1.0]),
        DesignBoxEdge::Left => ([0.0, a], [0.0, b]),
        DesignBoxEdge::Right => ([1.0, a], [1.0, b]),
    };
    // The existing bilinear enclosure visits ALL crossed grid intervals, not
    // just the band's endpoints or a few samples. Uncertain/cut support refuses
    // under the same strict-negative convention as canonical traction assembly.
    let enclosure = phi.enclose(lo, hi);
    if !(enclosure.lo().is_finite() && enclosure.hi().is_finite() && enclosure.hi() < 0.0) {
        return Err(invalid(format!(
            "multi-load case {case} requires its complete {:?} traction band {a}..={b} inside material; absent, cut or uncertain load support cannot be discarded",
            support.edge(),
        )));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn full_band_is_checked_on_each_named_edge() {
        for edge in [DesignBoxEdge::Bottom, DesignBoxEdge::Top,
            DesignBoxEdge::Left, DesignBoxEdge::Right]
        {
            let support = EdgeBand::new(edge, 0.25, 0.75).unwrap();
            assert!(require(&GridSdf::from_fn(8, &|_, _| -1.0), support, 0).is_ok());
            let absent = GridSdf::from_fn(8, &|x, y| match edge {
                DesignBoxEdge::Bottom => 0.1 - y,
                DesignBoxEdge::Top => y - 0.9,
                DesignBoxEdge::Left => 0.1 - x,
                DesignBoxEdge::Right => x - 0.9,
            });
            assert!(matches!(require(&absent, support, 2),
                Err(CutFemError::InvalidElasticityInput { .. })));
            // A boundary exactly on phi=0 is not a certified loaded trace.
            assert!(require(&GridSdf::from_fn(8, &|_, _| 0.0), support, 0).is_err());
        }
    }

    #[test]
    fn missing_interior_support_cannot_hide_between_endpoint_and_midpoint_probes() {
        let support = EdgeBand::new(DesignBoxEdge::Right, 0.125, 0.875).unwrap();
        let mut field = GridSdf::from_fn(8, &|_, _| -1.0);
        field.nodes_mut()[2 * 9 + 8] = 1.0;
        for y in [0.125, 0.5, 0.875] { assert!(field.value_at([1.0, y]) < 0.0); }
        assert!(require(&field, support, 0).is_err());
        // Do not widen the policy beyond the actually declared load band.
        let unrelated = EdgeBand::new(DesignBoxEdge::Right, 0.5, 0.875).unwrap();
        assert!(require(&field, unrelated, 0).is_ok());
    }

    #[test]
    fn polled_late_zero_weight_load_refuses_before_a_partial_family_is_returned() {
        let field = GridSdf::from_fn(8, &|_, y| (y - 0.5).abs() - 0.35);
        let cases = [
            RobustLoadCase::new(DesignBoxEdge::Right, 0.375, 0.625, [0.0, -1.0], 1.0).unwrap(),
            RobustLoadCase::new(DesignBoxEdge::Top, 0.375, 0.625, [1.0, 0.0], 0.0).unwrap(),
        ];
        let settings = OptimizeSettings { level: 3, ..OptimizeSettings::default() };
        let kernel = Kernel::new(&field, &cases, settings, RobustAggregate::WeightedSum).unwrap();
        let mut boundaries = Vec::new();
        let result = kernel.evaluate_scheduled(field, Some(1), |case, stage| {
            match stage {
                CaseProgress::Start => boundaries.push((case, false)),
                CaseProgress::Complete => boundaries.push((case, true)),
                CaseProgress::Iterations(_) => {},
            }
            ControlFlow::<()>::Continue(())
        });
        assert!(matches!(result, Err(CutFemError::InvalidElasticityInput { ref what })
            if what.contains("case 1") && what.contains("load support")));
        assert_eq!(boundaries, [(0, false), (0, true), (1, false)]);
    }
}
