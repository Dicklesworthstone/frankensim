//! Allocation-free arithmetic for the admitted bare-tube bank owner.
//! Inputs are already checked by `tube_bank`; no extrapolation policy lives here.
use fs_math::det;

/// Arrangement is explicit: equal longitudinal/transverse pitches do not imply inline tubes.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TubeBankLayout {
    Inline,
    Staggered,
}

pub(crate) fn gap_velocity_ratio(
    layout: TubeBankLayout,
    transverse: f64,
    longitudinal: f64,
) -> f64 {
    let transverse_gap = transverse - 1.0;
    // Two diagonal passages share each transverse-pitch inlet streamtube.
    let gap = if layout == TubeBankLayout::Staggered {
        let diagonal = det::sqrt(longitudinal * longitudinal + 0.25 * transverse * transverse);
        transverse_gap.min(2.0 * (diagonal - 1.0))
    } else {
        transverse_gap
    };
    transverse / gap
}

// Original table with declared linear interpolation between retained knots.
// At 16 or more rows the source uses the uncorrected mean.
fn row_correction(layout: TubeBankLayout, rows: u32) -> f64 {
    const ROWS: [u32; 9] = [1, 2, 3, 4, 5, 7, 10, 13, 16];
    let factors = match layout {
        TubeBankLayout::Inline => [0.70, 0.80, 0.86, 0.90, 0.93, 0.96, 0.98, 0.99, 1.0],
        TubeBankLayout::Staggered => [0.64, 0.76, 0.84, 0.89, 0.93, 0.96, 0.98, 0.99, 1.0],
    };
    for i in 1..ROWS.len() {
        if rows <= ROWS[i] {
            let fraction = f64::from(rows - ROWS[i - 1]) / f64::from(ROWS[i] - ROWS[i - 1]);
            return factors[i - 1] + fraction * (factors[i] - factors[i - 1]);
        }
    }
    1.0
}

pub(crate) fn mean_nusselt(
    layout: TubeBankLayout,
    rows: u32,
    reynolds: f64,
    prandtl: f64,
    surface_prandtl: f64,
    transverse_over_longitudinal: f64,
) -> f64 {
    let (coefficient, exponent) = match layout {
        TubeBankLayout::Inline => (0.27, 0.63),
        TubeBankLayout::Staggered => (0.35 * det::pow(transverse_over_longitudinal, 0.2), 0.6),
    };
    row_correction(layout, rows)
        * coefficient
        * det::pow(reynolds, exponent)
        * det::pow(prandtl, 0.36)
        * det::pow(prandtl / surface_prandtl, 0.25)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn near(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() < 2e-12 * expected.abs().max(1.0),
            "{actual} != {expected}"
        );
    }
    #[test]
    fn g1_mean_matches_independent_sixty_digit_source_arithmetic() {
        near(
            mean_nusselt(TubeBankLayout::Inline, 4, 20_000.0, 7.0, 3.5, 1.0),
            298.36506027896674,
        );
        near(
            mean_nusselt(TubeBankLayout::Staggered, 4, 20_000.0, 7.0, 3.5, 1.0),
            284.16385885169537,
        );
        for layout in [TubeBankLayout::Inline, TubeBankLayout::Staggered] {
            let base = mean_nusselt(layout, 4, 20_000.0, 7.0, 7.0, 1.0);
            near(
                mean_nusselt(layout, 4, 20_000.0, 7.0, 3.5, 1.0) / base,
                2_f64.sqrt().sqrt(),
            );
        }
    }
    #[test]
    fn g0_gap_velocity_follows_transverse_and_diagonal_mass_continuity() {
        near(gap_velocity_ratio(TubeBankLayout::Inline, 2.0, 2.0), 2.0);
        near(gap_velocity_ratio(TubeBankLayout::Staggered, 2.0, 2.0), 2.0);
        near(
            gap_velocity_ratio(TubeBankLayout::Staggered, 2.0, 1.0),
            1.0 / (2_f64.sqrt() - 1.0),
        );
    }
    #[test]
    fn g3_rows_and_flow_changes_keep_the_original_source_scaling() {
        for layout in [TubeBankLayout::Inline, TubeBankLayout::Staggered] {
            let full = mean_nusselt(layout, 16, 10_000.0, 0.71, 0.71, 1.0);
            for (rows, factor) in [
                (5, 0.93),
                (6, 0.945),
                (7, 0.96),
                (10, 0.98),
                (13, 0.99),
                (16, 1.0),
                (20, 1.0),
            ] {
                near(
                    mean_nusselt(layout, rows, 10_000.0, 0.71, 0.71, 1.0) / full,
                    factor,
                );
            }
            let exponent = if layout == TubeBankLayout::Inline {
                0.63
            } else {
                0.6
            };
            near(
                mean_nusselt(layout, 16, 20_000.0, 0.71, 0.71, 1.0) / full,
                2_f64.powf(exponent),
            );
        }
    }
}
