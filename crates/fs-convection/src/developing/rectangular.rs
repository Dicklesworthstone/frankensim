//! Arithmetic shared by the original rectangular mean and its finite sections.
//!
//! Each linear-in-Graetz segment `Nu(Gz) = m Gz + c` gives cumulative
//! transfer `x Nu(P/x) = m P + c x`. Its local transfer is therefore the
//! intercept `c`. Integrating positive segment lengths avoids subtracting
//! nearly equal cumulative values, even for adjacent floating-point endpoints.

const GRAETZ: [f64; 14] = [
    10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0, 200.0, 220.0,
];
const NUSSELT: [f64; 14] = [
    4.20, 4.79, 5.23, 5.61, 5.95, 6.27, 6.88, 7.42, 7.91, 8.37, 8.80, 9.20, 9.60, 10.00,
];

/// Preserve the original table interpolation and explicitly declared bridge.
pub(crate) fn mean_nusselt(graetz: f64, developed: f64) -> f64 {
    if graetz <= GRAETZ[0] {
        return developed + (NUSSELT[0] - developed) * graetz / GRAETZ[0];
    }
    for index in 1..GRAETZ.len() {
        if graetz <= GRAETZ[index] {
            let fraction = (graetz - GRAETZ[index - 1]) / (GRAETZ[index] - GRAETZ[index - 1]);
            return NUSSELT[index - 1] + fraction * (NUSSELT[index] - NUSSELT[index - 1]);
        }
    }
    NUSSELT[NUSSELT.len() - 1]
}

/// Endpoints have already passed the unchanged original source domain.
/// Exact-entrance sections use the original mean before reaching this kernel.
pub(crate) fn section_nusselt(peclet: f64, start: f64, end: f64, developed: f64) -> f64 {
    let width = end - start;
    // The Gz<10 engineering bridge has constant section transfer: its entire
    // entrance correction is already carried by the upstream cumulative mean.
    let bridge_start = start.max(peclet / GRAETZ[0]);
    let mut mean = developed * ((end - bridge_start).max(0.0) / width);
    for index in 1..GRAETZ.len() {
        let a = start.max(peclet / GRAETZ[index]);
        let b = end.min(peclet / GRAETZ[index - 1]);
        if b > a {
            let slope = (NUSSELT[index] - NUSSELT[index - 1]) / (GRAETZ[index] - GRAETZ[index - 1]);
            let intercept = NUSSELT[index - 1] - slope * GRAETZ[index - 1];
            mean += intercept * ((b - a) / width);
        }
    }
    mean
}

#[cfg(test)]
mod tests {
    use super::*;

    // Independently evaluated original Shah-London CWT polynomial at a=0.5.
    const DEVELOPED: f64 = 3.388736875;

    fn near(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() <= 3e-13 * expected.abs().max(1.0),
            "{actual} != {expected}"
        );
    }

    #[test]
    fn unequal_sections_preserve_cumulative_heat_across_table_and_bridge() {
        let peclet = 72.0;
        let positions = [0.4, 0.63, 0.9, 1.4, 2.6, 4.3, 7.8, 12.0];
        let start = positions[0];
        let end = positions[positions.len() - 1];
        let expected = end * mean_nusselt(peclet / end, DEVELOPED)
            - start * mean_nusselt(peclet / start, DEVELOPED);
        let sum: f64 = positions
            .windows(2)
            .map(|pair| section_nusselt(peclet, pair[0], pair[1], DEVELOPED) * (pair[1] - pair[0]))
            .sum();
        near(sum, expected);
        near(
            section_nusselt(peclet, start, end, DEVELOPED),
            expected / (end - start),
        );
        // Treating each section as a new entrance increases the heat spuriously.
        let reset: f64 = positions
            .windows(2)
            .map(|pair| {
                let length = pair[1] - pair[0];
                mean_nusselt(peclet / length, DEVELOPED) * length
            })
            .sum();
        assert!(reset > 1.2 * sum);
    }

    #[test]
    fn short_sections_retain_table_segment_limits_and_the_bridge_limit() {
        // Independent intercepts: (Gz, Nu)=(10,4.20),(20,4.79) -> 3.61;
        // (40,5.61),(50,5.95) -> 4.25; (180,9.20),(200,9.60) -> 5.60.
        for (start, expected) in [
            (5.0_f64, 3.61),
            (1.6, 4.25),
            (0.375, 5.6),
            (10.0, DEVELOPED),
        ] {
            near(
                section_nusselt(72.0, start, start.next_up(), DEVELOPED),
                expected,
            );
        }
        let a = 5.0_f64;
        let b = a.next_up();
        let naive = (b * mean_nusselt(72.0 / b, DEVELOPED) - a * mean_nusselt(72.0 / a, DEVELOPED))
            / (b - a);
        assert!(
            (naive - 3.61).abs() > 0.1,
            "fixture must expose cancellation"
        );
    }

    #[test]
    fn knot_crossing_uses_both_adjacent_intercepts_without_smoothing() {
        // Gz=20 is x=3.6. The preceding intercept is 3.91 and the following
        // one is 3.61; a symmetric spatial interval has their arithmetic mean.
        near(section_nusselt(72.0, 3.5, 3.7, DEVELOPED), 3.76);
        let knot = 3.6_f64;
        near(
            section_nusselt(72.0, knot.next_down(), knot, DEVELOPED),
            3.91,
        );
        near(section_nusselt(72.0, knot, knot.next_up(), DEVELOPED), 3.61);
        near(section_nusselt(72.0, 10.0, 100.0, DEVELOPED), DEVELOPED);
    }
}
