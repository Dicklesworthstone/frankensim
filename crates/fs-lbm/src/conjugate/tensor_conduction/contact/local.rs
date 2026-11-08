//! Conservative two-trace heat and the exact resistance chain rule.

/// Heat is positive from the first trace to the second. The second equation
/// receives the exact opposite row; no contact storage is introduced.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct ContactHeat {
    pub heat_w: f64,
    pub derivative_w_k: [f64; 2],
    pub arithmetic_w: f64,
}

impl ContactHeat {
    pub(crate) fn reversed(self) -> Self {
        Self {
            heat_w: -self.heat_w,
            derivative_w_k: [-self.derivative_w_k[1], -self.derivative_w_k[0]],
            arithmetic_w: self.arithmetic_w,
        }
    }
}

pub(super) fn linear(conductance_w_k: f64, temperatures_k: [f64; 2]) -> Option<ContactHeat> {
    if !conductance_w_k.is_finite()
        || conductance_w_k <= 0.0
        || temperatures_k.iter().any(|t| !t.is_finite() || *t <= 0.0)
    {
        return None;
    }
    let [a, b] = temperatures_k;
    let heat_w = conductance_w_k * (a - b);
    let arithmetic_w = conductance_w_k * (a.abs() + b.abs());
    (heat_w.is_finite() && arithmetic_w.is_finite()).then_some(ContactHeat {
        heat_w,
        derivative_w_k: [conductance_w_k, -conductance_w_k],
        arithmetic_w,
    })
}

/// q=A*(Ta-Tb)/R and dq/dTi=(+/- A/R)-q*(dR/dTi)/R.
/// R is area-specific; its partials are with respect to the TWO actual traces.
pub(super) fn original_resistance(
    area_m2: f64,
    temperatures_k: [f64; 2],
    resistance_m2_k_w: f64,
    resistance_partials_m2_w: [f64; 2],
) -> Option<ContactHeat> {
    if !area_m2.is_finite()
        || area_m2 <= 0.0
        || !resistance_m2_k_w.is_finite()
        || resistance_m2_k_w <= 0.0
        || resistance_partials_m2_w.iter().any(|d| !d.is_finite())
    {
        return None;
    }
    let mut result = linear(area_m2 / resistance_m2_k_w, temperatures_k)?;
    for (derivative, partial) in result
        .derivative_w_k
        .iter_mut()
        .zip(resistance_partials_m2_w)
    {
        // At equilibrium the feedback is exactly zero, even if a separately
        // formed logarithmic resistance slope would be unrepresentable.
        if result.heat_w != 0.0 && partial != 0.0 {
            *derivative = (-result.heat_w).mul_add(partial / resistance_m2_k_w, *derivative);
        }
        if !derivative.is_finite() {
            return None;
        }
    }
    Some(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn resistance([a, b]: [f64; 2]) -> f64 {
        0.5 + (a - 320.0) / 1024.0 + 3.0 * (b - 320.0) / 1024.0
    }

    #[test]
    fn original_two_trace_tangent_matches_independent_quotient_difference() {
        let area = 0.06;
        let t = [350.0, 320.0];
        let heat =
            original_resistance(area, t, resistance(t), [1.0 / 1024.0, 3.0 / 1024.0]).unwrap();
        assert!((heat.heat_w - area * (t[0] - t[1]) / resistance(t)).abs() < 1e-14);
        for side in 0..2 {
            let mut a = t;
            let mut b = t;
            let epsilon = 1e-3;
            a[side] += epsilon;
            b[side] -= epsilon;
            let q = |s: [f64; 2]| area * (s[0] - s[1]) / resistance(s);
            let numerical = (q(a) - q(b)) / (2.0 * epsilon);
            assert!((heat.derivative_w_k[side] - numerical).abs() < 2e-10);
        }
        assert_ne!(heat.derivative_w_k[0], -heat.derivative_w_k[1]);
    }

    #[test]
    fn source_relabeling_preserves_heat_and_the_complete_tangent() {
        let t = [350.0, 320.0];
        let forward =
            original_resistance(0.06, t, resistance(t), [1.0 / 1024.0, 3.0 / 1024.0]).unwrap();
        let reverse = original_resistance(
            0.06,
            [t[1], t[0]],
            resistance(t),
            [3.0 / 1024.0, 1.0 / 1024.0],
        )
        .unwrap();
        assert_eq!(reverse, forward.reversed());
        let direction = [0.75, -0.25];
        let dq =
            forward.derivative_w_k[0] * direction[0] + forward.derivative_w_k[1] * direction[1];
        assert_eq!(dq + -dq, 0.0);
        assert!(forward.heat_w * (t[1].recip() - t[0].recip()) > 0.0);
    }

    #[test]
    fn a_constant_source_retains_the_exact_linear_path() {
        let t = [350.0, 320.0];
        assert_eq!(
            original_resistance(0.125, t, 0.5, [0.0; 2]),
            linear(0.25, t)
        );
        let source = original_resistance(1e-200, [320.0; 2], 1e-200, [1e200; 2]).unwrap();
        assert_eq!(source.heat_w, 0.0);
        assert_eq!(source.derivative_w_k, [1.0, -1.0]);
    }

    #[test]
    fn unresolved_or_nonphysical_interface_data_refuse() {
        assert!(original_resistance(0.0, [300.0, 400.0], 1.0, [0.0; 2]).is_none());
        assert!(original_resistance(1.0, [0.0, 400.0], 1.0, [0.0; 2]).is_none());
        assert!(original_resistance(1.0, [300.0, 400.0], 0.0, [0.0; 2]).is_none());
        assert!(original_resistance(1.0, [300.0, 400.0], 1.0, [f64::NAN, 0.0]).is_none());
        assert!(
            original_resistance(f64::MIN_POSITIVE, [300.0, 400.0], f64::MAX, [0.0; 2]).is_none()
        );
    }
}
