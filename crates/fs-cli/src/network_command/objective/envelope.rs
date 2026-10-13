//! Design-only smooth upper bound for an explicitly selected discrete maximum.
//! The unnormalized log-sum-exp is essential: subtracting ln(n) would produce
//! an average-like surrogate below the peak and lose conservative feasibility.

use super::*;
use fs_math::det;

#[derive(Debug)]
pub(super) struct Envelope {
    pub upper_k: f64,
    pub maximum_k: f64,
    pub smoothing_k: f64,
    pub bias_bound_k: f64,
    pub weights: Vec<(usize, f64)>,
}

impl Envelope {
    pub fn new(cx: &Cx<'_>, temperatures: &[f64], vertices: &[usize], maximum: f64, tau: f64) -> Result<Self> {
        poll(cx)?;
        if vertices.is_empty() || !maximum.is_finite() || !tau.is_finite() || tau <= 0.0 {
            return Err(bad("upper envelope requires a finite maximum, positive smoothing and nonempty selection"));
        }
        let bias_bound_k = tau * det::ln(vertices.len() as f64);
        if !bias_bound_k.is_finite() { return Err(producer("upper-envelope bias bound overflow")); }
        let mut weights = Vec::with_capacity(vertices.len());
        let mut sum = 0.0_f64;
        let mut correction = 0.0_f64;
        let mut seen_maximum = false;
        for &vertex in vertices {
            poll(cx)?;
            let value = *temperatures.get(vertex).ok_or_else(|| bad("upper-envelope vertex outside the solved field"))?;
            if !value.is_finite() || value > maximum { return Err(producer("invalid upper-envelope maximum")); }
            // A very cold tail can legitimately underflow to zero. In contrast,
            // every maximizing entry contributes exactly one, even for tiny tau.
            let term = if value == maximum { seen_maximum = true; 1.0 } else {
                let shifted = (value - maximum) / tau;
                if shifted == f64::NEG_INFINITY { 0.0 } else { det::exp(shifted) }
            };
            if !term.is_finite() || term < 0.0 { return Err(producer("nonfinite upper-envelope weight")); }
            let next = sum + term;
            correction += if sum >= term { (sum - next) + term } else { (term - next) + sum };
            sum = next;
            weights.push((vertex, term));
        }
        sum += correction;
        if !seen_maximum || !sum.is_finite() || sum < 1.0 { return Err(producer("invalid upper-envelope normalization")); }
        let upper_k = maximum + tau * det::ln(sum);
        if !upper_k.is_finite() || upper_k < maximum { return Err(producer("unrepresentable upper-envelope temperature")); }
        for (_, weight) in &mut weights { poll(cx)?; *weight /= sum; }
        poll(cx)?;
        Ok(Self { upper_k, maximum_k: maximum, smoothing_k: tau, bias_bound_k, weights })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::network_command::tests::{request, with_cx};

    fn close(a: f64, b: f64, tol: f64) {
        assert!((a-b).abs() <= tol * b.abs().max(1.0), "{a:.16e} != {b:.16e}");
    }
    fn selection(r: &Request, json: &str, tau: f64) -> Objective {
        let mut objective = Objective::parse(&J::parse(json).unwrap(), &r.surfaces, &r.mesh).unwrap();
        objective.enable_upper_envelope(tau).unwrap();
        objective
    }

    #[test]
    fn ties_use_all_selected_nodes_without_underestimating_the_peak() {
        let r = request();
        let objective = selection(&r, r#"{"max_vertices":[1,3]}"#, 0.2);
        let mut temperatures = vec![900.0; r.mesh.vertex_count()];
        temperatures[1] = 350.0; temperatures[3] = 350.0;
        with_cx(|cx| {
            let state = objective.evaluate(cx, &temperatures, &[]).unwrap();
            close(state.value, 350.0 + 0.2 * 2.0_f64.ln(), 1e-14);
            assert_eq!(state.physical_value(), 350.0);
            assert_eq!(state.vertex, Some(1)); assert_eq!(state.ties, 2);
            assert_eq!(state.envelope.as_ref().unwrap().weights, vec![(1, 0.5), (3, 0.5)]);
            let document = J::parse(&objective.render(&state, &r.mesh).unwrap()).unwrap();
            assert_eq!(document.str_field("kind"), Some("smooth_upper_envelope"));
            assert_eq!(document.f64_field("actual_maximum_k"), Some(350.0));
            assert!(state.value > state.physical_value());
        });
    }

    #[test]
    fn smooth_weights_match_each_temperature_partial_and_shift_identity() {
        let r = request();
        let objective = selection(&r, r#"{"max_vertices":[0,1,2]}"#, 0.1);
        let mut temperatures = vec![300.0; r.mesh.vertex_count()];
        temperatures[..3].copy_from_slice(&[315.0, 315.01, 314.97]);
        with_cx(|cx| {
            let state = objective.evaluate(cx, &temperatures, &[]).unwrap();
            let envelope = state.envelope.as_ref().unwrap();
            close(envelope.weights.iter().map(|(_, w)| w).sum(), 1.0, 1e-14);
            for &(vertex, weight) in &envelope.weights {
                let h = 1e-5;
                let mut plus = temperatures.clone(); plus[vertex] += h;
                let mut minus = temperatures.clone(); minus[vertex] -= h;
                let difference = (objective.evaluate(cx, &plus, &[]).unwrap().value
                    - objective.evaluate(cx, &minus, &[]).unwrap().value) / (2.0*h);
                close(weight, difference, 1e-8);
            }
            for value in &mut temperatures { *value += 7.0; }
            close(objective.evaluate(cx, &temperatures, &[]).unwrap().value, state.value + 7.0, 1e-14);
            assert!(envelope.upper_k >= envelope.maximum_k);
            assert!(envelope.upper_k - envelope.maximum_k <= envelope.bias_bound_k + 1e-12);
        });
    }

    #[test]
    fn singletons_and_widely_separated_temperatures_remain_representable() {
        with_cx(|cx| {
            let one = Envelope::new(cx, &[1e300], &[0], 1e300, 1e-300).unwrap();
            assert_eq!(one.upper_k, one.maximum_k); assert_eq!(one.weights, vec![(0, 1.0)]);
            let wide = Envelope::new(cx, &[-1e300, 1e300], &[0,1], 1e300, 1e-300).unwrap();
            assert_eq!(wide.upper_k, 1e300); assert_eq!(wide.weights, vec![(0, 0.0), (1, 1.0)]);
        });
    }

    #[test]
    fn nonfinite_inputs_and_wrong_selector_refuse() {
        let r = request();
        for tau in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            let mut maximum = Objective::parse(&J::parse(r#"{"max_vertices":[0]}"#).unwrap(), &r.surfaces, &r.mesh).unwrap();
            assert!(maximum.enable_upper_envelope(tau).is_err());
        }
        let mut mean = Objective::MeanWall("first-face".into());
        assert!(mean.enable_upper_envelope(0.1).is_err());
        with_cx(|cx| {
            assert!(Envelope::new(cx, &[f64::NAN], &[0], 300.0, 0.1).is_err());
            assert!(Envelope::new(cx, &[300.0], &[1], 300.0, 0.1).is_err());
            assert!(Envelope::new(cx, &[300.0], &[], 300.0, 0.1).is_err());
            assert!(Envelope::new(cx, &[300.0], &[0], 301.0, 0.1).is_err());
            assert!(Envelope::new(cx, &[300.0; 8], &[0,1,2,3,4,5,6,7], 300.0, f64::MAX).is_err());
        });
    }
}
