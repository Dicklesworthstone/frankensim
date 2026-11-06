//! G1/G3: joint spectral resolution against closed triangular equations.

use fs_math::c64::C64;
use fs_orbit::sideband::{
    PreparedPeriodicLinearization, SidebandControls, SidebandError, SidebandScales,
    refinement::{SidebandRefinementControls, solve_sideband_adaptive},
    sensitivity::{
        SidebandCarrierDerivativeControls, SidebandCarrierDerivativeResponse,
        refinement::{
            SidebandCarrierRefinementControls, SidebandCarrierRefinementTolerances,
            solve_sideband_carrier_derivative_adaptive,
        },
    },
};
use fs_orbit::{HbOrbit, IslandStateJacobian, OrbitProblem};
use std::cell::Cell;

const D: usize = 3;

fn controls(harmonics: usize) -> SidebandControls {
    SidebandControls {
        harmonics,
        samples_per_period: 128.max(2 * (3 * harmonics + 2) + 1),
        maximum_unknowns: 512,
        maximum_payload_bytes: 64 << 20,
        maximum_work: 2_000_000_000,
        maximum_base_residual: 1e-10,
        maximum_tangent_sampling_defect: 1e-9,
        maximum_relative_residual: 1e-10,
        maximum_omitted_residual: 1e-8,
    }
}

fn derivative_controls() -> SidebandCarrierDerivativeControls {
    SidebandCarrierDerivativeControls {
        maximum_payload_bytes: 64 << 20,
        maximum_work: 10_000_000_000,
        maximum_relative_residual: 1e-10,
        maximum_omitted_residual: 1e-8,
    }
}

/// Three ORIGINAL states with one-way periodic conversion x -> y -> z.
/// Two conversions can return to band zero, which distinguishes an actual
/// reduced-order solve from clipping a larger solution. Optional rational
/// diagonal ports have exact analytic derivatives and are never time-realized.
struct Chain {
    omega: f64,
    modulation: isize,
    rates: [f64; D],
    coupling: [f64; 2],
    phase: [f64; 2],
    memory: [f64; D],
    tau: [f64; D],
    scales: SidebandScales,
    source_calls: Cell<usize>,
    tangent_calls: Cell<usize>,
    derivative_calls: Cell<usize>,
    refuse_tangent: Cell<bool>,
}

impl Chain {
    fn ordinary() -> Self {
        Self {
            omega: 1.7,
            modulation: 2,
            rates: [0.8, 1.3, 0.7],
            coupling: [0.8, 0.6],
            phase: [0.37, -0.23],
            memory: [0.0; D],
            tau: [0.15, 0.7, 0.3],
            scales: SidebandScales {
                states: vec![1.0; D],
                equations: vec![1.0; D],
            },
            source_calls: Cell::new(0),
            tangent_calls: Cell::new(0),
            derivative_calls: Cell::new(0),
            refuse_tangent: Cell::new(false),
        }
    }

    fn derivative_sensitive() -> Self {
        Self {
            omega: 1.0,
            modulation: 9,
            rates: [1e-3, 1.0, 2.0],
            coupling: [1e-6, 0.0],
            ..Self::ordinary()
        }
    }

    fn base(&self) -> HbOrbit {
        HbOrbit {
            omega: self.omega,
            coeffs: vec![vec![C64::ZERO]; D],
            residual_trace: vec![],
            residual: f64::NAN,
            aft_samples: 0,
            sampling_defect: f64::NAN,
        }
    }

    fn tolerances(&self, value: f64) -> SidebandCarrierRefinementTolerances {
        SidebandCarrierRefinementTolerances {
            equation_tolerances: self.scales.equations.iter().map(|s| value * s).collect(),
            state_tolerances: self.scales.states.iter().map(|s| value * s).collect(),
            derivative_equation_tolerances: self
                .scales
                .equations
                .iter()
                .map(|s| value * s)
                .collect(),
            derivative_state_tolerances: self.scales.states.iter().map(|s| value * s).collect(),
        }
    }

    fn central(&self, harmonics: usize, amplitude: f64) -> Vec<C64> {
        let mut forcing = vec![C64::ZERO; D * (2 * harmonics + 1)];
        forcing[D * harmonics] = C64::from_re(amplitude * self.scales.equations[0]);
        forcing
    }

    fn gain(&self, link: usize, time: f64) -> f64 {
        self.coupling[link] * (self.modulation as f64 * self.omega * time + self.phase[link]).cos()
    }

    fn tangent(&self, time: f64, _: &[f64], out: &mut [f64]) -> bool {
        self.tangent_calls.set(self.tangent_calls.get() + 1);
        if self.refuse_tangent.get() {
            return false;
        }
        out.fill(0.0);
        for row in 0..D {
            out[row * D + row] =
                -self.rates[row] * self.scales.equations[row] / self.scales.states[row];
            if row > 0 {
                out[row * D + row - 1] = self.gain(row - 1, time) * self.scales.equations[row]
                    / self.scales.states[row - 1];
            }
        }
        true
    }

    fn port_diagonal(&self, row: usize, s: C64) -> C64 {
        s + s.scale(self.memory[row]) / (C64::ONE + s.scale(self.tau[row]))
    }

    fn port_slope(&self, row: usize, frequency: f64) -> C64 {
        let relaxation = C64::new(1.0, self.tau[row] * frequency);
        C64::new(0.0, 1.0) * (C64::ONE + C64::from_re(self.memory[row]) / (relaxation * relaxation))
    }

    fn derivative(&self, frequency: f64, out: &mut [C64]) -> bool {
        self.derivative_calls.set(self.derivative_calls.get() + 1);
        out.fill(C64::ZERO);
        for row in 0..D {
            out[row * D + row] = self
                .port_slope(row, frequency)
                .scale(self.scales.equations[row] / self.scales.states[row]);
        }
        true
    }

    fn force_at(&self, forcing: &[C64], h: usize, band: isize, row: usize) -> C64 {
        if band < -(h as isize) || band > h as isize {
            C64::ZERO
        } else {
            forcing[D * (band + h as isize) as usize + row].scale(1.0 / self.scales.equations[row])
        }
    }

    /// Closed triangular solution and exact chain rule, independent of sampled
    /// Fourier coefficients, numerical factors and frequency differencing.
    #[allow(clippy::too_many_arguments)]
    fn exact(
        &self,
        forcing: &[C64],
        slope: &[C64],
        h: usize,
        carrier: f64,
        band: isize,
        row: usize,
    ) -> [C64; 2] {
        let mut rhs = self.force_at(forcing, h, band, row);
        let mut rhs_slope = self.force_at(slope, h, band, row);
        if row > 0 {
            for sign in [-1_isize, 1] {
                let previous = self.exact(
                    forcing,
                    slope,
                    h,
                    carrier,
                    band - sign * self.modulation,
                    row - 1,
                );
                let phase = sign as f64 * self.phase[row - 1];
                let gain = C64::new(phase.cos(), phase.sin()).scale(0.5 * self.coupling[row - 1]);
                rhs = rhs + gain * previous[0];
                rhs_slope = rhs_slope + gain * previous[1];
            }
        }
        let frequency = carrier + band as f64 * self.omega;
        let denominator =
            C64::from_re(self.rates[row]) + self.port_diagonal(row, C64::new(0.0, frequency));
        let value = rhs / denominator;
        [
            value,
            (rhs_slope - self.port_slope(row, frequency) * value) / denominator,
        ]
    }
}

impl OrbitProblem for Chain {
    fn dim(&self) -> usize {
        D
    }
    fn autonomous(&self) -> bool {
        false
    }
    fn island(&self, time: f64, state: &[f64], out: &mut [f64]) {
        self.source_calls.set(self.source_calls.get() + 1);
        for row in 0..D {
            let mut rhs = -self.rates[row] * state[row] / self.scales.states[row];
            if row > 0 {
                rhs += self.gain(row - 1, time) * state[row - 1] / self.scales.states[row - 1];
            }
            out[row] = self.scales.equations[row] * rhs;
        }
    }
    fn port(&self, s: C64) -> Vec<C64> {
        let mut result = vec![C64::ZERO; D * D];
        for row in 0..D {
            result[row * D + row] = self
                .port_diagonal(row, s)
                .scale(self.scales.equations[row] / self.scales.states[row]);
        }
        result
    }
    fn island_state_jacobian(&self) -> Option<IslandStateJacobian<Self>> {
        Some(Self::tangent)
    }
}

fn at(values: &[C64], h: usize, band: isize, row: usize) -> C64 {
    values[D * (band + h as isize) as usize + row]
}

fn near(actual: C64, expected: C64) {
    assert!(
        (actual - expected).abs() <= 2e-8 * (1.0 + expected.abs()),
        "actual {actual:?}, expected {expected:?}"
    );
}

fn check_exact(
    source: &Chain,
    pair: &SidebandCarrierDerivativeResponse,
    h: usize,
    carrier: f64,
    forcing: &[C64],
    slope: &[C64],
) {
    for band in -(h as isize)..=h as isize {
        for row in 0..D {
            let expected = source.exact(forcing, slope, 1, carrier, band, row);
            near(
                at(&pair.response.coefficients, h, band, row)
                    .scale(1.0 / source.scales.states[row]),
                expected[0],
            );
            near(
                at(&pair.angular_derivative.coefficients, h, band, row)
                    .scale(1.0 / source.scales.states[row]),
                expected[1],
            );
        }
    }
}

#[test]
fn derivative_only_distant_equation_defect_prevents_false_response_adequacy() {
    let source = Chain::derivative_sensitive();
    let base = source.base();
    let forcing = source.central(1, 1e-3);
    let slope = vec![C64::ZERO; forcing.len()];
    let carrier = 1e-3;
    let mut c = controls(1);
    c.maximum_omitted_residual = 1e-5;
    let mut tolerances = source.tolerances(1e-5);
    tolerances.derivative_state_tolerances.fill(1e-6);
    let primal = solve_sideband_adaptive(
        &source,
        &base,
        carrier,
        &forcing,
        c,
        &source.scales,
        &SidebandRefinementControls {
            maximum_harmonics: 31,
            equation_tolerances: tolerances.equation_tolerances.clone(),
            state_tolerances: tolerances.state_tolerances.clone(),
        },
        || false,
    )
    .unwrap();
    assert_eq!(primal.harmonics, 3);
    // |X|=1/sqrt(2), |X_nu|=500: the original distant residual is about
    // 5e-7, but its derivative residual is about 3.5e-4 in the same row chart.
    let refinement = SidebandCarrierRefinementControls {
        maximum_harmonics: 31,
        tolerances,
    };
    let result = solve_sideband_carrier_derivative_adaptive(
        &source,
        &base,
        carrier,
        &forcing,
        &slope,
        c,
        &source.scales,
        derivative_controls(),
        &refinement,
        |w, out| source.derivative(w, out),
        || false,
    )
    .unwrap();
    assert_eq!(result.harmonics, 31);
    assert_eq!(
        result
            .refinement_trace
            .iter()
            .map(|x| x.response.harmonics)
            .collect::<Vec<_>>(),
        [1, 3, 7, 15, 31]
    );
    let plateau = result.refinement_trace[1];
    assert!(plateau.response.equation_residual_ratio < 1.0);
    assert!(plateau.response.state_change_ratio.unwrap() < 1.0);
    assert!(plateau.angular_derivative.omitted_residual_ratio < 1.0);
    assert!(plateau.angular_derivative.equation_residual_ratio > 1.0);
    assert!(!plateau.equations_resolved());
    assert!(result.refinement_trace.last().unwrap().is_resolved());
    check_exact(
        &source,
        &result.response,
        result.harmonics,
        carrier,
        &forcing,
        &slope,
    );

    let limited = SidebandCarrierRefinementControls {
        maximum_harmonics: 3,
        ..refinement
    };
    assert!(
        matches!(solve_sideband_carrier_derivative_adaptive(&source, &base, carrier,
        &forcing, &slope, c, &source.scales, derivative_controls(), &limited,
        |w, out| source.derivative(w, out), || false),
        Err(SidebandError::CarrierHarmonicsUnresolved { last })
            if last.response.equation_residual_ratio < 1.0
                && last.angular_derivative.equation_residual_ratio > 1.0)
    );
}

#[test]
fn rational_ports_and_authored_forcing_slopes_preserve_signed_phase_and_units() {
    let source = Chain {
        modulation: 1,
        memory: [0.6, 0.2, 0.4],
        scales: SidebandScales {
            states: vec![1e-9, 1e4, 1e-2],
            equations: vec![1e6, 1e-5, 1e3],
        },
        ..Chain::ordinary()
    };
    let normalized = [
        C64::new(-0.3, 0.7),
        C64::new(0.4, -0.1),
        C64::new(-0.2, 0.3),
        C64::new(0.9, -0.2),
        C64::new(-0.6, 0.5),
        C64::new(0.1, 0.2),
        C64::new(0.25, 0.35),
        C64::new(0.1, 0.8),
        C64::new(-0.4, 0.1),
    ];
    let forcing: Vec<_> = normalized
        .iter()
        .enumerate()
        .map(|(i, z)| z.scale(source.scales.equations[i % D]))
        .collect();
    let slope: Vec<_> = normalized
        .iter()
        .enumerate()
        .map(|(i, z)| {
            (*z * C64::new(0.2 - 0.03 * i as f64, -0.13)).scale(source.scales.equations[i % D])
        })
        .collect();
    let carrier = -0.37;
    let result = solve_sideband_carrier_derivative_adaptive(
        &source,
        &source.base(),
        carrier,
        &forcing,
        &slope,
        controls(1),
        &source.scales,
        derivative_controls(),
        &SidebandCarrierRefinementControls {
            maximum_harmonics: 7,
            tolerances: source.tolerances(1e-7),
        },
        |w, out| source.derivative(w, out),
        || false,
    )
    .unwrap();
    assert_eq!(result.harmonics, 7);
    assert!(result.refinement_trace.last().unwrap().is_resolved());
    check_exact(
        &source,
        &result.response,
        result.harmonics,
        carrier,
        &forcing,
        &slope,
    );
}

#[test]
fn one_cached_tangent_serves_actual_lower_and_upper_solves_and_shared_assessment_samples() {
    let source = Chain::ordinary();
    let base = source.base();
    let c = controls(15);
    let prepared =
        PreparedPeriodicLinearization::new(&source, &base, c, &source.scales, || false).unwrap();
    let preparation_calls = 3 * c.samples_per_period + 1;
    assert_eq!(source.source_calls.get(), preparation_calls);
    assert_eq!(source.tangent_calls.get(), preparation_calls);
    let tolerances = source.tolerances(1e-7);
    let carrier = 0.37;
    let mut candidates = Vec::new();
    for h in [1, 7, 15] {
        let forcing = source.central(h, 0.9);
        let slope = vec![C64::ZERO; forcing.len()];
        let tangents_before = source.tangent_calls.get();
        let derivatives_before = source.derivative_calls.get();
        let result = prepared
            .solve_carrier_derivative_assessed_at_order(
                h,
                carrier,
                &forcing,
                &slope,
                derivative_controls(),
                &tolerances,
                |w, out| source.derivative(w, out),
                || false,
            )
            .unwrap();
        assert_eq!(source.source_calls.get(), preparation_calls);
        assert_eq!(
            source.tangent_calls.get() - tangents_before,
            6 * c.samples_per_period + 4
        );
        assert_eq!(
            source.derivative_calls.get() - derivatives_before,
            2 * h + 1
        );
        assert_eq!(
            result.check.response.assessment_samples,
            result.check.angular_derivative.assessment_samples
        );
        candidates.push(result);
    }
    let low = &candidates[0];
    let high = &candidates[1];
    assert!(at(&low.response.response.coefficients, 1, 0, 2).abs() < 1e-12);
    assert!(at(&high.response.response.coefficients, 7, 0, 2).abs() > 1e-3);
    // Merely clipping a high-order response would retain its nonzero central z
    // coefficient and fail the low-order assertion above.
    assert!(!low.check.equations_resolved());
    assert!(high.check.equations_resolved());
    assert!(!high.check.is_resolved());
    let mut compared = high.clone();
    compared.compare_with(low, &tolerances).unwrap();
    assert!(compared.check.response.state_change_ratio.unwrap() > 1.0);
    let mut final_order = candidates[2].clone();
    final_order.compare_with(high, &tolerances).unwrap();
    assert!(final_order.check.is_resolved());
    let forcing = source.central(1, 0.9);
    check_exact(
        &source,
        &final_order.response,
        final_order.harmonics,
        carrier,
        &forcing,
        &vec![C64::ZERO; forcing.len()],
    );
    let mut wrong_carrier = final_order.clone();
    wrong_carrier.carrier_angular_frequency = carrier + 0.1;
    assert!(matches!(
        wrong_carrier.compare_with(high, &tolerances),
        Err(SidebandError::Invalid(_))
    ));
}

#[test]
fn malformed_data_source_refusals_and_assessment_cancellation_do_not_publish_a_pair() {
    let source = Chain::ordinary();
    let base = source.base();
    let prepared =
        PreparedPeriodicLinearization::new(&source, &base, controls(7), &source.scales, || false)
            .unwrap();
    let forcing = source.central(1, 0.9);
    let slope = vec![C64::ZERO; forcing.len()];
    let tolerances = source.tolerances(1e-7);
    assert!(matches!(
        prepared.solve_carrier_derivative_assessed_at_order(
            1,
            0.37,
            &forcing[..forcing.len() - 1],
            &slope,
            derivative_controls(),
            &tolerances,
            |w, out| source.derivative(w, out),
            || false
        ),
        Err(SidebandError::Invalid(_))
    ));
    assert!(matches!(
        prepared.solve_carrier_derivative_assessed_at_order(
            8,
            0.37,
            &forcing,
            &slope,
            derivative_controls(),
            &tolerances,
            |w, out| source.derivative(w, out),
            || false
        ),
        Err(SidebandError::Invalid(_))
    ));
    assert!(matches!(
        prepared.solve_carrier_derivative_assessed_at_order(
            1,
            0.37,
            &forcing,
            &slope,
            derivative_controls(),
            &tolerances,
            |_, out| {
                out[0] = C64::new(0.0, 1.0);
                true
            },
            || false
        ),
        Err(SidebandError::CarrierDerivativeRefused { .. })
    ));
    let before = source.tangent_calls.get();
    assert!(matches!(
        prepared.solve_carrier_derivative_assessed_at_order(
            1,
            0.37,
            &forcing,
            &slope,
            derivative_controls(),
            &tolerances,
            |w, out| source.derivative(w, out),
            || source.tangent_calls.get() >= before + 2
        ),
        Err(SidebandError::Cancelled)
    ));
    assert_eq!(source.tangent_calls.get(), before + 2);
    source.refuse_tangent.set(true);
    assert!(matches!(
        prepared.solve_carrier_derivative_assessed_at_order(
            1,
            0.37,
            &forcing,
            &slope,
            derivative_controls(),
            &tolerances,
            |w, out| source.derivative(w, out),
            || false
        ),
        Err(SidebandError::TangentRefused { .. })
    ));
}

#[test]
fn added_derivative_work_is_cumulative_across_actual_orders() {
    let source = Chain::derivative_sensitive();
    let forcing = source.central(1, 1e-3);
    let slope = vec![C64::ZERO; forcing.len()];
    let mut c = controls(1);
    c.maximum_omitted_residual = 1e-5;
    let mut extra = derivative_controls();
    // Each individual round through N=15 fits 250M under the declared score;
    // their combined assessments exceed it. The core has its own larger cap.
    extra.maximum_work = 250_000_000;
    let mut tolerances = source.tolerances(1e-5);
    tolerances.derivative_state_tolerances.fill(1e-6);
    let result = solve_sideband_carrier_derivative_adaptive(
        &source,
        &source.base(),
        1e-3,
        &forcing,
        &slope,
        c,
        &source.scales,
        extra,
        &SidebandCarrierRefinementControls {
            maximum_harmonics: 15,
            tolerances,
        },
        |w, out| source.derivative(w, out),
        || false,
    );
    assert!(
        matches!(result, Err(SidebandError::Budget { quantity, requested, cap })
        if quantity == "adaptive carrier derivative work"
            && cap == extra.maximum_work as u128 && requested > cap)
    );
    assert!(source.source_calls.get() > 3 * c.samples_per_period + 1);
}
