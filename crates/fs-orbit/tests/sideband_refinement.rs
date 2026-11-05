//! G1/G3: original-equation refinement against an exact distant-band cascade.

use fs_math::c64::C64;
use fs_orbit::sideband::refinement::{SidebandRefinementControls, solve_sideband_adaptive};
use fs_orbit::sideband::{
    PreparedPeriodicLinearization, SidebandControls, SidebandError, SidebandScales,
};
use fs_orbit::{HbOrbit, IslandStateJacobian, OrbitProblem};
use std::cell::Cell;

const INITIAL_SAMPLES: usize = 64;

fn controls() -> SidebandControls {
    SidebandControls {
        harmonics: 1,
        samples_per_period: INITIAL_SAMPLES,
        maximum_unknowns: 256,
        maximum_payload_bytes: 64_000_000,
        maximum_work: 1_000_000_000,
        maximum_base_residual: 1e-10,
        maximum_tangent_sampling_defect: 1e-9,
        maximum_relative_residual: 1e-10,
        maximum_omitted_residual: 1e-8,
    }
}

/// In normalized physical coordinates, x' = -a*x and
/// y' = -b*y + c*cos(9*omega*t + phase)*x. The base orbit is exactly zero.
/// The probe is external equation forcing, never part of this base source.
struct DistantCascade {
    units: SidebandScales,
    phase: f64,
    tangent_kind: u8,
    partial_source: bool,
    source_calls: Cell<usize>,
    tangent_calls: Cell<usize>,
}

impl DistantCascade {
    const OMEGA: f64 = 1.7;
    const A: f64 = 0.8;
    const B: f64 = 1.3;
    const C: f64 = 0.7;
    const MODULATION: isize = 9;

    fn ordinary() -> Self {
        Self {
            units: SidebandScales {
                states: vec![1.0; 2],
                equations: vec![1.0; 2],
            },
            phase: 0.37,
            tangent_kind: 0,
            partial_source: false,
            source_calls: Cell::new(0),
            tangent_calls: Cell::new(0),
        }
    }

    fn base(&self) -> HbOrbit {
        HbOrbit {
            omega: Self::OMEGA,
            coeffs: vec![vec![C64::ZERO]; 2],
            residual_trace: vec![],
            // Fresh original-equation checks must ignore stale diagnostics.
            residual: f64::NAN,
            aft_samples: 0,
            sampling_defect: f64::NAN,
        }
    }

    fn refinement(&self, maximum_harmonics: usize) -> SidebandRefinementControls {
        SidebandRefinementControls {
            maximum_harmonics,
            equation_tolerances: self.units.equations.iter().map(|s| 1e-8 * s).collect(),
            state_tolerances: self.units.states.iter().map(|s| 1e-8 * s).collect(),
        }
    }

    fn derivative(&self, time: f64, _: &[f64], out: &mut [f64]) -> bool {
        self.tangent_calls.set(self.tangent_calls.get() + 1);
        if self.tangent_kind == 2
            || (self.tangent_kind == 4 && self.tangent_calls.get() > 3 * INITIAL_SAMPLES + 1)
        {
            return false;
        }
        out[0] = -Self::A * self.units.equations[0] / self.units.states[0];
        out[1] = 0.0;
        out[2] = Self::C
            * (Self::MODULATION as f64 * Self::OMEGA * time + self.phase).cos()
            * self.units.equations[1]
            / self.units.states[0];
        if self.tangent_kind != 3 {
            out[3] = -Self::B * self.units.equations[1] / self.units.states[1];
        }
        true
    }

    fn central_forcing(&self) -> Vec<C64> {
        let mut forcing = vec![C64::ZERO; 6];
        forcing[2] = C64::new(0.9, -0.2).scale(self.units.equations[0]);
        forcing
    }

    fn force_at(&self, forcing: &[C64], band: isize, row: usize) -> C64 {
        if (-1..=1).contains(&band) {
            forcing[2 * (band + 1) as usize + row].scale(1.0 / self.units.equations[row])
        } else {
            C64::ZERO
        }
    }

    /// Closed triangular equations, independently solved without any Fourier
    /// projection, generic matrix factorization, or time integration.
    fn expected(&self, forcing: &[C64], carrier: f64, band: isize) -> [C64; 2] {
        let x =
            |k| self.force_at(forcing, k, 0) / C64::new(Self::A, carrier + k as f64 * Self::OMEGA);
        let positive = C64::new(self.phase.cos(), self.phase.sin()).scale(0.5 * Self::C);
        let negative = C64::new(self.phase.cos(), -self.phase.sin()).scale(0.5 * Self::C);
        let y = (self.force_at(forcing, band, 1)
            + positive * x(band - Self::MODULATION)
            + negative * x(band + Self::MODULATION))
            / C64::new(Self::B, carrier + band as f64 * Self::OMEGA);
        [x(band), y]
    }
}

impl OrbitProblem for DistantCascade {
    fn dim(&self) -> usize {
        2
    }

    fn autonomous(&self) -> bool {
        false
    }

    fn island(&self, time: f64, state: &[f64], out: &mut [f64]) {
        self.source_calls.set(self.source_calls.get() + 1);
        let x = state[0] / self.units.states[0];
        out[0] = -self.units.equations[0] * Self::A * x;
        if !self.partial_source {
            out[1] = self.units.equations[1]
                * (-Self::B * state[1] / self.units.states[1]
                    + Self::C
                        * (Self::MODULATION as f64 * Self::OMEGA * time + self.phase).cos()
                        * x);
        }
    }

    fn port(&self, s: C64) -> Vec<C64> {
        vec![
            s.scale(self.units.equations[0] / self.units.states[0]),
            C64::ZERO,
            C64::ZERO,
            s.scale(self.units.equations[1] / self.units.states[1]),
        ]
    }

    fn island_state_jacobian(&self) -> Option<IslandStateJacobian<Self>> {
        if self.tangent_kind == 1 {
            None
        } else {
            Some(Self::derivative)
        }
    }
}

fn coefficient(values: &[C64], harmonics: usize, band: isize, row: usize) -> C64 {
    values[2 * (band + harmonics as isize) as usize + row]
}

fn close(actual: C64, expected: C64) {
    assert!(
        (actual - expected).abs() <= 5e-9 * (1.0 + expected.abs()),
        "actual {actual:?}, expected {expected:?}"
    );
}

#[test]
fn full_equation_refinement_escapes_a_false_low_order_plateau() {
    let source = DistantCascade::ordinary();
    let base = source.base();
    let carrier = 0.43;
    let forcing = source.central_forcing();

    // Both independently solved low orders appear converged under the existing
    // finite omitted-band check. The missing physical response is at +/-9.
    for harmonics in [1, 3] {
        let mut fixed = controls();
        fixed.harmonics = harmonics;
        let prepared =
            PreparedPeriodicLinearization::new(&source, &base, fixed, &source.units, || false)
                .unwrap();
        let mut padded = vec![C64::ZERO; 2 * (2 * harmonics + 1)];
        padded[2 * harmonics] = forcing[2];
        let result = prepared
            .factor(carrier, || false)
            .unwrap()
            .solve(&padded, || false)
            .unwrap();
        assert!(result.omitted_residual < 1e-10);
        close(
            coefficient(&result.coefficients, harmonics, 0, 0),
            source.expected(&forcing, carrier, 0)[0],
        );
        for band in -(harmonics as isize)..=harmonics as isize {
            close(
                coefficient(&result.coefficients, harmonics, band, 1),
                C64::ZERO,
            );
        }
    }

    let result = solve_sideband_adaptive(
        &source,
        &base,
        carrier,
        &forcing,
        controls(),
        &source.units,
        &source.refinement(31),
        || false,
    )
    .unwrap();
    assert_eq!(result.harmonics, 31);
    assert_eq!(
        result
            .refinement_trace
            .iter()
            .map(|check| check.harmonics)
            .collect::<Vec<_>>(),
        [1, 3, 7, 15, 31]
    );
    let plateau = &result.refinement_trace[1];
    assert!(plateau.state_change_ratio.unwrap() < 1.0);
    assert!(plateau.omitted_residual_ratio < 1.0);
    assert!(plateau.equation_residual_ratio > 1.0);
    assert!(result.refinement_trace[3].equation_residual_ratio < 1.0);
    assert!(result.refinement_trace[3].state_change_ratio.unwrap() > 1.0);
    let accepted = result.refinement_trace.last().unwrap();
    assert!(accepted.assessment_samples > 2 * result.harmonics);
    assert!(accepted.equation_residual_ratio <= 1.0);
    assert!(accepted.sampling_defect_ratio <= 1.0);
    assert!(accepted.omitted_residual_ratio <= 1.0);
    assert!(accepted.state_change_ratio.unwrap() <= 1.0);
    assert!(result.base_residual <= controls().maximum_base_residual);
    assert!(result.tangent_sampling_defect <= controls().maximum_tangent_sampling_defect);
    for band in -(result.harmonics as isize)..=result.harmonics as isize {
        for (row, expected) in source
            .expected(&forcing, carrier, band)
            .into_iter()
            .enumerate()
        {
            close(
                coefficient(&result.response.coefficients, result.harmonics, band, row),
                expected,
            );
        }
    }
    assert!(coefficient(&result.response.coefficients, result.harmonics, -9, 1).abs() > 1e-3);
    assert!(coefficient(&result.response.coefficients, result.harmonics, 9, 1).abs() > 1e-3);
}

#[test]
fn signed_complex_forcing_and_physical_units_survive_order_changes() {
    let mut source = DistantCascade::ordinary();
    source.phase = -0.73;
    source.units = SidebandScales {
        states: vec![1e-8, 1e4],
        equations: vec![1e3, 1e-5],
    };
    let normalized = [
        C64::new(-0.3, 0.7),
        C64::new(0.4, -0.1),
        C64::new(0.9, -0.2),
        C64::new(-0.6, 0.5),
        C64::new(0.25, 0.35),
        C64::new(0.1, 0.8),
    ];
    let forcing: Vec<_> = normalized
        .iter()
        .enumerate()
        .map(|(index, value)| value.scale(source.units.equations[index % 2]))
        .collect();
    let carrier = -0.43;
    let result = solve_sideband_adaptive(
        &source,
        &source.base(),
        carrier,
        &forcing,
        controls(),
        &source.units,
        &source.refinement(31),
        || false,
    )
    .unwrap();
    assert_eq!(result.harmonics, 31);
    for band in -(result.harmonics as isize)..=result.harmonics as isize {
        for (row, expected) in source
            .expected(&forcing, carrier, band)
            .into_iter()
            .enumerate()
        {
            let physical = coefficient(&result.response.coefficients, result.harmonics, band, row);
            close(physical.scale(1.0 / source.units.states[row]), expected);
        }
    }
    // Authored input at both +/-1 produces distinct physical lines at +/-10;
    // centering, conjugating, folding or doubling the forcing cannot pass.
    for band in [-10, -9, -8, 8, 9, 10] {
        assert!(coefficient(&result.response.coefficients, result.harmonics, band, 1).abs() > 1.0);
    }
}

#[test]
fn maximum_order_refuses_both_a_false_plateau_and_one_newly_resolved_candidate() {
    let source = DistantCascade::ordinary();
    for maximum in [3, 15] {
        let error = solve_sideband_adaptive(
            &source,
            &source.base(),
            0.43,
            &source.central_forcing(),
            controls(),
            &source.units,
            &source.refinement(maximum),
            || false,
        )
        .unwrap_err();
        match error {
            SidebandError::HarmonicsUnresolved { last } => {
                assert_eq!(last.harmonics, maximum);
                if maximum == 3 {
                    assert!(last.equation_residual_ratio > 1.0);
                    assert!(last.state_change_ratio.unwrap() < 1.0);
                } else {
                    assert!(last.equation_residual_ratio < 1.0);
                    assert!(last.state_change_ratio.unwrap() > 1.0);
                }
            }
            other => panic!("expected unresolved response order, got {other:?}"),
        }
    }
}

#[test]
fn cancellation_and_declared_source_refusals_propagate_without_refinement() {
    let source = DistantCascade::ordinary();
    let cancelled = solve_sideband_adaptive(
        &source,
        &source.base(),
        0.43,
        &source.central_forcing(),
        controls(),
        &source.units,
        &source.refinement(31),
        || true,
    );
    assert!(matches!(cancelled, Err(SidebandError::Cancelled)));
    assert_eq!(source.source_calls.get(), 0);
    assert_eq!(source.tangent_calls.get(), 0);
    let cancelled_during_source = solve_sideband_adaptive(
        &source,
        &source.base(),
        0.43,
        &source.central_forcing(),
        controls(),
        &source.units,
        &source.refinement(31),
        || source.source_calls.get() >= 16,
    );
    assert!(matches!(
        cancelled_during_source,
        Err(SidebandError::Cancelled)
    ));
    assert!(source.source_calls.get() >= 16);
    assert!(source.source_calls.get() < 64);

    for tangent_kind in [1, 2, 3] {
        let mut source = DistantCascade::ordinary();
        source.tangent_kind = tangent_kind;
        let result = solve_sideband_adaptive(
            &source,
            &source.base(),
            0.43,
            &source.central_forcing(),
            controls(),
            &source.units,
            &source.refinement(31),
            || false,
        );
        if tangent_kind == 1 {
            assert!(matches!(result, Err(SidebandError::TangentUnavailable)));
        } else {
            assert!(matches!(result, Err(SidebandError::TangentRefused { .. })));
            assert_eq!(source.tangent_calls.get(), 1);
        }
    }
    // Let both preparation grids finish, then refuse the first independent
    // original-equation assessment. A wrapper must not suppress this refusal
    // or reuse only the already truncated Fourier tangent for its check.
    let mut source = DistantCascade::ordinary();
    source.tangent_kind = 4;
    assert!(matches!(
        solve_sideband_adaptive(
            &source,
            &source.base(),
            0.43,
            &source.central_forcing(),
            controls(),
            &source.units,
            &source.refinement(31),
            || false,
        ),
        Err(SidebandError::TangentRefused { sample: 0, .. })
    ));
    assert_eq!(source.tangent_calls.get(), 3 * INITIAL_SAMPLES + 2);
    let mut source = DistantCascade::ordinary();
    source.partial_source = true;
    assert!(matches!(
        solve_sideband_adaptive(
            &source,
            &source.base(),
            0.43,
            &source.central_forcing(),
            controls(),
            &source.units,
            &source.refinement(31),
            || false,
        ),
        Err(SidebandError::SourceRefused { .. })
    ));
}

#[test]
fn adaptive_work_is_cumulative_across_the_actual_refinement_chain() {
    let source = DistantCascade::ordinary();
    let mut limited = controls();
    // Under the declared numerical work model, every individual round through
    // N=15 fits 40M; the 1,3,7,15 chain costs about 46.95M. Reapplying this cap
    // independently to each round would let the complete chain through.
    limited.maximum_work = 40_000_000;
    let result = solve_sideband_adaptive(
        &source,
        &source.base(),
        0.43,
        &source.central_forcing(),
        limited,
        &source.units,
        &source.refinement(15),
        || false,
    );
    match result {
        Err(SidebandError::Budget {
            quantity,
            requested,
            cap,
        }) => {
            assert_eq!(quantity, "adaptive numerical work");
            assert_eq!(cap, limited.maximum_work as u128);
            assert!(requested > cap);
        }
        other => panic!("expected cumulative work refusal, got {other:?}"),
    }
    assert!(source.source_calls.get() > 3 * INITIAL_SAMPLES + 1);
}

#[test]
fn malformed_forcing_and_missing_second_order_are_rejected() {
    let source = DistantCascade::ordinary();
    for forcing in [vec![C64::ZERO; 4], vec![C64::ZERO; 8]] {
        assert!(matches!(
            solve_sideband_adaptive(
                &source,
                &source.base(),
                0.43,
                &forcing,
                controls(),
                &source.units,
                &source.refinement(31),
                || false,
            ),
            Err(SidebandError::Invalid(_))
        ));
    }
    assert!(matches!(
        solve_sideband_adaptive(
            &source,
            &source.base(),
            0.43,
            &source.central_forcing(),
            controls(),
            &source.units,
            &source.refinement(1),
            || false,
        ),
        Err(SidebandError::Invalid(_))
    ));
}
