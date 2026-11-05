//! G1/G3: complete periodic tangents versus analytic sidebands and real ODEs.

use fs_math::c64::C64;
use fs_orbit::sideband::{
    PreparedPeriodicLinearization, SidebandControls, SidebandError, SidebandScales,
};
use fs_orbit::{HbOrbit, IslandStateJacobian, OrbitProblem};
use std::cell::Cell;

fn controls(harmonics: usize) -> SidebandControls {
    SidebandControls {
        harmonics,
        samples_per_period: 64,
        maximum_unknowns: 256,
        maximum_payload_bytes: 64_000_000,
        maximum_work: 1_000_000_000,
        maximum_base_residual: 1e-10,
        maximum_tangent_sampling_defect: 1e-10,
        maximum_relative_residual: 1e-12,
        maximum_omitted_residual: 1e-9,
    }
}
fn units(d: usize) -> SidebandScales {
    SidebandScales {
        states: vec![1.0; d],
        equations: vec![1.0; d],
    }
}
fn orbit(omega: f64, coeffs: Vec<Vec<C64>>) -> HbOrbit {
    // Deliberately invalid old diagnostics: preparation must recompute physics.
    HbOrbit {
        omega,
        coeffs,
        residual_trace: vec![],
        residual: f64::NAN,
        aft_samples: 0,
        sampling_defect: f64::NAN,
    }
}
fn at(values: &[C64], d: usize, h: usize, k: isize, row: usize) -> C64 {
    values[(k + h as isize) as usize * d + row]
}
fn close(a: C64, b: C64, tolerance: f64) {
    assert!((a - b).abs() < tolerance, "actual {a:?}, expected {b:?}");
}

struct CoupledMechanicalPort;
impl CoupledMechanicalPort {
    fn tangent(&self, _: f64, _: &[f64], out: &mut [f64]) -> bool {
        out.fill(0.0);
        true
    }
}
impl OrbitProblem for CoupledMechanicalPort {
    fn dim(&self) -> usize {
        2
    }
    fn autonomous(&self) -> bool {
        false
    }
    fn island(&self, _: f64, _: &[f64], out: &mut [f64]) {
        out.fill(0.0);
    }
    fn port(&self, s: C64) -> Vec<C64> {
        // Full symmetric mass/stiffness, nonproportional positive loss,
        // skew gyroscopic coupling and one rank-one Maxwell impedance.
        let relaxation = s.scale(1.2) / (C64::ONE + s.scale(0.4));
        let mass = [2.0, 0.3, 0.3, 1.5];
        let damping = [0.3, 0.47, -0.33, 0.5];
        let stiffness = [4.0, 0.2, 0.2, 2.0];
        let memory = [1.0, -0.5, -0.5, 0.25];
        (0..4)
            .map(|j| {
                (s * s).scale(mass[j])
                    + s.scale(damping[j])
                    + C64::from_re(stiffness[j])
                    + relaxation.scale(memory[j])
            })
            .collect()
    }
    fn island_state_jacobian(&self) -> Option<IslandStateJacobian<Self>> {
        Some(Self::tangent)
    }
}

#[test]
fn full_coupled_frequency_dependent_port_preserves_signed_gyroscopic_and_memory_response() {
    let source = CoupledMechanicalPort;
    let base = orbit(1.7, vec![vec![C64::ZERO; 2]; 2]);
    let scales = SidebandScales {
        states: vec![1e-6, 1e3],
        equations: vec![1e2, 1e-4],
    };
    let prepared =
        PreparedPeriodicLinearization::new(&source, &base, controls(1), &scales, || false).unwrap();
    let forcing = [
        C64::new(0.4, -0.1),
        C64::new(-0.2, 0.3),
        C64::new(0.7, 0.2),
        C64::new(0.1, -0.4),
        C64::new(-0.1, 0.5),
        C64::new(0.3, 0.2),
    ];
    for carrier in [-2.1, 0.0, 0.9] {
        let result = prepared
            .factor(carrier, || false)
            .unwrap()
            .solve(&forcing, || false)
            .unwrap();
        for band in 0..3 {
            let w = carrier + (band as f64 - 1.0) * base.omega;
            let a = source.port(C64::new(0.0, w));
            // Independent closed two-coordinate response, with no fitted Q,
            // diagonal modal assumption or conversion to a first-order port.
            let determinant = a[0] * a[3] - a[1] * a[2];
            let f = &forcing[2 * band..2 * band + 2];
            let expected = [
                (a[3] * f[0] - a[1] * f[1]) / determinant,
                (a[0] * f[1] - a[2] * f[0]) / determinant,
            ];
            for (row, &value) in expected.iter().enumerate() {
                close(result.coefficients[2 * band + row], value, 2e-12);
            }
        }
        assert_eq!(result.omitted_residual, 0.0);
    }
    let tiny = orbit(1.7, vec![vec![C64::from_re(1e-100)], vec![C64::ZERO]]);
    let erased = SidebandScales {
        states: vec![1e300; 2],
        equations: vec![1e300; 2],
    };
    assert!(matches!(
        PreparedPeriodicLinearization::new(&source, &tiny, controls(1), &erased, || false),
        Err(SidebandError::Invalid(
            "normalized base coefficient is not representable"
        ))
    ));
}

struct Cascade {
    scales: SidebandScales,
    phase: f64,
    calls: Cell<usize>,
    tangent: u8,
}
impl Cascade {
    const OMEGA: f64 = 1.7;
    fn derivative(&self, _: f64, x: &[f64], out: &mut [f64]) -> bool {
        if self.tangent == 2 {
            return false;
        }
        out[0] = -self.scales.equations[0] / self.scales.states[0];
        out[1] = 0.0;
        out[2] =
            6.0 * (x[0] / self.scales.states[0]) * self.scales.equations[1] / self.scales.states[0];
        if self.tangent != 3 {
            out[3] = -2.0 * self.scales.equations[1] / self.scales.states[1];
        }
        true
    }
    fn base(&self) -> HbOrbit {
        let x =
            C64::new(self.phase.cos(), self.phase.sin()).scale(0.5) / C64::new(1.0, Self::OMEGA);
        let y0 = C64::from_re(3.0 * x.norm_sq());
        let y2 = x * x * C64::from_re(3.0) / C64::new(2.0, 2.0 * Self::OMEGA);
        orbit(
            Self::OMEGA,
            vec![
                vec![C64::ZERO, x.scale(self.scales.states[0]), C64::ZERO],
                vec![
                    y0.scale(self.scales.states[1]),
                    C64::ZERO,
                    y2.scale(self.scales.states[1]),
                ],
            ],
        )
    }
}
impl OrbitProblem for Cascade {
    fn dim(&self) -> usize {
        2
    }
    fn autonomous(&self) -> bool {
        false
    }
    fn island(&self, t: f64, x: &[f64], out: &mut [f64]) {
        self.calls.set(self.calls.get() + 1);
        let x0 = x[0] / self.scales.states[0];
        out[0] = self.scales.equations[0] * ((Self::OMEGA * t + self.phase).cos() - x0);
        out[1] = self.scales.equations[1] * (3.0 * x0 * x0 - 2.0 * x[1] / self.scales.states[1]);
    }
    fn port(&self, s: C64) -> Vec<C64> {
        vec![
            s.scale(self.scales.equations[0] / self.scales.states[0]),
            C64::ZERO,
            C64::ZERO,
            s.scale(self.scales.equations[1] / self.scales.states[1]),
        ]
    }
    fn island_state_jacobian(&self) -> Option<IslandStateJacobian<Self>> {
        if self.tangent == 1 {
            None
        } else {
            Some(Self::derivative)
        }
    }
}

#[test]
fn nonlinear_orbit_generates_both_signed_sidebands_with_phase_units_and_reused_factors() {
    for scales in [
        units(2),
        SidebandScales {
            states: vec![1e-12, 1e6],
            equations: vec![1e8, 1e-7],
        },
    ] {
        for phase in [0.0, 0.63] {
            let source = Cascade {
                scales: scales.clone(),
                phase,
                calls: Cell::new(0),
                tangent: 0,
            };
            let prepared = PreparedPeriodicLinearization::new(
                &source,
                &source.base(),
                controls(2),
                &scales,
                || false,
            )
            .unwrap();
            let calls = source.calls.get();
            assert_eq!(calls, 64 + 129);
            assert!(prepared.base_residual < 1e-13);
            for carrier in [-0.4, 0.0, 0.6, Cascade::OMEGA] {
                let factor = prepared.factor(carrier, || false).unwrap();
                for drive in [C64::new(0.7, -0.2), C64::new(-0.1, 0.8)] {
                    let mut forcing = vec![C64::ZERO; 10];
                    forcing[4] = drive.scale(scales.equations[0]);
                    let response = factor.solve(&forcing, || false).unwrap();
                    let dx = drive / C64::new(1.0, carrier);
                    let base_x = source.base().coeffs[0][1].scale(1.0 / scales.states[0]);
                    for k in -2_isize..=2 {
                        close(
                            at(&response.coefficients, 2, 2, k, 0).scale(1.0 / scales.states[0]),
                            if k == 0 { dx } else { C64::ZERO },
                            2e-13,
                        );
                        let expected = match k {
                            1 => {
                                base_x * dx * C64::from_re(6.0)
                                    / C64::new(2.0, carrier + Cascade::OMEGA)
                            }
                            -1 => {
                                base_x.conj() * dx * C64::from_re(6.0)
                                    / C64::new(2.0, carrier - Cascade::OMEGA)
                            }
                            _ => C64::ZERO,
                        };
                        close(
                            at(&response.coefficients, 2, 2, k, 1).scale(1.0 / scales.states[1]),
                            expected,
                            2e-13,
                        );
                    }
                    assert!(response.relative_residual < 1e-14);
                    assert!(response.omitted_residual < 1e-13);
                }
            }
            assert_eq!(
                source.calls.get(),
                calls,
                "all frequency factors and RHS reuse the original tangent"
            );
        }
    }
}

struct Parametric {
    modulation: f64,
    damping: f64,
    modulation_order: f64,
}
impl OrbitProblem for Parametric {
    fn dim(&self) -> usize {
        2
    }
    fn autonomous(&self) -> bool {
        false
    }
    fn island(&self, t: f64, x: &[f64], out: &mut [f64]) {
        out[0] = x[1];
        out[1] = -(1.0 + self.modulation * (2.0 * self.modulation_order * t).cos()) * x[0]
            - self.damping * x[1];
    }
    fn island_state_jacobian(&self) -> Option<IslandStateJacobian<Self>> {
        Some(|source, t, _, out| {
            out.copy_from_slice(&[
                0.0,
                1.0,
                -(1.0 + source.modulation * (2.0 * source.modulation_order * t).cos()),
                -source.damping,
            ]);
            true
        })
    }
}
fn parametric_response(
    source: &Parametric,
    h: usize,
    omitted: f64,
) -> Result<Vec<C64>, SidebandError> {
    let mut c = controls(h);
    c.maximum_omitted_residual = omitted;
    let prepared = PreparedPeriodicLinearization::new(
        source,
        &orbit(2.0, vec![vec![C64::ZERO]; 2]),
        c,
        &units(2),
        || false,
    )?;
    let factor = prepared.factor(1.0, || false)?;
    let mut force = vec![C64::ZERO; 2 * (2 * h + 1)];
    force[2 * h + 1] = C64::ONE;
    Ok(factor.solve(&force, || false)?.coefficients)
}

fn time_oracle(source: &Parametric) -> (C64, C64) {
    // Independent real physical q/v RK4, driven by cos(t), then demodulation.
    let per_period = 1024;
    let dt = core::f64::consts::TAU / per_period as f64;
    let mut state = [0.0; 2];
    let rhs = |t: f64, x: [f64; 2]| {
        [
            x[1],
            t.cos() - (1.0 + source.modulation * (2.0 * t).cos()) * x[0] - source.damping * x[1],
        ]
    };
    let mut observed = [C64::ZERO; 2];
    for step in 0..100 * per_period {
        let t = step as f64 * dt;
        if step >= 99 * per_period {
            observed[0] =
                observed[0] + C64::new(t.cos(), -t.sin()).scale(2.0 * state[0] / per_period as f64);
            observed[1] = observed[1]
                + C64::new((3.0 * t).cos(), -(3.0 * t).sin())
                    .scale(2.0 * state[0] / per_period as f64);
        }
        let k1 = rhs(t, state);
        let k2 = rhs(
            t + 0.5 * dt,
            std::array::from_fn(|i| state[i] + 0.5 * dt * k1[i]),
        );
        let k3 = rhs(
            t + 0.5 * dt,
            std::array::from_fn(|i| state[i] + 0.5 * dt * k2[i]),
        );
        let k4 = rhs(t + dt, std::array::from_fn(|i| state[i] + dt * k3[i]));
        state = std::array::from_fn(|i| {
            state[i] + dt * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0
        });
    }
    (observed[0], observed[1])
}

#[test]
fn parametric_resonance_and_frequency_conversion_converge_to_the_original_ode() {
    let source = Parametric {
        modulation: 0.2,
        damping: 0.2,
        modulation_order: 1.0,
    };
    assert!(matches!(
        parametric_response(&source, 0, 1e-6),
        Err(SidebandError::OmittedSidebands { .. })
    ));
    let coarse = parametric_response(&source, 2, 1.0).unwrap();
    let medium = parametric_response(&source, 4, 1.0).unwrap();
    let fine = parametric_response(&source, 6, 1e-9).unwrap();
    let aggregate = |a: &[C64], h: usize| {
        (
            at(a, 2, h, 0, 0) + at(a, 2, h, -1, 0).conj(),
            at(a, 2, h, 1, 0) + at(a, 2, h, -2, 0).conj(),
        )
    };
    let (a, b, c) = (
        aggregate(&coarse, 2),
        aggregate(&medium, 4),
        aggregate(&fine, 6),
    );
    assert!((b.0 - c.0).abs() < (a.0 - c.0).abs() * 0.001);
    let oracle = time_oracle(&source);
    close(c.0, oracle.0, 3e-7);
    close(c.1, oracle.1, 3e-7);
    assert!(
        c.0.abs() > 6.0,
        "mean-tangent oscillator would have amplitude 5"
    );
    assert!(
        c.1.abs() > 0.01,
        "actual stiffness modulation generates a third harmonic"
    );
}

#[test]
fn stale_base_missing_harmonics_aliases_and_singular_poles_are_not_finite_successes() {
    let source = Cascade {
        scales: units(2),
        phase: 0.0,
        calls: Cell::new(0),
        tangent: 0,
    };
    let mut base = source.base();
    for row in &mut base.coeffs {
        row.truncate(2);
    }
    assert!(
        matches!(
            PreparedPeriodicLinearization::new(&source, &base, controls(1), &units(2), || false),
            Err(SidebandError::BaseUnresolved { .. })
        ),
        "retained base equations alone would miss generated second harmonic"
    );
    let alias = Parametric {
        modulation: 0.2,
        damping: 0.2,
        modulation_order: 63.0,
    };
    assert!(matches!(
        PreparedPeriodicLinearization::new(
            &alias,
            &orbit(2.0, vec![vec![C64::ZERO]; 2]),
            controls(2),
            &units(2),
            || false
        ),
        Err(SidebandError::SamplingUnresolved { .. })
    ));
    let lossless = Parametric {
        modulation: 0.0,
        damping: 0.0,
        modulation_order: 1.0,
    };
    let prepared = PreparedPeriodicLinearization::new(
        &lossless,
        &orbit(2.0, vec![vec![C64::ZERO]; 2]),
        controls(1),
        &units(2),
        || false,
    )
    .unwrap();
    assert!(matches!(
        prepared.factor(1.0, || false),
        Err(SidebandError::Singular { .. })
    ));
}

#[test]
fn exact_tangent_resources_and_cancellation_are_admitted_without_partial_results() {
    for tangent in [1, 2, 3] {
        let source = Cascade {
            scales: units(2),
            phase: 0.0,
            calls: Cell::new(0),
            tangent,
        };
        let result = PreparedPeriodicLinearization::new(
            &source,
            &source.base(),
            controls(2),
            &units(2),
            || false,
        );
        if tangent == 1 {
            assert!(matches!(result, Err(SidebandError::TangentUnavailable)));
        } else {
            assert!(matches!(
                result,
                Err(SidebandError::TangentRefused { sample: 0, .. })
            ));
        }
    }
    let source = Cascade {
        scales: units(2),
        phase: 0.0,
        calls: Cell::new(0),
        tangent: 0,
    };
    for resource in 0..3 {
        let mut c = controls(2);
        match resource {
            0 => c.maximum_unknowns = 1,
            1 => c.maximum_payload_bytes = 1,
            _ => c.maximum_work = 1,
        }
        assert!(matches!(
            PreparedPeriodicLinearization::new(&source, &source.base(), c, &units(2), || false),
            Err(SidebandError::Budget { .. })
        ));
        assert_eq!(source.calls.get(), 0);
    }
    assert!(matches!(
        PreparedPeriodicLinearization::new(&source, &source.base(), controls(2), &units(2), || {
            source.calls.get() >= 3
        }),
        Err(SidebandError::Cancelled)
    ));
    let prepared =
        PreparedPeriodicLinearization::new(&source, &source.base(), controls(2), &units(2), || {
            false
        })
        .unwrap();
    let factor = prepared.factor(0.6, || false).unwrap();
    let mut forcing = vec![C64::ZERO; 10];
    forcing[4] = C64::ONE;
    let a = factor.solve(&forcing, || false).unwrap();
    let mut polls = 0;
    assert!(matches!(
        factor.solve(&forcing, || {
            polls += 1;
            polls > 4
        }),
        Err(SidebandError::Cancelled)
    ));
    let b = factor.solve(&forcing, || false).unwrap();
    assert_eq!(a.coefficients, b.coefficients);
    assert_eq!(a.relative_residual.to_bits(), b.relative_residual.to_bits());
    assert_eq!(forcing[4], C64::ONE);
}
