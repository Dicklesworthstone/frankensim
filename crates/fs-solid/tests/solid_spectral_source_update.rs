//! Continuous source clocks survive physical stiffness/viscosity replacement.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_modal::time::{ModalToneBudget, PreparedModalTones};
use fs_solid::dynamics::spectral::{
    CompleteSpectralSolidControls, CompleteSpectralSolidError, DampedSpectralSolid,
    HarmonicDampedSpectralSolid, MultiToneDampedSpectralSolid,
};
use fs_solid::linear3::TetElasticAssembly;
use fs_sparse::{Coo, Csr};

fn with_cx(cancelled: bool, run: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    if cancelled {
        gate.request();
    }
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        run(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 65,
                kernel_id: 35,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ));
    });
}
fn controls() -> CompleteSpectralSolidControls {
    CompleteSpectralSolidControls {
        maximum_free_dofs: 16,
        maximum_preparation_bytes: 8 * 1024 * 1024,
        maximum_preparation_multiply_terms: 100_000,
        maximum_step_multiply_terms: 1_000_000,
        maximum_relative_operator_error: 1e-12,
        maximum_absolute_eigen_residual: 0.0,
        maximum_relative_state_error: 1e-12,
        maximum_work_error_j: 1e-9,
        maximum_steps: 1000,
    }
}
fn matrix(value: f64) -> Csr {
    let mut coo = Coo::new(1, 1);
    coo.push(0, 0, value);
    coo.assemble()
}
fn assembly(k: f64) -> TetElasticAssembly {
    TetElasticAssembly {
        stiffness: matrix(k),
        mass: matrix(4.0),
        free_dofs: vec![0],
        total_mass_kg: 4.0,
        minimum_scaled_jacobian: 1.0,
        element_volumes_m3: vec![1.0],
    }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e} +/- {tolerance:e}"
    );
}
fn harmonic_bits(s: &HarmonicDampedSpectralSolid) -> Vec<u64> {
    [
        s.displacement_m()[0],
        s.velocity_mps()[0],
        s.time_s(),
        s.stored_energy_j(),
        s.kinetic_energy_j(),
        s.cumulative_work_j(),
        s.cumulative_dissipated_j(),
        s.cumulative_parameter_work_j(),
        s.damping_rates_per_s()[0],
    ]
    .into_iter()
    .map(f64::to_bits)
    .chain([s.steps()])
    .collect()
}
fn solid_bits(s: &DampedSpectralSolid) -> Vec<u64> {
    [
        s.displacement_m()[0],
        s.velocity_mps()[0],
        s.time_s(),
        s.stored_energy_j(),
        s.kinetic_energy_j(),
        s.cumulative_work_j(),
        s.cumulative_dissipated_j(),
        s.cumulative_parameter_work_j(),
        s.damping_rates_per_s()[0],
    ]
    .into_iter()
    .map(f64::to_bits)
    .chain([s.steps()])
    .collect()
}
fn simpson(mut value: impl FnMut(f64) -> f64, h: f64) -> f64 {
    let mut sum = 0.0;
    for i in 0..=256 {
        let weight = if i == 0 || i == 256 {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        sum += weight * value(h * i as f64 / 256.0);
    }
    sum * h / 768.0
}

#[test]
fn retuning_into_resonance_keeps_phase_heat_and_clock_and_late_cache_refusal_is_atomic() {
    with_cx(false, |cx| {
        let mut solid = HarmonicDampedSpectralSolid::prepare_with_basis(
            &assembly(16.0),
            &matrix(0.8),
            &[4.0],
            &[0.5],
            3.0,
            0.2,
            &[0.1],
            &[-0.05],
            controls(),
            cx,
        )
        .unwrap();
        let (mut cosine, mut sine) = ([2.0], [-0.3]);
        for _ in 0..3 {
            solid
                .advance_continuous(&mut cosine, &mut sine, cx)
                .unwrap();
        }
        let before = harmonic_bits(&solid);
        let phase = solid.advance_force_phase(cosine[0], sine[0]).unwrap();
        // Complete physical preparation passes this budget; its new harmonic
        // scaling needs more than the admitted 100000 preparation terms.
        assert!(matches!(
            solid.replace_parameters(
                &matrix(4.0 * 1024.0 * 1024.0),
                &matrix(0.0),
                &[1024.0 * 1024.0],
                &[-0.5],
                1e9,
                cx
            ),
            Err(CompleteSpectralSolidError::Budget { .. })
        ));
        assert_eq!(harmonic_bits(&solid), before);
        assert_eq!(
            solid.advance_force_phase(cosine[0], sine[0]).unwrap(),
            phase
        );
        with_cx(true, |cancelled| {
            assert!(matches!(
                solid.replace_parameters(
                    &matrix(36.0),
                    &matrix(0.0),
                    &[9.0],
                    &[-0.5],
                    10.0,
                    cancelled
                ),
                Err(CompleteSpectralSolidError::Cancelled)
            ))
        });
        assert_eq!(harmonic_bits(&solid), before);
        let (u, v, old_heat, old_work) = (
            solid.displacement_m()[0],
            solid.velocity_mps()[0],
            solid.cumulative_dissipated_j(),
            solid.cumulative_work_j(),
        );
        let report = solid
            .replace_parameters(&matrix(36.0), &matrix(0.0), &[9.0], &[-0.5], 10.0, cx)
            .unwrap();
        near(report.parameter_work_j, 10.0 * u * u, 2e-15);
        assert_eq!(solid.displacement_m()[0].to_bits(), u.to_bits());
        assert_eq!(solid.velocity_mps()[0].to_bits(), v.to_bits());
        assert_eq!(
            solid.cumulative_dissipated_j().to_bits(),
            old_heat.to_bits()
        );
        assert_eq!(solid.cumulative_work_j().to_bits(), old_work.to_bits());
        assert_eq!(solid.steps(), 3);
        assert_eq!(
            solid.advance_force_phase(cosine[0], sine[0]).unwrap(),
            phase
        );
        let (fc, fs) = (cosine[0], sine[0]);
        // Independent undamped resonant convolution, with both quadratures.
        let state = |t: f64| {
            let (s, c) = (3.0 * t).sin_cos();
            let (a, b) = (fc / 4.0, fs / 4.0);
            (
                u * c + v * s / 3.0 + a * t * s / 6.0 + b * (s - 3.0 * t * c) / 18.0,
                -3.0 * u * s + v * c + a * (s / 6.0 + 0.5 * t * c) + 0.5 * b * t * s,
            )
        };
        let expected = state(0.2);
        let work = simpson(
            |t| {
                let (s, c) = (3.0 * t).sin_cos();
                (fc * c + fs * s) * state(t).1
            },
            0.2,
        );
        let frame = solid
            .advance_continuous(&mut cosine, &mut sine, cx)
            .unwrap();
        near(solid.displacement_m()[0], expected.0, 3e-14);
        near(solid.velocity_mps()[0], expected.1, 3e-14);
        near(frame.energy.supplied_work_j, work, 3e-13);
        assert_eq!(frame.energy.dissipated_j, 0.0);
        assert_eq!(
            frame.energy.cumulative_dissipated_j.to_bits(),
            old_heat.to_bits()
        );
        near(frame.energy.cumulative_balance_error_j, 0.0, 2e-13);
        assert_eq!((cosine[0], sine[0]), phase);
        assert_eq!(solid.steps(), 4);
        assert_eq!(solid.forcing_angular_frequency_rad_s(), 3.0);

        // DC keeps the lower-cost held-force path even after a basis update.
        let mut dc_controls = controls();
        dc_controls.maximum_preparation_multiply_terms = 20_000;
        let mut dc = HarmonicDampedSpectralSolid::prepare_with_basis(
            &assembly(16.0),
            &matrix(0.0),
            &[4.0],
            &[0.5],
            0.0,
            0.1,
            &[0.2],
            &[0.3],
            dc_controls,
            cx,
        )
        .unwrap();
        dc.replace_parameters(&matrix(0.0), &matrix(0.0), &[0.0], &[-0.5], 10.0, cx)
            .unwrap();
        let (mut c, mut s) = ([1.0], [27.0]);
        dc.advance_continuous(&mut c, &mut s, cx).unwrap();
        near(dc.displacement_m()[0], 0.2 + 0.03 + 0.01 / 8.0, 2e-15);
        near(dc.velocity_mps()[0], 0.3 + 0.1 / 4.0, 2e-15);
        assert_eq!(c, [1.0]);
        assert_eq!(s, [27.0]);
    });
}

// Exact affine-quadrature particular solution plus the independent homogeneous
// underdamped poles. This integrates physical F/M, without calling a production
// force transition or using endpoint energy to infer work/heat.
struct AffineOracle<'a> {
    frequencies: &'a [f64],
    c0: &'a [f64],
    s0: &'a [f64],
    c1: &'a [f64],
    s1: &'a [f64],
    initial: [f64; 2],
    h: f64,
    lambda: f64,
    rate: f64,
}
impl AffineOracle<'_> {
    fn state(&self, t: f64) -> (f64, f64) {
        let (mut q, mut v, mut q0, mut v0) = (0.0, 0.0, 0.0, 0.0);
        for j in 0..self.frequencies.len() {
            let w = self.frequencies[j];
            let (d, e) = (self.lambda - w * w, self.rate * w);
            let denom = d * d + e * e;
            let (fc, fs) = (self.c0[j] / 4.0, self.s0[j] / 4.0);
            let (dfc, dfs) = (
                (self.c1[j] - self.c0[j]) / (4.0 * self.h),
                (self.s1[j] - self.s0[j]) / (4.0 * self.h),
            );
            let (b, f) = ((d * dfc - e * dfs) / denom, (e * dfc + d * dfs) / denom);
            let (rc, rs) = (
                fc - self.rate * b - 2.0 * w * f,
                fs - self.rate * f + 2.0 * w * b,
            );
            let (a, g) = ((d * rc - e * rs) / denom, (e * rc + d * rs) / denom);
            let (sin, cos) = (w * t).sin_cos();
            q += (a + b * t) * cos + (g + f * t) * sin;
            v += (b + w * (g + f * t)) * cos + (f - w * (a + b * t)) * sin;
            q0 += a;
            v0 += b + w * g;
        }
        let (x, y) = (self.initial[0] - q0, self.initial[1] - v0);
        let decay = 0.5 * self.rate;
        let w = (self.lambda - decay * decay).sqrt();
        let (sin, cos) = (w * t).sin_cos();
        let exp = (-decay * t).exp();
        (
            q + exp * (x * cos + (y + decay * x) * sin / w),
            v + exp * (y * cos - (decay * y + self.lambda * x) * sin / w),
        )
    }
    fn force(&self, t: f64) -> f64 {
        self.frequencies
            .iter()
            .enumerate()
            .map(|(j, &w)| {
                let (s, c) = (w * t).sin_cos();
                (self.c0[j] + (self.c1[j] - self.c0[j]) * t / self.h) * c
                    + (self.s0[j] + (self.s1[j] - self.s0[j]) * t / self.h) * s
            })
            .sum()
    }
}

#[test]
fn affine_multitone_replacement_keeps_complete_interference_power_and_counts_both_banks() {
    with_cx(false, |cx| {
        let frequencies = [0.0, 2.0, 3.0];
        let mut solid = MultiToneDampedSpectralSolid::prepare_with_basis(
            &assembly(16.0),
            &matrix(1.6),
            &[4.0],
            &[0.5],
            &frequencies,
            0.11,
            &[0.2],
            &[-0.1],
            true,
            controls(),
            ModalToneBudget::realtime_reference(),
            cx,
        )
        .unwrap();
        let (mut c, mut s) = ([0.3, 0.7, -0.2], [0.0, -0.1, 0.5]);
        solid
            .advance_affine_continuous(&mut c, &mut s, &[0.5, 0.8, -0.1], &[0.0, 0.2, 0.4], cx)
            .unwrap();
        let before = solid_bits(solid.solid());
        let phases: Vec<_> = (0..3)
            .map(|i| solid.advance_force_phase(i, c[i], s[i]).unwrap())
            .collect();
        let mut too_little = ModalToneBudget::realtime_reference();
        too_little.maximum_preparation_work = 0;
        // Rejection is in the new native source bank, after physical K/C/state
        // preparation. The published owner must keep its old kernels as well.
        assert!(
            solid
                .replace_parameters(
                    &matrix(64.0),
                    &matrix(4.8),
                    &[16.0],
                    &[-0.5],
                    10.0,
                    too_little,
                    cx
                )
                .is_err()
        );
        assert_eq!(solid_bits(solid.solid()), before);
        for i in 0..3 {
            assert_eq!(solid.advance_force_phase(i, c[i], s[i]).unwrap(), phases[i]);
        }
        let initial = [
            solid.solid().displacement_m()[0],
            solid.solid().velocity_mps()[0],
        ];
        let (heat, work) = (
            solid.solid().cumulative_dissipated_j(),
            solid.solid().cumulative_work_j(),
        );
        let report = solid
            .replace_parameters(
                &matrix(64.0),
                &matrix(4.8),
                &[16.0],
                &[-0.5],
                10.0,
                ModalToneBudget::realtime_reference(),
                cx,
            )
            .unwrap();
        near(
            report.parameter_work_j,
            24.0 * initial[0] * initial[0],
            3e-15,
        );
        assert_eq!(
            solid.solid().displacement_m()[0].to_bits(),
            initial[0].to_bits()
        );
        assert_eq!(
            solid.solid().velocity_mps()[0].to_bits(),
            initial[1].to_bits()
        );
        assert_eq!(
            solid.solid().cumulative_dissipated_j().to_bits(),
            heat.to_bits()
        );
        assert_eq!(solid.solid().cumulative_work_j().to_bits(), work.to_bits());
        assert_eq!(solid.solid().steps(), 1);
        assert!(solid.is_affine());
        assert_eq!(solid.angular_frequencies_rad_s(), frequencies.as_slice());
        for i in 0..3 {
            assert_eq!(solid.advance_force_phase(i, c[i], s[i]).unwrap(), phases[i]);
        }
        let (c0, s0) = (c, s);
        let (ce, se) = (
            [c[0] + 0.05, c[1] - 0.04, c[2] + 0.02],
            [s[0], s[1] + 0.03, s[2] - 0.06],
        );
        let oracle = AffineOracle {
            frequencies: &frequencies,
            c0: &c0,
            s0: &s0,
            c1: &ce,
            s1: &se,
            initial,
            h: 0.11,
            lambda: 16.0,
            rate: 1.2,
        };
        let expected = oracle.state(0.11);
        let expected_work = simpson(|t| oracle.force(t) * oracle.state(t).1, 0.11);
        let expected_heat = simpson(|t| 4.8 * oracle.state(t).1.powi(2), 0.11);
        let frame = solid
            .advance_affine_continuous(&mut c, &mut s, &ce, &se, cx)
            .unwrap();
        near(solid.solid().displacement_m()[0], expected.0, 3e-14);
        near(solid.solid().velocity_mps()[0], expected.1, 3e-14);
        near(frame.energy.supplied_work_j, expected_work, 5e-13);
        near(frame.energy.dissipated_j, expected_heat, 5e-13);
        near(frame.energy.cumulative_balance_error_j, 0.0, 3e-13);
        for i in 0..3 {
            let (sin, cos) = (frequencies[i] * 0.11).sin_cos();
            near(c[i], ce[i] * cos + se[i] * sin, 2e-15);
            near(s[i], se[i] * cos - ce[i] * sin, 2e-15);
        }
        // The shared native envelope exactly admits initial construction.
        // Retaining that bank during a physical replacement needs more space.
        let bank = PreparedModalTones::try_new_affine_with_damping_rates(
            &[2.0],
            &[0.4],
            &frequencies,
            0.11,
            ModalToneBudget::realtime_reference(),
        )
        .unwrap();
        let mut tight = controls();
        tight.maximum_preparation_bytes =
            (96 + 4096 + 32 * 3 + bank.preparation_storage_bound_bytes()) as usize;
        let mut bounded = MultiToneDampedSpectralSolid::prepare_with_basis(
            &assembly(16.0),
            &matrix(1.6),
            &[4.0],
            &[0.5],
            &frequencies,
            0.11,
            &[0.2],
            &[-0.1],
            true,
            tight,
            ModalToneBudget::realtime_reference(),
            cx,
        )
        .unwrap();
        let before = solid_bits(bounded.solid());
        assert!(matches!(
            bounded.replace_parameters(
                &matrix(64.0),
                &matrix(4.8),
                &[16.0],
                &[-0.5],
                10.0,
                ModalToneBudget::realtime_reference(),
                cx
            ),
            Err(CompleteSpectralSolidError::Budget { .. })
        ));
        assert_eq!(solid_bits(bounded.solid()), before);
    });
}
