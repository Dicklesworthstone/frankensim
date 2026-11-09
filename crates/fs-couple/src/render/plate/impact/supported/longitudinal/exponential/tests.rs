use super::super::relaxation::SupportedRelaxationArm;
use super::*;
use crate::modal_acoustic_time::ModalAcousticState;
use crate::render::plate::impact::{
    ImpactConfig, VolumeSpring, exponential::ExponentialImpactSystem, relaxation::InitialMemory,
};
use fs_exec::CancelGate;

fn support(ends: StringEndMotion) -> SupportedStrings {
    let wire = WireSpan {
        endpoints_m: [[0.0, 0.0], [0.2, 0.0]],
        linear_density_kg_m: 0.012,
        tension_n: 8.0,
        bending_stiffness_n_m2: 1e-7,
        damping_per_s: vec![0.12],
    };
    let geometry = StringSupportGeometry::try_new(
        [[0.0; 3], [0.2, 0.0, 0.0]],
        [0.0, 1.0, 0.0],
        [0.6, 0.64, 0.48],
    )
    .unwrap();
    SupportedStrings::with_longitudinal_motion(
        &[wire],
        &[Some(StringStretching {
            axial_rigidity_n: 100.0,
            maximum_slope: 0.2,
        })],
        &[ends],
        &[geometry],
        &[LongitudinalSupportModes {
            solid_linear_density_kg_m: 0.01,
            modes: vec![1, 2],
            maximum_strain: 0.02,
        }],
        TranslatingSupport {
            mass_kg: 0.01,
            stiffness_n_m: 12.0,
            damping_n_s_m: 0.02,
            initial_position_m: 2e-5,
            initial_velocity_m_s: 0.0,
            maximum_travel_m: 0.01,
            maximum_slope: 0.2,
        },
    )
    .unwrap()
}
fn config(dt_s: f64) -> ImpactConfig {
    ImpactConfig {
        dt_s,
        max_steps: 512,
        maximum_energy_j: 1.0,
        energy_absolute_tolerance_j: 1e-12,
        energy_relative_tolerance: 1e-8,
        maximum_generalized_force: 1e3,
    }
}
fn initial() -> Vec<ModalAcousticState> {
    [1e-5, -7e-6, 2e-8, -1e-8]
        .into_iter()
        .map(|q| ModalAcousticState {
            displacement_m_sqrt_kg: q,
            velocity_m_sqrt_kg_per_s: 0.0,
        })
        .collect()
}
fn arm() -> SupportedRelaxationArm {
    SupportedRelaxationArm {
        string: 0,
        axial_rigidity_n: 50.0,
        relaxation_time_s: 0.001,
    }
}

// Independent two-word polynomial/automatic derivatives for this physical
// one-sine-per-plane, two-axial-mode fixture. Neither the production Fourier
// strain, support derivative, quadratic matrix nor work helper is called.
#[derive(Clone, Copy)]
struct Polynomial {
    value: Dd,
    derivative: [Dd; 5],
}
impl Polynomial {
    fn constant(v: f64) -> Self {
        Self {
            value: Dd::from_f64(v),
            derivative: [Dd::from_f64(0.0); 5],
        }
    }
    fn variable(v: f64, i: usize) -> Self {
        let mut result = Self::constant(v);
        result.derivative[i] = Dd::from_f64(1.0);
        result
    }
}
impl core::ops::Add for Polynomial {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self {
            value: self.value + rhs.value,
            derivative: core::array::from_fn(|i| self.derivative[i] + rhs.derivative[i]),
        }
    }
}
impl core::ops::Sub for Polynomial {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self {
            value: self.value - rhs.value,
            derivative: core::array::from_fn(|i| self.derivative[i] - rhs.derivative[i]),
        }
    }
}
impl core::ops::Mul for Polynomial {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        Self {
            value: self.value * rhs.value,
            derivative: core::array::from_fn(|i| {
                self.derivative[i] * rhs.value + self.value * rhs.derivative[i]
            }),
        }
    }
}
fn polynomial(s: &SupportedStrings, z: &[f64], nonlinear_only: bool) -> Polynomial {
    let c = Polynomial::constant;
    let SupportedPotential::Longitudinal(p) = &s.0.potentials[0] else {
        unreachable!()
    };
    let x = Polynomial::variable(z[4], 4) * c(s.0.weight);
    let q: Vec<_> = (0..4)
        .map(|i| Polynomial::variable(z[i], i) - x * c(s.0.translation[i]))
        .collect();
    let y0 = x * c(p.chord[0]);
    let y1 = q[0] * c(p.scales[0]);
    let z0 = x * c(p.chord[1]);
    let z1 = q[1] * c(p.scales[0]);
    let linear = [
        x * c(p.chord[2]),
        q[2] * c(p.model.axial_slope_scales()[0]),
        q[3] * c(p.model.axial_slope_scales()[1]),
    ];
    let geometric = [
        c(0.5) * (y0 * y0 + z0 * z0) + c(0.25) * (y1 * y1 + z1 * z1),
        y0 * y1 + z0 * z1,
        c(0.25) * (y1 * y1 + z1 * z1),
    ];
    let coefficient = p.model.strain_energy_coefficient();
    let mut energy = c(0.0);
    for i in 0..3 {
        let k = c(coefficient * if i == 0 { 2.0 } else { 1.0 });
        energy = energy + k * (linear[i] * geometric[i] + c(0.5) * geometric[i] * geometric[i]);
        if !nonlinear_only {
            energy = energy + c(0.5) * k * linear[i] * linear[i];
        }
    }
    if !nonlinear_only {
        energy =
            energy + c(0.5 * s.0.support.stiffness_n_m) * x * x + c(0.5 * p.lift_stiffness) * x * x;
        for i in 0..2 {
            energy = energy + c(0.5 * p.omegas[i] * p.omegas[i]) * q[i] * q[i];
        }
    }
    energy
}

#[test]
fn original_signed_support_q_and_full_strain_work_match_independent_polynomials() {
    for ends in [
        StringEndMotion::First,
        StringEndMotion::Second,
        StringEndMotion::Both,
    ] {
        let s = support(ends);
        let split = s.exponential_split().unwrap();
        assert_eq!(split.mode_count(), 5);
        assert!(
            split
                .stiffness()
                .iter()
                .enumerate()
                .any(|(i, k)| i / 5 != i % 5 && *k != 0.0)
        );
        for j in 0..5 {
            let mut direction = [0.0; 5];
            direction[j] = 1.0;
            let mut tangent = [0.0; 5];
            s.hessian_vector(&[0.0; 5], &direction, &mut tangent);
            for i in 0..5 {
                let k = split.stiffness()[5 * i + j];
                assert!((k - tangent[i]).abs() < 3e-12 * k.abs().max(tangent[i].abs()).max(1.0));
            }
        }
        let a = [3e-6, -2e-6, 1e-8, -5e-9, 2e-6];
        let reference = polynomial(&s, &a, false);
        assert!(
            (s.potential(&a) - reference.value.to_f64()).abs() < 1e-13 * reference.value.to_f64()
        );
        let mut gradient = [0.0; 5];
        s.gradient(&a, &mut gradient);
        for i in 0..5 {
            let expected = reference.derivative[i].to_f64();
            assert!((gradient[i] - expected).abs() < 2e-12 * expected.abs().max(1e-9));
        }
        let direction = [1e-6, 2e-6, -1e-8, 2e-8, -5e-7];
        for exponent in [0, 8, 16, 22] {
            let step = 2.0_f64.powi(-exponent);
            let b: Vec<_> = a.iter().zip(direction).map(|(a, d)| a + step * d).collect();
            let mid: Vec<_> = a.iter().zip(&b).map(|(a, b)| (a + b) * 0.5).collect();
            let full = |q: &[f64]| q.iter().flat_map(|q| [*q, 0.0]).collect::<Vec<_>>();
            let actual = split
                .midpoint_nonlinear_remainder(&full(&a), &full(&b), &full(&mid))
                .unwrap();
            let mut expected = polynomial(&s, &b, true).value - polynomial(&s, &a, true).value;
            let g = polynomial(&s, &mid, true).derivative;
            for i in 0..5 {
                expected = expected - g[i] * (Dd::from_f64(b[i]) - Dd::from_f64(a[i]));
            }
            let expected = expected.to_f64();
            assert!(
                (actual - expected).abs() < 2e-7 * expected.abs() + 2e-36,
                "{ends:?}, 2^-{exponent}: {actual:e} != {expected:e}"
            );
        }
    }
}

#[test]
fn numerical_envelope_keeps_original_storage_damping_and_every_axial_history() {
    let original = support(StringEndMotion::Second)
        .with_axial_viscosity(vec![0.003])
        .unwrap()
        .with_longitudinal_relaxation(
            vec![
                arm(),
                SupportedRelaxationArm {
                    relaxation_time_s: 0.003,
                    ..arm()
                },
            ],
            InitialMemory::Relaxed,
            6,
        )
        .unwrap();
    let prepared = original.clone().with_exponential_envelope().unwrap();
    assert!(!original.has_exponential_envelope());
    assert!(prepared.has_exponential_envelope());
    assert_eq!(prepared.relaxation_state_count(), 6);
    assert_eq!(prepared.mode_count(), original.mode_count());
    assert!(prepared.clock_bounds()[0] > original.clock_bounds()[0]);
    assert_eq!(
        prepared.clock_bounds(),
        prepared
            .clone()
            .with_exponential_envelope()
            .unwrap()
            .clock_bounds()
    );
    let x = [3e-6, -2e-6, 1e-8, -5e-9, 2e-6];
    assert_eq!(
        prepared.potential(&x).to_bits(),
        original.potential(&x).to_bits()
    );
    let (mut a, mut b) = ([0.0; 5], [0.0; 5]);
    original.gradient(&x, &mut a);
    prepared.gradient(&x, &mut b);
    assert_eq!(a, b);
    let (mut a, mut b) = (vec![0.0; 100], vec![0.0; 100]);
    original.add_resistance(0, 10, &mut a).unwrap();
    prepared.add_resistance(0, 10, &mut b).unwrap();
    assert_eq!(a, b);
    assert!(
        a.iter()
            .enumerate()
            .any(|(i, v)| i / 10 != i % 10 && *v != 0.0)
    );
    let branches = prepared.memory_branches(0, 10).unwrap();
    assert_eq!(branches.len(), 6);
    assert!(
        branches
            .iter()
            .all(|b| b.strain.zero_constant_linear_strain().is_some())
    );
    assert!(
        support(StringEndMotion::Second)
            .with_finite_strain_axial_viscosity(vec![0.003])
            .unwrap()
            .with_exponential_envelope()
            .is_err()
    );
}

#[test]
fn complete_supported_material_motion_and_heat_converge_without_separate_history_steps() {
    let gate = CancelGate::new_clock_free();
    let loaded = support(StringEndMotion::Second)
        .with_axial_viscosity(vec![0.003])
        .unwrap()
        .with_longitudinal_relaxation(vec![arm()], InitialMemory::Relaxed, 3)
        .unwrap();
    let run = |dt: f64, steps: usize| {
        let mut run = ExponentialImpactSystem::new(
            vec![loaded.body_with_relative_states(&initial()).unwrap()],
            vec![],
            vec![],
            config(dt),
            &gate,
        )
        .unwrap();
        assert_eq!(run.state().len(), 13);
        assert_eq!(run.relaxation_branch_count(), 3);
        let before = run.stored_energy_j();
        let (mut heat, mut work) = (0.0, 0.0);
        let mut force = [0.0; 5];
        force[4] = 0.02 * loaded.force_weight();
        for _ in 0..steps {
            let frame = run.step(&force, &gate).unwrap();
            assert!(frame.dissipated_energy_j >= 0.0);
            heat += frame.dissipated_energy_j;
            work += frame.supplied_work_j;
        }
        assert!(heat > 0.0);
        assert!((run.stored_energy_j() + heat - before - work).abs() < 1e-10);
        (run.state().to_vec(), heat, work)
    };
    let reference = run(2.5e-7, 256);
    let coarse = run(2e-6, 32);
    let fine = run(1e-6, 64);
    let error = |result: &(Vec<f64>, f64, f64)| {
        result
            .0
            .iter()
            .zip(&reference.0)
            .map(|(a, b)| (a - b).abs())
            .chain([
                (result.1 - reference.1).abs(),
                (result.2 - reference.2).abs(),
            ])
            .fold(0.0_f64, f64::max)
    };
    assert!(error(&fine) < 0.4 * error(&coarse) + 1e-12);
    assert!(error(&fine) < 2e-7);
}

#[test]
fn complete_source_capacity_lossless_envelope_and_atomic_refusal_remain_binding() {
    let gate = CancelGate::new_clock_free();
    let dt = 1e-6;
    let base = support(StringEndMotion::Second);
    let mut run = ExponentialImpactSystem::new(
        vec![base.body_with_relative_states(&initial()).unwrap()],
        vec![],
        vec![],
        config(dt),
        &gate,
    )
    .unwrap();
    assert_eq!(run.state().len(), 10);
    assert_eq!(run.relaxation_branch_count(), 0);
    assert!(run.relaxation.is_some()); // envelope metadata, zero extra states.
    let state = run.state().to_vec();
    let energy = run.stored_energy_j();
    assert!(run.step(&[2e3; 5], &gate).is_err());
    assert_eq!(run.state(), state);
    assert_eq!(run.stored_energy_j().to_bits(), energy.to_bits());
    assert_eq!(run.samples(), 0);
    let full = base
        .clone()
        .with_longitudinal_relaxation(vec![arm(); 40], InitialMemory::Relaxed, 120)
        .unwrap();
    assert_eq!(full.relaxation_state_count(), 120);
    assert!(
        ExponentialImpactSystem::new(vec![full.body()], vec![], vec![], config(dt), &gate).is_err()
    );
    let upper = base
        .clone()
        .with_exponential_envelope()
        .unwrap()
        .clock_bounds()[0]
        * 1.01;
    let sourced = base
        .with_longitudinal_relaxation(vec![], InitialMemory::Relaxed, 1)
        .unwrap()
        .bind_longitudinal_relaxation_domain(upper, dt)
        .unwrap();
    assert!(
        ExponentialImpactSystem::new(
            vec![sourced.body()],
            vec![],
            vec![],
            config(dt * 0.5),
            &gate
        )
        .is_err()
    );
    let volume = VolumeSpring {
        bulk_modulus_pa: 4.0 * upper * upper,
        volume_m3: 1.0,
        areas: vec![0.0, 0.0, 0.0, 0.0, 1.0],
    };
    assert!(
        ExponentialImpactSystem::new(
            vec![sourced.body()],
            vec![],
            vec![volume],
            config(dt),
            &gate
        )
        .is_err()
    );
    let valid =
        ExponentialImpactSystem::new(vec![sourced.body()], vec![], vec![], config(dt), &gate)
            .unwrap();
    let radiation = crate::render::plate::impact::radiation::Model {
        ports: 1,
        poles: vec![crate::render::plate::impact::radiation::Pole {
            omega: 2.0 * upper,
            zeta: 0.2,
            coupling: vec![1.0],
        }],
    };
    assert!(
        valid
            .with_radiation_load(&radiation, &[4], 1, &gate)
            .is_err()
    );
}
