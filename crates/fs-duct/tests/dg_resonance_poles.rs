//! G0/G1/G3: free resonance and decay from original spatial/material operators.
use fs_duct::dg::continuous::poles::{DgPoleError, PoleOptions};
use fs_duct::dg::{DgBoundary, DgDuct, DgMesh, PortId};
use fs_duct::{Duct, Segment};
use fs_material::gas::{GasSpec, GasState};
use fs_math::c64::C64;
use fs_phs::WallPin;

const WORK: usize = 1_000_000_000;
fn air() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap()
}
fn cylinder(
    elements: usize,
    order: usize,
    radius: f64,
    length: f64,
    inlet: DgBoundary,
    outlet: DgBoundary,
) -> DgDuct {
    DgDuct::new(
        &Duct {
            segments: vec![Segment::Cylinder { radius, length }],
        },
        &air(),
        DgMesh { elements, order },
        inlet,
        outlet,
    )
    .unwrap()
}
fn options(scale: f64) -> PoleOptions {
    PoleOptions {
        spectral_scale: scale,
        step_tolerance: 1e-11,
        residual_tolerance: 1e-10,
        max_iterations: 24,
        max_dimension: 1024,
        max_dense_bytes: 128 * 1024 * 1024,
        max_work: 2_000_000_000,
    }
}
fn seed(duct: &DgDuct, omega: f64) -> Vec<C64> {
    duct.harmonic_response(omega, &[(PortId(7), C64::ONE)], WORK, 1e-10, &mut || false)
        .unwrap()
        .energy_coordinates
}
fn near(actual: f64, expected: f64, relative: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= relative * scale.max(expected.abs()),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn g1_lossless_pole_and_complete_standing_wave_converge_to_cylinder() {
    let (radius, length) = (0.01, 0.5);
    let gas = air();
    let omega = core::f64::consts::PI * gas.sound_speed / (2.0 * length);
    let mut errors = Vec::new();
    for elements in [1, 2, 4] {
        let duct = cylinder(
            elements,
            2,
            radius,
            length,
            DgBoundary::Port(PortId(7)),
            DgBoundary::IdealOpen,
        );
        let guess = seed(&duct, 0.98 * omega);
        let pole = duct
            .resonance_pole(C64::new(0.0, -omega), &guess, &options(omega), || false)
            .unwrap();
        assert_eq!(pole.energy_coordinates.len(), duct.total_dofs());
        assert_eq!(pole.ports[0].id, PortId(7));
        assert_eq!(pole.wall_dissipation_w, 0.0);
        assert_eq!(pole.molecular_dissipation_w, 0.0);
        assert_eq!(pole.matched_absorption_w, 0.0);
        assert_eq!(pole.auxiliary_stored_energy_j, 0.0);
        assert!(pole.amplitude_decay_rate_per_s.abs() < 1e-10 * omega);
        near(pole.stored_energy_j, 0.5, 2e-14, 1.0);
        assert!(pole.relative_residual < 1e-10);
        let k = core::f64::consts::PI / (2.0 * length);
        let z = gas.density * gas.sound_speed / (core::f64::consts::PI * radius * radius);
        let amplitude = pole.pressure_pa[0];
        let mut field_error = 0.0_f64;
        for (i, x) in duct.node_positions().into_iter().enumerate() {
            field_error = field_error
                .max((pole.pressure_pa[i] / amplitude - C64::from_re((k * x).cos())).abs());
            field_error = field_error.max(
                (pole.volume_flow_m3_s[i].scale(z) / amplitude - C64::new(0.0, (k * x).sin()))
                    .abs(),
            );
        }
        errors.push((
            (pole.oscillation_frequency_rad_s / omega - 1.0).abs(),
            field_error,
        ));
    }
    assert!(
        errors[1].0 < errors[0].0 * 0.2 && errors[2].0 < errors[1].0 * 0.2,
        "{errors:?}"
    );
    assert!(errors[2].0 < 1e-5 && errors[2].1 < 0.002, "{errors:?}");
}

#[test]
fn g1_original_compliant_wall_gives_analytic_resonance_decay_and_heat() {
    let (radius, length) = (0.01, 0.3);
    let gas = air();
    let bulk = gas.density * gas.sound_speed.powi(2);
    let beta = 2.0 * bulk / radius;
    for resistance in [0.0, 300.0] {
        let pin = WallPin {
            surface_density: 0.2,
            stiffness_per_area: 1e7,
            resistance,
        };
        let duct = cylinder(1, 2, radius, length, DgBoundary::Closed, DgBoundary::Closed)
            .with_compliant_wall(pin)
            .unwrap();
        let decay = resistance / (2.0 * pin.surface_density);
        let omega = ((pin.stiffness_per_area + beta) / pin.surface_density - decay * decay).sqrt();
        let expected = C64::new(-decay, -omega);
        let coupling = (beta / pin.surface_density).sqrt();
        let wall_omega = (pin.stiffness_per_area / pin.surface_density).sqrt();
        let m = duct.field_dofs();
        let base = duct.to_energy_coords(&vec![1.0; m], &vec![0.0; m]);
        let mut guess: Vec<_> = base.into_iter().map(C64::from_re).collect();
        for i in 0..m {
            guess[2 * m + i] = guess[i].scale(-wall_omega / coupling);
            guess[3 * m + i] = -expected * guess[i].scale(1.0 / coupling);
        }
        let pole = duct
            .resonance_pole(expected.scale(0.99), &guess, &options(omega), || false)
            .unwrap();
        assert!((pole.laplace_pole_per_s - expected).abs() < 1e-9 * omega);
        let p = pole.pressure_pa[0];
        for &other in &pole.pressure_pa {
            assert!((other / p - C64::ONE).abs() < 1e-9);
        }
        let displacement = p.scale(-1.0 / beta);
        let velocity = expected * displacement;
        let surface = core::f64::consts::TAU * radius * length;
        let wall_energy = 0.5
            * surface
            * (pin.surface_density * velocity.norm_sq()
                + pin.stiffness_per_area * displacement.norm_sq());
        let heat = surface * resistance * velocity.norm_sq();
        near(pole.auxiliary_stored_energy_j, wall_energy, 1e-9, 1.0);
        near(pole.wall_dissipation_w, heat, 1e-9, 1.0);
        near(
            pole.decay_power_residual_w,
            0.0,
            1e-9,
            omega * pole.stored_energy_j,
        );
        near(pole.amplitude_decay_rate_per_s, decay, 1e-9, omega);
        assert_eq!(pole.molecular_dissipation_w, 0.0);
        assert_eq!(pole.matched_absorption_w, 0.0);
    }
}

// Independent uniform telegraph dispersion from the ORIGINAL scalar material
// branches: (s+thermal(s))*(s+drag+viscous(s)) + (c*k)^2 = 0.
// This uses neither the DG generator nor its vector pole solver.
fn wall_dispersion(model: &fs_phs::ZkDiffusive, omega: f64) -> C64 {
    let mut s = C64::new(-10.0, -omega);
    for _ in 0..40 {
        let mut p = s;
        let mut q = s + C64::from_re(model.series_drag_per_s());
        let (mut dp, mut dq) = (C64::ONE, C64::ONE);
        for (terms, value, derivative) in [
            (model.shunt_terms(), &mut p, &mut dp),
            (model.series_terms(), &mut q, &mut dq),
        ] {
            for term in terms {
                let denominator = s + C64::from_re(term.rate_per_s);
                *value = *value + (s / denominator).scale(term.gain_per_s);
                *derivative = *derivative
                    + C64::from_re(term.gain_per_s * term.rate_per_s) / (denominator * denominator);
            }
        }
        let step = (p * q + C64::from_re(omega * omega)) / (dp * q + p * dq);
        s = s - step;
        if step.abs() < 1e-13 * omega {
            break;
        }
    }
    s
}

#[test]
fn g1_viscothermal_poles_retain_all_memories_and_converge_to_telegraph_dispersion() {
    let (radius, length) = (0.001, 0.3);
    let omega = core::f64::consts::PI * air().sound_speed / (2.0 * length);
    let mut errors = Vec::new();
    for elements in [1, 2, 4] {
        let duct = cylinder(
            elements,
            2,
            radius,
            length,
            DgBoundary::Port(PortId(7)),
            DgBoundary::IdealOpen,
        )
        .with_viscothermal(4, [80.0, 1200.0], 0.2)
        .unwrap();
        let expected = wall_dispersion(duct.viscothermal_model().unwrap(), omega);
        let guess = seed(&duct, -expected.im * 0.98);
        let pole = duct
            .resonance_pole(expected, &guess, &options(omega), || false)
            .unwrap();
        errors.push((pole.laplace_pole_per_s - expected).abs() / omega);
        assert_eq!(pole.energy_coordinates.len(), duct.total_dofs());
        assert!(pole.energy_coordinates.len() > 2 * duct.field_dofs());
        assert!(pole.auxiliary_stored_energy_j > 0.0);
        assert!(pole.amplitude_decay_rate_per_s > 0.0);
        assert!(pole.oscillation_frequency_rad_s < omega);
        assert!(pole.wall_dissipation_w > 0.0);
        near(
            pole.decay_power_residual_w,
            0.0,
            1e-9,
            pole.wall_dissipation_w,
        );
        near(
            pole.wall_dissipation_w,
            2.0 * pole.amplitude_decay_rate_per_s * pole.stored_energy_j,
            1e-9,
            1.0,
        );
    }
    assert!(
        errors[1] < errors[0] * 0.2 && errors[2] < errors[1] * 0.2,
        "{errors:?}"
    );
    assert!(errors[2] < 1e-5, "{errors:?}");
}

#[test]
fn g1_real_material_relaxation_has_zero_frequency_and_original_positive_loss() {
    let duct = cylinder(1, 1, 0.001, 0.05, DgBoundary::Closed, DgBoundary::Closed)
        .with_viscothermal(1, [10.0, 1000.0], 0.95)
        .unwrap();
    let model = duct.viscothermal_model().unwrap();
    assert_eq!(model.shunt_terms().len(), 1);
    let term = model.shunt_terms()[0];
    let decay = term.gain_per_s + term.rate_per_s;
    let m = duct.field_dofs();
    let mut guess: Vec<_> = duct
        .to_energy_coords(&vec![1.0; m], &vec![0.0; m])
        .into_iter()
        .map(C64::from_re)
        .collect();
    let offset = 2 * m + model.series_terms().len() * m;
    for i in 0..m {
        guess[offset + i] = guess[i].scale(-(term.rate_per_s / term.gain_per_s).sqrt());
    }
    let pole = duct
        .resonance_pole(C64::from_re(-0.99 * decay), &guess, &options(decay), || {
            false
        })
        .unwrap();
    assert_eq!(pole.oscillation_frequency_rad_s, 0.0);
    near(pole.amplitude_decay_rate_per_s, decay, 1e-10, 1.0);
    assert!(pole.wall_dissipation_w > 0.0);
    near(
        pole.decay_power_residual_w,
        0.0,
        1e-10,
        pole.wall_dissipation_w,
    );
}

#[test]
fn g0_g4_band_capacity_seed_and_cancellation_refuse_without_changing_source() {
    let duct = cylinder(
        1,
        2,
        0.001,
        0.3,
        DgBoundary::Port(PortId(7)),
        DgBoundary::IdealOpen,
    )
    .with_viscothermal(4, [80.0, 1200.0], 0.2)
    .unwrap();
    let omega = 1000.0;
    let guess = seed(&duct, omega);
    let before = duct.continuous_operator().generator_per_s().to_vec();
    assert!(matches!(
        duct.resonance_pole(C64::new(0.0, -omega), &guess, &options(omega), || true),
        Err(DgPoleError::Source(_))
    ));
    let mut limit = options(omega);
    limit.max_dimension = duct.total_dofs() - 1;
    assert!(
        duct.resonance_pole(C64::new(0.0, -omega), &guess, &limit, || false)
            .is_err()
    );
    limit = options(omega);
    limit.max_dense_bytes = 1;
    assert!(
        duct.resonance_pole(C64::new(0.0, -omega), &guess, &limit, || false)
            .is_err()
    );
    limit = options(omega);
    limit.max_work = 1;
    assert!(
        duct.resonance_pole(C64::new(0.0, -omega), &guess, &limit, || false)
            .is_err()
    );
    assert!(
        duct.resonance_pole(C64::new(0.0, -1.0), &guess, &options(omega), || false)
            .is_err()
    );
    assert!(
        duct.resonance_pole(C64::new(f64::NAN, -omega), &guess, &options(omega), || {
            false
        })
        .is_err()
    );
    assert!(
        duct.resonance_pole(C64::new(0.0, -omega), &guess[..1], &options(omega), || {
            false
        })
        .is_err()
    );
    assert_eq!(before, duct.continuous_operator().generator_per_s());
}
