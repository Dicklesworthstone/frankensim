//! G0/G1/G3/G4: full two-fluid waves, physical ports, dissipation and rollback.

use fs_duct::superfluid::{
    BoundaryVelocities, ChannelError, EndVelocities, LinearizationLimits, MutualFriction,
    TwoFluidChannel,
};
use fs_material::superfluid::{TwoFluidProperties, TwoFluidState};

fn material(alpha: f64) -> TwoFluidState {
    // Manufactured constitutive data for independent analytic solutions.
    TwoFluidState::try_new(TwoFluidProperties {
        temperature_k: 2.0,
        density_kg_m3: 4.0,
        superfluid_fraction: 0.75,
        entropy_j_kg_k: 4.0,
        specific_heat_cp_j_kg_k: 16.0,
        adiabatic_sound_speed_m_s: 10.0,
        thermal_expansion_per_k: alpha,
        normal_viscosity_pa_s: 0.03,
    })
    .unwrap()
}

fn channel(n: usize, source: TwoFluidState) -> TwoFluidChannel {
    TwoFluidChannel::try_new(
        source,
        1.0,
        0.01,
        n,
        LinearizationLimits {
            max_pressure_perturbation_pa: 1.0,
            max_temperature_perturbation_k: 0.02,
            max_component_speed_m_s: 0.1,
        },
    )
    .unwrap()
}

fn standing(n: usize, source: TwoFluidState, p: f64, t: f64) -> TwoFluidChannel {
    let mut duct = channel(n, source);
    let shape: Vec<_> = (0..n)
        .map(|i| (std::f64::consts::PI * (i as f64 + 0.5) / n as f64).cos())
        .collect();
    duct.initialize(
        &shape.iter().map(|x| p * x).collect::<Vec<_>>(),
        &shape.iter().map(|x| t * x).collect::<Vec<_>>(),
        &vec![0.0; n + 1],
        &vec![0.0; n + 1],
        || false,
    )
    .unwrap();
    duct
}

fn advance_to(duct: &mut TwoFluidChannel, end: f64, friction: MutualFriction) {
    let steps = (end / (0.7 * duct.maximum_step_s())).ceil() as usize;
    let dt = end / steps as f64;
    for _ in 0..steps {
        let report = duct
            .step(dt, BoundaryVelocities::CLOSED, friction, || false)
            .unwrap();
        assert!(report.maximum_relative_linear_residual < 1e-13);
        assert!(
            report.wave_balance_defect_j.abs()
                < 1e-10 * (report.wave_energy_before_j + report.wave_energy_after_j)
        );
    }
}

fn relative_error(duct: &TwoFluidChannel, p: f64, t: f64, speed: f64, time: f64) -> f64 {
    let phase = std::f64::consts::PI * speed * time;
    let mut sum = 0.0;
    for i in 0..duct.cells() {
        let shape = (std::f64::consts::PI * (i as f64 + 0.5) / duct.cells() as f64).cos();
        let cell = duct.cell(i).unwrap();
        if p != 0.0 {
            sum += (cell.pressure_perturbation_pa / p - shape * phase.cos()).powi(2);
        }
        if t != 0.0 {
            sum += (cell.temperature_perturbation_k / t - shape * phase.cos()).powi(2);
        }
    }
    (sum / duct.cells() as f64).sqrt()
}

#[test]
fn first_and_second_sound_have_independent_speeds_and_second_order_convergence() {
    let source = material(0.0);
    for (p, t, speed) in [(0.01, 0.0, 10.0), (0.0, 0.001, 6.0_f64.sqrt())] {
        let time = 0.37 / speed;
        let mut errors = Vec::new();
        for n in [12, 24, 48] {
            let mut duct = standing(n, source, p, t);
            let energy = duct.wave_energy_j().unwrap();
            advance_to(&mut duct, time, MutualFriction::None);
            errors.push(relative_error(&duct, p, t, speed, time));
            assert!((duct.wave_energy_j().unwrap() / energy - 1.0).abs() < 2e-12);
            for i in 0..n {
                let cell = duct.cell(i).unwrap();
                if p == 0.0 {
                    assert_eq!(cell.pressure_perturbation_pa, 0.0);
                }
                if t == 0.0 {
                    assert_eq!(cell.temperature_perturbation_k, 0.0);
                }
            }
            let middle = duct.face(n / 2).unwrap();
            if p == 0.0 {
                assert_eq!(middle.mass_flux_kg_m2_s, 0.0);
                assert!(middle.normal_velocity_m_s * middle.superfluid_velocity_m_s < 0.0);
                assert!(middle.counterflow_heat_flux_w_m2 > 0.0);
            } else {
                assert_eq!(middle.normal_velocity_m_s, middle.superfluid_velocity_m_s);
                assert_eq!(middle.counterflow_heat_flux_w_m2, 0.0);
            }
        }
        assert!(errors[0] / errors[1] > 3.6, "{errors:?}");
        assert!(errors[1] / errors[2] > 3.6, "{errors:?}");
        assert!(errors[2] < 0.0007, "{errors:?}");
    }
}

#[test]
fn finite_expansion_propagates_coupled_eigenmodes_with_correct_polarization() {
    for alpha in [-0.3, 0.3] {
        let source = material(alpha);
        let modes = source.sound_modes();
        // Eigenvectors derived directly from the physical pressure equation:
        // (c_mode²-c_ad²)*p = rho*c_ad²*alpha*(fs/fn)*T*s²/cp * theta.
        let coupling = 4.0 * 100.0 * alpha * 3.0 * 2.0 * 16.0 / 16.0;
        for speed in [modes.fast_m_s, modes.slow_m_s] {
            let temperature = 0.0001;
            let pressure = coupling * temperature / (speed * speed - 100.0);
            let mut duct = standing(64, source, pressure, temperature);
            let initial = duct.wave_energy_j().unwrap();
            let time = 0.37 / speed;
            advance_to(&mut duct, time, MutualFriction::None);
            assert!(relative_error(&duct, pressure, temperature, speed, time) < 0.0007);
            assert!((duct.wave_energy_j().unwrap() / initial - 1.0).abs() < 3e-12);
            let face = duct.face(32).unwrap();
            let factor = (std::f64::consts::PI * speed * time).sin();
            let expected_bulk = pressure / (source.density_kg_m3() * speed) * factor;
            let expected_relative =
                source.entropy_j_kg_k() / source.normal_fraction() * temperature / speed * factor;
            assert!((face.barycentric_velocity_m_s / expected_bulk - 1.0).abs() < 0.001);
            assert!((face.counterflow_velocity_m_s / expected_relative - 1.0).abs() < 0.001);
        }
    }
}

#[test]
fn thermal_port_conserves_mass_and_entropy_and_counts_wave_work_separately() {
    let source = material(0.3);
    let mut duct = channel(64, source);
    let boundary = BoundaryVelocities {
        left: EndVelocities {
            normal_m_s: 3e-4,
            superfluid_m_s: -1e-4,
        },
        right: EndVelocities::default(),
    };
    let mut work = 0.0;
    let mut mass = 0.0;
    let mut entropy = 0.0;
    let mut heat = 0.0;
    for _ in 0..12 {
        let report = duct
            .step(
                duct.maximum_step_s() * 0.5,
                boundary,
                MutualFriction::None,
                || false,
            )
            .unwrap();
        work += report.boundary_wave_work_j;
        mass += report.boundary_mass_kg;
        entropy += report.boundary_entropy_j_k;
        heat += report.boundary_counterflow_heat_j;
    }
    let volume = duct.length_m() * duct.area_m2() / duct.cells() as f64;
    let mut stored_mass = 0.0;
    let mut stored_entropy = 0.0;
    for i in 0..duct.cells() {
        let state = duct.cell(i).unwrap();
        stored_mass += volume * state.density_perturbation_kg_m3;
        stored_entropy += volume
            * (source.density_kg_m3() * state.specific_entropy_perturbation_j_kg_k
                + source.entropy_j_kg_k() * state.density_perturbation_kg_m3);
    }
    assert!(stored_mass.abs() < 1e-17);
    assert!(mass.abs() < 1e-17);
    assert!((stored_entropy / entropy - 1.0).abs() < 2e-13);
    assert!((heat / (source.temperature_k() * entropy) - 1.0).abs() < 2e-13);
    assert!((duct.wave_energy_j().unwrap() / work - 1.0).abs() < 2e-12);
    assert!(
        heat > 1000.0 * work,
        "transported heat is not wave availability"
    );
}

#[test]
fn friction_laws_export_only_lost_kinetic_energy_to_thermal_consumers() {
    for law in [
        MutualFriction::Linear {
            gamma_kg_m3_s: 15.0,
        },
        MutualFriction::GorterMellink { a_m_s_kg: 2e6 },
        MutualFriction::IsotropicVortexTangle {
            circulation_quantum_m2_s: 1e-7,
            line_density_per_m2: 6e8,
            dissipative_b: 1.0,
        },
    ] {
        let mut duct = standing(64, material(-0.3), 0.001, 0.001);
        let initial = duct.wave_energy_j().unwrap();
        let mut exported = 0.0;
        for _ in 0..150 {
            let report = duct
                .step(
                    duct.maximum_step_s() * 0.8,
                    BoundaryVelocities::CLOSED,
                    law,
                    || false,
                )
                .unwrap();
            let cells: f64 = duct.last_mutual_friction_heat_j().iter().sum();
            assert!((cells - report.mutual_friction_heat_j).abs() < 1e-13 * initial);
            assert!(duct.last_mutual_friction_heat_j().iter().all(|x| *x >= 0.0));
            exported += cells;
        }
        assert!(exported > 0.1 * initial);
        assert!((duct.wave_energy_j().unwrap() + exported - initial).abs() < 1e-11 * initial);
        assert!((duct.cumulative_mutual_friction_heat_j() - exported).abs() < 1e-12 * initial);
    }
}

#[test]
fn isotropic_vortex_tangle_matches_scalar_projection_and_refuses_bad_inputs_atomically() {
    let mut tangle = standing(24, material(-0.3), 0.001, 0.001);
    let mut scalar = standing(24, material(-0.3), 0.001, 0.001);
    let law = MutualFriction::IsotropicVortexTangle {
        circulation_quantum_m2_s: 1e-7,
        line_density_per_m2: 6e8,
        dissipative_b: 1.0,
    };
    // Independently: gamma=(rho_n*rho_s/rho)*B*kappa*L/3=0.75*60/3.
    for _ in 0..20 {
        let dt = tangle.maximum_step_s() * 0.8;
        let a = tangle
            .step(dt, BoundaryVelocities::CLOSED, law, || false)
            .unwrap();
        let b = scalar
            .step(
                dt,
                BoundaryVelocities::CLOSED,
                MutualFriction::Linear {
                    gamma_kg_m3_s: 15.0,
                },
                || false,
            )
            .unwrap();
        assert_eq!(a, b);
    }
    assert_eq!(snapshot(&tangle), snapshot(&scalar));
    for (quantum, density, b) in [
        (0.0, 1.0, 1.0),
        (1.0, -1.0, 1.0),
        (1.0, 1.0, -1.0),
        (1.0, 1.0, f64::NAN),
        (f64::MAX, 2.0, 1.0),
    ] {
        let before = snapshot(&tangle);
        assert!(
            tangle
                .step(
                    tangle.maximum_step_s(),
                    BoundaryVelocities::CLOSED,
                    MutualFriction::IsotropicVortexTangle {
                        circulation_quantum_m2_s: quantum,
                        line_density_per_m2: density,
                        dissipative_b: b,
                    },
                    || false
                )
                .is_err()
        );
        assert_eq!(snapshot(&tangle), before);
    }
}

fn snapshot(duct: &TwoFluidChannel) -> (Vec<f64>, Vec<f64>) {
    let mut cells = vec![duct.time_s(), duct.cumulative_mutual_friction_heat_j()];
    for i in 0..duct.cells() {
        let c = duct.cell(i).unwrap();
        cells.extend([
            c.pressure_perturbation_pa,
            c.temperature_perturbation_k,
            c.density_perturbation_kg_m3,
            c.specific_entropy_perturbation_j_kg_k,
        ]);
    }
    cells.extend(duct.last_mutual_friction_heat_j());
    let mut faces = Vec::new();
    for i in 0..=duct.cells() {
        let f = duct.face(i).unwrap();
        faces.extend([
            f.normal_velocity_m_s,
            f.superfluid_velocity_m_s,
            f.counterflow_heat_flux_w_m2,
        ]);
    }
    (cells, faces)
}

#[test]
fn cancellation_is_atomic_at_early_factor_and_final_publication_boundaries() {
    let source = material(0.3);
    let law = MutualFriction::Linear { gamma_kg_m3_s: 2.0 };
    let mut baseline = standing(256, source, 0.001, 0.001);
    let dt = baseline.maximum_step_s() * 0.5;
    let mut total_polls = 0;
    baseline
        .step(dt, BoundaryVelocities::CLOSED, law, || {
            total_polls += 1;
            false
        })
        .unwrap();
    assert!(total_polls > 20);
    for stop in [1, total_polls / 2, total_polls] {
        let mut duct = standing(256, source, 0.001, 0.001);
        let before = snapshot(&duct);
        let mut calls = 0;
        assert_eq!(
            duct.step(dt, BoundaryVelocities::CLOSED, law, || {
                calls += 1;
                calls == stop
            }),
            Err(ChannelError::Cancelled)
        );
        assert_eq!(snapshot(&duct), before);
        // A cancelled cache/factor replacement cannot poison a retry.
        duct.step(dt, BoundaryVelocities::CLOSED, law, || false)
            .unwrap();
        assert_eq!(snapshot(&duct), snapshot(&baseline));
    }
}

#[test]
fn domain_refusal_leaves_the_last_physical_state_and_heat_untouched() {
    let mut duct = standing(16, material(0.0), 0.01, 0.001);
    let before = snapshot(&duct);
    for dt in [
        0.0,
        -1.0,
        f64::NAN,
        f64::INFINITY,
        2.0 * duct.maximum_step_s(),
    ] {
        assert!(
            duct.step(dt, BoundaryVelocities::CLOSED, MutualFriction::None, || {
                false
            })
            .is_err()
        );
        assert_eq!(snapshot(&duct), before);
    }
    let boundary = BoundaryVelocities {
        left: EndVelocities {
            normal_m_s: f64::NAN,
            superfluid_m_s: 0.0,
        },
        right: EndVelocities::default(),
    };
    assert!(
        duct.step(
            duct.maximum_step_s(),
            boundary,
            MutualFriction::None,
            || false
        )
        .is_err()
    );
    assert!(
        duct.initialize(&[0.0; 16], &[1.0; 16], &[0.0; 17], &[0.0; 17], || false)
            .is_err()
    );
    assert_eq!(snapshot(&duct), before);
}

#[test]
fn sourced_helium_reference_drives_both_physical_wave_modes() {
    use fs_material::superfluid::helium::HeliumSvpReference;
    let helium = HeliumSvpReference::at_1_8_k().unwrap();
    let source = *helium.state();
    let mut duct = TwoFluidChannel::try_new(
        source,
        0.1,
        1e-4,
        64,
        LinearizationLimits {
            max_pressure_perturbation_pa: 10.0,
            max_temperature_perturbation_k: 1e-4,
            max_component_speed_m_s: 0.01,
        },
    )
    .unwrap();
    let theta: Vec<_> = (0..64)
        .map(|i| 1e-6 * (std::f64::consts::PI * (i as f64 + 0.5) / 64.0).cos())
        .collect();
    duct.initialize(&[0.0; 64], &theta, &[0.0; 65], &[0.0; 65], || false)
        .unwrap();
    let before = duct.wave_energy_j().unwrap();
    let time = 0.1 / (4.0 * source.sound_modes().slow_m_s);
    advance_to(&mut duct, time, MutualFriction::None);
    assert!((duct.wave_energy_j().unwrap() / before - 1.0).abs() < 5e-12);
    assert!(duct.face(32).unwrap().counterflow_heat_flux_w_m2 > 1.0);
    assert!(
        duct.cell(0).unwrap().pressure_perturbation_pa.abs() > 1e-7,
        "the finite expansion must produce the physical pressure coupling"
    );
}

#[test]
fn heat_owner_refusal_and_post_candidate_cancellation_do_not_publish_wave_state() {
    use fs_duct::superfluid::HeatStepError;
    use std::cell::Cell;

    let mut duct = standing(32, material(0.3), 0.001, 0.001);
    let before = snapshot(&duct);
    let dt = 0.5 * duct.maximum_step_s();
    let friction = MutualFriction::Linear {
        gamma_kg_m3_s: 20.0,
    };
    let refused = duct.step_with_heat(
        dt,
        BoundaryVelocities::CLOSED,
        friction,
        || false,
        |report, heat| {
            assert!(report.mutual_friction_heat_j > 0.0);
            assert_eq!(heat.len(), 32);
            Err::<(), _>("thermal material domain")
        },
    );
    assert_eq!(refused, Err(HeatStepError::Heat("thermal material domain")));
    assert_eq!(snapshot(&duct), before);

    let prepared = Cell::new(false);
    let cancelled = duct.step_with_heat(
        dt,
        BoundaryVelocities::CLOSED,
        friction,
        || prepared.get(),
        |_, heat| {
            prepared.set(true);
            Ok::<_, &'static str>(heat.iter().sum::<f64>())
        },
    );
    assert!(
        prepared.get(),
        "cancellation must happen after preparing the external candidate"
    );
    assert_eq!(
        cancelled,
        Err(HeatStepError::Channel(ChannelError::Cancelled))
    );
    assert_eq!(snapshot(&duct), before);

    let (report, heat_candidate) = duct
        .step_with_heat(
            dt,
            BoundaryVelocities::CLOSED,
            friction,
            || false,
            |_, heat| Ok::<_, &'static str>(heat.iter().sum::<f64>()),
        )
        .unwrap();
    assert_eq!(duct.time_s(), dt);
    assert!(
        (heat_candidate - report.mutual_friction_heat_j).abs()
            < 1e-13 * report.wave_energy_before_j
    );
    assert_eq!(
        duct.cumulative_mutual_friction_heat_j(),
        report.mutual_friction_heat_j
    );
}
