//! G1/G3/G4: fixed-geometry material substitution preserves canonical motion,
//! reports external work, and updates forces, losses and the pressure observer.
use super::*;

fn specimen(
    name: &str,
    density: f64,
    young: f64,
    viscosity: f64,
    modes: usize,
) -> fs_couple::string_specimen::ResolvedStringSpecimen {
    let material = state_with_domain(
        name,
        &[
            ("density", QuantitySpec::dimensional(Density::DIMS), density),
            (
                "young_modulus",
                QuantitySpec::dimensional(Pressure::DIMS),
                young,
            ),
            (
                KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY,
                QuantitySpec::dimensional(DynViscosity::DIMS),
                viscosity,
            ),
        ],
        ValidityDomain::unconstrained().with("omega", 1.0, 20_000.0),
        QueryPoint::new().with("omega", 1.0).unwrap(),
        None,
    );
    with_uniform_circular_material_state(
        PrestressedString {
            n_modes: modes,
            ..loss_template()
        },
        0.0005,
        &material,
    )
    .unwrap()
    .with_kelvin_voigt_bending_loss()
    .unwrap()
}

fn original() -> fs_couple::string_specimen::ResolvedStringSpecimen {
    specimen("original synthetic wire", 7800.0, 200e9, 1e6, 2)
}

fn replacement() -> fs_couple::string_specimen::ResolvedStringSpecimen {
    specimen("replacement synthetic wire", 15600.0, 80e9, 4e6, 2)
}

fn running(cx: &fs_exec::Cx<'_>) -> (LinearMaterialStringRuntime, Vec<f64>) {
    let mut runtime = LinearMaterialStringRuntime::try_new(
        cx,
        original(),
        Some(Pluck {
            station_frac: 0.4,
            height_m: 1e-5,
        }),
        incremental_ambient(),
        1.0,
        48_000,
    )
    .unwrap();
    let forces = runtime.project_point_forces(cx, &loads()).unwrap();
    for _ in 0..8 {
        runtime.step(cx, &forces).unwrap();
    }
    (runtime, forces)
}

fn loads() -> [StringPointForce; 1] {
    [StringPointForce {
        station_frac: 0.3,
        force_n: 0.25,
    }]
}

#[test]
fn g1_substitution_preserves_displacement_momentum_and_rebuilds_actual_sound() {
    with_string_cx(|cx, _| {
        let (mut runtime, old_force) = running(cx);
        let old = runtime.specimen().string();
        let before = runtime.states().to_vec();
        let old_modes = runtime.modes().to_vec();
        let old_point = runtime.point_motion(cx, 0.3).unwrap();
        let old_mass = old.lin_density_kg_m * old.length_m / 2.0;
        let source = replacement();
        let new = source.string();
        let new_mass = new.lin_density_kg_m * new.length_m / 2.0;
        let scale = (new_mass / old_mass).sqrt();
        // The ordinary thermal/material update keeps its stricter source gate.
        assert!(runtime.rebind(cx, source.clone(), 1.0).is_err());
        let update = runtime
            .substitute_material(cx, source.clone(), 1.0)
            .unwrap();
        assert_eq!(runtime.specimen(), &source);
        assert_eq!(runtime.accepted_samples(), 8);
        assert_eq!(runtime.epoch(), 9);
        close(
            runtime.point_motion(cx, 0.3).unwrap().displacement_m,
            old_point.displacement_m,
        );
        close(
            runtime.point_motion(cx, 0.3).unwrap().velocity_m_s * new_mass,
            old_point.velocity_m_s * old_mass,
        );
        let mut expected_energy = 0.0;
        for (k, (old_state, new_state)) in before.iter().zip(runtime.states()).enumerate() {
            let displacement = old_state.displacement_m_sqrt_kg / old_mass.sqrt();
            let momentum = old_state.velocity_m_sqrt_kg_per_s * old_mass.sqrt();
            close(
                new_state.displacement_m_sqrt_kg / new_mass.sqrt(),
                displacement,
            );
            close(
                new_state.velocity_m_sqrt_kg_per_s * new_mass.sqrt(),
                momentum,
            );
            let wave = (k + 1) as f64 * core::f64::consts::PI / new.length_m;
            let stiffness = (new.tension_n * wave.powi(2)
                + new.bending_stiffness_n_m2 * wave.powi(4))
                * new.length_m
                / 2.0;
            expected_energy +=
                momentum.powi(2) / (2.0 * new_mass) + 0.5 * stiffness * displacement.powi(2);
        }
        close(update.energy_after_j, expected_energy);
        close(
            update.parameter_work_j,
            expected_energy - update.energy_before_j,
        );
        assert!(update.parameter_work_j < 0.0);
        assert_ne!(runtime.modes(), old_modes);

        // Same physical force, new mass normalization. Independently solve the
        // forced damped oscillator and compact strip pressure, including the
        // even mode's accepted physical acceleration history across the swap.
        let force = runtime.project_point_forces(cx, &loads()).unwrap();
        for (new_force, old_force) in force.iter().zip(&old_force) {
            close(new_force * scale, *old_force);
        }
        let dt: f64 = 1.0 / 48_000.0;
        let gas = *runtime.ambient();
        let factor = gas.density * new.width_m / (4.0 * core::f64::consts::PI * new_mass.sqrt());
        let mut expected_pressure = 0.0;
        let mut expected_states = Vec::new();
        for (k, ((mode, state), &input)) in runtime
            .modes()
            .iter()
            .zip(runtime.states())
            .zip(&force)
            .enumerate()
        {
            let omega = mode.angular_frequency_rad_s;
            let decay = mode.damping_ratio * omega;
            let wd = (omega.powi(2) - decay.powi(2)).sqrt();
            let qeq = input / omega.powi(2);
            let q = state.displacement_m_sqrt_kg - qeq;
            let v = state.velocity_m_sqrt_kg_per_s;
            let next_q = qeq
                + (-decay * dt).exp()
                    * (q * (wd * dt).cos() + (v + decay * q) / wd * (wd * dt).sin());
            let next_v = (-decay * dt).exp()
                * (v * (wd * dt).cos() - (decay * v + omega.powi(2) * q) / wd * (wd * dt).sin());
            expected_states.push((next_q, next_v));
            let acceleration = input - omega.powi(2) * next_q - 2.0 * decay * next_v;
            let previous = (old_force[k]
                - old_modes[k].angular_frequency_rad_s.powi(2) * before[k].displacement_m_sqrt_kg
                - 2.0
                    * old_modes[k].damping_ratio
                    * old_modes[k].angular_frequency_rad_s
                    * before[k].velocity_m_sqrt_kg_per_s)
                * scale;
            let kpi = (k + 1) as f64 * core::f64::consts::PI;
            expected_pressure += if k % 2 == 0 {
                factor * 2.0 * new.length_m / kpi * acceleration
            } else {
                -factor * new.length_m.powi(2) / (kpi * gas.sound_speed) * (acceleration - previous)
                    / dt
            };
        }
        let frame = runtime.step(cx, &force).unwrap();
        for (state, (q, v)) in runtime.states().iter().zip(expected_states) {
            close(state.displacement_m_sqrt_kg, q);
            close(state.velocity_m_sqrt_kg_per_s, v);
        }
        close(frame.acoustic.observer_pressure_pa, expected_pressure);
        assert!(frame.dissipation.material_heat_j > 0.0);
        assert!(frame.dissipation.air_loss_j > 0.0);
        assert_eq!(frame.dissipation.authored_loss_j, 0.0);
    });
}

#[test]
fn g3_substitution_round_trip_and_identical_source_preserve_motion() {
    with_string_cx(|cx, _| {
        let (mut runtime, force) = running(cx);
        let mut reference = runtime.clone();
        let before = runtime.states().to_vec();
        let forward = runtime.substitute_material(cx, replacement(), 1.0).unwrap();
        let reverse = runtime.substitute_material(cx, original(), 1.0).unwrap();
        close(forward.parameter_work_j, -reverse.parameter_work_j);
        assert_eq!(runtime.specimen(), reference.specimen());
        assert_eq!(runtime.modes(), reference.modes());
        for (returned, initial) in runtime.states().iter().zip(&before) {
            close(
                returned.displacement_m_sqrt_kg,
                initial.displacement_m_sqrt_kg,
            );
            close(
                returned.velocity_m_sqrt_kg_per_s,
                initial.velocity_m_sqrt_kg_per_s,
            );
        }
        close(
            runtime
                .step(cx, &force)
                .unwrap()
                .acoustic
                .observer_pressure_pa,
            reference
                .step(cx, &force)
                .unwrap()
                .acoustic
                .observer_pressure_pa,
        );
        let mut same = reference.clone();
        let update = same.substitute_material(cx, original(), 0.0).unwrap();
        assert_eq!(update.parameter_work_j, 0.0);
        assert_eq!(same.states(), reference.states());
        for _ in 0..8 {
            // Epoch is the only difference; no sample/history is reset.
            let a = same.step(cx, &force).unwrap();
            let b = reference.step(cx, &force).unwrap();
            assert_eq!(a.acoustic, b.acoustic);
            assert_eq!(a.dissipation, b.dissipation);
        }
    });
}

#[test]
fn g4_substitution_refusal_preserves_complete_pressure_continuation() {
    with_string_cx(|cx, _| {
        let (mut runtime, force) = running(cx);
        let mut reference = runtime.clone();
        for budget in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(
                runtime
                    .substitute_material(cx, replacement(), budget)
                    .is_err()
            );
        }
        for (length, radius, modes) in [(0.501, 0.0005, 2), (0.5, 0.0006, 2), (0.5, 0.0005, 1)] {
            let source = replacement();
            let changed = with_uniform_circular_material_state(
                PrestressedString {
                    length_m: length,
                    n_modes: modes,
                    ..source.string()
                },
                radius,
                source.material(),
            )
            .unwrap();
            assert!(runtime.substitute_material(cx, changed, 1.0).is_err());
        }
        with_string_cx(|cancelled, gate| {
            gate.request();
            assert_eq!(
                runtime.substitute_material(cancelled, replacement(), 1.0),
                Err(AcousticRealizeError::Cancelled)
            );
        });
        assert_eq!(runtime.specimen(), reference.specimen());
        assert_eq!(runtime.states(), reference.states());
        assert_eq!(runtime.modes(), reference.modes());
        assert_eq!(runtime.epoch(), reference.epoch());
        assert_eq!(runtime.accepted_samples(), reference.accepted_samples());
        for _ in 0..8 {
            assert_eq!(
                runtime.step(cx, &force).unwrap(),
                reference.step(cx, &force).unwrap()
            );
        }
    });
}
