use super::*;
use fs_math::det;
use std::cell::Cell;

const VOLUME: f64 = LENGTH * HEIGHT * WIDTH;

fn sealed<'c>(
    mesh: &'c TetComplex,
    xyz: &[[f64; 3]],
    medium: ThermoviscousMedium,
    thermodynamics: ThermoviscousThermodynamics,
) -> ThermoviscousFe<'c> {
    let beta_s = thermodynamics
        .adiabatic_compressibility_pa_inv(medium)
        .unwrap();
    let geometry = HelmholtzFe::new(
        mesh,
        xyz,
        2,
        Medium {
            density_kg_m3: medium.density_kg_m3,
            sound_speed_m_s: (medium.density_kg_m3 * beta_s).sqrt().recip(),
        },
    )
    .unwrap();
    ThermoviscousFe::new_with_thermodynamics(
        mesh,
        xyz,
        medium,
        thermodynamics,
        &ThermoviscousBoundary {
            no_slip_faces: geometry.boundary_faces(),
            isothermal_faces: vec![],
            ports: vec![],
        },
        ThermoviscousLimits::default(),
        &mut || false,
    )
    .unwrap()
}

fn closed_gas<'c>(
    mesh: &'c TetComplex,
    xyz: &[[f64; 3]],
    relaxing: bool,
) -> ThermoviscousInterior<'c> {
    let gas = gas();
    let medium = source_medium(&gas, 0.0).unwrap();
    ThermoviscousInterior::from_relaxation(
        sealed(
            mesh,
            xyz,
            medium,
            ThermoviscousThermodynamics::ideal_gas(medium).unwrap(),
        ),
        source(gas, relaxing, 1.0),
        BAND,
    )
    .unwrap()
}

fn scalar_close(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() < 2e-8 * scale,
        "{actual:e} != {expected:e}, scale={scale:e}",
    );
}

fn same_fields(actual: &ThermoviscousField, expected: &ThermoviscousField) {
    assert_eq!(actual.pressure_drive_pa, expected.pressure_drive_pa);
    assert_eq!(actual.inward_flow_m3_s, expected.inward_flow_m3_s);
    assert_eq!(actual.pressure_pa, expected.pressure_pa);
    assert_eq!(actual.velocity_m_s, expected.velocity_m_s);
    assert_eq!(actual.temperature_k, expected.temperature_k);
    for (a, b) in [
        (actual.power.input_w, expected.power.input_w),
        (
            actual.power.heat_source_work_w,
            expected.power.heat_source_work_w,
        ),
        (actual.power.viscous_w, expected.power.viscous_w),
        (actual.power.bulk_w, expected.power.bulk_w),
        (actual.power.conductive_w, expected.power.conductive_w),
        (
            actual.power.solid_conductive_w,
            expected.power.solid_conductive_w,
        ),
        (actual.power.caloric_w, expected.power.caloric_w),
        (
            actual.power.structural_damping_w,
            expected.power.structural_damping_w,
        ),
    ] {
        assert_eq!(a.to_bits(), b.to_bits());
    }
}

#[test]
fn original_molecular_heat_release_has_absolute_closed_volume_phase_and_work() {
    let (mesh, unit) = kuhn_cube(1);
    let xyz = points(&unit);
    let shape_integral = -2.5;
    let amplitude = C64::new(0.001, -0.0003);
    for relaxing in [false, true] {
        let model = closed_gas(&mesh, &xyz, relaxing);
        let calls = Cell::new(0);
        let heat = model
            .volume_heat_source(
                4,
                &mut |_, _| {
                    calls.set(calls.get() + 1);
                    Ok(shape_integral / VOLUME)
                },
                &mut || false,
            )
            .unwrap();
        let prepared_calls = calls.get();
        assert!(prepared_calls > 0);
        scalar_close(
            heat.integrated_shape(),
            shape_integral,
            shape_integral.abs(),
        );
        for frequency in [80.0, 1300.0] {
            let omega = std::f64::consts::TAU * frequency;
            let cp = model
                .relaxation()
                .unwrap()
                .caloric_response(omega)
                .unwrap()
                .cp_j_kg_k
                .conj();
            let gas = model.gas();
            let cv = cp - C64::from_re(gas.specific_gas_constant);
            let total_heat = amplitude.scale(shape_integral);
            // Independent sealed-volume EOS: p'=rho R T', Q'=i w rho V Cv T'.
            let temperature = total_heat / (cv * C64::new(0.0, omega * gas.density * VOLUME));
            let pressure = temperature.scale(gas.pressure / gas.temperature);
            let prepared = model.prepare_response(frequency, &mut || false).unwrap();
            let actual = prepared
                .response_with_heat_sources(&[], &[], &[(&heat, amplitude)], &mut || false)
                .unwrap();
            let original = model
                .inner()
                .prepare_with_heat_capacity(frequency, cp, &mut || false)
                .unwrap()
                .response_with_heat_sources(&[], &[], &[(&heat, amplitude)], &mut || false)
                .unwrap();
            same_fields(&actual, &original);
            same_fields(
                &actual,
                &model
                    .response_with_heat_sources(
                        frequency,
                        &[],
                        &[],
                        &[(&heat, amplitude)],
                        &mut || false,
                    )
                    .unwrap(),
            );
            for &p in &actual.pressure_pa {
                close(p, pressure, pressure.abs(), 2e-8);
            }
            for &t in &actual.temperature_k[..mesh.vertex_count] {
                close(t, temperature, temperature.abs(), 2e-8);
            }
            for &t in &actual.temperature_k[mesh.vertex_count..] {
                assert!(t.abs() < 2e-8 * temperature.abs());
            }
            let work = 0.5 * (total_heat * temperature.conj()).re / gas.temperature;
            let scale = total_heat.abs() * temperature.abs() / gas.temperature;
            scalar_close(actual.power.heat_source_work_w, work, scale);
            scalar_close(actual.power.caloric_w, work, scale);
            assert!(actual.power.balance_residual_w().abs() < 2e-8 * scale);
            if relaxing {
                assert!(work > 0.0);
            }
            // Coherent signed sources cancel before the original field solve.
            let cancelled = prepared
                .response_with_heat_sources(
                    &[],
                    &[],
                    &[(&heat, amplitude), (&heat, -amplitude)],
                    &mut || false,
                )
                .unwrap();
            assert!(cancelled.pressure_pa.iter().all(|p| *p == C64::ZERO));
            assert!(cancelled.temperature_k.iter().all(|t| *t == C64::ZERO));
        }
        assert_eq!(calls.get(), prepared_calls);
    }
}

#[test]
fn material_heat_release_retains_signed_expansion_and_original_card() {
    let (mesh, unit) = kuhn_cube(1);
    let xyz = points(&unit);
    let frequency = 900.0;
    let amplitude = C64::new(0.001, 0.0003);
    for alpha in [-2e-4, 0.0, 2e-4] {
        let (card, source) = fluid::liquid(Some(alpha), 0.004, 1.0);
        let (medium, thermo) = source_fluid_medium(&source, 101325.0, 0.002).unwrap();
        let model =
            ThermoviscousInterior::from_material(sealed(&mesh, &xyz, medium, thermo), source, BAND)
                .unwrap();
        assert_eq!(
            model.material_source().material().card_identity(),
            card.content_hash()
        );
        let heat = model
            .volume_heat_source(4, &mut |_, _| Ok(1.0 / VOLUME), &mut || false)
            .unwrap();
        let actual = model
            .response_with_heat_sources(frequency, &[], &[], &[(&heat, amplitude)], &mut || false)
            .unwrap();
        let beta_t = thermo.isothermal_compressibility_pa_inv;
        let cv_volume =
            medium.density_kg_m3 * medium.cp_j_kg_k - medium.temperature_k * alpha * alpha / beta_t;
        let temperature =
            amplitude / C64::new(0.0, std::f64::consts::TAU * frequency * VOLUME * cv_volume);
        let pressure = temperature.scale(alpha / beta_t);
        for &t in &actual.temperature_k[..mesh.vertex_count] {
            close(t, temperature, temperature.abs(), 2e-8);
        }
        for &p in &actual.pressure_pa {
            close(p, pressure, pressure.abs().max(1e-15), 2e-8);
        }
        // The same material closure owns causal heating, including alpha=0.
        let dt = 1e-7;
        let causal = model
            .prepare_time_with_heat_sources(
                dt,
                &[&heat],
                ThermoviscousTimeLimits::default(),
                &mut || false,
            )
            .unwrap();
        let step = causal
            .step(&causal.zero_state(), &[], &[], &[0.001], &mut || false)
            .unwrap();
        let temperature = 0.001 * dt / (VOLUME * cv_volume);
        for &t in &step.field.state.temperature_k[..mesh.vertex_count] {
            scalar_close(t, temperature, temperature.abs());
        }
        let pressure = alpha / beta_t * temperature;
        for &p in &step.field.state.pressure_pa {
            scalar_close(p, pressure, pressure.abs().max(1e-15));
        }
    }
}

#[test]
fn loaded_heat_and_incident_waves_share_original_calorics_and_affine_emission() {
    let (mesh, unit) = kuhn_cube(1);
    let xyz = points(&unit);
    let gas = gas();
    let model = ThermoviscousInterior::from_relaxation(
        build(&mesh, &xyz, &gas, true),
        source(gas, true, 1.0),
        BAND,
    )
    .unwrap();
    let heat = model
        .volume_heat_source(4, &mut |_, p| Ok(2.0 * p[0] / LENGTH / VOLUME), &mut || {
            false
        })
        .unwrap();
    let frequency = 700.0;
    let reference: [f64; 2] = [1.3e8, 0.7e8];
    let incident = [C64::new(0.0001, 0.0002), C64::new(-0.0002, 0.00005)];
    let sources = [(&heat, C64::new(0.0003, -0.0002))];
    let loaded = model
        .prepare_scattering(frequency, &reference, &mut || false)
        .unwrap();
    let actual = loaded
        .response_with_sources(&incident, &[], &sources, &mut || false)
        .unwrap();
    let cp = model
        .relaxation()
        .unwrap()
        .caloric_response(std::f64::consts::TAU * frequency)
        .unwrap()
        .cp_j_kg_k
        .conj();
    let original = model
        .inner()
        .prepare_scattering(frequency, cp, &reference, &mut || false)
        .unwrap()
        .response_with_sources(&incident, &[], &sources, &mut || false)
        .unwrap();
    same_fields(&actual, &original);
    let offset = loaded
        .source_outgoing(&[], &sources, &mut || false)
        .unwrap();
    let scattering = loaded.scattering(&mut || false).unwrap();
    let duct = model
        .duct_scattering(frequency, &reference, &mut || false)
        .unwrap();
    assert_eq!(
        duct,
        scattering.iter().map(|s| s.conj()).collect::<Vec<_>>()
    );
    let mut boundary_work = 0.0;
    for i in 0..2 {
        let outgoing = (actual.pressure_drive_pa[i]
            - actual.inward_flow_m3_s[i].scale(reference[i]))
        .scale(0.5 / reference[i].sqrt());
        let expected = (0..2).fold(offset[i], |sum, j| {
            sum + scattering[2 * i + j] * incident[j]
        });
        close(outgoing, expected, outgoing.abs() + expected.abs(), 2e-8);
        boundary_work += 0.5 * (incident[i].norm_sq() - outgoing.norm_sq());
    }
    let scale = actual.power.dissipated_w().abs()
        + actual.power.heat_source_work_w.abs()
        + boundary_work.abs();
    scalar_close(
        actual.power.input_w,
        boundary_work + actual.power.heat_source_work_w,
        scale,
    );
    assert!(actual.power.balance_residual_w().abs() < 2e-8 * scale);
}

#[test]
fn held_heat_preserves_original_molecular_lag_and_integrated_signed_work() {
    let (mesh, unit) = kuhn_cube(1);
    let xyz = points(&unit);
    let dt = 1e-4;
    let drive = -0.001;
    for relaxing in [false, true] {
        let model = closed_gas(&mesh, &xyz, relaxing);
        let heat = model
            .volume_heat_source(4, &mut |_, _| Ok(1.0 / VOLUME), &mut || false)
            .unwrap();
        let prepared = model
            .prepare_time_with_heat_sources(
                dt,
                &[&heat],
                ThermoviscousTimeLimits::default(),
                &mut || false,
            )
            .unwrap();
        let source = model.relaxation().unwrap();
        assert_eq!(
            prepared.frozen_cp_j_kg_k(),
            source.frozen_cv_j_kg_k() + model.gas().specific_gas_constant
        );
        assert_eq!(
            prepared.caloric_reservoirs().len(),
            source.processes().len()
        );
        for (actual, original) in prepared.caloric_reservoirs().iter().zip(source.processes()) {
            assert_eq!(actual.heat_capacity_j_kg_k, original.heat_capacity_j_kg_k);
            assert_eq!(actual.relaxation_time_s, original.relaxation_time_s);
        }
        let initial = prepared.zero_state();
        let original_initial = initial.clone();
        let step = prepared
            .step(&initial, &[], &[], &[drive], &mut || false)
            .unwrap();
        assert_eq!(initial, original_initial);
        let original = model
            .inner()
            .prepare_time_with_heat_sources(
                dt,
                prepared.frozen_cp_j_kg_k(),
                prepared.caloric_reservoirs(),
                None,
                None,
                &[&heat],
                ThermoviscousTimeLimits::default(),
                &mut || false,
            )
            .unwrap()
            .step(&initial, &[], &[], &[drive], &mut || false)
            .unwrap();
        assert_eq!(step.field.state, original.field.state);
        assert_eq!(
            step.heat_source_work_by_source_j,
            original.heat_source_work_by_source_j
        );
        let mass = model.gas().density * VOLUME;
        let cv = source.frozen_cv_j_kg_k();
        let (temperature, temperature_integral) = if relaxing {
            assert_eq!(source.processes().len(), 1);
            let process = &source.processes()[0];
            let capacity = process.heat_capacity_j_kg_k;
            let equilibrium_cv = cv + capacity;
            let rate = (1.0 + capacity / cv) / process.relaxation_time_s;
            let decay = det::exp(-rate * dt);
            let slope = drive / (mass * cv);
            let difference = slope / rate * (1.0 - decay);
            let equilibrium = drive * dt / (mass * equilibrium_cv);
            let temperature = equilibrium + capacity / equilibrium_cv * difference;
            let molecular = equilibrium - cv / equilibrium_cv * difference;
            for &t in &step.field.state.molecular_temperature_k[0][..mesh.vertex_count] {
                scalar_close(t, molecular, molecular.abs());
            }
            let integral = drive * dt * dt / (2.0 * mass * equilibrium_cv)
                + capacity / equilibrium_cv * slope / rate * (dt - (1.0 - decay) / rate);
            let dissipation = mass * capacity
                / (process.relaxation_time_s * model.gas().temperature)
                * (slope / rate).powi(2)
                * (dt - 2.0 * (1.0 - decay) / rate + (1.0 - decay * decay) / (2.0 * rate));
            scalar_close(
                step.field.work.molecular_relaxation_heat_j,
                dissipation,
                dissipation.abs(),
            );
            (temperature, integral)
        } else {
            let temperature = drive * dt / (mass * cv);
            (temperature, 0.5 * dt * temperature)
        };
        for &t in &step.field.state.temperature_k[..mesh.vertex_count] {
            scalar_close(t, temperature, temperature.abs());
        }
        let pressure = model.gas().pressure / model.gas().temperature * temperature;
        for &p in &step.field.state.pressure_pa {
            scalar_close(p, pressure, pressure.abs());
        }
        let work = drive * temperature_integral / model.gas().temperature;
        scalar_close(step.field.work.heat_source_work_j, work, work.abs());
        scalar_close(step.heat_source_work_by_source_j[0], work, work.abs());
        assert!(step.field.work.balance_residual_j().abs() < 2e-8 * work.abs());
        if !relaxing {
            let returned = prepared
                .step(&step.field.state, &[], &[], &[-drive], &mut || false)
                .unwrap();
            scalar_close(returned.field.work.heat_source_work_j, -work, work.abs());
            assert!(returned.field.work.final_energy_j.abs() < 2e-8 * work.abs());
        }
    }
}

#[test]
fn causal_heat_keeps_original_material_wall_depth_and_tangential_transport() {
    use super::super::wall::{MaterialTangentialWallRegion, MaterialWallRegion};
    use fs_feec::acoustics::helmholtz::boundary_layer::BoundaryLayerGeometry;
    use fs_feec::acoustics::helmholtz::boundary_layer::thermal_wall::{
        ThermalBacking, ThermalWallDepth,
    };

    let (mesh, unit) = kuhn_cube(1);
    let xyz = points(&unit);
    let gas = gas();
    let geometry = HelmholtzFe::new(
        &mesh,
        &xyz,
        2,
        Medium {
            density_kg_m3: gas.density,
            sound_speed_m_s: gas.sound_speed,
        },
    )
    .unwrap();
    let face = geometry
        .boundary_faces()
        .into_iter()
        .find(|f| f.iter().all(|&v| xyz[v as usize][1] == 0.0))
        .unwrap();
    let material = time::wall_transport(1000.0, 2.0, 0.003);
    let sample = material.at_temperature(gas.temperature).unwrap();
    for tangential in [false, true] {
        let model = ThermoviscousInterior::from_gas(build(&mesh, &xyz, &gas, false), gas, BAND)
            .unwrap()
            .with_conductive_walls(
                &[MaterialWallRegion {
                    faces: &[face],
                    sample: &sample,
                    depth: ThermalWallDepth::Slab {
                        thickness_m: 0.00002,
                        backing: ThermalBacking::Adiabatic,
                    },
                    minimum_tangential_scale_m: 0.004,
                    geometry: BoundaryLayerGeometry {
                        minimum_curvature_radius_m: f64::INFINITY,
                        minimum_half_gap_m: 0.0001,
                        maximum_layer_ratio: 0.05,
                    },
                }],
                &mut || false,
            )
            .unwrap();
        let model = if tangential {
            model
                .with_tangentially_isotropic_wall_transport(
                    &[MaterialTangentialWallRegion {
                        normal: &material,
                        tangential: &material,
                    }],
                    &mut || false,
                )
                .unwrap()
        } else {
            model
        };
        let heat = model
            .volume_heat_source(4, &mut |_, p| Ok(2.0 * p[0] / LENGTH / VOLUME), &mut || {
                false
            })
            .unwrap();
        assert!(
            model
                .prepare_time_with_heat_sources(
                    1e-7,
                    &[&heat],
                    ThermoviscousTimeLimits::default(),
                    &mut || false,
                )
                .is_err()
        );
        let prepared = model
            .prepare_time_with_wall_diffusion_and_heat_sources(
                1e-7,
                &[2],
                &[&heat],
                ThermoviscousTimeLimits::default(),
                &mut || false,
            )
            .unwrap();
        let transport = [TangentialWallConductivity {
            isotropic_w_m_k: sample.conductivity_w_m_k(),
        }];
        let original = model
            .inner()
            .prepare_time_with_heat_sources(
                1e-7,
                gas.specific_heat_cp,
                &[],
                Some(ThermoviscousWallDiscretization {
                    normal_elements: &[2],
                    maximum_frequency_hz: BAND[1],
                }),
                if tangential { Some(&transport) } else { None },
                &[&heat],
                ThermoviscousTimeLimits::default(),
                &mut || false,
            )
            .unwrap();
        let initial = prepared.zero_state();
        assert_eq!(initial, original.zero_state());
        let actual = prepared
            .step(&initial, &[0.0; 2], &[], &[0.001], &mut || false)
            .unwrap();
        let expected = original
            .step(&initial, &[0.0; 2], &[], &[0.001], &mut || false)
            .unwrap();
        assert_eq!(actual.field.state, expected.field.state);
        assert_eq!(
            actual.heat_source_work_by_source_j,
            expected.heat_source_work_by_source_j
        );
        assert_eq!(
            actual.field.work.solid_conductive_heat_j,
            expected.field.work.solid_conductive_heat_j
        );
        assert!(actual.field.work.solid_conductive_heat_j > 0.0);
        assert!(
            actual.field.state.conductive_wall_temperature_k[0]
                .iter()
                .any(|t| *t != 0.0)
        );
    }
}

#[test]
fn heat_source_refusals_preserve_owner_band_clock_and_candidate_state() {
    let (mesh, unit) = kuhn_cube(1);
    let xyz = points(&unit);
    let model = closed_gas(&mesh, &xyz, true);
    let foreign = closed_gas(&mesh, &xyz, true);
    let heat = model
        .volume_heat_source(4, &mut |_, _| Ok(1.0 / VOLUME), &mut || false)
        .unwrap();
    let other = foreign
        .volume_heat_source(4, &mut |_, _| Ok(1.0 / VOLUME), &mut || false)
        .unwrap();
    assert!(
        model
            .response_with_heat_sources(300.0, &[], &[], &[(&other, C64::ONE)], &mut || false)
            .is_err()
    );
    assert!(
        model
            .prepare_response(BAND[0] * 0.9, &mut || false)
            .is_err()
    );
    assert!(
        model
            .prepare_scattering(BAND[1] * 1.1, &[], &mut || false)
            .is_err()
    );
    assert!(model.prepare_response(300.0, &mut || true).is_err());
    for (dt, sources) in [
        (1e-4, vec![&other]),
        (1.0 / BAND[1], vec![&heat]),
        (1e-4, vec![]),
    ] {
        assert!(
            model
                .prepare_time_with_heat_sources(
                    dt,
                    &sources,
                    ThermoviscousTimeLimits::default(),
                    &mut || false
                )
                .is_err()
        );
    }
    let mut limits = ThermoviscousTimeLimits::default();
    limits.max_states = 1;
    assert!(
        model
            .prepare_time_with_heat_sources(1e-4, &[&heat], limits, &mut || false)
            .is_err()
    );
    let prepared = model
        .prepare_time_with_heat_sources(
            1e-4,
            &[&heat],
            ThermoviscousTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let initial = prepared.zero_state();
    let before = initial.clone();
    assert!(
        prepared
            .step(&initial, &[], &[], &[f64::NAN], &mut || false)
            .is_err()
    );
    assert!(
        prepared
            .step(&initial, &[], &[], &[1.0], &mut || true)
            .is_err()
    );
    assert_eq!(initial, before);
}
