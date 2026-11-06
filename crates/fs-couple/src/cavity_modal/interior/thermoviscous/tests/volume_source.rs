use super::super::wall::{MaterialTangentialWallRegion, MaterialWallRegion};
use super::*;
use fs_feec::acoustics::helmholtz::boundary_layer::BoundaryLayerGeometry;
use fs_feec::acoustics::helmholtz::boundary_layer::thermal_wall::{
    ThermalBacking, ThermalWallDepth,
};
use fs_math::det;
use std::cell::Cell;

const VOLUME: f64 = LENGTH * HEIGHT * WIDTH;

#[test]
fn original_material_forces_share_molecular_calorics_and_signed_harmonic_work() {
    let (mesh, unit) = kuhn_cube(1);
    let xyz = points(&unit);
    let gas = gas();
    let model = ThermoviscousInterior::from_relaxation(
        build(&mesh, &xyz, &gas, false),
        source(gas, true, 1.0),
        BAND,
    )
    .unwrap();
    let calls = Cell::new(0usize);
    let force = model
        .volume_force_source(
            4,
            &mut |_, x| {
                calls.set(calls.get() + 1);
                Ok([(1.0 + x[0] / LENGTH) / VOLUME, -0.5 / VOLUME, 0.0])
            },
            &mut || false,
        )
        .unwrap();
    let heat = model
        .volume_heat_source(4, &mut |_, _| Ok(1.0 / VOLUME), &mut || false)
        .unwrap();
    let sampled = calls.get();
    let frequency = 733.0;
    let traction = [C64::new(0.5, 0.1), C64::new(-0.2, 0.3)];
    let forces = [(&force, C64::new(2e-6, -3e-7))];
    let heating = [(&heat, C64::new(2e-4, 5e-5))];
    let cp = model
        .relaxation()
        .unwrap()
        .caloric_response(std::f64::consts::TAU * frequency)
        .unwrap()
        .cp_j_kg_k
        .conj();
    let original = model
        .inner()
        .prepare_with_heat_capacity(frequency, cp, &mut || false)
        .unwrap()
        .response_with_volume_sources(&traction, &[], &forces, &heating, &mut || false)
        .unwrap();
    let bound = model
        .response_with_volume_sources(frequency, &traction, &[], &forces, &heating, &mut || false)
        .unwrap();
    assert_eq!(bound.pressure_pa, original.pressure_pa);
    assert_eq!(bound.velocity_m_s, original.velocity_m_s);
    assert_eq!(bound.temperature_k, original.temperature_k);
    assert_eq!(bound.inward_flow_m3_s, original.inward_flow_m3_s);
    assert_eq!(
        bound.power.body_force_work_w,
        original.power.body_force_work_w
    );
    assert_eq!(
        bound.power.heat_source_work_w,
        original.power.heat_source_work_w
    );
    assert_eq!(bound.power.caloric_w, original.power.caloric_w);
    assert!(bound.power.body_force_work_w.abs() > 0.0);
    assert!(bound.power.caloric_w > 0.0);
    let scale = bound
        .power
        .input_w
        .abs()
        .max(bound.power.dissipated_w().abs());
    assert!(bound.power.balance_residual_w().abs() < 2e-8 * scale);
    assert_eq!(calls.get(), sampled);
    assert!(
        model
            .response_with_volume_sources(
                BAND[1] + 1.0,
                &traction,
                &[],
                &forces,
                &heating,
                &mut || false,
            )
            .is_err()
    );
    let foreign_model =
        ThermoviscousInterior::from_gas(build(&mesh, &xyz, &gas, false), gas, BAND).unwrap();
    let foreign = foreign_model
        .volume_force_source(4, &mut |_, _| Ok([1.0; 3]), &mut || false)
        .unwrap();
    assert!(
        model
            .response_with_volume_sources(
                frequency,
                &traction,
                &[],
                &[(&foreign, C64::ONE)],
                &[],
                &mut || false,
            )
            .is_err()
    );
    assert!(matches!(
        model.volume_force_source(4, &mut |_, _| Ok([0.0; 3]), &mut || true,),
        Err(ThermoviscousSourceError::Interior(
            ThermoviscousError::Cancelled
        ))
    ));
}

#[test]
fn duct_affine_source_scattering_converts_every_phasor_and_closes_matched_power() {
    let (mesh, unit) = kuhn_cube(1);
    let xyz = points(&unit);
    let gas = gas();
    let model = ThermoviscousInterior::from_relaxation(
        build(&mesh, &xyz, &gas, false),
        source(gas, true, 1.0),
        BAND,
    )
    .unwrap();
    let force = model
        .volume_force_source(4, &mut |_, _| Ok([1.0 / VOLUME, 0.0, 0.0]), &mut || false)
        .unwrap();
    let heat = model
        .volume_heat_source(4, &mut |_, x| Ok(2.0 * x[0] / LENGTH / VOLUME), &mut || {
            false
        })
        .unwrap();
    let frequency = 600.0;
    let z = gas.density * gas.sound_speed / (HEIGHT * WIDTH);
    let references = [z, 1.7 * z];
    let force_minus = C64::new(2e-6, 8e-7);
    let heat_minus = C64::new(2e-4, -6e-5);
    let (scattering, outgoing) = model
        .duct_scattering_with_volume_sources(
            frequency,
            &references,
            &[],
            &[(&force, force_minus)],
            &[(&heat, heat_minus)],
            &mut || false,
        )
        .unwrap();
    assert_eq!(
        scattering,
        model
            .duct_scattering(frequency, &references, &mut || false)
            .unwrap()
    );
    assert!(outgoing.iter().any(|z| z.abs() > 0.0));
    let prepared = model
        .prepare_scattering(frequency, &references, &mut || false)
        .unwrap();
    let incident_minus = [C64::new(2e-5, -1e-5), C64::new(-1e-5, 3e-5)];
    let incident_plus = incident_minus.map(C64::conj);
    let field = prepared
        .response_with_volume_sources(
            &incident_plus,
            &[],
            &[(&force, force_minus.conj())],
            &[(&heat, heat_minus.conj())],
            &mut || false,
        )
        .unwrap();
    for i in 0..2 {
        let predicted = outgoing[i]
            + (0..2).fold(C64::ZERO, |sum, j| {
                sum + scattering[i * 2 + j] * incident_minus[j]
            });
        // b = a - sqrt(R) Q for inward-positive Q, converted back to exp(-i wt).
        let original = incident_minus[i]
            - field.inward_flow_m3_s[i]
                .conj()
                .scale(det::sqrt(references[i]));
        close(
            predicted,
            original,
            predicted.abs().max(original.abs()),
            3e-9,
        );
    }
    let matched = prepared
        .response_with_volume_sources(
            &[C64::ZERO; 2],
            &[],
            &[(&force, force_minus.conj())],
            &[(&heat, heat_minus.conj())],
            &mut || false,
        )
        .unwrap();
    let radiated = 0.5
        * outgoing
            .iter()
            .map(|z| z.re * z.re + z.im * z.im)
            .sum::<f64>();
    let supplied = matched.power.body_force_work_w + matched.power.heat_source_work_w;
    let expected = matched.power.dissipated_w() + radiated;
    assert!((supplied - expected).abs() < 2e-8 * supplied.abs().max(expected.abs()));
    let (passive, silent) = model
        .duct_scattering_with_volume_sources(frequency, &references, &[], &[], &[], &mut || false)
        .unwrap();
    assert_eq!(passive, scattering);
    assert!(silent.iter().all(|z| *z == C64::ZERO));
}

#[test]
fn source_bound_causal_forces_keep_every_original_molecular_reservoir() {
    let (mesh, unit) = kuhn_cube(1);
    let xyz = points(&unit);
    let gas = gas();
    let model = ThermoviscousInterior::from_relaxation(
        build(&mesh, &xyz, &gas, false),
        source(gas, true, 1.0),
        BAND,
    )
    .unwrap();
    let force = model
        .volume_force_source(4, &mut |_, _| Ok([1.0 / VOLUME, 0.0, 0.0]), &mut || false)
        .unwrap();
    let heat = model
        .volume_heat_source(4, &mut |_, _| Ok(1.0 / VOLUME), &mut || false)
        .unwrap();
    let duration = 1.0 / 48_000.0;
    let bound = model
        .prepare_time_with_volume_sources(
            duration,
            &[&force],
            &[&heat],
            ThermoviscousTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let original_source = model.relaxation().unwrap();
    let frozen = original_source.frozen_cv_j_kg_k() + gas.specific_gas_constant;
    let reservoirs = original_source
        .processes()
        .iter()
        .map(|p| ThermoviscousCaloricReservoir {
            heat_capacity_j_kg_k: p.heat_capacity_j_kg_k,
            relaxation_time_s: p.relaxation_time_s,
        })
        .collect::<Vec<_>>();
    let original = model
        .inner()
        .prepare_time_with_volume_sources(
            duration,
            frozen,
            &reservoirs,
            None,
            None,
            &[&force],
            &[&heat],
            ThermoviscousTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(bound.frozen_cp_j_kg_k(), frozen);
    assert_eq!(
        bound.caloric_reservoirs().len(),
        original_source.processes().len()
    );
    let a = bound
        .step(
            &bound.zero_state(),
            &[0.2, -0.1],
            &[],
            &[2e-6],
            &[2e-4],
            &mut || false,
        )
        .unwrap();
    let b = original
        .step(
            &original.zero_state(),
            &[0.2, -0.1],
            &[],
            &[2e-6],
            &[2e-4],
            &mut || false,
        )
        .unwrap();
    assert_eq!(a.field.state, b.field.state);
    assert_eq!(a.body_force_work_by_source_j, b.body_force_work_by_source_j);
    assert_eq!(
        a.heat_source_work_by_source_j,
        b.heat_source_work_by_source_j
    );
    assert_eq!(a.field.work.supplied_j(), b.field.work.supplied_j());
    assert!(
        a.field
            .state
            .molecular_temperature_k
            .iter()
            .flatten()
            .any(|v| *v != 0.0)
    );
    assert!(
        model
            .prepare_time_with_volume_sources(
                1.0 / BAND[1],
                &[&force],
                &[&heat],
                ThermoviscousTimeLimits::default(),
                &mut || false,
            )
            .is_err()
    );
}

#[test]
fn material_volume_drives_retain_finite_wall_depth_and_in_plane_transport() {
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
        .unwrap()
        .with_tangentially_isotropic_wall_transport(
            &[MaterialTangentialWallRegion {
                normal: &material,
                tangential: &material,
            }],
            &mut || false,
        )
        .unwrap();
    let force = model
        .volume_force_source(4, &mut |_, _| Ok([1.0 / VOLUME, 0.0, 0.0]), &mut || false)
        .unwrap();
    let heat = model
        .volume_heat_source(4, &mut |_, x| Ok(2.0 * x[0] / LENGTH / VOLUME), &mut || {
            false
        })
        .unwrap();
    assert!(
        model
            .prepare_time_with_volume_sources(
                1e-7,
                &[&force],
                &[&heat],
                ThermoviscousTimeLimits::default(),
                &mut || false,
            )
            .is_err()
    );
    let prepared = model
        .prepare_time_with_wall_diffusion_and_volume_sources(
            1e-7,
            &[2],
            &[&force],
            &[&heat],
            ThermoviscousTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(prepared.conductive_wall_layouts().len(), 1);
    assert_eq!(prepared.maximum_wall_frequency_hz(), Some(BAND[1]));
    let transport = [TangentialWallConductivity {
        isotropic_w_m_k: sample.conductivity_w_m_k(),
    }];
    assert_eq!(
        prepared.tangential_wall_conductivity(),
        Some(transport.as_slice())
    );
    let original = model
        .inner()
        .prepare_time_with_volume_sources(
            1e-7,
            gas.specific_heat_cp,
            &[],
            Some(ThermoviscousWallDiscretization {
                normal_elements: &[2],
                maximum_frequency_hz: BAND[1],
            }),
            Some(&transport),
            &[&force],
            &[&heat],
            ThermoviscousTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let a = prepared
        .step(
            &prepared.zero_state(),
            &[0.0; 2],
            &[],
            &[2e-6],
            &[2e-4],
            &mut || false,
        )
        .unwrap();
    let b = original
        .step(
            &original.zero_state(),
            &[0.0; 2],
            &[],
            &[2e-6],
            &[2e-4],
            &mut || false,
        )
        .unwrap();
    assert_eq!(a.field.state, b.field.state);
    assert_eq!(
        a.field.work.solid_conductive_heat_j,
        b.field.work.solid_conductive_heat_j
    );
    assert!(
        a.field
            .state
            .conductive_wall_temperature_k
            .iter()
            .flatten()
            .any(|v| *v != 0.0)
    );
    assert!(a.field.work.solid_conductive_heat_j > 0.0);
}
