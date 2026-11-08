//! The original reduced material, its interface and resolved fluid share one state.
use super::*;
use crate::shear_shell::coating::{
    CoatingSide, ShellCoating, relaxation::CoatingRelaxationOptions,
};
use crate::shear_shell::constant_loss::relaxation::{
    SubstrateRelaxationOptions, SubstrateStorageReference,
};
use fs_material::coating::ResolvedCoating;

fn film() -> ResolvedCoating {
    // A manufactured one-arm Zener source. Its relaxation time is the exact
    // midpoint of the original one-arm lowering ladder, so this test needs no
    // approximate fractional fit or measured-material claim.
    let mut properties = ClaimSet::new();
    for (name, dims, value) in [
        ("density", Density::DIMS, 1100.0),
        ("varnish_e0", Pressure::DIMS, 1e6),
        ("varnish_einf", Pressure::DIMS, 2e6),
        ("varnish_alpha", Dims::NONE, 1.0),
        ("varnish_beta", Dims::NONE, 1.0),
        (
            "varnish_tau",
            Dims([0, 0, 1, 0, 0, 0]),
            1.0 / 1e7_f64.sqrt(),
        ),
    ] {
        properties
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained()
                    .with("T", 293.15, 293.15)
                    .with("omega", 100.0, 100_000.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic one-arm shell film regression".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic Zener film".into(),
            phase: "solid".into(),
            process: "physical regression".into(),
            revision: 0,
        },
        properties,
        vec![],
    )
    .unwrap();
    let state = resolve_material_state_point(
        &card,
        &QueryPoint::new()
            .with("T", 293.15)
            .unwrap()
            .with("omega", 1000.0)
            .unwrap(),
        &ResolvedCoating::requirements().unwrap(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    ResolvedCoating::from_material(state, 0.3, 0.001).unwrap()
}

fn memory_image(permutation: [usize; 5], coated: bool, lower_film: bool) -> MaterialShellImage {
    let gate = CancelGate::new_clock_free();
    let (mut shell, directors, options, ports) =
        source(permutation, ShellFace::Negative, true, &gate);
    if coated {
        shell = shell
            .with_coatings(
                vec![ShellCoating {
                    region: "panel".into(),
                    material: film(),
                    side: CoatingSide::Positive,
                    thickness_m: 1e-6,
                    reference_omega_rad_s: 1000.0,
                    max_storage_relative_drift: 0.9,
                }],
                &gate,
            )
            .unwrap();
        if lower_film {
            shell = shell
                .with_coating_relaxation(
                    CoatingRelaxationOptions {
                        terms: 1,
                        storage_relative_tolerance: 1e-6,
                        loss_relative_tolerance: 1e-6,
                        max_memory_coordinates: 8,
                    },
                    &gate,
                )
                .unwrap();
        }
    }
    shell
        .with_substrate_relaxation(
            SubstrateRelaxationOptions {
                storage_reference: SubstrateStorageReference::Relaxed,
                max_memory_coordinates: 16,
            },
            &gate,
        )
        .unwrap()
        .reduce(&directors, &options, &ports, &gate)
        .unwrap()
}

fn original_channels(image: &MaterialShellImage) -> (Vec<Vec<f64>>, Vec<f64>) {
    let coating = image.coating_relaxation();
    let substrate = image.substrate_relaxation().unwrap();
    let rows = coating
        .into_iter()
        .flat_map(|m| m.coupling_rows_per_s())
        .chain(substrate.coupling_rows_per_s())
        .cloned()
        .collect();
    let times = coating
        .into_iter()
        .flat_map(|m| m.relaxation_times_s())
        .chain(substrate.relaxation_times_s())
        .copied()
        .collect();
    (rows, times)
}

fn scalar_close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 2e-9 * expected.abs().max(1e-30),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn source_coating_and_substrate_histories_match_native_motion_and_preload_heat() {
    let gate = CancelGate::new_clock_free();
    let image = memory_image([0, 1, 2, 3, 4], true, true);
    let (mesh, positions) = fluid([0, 1, 2, 3, 4, 5]);
    let (rows, times) = original_channels(&image);
    assert!(
        !image
            .coating_relaxation()
            .unwrap()
            .coupling_rows_per_s()
            .is_empty()
    );
    assert!(
        !image
            .substrate_relaxation()
            .unwrap()
            .coupling_rows_per_s()
            .is_empty()
    );
    let source = assembly(&image, &mesh, &positions)
        .bind_gas(gas(), BAND, &gate)
        .unwrap();
    let actual = source.interior().inner();
    assert_eq!(
        actual.moving_wall_relaxation_channels_per_s(),
        Some(rows.as_slice())
    );
    assert_eq!(
        actual.moving_wall_relaxation_times_s(),
        Some(times.as_slice())
    );
    for (native, original) in actual
        .moving_wall()
        .unwrap()
        .structure
        .modes()
        .iter()
        .zip(image.radiators())
    {
        assert_eq!(native.angular_frequency_rad_s, original.omega);
    }
    let c = actual.moving_wall_damping_matrix_per_s().unwrap();
    let interface = image.source().interfaces()[0].law().damping_n_s_m;
    let column = &image.interface_columns()[0];
    for i in 0..2 {
        for j in 0..2 {
            scalar_close(
                c[2 * i + j],
                image.damping().matrix_per_s()[2 * i + j] + interface * column[i] * column[j],
            );
        }
    }
    let native = manual(&image, &mesh, &positions)
        .with_moving_wall_relaxation_channels(
            rows.clone(),
            times,
            controls().relaxation,
            &mut || false,
        )
        .unwrap();
    let (field, force) = response(&source, &gate);
    let expected = native
        .prepare(1000.0, &mut || false)
        .unwrap()
        .response_with_wall_forces(&[], &force, &mut || false)
        .unwrap();
    close(&field.pressure_pa, &expected.pressure_pa, 2e-9);
    close(
        &field.moving_wall.modal_velocity_m_sqrt_kg_per_s,
        &expected.moving_wall.modal_velocity_m_sqrt_kg_per_s,
        2e-9,
    );
    close(
        &field.moving_wall.relaxation_coordinates_sqrt_j,
        &expected.moving_wall.relaxation_coordinates_sqrt_j,
        2e-9,
    );
    scalar_close(
        field.power.structural_relaxation_w,
        expected.power.structural_relaxation_w,
    );
    assert!(field.power.structural_relaxation_w > 0.0 && field.power.viscous_w > 0.0);
    assert!(field.power.balance_residual_w().abs() < 2e-8 * field.power.input_w);

    let time = source
        .interior()
        .prepare_time(1e-5, ThermoviscousTimeLimits::default(), &mut || false)
        .unwrap();
    let reference = native
        .prepare_time(1e-5, ThermoviscousTimeLimits::default(), &mut || false)
        .unwrap();
    let mut preload = time.zero_state();
    assert_eq!(preload.wall_relaxation_coordinates_sqrt_j.len(), rows.len());
    for (i, value) in preload
        .wall_relaxation_coordinates_sqrt_j
        .iter_mut()
        .enumerate()
    {
        *value = if i % 2 == 0 { 1.0 } else { -1.0 } * (i + 1) as f64 * 1e-8;
    }
    let initial_energy: f64 = preload
        .wall_relaxation_coordinates_sqrt_j
        .iter()
        .map(|x| 0.5 * x * x)
        .sum();
    scalar_close(time.stored_energy_j(&preload).unwrap(), initial_energy);
    let original = preload.clone();
    let step = time
        .step(&preload, &[], &[0.0, 0.0], &mut || false)
        .unwrap();
    let expected = reference
        .step(&preload, &[], &[0.0, 0.0], &mut || false)
        .unwrap();
    assert_eq!(preload, original);
    for (&actual, &expected) in step
        .state
        .wall_relaxation_coordinates_sqrt_j
        .iter()
        .zip(&expected.state.wall_relaxation_coordinates_sqrt_j)
    {
        scalar_close(actual, expected);
    }
    scalar_close(
        step.work.structural_relaxation_heat_j,
        expected.work.structural_relaxation_heat_j,
    );
    assert!(step.work.structural_relaxation_heat_j > 0.0);
    assert!(step.work.balance_residual_j().abs() < 2e-8 * initial_energy);
}

#[test]
fn signed_material_histories_preserve_relabelled_physical_response() {
    let gate = CancelGate::new_clock_free();
    let original = memory_image([0, 1, 2, 3, 4], true, true);
    let relabelled = memory_image([3, 0, 4, 1, 2], true, true);
    let (mesh, positions) = fluid([0, 1, 2, 3, 4, 5]);
    let a = assembly(&original, &mesh, &positions)
        .bind_gas(gas(), BAND, &gate)
        .unwrap();
    let b = assembly(&relabelled, &mesh, &positions)
        .bind_gas(gas(), BAND, &gate)
        .unwrap();
    let (first, _) = response(&a, &gate);
    let (second, _) = response(&b, &gate);
    close(&first.pressure_pa, &second.pressure_pa, 2e-8);
    scalar_close(first.power.input_w, second.power.input_w);
    scalar_close(
        first.power.structural_relaxation_w,
        second.power.structural_relaxation_w,
    );
}

#[test]
fn source_memory_support_and_complete_resource_limits_remain_authoritative() {
    let gate = CancelGate::new_clock_free();
    let image = memory_image([0, 1, 2, 3, 4], true, true);
    let (mesh, positions) = fluid([0, 1, 2, 3, 4, 5]);
    let medium = source_medium(&gas(), 0.0).unwrap();
    let thermo = ThermoviscousThermodynamics::ideal_gas(medium).unwrap();
    let build = |image: &MaterialShellImage, controls| {
        ShellWallAssembly::new(
            image,
            geometry(&mesh, &positions, medium, thermo),
            &map(image, &positions),
            medium,
            thermo,
            &boundary(&mesh),
            controls,
            &gate,
        )
        .map(|_| ())
    };
    for quantity in 0..3 {
        let mut bounded = controls();
        match quantity {
            0 => bounded.relaxation.maximum_channels = original_channels(&image).0.len() - 1,
            1 => bounded.relaxation.maximum_workspace_bytes = 1,
            _ => bounded.relaxation.maximum_projection_terms = 0,
        }
        assert!(matches!(
            build(&image, bounded),
            Err(ShellWallError::Budget { .. })
        ));
    }
    // 1 Hz is admitted by substrate/interface data, but not by this film's
    // explicit 100 rad/s lower support. A passive continuation does not widen it.
    assert!(
        assembly(&image, &mesh, &positions)
            .bind_gas(gas(), [1.0, 1000.0], &gate)
            .is_err()
    );
    let slice = memory_image([0, 1, 2, 3, 4], true, false);
    assert!(build(&slice, controls()).is_err());
    let substrate_only = memory_image([0, 1, 2, 3, 4], false, false);
    let source = assembly(&substrate_only, &mesh, &positions)
        .bind_gas(gas(), BAND, &gate)
        .unwrap();
    assert_eq!(
        source
            .interior()
            .inner()
            .moving_wall_relaxation_channels_per_s()
            .unwrap(),
        substrate_only
            .substrate_relaxation()
            .unwrap()
            .coupling_rows_per_s()
    );
}
