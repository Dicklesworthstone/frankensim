//! G1/G3: original tetrahedral material data drives the ordinary DG microphone.
//! Manufactured geometry/laws do not establish measured instrument realism.
mod surface;
use super::*;
use crate::acoustic_realize::AcousticRealizeError;
use crate::bernoulli_aperture::dg::loaded::performance::LoadedDgPerformance;
use crate::bernoulli_aperture::performance::file::material_modal::dg::MaterialModalDgValvePerformance;

const INPUT: &str =
    include_str!("../../../../../../examples/material_modal_valve/solid-dg.performance");

fn source() -> String {
    INPUT.replace("audio 96000 96000 2048 1", "audio 96000 96000 512 1")
}

fn solid_pack(
    rho: f64,
    young: f64,
    c33: f64,
    missing: &str,
    high: f64,
) -> NormalizedMaterialCardPack {
    solid_pack_with_strain(rho, young, c33, missing, high, 0.02)
}
fn solid_pack_with_strain(
    rho: f64,
    young: f64,
    c33: f64,
    missing: &str,
    high: f64,
    strain_limit: f64,
) -> NormalizedMaterialCardPack {
    solid_pack_with_memory(rho, young, c33, missing, high, strain_limit, 0.0001)
}
fn solid_pack_with_memory(
    rho: f64,
    young: f64,
    c33: f64,
    missing: &str,
    high: f64,
    strain_limit: f64,
    time: f64,
) -> NormalizedMaterialCardPack {
    let mut claims = ClaimSet::new();
    let mut identity = [
        rho.to_le_bytes(),
        young.to_le_bytes(),
        c33.to_le_bytes(),
        high.to_le_bytes(),
        strain_limit.to_le_bytes(),
        time.to_le_bytes(),
    ]
    .concat();
    identity.extend_from_slice(missing.as_bytes());
    let artifact = hash_domain("org.frankensim.test.native-solid-valve", &identity);
    let provenance = Provenance {
        source: "manufactured 3D material valve".into(),
        license: "CC0-1.0".into(),
        artifact: None,
    };
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: "manufactured two-coordinate tetrahedron".into(),
            method: "declared physical regression".into(),
            artifact,
            caveats: "not measured cane or tissue".into(),
            provenance: provenance.clone(),
        })
        .unwrap();
    for (name, value, dims) in [
        ("density", rho, Density::DIMS),
        ("young_modulus", young, Pressure::DIMS),
        ("poisson_ratio", 0.3, Dims::NONE),
        ("strain_limit", strain_limit, Dims::NONE),
        ("relaxation_time", time, Dims([0, 0, 1, 0, 0, 0])),
        ("branch_11", 100.0, Pressure::DIMS),
        ("branch_22", 200.0, Pressure::DIMS),
        ("branch_33", c33, Pressure::DIMS),
        ("zero_modulus", 0.0, Pressure::DIMS),
    ] {
        if name == missing {
            continue;
        }
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained()
                    .with("T", 280.0, 320.0)
                    .with("omega", 0.0, high),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                provenance: provenance.clone(),
                observations: vec![observation.clone()],
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured".into(),
            phase: "soft-solid".into(),
            process: "reference".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let normalized = NormalizedPack::new(
        "native-solid-valve",
        "synthetic-v1",
        artifact,
        "CC0-1.0",
        card.claims().clone(),
        vec![],
        vec![],
    )
    .unwrap();
    NormalizedMaterialCardPack::new(card.id().clone(), normalized).unwrap()
}

fn build(
    text: &str,
    block: usize,
    rho: f64,
    young: f64,
    c33: f64,
) -> MaterialModalDgValvePerformance {
    MaterialModalDgValvePerformance::from_bytes_with_sources(
        text.as_bytes(),
        block,
        &CancelGate::new_clock_free(),
        None,
        |_| Ok(solid_pack(rho, young, c33, "", 100_000.0)),
    )
    .unwrap()
}
fn valve(source: &MaterialModalDgValvePerformance) -> &ModalAperture {
    let LoadedAperture::Modal(valve) = source.renderer().system().aperture() else {
        panic!("the original solid must keep the ordinary modal aperture");
    };
    valve
}
fn render(source: MaterialModalDgValvePerformance) -> (Vec<f64>, LoadedDgPerformance) {
    let samples = source.info().samples as usize;
    let block = source.renderer().config().max_block;
    let mut renderer = source.into_renderer();
    let mut pressure = vec![0.0; samples];
    for chunk in pressure.chunks_mut(block) {
        renderer.block(chunk).unwrap();
    }
    (pressure, renderer)
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() < tolerance * b.abs().max(1e-25),
        "{a:e} != {b:e}"
    );
}

#[test]
fn complete_solid_native_source_retains_tensor_data_and_drives_the_original_microphone() {
    let text = source();
    let mut calls = 0;
    let built = MaterialModalDgValvePerformance::from_bytes_with_sources(
        text.as_bytes(),
        37,
        &CancelGate::new_clock_free(),
        None,
        |path| {
            assert_eq!(path, "solid-valve.matpack");
            calls += 1;
            Ok(solid_pack(1000.0, 20000.0, 400.0, "", 100000.0))
        },
    )
    .unwrap();
    assert_eq!(calls, 1);
    assert_eq!(
        (
            built.info().nodes,
            built.info().modes,
            built.info().triangles
        ),
        (4, 2, 0)
    );
    assert!(built.material_records().is_empty());
    let record = built.solid_record().unwrap();
    assert_eq!((record.element_order(), record.tetrahedra()), ("tet4", 1));
    assert_eq!(record.materials().len(), 1);
    assert_eq!(
        record.materials()[0]
            .resolved()
            .property("branch_33")
            .unwrap()
            .value_si(),
        400.0
    );
    let mechanics = valve(&built);
    assert!(mechanics.material_plates().is_empty());
    assert!(mechanics.material_relaxation().is_none());
    assert_eq!(mechanics.blade_propagator(), BladePropagator::ExactZoh);
    let solid = mechanics.material_solid().unwrap();
    assert_eq!(solid.source().modes().len(), 2);
    assert_eq!(solid.source().source_id(), "manufactured-source-soft-solid");
    assert_eq!(solid.source().approximation().reference_omega_rad_s, 2000.0);
    let full_c = mechanics.viscosity_matrix_n_s_m().unwrap();
    assert_eq!(full_c.len(), 4);
    assert!(
        full_c[1].abs() > 1e-12,
        "source tensor coupling must not be diagonalized away"
    );
    assert_eq!(
        mechanics.gas(),
        Some(built.renderer().system().load().duct().inlet_gas())
    );

    let (reference, run) = render(built);
    assert!(reference.iter().all(|v| v.is_finite()));
    assert!(reference.iter().any(|v| v.abs() > 1e-15));
    assert_eq!(run.system().load().state().step_count(), 512);
    assert_eq!(run.applied_plate_force_controls(0).len(), 2);
    assert_eq!(run.held_generalized_forces_n(), [0.0, 0.0]);
    assert_eq!(run.samples_rendered(), 512);
    for (rho, young, c33) in [
        (1150.0, 20000.0, 400.0),
        (1000.0, 26000.0, 400.0),
        (1000.0, 20000.0, 2000.0),
    ] {
        let (changed, changed_run) = render(build(&text, 37, rho, young, c33));
        assert!(
            reference
                .iter()
                .zip(&changed)
                .any(|(a, b)| (a - b).abs() > 1e-15)
        );
        assert_ne!(
            run.system().load().state(),
            changed_run.system().load().state()
        );
        assert_ne!(
            run.last_frame().unwrap().aperture.bore_pressure_pa,
            changed_run.last_frame().unwrap().aperture.bore_pressure_pa
        );
    }
    for block in [1, 512] {
        let (other, other_run) = render(build(&text, block, 1000.0, 20000.0, 400.0));
        assert_eq!(reference, other);
        assert_eq!(
            run.system().load().state(),
            other_run.system().load().state()
        );
        assert_eq!(
            encode_pcm16_wav(&reference, 96000, 1.0).unwrap(),
            encode_pcm16_wav(&other, 96000, 1.0).unwrap()
        );
    }
}

#[test]
fn source_pressure_load_and_complete_volume_image_match_an_independent_tetrahedral_limit() {
    let mut text = source();
    let start = text.find("solid_branches 1\n").unwrap();
    let end = text[start..].find("\nsolid_nodes").unwrap() + start;
    text.replace_range(start..end, "solid_branches 0");
    text = text
        .replace("solid_node 0 0.01 0 1 0 1", "solid_node 0 0.01 0 1 1 1")
        .replace("initial 2\nstate 0 0\nstate 0 0", "initial 1\nstate 0 0")
        .replace(
            "source_forces 2\nsource_force 0 0 0.00001\nsource_force 192 0 0\n",
            "",
        );
    let built = build(&text, 37, 1000.0, 20000.0, 400.0);
    let source = valve(&built).material_solid().unwrap().source();
    assert_eq!(source.modes().len(), 1);
    // Single free u_z at vertex (0,0,a). Its affine strain is epsilon_zz=u/a.
    // Exact consistent diagonal mass is rho*V/10; K=C_3333*V/a^2.
    // The opposing exterior face pressure has signed nodal area a^2/6.
    let a = 0.01_f64;
    let volume = a.powi(3) / 6.0;
    let free_mass = 1000.0 * volume / 10.0;
    let c33 = 20000.0 * (1.0 - 0.3) / ((1.0 + 0.3) * (1.0 - 2.0 * 0.3));
    let stiffness = c33 * volume / a.powi(2);
    let pressure_area = a.powi(2) / 6.0;
    close(source.total_mass_kg(), 1000.0 * volume, 1e-12);
    close(source.modes()[0].lambda, stiffness / free_mass, 1e-10);
    close(
        source.pressure_projection_m2_per_sqrt_kg()[0].abs(),
        pressure_area / free_mass.sqrt(),
        1e-10,
    );
    close(
        source.opening_projection_per_sqrt_kg()[0].abs(),
        1.0 / (3.0 * free_mass.sqrt()),
        1e-10,
    );
    close(
        source.compliance().exact_volume_m3_per_pa,
        pressure_area.powi(2) / stiffness,
        1e-10,
    );
    close(
        source.compliance().retained_volume_m3_per_pa,
        pressure_area.powi(2) / stiffness,
        1e-10,
    );
    // The source map is the opening-positive conjugate reference. Positive
    // upstream-minus-bore pressure applies its NEGATIVE in the aperture.
    let mechanics = valve(&built);
    let static_coordinate_per_pa: Vec<_> = mechanics
        .spec()
        .modes
        .iter()
        .map(|mode| -mode.pressure_area_m2 / mode.stiffness_n_m)
        .collect();
    let nodal = mechanics
        .material_solid()
        .unwrap()
        .expand_nodal(&static_coordinate_per_pa, &CancelGate::new_clock_free())
        .unwrap();
    close(nodal[3][2], -pressure_area / stiffness, 1e-10);
    assert!(
        nodal[3][2] < 0.0,
        "positive pressure must close this opening-positive source"
    );
    let conjugate_volume: f64 = mechanics
        .spec()
        .modes
        .iter()
        .zip(&static_coordinate_per_pa)
        .map(|(mode, q)| mode.pressure_area_m2 * q)
        .sum();
    close(conjugate_volume, pressure_area * nodal[3][2], 1e-10);
}

#[test]
fn native_solid_refuses_missing_sources_unresolved_images_and_nonphysical_pressure_faces() {
    let text = source();
    for (missing, band) in [("density", 100000.0), ("branch_33", 100000.0), ("", 1000.0)] {
        assert!(
            MaterialModalDgValvePerformance::from_bytes_with_sources(
                text.as_bytes(),
                37,
                &CancelGate::new_clock_free(),
                None,
                |_| Ok(solid_pack(1000.0, 20000.0, 400.0, missing, band)),
            )
            .is_err()
        );
    }
    for changed in [
        text.replace("solid tet4", "solid tet10"),
        text.replace(
            "solid_axes 1 0 0 0 1 0 0 0 1",
            "solid_axes 2 0 0 0 1 0 0 0 1",
        ),
        text.replace("solid_limits 4 1 2 1e-6 6", "solid_limits 4 1 1 1e-6 6"),
        text.replace("solid_limits 4 1 2 1e-6 6", "solid_limits 4 1 2 1e-6 5"),
        text.replace(
            "solid_pressure_face 1 2 3 -1",
            "solid_pressure_face 0 0 3 -1",
        ),
        text.replace(
            "solid_pressure_face 1 2 3 -1",
            "solid_pressure_face 1 2 3 0.5",
        ),
        text.replace(
            "solid_pressure_faces 1\nsolid_pressure_face 1 2 3 -1",
            "solid_pressure_faces 2\nsolid_pressure_face 1 2 3 -1\nsolid_pressure_face 3 1 2 -1",
        ),
        text.replace("valve 0.0004", "valve 0.0005"),
        text.replace(
            "solid_approximation 2000 0.02 0.02",
            "solid_approximation 2000 0 0",
        ),
        text.replace("initial 2\nstate 0 0\nstate 0 0", "initial 1\nstate 0 0"),
        text.replace("source_forces 2", "wetted_surfaces 0\nsource_forces 2"),
        text.replace(
            "solid_material solid-valve.matpack isotropic",
            "solid_material solid-valve.matpack orthotropic",
        ),
    ] {
        assert!(
            MaterialModalDgValvePerformance::from_bytes_with_sources(
                changed.as_bytes(),
                37,
                &CancelGate::new_clock_free(),
                None,
                |_| Ok(solid_pack(1000.0, 20000.0, 400.0, "", 100000.0)),
            )
            .is_err(),
            "must refuse {changed}"
        );
    }
    let gate = CancelGate::new_clock_free();
    gate.request();
    assert!(matches!(
        MaterialModalDgValvePerformance::from_bytes_with_sources(
            text.as_bytes(),
            37,
            &gate,
            None,
            |_| -> Result<NormalizedMaterialCardPack, String> {
                panic!("cancel before loading source")
            },
        ),
        Err(PlateValveInputError::Physics(
            AcousticRealizeError::Cancelled
        ))
    ));
}

#[test]
fn source_small_strain_domain_rejects_an_invalid_native_initial_shape() {
    let text = source().replacen("state 0 0", "state 0.000001 0", 1);
    assert!(
        MaterialModalDgValvePerformance::from_bytes_with_sources(
            text.as_bytes(),
            37,
            &CancelGate::new_clock_free(),
            None,
            |_| Ok(solid_pack_with_strain(
                1000.0, 20000.0, 400.0, "", 100000.0, 1e-8
            )),
        )
        .is_err()
    );
    // The very same coordinates are admitted when the original material
    // explicitly supports their complete three-dimensional strain.
    build(&text, 37, 1000.0, 20000.0, 400.0);
}

const HEREDITARY: &str =
    include_str!("../../../../../../examples/material_modal_valve/solid-hereditary-dg.performance");
fn hereditary_source() -> String {
    HEREDITARY.replace("audio 96000 96000 2048 1", "audio 96000 96000 512 1")
}

#[test]
fn full_volume_tensor_histories_reach_the_same_native_pressure_and_force_trajectory() {
    let text = hereditary_source();
    let built = build(&text, 37, 1000.0, 20000.0, 400.0);
    let mechanics = valve(&built);
    assert!(mechanics.material_solid().is_none());
    let source = mechanics.hereditary_solid().unwrap().source();
    assert_eq!(source.modes().len(), 2);
    assert_eq!(source.material_points().len(), 1);
    assert_eq!(source.point_elements(), &[0]);
    assert_eq!(source.memory_channel_count(), 2);
    assert_eq!(source.branches().len(), 1);
    assert_eq!(source.branches()[0].relaxation_time_s(), 0.0001);
    assert!(source.basis_error() <= 1e-10);
    assert!(source.instantaneous_frequency_bound_rad_s() / 96000.0 < 0.5);
    assert!(
        mechanics
            .viscosity_matrix_n_s_m()
            .unwrap()
            .iter()
            .all(|c| *c == 0.0),
        "Maxwell histories must not also become frozen Kelvin loss"
    );
    let (reference, run) = render(built);
    assert!(reference.iter().any(|p| p.abs() > 1e-15));
    let LoadedAperture::Modal(final_mechanics) = run.system().aperture() else {
        panic!()
    };
    assert!(
        final_mechanics
            .material_relaxation()
            .unwrap()
            .memory_sqrt_j()
            .iter()
            .flatten()
            .any(|q| *q != 0.0)
    );
    assert!(run.last_frame().unwrap().aperture.relaxation_loss_j > 0.0);
    assert_eq!(run.applied_plate_force_controls(0).len(), 2);
    for (c33, time) in [(2000.0, 0.0001), (400.0, 0.0002)] {
        let changed = MaterialModalDgValvePerformance::from_bytes_with_sources(
            text.as_bytes(),
            37,
            &CancelGate::new_clock_free(),
            None,
            |_| {
                Ok(solid_pack_with_memory(
                    1000.0, 20000.0, c33, "", 100000.0, 0.02, time,
                ))
            },
        )
        .unwrap();
        let (pressure, other) = render(changed);
        assert!(
            reference
                .iter()
                .zip(&pressure)
                .any(|(a, b)| (a - b).abs() > 1e-15)
        );
        assert_ne!(run.system().load().state(), other.system().load().state());
    }
    for block in [1, 512] {
        let (other, other_run) = render(build(&text, block, 1000.0, 20000.0, 400.0));
        assert_eq!(reference, other);
        assert_eq!(
            run.system().load().state(),
            other_run.system().load().state()
        );
        assert_eq!(
            encode_pcm16_wav(&reference, 96000, 1.0).unwrap(),
            encode_pcm16_wav(&other, 96000, 1.0).unwrap()
        );
    }
}

#[test]
fn native_volume_memory_initialization_and_original_projection_limits_are_explicit() {
    let text = hereditary_source().replacen("state 0 0", "state 0.000001 0", 1);
    let relaxed = build(&text, 37, 1000.0, 20000.0, 400.0);
    let unrelaxed = build(
        &text.replace("solid_memory relaxed", "solid_memory unrelaxed"),
        37,
        1000.0,
        20000.0,
        400.0,
    );
    assert_eq!(
        valve(&relaxed).modal_states(),
        valve(&unrelaxed).modal_states()
    );
    assert!(valve(&unrelaxed).stored_energy_j() > valve(&relaxed).stored_energy_j());
    assert_ne!(render(relaxed).0, render(unrelaxed).0);
    let selection = "solid_memory relaxed 0.5 4 1 1048576 1048576 1e-10";
    for change in [
        "solid_memory relaxed 0.0001 4 1 1048576 1048576 1e-10",
        "solid_memory relaxed 0.5 3 1 1048576 1048576 1e-10",
        "solid_memory relaxed 0.5 4 0 1048576 1048576 1e-10",
        "solid_memory relaxed 0.5 4 1 0 1048576 1e-10",
        "solid_memory relaxed 0.5 4 1 1048576 0 1e-10",
    ] {
        assert!(
            MaterialModalDgValvePerformance::from_bytes_with_sources(
                text.replace(selection, change).as_bytes(),
                37,
                &CancelGate::new_clock_free(),
                None,
                |_| Ok(solid_pack(1000.0, 20000.0, 400.0, "", 100000.0)),
            )
            .is_err()
        );
    }
    for (strain_limit, time) in [(1e-8, 0.0001), (0.02, 1e-7)] {
        assert!(
            MaterialModalDgValvePerformance::from_bytes_with_sources(
                text.as_bytes(),
                37,
                &CancelGate::new_clock_free(),
                None,
                |_| Ok(solid_pack_with_memory(
                    1000.0,
                    20000.0,
                    400.0,
                    "",
                    100000.0,
                    strain_limit,
                    time
                )),
            )
            .is_err()
        );
    }
}
