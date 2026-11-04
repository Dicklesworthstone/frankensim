//! G1/G3: actual LC outlet flow reaches the existing physical microphone.
use fs_couple::acoustic_realize::{
    AcousticRealizeError, BoreRealization, BoreRealizationReceipt, RealizeOptions,
    duct_observer::BaffledDuctObserver, realize_assembly_with_options,
};
use fs_couple::air_path::GasAbsorptionModel;
use fs_couple::bernoulli_aperture::wall::material::{
    AxialConstraint, LocalMaterialWall, WallMaterialLoss, WallRelaxation,
};
use fs_couple::driving_point::regional::material_wall::MaterialWallLayout;
use fs_couple::pcm_wav::baffled::CircularOutletReceiver;
use fs_material::gas::{GasSpec, GasState};
use fs_scenario::{
    AcousticAssembly, AmbientGas, BeatingReed, CylinderSegment, Listener, ViscothermalDuct,
    WaveguideEnd,
};

const RATE: u32 = 16_000;

// Source-bearing numerical fixture; these coefficients are not measured data.
fn sourced_wall(relaxation_time_s: Option<f64>, band: [f64; 2]) -> LocalMaterialWall {
    use fs_evidence::ValidityDomain;
    use fs_matdb::{
        ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
        PropertyValue, Provenance, QueryPoint, UncertaintyModel,
    };
    use fs_material::state_point::{
        MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
        resolve_material_state_point,
    };
    use fs_qty::{Density, Dims, DynViscosity, Pressure, QuantitySpec};
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, dims, value) in [
        ("density", Density::DIMS, 1000.0),
        ("young_modulus", Pressure::DIMS, 3e6),
        ("poisson_ratio", Dims::NONE, 0.3),
        ("kelvin_voigt_young_viscosity", DynViscosity::DIMS, 100.0),
        ("wall_delta_e", Pressure::DIMS, 2e6),
        (
            "wall_tau",
            Dims([0, 0, 1, 0, 0, 0]),
            relaxation_time_s.unwrap_or(0.001),
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained()
                    .with("T", 290.0, 300.0)
                    .with(
                        "omega",
                        core::f64::consts::TAU * band[0],
                        core::f64::consts::TAU * band[1],
                    ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic local wall for native pressure regression".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_with_quantity(
                name,
                QuantitySpec::dimensional(dims),
                ScalarAdmissibility::Finite,
            )
            .unwrap(),
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic pressure-test wall".into(),
            phase: "solid".into(),
            process: "manufactured test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let material = resolve_material_state_point(
        &card,
        &QueryPoint::new()
            .with("T", 293.15)
            .unwrap()
            .with(
                "omega",
                core::f64::consts::TAU * f64::midpoint(band[0], band[1]),
            )
            .unwrap(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let loss = if relaxation_time_s.is_some() {
        WallMaterialLoss::Maxwell(
            WallRelaxation::from_state_point(
                &material,
                &[("wall_delta_e".into(), "wall_tau".into())],
                band,
                0.5,
                0.5,
            )
            .unwrap(),
        )
    } else {
        WallMaterialLoss::KelvinVoigt
    };
    LocalMaterialWall::new(material, 0.0001, AxialConstraint::StressFree, loss, band).unwrap()
}

fn gas(pressure_pa: f64) -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, pressure_pa).unwrap()
}

fn input(mass_kg: f64) -> AcousticAssembly {
    AcousticAssembly {
        ambient: AmbientGas::sea_level(),
        string: None,
        pluck: None,
        bow: None,
        blow: None,
        soundboard: None,
        body_modes: vec![],
        plate: None,
        cavity: None,
        obstacles: vec![],
        contact_texture: None,
        duct: Some(ViscothermalDuct {
            segments: vec![CylinderSegment {
                radius_m: 0.004,
                outlet_radius_m: 0.004,
                length_m: 0.3,
            }],
            tone_holes: vec![],
            termination: WaveguideEnd::FlangedOpen,
            wall: None,
        }),
        reed: Some(BeatingReed {
            rest_opening_m: 0.0004,
            width_m: 0.01,
            closing_pressure_pa: 10_000.0,
            blowing_pressure_pa: 100.0,
            attack_s: 0.001,
            mass_kg,
            stiffness_n_m: 100.0,
            damping_ratio: 0.1,
        }),
        listener: Listener {
            distance_m: f64::NAN,
        },
        sample_rate_hz: RATE,
        duration_s: 128.0 / f64::from(RATE),
    }
}

fn receiver(gas: GasState, path_samples: usize, radius_m: f64) -> CircularOutletReceiver {
    // With one equal-area annulus every ray has this exact travel length,
    // even though the original finite-aperture receiver remains in use.
    let path = path_samples as f64 * gas.sound_speed / f64::from(RATE);
    CircularOutletReceiver {
        position_m: [0.0, 0.0, (path * path - radius_m * radius_m / 2.0).sqrt()],
        radial_rings: 1,
        angular_points: 8,
        maximum_frequency_hz: 500.0,
    }
}

fn options(gas: GasState, path_samples: usize) -> RealizeOptions<'static> {
    RealizeOptions {
        // Deliberately coarse numerical fixture; its original passive-pitch
        // admission still runs. This test is not a musical accuracy claim.
        bore: BoreRealization::lumped_ladder_with_tolerance(4, 500.0),
        gas: Some(gas),
        observer_absorption: Some(GasAbsorptionModel::Lossless),
        baffled_duct: Some(BaffledDuctObserver {
            segment_gases: None,
            exterior_gas: None,
            material_wall: None,
            receiver: receiver(gas, path_samples, 0.004),
        }),
        ..RealizeOptions::default()
    }
}

#[test]
fn cylinder_and_single_cone_reeds_emit_retarded_pressure_with_inverse_distance_amplitude() {
    for (mass, radius, cells) in [
        (0.0, 0.004, 4),
        (1e-4, 0.004, 4),
        (0.0, 0.006, 8),
        (1e-4, 0.006, 8),
    ] {
        let mut input = input(mass);
        input.duct.as_mut().unwrap().segments[0].outlet_radius_m = radius;
        let gas = gas(101_325.0);
        let observer_options = |samples| {
            let mut result = options(gas, samples);
            result.bore = BoreRealization::lumped_ladder_with_tolerance(cells, 500.0);
            result.baffled_duct.as_mut().unwrap().receiver = receiver(gas, samples, radius);
            result
        };
        let near = realize_assembly_with_options(&input, &observer_options(16)).unwrap();
        let far = realize_assembly_with_options(&input, &observer_options(32)).unwrap();
        assert!(matches!(
            near.bore_receipt,
            Some(BoreRealizationReceipt::LumpedLadder { cells: actual, .. }) if actual == cells
        ));
        let binding = near.duct_observer.as_ref().unwrap();
        assert_eq!(binding.axial_delay_samples, 0);
        assert!(near.pressure_pa[..15].iter().all(|p| *p == 0.0));
        let amplitude = near
            .pressure_pa
            .iter()
            .map(|p| p.abs())
            .fold(0.0_f64, f64::max);
        assert!(amplitude > 0.0 && amplitude.is_finite());
        for i in 0..near.pressure_pa.len() - 16 {
            assert!(
                (far.pressure_pa[i + 16] - 0.5 * near.pressure_pa[i]).abs() < 1e-10 * amplitude
            );
        }
        // The explicitly selected receiver never contains instantaneous local
        // mouth pressure, even though that diagnostic has the same Pa units.
        let mut diagnostic = observer_options(16);
        diagnostic.baffled_duct = None;
        let mut local = input.clone();
        local.listener.distance_m = 1.0;
        let source = realize_assembly_with_options(&local, &diagnostic).unwrap();
        assert!(source.pressure_pa[..15].iter().any(|p| *p != 0.0));
    }
}

#[test]
fn single_cone_pitch_admission_qualifies_the_same_frustum_that_radiates() {
    let gas = gas(101_325.0);
    let mut input = input(0.0);
    input.duct.as_mut().unwrap().segments[0].outlet_radius_m = 0.006;
    input.listener.distance_m = 1.0;
    input.duration_s = 1.0 / f64::from(RATE);
    let mut physical_options = options(gas, 16);
    physical_options.bore = BoreRealization::lumped_ladder_with_tolerance(8, 500.0);
    physical_options.baffled_duct.as_mut().unwrap().receiver = receiver(gas, 16, 0.006);
    let physical = realize_assembly_with_options(&input, &physical_options).unwrap();
    let mut legacy_options = physical_options;
    legacy_options.baffled_duct = None;
    let legacy = realize_assembly_with_options(&input, &legacy_options).unwrap();
    let (physical_error, dispersion) = match physical.bore_receipt.unwrap() {
        BoreRealizationReceipt::LumpedLadder {
            pitch_error_cents,
            dispersion_error_cents,
            ..
        } => (pitch_error_cents.abs(), dispersion_error_cents),
        other => panic!("unexpected physical receipt: {other:?}"),
    };
    let legacy_error = match legacy.bore_receipt.unwrap() {
        BoreRealizationReceipt::LumpedLadder {
            pitch_error_cents, ..
        } => pitch_error_cents.abs(),
        other => panic!("unexpected legacy receipt: {other:?}"),
    };
    // The coarse frustum has a different finite spatial error from the
    // legacy psi-line. A tolerance between the two must reject this physical
    // image even though the old source image passes for the same geometry.
    assert!(physical_error > legacy_error + 1.0);
    let tolerance = f64::midpoint(physical_error, legacy_error);
    assert!(dispersion < tolerance);
    physical_options.bore = BoreRealization::lumped_ladder_with_tolerance(8, tolerance);
    match realize_assembly_with_options(&input, &physical_options).unwrap_err() {
        AcousticRealizeError::LumpedLadderBand {
            pitch_error_cents, ..
        } => {
            assert!(pitch_error_cents.abs() > tolerance);
            assert!((pitch_error_cents.abs() - physical_error).abs() < 1e-9);
        }
        other => panic!("unexpected physical admission refusal: {other:?}"),
    }
    legacy_options.bore = physical_options.bore;
    assert!(realize_assembly_with_options(&input, &legacy_options).is_ok());
}

#[test]
fn pressure_changes_reach_the_same_source_and_microphone_without_pressure_normalization() {
    let input = input(0.0);
    let low = realize_assembly_with_options(&input, &options(gas(70_000.0), 16)).unwrap();
    let high = realize_assembly_with_options(&input, &options(gas(140_000.0), 16)).unwrap();
    assert_eq!(low.gas.pressure, 70_000.0);
    assert_eq!(high.gas.pressure, 140_000.0);
    let difference: f64 = low
        .pressure_pa
        .iter()
        .zip(&high.pressure_pa)
        .map(|(a, b)| (a - b).abs())
        .sum();
    assert!(difference > 1e-8);
}

#[test]
fn warm_bore_uses_its_own_source_gas_and_the_rooms_retarded_microphone() {
    let exterior = gas(101_325.0);
    let interior =
        GasState::try_new(&GasSpec::dry_air_ussa1976(), 310.0, exterior.pressure).unwrap();
    for mass in [0.0, 1e-4] {
        let input = input(mass);
        let observer_options = |path_samples| {
            let mut options = options(interior, path_samples);
            let binding = options.baffled_duct.as_mut().unwrap();
            binding.exterior_gas = Some(exterior);
            binding.receiver = receiver(exterior, path_samples, 0.004);
            options
        };
        let near = realize_assembly_with_options(&input, &observer_options(16)).unwrap();
        let far = realize_assembly_with_options(&input, &observer_options(32)).unwrap();
        let binding = near.duct_observer.as_ref().unwrap();
        assert_eq!(binding.segment_gases, vec![interior]);
        assert_eq!(binding.exterior_gas, exterior);
        assert_eq!(near.gas, exterior);
        assert!(near.pressure_pa[..15].iter().all(|p| *p == 0.0));
        let amplitude = near
            .pressure_pa
            .iter()
            .map(|p| p.abs())
            .fold(0.0_f64, f64::max);
        assert!(amplitude > 0.0 && amplitude.is_finite());
        for i in 0..near.pressure_pa.len() - 16 {
            assert!(
                (far.pressure_pa[i + 16] - 0.5 * near.pressure_pa[i]).abs() < 1e-10 * amplitude
            );
        }
        // Replacing the bore by room gas changes its compression, inertia,
        // losses and Bernoulli source. The exterior override must not perform
        // that replacement implicitly while leaving only a different receipt.
        let uniform_room = realize_assembly_with_options(&input, &options(exterior, 16)).unwrap();
        let difference = near
            .pressure_pa
            .iter()
            .zip(&uniform_room.pressure_pa)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0_f64, f64::max);
        assert!(difference > 1e-8 * amplitude);
    }
}

#[test]
fn explicit_same_exterior_preserves_the_complete_physical_output_bits() {
    let gas = gas(101_325.0);
    let input = input(0.0);
    let implicit = options(gas, 16);
    let mut explicit = implicit;
    explicit.baffled_duct.as_mut().unwrap().exterior_gas = Some(gas);
    let implicit = realize_assembly_with_options(&input, &implicit).unwrap();
    let explicit = realize_assembly_with_options(&input, &explicit).unwrap();
    let bits = |values: &[f64]| {
        values
            .iter()
            .map(|value| value.to_bits())
            .collect::<Vec<_>>()
    };
    assert_eq!(bits(&implicit.pressure_pa), bits(&explicit.pressure_pa));
    assert_eq!(implicit.bore_receipt, explicit.bore_receipt);
    assert_eq!(implicit.duct_observer, explicit.duct_observer);
}

#[test]
fn regional_cylinders_and_frusta_emit_their_loaded_flow_into_the_same_room() {
    let exterior = gas(95_000.0);
    let spec = GasSpec::dry_air_ussa1976();
    let regions = [
        GasState::try_new(&spec, 310.0, exterior.pressure).unwrap(),
        GasState::try_new(&spec, 275.0, exterior.pressure).unwrap(),
    ];
    for (mass, radius) in [(0.0, 0.004), (1e-4, 0.005)] {
        let mut input = input(mass);
        input.duct.as_mut().unwrap().segments = vec![
            CylinderSegment::cylinder(0.004, 0.12),
            CylinderSegment::taper(0.004, radius, 0.18),
        ];
        let mut uniform = options(regions[0], 16);
        uniform.bore = BoreRealization::lumped_ladder_with_tolerance(6, 500.0);
        let binding = uniform.baffled_duct.as_mut().unwrap();
        binding.exterior_gas = Some(exterior);
        binding.receiver = receiver(exterior, 16, radius);
        let regional = RealizeOptions {
            gas: None,
            baffled_duct: Some(BaffledDuctObserver {
                segment_gases: Some(&regions),
                ..uniform.baffled_duct.unwrap()
            }),
            ..uniform
        };
        let local = realize_assembly_with_options(&input, &regional).unwrap();
        let uniform = realize_assembly_with_options(&input, &uniform).unwrap();
        let observed = local.duct_observer.as_ref().unwrap();
        assert_eq!(observed.segment_gases, regions.to_vec());
        assert_eq!(observed.exterior_gas, exterior);
        assert_eq!(observed.axial_delay_samples, 0);
        assert_eq!(local.gas, exterior);
        assert!(
            local.pressure_pa[..15]
                .iter()
                .all(|pressure| *pressure == 0.0)
        );
        let amplitude = local
            .pressure_pa
            .iter()
            .map(|pressure| pressure.abs())
            .fold(0.0_f64, f64::max);
        assert!(amplitude.is_finite() && amplitude > 0.0);
        let difference = local
            .pressure_pa
            .iter()
            .zip(&uniform.pressure_pa)
            .map(|(local, uniform)| (local - uniform).abs())
            .fold(0.0_f64, f64::max);
        assert!(difference > 1e-8 * amplitude);
        match local.bore_receipt.unwrap() {
            BoreRealizationReceipt::LumpedLadder {
                cells,
                pitch_error_cents,
                ..
            } => {
                assert_eq!(cells, 6);
                assert!(pitch_error_cents.is_finite() && pitch_error_cents.abs() <= 500.0);
            }
            other => panic!("unexpected regional receipt: {other:?}"),
        }
    }
}

#[test]
fn identical_regional_profile_preserves_uniform_output_and_receipt_bits() {
    let gas = gas(101_325.0);
    let mut input = input(0.0);
    input.duct.as_mut().unwrap().segments = vec![CylinderSegment::cylinder(0.004, 0.15); 2];
    let regions = [gas; 2];
    let uniform = options(gas, 16);
    let regional = RealizeOptions {
        gas: None,
        baffled_duct: Some(BaffledDuctObserver {
            segment_gases: Some(&regions),
            ..uniform.baffled_duct.unwrap()
        }),
        ..uniform
    };
    let uniform = realize_assembly_with_options(&input, &uniform).unwrap();
    let regional = realize_assembly_with_options(&input, &regional).unwrap();
    assert_eq!(
        uniform
            .pressure_pa
            .iter()
            .map(|value| value.to_bits())
            .collect::<Vec<_>>(),
        regional
            .pressure_pa
            .iter()
            .map(|value| value.to_bits())
            .collect::<Vec<_>>()
    );
    assert_eq!(uniform.bore_receipt, regional.bore_receipt);
    assert_eq!(uniform.duct_observer, regional.duct_observer);
}

#[test]
fn sourced_kelvin_wall_preserves_the_equivalent_original_cylinder_pressure_and_pitch() {
    let gas = gas(101_325.0);
    let source = sourced_wall(None, [1.0, 8000.0]);
    let pin = source.base_patch(0.004, 1.0).unwrap().wall;
    let layout = MaterialWallLayout::Uniform(source);
    let mut numerical_input = input(0.0);
    numerical_input.duct.as_mut().unwrap().wall = Some(fs_scenario::LocallyReactingWall {
        surface_density_kg_m2: pin.surface_density,
        stiffness_pa_per_m: pin.stiffness_per_area,
        resistance_pa_s_per_m: pin.resistance,
    });
    let numerical_options = options(gas, 16);
    let numerical = realize_assembly_with_options(&numerical_input, &numerical_options).unwrap();
    let material_options = RealizeOptions {
        baffled_duct: Some(BaffledDuctObserver {
            material_wall: Some(&layout),
            ..numerical_options.baffled_duct.unwrap()
        }),
        ..numerical_options
    };
    let sourced = realize_assembly_with_options(&input(0.0), &material_options).unwrap();
    let bits = |values: &[f64]| {
        values
            .iter()
            .map(|value| value.to_bits())
            .collect::<Vec<_>>()
    };
    assert_eq!(bits(&numerical.pressure_pa), bits(&sourced.pressure_pa));
    assert_eq!(numerical.bore_receipt, sourced.bore_receipt);
    assert_eq!(
        sourced.duct_observer.unwrap().material_wall,
        Some(layout.clone())
    );
    assert!(
        realize_assembly_with_options(&numerical_input, &material_options).is_err(),
        "the sourced wall cannot be counted again as an authored numerical wall"
    );
}

#[test]
fn regional_maxwell_wall_sources_change_the_same_loaded_frustum_microphone() {
    let exterior = gas(101_325.0);
    let regions = [
        exterior,
        GasState::try_new(&GasSpec::dry_air_ussa1976(), 310.0, exterior.pressure).unwrap(),
    ];
    let walls = MaterialWallLayout::Segments(vec![
        sourced_wall(Some(0.0002), [1.0, 8000.0]),
        sourced_wall(Some(0.002), [1.0, 8000.0]),
    ]);
    let changed = MaterialWallLayout::Segments(vec![
        sourced_wall(Some(0.0002), [1.0, 8000.0]),
        sourced_wall(Some(0.01), [1.0, 8000.0]),
    ]);
    for mass in [0.0, 1e-4] {
        let mut input = input(mass);
        input.duct.as_mut().unwrap().segments = vec![
            CylinderSegment::cylinder(0.004, 0.12),
            CylinderSegment::taper(0.004, 0.005, 0.18),
        ];
        let render = |wall| {
            let settings = RealizeOptions {
                bore: BoreRealization::lumped_ladder_with_tolerance(6, 500.0),
                observer_absorption: Some(GasAbsorptionModel::Lossless),
                baffled_duct: Some(BaffledDuctObserver {
                    segment_gases: Some(&regions),
                    exterior_gas: Some(exterior),
                    material_wall: Some(wall),
                    receiver: receiver(exterior, 16, 0.005),
                }),
                ..RealizeOptions::default()
            };
            realize_assembly_with_options(&input, &settings).unwrap()
        };
        let original = render(&walls);
        let changed = render(&changed);
        let receipt = original.duct_observer.as_ref().unwrap();
        assert_eq!(receipt.material_wall.as_ref(), Some(&walls));
        assert_eq!(receipt.segment_gases, regions.to_vec());
        assert_eq!(receipt.exterior_gas, exterior);
        assert_eq!(receipt.axial_delay_samples, 0);
        assert!(
            original.pressure_pa[..15]
                .iter()
                .all(|pressure| *pressure == 0.0)
        );
        let amplitude = original
            .pressure_pa
            .iter()
            .map(|p| p.abs())
            .fold(0.0_f64, f64::max);
        assert!(amplitude > 0.0 && amplitude.is_finite());
        let difference = original
            .pressure_pa
            .iter()
            .zip(&changed.pressure_pa)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0_f64, f64::max);
        assert!(difference > 1e-8 * amplitude);
    }
}

#[test]
fn original_material_frequency_refusal_survives_physical_lumped_admission() {
    let layout = MaterialWallLayout::Uniform(sourced_wall(Some(0.001), [1.0, 50.0]));
    let uniform = options(gas(101_325.0), 16);
    let settings = RealizeOptions {
        baffled_duct: Some(BaffledDuctObserver {
            material_wall: Some(&layout),
            ..uniform.baffled_duct.unwrap()
        }),
        ..uniform
    };
    match realize_assembly_with_options(&input(0.0), &settings).unwrap_err() {
        AcousticRealizeError::Duct(fs_duct::DuctError::BadParameter { what }) => {
            assert_eq!(
                what,
                "material wall frequency lies outside its declared source band"
            );
        }
        other => panic!("source refusal was replaced: {other:?}"),
    }
}

#[test]
fn unsupported_boundary_and_region_physics_are_not_silently_omitted() {
    let gas = gas(101_325.0);
    let mut input = input(0.0);
    let mut options = options(gas, 16);
    options.baffled_duct.as_mut().unwrap().exterior_gas =
        Some(GasState::try_new(&GasSpec::dry_air_ussa1976(), 310.0, 1.01 * gas.pressure).unwrap());
    assert!(realize_assembly_with_options(&input, &options).is_err());
    options.baffled_duct.as_mut().unwrap().exterior_gas = None;
    input.duct.as_mut().unwrap().termination = WaveguideEnd::Closed;
    assert!(realize_assembly_with_options(&input, &options).is_err());
    let duct = input.duct.as_mut().unwrap();
    duct.termination = WaveguideEnd::FlangedOpen;
    duct.segments = vec![CylinderSegment::cylinder(0.004, 0.15); 2];
    let regions = [
        gas,
        GasState::try_new(&GasSpec::dry_air_ussa1976(), 310.0, 1.01 * gas.pressure).unwrap(),
    ];
    let regional_options = RealizeOptions {
        gas: None,
        baffled_duct: Some(BaffledDuctObserver {
            segment_gases: Some(&regions),
            ..options.baffled_duct.unwrap()
        }),
        ..options
    };
    assert!(realize_assembly_with_options(&input, &regional_options).is_err());
    let matched = [
        gas,
        GasState::try_new(&GasSpec::dry_air_ussa1976(), 310.0, gas.pressure).unwrap(),
    ];
    let underresolved = RealizeOptions {
        bore: BoreRealization::lumped_ladder_with_tolerance(3, 500.0),
        baffled_duct: Some(BaffledDuctObserver {
            segment_gases: Some(&matched),
            ..regional_options.baffled_duct.unwrap()
        }),
        ..regional_options
    };
    match realize_assembly_with_options(&input, &underresolved).unwrap_err() {
        AcousticRealizeError::InvalidDescription { what } => assert!(what.contains("two cells")),
        other => panic!("unexpected regional cell-budget refusal: {other:?}"),
    }
}
