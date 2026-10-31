//! The ordinary assembly observes the same physical bore and material media.
use fs_couple::acoustic_realize::{
    BoreRealization, RealizeOptions, duct_observer::BaffledDuctObserver,
    realize_assembly_with_options,
};
use fs_couple::air_path::{GasAbsorptionModel, ObserverFilter};
use fs_couple::driving_point::regional::RegionalDuct;
use fs_couple::pcm_wav::baffled::CircularOutletReceiver;
use fs_couple::render::ReedBoreVoice;
use fs_couple::thin_plate::PlateBank;
use fs_duct::{Duct, Segment, Termination};
use fs_material::gas::{GasSpec, GasState};
use fs_material::gas_species::GasSpecies;
use fs_scenario::{
    AcousticAssembly, AmbientGas, BeatingReed, CylinderSegment, Listener, ToneHole,
    ViscothermalDuct, VolumeVelocityPulse, WaveguideEnd,
};

use fs_couple::bernoulli_aperture::wall::material::{
    AxialConstraint, LocalMaterialWall, WallMaterialLoss, relaxation::WallRelaxation,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use fs_qty::{Density, Dims, Pressure};

fn wall_source(young_pa: f64, relaxation_s: Option<f64>, band_hz: [f64; 2])
    -> (MaterialCard, QueryPoint, Vec<ScalarPropertyRequirement>)
{
    let mut properties = vec![
        ("density", Density::DIMS, 1000.0),
        ("young_modulus", Pressure::DIMS, young_pa),
        ("poisson_ratio", Dims([0; 6]), 0.3),
    ];
    if let Some(time) = relaxation_s {
        properties.extend([
            ("wall_delta_e", Pressure::DIMS, 2e7),
            ("wall_tau", Dims([0, 0, 1, 0, 0, 0]), time),
        ]);
    }
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, dims, value) in properties {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained().with("T", 290.0, 310.0)
                .with("omega", core::f64::consts::TAU * band_hz[0],
                    core::f64::consts::TAU * band_hz[1]),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured wall for integration verification, not measured material".into(),
                license: "CC0-1.0".into(), artifact: None,
            },
        }).unwrap();
        requirements.push(ScalarPropertyRequirement::try_new(
            name, dims, ScalarAdmissibility::Finite,
        ).unwrap());
    }
    let card = MaterialCard::assemble(MaterialStateId {
        chemistry: "manufactured wall".into(), phase: "solid".into(),
        process: "numerical fixture".into(), revision: 0,
    }, claims, vec![]).unwrap();
    (card, QueryPoint::new().with("T", 300.0).unwrap()
        .with("omega", 1000.0).unwrap(), requirements)
}

fn material_wall(young_pa: f64, relaxation_s: Option<f64>, band_hz: [f64; 2]) -> LocalMaterialWall {
    let (card, point, requirements) = wall_source(young_pa, relaxation_s, band_hz);
    let source = resolve_material_state_point(&card, &point,
        &requirements, MaterialPropertySelection::SingleClaimOnly).unwrap();
    let loss = if relaxation_s.is_some() {
        WallMaterialLoss::Maxwell(WallRelaxation::from_state_point(&source,
            &[("wall_delta_e".into(), "wall_tau".into())], band_hz, 0.5, 0.5).unwrap())
    } else { WallMaterialLoss::Elastic };
    LocalMaterialWall::new(source, 0.00005, AxialConstraint::StressFree, loss, band_hz).unwrap()
}

/// G1/G3: the common compiler resolves the same original wall source as the
/// direct physical path. Changing only card values changes the actual output.
#[test]
fn common_material_compiler_carries_elastic_and_relaxing_walls_into_pressure() {
    use fs_couple::material_assembly::{
        AcousticMaterialBindings, MaterialSource, compile_material_assembly,
        duct::{DuctWallMaterialBinding, WallConstitutiveSource, WallMaterialLossSource},
    };
    let template = input();
    let gases = [gas(293.15), gas(310.15)];
    let observer = binding(Some(&gases), Some(gas(283.15)));
    let band = [1e-6, 4500.0];
    for variants in [
        [(3e7, None), (6e7, None)],
        [(3e7, Some(0.0003)), (3e7, Some(0.003))],
    ] {
        let mut pressures = Vec::new();
        for (young, relaxation) in variants {
            let (card, point, requirements) = wall_source(young, relaxation, band);
            let selection = MaterialPropertySelection::PinnedByProperty(
                requirements.iter().map(|requirement| {
                    (requirement.name().to_owned(), card.claims_for(requirement.name())[0].0)
                }).collect(),
            );
            let loss = if relaxation.is_some() {
                WallMaterialLossSource::Maxwell {
                    properties: vec![("wall_delta_e".into(), "wall_tau".into())],
                    max_dt_over_tau: 0.5,
                    max_angular_step: 0.5,
                }
            } else { WallMaterialLossSource::Elastic };
            let compiled = compile_material_assembly(&template, &AcousticMaterialBindings {
                duct_wall: Some(DuctWallMaterialBinding {
                    source: MaterialSource { card: &card, point: &point, selection },
                    constitutive: WallConstitutiveSource::Isotropic(loss),
                    thickness_m: 0.00005,
                    axial: AxialConstraint::StressFree,
                    band_hz: band,
                }),
                ..AcousticMaterialBindings::default()
            }).unwrap();
            assert_eq!(compiled.assembly(), &template);
            let wall = compiled.duct_wall().unwrap();
            assert_eq!(wall.material().query(), point);
            assert_eq!(wall.material().properties().len(), requirements.len());
            let actual = compiled.realize_with_options(&options(observer)).unwrap();
            assert_eq!(actual.duct_observer.as_ref().unwrap().material_wall.as_ref(), Some(wall));
            assert_eq!(actual.duct_observer.as_ref().unwrap().segment_gases, gases);
            assert_eq!(actual.gas, observer.exterior_gas.unwrap());
            let direct_wall = material_wall(young, relaxation, band);
            let direct = realize_assembly_with_options(&template, &options(BaffledDuctObserver {
                material_wall: Some(&direct_wall), ..observer
            })).unwrap();
            assert_eq!(actual.pressure_pa, direct.pressure_pa);
            assert_eq!(actual, compiled.realize_with_options(&options(observer)).unwrap());
            pressures.push(actual.pressure_pa);
        }
        assert_ne!(pressures[0], pressures[1]);
    }
}

/// G4: no authored override, unavailable source or unrouted material may be
/// silently substituted. Failed cold compilation leaves prior state untouched.
#[test]
fn common_material_compiler_refuses_ambiguous_unrouted_and_invalid_wall_sources() {
    use fs_couple::material_assembly::{
        AcousticMaterialBindings, MaterialSource, compile_material_assembly,
        duct::{DuctWallMaterialBinding, WallConstitutiveSource, WallMaterialLossSource},
    };
    let template = input();
    let band = [1e-6, 4500.0];
    let (card, point, _) = wall_source(3e7, None, band);
    let binding = DuctWallMaterialBinding {
        source: MaterialSource { card: &card, point: &point,
            selection: MaterialPropertySelection::SingleClaimOnly },
        constitutive: WallConstitutiveSource::Isotropic(WallMaterialLossSource::Elastic),
        thickness_m: 0.00005,
        axial: AxialConstraint::StressFree,
        band_hz: band,
    };
    let compile = |input: &AcousticAssembly, wall: DuctWallMaterialBinding<'_>| {
        compile_material_assembly(input, &AcousticMaterialBindings {
            duct_wall: Some(wall), ..AcousticMaterialBindings::default()
        })
    };
    let compiled = compile(&template, binding.clone()).unwrap();
    let original_wall = compiled.duct_wall().unwrap().clone();
    assert!(compiled.realize().is_err(), "material wall must reach a physical observer");
    let mut absent = template.clone();
    absent.duct = None;
    assert!(compile(&absent, binding.clone()).is_err());
    let mut authored = template.clone();
    authored.duct.as_mut().unwrap().wall = Some(fs_scenario::LocallyReactingWall {
        surface_density_kg_m2: 1.0, stiffness_pa_per_m: 1e8, resistance_pa_s_per_m: 0.0,
    });
    assert!(compile(&authored, binding.clone()).is_err());
    for wall in [
        DuctWallMaterialBinding { thickness_m: 0.001, ..binding.clone() },
        DuctWallMaterialBinding { band_hz: [1e-6, 9000.0], ..binding.clone() },
        DuctWallMaterialBinding {
            constitutive: WallConstitutiveSource::Isotropic(WallMaterialLossSource::Maxwell {
                properties: vec![("wall_delta_e".into(), "wall_tau".into())],
                max_dt_over_tau: 0.5, max_angular_step: 0.5,
            }), ..binding.clone()
        },
        DuctWallMaterialBinding {
            constitutive: WallConstitutiveSource::Isotropic(WallMaterialLossSource::Maxwell {
                properties: vec![("wall_delta_e".into(), "wall_tau".into()); 9],
                max_dt_over_tau: 0.5, max_angular_step: 0.5,
            }), ..binding.clone()
        },
    ] {
        assert!(compile(&template, wall).is_err());
    }
    let observer = BaffledDuctObserver {
        material_wall: Some(&original_wall), ..self::binding(None, None)
    };
    assert!(compiled.realize_with_options(&options(observer)).is_err());
    assert_eq!(compiled.duct_wall(), Some(&original_wall));
    assert_eq!(compiled.assembly(), &template);
}

fn gas(temperature: f64) -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), temperature, 100_000.0).unwrap()
}
fn input() -> AcousticAssembly {
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
            segments: [0.2, 0.3]
                .map(|length_m| CylinderSegment {
                    radius_m: 0.0015,
                    outlet_radius_m: 0.0015,
                    length_m,
                })
                .to_vec(),
            tone_holes: vec![],
            termination: WaveguideEnd::FlangedOpen,
            wall: None,
        }),
        reed: Some(BeatingReed {
            rest_opening_m: 0.0004,
            width_m: 0.013,
            closing_pressure_pa: 6000.0,
            blowing_pressure_pa: 2800.0,
            attack_s: 0.002,
            mass_kg: 0.0,
            stiffness_n_m: 0.0,
            damping_ratio: 0.35,
        }),
        listener: Listener { distance_m: 1.0 },
        sample_rate_hz: 8_000,
        duration_s: 0.064,
    }
}
fn binding(gases: Option<&[GasState]>, exterior: Option<GasState>) -> BaffledDuctObserver<'_> {
    BaffledDuctObserver {
        segment_gases: gases,
        exterior_gas: exterior,
        material_wall: None,
        receiver: CircularOutletReceiver {
            position_m: [0.0, 0.0, 0.5],
            radial_rings: 4,
            angular_points: 16,
            maximum_frequency_hz: 500.0,
        },
    }
}
fn options(binding: BaffledDuctObserver<'_>) -> RealizeOptions<'_> {
    RealizeOptions {
        baffled_duct: Some(binding),
        observer_absorption: Some(GasAbsorptionModel::Lossless),
        ..RealizeOptions::default()
    }
}

/// G1/G3: compare the public assembly against the existing pressure/flow and
/// spatial receiver owners, then independently check the earliest causal arrival.
#[test]
fn physical_assembly_matches_original_owners_and_retains_propagation_delay() {
    let mut input = input();
    // Explicit material and receiver inputs must not consult these unused defaults.
    input.ambient.temperature_k = -1.0;
    input.listener.distance_m = f64::NAN;
    let gases = [gas(293.15), gas(310.15)];
    let binding = binding(Some(&gases), Some(gas(283.15)));
    let actual = realize_assembly_with_options(&input, &options(binding)).unwrap();
    let info = actual.duct_observer.as_ref().unwrap();
    assert_eq!(info.segment_gases, gases);
    assert_eq!(actual.gas, binding.exterior_gas.unwrap());
    let physical = Duct {
        segments: [0.2, 0.3]
            .map(|length| Segment::Cylinder {
                radius: 0.0015,
                length,
            })
            .to_vec(),
    };
    let regional = RegionalDuct::new(&physical, &gases)
        .unwrap()
        .with_exterior_gas(actual.gas)
        .unwrap();
    let mut voice = ReedBoreVoice::new_with_regional(
        &regional,
        input.reed.unwrap(),
        Termination::FlangedOpen,
        PlateBank::default(),
        0.5,
        input.sample_rate_hz,
        actual.pressure_pa.len(),
        None,
    )
    .unwrap();
    let mut microphone = regional
        .baffled_outlet(
            input.sample_rate_hz,
            voice.line_realization().n_fft,
            binding.receiver,
            None,
        )
        .unwrap();
    microphone.step(voice.outgoing_pressure_pa()).unwrap();
    for &pressure in &actual.pressure_pa {
        voice.step_block(&mut [0.0]).unwrap();
        assert_eq!(
            pressure,
            microphone.step(voice.outgoing_pressure_pa()).unwrap()
        );
    }
    let first_possible = info.axial_delay_samples + info.exterior_delay_samples.0 - 1;
    assert!(
        actual.pressure_pa[..first_possible]
            .iter()
            .all(|p| *p == 0.0)
    );
    assert!(actual.pressure_pa.iter().all(|p| p.is_finite()));
    assert!(actual.pressure_pa.iter().any(|p| *p != 0.0));
}

/// G3/G5: exterior material changes cannot replace the authored interior or
/// its travel time. Interior temperature/composition reaches actual propagation.
#[test]
fn physical_assembly_keeps_interior_exterior_and_observer_law_independent() {
    let input = input();
    let gases = [gas(293.15); 2];
    let homogeneous =
        realize_assembly_with_options(&input, &options(binding(Some(&gases), None))).unwrap();
    let same =
        realize_assembly_with_options(&input, &options(binding(Some(&gases), Some(gases[1]))))
            .unwrap();
    assert_eq!(homogeneous, same);
    let exterior = gas(323.15);
    let outside =
        realize_assembly_with_options(&input, &options(binding(Some(&gases), Some(exterior))))
            .unwrap();
    assert_ne!(homogeneous.pressure_pa, outside.pressure_pa);
    assert_eq!(outside.duct_observer.as_ref().unwrap().segment_gases, gases);
    assert_eq!(
        homogeneous
            .duct_observer
            .as_ref()
            .unwrap()
            .axial_delay_samples,
        outside.duct_observer.as_ref().unwrap().axial_delay_samples
    );
    let helium = GasState::try_new(&GasSpecies::Helium.spec(), 293.15, 100_000.0).unwrap();
    let other_inside = [gases[0], helium];
    let inside = realize_assembly_with_options(
        &input,
        &options(binding(Some(&other_inside), Some(exterior))),
    )
    .unwrap();
    assert_eq!(inside.gas, outside.gas);
    assert_ne!(inside.pressure_pa, outside.pressure_pa);
    assert!(
        inside.duct_observer.as_ref().unwrap().axial_delay_samples
            < outside.duct_observer.as_ref().unwrap().axial_delay_samples
    );
    let mut selected = options(binding(Some(&other_inside), Some(exterior)));
    selected.observer_absorption = Some(GasAbsorptionModel::StokesKirchhoff);
    let absorbed = realize_assembly_with_options(&input, &selected).unwrap();
    let original = ObserverFilter::stokes_kirchhoff(&exterior, 0.5, input.sample_rate_hz).unwrap();
    assert_eq!(
        absorbed.duct_observer.as_ref().unwrap().filter,
        *original.receipt()
    );
    let receipt = absorbed.duct_observer.as_ref().unwrap().filter;
    assert!(
        (receipt.receiver_delay_s + receipt.low_frequency_delay_s - 0.5 / exterior.sound_speed)
            .abs()
            < 1e-15
    );
    assert_ne!(absorbed.pressure_pa, inside.pressure_pa);
}

/// G0: unsupported participants, incoherent snapshots and conflicting source
/// assignments refuse rather than becoming a partial or differently loaded render.
#[test]
fn physical_assembly_refuses_incomplete_or_inconsistent_physical_bindings() {
    let input = input();
    let gases = [gas(293.15); 2];
    let selected = options(binding(Some(&gases), None));
    for altered in [
        RealizeOptions {
            observer_absorption: None,
            ..selected
        },
        RealizeOptions {
            gas: Some(gases[0]),
            ..selected
        },
        RealizeOptions {
            bore: BoreRealization::lumped_ladder(16),
            ..selected
        },
        RealizeOptions {
            observer_absorption: Some(GasAbsorptionModel::MolecularRelaxation),
            ..selected
        },
    ] {
        assert!(realize_assembly_with_options(&input, &altered).is_err());
    }
    let mut pressure_mismatch = gases;
    pressure_mismatch[1] =
        GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 80_000.0).unwrap();
    let mut stale = gases;
    stale[1].sound_speed *= 1.1;
    for medium in [&gases[..1], &pressure_mismatch, &stale] {
        assert!(
            realize_assembly_with_options(&input, &options(binding(Some(medium), None))).is_err()
        );
    }
    let mut coupled = input.clone();
    coupled.blow = Some(VolumeVelocityPulse {
        peak_m3_s: 1e-5,
        duration_s: 0.01,
    });
    assert!(realize_assembly_with_options(&coupled, &selected).is_err());
    let mut unflanged = input.clone();
    unflanged.duct.as_mut().unwrap().termination = WaveguideEnd::UnflangedOpen;
    assert!(realize_assembly_with_options(&unflanged, &selected).is_err());
    let mut branched = input;
    branched.duct.as_mut().unwrap().tone_holes.push(ToneHole {
        after_segment: 0,
        radius_m: 0.0005,
        chimney_m: 0.001,
        open_fraction: 1.0,
    });
    assert!(realize_assembly_with_options(&branched, &selected).is_err());
}

/// G1/G3: the source material's thin-cylinder reduction matches an independently
/// calculated constant wall, while a data-only stiffness swap changes sound.
#[test]
fn material_wall_reaches_feedback_and_physical_observation() {
    let input = input();
    let gases = [gas(293.15); 2];
    let wall = material_wall(3e7, None, [1e-6, 4500.0]);
    let bound = BaffledDuctObserver { material_wall: Some(&wall),
        ..binding(Some(&gases), None) };
    let actual = realize_assembly_with_options(&input, &options(bound)).unwrap();
    assert_eq!(actual.duct_observer.as_ref().unwrap().material_wall.as_ref(), Some(&wall));
    assert_eq!(actual.duct_observer.as_ref().unwrap().segment_gases, gases);
    // Independent membrane limit: sigma = rho*h, K = E*h/r^2 at free axial stress.
    let mut reference = input.clone();
    reference.duct.as_mut().unwrap().wall = Some(fs_scenario::acoustic::LocallyReactingWall {
        surface_density_kg_m2: 1000.0 * 0.00005,
        stiffness_pa_per_m: 3e7 * 0.00005 / (0.0015_f64 * 0.0015),
        resistance_pa_s_per_m: 0.0,
    });
    let expected = realize_assembly_with_options(&reference,
        &options(binding(Some(&gases), None))).unwrap();
    let peak = expected.pressure_pa.iter().map(|p| p.abs()).fold(0.0, f64::max);
    assert!(peak > 0.0);
    let error = actual.pressure_pa.iter().zip(&expected.pressure_pa)
        .map(|(a, b)| (a-b).abs()).fold(0.0, f64::max);
    assert!(error / peak < 1e-8, "relative pressure discrepancy {}", error / peak);
    let stiff = material_wall(6e7, None, [1e-6, 4500.0]);
    let changed = realize_assembly_with_options(&input,
        &options(BaffledDuctObserver { material_wall: Some(&stiff), ..bound })).unwrap();
    assert_ne!(actual.pressure_pa, changed.pressure_pa);
    assert_eq!(actual.gas, changed.gas);
    assert_ne!(wall.material().identity(), stiff.material().identity());
}

/// G3/G5: a relaxation-time edit reaches the complete material law; it is not
/// reduced to the identical equilibrium E/rho or added as an output filter.
#[test]
fn material_wall_preserves_hereditary_spectrum_and_independent_gases() {
    let input = input();
    let gases = [gas(293.15), gas(305.15)];
    let exterior = gas(283.15);
    let fast = material_wall(3e7, Some(0.0003), [1e-6, 4500.0]);
    let slow = material_wall(3e7, Some(0.003), [1e-6, 4500.0]);
    let render = |wall: &LocalMaterialWall| {
        realize_assembly_with_options(&input, &options(BaffledDuctObserver {
            material_wall: Some(wall), ..binding(Some(&gases), Some(exterior))
        })).unwrap()
    };
    let first = render(&fast);
    assert_eq!(first, render(&fast));
    let second = render(&slow);
    assert_ne!(first.pressure_pa, second.pressure_pa);
    assert_eq!(first.gas, exterior);
    let retained = first.duct_observer.unwrap();
    assert_eq!(retained.segment_gases, gases);
    assert_eq!(retained.material_wall.as_ref(), Some(&fast));
    // A wall's original 300 K source is independent of both gas temperatures.
    assert_eq!(fast.material().query().axes()["T"], 300.0);
    assert_eq!(retained.exterior_gas, exterior);
}

/// G0: a material binding cannot silently replace another wall, extrapolate
/// its frequency support, or bypass thin-wall geometry admission.
#[test]
fn material_wall_refuses_conflicting_and_unsupported_bindings() {
    let mut input = input();
    let gases = [gas(293.15); 2];
    let wall = material_wall(3e7, None, [1e-6, 4500.0]);
    let bound = BaffledDuctObserver { material_wall: Some(&wall),
        ..binding(Some(&gases), None) };
    input.duct.as_mut().unwrap().wall = Some(fs_scenario::acoustic::LocallyReactingWall {
        surface_density_kg_m2: 1.0, stiffness_pa_per_m: 1e8, resistance_pa_s_per_m: 1.0,
    });
    assert!(realize_assembly_with_options(&input, &options(bound)).is_err());
    input.duct.as_mut().unwrap().wall = None;
    let narrow = material_wall(3e7, None, [1e-6, 1000.0]);
    assert!(realize_assembly_with_options(&input,
        &options(BaffledDuctObserver { material_wall: Some(&narrow), ..bound })).is_err());
    let duct = input.duct.as_mut().unwrap();
    duct.segments[0].radius_m = 0.0005;
    duct.segments[0].outlet_radius_m = 0.0005;
    assert!(realize_assembly_with_options(&input, &options(bound)).is_err());
}
