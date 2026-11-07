use super::*;
use crate::charts::TriMesh;
use crate::conductor::{
    CoatingIorSample, CoherentAbsorbingCoating, ConductorDataStatus, ConductorIorSample,
    ConductorOptics, ConductorSource,
};
use crate::dielectric::{BeerLambertAbsorption, CauchyIor, DielectricGlass, GlassProvenance};
use crate::instances::{GeometryInstance, InstanceScene, RigidTransform, SharedGeometry};
use crate::polarization::Stokes;
use crate::polarization::dielectric::scene::{
    Interaction, PathTermination, PolarizedScene, SurfaceBinding, TraceSettings,
};
use crate::polarization::dielectric::transport::PolarizationFrame;
use crate::thermal_surface::heating::{
    SubstrateOpticsModel, feedback::SurfaceConductorLaw, mirror_state::StatefulMirrorReceiver,
};
use fs_blake3::ContentHash;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, Event, MaterialId, StepConfig,
    ThermalLaw, TimedEvent,
};
use std::sync::Arc;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 31,
                kernel_id: 7,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn metal(eta: f64, coated: bool, source_tag: u8) -> ConductorOptics {
    let source = ConductorSource::try_new(
        ContentHash([source_tag; 32]),
        ConductorDataStatus::Representative,
    )
    .unwrap();
    let bare = ConductorOptics::try_new(
        core::array::from_fn(|i| {
            ConductorIorSample::try_new(380.0 + 50.0 * i as f64, eta, 3.0).unwrap()
        }),
        source,
    )
    .unwrap();
    if !coated {
        return bare;
    }
    bare.with_absorbing_coating(
        CoherentAbsorbingCoating::try_new(
            core::array::from_fn(|i| {
                CoatingIorSample::try_new(380.0 + 50.0 * i as f64, 1.45, 0.07).unwrap()
            }),
            100e-9,
            source,
        )
        .unwrap(),
    )
}

fn plane(id: u64, z: f64) -> GeometryInstance {
    GeometryInstance::try_new(
        id,
        fs_blake3::hash_domain(
            "org.frankensim.test.mirror-heating-plane",
            &id.to_le_bytes(),
        ),
        SharedGeometry::mesh(TriMesh::new(
            vec![[-10.0, -10.0, z], [10.0, -10.0, z], [0.0, 10.0, z]],
            vec![[0, 1, 2]],
        )),
        RigidTransform::identity(),
    )
    .unwrap()
}

fn trace_path(
    extinction: f64,
    eta: f64,
    coated: bool,
    launch_power: f64,
    cavity: bool,
    source_tag: u8,
    path_id: u64,
) -> MirrorPathResult {
    trace_path_at_wavelength(
        extinction,
        eta,
        coated,
        launch_power,
        cavity,
        source_tag,
        path_id,
        550.0,
    )
}

#[allow(clippy::too_many_arguments)]
fn trace_path_at_wavelength(
    extinction: f64,
    eta: f64,
    coated: bool,
    launch_power: f64,
    cavity: bool,
    source_tag: u8,
    path_id: u64,
    wavelength_nm: f64,
) -> MirrorPathResult {
    with_cx(|cx| {
        let geometry = InstanceScene::try_new(vec![plane(1, 1.0), plane(2, -1.0)]).unwrap();
        let media = [DielectricGlass::new(
            CauchyIor::try_constant(1.0).unwrap(),
            BeerLambertAbsorption::try_constant(extinction).unwrap(),
            GlassProvenance::Custom,
        )];
        let conductors = [metal(eta, coated, source_tag)];
        let mirror = Interaction::SmoothConductor {
            optics: 0,
            medium: 0,
        };
        let bindings = vec![
            SurfaceBinding {
                object_id: 1,
                interaction: mirror,
            },
            SurfaceBinding {
                object_id: 2,
                interaction: if cavity {
                    mirror
                } else {
                    Interaction::Detector {
                        reference_axis: Vec3::new(1.0, 0.0, 0.0),
                    }
                },
            },
        ];
        let scene =
            PolarizedScene::try_new_with_conductors(&geometry, &media, &conductors, bindings)
                .unwrap();
        let beam = PolarizedBeam::try_new(
            Point3::new(0.0, 0.0, 0.0),
            PolarizationFrame::try_new(Vec3::new(0.0, 0.0, 1.0), Vec3::new(1.0, 0.0, 0.0)).unwrap(),
            wavelength_nm,
            media[0],
            Stokes::unpolarized(launch_power),
        )
        .unwrap();
        scene
            .trace_with_mirror_balances(
                cx,
                beam,
                0,
                TraceSettings {
                    max_surface_events: 3,
                    max_distance_m: 10.0,
                    surface_epsilon_m: 1e-9,
                    query_offset_m: 1e-8,
                },
                31,
                path_id,
            )
            .unwrap()
    })
}

struct CaloricLaw {
    latent: bool,
}
impl ThermalLaw for CaloricLaw {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if !h.is_finite() || !(0.0..=1000.0).contains(&h) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(if !self.latent {
            300.0 + h / 4.0
        } else if h <= 10.0 {
            300.0 + h
        } else if h <= 110.0 {
            310.0
        } else {
            310.0 + h - 110.0
        })
    }
    fn maximum_temperature_slope(&self) -> f64 {
        if self.latent { 1.0 } else { 0.25 }
    }
}
static SENSIBLE: CaloricLaw = CaloricLaw { latent: false };
static LATENT: CaloricLaw = CaloricLaw { latent: true };

fn config() -> StepConfig {
    StepConfig {
        max_dt_s: 3.0,
        energy_tolerance_j: 1e-10,
    }
}

fn thermal(mass_kg: f64, law: &'static dyn ThermalLaw) -> DepositionSimulation<'static> {
    let cell = DepositedCell {
        id: CellId(11),
        segment_id: 11,
        body_id: 1,
        material_id: MaterialId(7),
        mass_kg,
        specific_h_j_kg: 0.0,
        geometry: CellGeometry {
            start_m: [0.0; 3],
            end_m: [1.0, 0.0, 0.0],
            volume_m3: 1.0,
            surface_area_m2: 6.0,
        },
        boundary: Boundary {
            ambient_k: 300.0,
            htc_w_m2_k: 0.0,
            reservoir: None,
            source_w: 0.0,
        },
    };
    DepositionSimulation::new(
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Birth(cell),
        }],
        &[(MaterialId(7), law)],
        DepositionBudget {
            max_cells: 2,
            max_links: 0,
            max_events: 2,
            max_substeps: 64,
        },
    )
    .unwrap()
    .advance_to(0.0, config(), || Ok(()))
    .unwrap()
    .0
}

fn bindings() -> BTreeMap<u64, SurfaceHeatingReceiver> {
    let receiver = SurfaceHeatingReceiver {
        thermal_cell: CellId(11),
        material_id: MaterialId(7),
        temperature_bounds_k: (300.0, 600.0),
    };
    BTreeMap::from([(1, receiver), (2, receiver)])
}

fn limits() -> SurfaceHeatingBudget {
    SurfaceHeatingBudget {
        max_sources: 4,
        max_total_power_segments: 16,
        max_source_evaluations: 64,
    }
}

fn constant_power(watts: f64) -> SpectralPowerSchedule<1> {
    SpectralPowerSchedule::try_new(
        &[SpectralPowerSegment {
            start_s: 0.0,
            end_s: 1.0,
            start_power_w: [watts],
            end_power_w: [watts],
        }],
        1,
    )
    .unwrap()
}

fn heating(
    path: &MirrorPathResult,
    power: &SpectralPowerSchedule<1>,
    thermal: &DepositionSimulation<'_>,
) -> SurfaceOpticalHeating<1> {
    SurfaceOpticalHeating::try_new_from_mirror_path(
        path,
        500,
        power,
        &bindings(),
        thermal,
        limits(),
        SurfaceHeatingTolerance {
            power_roundoff_w: 1e-8,
        },
        || Ok(()),
    )
    .unwrap()
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-12 * expected.abs().max(1e-30),
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

#[test]
fn g1_coated_mirror_path_heats_native_cell_after_beer_attenuation_and_exact_schedule() {
    let power = SpectralPowerSchedule::try_new(
        &[
            SpectralPowerSegment {
                start_s: 0.25,
                end_s: 1.25,
                start_power_w: [4.0],
                end_power_w: [12.0],
            },
            SpectralPowerSegment {
                start_s: 1.25,
                end_s: 1.5,
                start_power_w: [4.0],
                end_power_w: [4.0],
            },
        ],
        2,
    )
    .unwrap();
    // Independent Airy reflectance for n_film=1.45+0.07i, d=100 nm,
    // n_substrate=0.2+3i at 550 nm and normal incidence.
    let r = 0.515_088_324_939_222_6;
    let intercepted_j = 9.0 * fs_math::det::exp(-0.2);
    let heat_j = intercepted_j * (1.0 - r);
    for reference_w in [1.0, 3.0] {
        let path = trace_path(0.2, 0.2, true, reference_w, false, 77, 7);
        assert_eq!(path.source_beam().power().i, reference_w);
        assert_eq!(path.path.termination(), PathTermination::Detected);
        let original = thermal(2.0, &SENSIBLE);
        let input = heating(&path, &power, &original);
        let (heated, accepted, audit) = original
            .advance_observed(&input, 2.0, config(), || Ok(()))
            .unwrap();
        close(accepted.energy().incident_j[0], intercepted_j);
        close(accepted.energy().absorbed_j[0], heat_j);
        close(accepted.energy().thermal_input_j, heat_j);
        close(audit.external_heat_j, heat_j);
        close(heated.cells()[&CellId(11)].specific_h_j_kg, heat_j / 2.0);
        close(
            heated.temperature_k(CellId(11)).unwrap(),
            300.0 + heat_j / 8.0,
        );
        assert_eq!(accepted.source_evaluations(), 4);
        assert_eq!(original.temperature_k(CellId(11)).unwrap(), 300.0);
    }
}

#[test]
fn g1_repeated_mirror_visits_heat_one_cell_into_its_native_latent_plateau() {
    let path = trace_path(0.0, 0.2, false, 4.0, true, 77, 7);
    assert_eq!(path.path.termination(), PathTermination::SurfaceLimit);
    assert_eq!(path.mirrors.len(), 3);
    let original = thermal(1.0, &LATENT);
    let input = heating(&path, &constant_power(250.0), &original);
    assert_eq!(input.source_count(), 3);
    let (heated, accepted, audit) = original
        .advance_observed(&input, 1.0, config(), || Ok(()))
        .unwrap();
    let r = (0.8 * 0.8 + 9.0) / (1.2 * 1.2 + 9.0);
    let heat_j = 250.0 * (1.0 - r * r * r);
    close(heated.cells()[&CellId(11)].specific_h_j_kg, heat_j);
    assert_eq!(heated.temperature_k(CellId(11)).unwrap(), 310.0);
    close(accepted.energy().thermal_input_j, heat_j);
    close(audit.external_heat_j, heat_j);
    assert_eq!(accepted.sources_w().len(), 1);
    close(accepted.sources_w()[&CellId(11)], heat_j);
    // Local incidences include each return; the surviving capped ray stays light.
    close(accepted.energy().incident_j[0], 250.0 * (1.0 + r + r * r));
    close(path.path.beam().power().i, 4.0 * r * r * r);
}

#[test]
fn g1_rounded_unit_mirror_reflection_still_increases_native_enthalpy() {
    let path = trace_path(0.0, 1e-20, false, 1.0, false, 77, 7);
    assert_eq!(path.mirrors[0].balance.reflection().beam.power().i, 1.0);
    let original = thermal(1.0, &SENSIBLE);
    let input = heating(&path, &constant_power(1e12), &original);
    let (heated, accepted, audit) = original
        .advance_observed(&input, 1.0, config(), || Ok(()))
        .unwrap();
    // A=4 eta/((1+eta)^2+k^2)=4e-21, independently of the rounded R.
    let heat_j = 4e-9;
    close(heated.cells()[&CellId(11)].specific_h_j_kg, heat_j);
    assert!(heated.temperature_k(CellId(11)).unwrap() > 300.0);
    close(accepted.energy().absorbed_j[0], heat_j);
    close(audit.external_heat_j, heat_j);
    close(accepted.energy().residual_j[0], -heat_j);
}

#[test]
fn g3_mirror_heating_retains_original_optical_source_and_path_identity() {
    let power = constant_power(2.0);
    let mut identities = Vec::new();
    for (source_tag, path_id) in [(77, 7), (78, 7), (77, 8)] {
        let path = trace_path(0.0, 0.2, false, 3.0, false, source_tag, path_id);
        assert_eq!(path.seed(), 31);
        assert_eq!(path.path_id(), path_id);
        assert_eq!(
            path.mirrors[0].balance.optics_identity(),
            metal(0.2, false, source_tag).content_identity()
        );
        let sources = SurfaceIlluminationSource::from_mirror_path(
            &path,
            500,
            &power,
            &bindings(),
            limits(),
            || Ok(()),
        )
        .unwrap();
        identities.push(sources[0].response.identity());
        close(sources[0].power.peak_power_w()[0], 2.0);
    }
    assert_ne!(identities[0], identities[1]);
    assert_ne!(identities[0], identities[2]);
}

#[test]
fn g0_mirror_thermal_preparation_refuses_zero_normalization_unmapped_cells_and_cancel() {
    let original = thermal(1.0, &SENSIBLE);
    let power = constant_power(2.0);
    let zero = trace_path(0.0, 0.2, false, 0.0, false, 77, 7);
    assert!(
        SurfaceIlluminationSource::from_mirror_path(
            &zero,
            500,
            &power,
            &bindings(),
            limits(),
            || Ok(()),
        )
        .is_err()
    );
    // A later macroscopic schedule must not amplify reference-ray underflow
    // into lost or severely quantized heat. Rescale and retrace these rays.
    for (eta, reference) in [(0.2, f64::from_bits(1)), (1e-20, 1e-300)] {
        let unresolved = trace_path(0.0, eta, false, reference, false, 77, 7);
        assert!(
            SurfaceIlluminationSource::from_mirror_path(
                &unresolved,
                500,
                &power,
                &bindings(),
                limits(),
                || Ok(()),
            )
            .is_err()
        );
    }
    let path = trace_path(0.0, 0.2, false, 1.0, false, 77, 7);
    assert!(
        SurfaceIlluminationSource::from_mirror_path(
            &path,
            500,
            &power,
            &BTreeMap::new(),
            limits(),
            || Ok(()),
        )
        .is_err()
    );
    assert!(matches!(
        SurfaceIlluminationSource::from_mirror_path(
            &path,
            500,
            &power,
            &bindings(),
            limits(),
            || Err(DepositionError::Cancelled),
        ),
        Err(SurfaceHeatingError::Thermal(DepositionError::Cancelled))
    ));
    let mut wrong = bindings();
    wrong.get_mut(&1).unwrap().material_id = MaterialId(999);
    assert!(
        SurfaceOpticalHeating::try_new_from_mirror_path(
            &path,
            500,
            &power,
            &wrong,
            &original,
            limits(),
            SurfaceHeatingTolerance {
                power_roundoff_w: 1e-8
            },
            || Ok(()),
        )
        .is_err()
    );
    assert_eq!(original.time_s(), 0.0);
    assert_eq!(original.cells()[&CellId(11)].specific_h_j_kg, 0.0);
}

#[test]
fn g1_spectral_mirror_histories_share_one_native_temperature_and_preserve_band_fluence() {
    let blue = trace_path_at_wavelength(0.0, 0.2, true, 3.0, true, 77, 10, 450.0);
    let red = trace_path_at_wavelength(0.2, 0.2, true, 7.0, false, 77, 11, 650.0);
    let blue_power = SpectralPowerSchedule::try_new(
        &[SpectralPowerSegment {
            start_s: 0.25,
            end_s: 1.25,
            start_power_w: [0.0],
            end_power_w: [16.0],
        }],
        1,
    )
    .unwrap();
    let red_power = SpectralPowerSchedule::try_new(
        &[
            SpectralPowerSegment {
                start_s: 0.0,
                end_s: 0.5,
                start_power_w: [6.0],
                end_power_w: [6.0],
            },
            SpectralPowerSegment {
                start_s: 1.5,
                end_s: 2.0,
                start_power_w: [10.0],
                end_power_w: [10.0],
            },
        ],
        2,
    )
    .unwrap();
    let original = thermal(2.0, &SENSIBLE);
    // Path order deliberately differs from wavelength order. Blue returns to
    // the same receiver three times; red contributes one attenuated incidence.
    let input = SurfaceOpticalHeating::try_new_from_mirror_paths(
        [450.0, 650.0],
        &[
            MirrorPathIllumination {
                path: &red,
                power: &red_power,
            },
            MirrorPathIllumination {
                path: &blue,
                power: &blue_power,
            },
        ],
        500,
        &bindings(),
        &original,
        limits(),
        SurfaceHeatingTolerance {
            power_roundoff_w: 1e-8,
        },
        || Ok(()),
    )
    .unwrap();
    assert_eq!(input.source_count(), 4);
    let (heated, accepted, audit) = original
        .advance_observed(&input, 2.0, config(), || Ok(()))
        .unwrap();
    // Independent complex Airy formula, n_film=1.45+0.07i, d=100 nm,
    // n_substrate=0.2+3i; phase is exp(4 pi i n_film d/lambda).
    let r_blue = 0.563_260_630_633_973_4;
    let r_red = 0.508_094_446_671_066_5;
    let red_incident = 8.0 * fs_math::det::exp(-0.2);
    let heat = [
        8.0 * (1.0 - r_blue * r_blue * r_blue),
        red_incident * (1.0 - r_red),
    ];
    close(
        accepted.energy().incident_j[0],
        8.0 * (1.0 + r_blue + r_blue * r_blue),
    );
    close(accepted.energy().incident_j[1], red_incident);
    close(
        accepted.energy().reflected_j[0],
        8.0 * r_blue * (1.0 + r_blue + r_blue * r_blue),
    );
    close(accepted.energy().reflected_j[1], red_incident * r_red);
    for (lane, expected) in heat.into_iter().enumerate() {
        close(accepted.energy().absorbed_j[lane], expected);
    }
    let total_heat = heat[0] + heat[1];
    close(accepted.energy().thermal_input_j, total_heat);
    close(audit.external_heat_j, total_heat);
    close(
        heated.cells()[&CellId(11)].specific_h_j_kg,
        total_heat / 2.0,
    );
    close(
        heated.temperature_k(CellId(11)).unwrap(),
        300.0 + total_heat / 8.0,
    );
}

#[test]
fn g3_spectral_path_splitting_preserves_native_heat_without_implicit_averaging() {
    let blue = trace_path_at_wavelength(0.3, 0.2, true, 1.0, false, 77, 10, 450.0);
    let blue_copy = trace_path_at_wavelength(0.3, 0.2, true, 7.0, false, 77, 12, 450.0);
    let red = trace_path_at_wavelength(0.2, 0.2, true, 3.0, false, 77, 11, 650.0);
    let full = constant_power(12.0);
    let one_third = constant_power(4.0);
    let two_thirds = constant_power(8.0);
    let red_power = constant_power(3.0);
    let original = thermal(1.0, &SENSIBLE);
    let run = |paths: &[MirrorPathIllumination<'_>]| {
        let input = SurfaceOpticalHeating::try_new_from_mirror_paths(
            [450.0, 650.0],
            paths,
            500,
            &bindings(),
            &original,
            limits(),
            SurfaceHeatingTolerance {
                power_roundoff_w: 1e-8,
            },
            || Ok(()),
        )
        .unwrap();
        original
            .advance_observed(&input, 1.0, config(), || Ok(()))
            .unwrap()
    };
    let (unsplit, single, _) = run(&[
        MirrorPathIllumination {
            path: &blue,
            power: &full,
        },
        MirrorPathIllumination {
            path: &red,
            power: &red_power,
        },
    ]);
    let (split, divided, _) = run(&[
        MirrorPathIllumination {
            path: &blue,
            power: &one_third,
        },
        MirrorPathIllumination {
            path: &red,
            power: &red_power,
        },
        MirrorPathIllumination {
            path: &blue_copy,
            power: &two_thirds,
        },
    ]);
    for lane in 0..2 {
        close(
            divided.energy().incident_j[lane],
            single.energy().incident_j[lane],
        );
        close(
            divided.energy().absorbed_j[lane],
            single.energy().absorbed_j[lane],
        );
    }
    close(
        split.cells()[&CellId(11)].specific_h_j_kg,
        unsplit.cells()[&CellId(11)].specific_h_j_kg,
    );
    close(
        split.temperature_k(CellId(11)).unwrap(),
        unsplit.temperature_k(CellId(11)).unwrap(),
    );
}

#[test]
fn g0_spectral_mirror_batches_require_exact_lanes_and_aggregate_budgets() {
    let path = trace_path(0.0, 0.2, false, 1.0, false, 77, 7);
    let power = constant_power(2.0);
    let input = MirrorPathIllumination {
        path: &path,
        power: &power,
    };
    for packet in [[450.0, 650.0], [550.0, 550.0]] {
        assert!(
            SurfaceIlluminationSource::from_mirror_paths(
                packet,
                &[input],
                500,
                &bindings(),
                limits(),
                || Ok(()),
            )
            .is_err()
        );
    }
    for budget in [
        SurfaceHeatingBudget {
            max_sources: 1,
            ..limits()
        },
        SurfaceHeatingBudget {
            max_total_power_segments: 1,
            ..limits()
        },
    ] {
        assert!(matches!(
            SurfaceIlluminationSource::from_mirror_paths(
                [550.0],
                &[input, input],
                500,
                &bindings(),
                budget,
                || Ok(()),
            ),
            Err(SurfaceHeatingError::Budget)
        ));
    }
    assert!(
        SurfaceIlluminationSource::from_mirror_paths(
            [550.0],
            &[input, input],
            u64::MAX,
            &bindings(),
            limits(),
            || Ok(()),
        )
        .is_err()
    );
    let admitted = SurfaceIlluminationSource::from_mirror_paths(
        [450.0, 550.0],
        &[input],
        u64::MAX,
        &bindings(),
        limits(),
        || Ok(()),
    )
    .unwrap();
    assert_eq!(admitted[0].source_id, u64::MAX);
    assert_eq!(admitted[0].power.peak_power_w(), [0.0, 2.0]);
    assert_eq!(admitted[0].response.absorbed_fraction()[0], 0.0);
}

fn temperature_mirror_law(high_k: f64, high_eta: f64) -> Arc<SurfaceConductorLaw> {
    Arc::new(
        SurfaceConductorLaw::try_temperature(
            ConductorSource::try_new(ContentHash([91; 32]), ConductorDataStatus::Representative)
                .unwrap(),
            &[
                (300.0, metal(0.2, false, 77)),
                (high_k, metal(high_eta, false, 77)),
            ],
            2,
        )
        .unwrap(),
    )
}

fn law_receivers(law: &Arc<SurfaceConductorLaw>) -> BTreeMap<u64, StatefulMirrorReceiver> {
    bindings()
        .into_iter()
        .map(|(id, receiver)| {
            (
                id,
                StatefulMirrorReceiver {
                    receiver,
                    optics: SubstrateOpticsModel::Law(Arc::clone(law)),
                },
            )
        })
        .collect()
}

fn thermal_pair() -> DepositionSimulation<'static> {
    let first = thermal(1.0, &SENSIBLE).cells()[&CellId(11)].clone();
    let mut second = first.clone();
    second.id = CellId(12);
    second.segment_id = 12;
    DepositionSimulation::new(
        vec![
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(first),
            },
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(second),
            },
        ],
        &[(MaterialId(7), &SENSIBLE as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 2,
            max_links: 0,
            max_events: 2,
            max_substeps: 64,
        },
    )
    .unwrap()
    .advance_to(0.0, config(), || Ok(()))
    .unwrap()
    .0
}

fn normal_metal_r(eta: f64) -> f64 {
    ((eta - 1.0) * (eta - 1.0) + 9.0) / ((eta + 1.0) * (eta + 1.0) + 9.0)
}

#[test]
fn g1_stateful_upstream_mirror_changes_downstream_heat_and_converges_to_coupled_fresnel_ode() {
    let path = trace_path(0.0, 0.2, false, 1.0, true, 77, 7);
    let power = constant_power(250.0);
    let law = temperature_mirror_law(500.0, 6.2);
    let mut receivers = law_receivers(&law);
    receivers.get_mut(&2).unwrap().receiver.thermal_cell = CellId(12);
    receivers.get_mut(&2).unwrap().optics = SubstrateOpticsModel::Fixed(metal(0.2, false, 77));
    // Independent normal-incidence Fresnel ODE, with first/third hits at cell 11
    // and the intervening reflection at cell 12. Both masses are 1 kg and cp=4.
    let rates = |h: [f64; 2]| {
        let r1 = normal_metal_r(0.2 + 0.03 * h[0] / 4.0);
        let r2 = normal_metal_r(0.2);
        [
            250.0 * (1.0 - r1) * (1.0 + r1 * r2),
            250.0 * r1 * (1.0 - r2),
        ]
    };
    let mut reference = [0.0; 2];
    let fine_dt = 1.0 / 2048.0;
    for _ in 0..2048 {
        let a = rates(reference);
        let b = rates(core::array::from_fn(|i| {
            reference[i] + 0.5 * fine_dt * a[i]
        }));
        let c = rates(core::array::from_fn(|i| {
            reference[i] + 0.5 * fine_dt * b[i]
        }));
        let d = rates(core::array::from_fn(|i| reference[i] + fine_dt * c[i]));
        for i in 0..2 {
            reference[i] += fine_dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    let mut errors = Vec::new();
    for steps in [8, 16, 32] {
        let original = thermal_pair();
        let input = SurfaceOpticalHeating::try_new_from_stateful_mirror_paths(
            [550.0],
            &[MirrorPathIllumination {
                path: &path,
                power: &power,
            }],
            500,
            &receivers,
            &original,
            limits(),
            2,
            SurfaceHeatingTolerance {
                power_roundoff_w: 1e-8,
            },
            || Ok(()),
        )
        .unwrap();
        let dt = 1.0 / steps as f64;
        let (heated, accepted, audit) = original
            .advance_observed(
                &input,
                1.0,
                StepConfig {
                    max_dt_s: dt,
                    ..config()
                },
                || Ok(()),
            )
            .unwrap();
        let mut explicit = [0.0; 2];
        for _ in 0..steps {
            let rate = rates(explicit);
            for i in 0..2 {
                explicit[i] += dt * rate[i];
            }
        }
        for (i, id) in [CellId(11), CellId(12)].into_iter().enumerate() {
            close(heated.cells()[&id].specific_h_j_kg, explicit[i]);
        }
        close(audit.external_heat_j, explicit[0] + explicit[1]);
        close(accepted.energy().thermal_input_j, explicit[0] + explicit[1]);
        assert_eq!(accepted.response_evaluations(), 3 * (steps - 1));
        assert!(explicit[1] < 250.0 * normal_metal_r(0.2) * (1.0 - normal_metal_r(0.2)));
        errors.push((explicit[0] - reference[0]).abs() + (explicit[1] - reference[1]).abs());
    }
    assert!(
        errors[0] > 1.8 * errors[1] && errors[1] > 1.8 * errors[2],
        "{errors:?}"
    );
}

struct MirrorPlateau(u8);
impl ThermalLaw for MirrorPlateau {
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        Some([self.0; 32])
    }
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if !h.is_finite() || !(0.0..=1000.0).contains(&h) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(310.0)
    }
    fn maximum_temperature_slope(&self) -> f64 {
        0.0
    }
}
static MIRROR_PLATEAU: MirrorPlateau = MirrorPlateau(61);
static WRONG_MIRROR_PLATEAU: MirrorPlateau = MirrorPlateau(62);

#[test]
fn g1_full_mirror_itinerary_follows_latent_enthalpy_at_unchanged_temperature() {
    let path = trace_path(0.0, 0.2, false, 1.0, true, 77, 7);
    let power = constant_power(250.0);
    let law = Arc::new(
        SurfaceConductorLaw::try_enthalpy(
            ConductorSource::try_new(ContentHash([91; 32]), ConductorDataStatus::Representative)
                .unwrap(),
            ContentHash([61; 32]),
            (310.0, 310.0),
            &[(0.0, metal(0.2, false, 77)), (200.0, metal(2.2, false, 77))],
            2,
        )
        .unwrap(),
    );
    let receivers = law_receivers(&law);
    let original = thermal(1.0, &MIRROR_PLATEAU);
    let input = SurfaceOpticalHeating::try_new_from_stateful_mirror_paths(
        [550.0],
        &[MirrorPathIllumination {
            path: &path,
            power: &power,
        }],
        500,
        &receivers,
        &original,
        limits(),
        2,
        SurfaceHeatingTolerance {
            power_roundoff_w: 1e-8,
        },
        || Ok(()),
    )
    .unwrap();
    let mut state = original.clone();
    let mut observer = input.clone();
    let mut expected_h = 0.0;
    for step in 1..=3 {
        let r = normal_metal_r(0.2 + 0.01 * expected_h);
        expected_h += 0.125 * 250.0 * (1.0 - r * r * r);
        let (next, accepted, _) = state
            .advance_observed(
                &observer,
                step as f64 * 0.125,
                StepConfig {
                    max_dt_s: 0.125,
                    ..config()
                },
                || Ok(()),
            )
            .unwrap();
        close(next.cells()[&CellId(11)].specific_h_j_kg, expected_h);
        assert_eq!(next.temperature_k(CellId(11)).unwrap(), 310.0);
        assert_eq!(
            accepted.response_evaluations(),
            if step == 1 { 0 } else { 3 }
        );
        state = next;
        observer = accepted;
    }
    let resumed = SurfaceOpticalHeating::try_new_from_stateful_mirror_paths(
        [550.0],
        &[MirrorPathIllumination {
            path: &path,
            power: &power,
        }],
        500,
        &receivers,
        &state,
        limits(),
        2,
        SurfaceHeatingTolerance {
            power_roundoff_w: 1e-8,
        },
        || Ok(()),
    )
    .unwrap();
    assert_eq!(input.model_identity(), resumed.model_identity());
    assert!(
        SurfaceOpticalHeating::try_new_from_stateful_mirror_paths(
            [550.0],
            &[MirrorPathIllumination {
                path: &path,
                power: &power
            }],
            500,
            &receivers,
            &thermal(1.0, &WRONG_MIRROR_PLATEAU),
            limits(),
            2,
            SurfaceHeatingTolerance {
                power_roundoff_w: 1e-8
            },
            || Ok(()),
        )
        .is_err()
    );
    assert_eq!(original.cells()[&CellId(11)].specific_h_j_kg, 0.0);
}

#[test]
fn g0_stateful_mirror_paths_refuse_missing_events_and_roll_back_unsupported_native_endpoints() {
    let path = trace_path(0.0, 0.2, false, 1.0, true, 77, 7);
    let power = constant_power(250.0);
    let law = temperature_mirror_law(300.01, 0.21);
    let receivers = law_receivers(&law);
    let original = thermal(1.0, &SENSIBLE);
    let input = SurfaceOpticalHeating::try_new_from_stateful_mirror_paths(
        [550.0],
        &[MirrorPathIllumination {
            path: &path,
            power: &power,
        }],
        500,
        &receivers,
        &original,
        limits(),
        2,
        SurfaceHeatingTolerance {
            power_roundoff_w: 1e-8,
        },
        || Ok(()),
    )
    .unwrap();
    assert!(
        original
            .advance_observed(
                &input,
                0.125,
                StepConfig {
                    max_dt_s: 0.125,
                    ..config()
                },
                || Ok(())
            )
            .is_err()
    );
    assert_eq!(original.time_s(), 0.0);
    assert_eq!(original.cells()[&CellId(11)].specific_h_j_kg, 0.0);
    assert_eq!(input.energy().thermal_input_j, 0.0);
    assert_eq!(input.response_evaluations(), 0);
    let mut missing_event = path.clone();
    missing_event.mirrors.remove(1);
    assert!(
        SurfaceOpticalHeating::try_new_from_stateful_mirror_paths(
            [550.0],
            &[MirrorPathIllumination {
                path: &missing_event,
                power: &power
            }],
            500,
            &receivers,
            &original,
            limits(),
            2,
            SurfaceHeatingTolerance {
                power_roundoff_w: 1e-8
            },
            || Ok(()),
        )
        .is_err()
    );
}

fn folded_normals() -> [Vec3; 3] {
    [
        Vec3::new(0.6, 0.0, -0.8),
        Vec3::new(0.1, 0.8, -0.59),
        Vec3::new(-0.35, 0.2, 0.91),
    ]
    .map(|n| n.scale(1.0 / fs_math::det::sqrt(n.dot(n))))
}

fn trace_folded_path(optics: &[ConductorOptics; 3], polarization: Stokes) -> MirrorPathResult {
    with_cx(|cx| {
        let cross = |a: Vec3, b: Vec3| {
            Vec3::new(
                a.y * b.z - a.z * b.y,
                a.z * b.x - a.x * b.z,
                a.x * b.y - a.y * b.x,
            )
        };
        let mut origin = Point3::new(0.0, 0.0, 0.0);
        let mut direction = Vec3::new(0.0, 0.0, 1.0);
        let mut geometry = Vec::new();
        for (index, normal) in folded_normals().into_iter().enumerate() {
            origin = Point3::new(
                origin.x + direction.x,
                origin.y + direction.y,
                origin.z + direction.z,
            );
            let u = cross(Vec3::new(0.0, 0.0, 1.0), normal);
            let u = u.scale(1.0 / fs_math::det::sqrt(u.dot(u)));
            let v = cross(normal, u);
            let vertices = [[0.1, 0.0], [-0.1, 0.1], [-0.1, -0.1]]
                .map(|[a, b]| {
                    [
                        origin.x + a * u.x + b * v.x,
                        origin.y + a * u.y + b * v.y,
                        origin.z + a * u.z + b * v.z,
                    ]
                })
                .to_vec();
            let object_id = index as u64 + 1;
            geometry.push(
                GeometryInstance::try_new(
                    object_id,
                    fs_blake3::hash_domain(
                        "org.frankensim.test.folded-mirror",
                        &object_id.to_le_bytes(),
                    ),
                    SharedGeometry::mesh(TriMesh::new(vertices, vec![[0, 1, 2]])),
                    RigidTransform::identity(),
                )
                .unwrap(),
            );
            let d = 2.0 * direction.dot(normal);
            direction = Vec3::new(
                direction.x - d * normal.x,
                direction.y - d * normal.y,
                direction.z - d * normal.z,
            );
        }
        let geometry = InstanceScene::try_new(geometry).unwrap();
        let media = [DielectricGlass::new(
            CauchyIor::try_constant(1.0).unwrap(),
            BeerLambertAbsorption::try_constant(0.0).unwrap(),
            GlassProvenance::Custom,
        )];
        let scene = PolarizedScene::try_new_with_conductors(
            &geometry,
            &media,
            optics,
            (0..3)
                .map(|index| SurfaceBinding {
                    object_id: index as u64 + 1,
                    interaction: Interaction::SmoothConductor {
                        optics: index,
                        medium: 0,
                    },
                })
                .collect(),
        )
        .unwrap();
        let beam = PolarizedBeam::try_new(
            Point3::new(0.0, 0.0, 0.0),
            PolarizationFrame::try_new(Vec3::new(0.0, 0.0, 1.0), Vec3::new(1.0, 0.0, 0.0)).unwrap(),
            550.0,
            media[0],
            polarization,
        )
        .unwrap();
        scene
            .trace_with_mirror_balances(
                cx,
                beam,
                0,
                TraceSettings {
                    max_surface_events: 3,
                    max_distance_m: 10.0,
                    surface_epsilon_m: 1e-9,
                    query_offset_m: 1e-8,
                },
                31,
                41,
            )
            .unwrap()
    })
}

#[test]
fn g3_stateful_noncoplanar_mirrors_replay_the_complete_polarized_native_path() {
    let polarization = Stokes::try_new(1.0, 0.1, 0.3, 0.7).unwrap();
    let path = trace_folded_path(&[metal(0.2, false, 77); 3], polarization);
    assert_eq!(path.mirrors.len(), 3);
    let power = constant_power(250.0);
    let law = temperature_mirror_law(500.0, 6.2);
    let mut receivers = law_receivers(&law);
    receivers.insert(3, receivers[&1].clone());
    let mut thermal = thermal(1.0, &SENSIBLE);
    let mut observer = SurfaceOpticalHeating::try_new_from_stateful_mirror_paths(
        [550.0],
        &[MirrorPathIllumination {
            path: &path,
            power: &power,
        }],
        500,
        &receivers,
        &thermal,
        limits(),
        2,
        SurfaceHeatingTolerance {
            power_roundoff_w: 1e-8,
        },
        || Ok(()),
    )
    .unwrap();
    for index in 0..3 {
        let current = law.resolve(&thermal, CellId(11), MaterialId(7)).unwrap();
        let fresh = trace_folded_path(&[current; 3], polarization);
        let expected_heat = 0.05
            * 250.0
            * fresh
                .mirrors
                .iter()
                .map(|hit| hit.balance.absorbed_power_w())
                .sum::<f64>();
        let expected_incident = 0.05
            * 250.0
            * fresh
                .mirrors
                .iter()
                .map(|hit| hit.balance.incident_power_w())
                .sum::<f64>();
        let (next, accepted, _) = thermal
            .advance_observed(
                &observer,
                thermal.time_s() + 0.05,
                StepConfig {
                    max_dt_s: 0.05,
                    ..config()
                },
                || Ok(()),
            )
            .unwrap();
        close(
            next.cells()[&CellId(11)].specific_h_j_kg
                - thermal.cells()[&CellId(11)].specific_h_j_kg,
            expected_heat,
        );
        close(
            accepted.energy().incident_j[0] - observer.energy().incident_j[0],
            expected_incident,
        );
        assert_eq!(
            accepted.response_evaluations(),
            if index == 0 { 0 } else { 3 }
        );
        thermal = next;
        observer = accepted;
    }
}
