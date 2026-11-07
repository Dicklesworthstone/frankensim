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
use fs_blake3::ContentHash;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, Event, MaterialId, StepConfig,
    ThermalLaw, TimedEvent,
};

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
