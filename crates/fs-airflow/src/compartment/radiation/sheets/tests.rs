use super::super::{AbsorptionCoefficient, OpticalPath, OpticalSegment, STEFAN_BOLTZMANN_W_M2_K4};
use super::*;
use fs_conduction::radiation::{
    EMISSIVITY_DIMS, SURFACE_EMISSIVITY_PROPERTY, ViewFactorEvidence, ViewFactorMatrix,
    ViewFactorTolerance,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, SelectionPolicy, UncertaintyModel,
};

use super::super::soot::{RayleighSoot, SootOpticalDomain};
use super::super::spectral::{
    BandSurfaceEmissivity, SootPopulation, SpectralBand, SpectralOptics, SpectralSolve,
};
use fs_conduction::radiation::BlackbodyBand;

use fs_qty::{Area, Density, Length};

fn emissivity(value: f64) -> SurfaceEmissivity {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS),
            value: PropertyValue::Scalar {
                value,
                dims: EMISSIVITY_DIMS,
            },
            validity: ValidityDomain::unconstrained().with("T", 200.0, 1500.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: Vec::new(),
            provenance: Provenance {
                source: "analytic opaque bath".into(),
                license: "internal-test-use".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "fixture".into(),
            phase: "solid".into(),
            process: "opaque bath".into(),
            revision: 0,
        },
        claims,
        Vec::new(),
    )
    .unwrap();
    SurfaceEmissivity::from_card("bath", &card, 300.0, SelectionPolicy::SingleClaimOnly).unwrap()
}

fn sheet(a: [f64; 2], t: f64) -> SheetPair {
    SheetPair {
        surfaces: [1, 2],
        optics: ReciprocalSheetOptics::new(a, t, [200.0, 1500.0], "analytic diffuse sheet")
            .unwrap(),
    }
}

fn absorption(k: [f64; 2]) -> [LayerAbsorption; 2] {
    k.map(|value| {
        LayerAbsorption::new(
            AbsorptionCoefficient::new(value),
            "explicit uniform gas absorption",
        )
        .unwrap()
    })
}

fn geometry(areas: [f64; 2], gas: bool) -> ParticipatingEnclosure {
    let factors = vec![
        vec![0.0, 1.0, 0.0, 0.0],
        vec![1.0, 0.0, 0.0, 0.0],
        vec![0.0, 0.0, 0.0, 1.0],
        vec![0.0, 0.0, 1.0, 0.0],
    ];
    let matrix = ViewFactorMatrix::admit(
        vec![areas[0], areas[0], areas[1], areas[1]],
        factors,
        ViewFactorEvidence::Analytic {
            geometry: "two opposed black-bath/sheet pairs".into(),
        },
        ViewFactorTolerance::default(),
    )
    .unwrap();
    let paths = [[0, 1], [2, 3]]
        .into_iter()
        .enumerate()
        .map(|(side, surfaces)| OpticalPath {
            surfaces,
            exchange_area: Area::new(areas[side]),
            segments: if gas {
                vec![OpticalSegment {
                    layer: if side == 0 {
                        Layer::Lower
                    } else {
                        Layer::Upper
                    },
                    length: Length::new(if side == 0 { 1.0 } else { 2.0 }),
                }]
            } else {
                Vec::new()
            },
        })
        .collect();
    ParticipatingEnclosure::new(
        matrix,
        paths,
        "declared paired-chamber paths",
        1e-12,
        || false,
    )
    .unwrap()
}

fn blackbody(t: f64) -> f64 {
    STEFAN_BOLTZMANN_W_M2_K4 * t.powi(4)
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-11 * expected.abs().max(1.0),
        "{actual:e} != {expected:e}"
    );
}

fn run(model: &ParticipatingEnclosure, pair: SheetPair, temperatures: [f64; 4]) -> RadiationReport {
    model
        .evaluate_with_sheets(
            &temperatures.map(Temperature::new),
            &[emissivity(1.0), emissivity(1.0)],
            &[pair],
            [Temperature::new(300.0); 2],
            &absorption([0.0; 2]),
            RadiationSolve::default(),
            &mut model.workspace().unwrap(),
            || false,
        )
        .unwrap()
}

fn conserved(report: &RadiationReport) {
    let total: f64 = report
        .surface_heat_power
        .iter()
        .chain(report.gas_heat_power.iter())
        .map(|q| q.value())
        .sum();
    assert!(total.abs() < 1e-7, "net power {total:e}");
    assert_eq!(total, report.energy_residual_w);
}

#[test]
fn transparent_sheet_transfers_between_baths_without_heating() {
    let model = geometry([2.0; 2], false);
    let source = sheet([0.0; 2], 1.0);
    let report = run(&model, source.clone(), [900.0, 500.0, 500.0, 300.0]);
    let heat = 2.0 * (blackbody(900.0) - blackbody(300.0));
    close(report.surface_heat_power[0].value(), -heat);
    close(report.surface_heat_power[3].value(), heat);
    assert_eq!(report.surface_heat_power[1].value(), 0.0);
    assert_eq!(report.surface_heat_power[2].value(), 0.0);
    close(report.surface_radiosity_w_m2[1], blackbody(300.0));
    close(report.surface_radiosity_w_m2[2], blackbody(900.0));
    assert_eq!(report.sheets, vec![source]);
    assert_eq!(report.surface_emissivities.len(), 2);
    conserved(&report);
}

#[test]
fn asymmetric_sheet_matches_face_absorption_and_outgoing_energy() {
    let model = geometry([1.7; 2], false);
    let a = [0.2, 0.4];
    let t = 0.3;
    let source = sheet(a, t);
    let r = source.optics.reflectance();
    let report = run(&model, source, [900.0, 500.0, 500.0, 300.0]);
    let [hot, wall, cold] = [900.0, 500.0, 300.0].map(blackbody);
    let j = [
        a[0] * wall + r[0] * hot + t * cold,
        a[1] * wall + r[1] * cold + t * hot,
    ];
    close(report.surface_radiosity_w_m2[1], j[0]);
    close(report.surface_radiosity_w_m2[2], j[1]);
    close(
        report.surface_heat_power[1].value(),
        1.7 * a[0] * (hot - wall),
    );
    close(
        report.surface_heat_power[2].value(),
        1.7 * a[1] * (cold - wall),
    );
    close(report.surface_heat_power[0].value(), 1.7 * (j[0] - hot));
    close(report.surface_heat_power[3].value(), 1.7 * (j[1] - cold));
    conserved(&report);
}

#[test]
fn asymmetric_sheet_at_equilibrium_has_no_net_heat() {
    let model = geometry([1.0; 2], false);
    let report = run(&model, sheet([0.1, 0.65], 0.25), [700.0; 4]);
    for power in report
        .surface_heat_power
        .iter()
        .chain(report.gas_heat_power.iter())
    {
        assert_eq!(power.value(), 0.0);
    }
    for &j in &report.surface_radiosity_w_m2 {
        close(j, blackbody(700.0));
    }
    conserved(&report);
}

#[test]
fn zero_transmission_recovers_original_opaque_solver() {
    let model = geometry([1.0; 2], false);
    let temperatures = [900.0, 500.0, 500.0, 300.0].map(Temperature::new);
    let cards = [0.65, 0.2, 0.6, 0.8].map(emissivity);
    let old = model
        .evaluate(
            &temperatures,
            &cards,
            [Temperature::new(300.0); 2],
            &absorption([0.0; 2]),
            RadiationSolve::default(),
            &mut model.workspace().unwrap(),
            || false,
        )
        .unwrap();
    let new = model
        .evaluate_with_sheets(
            &temperatures,
            &[cards[0].clone(), cards[3].clone()],
            &[sheet([0.2, 0.6], 0.0)],
            [Temperature::new(300.0); 2],
            &absorption([0.0; 2]),
            RadiationSolve::default(),
            &mut model.workspace().unwrap(),
            || false,
        )
        .unwrap();
    for (a, b) in new.surface_heat_power.iter().zip(&old.surface_heat_power) {
        close(a.value(), b.value());
    }
    for (&a, &b) in new
        .surface_radiosity_w_m2
        .iter()
        .zip(&old.surface_radiosity_w_m2)
    {
        close(a, b);
    }
    assert!(old.sheets.is_empty());
    conserved(&new);
}

#[test]
fn transparent_sheet_carries_radiation_into_the_opposite_gas() {
    let model = geometry([1.3; 2], true);
    let report = model
        .evaluate_with_sheets(
            &[900.0, 500.0, 500.0, 300.0].map(Temperature::new),
            &[emissivity(1.0), emissivity(1.0)],
            &[sheet([0.0; 2], 1.0)],
            [600.0, 350.0].map(Temperature::new),
            &absorption([0.7, 0.9]),
            RadiationSolve::default(),
            &mut model.workspace().unwrap(),
            || false,
        )
        .unwrap();
    let [hot, cold, lower, upper] = [900.0, 300.0, 600.0, 350.0].map(blackbody);
    let [u, v] = [(-0.7_f64).exp(), (-1.8_f64).exp()];
    let j_left = v * cold + (1.0 - v) * upper;
    let j_right = u * hot + (1.0 - u) * lower;
    close(report.surface_radiosity_w_m2[1], j_left);
    close(report.surface_radiosity_w_m2[2], j_right);
    close(
        report.gas_heat_power[0].value(),
        1.3 * (1.0 - u) * (hot + j_left - 2.0 * lower),
    );
    close(
        report.gas_heat_power[1].value(),
        1.3 * (1.0 - v) * (j_right + cold - 2.0 * upper),
    );
    close(
        report.surface_heat_power[0].value(),
        1.3 * (u * j_left + (1.0 - u) * lower - hot),
    );
    close(
        report.surface_heat_power[3].value(),
        1.3 * (v * j_right + (1.0 - v) * upper - cold),
    );
    assert_eq!(report.surface_heat_power[1].value(), 0.0);
    assert_eq!(report.surface_heat_power[2].value(), 0.0);
    conserved(&report);
}

#[test]
fn incompatible_pair_geometry_state_and_source_axes_refuse() {
    for (areas, temperatures, pairs, cards) in [
        ([1.0, 2.0], [500.0; 4], vec![sheet([0.2; 2], 0.5)], 2),
        (
            [1.0; 2],
            [500.0, 500.0, 501.0, 500.0],
            vec![sheet([0.2; 2], 0.5)],
            2,
        ),
        (
            [1.0; 2],
            [500.0, 190.0, 190.0, 500.0],
            vec![sheet([0.2; 2], 0.5)],
            2,
        ),
        ([1.0; 2], [500.0; 4], vec![sheet([0.2; 2], 0.5); 2], 0),
        ([1.0; 2], [500.0; 4], vec![sheet([0.2; 2], 0.5)], 4),
    ] {
        let model = geometry(areas, false);
        assert!(
            model
                .evaluate_with_sheets(
                    &temperatures.map(Temperature::new),
                    &vec![emissivity(1.0); cards],
                    &pairs,
                    [Temperature::new(300.0); 2],
                    &absorption([0.0; 2]),
                    RadiationSolve::default(),
                    &mut model.workspace().unwrap(),
                    || false,
                )
                .is_err()
        );
    }
}

#[test]
fn lossless_transmitting_or_reflecting_cycle_refuses() {
    let model = ParticipatingEnclosure::new(
        ViewFactorMatrix::infinite_parallel_plates(1.0).unwrap(),
        vec![OpticalPath {
            surfaces: [0, 1],
            exchange_area: Area::new(1.0),
            segments: Vec::new(),
        }],
        "closed sheet-only cycle",
        1e-12,
        || false,
    )
    .unwrap();
    for t in [0.0, 0.6, 1.0] {
        let mut pair = sheet([0.0; 2], t);
        pair.surfaces = [0, 1];
        let result = model.evaluate_with_sheets(
            &[Temperature::new(500.0); 2],
            &[],
            &[pair],
            [Temperature::new(300.0); 2],
            &absorption([0.0; 2]),
            RadiationSolve::default(),
            &mut model.workspace().unwrap(),
            || false,
        );
        assert!(matches!(
            result,
            Err(RadiationError::Invalid(
                "closed lossless sheet/surface component has indeterminate radiosity"
            ))
        ));
    }
}

#[test]
fn transmission_uses_the_opposite_faces_implicit_self_return() {
    let s = 0.25;
    let factors = vec![
        vec![s, 1.0 - s, 0.0, 0.0],
        vec![1.0 - s, s, 0.0, 0.0],
        vec![0.0, 0.0, s, 1.0 - s],
        vec![0.0, 0.0, 1.0 - s, s],
    ];
    let matrix = ViewFactorMatrix::admit(
        vec![1.0; 4],
        factors,
        ViewFactorEvidence::Analytic {
            geometry: "paired diffuse surfaces with explicit self view".into(),
        },
        ViewFactorTolerance::default(),
    )
    .unwrap();
    let mut paths: Vec<_> = (0..4)
        .map(|i| OpticalPath {
            surfaces: [i, i],
            exchange_area: Area::new(0.5 * s),
            segments: Vec::new(),
        })
        .collect();
    for surfaces in [[0, 1], [2, 3]] {
        paths.push(OpticalPath {
            surfaces,
            exchange_area: Area::new(1.0 - s),
            segments: Vec::new(),
        });
    }
    let model =
        ParticipatingEnclosure::new(matrix, paths, "declared self-return chords", 1e-12, || {
            false
        })
        .unwrap();
    let report = run(&model, sheet([0.0; 2], 1.0), [900.0, 500.0, 500.0, 300.0]);
    let hot = blackbody(900.0);
    let cold = blackbody(300.0);
    close(
        report.surface_radiosity_w_m2[1],
        (s * hot + cold) / (1.0 + s),
    );
    close(
        report.surface_radiosity_w_m2[2],
        (hot + s * cold) / (1.0 + s),
    );
    let heat = (1.0 - s) / (1.0 + s) * (hot - cold);
    close(report.surface_heat_power[0].value(), -heat);
    close(report.surface_heat_power[3].value(), heat);
    assert_eq!(report.surface_heat_power[1].value(), 0.0);
    assert_eq!(report.surface_heat_power[2].value(), 0.0);
    conserved(&report);
}
fn band_source(value: f64) -> BandSurfaceEmissivity {
    BandSurfaceEmissivity::new(
        value,
        [200.0, 1500.0].map(Temperature::new),
        "analytic constant spectral absorptance",
    )
    .unwrap()
}

fn band(lo: f64, hi: f64, k: [f64; 2]) -> SpectralBand {
    SpectralBand {
        wavelengths: [lo, hi].map(Length::new),
        absorption: absorption(k),
    }
}

#[test]
fn selective_spectral_sheet_uses_one_solve_per_band_and_retains_sources() {
    let model = geometry([1.7; 2], false);
    let pairs = vec![vec![sheet([0.0; 2], 1.0)], vec![sheet([0.4, 0.2], 0.1)]];
    let optics = SpectralOptics::new_with_sheet_bands(
        vec![band(1e-7, 5e-6, [0.0; 2]), band(5e-6, 1e-3, [0.0; 2])],
        [None; 2],
        vec![
            vec![1.0, 0.0, 0.0, 1.0]
                .into_iter()
                .map(band_source)
                .collect(),
            vec![1.0, 0.4, 0.2, 1.0]
                .into_iter()
                .map(band_source)
                .collect(),
        ],
        pairs.clone(),
    )
    .unwrap();
    let report = model
        .evaluate_spectral(
            &[900.0, 500.0, 500.0, 300.0].map(Temperature::new),
            &[],
            [Temperature::new(300.0); 2],
            &optics,
            SpectralSolve {
                maximum_evaluations: 0,
                ..SpectralSolve::default()
            },
            &mut model.spectral_workspace().unwrap(),
            || false,
        )
        .unwrap();
    let mut heat = [0.0; 4];
    let mut radiosity = [0.0; 4];
    for (i, interval) in [[1e-7, 5e-6], [5e-6, 1e-3]].into_iter().enumerate() {
        let window = BlackbodyBand::new(interval[0], interval[1]).unwrap();
        let [hot, wall, cold] =
            [900.0, 500.0, 300.0].map(|t| window.evaluate(t).unwrap().power_w_m2);
        let source = &pairs[i][0].optics;
        let a = source.absorptance();
        let r = source.reflectance();
        let t = source.transmittance();
        let j = [
            a[0] * wall + r[0] * hot + t * cold,
            a[1] * wall + r[1] * cold + t * hot,
        ];
        for (total, value) in radiosity.iter_mut().zip([hot, j[0], j[1], cold]) {
            *total += value;
        }
        for (total, value) in heat.iter_mut().zip([
            1.7 * (j[0] - hot),
            1.7 * a[0] * (hot - wall),
            1.7 * a[1] * (cold - wall),
            1.7 * (j[1] - cold),
        ]) {
            *total += value;
        }
    }
    for (actual, expected) in report.surface_heat_power.iter().zip(heat) {
        close(actual.value(), expected);
    }
    for (&actual, expected) in report.surface_radiosity_w_m2.iter().zip(radiosity) {
        close(actual, expected);
    }
    assert_eq!(report.factorizations, 2);
    assert_eq!(report.wavelength_evaluations, 0);
    assert_eq!(report.accepted_panels, 0);
    assert_eq!(report.optics.sheet_bands().unwrap(), pairs.as_slice());
    assert!(report.surface_emissivities.is_empty());
    assert!(report.energy_residual_w.abs() < 1e-7);
}

#[test]
fn spectral_sheet_rejects_conflicting_absorptance_and_changing_pair_axis() {
    let bands = vec![band(1e-7, 5e-6, [0.0; 2]), band(5e-6, 1e-3, [0.0; 2])];
    let sources = vec![
        vec![1.0, 0.2, 0.4, 1.0]
            .into_iter()
            .map(band_source)
            .collect();
        2
    ];
    let pairs = vec![vec![sheet([0.2, 0.4], 0.1)]; 2];
    let mut conflicting = pairs.clone();
    conflicting[1][0].optics = sheet([0.3, 0.4], 0.1).optics;
    assert!(
        SpectralOptics::new_with_sheet_bands(
            bands.clone(),
            [None; 2],
            sources.clone(),
            conflicting,
        )
        .is_err()
    );
    let mut changed = pairs.clone();
    changed[1][0].surfaces = [2, 1];
    assert!(matches!(
        SpectralOptics::new_with_sheet_bands(bands, [None; 2], sources, changed,),
        Err(RadiationError::Invalid(
            "spectral sheet pairing and side order must be identical in every band"
        ))
    ));
}

#[test]
fn spectral_sheet_checks_actual_isothermal_source_temperature() {
    let model = geometry([1.0; 2], false);
    let mut pair = sheet([0.2, 0.4], 0.1);
    pair.optics =
        ReciprocalSheetOptics::new([0.2, 0.4], 0.1, [600.0, 800.0], "limited-temperature sheet")
            .unwrap();
    let optics = SpectralOptics::new_with_sheet_bands(
        vec![band(1e-7, 1e-3, [0.0; 2])],
        [None; 2],
        vec![
            vec![1.0, 0.2, 0.4, 1.0]
                .into_iter()
                .map(band_source)
                .collect(),
        ],
        vec![vec![pair]],
    )
    .unwrap();
    assert!(matches!(
        model.evaluate_spectral(
            &[900.0, 500.0, 500.0, 300.0].map(Temperature::new),
            &[],
            [Temperature::new(300.0); 2],
            &optics,
            SpectralSolve::default(),
            &mut model.spectral_workspace().unwrap(),
            || false,
        ),
        Err(RadiationError::Invalid(
            "sheet temperature is outside original optical support"
        ))
    ));
}

#[test]
fn soot_quadrature_across_a_transparent_sheet_matches_the_original_two_segment_path() {
    let model = geometry([1.3; 2], true);
    let soot = RayleighSoot::new(
        1.75,
        1.03,
        Density::new(1800.0),
        SootOpticalDomain {
            spectral_temperature: [200.0, 1500.0].map(Temperature::new),
            maximum_volume_fraction: 0.001,
        },
        "explicit analytic constant-index soot",
    )
    .unwrap();
    let population = [
        Some(SootPopulation {
            model: &soot,
            volume_fraction: 1e-6,
        }),
        None,
    ];
    let bands = vec![band(1e-7, 1e-3, [0.1, 0.8])];
    let optics = SpectralOptics::new_with_sheet_bands(
        bands.clone(),
        population,
        vec![
            vec![1.0, 0.0, 0.0, 1.0]
                .into_iter()
                .map(band_source)
                .collect(),
        ],
        vec![vec![sheet([0.0; 2], 1.0)]],
    )
    .unwrap();
    let unsplit = ParticipatingEnclosure::new(
        ViewFactorMatrix::infinite_parallel_plates(1.3).unwrap(),
        vec![OpticalPath {
            surfaces: [0, 1],
            exchange_area: Area::new(1.3),
            segments: vec![
                OpticalSegment {
                    layer: Layer::Lower,
                    length: Length::new(1.0),
                },
                OpticalSegment {
                    layer: Layer::Upper,
                    length: Length::new(2.0),
                },
            ],
        }],
        "the same two gas segments with no intermediate boundary",
        1e-12,
        || false,
    )
    .unwrap();
    let original =
        SpectralOptics::new_with_surface_bands(bands, population, vec![vec![band_source(1.0); 2]])
            .unwrap();
    let controls = SpectralSolve {
        heat_tolerance_w: 1e-7,
        radiosity_tolerance_w_m2: 1e-7,
        ..SpectralSolve::default()
    };
    let gas = [600.0, 350.0].map(Temperature::new);
    let report = model
        .evaluate_spectral(
            &[900.0, 500.0, 500.0, 300.0].map(Temperature::new),
            &[],
            gas,
            &optics,
            controls,
            &mut model.spectral_workspace().unwrap(),
            || false,
        )
        .unwrap();
    let old = unsplit
        .evaluate_spectral(
            &[900.0, 300.0].map(Temperature::new),
            &[],
            gas,
            &original,
            controls,
            &mut unsplit.spectral_workspace().unwrap(),
            || false,
        )
        .unwrap();
    for (actual, expected) in [report.surface_heat_power[0], report.surface_heat_power[3]]
        .iter()
        .zip(&old.surface_heat_power)
        .chain(report.gas_heat_power.iter().zip(&old.gas_heat_power))
    {
        assert!((actual.value() - expected.value()).abs() < 1e-6);
    }
    assert_eq!(report.surface_heat_power[1].value(), 0.0);
    assert_eq!(report.surface_heat_power[2].value(), 0.0);
    assert!(report.wavelength_evaluations > 0);
    assert_eq!(report.factorizations, report.wavelength_evaluations);
    assert!(report.energy_residual_w.abs() < 1e-7);
}
