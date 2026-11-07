use super::super::soot::{RayleighSoot, SootOpticalDomain};
use super::super::spectral::{
    BandSurfaceEmissivity, SootPopulation, SpectralBand, SpectralOptics, SpectralSolve,
};
use super::super::{AbsorptionCoefficient, OpticalPath, OpticalSegment, STEFAN_BOLTZMANN_W_M2_K4};
use super::*;
use fs_conduction::radiation::BlackbodyBand;
use fs_conduction::radiation::{
    EMISSIVITY_DIMS, SURFACE_EMISSIVITY_PROPERTY, ViewFactorEvidence, ViewFactorMatrix,
    ViewFactorTolerance,
    slab::{AbsorbingLayer, SlabQuadrature},
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, SelectionPolicy, UncertaintyModel,
};
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
                source: "analytic gray reservoir".into(),
                license: "internal-test-use".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "fixture".into(),
            phase: "solid".into(),
            process: "opaque reservoir".into(),
            revision: 0,
        },
        claims,
        Vec::new(),
    )
    .unwrap();
    SurfaceEmissivity::from_card("bath", &card, 300.0, SelectionPolicy::SingleClaimOnly).unwrap()
}

fn slab(depths: &[f64]) -> SlabPair {
    SlabPair {
        surfaces: [1, 2],
        optics: AbsorbingSlab::new(
            depths
                .iter()
                .map(|&depth| {
                    AbsorbingLayer::new(
                        1.0,
                        depth,
                        [200.0, 1500.0],
                        "explicit constant absorption coefficient",
                    )
                    .unwrap()
                })
                .collect(),
            SlabQuadrature::default(),
            || false,
        )
        .unwrap(),
    }
}

fn absorption(k: [f64; 2]) -> [LayerAbsorption; 2] {
    k.map(|k| {
        LayerAbsorption::new(
            AbsorptionCoefficient::new(k),
            "explicit uniform gas coefficient",
        )
        .unwrap()
    })
}

fn geometry(areas: [f64; 2], gas: bool) -> ParticipatingEnclosure {
    let matrix = ViewFactorMatrix::admit(
        vec![areas[0], areas[0], areas[1], areas[1]],
        vec![
            vec![0.0, 1.0, 0.0, 0.0],
            vec![1.0, 0.0, 0.0, 0.0],
            vec![0.0, 0.0, 0.0, 1.0],
            vec![0.0, 0.0, 1.0, 0.0],
        ],
        ViewFactorEvidence::Analytic {
            geometry: "two diffuse slab ports facing reservoirs".into(),
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
                    length: Length::new(1.0),
                }]
            } else {
                Vec::new()
            },
        })
        .collect();
    ParticipatingEnclosure::new(matrix, paths, "authored diffuse-port paths", 1e-12, || {
        false
    })
    .unwrap()
}

fn blackbody(t: f64) -> f64 {
    STEFAN_BOLTZMANN_W_M2_K4 * t.powi(4)
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-8 * expected.abs().max(1.0),
        "{actual:e} != {expected:e}"
    );
}
fn conserved(report: &RadiationReport) {
    let total: f64 = report
        .surface_heat_power
        .iter()
        .chain(report.gas_heat_power.iter())
        .chain(report.layer_heat_power.iter().flatten())
        .map(|q| q.value())
        .sum();
    assert!(total.abs() < 1e-7, "net heat {total:e}");
    assert_eq!(report.energy_residual_w, total);
}
fn run(model: &ParticipatingEnclosure, pair: SlabPair, temperatures: &[f64]) -> RadiationReport {
    model
        .evaluate_with_slabs(
            &[
                900.0,
                temperatures[0],
                temperatures[temperatures.len() - 1],
                250.0,
            ]
            .map(Temperature::new),
            &[1.0; 2].map(emissivity),
            &[pair],
            &[temperatures.iter().copied().map(Temperature::new).collect()],
            [Temperature::new(300.0); 2],
            &absorption([0.0; 2]),
            RadiationSolve::default(),
            &mut model.workspace().unwrap(),
            || false,
        )
        .unwrap()
}

#[test]
fn nonisothermal_two_layer_heat_matches_independent_exponential_integral_moments() {
    // Exact diffuse transmission: 2E3(0.5), 2E3(1), independently tabulated.
    let half = 0.443_208_728_550_356_9;
    let full = 0.219_383_934_395_520_27;
    let near = 1.0 - half;
    let far = half - full;
    let internal = 1.0 - 2.0 * half + full;
    for with_gas in [false, true] {
        let model = geometry([1.3; 2], with_gas);
        let pair = slab(&[0.5, 0.5]);
        let report = model
            .evaluate_with_slabs(
                &[900.0, 700.0, 400.0, 250.0].map(Temperature::new),
                &[1.0; 2].map(emissivity),
                &[pair.clone()],
                &[vec![Temperature::new(700.0), Temperature::new(400.0)]],
                [600.0, 350.0].map(Temperature::new),
                &absorption(if with_gas { [0.7, 0.9] } else { [0.0; 2] }),
                RadiationSolve::default(),
                &mut model.workspace().unwrap(),
                || false,
            )
            .unwrap();
        let [hot, first, second, cold, lower, upper] =
            [900.0, 700.0, 400.0, 250.0, 600.0, 350.0].map(blackbody);
        let u = if with_gas { (-0.7_f64).exp() } else { 1.0 };
        let v = if with_gas { (-0.9_f64).exp() } else { 1.0 };
        let h_left = u * hot + (1.0 - u) * lower;
        let h_right = v * cold + (1.0 - v) * upper;
        let j_left = full * h_right + near * first + far * second;
        let j_right = full * h_left + far * first + near * second;
        close(report.surface_radiosity_w_m2[1], j_left);
        close(report.surface_radiosity_w_m2[2], j_right);
        close(
            report.layer_heat_power[0][0].value(),
            1.3 * (near * (h_left - first) + far * (h_right - first) + internal * (second - first)),
        );
        close(
            report.layer_heat_power[0][1].value(),
            1.3 * (far * (h_left - second)
                + near * (h_right - second)
                + internal * (first - second)),
        );
        close(
            report.gas_heat_power[0].value(),
            1.3 * (1.0 - u) * (hot + j_left - 2.0 * lower),
        );
        close(
            report.gas_heat_power[1].value(),
            1.3 * (1.0 - v) * (cold + j_right - 2.0 * upper),
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
        assert_eq!(report.slabs, vec![pair]);
        conserved(&report);
    }
}

#[test]
fn transparent_layers_pass_radiation_and_receive_exactly_zero_heat() {
    let model = geometry([2.0; 2], false);
    let report = run(&model, slab(&[0.0; 3]), &[300.0, 800.0, 500.0]);
    for q in &report.layer_heat_power[0] {
        assert_eq!(q.value(), 0.0);
    }
    assert_eq!(report.surface_heat_power[1].value(), 0.0);
    assert_eq!(report.surface_heat_power[2].value(), 0.0);
    close(
        report.surface_heat_power[0].value(),
        2.0 * (blackbody(250.0) - blackbody(900.0)),
    );
    close(
        report.surface_heat_power[3].value(),
        -report.surface_heat_power[0].value(),
    );
    conserved(&report);
}

#[test]
fn g1_refractive_slab_reflection_reaches_enclosure_and_original_layer_heat() {
    let model = geometry([1.3; 2], false);
    for indices in [[1.5; 2], [1.4, 2.1]] {
        let mut pair = slab(&[0.2, 0.6]);
        pair.optics = AbsorbingSlab::new(
            pair.optics
                .layers()
                .iter()
                .cloned()
                .zip(indices)
                .map(|(layer, index)| {
                    layer
                        .with_refractive_index(index, "authored layer index")
                        .unwrap()
                })
                .collect(),
            SlabQuadrature::default(),
            || false,
        )
        .unwrap();
        let report = run(&model, pair.clone(), &[700.0, 400.0]);
        let h = [blackbody(900.0), blackbody(250.0)];
        let b = [blackbody(700.0), blackbody(400.0)];
        let emission = pair.optics.emission_to_faces(&b).unwrap();
        let mut heat = [0.0; 2];
        pair.optics.layer_net_fluxes(h, &b, &mut heat).unwrap();
        for side in 0..2 {
            close(
                report.surface_radiosity_w_m2[1 + side],
                emission[side]
                    + pair.optics.reflectances()[side] * h[side]
                    + pair.optics.transmittance() * h[1 - side],
            );
            close(report.layer_heat_power[0][side].value(), 1.3 * heat[side]);
            assert_eq!(report.surface_heat_power[1 + side].value(), 0.0);
        }
        assert!(pair.optics.reflectance() > 0.08);
        if indices[0] != indices[1] {
            assert!((pair.optics.reflectances()[0] - pair.optics.reflectances()[1]).abs() > 0.01);
        }
        conserved(&report);

        let equilibrium = model
            .evaluate_with_slabs(
                &[Temperature::new(700.0); 4],
                &[0.6; 2].map(emissivity),
                &[pair],
                &[vec![Temperature::new(700.0); 2]],
                [Temperature::new(700.0); 2],
                &absorption([0.0; 2]),
                RadiationSolve::default(),
                &mut model.workspace().unwrap(),
                || false,
            )
            .unwrap();
        for q in equilibrium
            .surface_heat_power
            .iter()
            .chain(equilibrium.layer_heat_power.iter().flatten())
        {
            assert_eq!(q.value(), 0.0);
        }
    }
}

#[test]
fn opaque_layers_emit_locally_and_exchange_heat_internally() {
    let model = geometry([1.0; 2], false);
    let report = run(&model, slab(&[1000.0; 2]), &[700.0, 400.0]);
    let [hot, first, second, cold] = [900.0, 700.0, 400.0, 250.0].map(blackbody);
    close(report.surface_radiosity_w_m2[1], first);
    close(report.surface_radiosity_w_m2[2], second);
    close(
        report.layer_heat_power[0][0].value(),
        hot + second - 2.0 * first,
    );
    close(
        report.layer_heat_power[0][1].value(),
        first + cold - 2.0 * second,
    );
    conserved(&report);
}

#[test]
fn isothermal_subdivision_preserves_full_stack_angular_memory() {
    let model = geometry([1.0; 2], false);
    let whole = run(&model, slab(&[1.0]), &[600.0]);
    let split = run(&model, slab(&[0.5, 0.5]), &[600.0, 600.0]);
    for (a, b) in split
        .surface_heat_power
        .iter()
        .zip(&whole.surface_heat_power)
    {
        close(a.value(), b.value());
    }
    for (&a, &b) in split
        .surface_radiosity_w_m2
        .iter()
        .zip(&whole.surface_radiosity_w_m2)
    {
        close(a, b);
    }
    close(
        split.layer_heat_power[0].iter().map(|q| q.value()).sum(),
        whole.layer_heat_power[0][0].value(),
    );
    conserved(&split);
}

#[test]
fn full_thermal_equilibrium_has_zero_layer_gas_and_surface_heat() {
    let model = geometry([1.0; 2], true);
    let report = model
        .evaluate_with_slabs(
            &[Temperature::new(700.0); 4],
            &[0.6; 2].map(emissivity),
            &[slab(&[0.2, 3.0])],
            &[vec![Temperature::new(700.0); 2]],
            [Temperature::new(700.0); 2],
            &absorption([0.3, 0.8]),
            RadiationSolve::default(),
            &mut model.workspace().unwrap(),
            || false,
        )
        .unwrap();
    for q in report
        .surface_heat_power
        .iter()
        .chain(report.gas_heat_power.iter())
        .chain(report.layer_heat_power.iter().flatten())
    {
        assert_eq!(q.value(), 0.0);
    }
    conserved(&report);
}

#[test]
fn slab_source_geometry_temperature_and_lossless_axes_are_enforced() {
    for (areas, face, temperatures) in [
        ([1.0, 2.0], [900.0, 700.0, 400.0, 250.0], vec![700.0, 400.0]),
        ([1.0; 2], [900.0, 701.0, 400.0, 250.0], vec![700.0, 400.0]),
        ([1.0; 2], [900.0, 700.0, 400.0, 250.0], vec![700.0]),
        ([1.0; 2], [900.0, 700.0, 190.0, 250.0], vec![700.0, 190.0]),
    ] {
        let model = geometry(areas, false);
        assert!(
            model
                .evaluate_with_slabs(
                    &face.map(Temperature::new),
                    &[1.0; 2].map(emissivity),
                    &[slab(&[0.5; 2])],
                    &[temperatures.into_iter().map(Temperature::new).collect()],
                    [Temperature::new(300.0); 2],
                    &absorption([0.0; 2]),
                    RadiationSolve::default(),
                    &mut model.workspace().unwrap(),
                    || false,
                )
                .is_err()
        );
    }
    let closed = ParticipatingEnclosure::new(
        ViewFactorMatrix::infinite_parallel_plates(1.0).unwrap(),
        vec![OpticalPath {
            surfaces: [0, 1],
            exchange_area: Area::new(1.0),
            segments: Vec::new(),
        }],
        "lossless closed ports",
        1e-12,
        || false,
    )
    .unwrap();
    let mut pair = slab(&[0.0]);
    pair.surfaces = [0, 1];
    assert!(
        closed
            .evaluate_with_slabs(
                &[Temperature::new(500.0); 2],
                &[],
                &[pair],
                &[vec![Temperature::new(500.0)]],
                [Temperature::new(300.0); 2],
                &absorption([0.0; 2]),
                RadiationSolve::default(),
                &mut closed.workspace().unwrap(),
                || false,
            )
            .is_err()
    );
}

fn spectral_band(lo: f64, hi: f64, k: [f64; 2]) -> SpectralBand {
    SpectralBand {
        wavelengths: [lo, hi].map(Length::new),
        absorption: absorption(k),
    }
}
fn band_emissivity(value: f64) -> BandSurfaceEmissivity {
    BandSurfaceEmissivity::new(
        value,
        [200.0, 1500.0].map(Temperature::new),
        "explicit constant band emissivity",
    )
    .unwrap()
}

#[test]
fn spectral_layers_match_reflecting_reservoir_solution_with_one_factorization_per_band() {
    let model = geometry([1.3; 2], false);
    let pairs = vec![vec![slab(&[0.0; 2])], vec![slab(&[0.5; 2])]];
    let optics = SpectralOptics::new_with_slab_bands(
        vec![
            spectral_band(1e-7, 5e-6, [0.0; 2]),
            spectral_band(5e-6, 1e-3, [0.0; 2]),
        ],
        [None; 2],
        vec![vec![band_emissivity(0.6), band_emissivity(0.8)]; 2],
        pairs.clone(),
    )
    .unwrap();
    let layer_t = vec![vec![Temperature::new(700.0), Temperature::new(400.0)]];
    let report = model
        .evaluate_spectral_with_slabs(
            &[900.0, 700.0, 400.0, 250.0].map(Temperature::new),
            &[],
            &layer_t,
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
    let mut heat = [0.0; 4]; // hot bath, cold bath, first material layer, second layer
    let mut radiosity = [0.0; 4];
    for (band, window) in [[1e-7, 5e-6], [5e-6, 1e-3]].into_iter().enumerate() {
        let planck = BlackbodyBand::new(window[0], window[1]).unwrap();
        let [hot, first, second, cold] =
            [900.0, 700.0, 400.0, 250.0].map(|t| planck.evaluate(t).unwrap().power_w_m2);
        let (t, near, far, exchange) = if band == 0 {
            (1.0, 0.0, 0.0, 0.0)
        } else {
            let half = 0.443_208_728_550_356_9;
            let full = 0.219_383_934_395_520_27;
            (full, 1.0 - half, half - full, 1.0 - 2.0 * half + full)
        };
        let emitted = [near * first + far * second, far * first + near * second];
        // Solve two independently derived reflecting-reservoir equations.
        let a = 0.6 * hot + 0.4 * emitted[0];
        let b = 0.8 * cold + 0.2 * emitted[1];
        let denominator = 1.0 - 0.4 * 0.2 * t * t;
        let j_hot = (a + 0.4 * t * b) / denominator;
        let j_cold = (b + 0.2 * t * a) / denominator;
        let j_left = t * j_cold + emitted[0];
        let j_right = t * j_hot + emitted[1];
        for (sum, value) in radiosity.iter_mut().zip([j_hot, j_left, j_right, j_cold]) {
            *sum += value;
        }
        for (sum, value) in heat.iter_mut().zip([
            1.3 * (j_left - j_hot),
            1.3 * (j_right - j_cold),
            1.3 * (near * (j_hot - first) + far * (j_cold - first) + exchange * (second - first)),
            1.3 * (far * (j_hot - second) + near * (j_cold - second) + exchange * (first - second)),
        ]) {
            *sum += value;
        }
    }
    for (&actual, expected) in report.surface_radiosity_w_m2.iter().zip(radiosity) {
        close(actual, expected);
    }
    for (actual, expected) in [
        report.surface_heat_power[0],
        report.surface_heat_power[3],
        report.layer_heat_power[0][0],
        report.layer_heat_power[0][1],
    ]
    .into_iter()
    .zip(heat)
    {
        close(actual.value(), expected);
    }
    assert_eq!(report.surface_heat_power[1].value(), 0.0);
    assert_eq!(report.surface_heat_power[2].value(), 0.0);
    assert_eq!(report.factorizations, 2);
    assert_eq!(report.wavelength_evaluations, 0);
    assert_eq!(report.layer_temperatures, layer_t);
    assert_eq!(report.optics.slab_bands().unwrap(), pairs.as_slice());
    assert!(report.layer_coverage[0].iter().all(|c| c.participating));
    assert!(
        report
            .layer_heat_error_estimate_w
            .iter()
            .flatten()
            .all(|&e| e == 0.0)
    );
    assert!(report.energy_residual_w.abs() < 1e-7);
}

#[test]
fn slab_spectra_require_owned_layer_temperatures_and_unchanged_physical_geometry() {
    let model = geometry([1.0; 2], false);
    let make = |pairs| {
        SpectralOptics::new_with_slab_bands(
            vec![
                spectral_band(1e-7, 5e-6, [0.0; 2]),
                spectral_band(5e-6, 1e-3, [0.0; 2]),
            ],
            [None; 2],
            vec![vec![band_emissivity(1.0); 2]; 2],
            pairs,
        )
    };
    assert!(make(vec![vec![slab(&[0.5; 2])], vec![slab(&[1.0])]]).is_err());
    let mut geometry_changed = slab(&[0.5; 2]);
    geometry_changed.optics = AbsorbingSlab::new(
        (0..2)
            .map(|_| {
                AbsorbingLayer::new(
                    2.0,
                    0.25,
                    [200.0, 1500.0],
                    "same optical depth with different material thickness",
                )
                .unwrap()
            })
            .collect(),
        SlabQuadrature::default(),
        || false,
    )
    .unwrap();
    assert!(make(vec![vec![slab(&[0.5; 2])], vec![geometry_changed]]).is_err());
    let optics = make(vec![vec![slab(&[0.5; 2])]; 2]).unwrap();
    assert!(matches!(
        model.evaluate_spectral(
            &[900.0, 700.0, 400.0, 250.0].map(Temperature::new),
            &[],
            [Temperature::new(300.0); 2],
            &optics,
            SpectralSolve::default(),
            &mut model.spectral_workspace().unwrap(),
            || false,
        ),
        Err(RadiationError::Invalid(
            "slab spectra require explicit layer temperatures through evaluate_spectral_with_slabs"
        ))
    ));
}

#[test]
fn hottest_interior_layer_controls_omitted_emission_even_when_both_ports_are_cold() {
    let model = geometry([1.0; 2], false);
    let optics = SpectralOptics::new_with_slab_bands(
        vec![spectral_band(5e-6, 1e-3, [0.0; 2])],
        [None; 2],
        vec![vec![band_emissivity(1.0); 2]],
        vec![vec![slab(&[0.5; 3])]],
    )
    .unwrap();
    let controls = SpectralSolve {
        maximum_omitted_blackbody_fraction: 0.05,
        ..SpectralSolve::default()
    };
    let mut work = model.spectral_workspace().unwrap();
    model
        .evaluate_spectral_with_slabs(
            &[Temperature::new(300.0); 4],
            &[],
            &[vec![Temperature::new(300.0); 3]],
            [Temperature::new(300.0); 2],
            &optics,
            controls,
            &mut work,
            || false,
        )
        .unwrap();
    assert!(matches!(model.evaluate_spectral_with_slabs(
        &[Temperature::new(300.0);4],&[],
        &[vec![Temperature::new(300.0),Temperature::new(1300.0),Temperature::new(300.0)]],
        [Temperature::new(300.0);2],&optics,controls,&mut work,||false,
    ),Err(RadiationError::SpectralTruncation{maximum_fraction,..}) if maximum_fraction>0.5));
}

#[test]
fn adaptive_soot_transport_resolves_each_nonisothermal_layer_against_independent_simpson_reference()
{
    let model = geometry([1.3; 2], true);
    let soot = RayleighSoot::new(
        1.75,
        1.03,
        Density::new(1800.0),
        SootOpticalDomain {
            spectral_temperature: [200.0, 1500.0].map(Temperature::new),
            maximum_volume_fraction: 0.001,
        },
        "explicit spectral absorption fixture",
    )
    .unwrap();
    let optics = SpectralOptics::new_with_slab_bands(
        vec![spectral_band(1e-7, 1e-3, [0.2, 0.7])],
        [
            Some(SootPopulation {
                model: &soot,
                volume_fraction: 1e-6,
            }),
            None,
        ],
        vec![vec![band_emissivity(1.0); 2]],
        vec![vec![slab(&[0.5; 2])]],
    )
    .unwrap();
    let report = model
        .evaluate_spectral_with_slabs(
            &[900.0, 700.0, 400.0, 250.0].map(Temperature::new),
            &[],
            &[vec![Temperature::new(700.0), Temperature::new(400.0)]],
            [600.0, 350.0].map(Temperature::new),
            &optics,
            SpectralSolve {
                heat_tolerance_w: 1e-6,
                radiosity_tolerance_w_m2: 1e-6,
                ..SpectralSolve::default()
            },
            &mut model.spectral_workspace().unwrap(),
            || false,
        )
        .unwrap();
    let panels = 8192;
    let left = 1e-7_f64.ln();
    let width = (1e-3_f64.ln() - left) / panels as f64;
    let mut reference = [0.0; 6];
    let t = 0.219_383_934_395_520_27;
    let half = 0.443_208_728_550_356_9;
    let near = 1.0 - half;
    let far = half - t;
    let exchange = 1.0 - 2.0 * half + t;
    for point in 0..=panels {
        let wavelength = (left + point as f64 * width).exp();
        let [hot, first, second, cold, lower, upper] = [900.0, 700.0, 400.0, 250.0, 600.0, 350.0]
            .map(|temperature| {
                BlackbodyBand::spectral_power_per_log_wavelength(wavelength, temperature).unwrap()
            });
        let u = (-(0.2
            + soot
                .spectral_coefficient(1e-6, Length::new(wavelength))
                .unwrap()
                .value()))
        .exp();
        let v = (-0.7_f64).exp();
        let h_left = u * hot + (1.0 - u) * lower;
        let h_right = v * cold + (1.0 - v) * upper;
        let j_left = t * h_right + near * first + far * second;
        let j_right = t * h_left + far * first + near * second;
        let weight = if point == 0 || point == panels {
            1.0
        } else if point % 2 == 0 {
            2.0
        } else {
            4.0
        };
        for (sum, value) in reference.iter_mut().zip([
            1.3 * (near * (h_left - first) + far * (h_right - first) + exchange * (second - first)),
            1.3 * (far * (h_left - second)
                + near * (h_right - second)
                + exchange * (first - second)),
            1.3 * (1.0 - u) * (hot + j_left - 2.0 * lower),
            1.3 * (1.0 - v) * (cold + j_right - 2.0 * upper),
            1.3 * (u * j_left + (1.0 - u) * lower - hot),
            1.3 * (v * j_right + (1.0 - v) * upper - cold),
        ]) {
            *sum += weight * value * width / 3.0;
        }
    }
    for (actual, expected) in [
        report.layer_heat_power[0][0],
        report.layer_heat_power[0][1],
        report.gas_heat_power[0],
        report.gas_heat_power[1],
        report.surface_heat_power[0],
        report.surface_heat_power[3],
    ]
    .into_iter()
    .zip(reference)
    {
        assert!((actual.value() - expected).abs() < 5e-5);
    }
    assert!(report.wavelength_evaluations > 0);
    assert_eq!(report.factorizations, report.wavelength_evaluations);
    assert!(
        report
            .layer_heat_error_estimate_w
            .iter()
            .flatten()
            .all(|&e| e <= 1e-6)
    );
    assert!(report.energy_residual_w.abs() < 1e-7);
}
