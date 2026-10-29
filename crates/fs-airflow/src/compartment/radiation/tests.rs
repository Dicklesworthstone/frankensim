use super::*;
use fs_conduction::radiation::{
    EMISSIVITY_DIMS, SURFACE_EMISSIVITY_PROPERTY, ViewFactorEvidence, ViewFactorTolerance,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, SelectionPolicy, UncertaintyModel,
};

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
                source: "analytic gray-boundary fixture".into(),
                license: "internal-test-use".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "fixture".into(),
            phase: "solid".into(),
            process: "gray surface".into(),
            revision: 0,
        },
        claims,
        Vec::new(),
    )
    .unwrap();
    SurfaceEmissivity::from_card("fixture", &card, 300.0, SelectionPolicy::SingleClaimOnly).unwrap()
}
fn segment(layer: Layer, length: f64) -> OpticalSegment {
    OpticalSegment {
        layer,
        length: Length::new(length),
    }
}
fn absorption(values: [f64; 2]) -> [LayerAbsorption; 2] {
    values.map(|value| {
        LayerAbsorption::new(
            AbsorptionCoefficient::new(value),
            "explicit analytic gray-medium fixture",
        )
        .unwrap()
    })
}
fn plates(area: f64, segments: Vec<OpticalSegment>) -> ParticipatingEnclosure {
    ParticipatingEnclosure::new(
        ViewFactorMatrix::infinite_parallel_plates(area).unwrap(),
        vec![OpticalPath {
            surfaces: [0, 1],
            exchange_area: Area::new(area),
            segments,
        }],
        "declared single-chord optical quadrature",
        1e-12,
        || false,
    )
    .unwrap()
}
fn run(
    model: &ParticipatingEnclosure,
    wall: [f64; 2],
    gas: [f64; 2],
    kappa: [f64; 2],
    e: [f64; 2],
) -> RadiationReport {
    model
        .evaluate(
            &wall.map(Temperature::new),
            &e.map(emissivity),
            gas.map(Temperature::new),
            &absorption(kappa),
            RadiationSolve::default(),
            &mut model.workspace().unwrap(),
            || false,
        )
        .unwrap()
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 5e-12 * expected.abs().max(1.0),
        "{actual:e} != {expected:e}"
    );
}
fn blackbody(t: f64) -> f64 {
    STEFAN_BOLTZMANN_W_M2_K4 * t.powi(4)
}
fn conserved(report: &RadiationReport) {
    let total: f64 = report
        .surface_heat_power
        .iter()
        .chain(report.gas_heat_power.iter())
        .map(|q| q.value())
        .sum();
    assert!(total.abs() < 1e-7, "net heat {total:e}");
    close(report.energy_residual_w, total);
}

#[test]
fn g1_transparent_medium_recovers_diffuse_parallel_plate_exchange() {
    let model = plates(
        2.0,
        vec![segment(Layer::Lower, 1.0), segment(Layer::Upper, 2.0)],
    );
    let report = run(
        &model,
        [900.0, 300.0],
        [1100.0, 220.0],
        [0.0, 0.0],
        [0.6, 0.8],
    );
    let heat = 2.0 * (blackbody(900.0) - blackbody(300.0)) / (1.0 / 0.6 + 1.0 / 0.8 - 1.0);
    close(report.surface_heat_power[0].value(), -heat);
    close(report.surface_heat_power[1].value(), heat);
    assert_eq!(report.gas_heat_power, [Power::new(0.0); 2]);
    conserved(&report);
}

#[test]
fn g1_black_boundaries_and_one_uniform_chord_match_formal_transfer() {
    let model = plates(1.3, vec![segment(Layer::Upper, 2.0)]);
    let report = run(
        &model,
        [1000.0, 300.0],
        [300.0, 700.0],
        [0.0, 0.7],
        [1.0, 1.0],
    );
    let transmission = (-1.4_f64).exp();
    let absorbed = -(-1.4_f64).exp_m1();
    let [b0, b1, bg] = [blackbody(1000.0), blackbody(300.0), blackbody(700.0)];
    close(
        report.surface_heat_power[0].value(),
        1.3 * (transmission * (b1 - b0) + absorbed * (bg - b0)),
    );
    close(
        report.surface_heat_power[1].value(),
        1.3 * (transmission * (b0 - b1) + absorbed * (bg - b1)),
    );
    close(
        report.gas_heat_power[1].value(),
        1.3 * absorbed * (b0 + b1 - 2.0 * bg),
    );
    assert_eq!(report.gas_heat_power[0].value(), 0.0);
    conserved(&report);
}

#[test]
fn g1_black_parallel_plate_slab_converges_to_hemispherical_exact_solution() {
    // For optical depth one, angular transmittance 2 E_3(1) = E_1(1).
    // Here E_n(x)=integral_1^infinity exp(-xt)t^(-n)dt; integration by parts
    // gives n E_(n+1)=exp(-x)-x E_n. This is an independent continuum value.
    let exact_transmission = 0.219_383_934_395_520_27;
    let mut errors = Vec::new();
    for count in [64, 256, 512] {
        let paths = (0..count)
            .map(|i| {
                let mu = (i as f64 + 0.5) / count as f64;
                OpticalPath {
                    surfaces: [0, 1],
                    exchange_area: Area::new(2.0 * mu / count as f64),
                    segments: vec![segment(Layer::Upper, 1.0 / mu)],
                }
            })
            .collect();
        let model = ParticipatingEnclosure::new(
            ViewFactorMatrix::infinite_parallel_plates(1.0).unwrap(),
            paths,
            "midpoint quadrature of exact parallel-plate diffuse angular measure 2 mu dmu",
            1e-12,
            || false,
        )
        .unwrap();
        let report = run(
            &model,
            [1000.0, 300.0],
            [300.0, 700.0],
            [0.0, 1.0],
            [1.0, 1.0],
        );
        let bg = blackbody(700.0);
        let inferred = (report.surface_heat_power[0].value() - (bg - blackbody(1000.0)))
            / (blackbody(300.0) - bg);
        errors.push((inferred - exact_transmission).abs());
        conserved(&report);
    }
    assert!(errors[1] < errors[0] / 15.0);
    assert!(errors[2] < errors[1] / 3.9);
    assert!(errors[2] < 3e-7);
}

#[test]
fn g0_two_layer_order_matches_both_direct_transfers_and_reversal() {
    let path = vec![segment(Layer::Lower, 0.4), segment(Layer::Upper, 0.8)];
    let model = plates(1.0, path.clone());
    let wall = [800.0, 350.0];
    let gas = [500.0, 650.0];
    let kappa = [0.2, 0.9];
    let report = run(&model, wall, gas, kappa, [0.55, 0.85]);
    let mut gas_heat = [0.0; 2];
    let mut incoming = report.surface_radiosity_w_m2[0];
    for part in &path {
        let i = part.layer.index();
        let absorbed = -(-kappa[i] * part.length.value()).exp_m1();
        let q = absorbed * (incoming - blackbody(gas[i]));
        gas_heat[i] += q;
        incoming -= q;
    }
    close(
        report.surface_heat_power[1].value(),
        incoming - report.surface_radiosity_w_m2[1],
    );
    incoming = report.surface_radiosity_w_m2[1];
    for part in path.iter().rev() {
        let i = part.layer.index();
        let absorbed = -(-kappa[i] * part.length.value()).exp_m1();
        let q = absorbed * (incoming - blackbody(gas[i]));
        gas_heat[i] += q;
        incoming -= q;
    }
    close(
        report.surface_heat_power[0].value(),
        incoming - report.surface_radiosity_w_m2[0],
    );
    for i in 0..2 {
        close(report.gas_heat_power[i].value(), gas_heat[i]);
    }
    let reversed = ParticipatingEnclosure::new(
        ViewFactorMatrix::infinite_parallel_plates(1.0).unwrap(),
        vec![OpticalPath {
            surfaces: [1, 0],
            exchange_area: Area::new(1.0),
            segments: path.into_iter().rev().collect(),
        }],
        "same reciprocal chord with reversed traversal",
        1e-12,
        || false,
    )
    .unwrap();
    let again = run(&reversed, wall, gas, kappa, [0.55, 0.85]);
    for (a, b) in report
        .surface_heat_power
        .iter()
        .chain(report.gas_heat_power.iter())
        .zip(
            again
                .surface_heat_power
                .iter()
                .chain(again.gas_heat_power.iter()),
        )
    {
        close(a.value(), b.value());
    }
    conserved(&report);
    conserved(&again);
}

#[test]
fn g0_equilibrium_is_exact_and_adjacent_temperatures_keep_heat() {
    let model = plates(
        1.0,
        vec![segment(Layer::Lower, 2.0), segment(Layer::Upper, 3.0)],
    );
    for kappa in [[0.0, 0.0], [1e-12, 0.7], [1e4, 2e4]] {
        let report = run(&model, [300.0; 2], [300.0; 2], kappa, [0.1, 0.8]);
        assert!(
            report
                .surface_heat_power
                .iter()
                .chain(report.gas_heat_power.iter())
                .all(|q| q.value() == 0.0)
        );
    }
    let hot = f64::from_bits(300.0_f64.to_bits() + 1);
    let report = run(&model, [hot, 300.0], [220.0, 1200.0], [0.0; 2], [1.0; 2]);
    let expected =
        STEFAN_BOLTZMANN_W_M2_K4 * (hot - 300.0) * (hot + 300.0) * (hot * hot + 300.0 * 300.0);
    assert!(report.surface_heat_power[1].value() > 0.0);
    assert!((report.surface_heat_power[1].value() / expected - 1.0).abs() < 1e-14);
}

#[test]
fn g1_opaque_layers_shield_remote_walls_and_exchange_at_the_interface() {
    let model = plates(
        1.0,
        vec![segment(Layer::Lower, 1.0), segment(Layer::Upper, 1.0)],
    );
    let report = run(
        &model,
        [1000.0, 300.0],
        [500.0, 700.0],
        [1e4, 1e4],
        [1.0, 1.0],
    );
    let b = [1000.0, 300.0, 500.0, 700.0].map(blackbody);
    close(report.surface_heat_power[0].value(), b[2] - b[0]);
    close(report.surface_heat_power[1].value(), b[3] - b[1]);
    close(report.gas_heat_power[0].value(), b[0] + b[3] - 2.0 * b[2]);
    close(report.gas_heat_power[1].value(), b[1] + b[2] - 2.0 * b[3]);
    conserved(&report);
}

#[test]
fn g3_subdividing_one_isothermal_segment_keeps_the_same_physics() {
    let whole = plates(1.0, vec![segment(Layer::Upper, 1.0)]);
    let split = plates(
        1.0,
        vec![segment(Layer::Upper, 0.25), segment(Layer::Upper, 0.75)],
    );
    for kappa in [1e-14, 0.7, 100.0] {
        let a = run(
            &whole,
            [900.0, 300.0],
            [300.0, 600.0],
            [0.0, kappa],
            [0.2, 0.7],
        );
        let b = run(
            &split,
            [900.0, 300.0],
            [300.0, 600.0],
            [0.0, kappa],
            [0.2, 0.7],
        );
        for (a, b) in a
            .surface_heat_power
            .iter()
            .chain(a.gas_heat_power.iter())
            .zip(b.surface_heat_power.iter().chain(b.gas_heat_power.iter()))
        {
            close(a.value(), b.value());
        }
        conserved(&a);
        conserved(&b);
    }
}

#[test]
fn g0_self_view_chords_count_both_directions_without_creating_heat() {
    let geometry = ViewFactorMatrix::admit(
        vec![1.0, 2.0],
        vec![vec![0.0, 1.0], vec![0.5, 0.5]],
        ViewFactorEvidence::Analytic {
            geometry: "reciprocal two-patch closure".into(),
        },
        ViewFactorTolerance::default(),
    )
    .unwrap();
    let paths = vec![
        OpticalPath {
            surfaces: [0, 1],
            exchange_area: Area::new(1.0),
            segments: vec![segment(Layer::Upper, 1.0)],
        },
        OpticalPath {
            surfaces: [1, 1],
            exchange_area: Area::new(0.5),
            segments: vec![segment(Layer::Upper, 1.0)],
        },
    ];
    let model =
        ParticipatingEnclosure::new(geometry, paths, "declared self-view fixture", 1e-12, || {
            false
        })
        .unwrap();
    let report = run(&model, [300.0; 2], [300.0, 600.0], [0.0, 0.5], [1.0; 2]);
    let q = -(-0.5_f64).exp_m1() * (blackbody(600.0) - blackbody(300.0));
    close(report.surface_heat_power[0].value(), q);
    close(report.surface_heat_power[1].value(), 2.0 * q);
    close(report.gas_heat_power[1].value(), -3.0 * q);
    conserved(&report);
}

#[test]
fn g4_invalid_geometry_source_domains_and_cancellation_publish_nothing() {
    assert!(LayerAbsorption::new(AbsorptionCoefficient::new(-1.0), "source").is_err());
    assert!(LayerAbsorption::new(AbsorptionCoefficient::new(0.0), " ").is_err());
    let broken = ParticipatingEnclosure::new(
        ViewFactorMatrix::infinite_parallel_plates(1.0).unwrap(),
        vec![OpticalPath {
            surfaces: [0, 1],
            exchange_area: Area::new(0.9),
            segments: vec![],
        }],
        "incomplete model",
        1e-12,
        || false,
    );
    assert!(broken.is_err());
    let model = plates(1.0, vec![segment(Layer::Upper, 1.0)]);
    let mut work = model.workspace().unwrap();
    let surfaces = [emissivity(0.8), emissivity(0.5)];
    for wall in [[1600.0, 300.0], [f64::NAN, 300.0], [0.0, 300.0]] {
        assert!(
            model
                .evaluate(
                    &wall.map(Temperature::new),
                    &surfaces,
                    [Temperature::new(600.0); 2],
                    &absorption([0.0, 0.5]),
                    RadiationSolve::default(),
                    &mut work,
                    || false
                )
                .is_err()
        );
    }
    let mut checks = 0;
    let result = model.evaluate(
        &[Temperature::new(300.0); 2],
        &surfaces,
        [Temperature::new(600.0); 2],
        &absorption([0.0, 0.5]),
        RadiationSolve::default(),
        &mut work,
        || {
            checks += 1;
            checks == 4
        },
    );
    assert_eq!(result, Err(RadiationError::Cancelled));
    let mut report = run(&model, [300.0; 2], [300.0, 600.0], [0.0, 0.5], [0.8, 0.5]);
    let mut rates = [LayerRates::default(); 2];
    report.add_gas_heat(&mut rates).unwrap();
    for i in 0..2 {
        close(
            rates[i].referenced_energy_w,
            report.gas_heat_power[i].value(),
        );
    }
    assert!(rates.iter().flat_map(|r| r.species_kg_s).all(|v| v == 0.0));
    let saved = rates;
    report.gas_heat_power = [Power::new(1.0), Power::new(f64::INFINITY)];
    assert!(report.add_gas_heat(&mut rates).is_err());
    assert_eq!(rates, saved);
}

// Reuse the existing original NASA9 air fixture instead of a second caloric bank.
#[path = "../../../examples/compartment_fire/air.rs"]
mod air;

#[test]
fn g1_original_gas_and_both_wall_stores_conserve_signed_radiation_energy() {
    use crate::compartment::wall::{
        GasContact, HeatTransferCoefficient, SpecificEnthalpy, ThermalConductivity, Wall1d,
        WallBoundary, WallGeometry, WallSolve,
    };
    use fs_gas::multicomponent::compartment::ZoneSolve;
    use fs_qty::{Density, Pressure, Time};
    use fs_thermal::enthalpy::{DepositionError, ThermalLaw};

    struct Sensible;
    impl ThermalLaw for Sensible {
        fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
            let t = 300.0 + h / 1000.0;
            if (200.0..=1500.0).contains(&t) {
                Ok(t)
            } else {
                Err(DepositionError::MaterialDomain)
            }
        }
        fn maximum_temperature_slope(&self) -> f64 {
            0.001
        }
    }
    let law = Sensible;
    let (gas, fractions) = air::source().unwrap();
    let mut room = TwoZone::new(
        gas,
        Area::new(2.0),
        Length::new(2.0),
        Length::new(1.0),
        Pressure::new(101_325.0),
        [Temperature::new(400.0), Temperature::new(600.0)],
        [fractions.clone(), fractions],
        || false,
    )
    .unwrap();
    let geometry = WallGeometry {
        area: Area::new(2.0),
        thickness: Length::new(0.04),
        density: Density::new(500.0),
        conductivity: ThermalConductivity::new(10.0),
        cells: 4,
    };
    let walls = [800.0, 300.0].map(|t| {
        Wall1d::new(
            &law,
            geometry,
            SpecificEnthalpy::new(1000.0 * (t - 300.0)),
            || false,
        )
        .unwrap()
    });
    let initial = room.physical_internal_energy().value()
        + walls.iter().map(|w| w.simulation().energy_j()).sum::<f64>();
    let initial_species =
        [Layer::Lower, Layer::Upper].map(|l| room.layer(l).species_masses().to_vec());
    let optics = plates(
        2.0,
        vec![segment(Layer::Lower, 1.0), segment(Layer::Upper, 1.0)],
    );
    // These authored coefficients exercise a gray absorbing numerical fixture.
    // They are not represented as measured nitrogen/oxygen absorption data.
    let radiation = optics
        .evaluate_room(
            &room,
            &[Temperature::new(800.0), Temperature::new(300.0)],
            &[emissivity(0.7), emissivity(0.8)],
            &absorption([0.4, 0.8]),
            RadiationSolve::default(),
            &mut optics.workspace().unwrap(),
            || false,
        )
        .unwrap();
    assert!(
        radiation.surface_heat_power[0].value() < 0.0,
        "hot wall must emit"
    );
    assert!(
        radiation.surface_heat_power[1].value() > 0.0,
        "cold wall must absorb"
    );
    let boundary = WallBoundary {
        gas: [GasContact {
            area: Area::new(0.0),
            temperature: Temperature::new(300.0),
            htc: HeatTransferCoefficient::new(0.0),
        }; 2],
        ambient_temperature: Temperature::new(300.0),
        ambient_htc: HeatTransferCoefficient::new(0.0),
        external_radiant_power: Power::new(0.0),
    };
    let duration = Time::new(0.01);
    let candidates = std::array::from_fn::<_, 2, _>(|i| {
        walls[i]
            .candidate_with_exchange(
                duration,
                boundary,
                radiation.surface_heat_power[i],
                WallSolve::default(),
                &mut walls[i].workspace().unwrap(),
                || false,
            )
            .unwrap()
    });
    let mut rates = [LayerRates::default(); 2];
    radiation.add_gas_heat(&mut rates).unwrap();
    room.advance(duration, rates, ZoneSolve::default(), || false)
        .unwrap();
    let final_energy = room.physical_internal_energy().value()
        + candidates
            .iter()
            .map(|c| c.wall.simulation().energy_j())
            .sum::<f64>();
    assert!(
        (final_energy - initial).abs() < 4e-5,
        "combined energy error {} J",
        final_energy - initial
    );
    for (i, candidate) in candidates.iter().enumerate() {
        assert_eq!(candidate.wall.time(), duration);
        assert_eq!(candidate.report.external_radiant_energy.value(), 0.0);
        close(
            candidate.report.coupled_radiant_energy.value(),
            duration.value() * radiation.surface_heat_power[i].value(),
        );
    }
    assert_eq!(room.time(), duration);
    for (i, layer) in [Layer::Lower, Layer::Upper].iter().enumerate() {
        assert_eq!(room.layer(*layer).species_masses(), initial_species[i]);
    }
}
