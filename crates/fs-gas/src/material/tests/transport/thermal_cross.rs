use super::*;
use crate::material::transport::{
    MaterialSoretDufourTransport, SpeciesHeatOfTransportMaterial, ThermalConductivityConvention,
};
use crate::multicomponent::NormalPrimitive;
use crate::multicomponent::dissipation::{MaxwellStefanSolveV1, PlanarTransport};
use fs_qty::Dims;
use fs_qty::semantic::CompositionBasis as CoordinateBasis;
use fs_thermochem::MassSpecificEnergyQuantityV1;

const MASS: f64 = 0.01;
const DIFFUSIVITY: f64 = 1e-4;

fn mole_quantity() -> QuantitySpec {
    quantity(QuantityKind::Composition(CoordinateBasis::MoleFraction))
}
fn reference() -> QueryPoint {
    context()
        .with("temperature", 400.0)
        .unwrap()
        .with("pressure", 100000.0)
        .unwrap()
        .with_quantity("x_A", mole_quantity(), 0.4)
        .unwrap()
        .with_quantity("x_B", mole_quantity(), 0.6)
        .unwrap()
}
fn axes() -> Vec<(SpeciesId, String)> {
    vec![(id("B"), "x_B".into()), (id("A"), "x_A".into())]
}
fn expected_heats(t: f64, p: f64, x: f64) -> [f64; 2] {
    [
        10000.0 + 100.0 * t + 40000.0 * x,
        -90000.0 + 0.2 * p - 20000.0 * (1.0 - x),
    ]
}
fn heat_claim(name: &str, maximum_temperature: f64) -> PropertyClaim {
    let (source_axes, dims, samples) = if name == "A" {
        (
            ["temperature".into(), "x_A".into()],
            [Temperature::DIMS, Dims::NONE],
            [
                (200.0, 0.0),
                (maximum_temperature, 0.0),
                (200.0, 1.0),
                (maximum_temperature, 1.0),
            ]
            .map(|(t, x)| ([t, x], expected_heats(t, 0.0, x)[0])),
        )
    } else {
        (
            ["pressure".into(), "x_B".into()],
            [Pressure::DIMS, Dims::NONE],
            [
                (50000.0, 0.0),
                (200000.0, 0.0),
                (50000.0, 1.0),
                (200000.0, 1.0),
            ]
            .map(|(p, xb)| ([p, xb], expected_heats(0.0, p, 1.0 - xb)[1])),
        )
    };
    PropertyClaim {
        key: PropertyKey::new(
            "kinetic_heat_of_transport",
            MassSpecificEnergyQuantityV1::DIMS,
        ),
        value: PropertyValue::Surface2D {
            axes: source_axes,
            axis_dims: dims,
            samples: samples.to_vec(),
            triangles: vec![[0, 1, 2], [1, 3, 2]],
            dims: MassSpecificEnergyQuantityV1::DIMS,
        },
        validity: ValidityDomain::unconstrained()
            .with_quantity(
                "temperature",
                quantity(QuantityKind::AbsoluteTemperature),
                200.0,
                maximum_temperature,
            )
            .with_quantity(
                "pressure",
                quantity(QuantityKind::Pressure),
                50000.0,
                200000.0,
            )
            .with_quantity("x_A", mole_quantity(), 0.0, 1.0)
            .with_quantity("x_B", mole_quantity(), 0.0, 1.0)
            .with_quantity(
                "humidity",
                quantity(QuantityKind::Moisture(MoistureBasis::RelativeHumidity)),
                0.0,
                1.0,
            ),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::LinearInside,
        observations: vec![],
        provenance: provenance(),
    }
}
fn heat_source(name: &str, claims: Vec<PropertyClaim>) -> SpeciesHeatOfTransportMaterial {
    let mut claim_set = ClaimSet::new();
    for claim in claims {
        claim_set.insert_claim(claim).unwrap();
    }
    SpeciesHeatOfTransportMaterial {
        species: id(name),
        material: MaterialCard::assemble(
            MaterialStateId {
                chemistry: name.into(),
                phase: "gas".into(),
                process: "manufactured-kinetic-cross-coefficient".into(),
                revision: 0,
            },
            claim_set,
            vec![],
        )
        .unwrap(),
        property: PropertyKey::new(
            "kinetic_heat_of_transport",
            MassSpecificEnergyQuantityV1::DIMS,
        ),
        selection: MaterialPropertySelection::SingleClaimOnly,
    }
}
fn sources() -> Vec<SpeciesHeatOfTransportMaterial> {
    vec![
        heat_source("B", vec![heat_claim("B", 800.0)]),
        heat_source("A", vec![heat_claim("A", 800.0)]),
    ]
}
fn sourced() -> MaterialSoretDufourTransport<MaterialTransport<ConstantTransport>> {
    MaterialSoretDufourTransport::new(
        adapter(model(0.01, 0.03, 2000.0)),
        ThermalConductivityConvention::ZeroDiffusiveSpeciesFlux,
        sources(),
        axes(),
        reference(),
        || false,
    )
    .unwrap()
}
fn relative(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-10 * expected.abs().max(1e-18),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn g1_material_soret_surfaces_drive_shared_faces_at_actual_temperature_pressure_and_composition() {
    let law = sourced();
    let gas = MixtureGas::new(
        law.source().clone(),
        CaloricDomain::new(Temperature::new(200.0), Temperature::new(1900.0)),
        || false,
    )
    .unwrap();
    let mut face = PlanarTransport::new(gas, law, || false).unwrap();
    // Independently vary T, pressure and composition. Equal-mass isomers have
    // no pressure diffusion or partial-enthalpy contrast, exposing ONLY the
    // source-card Soret/Dufour term alongside the material Fourier flux.
    for (t, p, x, gx, gt) in [
        (400.0, 100000.0, 0.4, 0.0, 200.0),
        (600.0, 100000.0, 0.4, 0.0, 200.0),
        (400.0, 170000.0, 0.4, 0.0, 200.0),
        (400.0, 100000.0, 0.7, 0.0, 200.0),
        (500.0, 140000.0, 0.6, 0.2, 0.0),
        (500.0, 140000.0, 0.6, 0.2, 200.0),
    ] {
        let distance = 0.1;
        let densities = |x: f64, t: f64| [p * x * MASS / (R * t), p * (1.0 - x) * MASS / (R * t)];
        let left = densities(x - 0.5 * distance * gx, t - 0.5 * distance * gt);
        let right = densities(x + 0.5 * distance * gx, t + 0.5 * distance * gt);
        let primitive = |densities| NormalPrimitive {
            partial_densities: densities,
            normal_velocity: 0.0,
            transverse_velocity: 0.0,
            pressure: p,
        };
        let mut flux = [0.0; 5];
        let report = face
            .face_from_primitives(
                primitive(&left),
                primitive(&right),
                distance,
                [0.0; 2],
                MaxwellStefanSolveV1::default(),
                &mut flux,
                || false,
            )
            .unwrap();
        let heat = expected_heats(t, p, x);
        let drive = gx + x * (1.0 - x) * MASS * (heat[0] - heat[1]) / (R * t) * (gt / t);
        let j = -p / (R * t) * DIFFUSIVITY * MASS * drive;
        relative(flux[0], j);
        relative(flux[1], -j);
        // Independent two-component WMS expression for equal molar masses.
        let phi_ab = (1.0 + (1.0_f64 / 3.0).sqrt()).powi(2) / 4.0;
        let phi_ba = (1.0 + 3.0_f64.sqrt()).powi(2) / 4.0;
        let k = x * 0.03 / (x + (1.0 - x) * phi_ab) + (1.0 - x) * 0.15 / ((1.0 - x) + x * phi_ba);
        let reduced_heat = -k * gt + (heat[0] - heat[1]) * j;
        relative(flux[4], reduced_heat);
        let entropy = -R * j / MASS * gx / (x * (1.0 - x)) - reduced_heat * gt / (t * t);
        let pair_entropy = R * j * j / (p / (R * t) * DIFFUSIVITY * MASS * MASS * x * (1.0 - x));
        relative(entropy, pair_entropy + k * (gt / t).powi(2));
        assert!(entropy > 0.0);
        assert!(report.mass_frame_residual_kg_m2_s.abs() < 1e-16);
        let resolved = face.law().last_resolved();
        for (value, expected) in resolved.heats().iter().zip(heat) {
            relative(value.value(), expected);
        }
        for (source, material) in face.law().sources().iter().zip(resolved.components()) {
            assert_eq!(material.card_identity(), source.material.content_hash());
            let query = material.query();
            relative(query.axes()["temperature"], t);
            relative(query.axes()["pressure"], p);
            relative(query.axes()["x_A"], x);
            assert_eq!(query.axes()["sample_id"], 7.0);
            assert_eq!(query.axis_quantities()["x_A"], mole_quantity());
            assert_eq!(query.axes()["humidity"], 0.3);
        }
    }
}

#[test]
fn g0_material_soret_pins_selected_claims_and_refuses_missing_domains_without_publication() {
    let mut original = heat_claim("A", 500.0);
    // The immutable card has a tempting disjoint replacement at high T. The
    // original selected claim must refuse there, even though a fresh selection
    // against this exact card could succeed.
    let mut replacement = heat_claim("A", 800.0);
    replacement.validity = replacement.validity.with_quantity(
        "temperature",
        quantity(QuantityKind::AbsoluteTemperature),
        600.0,
        800.0,
    );
    if let PropertyValue::Surface2D { samples, .. } = &mut replacement.value {
        for (point, value) in samples {
            point[0] = point[0].max(600.0);
            *value = expected_heats(point[0], 0.0, point[1])[0];
        }
    }
    replacement.provenance.source = "different experimental cross-coefficient source".into();
    original.provenance.source = "original lower-temperature source".into();
    let mut bound = sources();
    bound[1] = heat_source("A", vec![original, replacement]);
    assert!(
        MaterialSoretDufourTransport::new(
            adapter(model(0.01, 0.03, 2000.0)),
            ThermalConductivityConvention::ZeroDiffusiveSpeciesFlux,
            bound.clone(),
            axes(),
            reference().with("temperature", 700.0).unwrap(),
            || false,
        )
        .is_ok()
    );
    let mut law = MaterialSoretDufourTransport::new(
        adapter(model(0.01, 0.03, 2000.0)),
        ThermalConductivityConvention::ZeroDiffusiveSpeciesFlux,
        bound,
        axes(),
        reference(),
        || false,
    )
    .unwrap();
    let original = law.last_resolved().clone();
    let mut output = [11.0; 2];
    for (t, p, x) in [
        (700.0, 100000.0, [0.4, 0.6]),
        (400.0, 300000.0, [0.4, 0.6]),
        (400.0, 100000.0, [0.3, 0.6]),
    ] {
        assert!(
            law.heats_of_transport(
                Temperature::new(t),
                Pressure::new(p),
                &x,
                &mut output,
                &mut || false
            )
            .is_err()
        );
        assert_eq!(output, [11.0; 2]);
        assert_eq!(law.last_resolved(), &original);
    }
    let mut binary = [9.0; 4];
    assert!(
        law.evaluate(
            Temperature::new(700.0),
            Pressure::new(100000.0),
            &[0.4, 0.6],
            &mut binary,
            &mut || false
        )
        .is_err()
    );
    assert_eq!(binary, [9.0; 4]);
    assert_eq!(law.last_resolved(), &original);
    law.heats_of_transport(
        Temperature::new(450.0),
        Pressure::new(130000.0),
        &[1.0, 0.0],
        &mut output,
        &mut || false,
    )
    .unwrap();
    assert_eq!(law.last_resolved().mole_fractions(), &[1.0, 0.0]);
    assert_eq!(law.last_resolved().components().len(), 2);
    for (previous, refreshed) in original
        .components()
        .iter()
        .zip(law.last_resolved().components())
    {
        assert_eq!(
            previous.properties()[0].answer().receipt.selected,
            refreshed.properties()[0].answer().receipt.selected
        );
    }
}

#[test]
fn g0_material_soret_requires_explicit_species_units_coordinates_and_complete_source_data() {
    let construct = |sources, coordinates, point| {
        MaterialSoretDufourTransport::new(
            adapter(model(0.01, 0.03, 2000.0)),
            ThermalConductivityConvention::ZeroDiffusiveSpeciesFlux,
            sources,
            coordinates,
            point,
            || false,
        )
    };
    let mut missing = sources();
    missing.pop();
    assert!(construct(missing, axes(), reference()).is_err());
    let mut duplicate = sources();
    duplicate[0].species = id("A");
    assert!(construct(duplicate, axes(), reference()).is_err());
    let mut wrong_units = sources();
    wrong_units[0].property = PropertyKey::new("kinetic_heat_of_transport", fs_qty::Energy::DIMS);
    assert!(construct(wrong_units, axes(), reference()).is_err());
    let mut missing_property = sources();
    missing_property[0].property = PropertyKey::new(
        "missing_kinetic_cross_coefficient",
        MassSpecificEnergyQuantityV1::DIMS,
    );
    assert!(construct(missing_property, axes(), reference()).is_err());
    assert!(
        construct(
            sources(),
            vec![(id("A"), "x_A".into()), (id("B"), "x_A".into())],
            reference()
        )
        .is_err()
    );
    assert!(
        construct(
            sources(),
            axes(),
            reference()
                .with_quantity(
                    "x_A",
                    quantity(QuantityKind::Composition(CoordinateBasis::MassFraction)),
                    0.4
                )
                .unwrap()
        )
        .is_err()
    );
    assert!(
        construct(
            sources(),
            axes(),
            reference()
                .with_quantity(
                    "stale_temperature",
                    quantity(QuantityKind::AbsoluteTemperature),
                    400.0
                )
                .unwrap()
        )
        .is_err()
    );
    let frequency = quantity(QuantityKind::Frequency(
        fs_qty::semantic::FrequencyConvention::Cyclic,
    ));
    let mut dynamic = heat_claim("A", 800.0);
    dynamic.validity = dynamic
        .validity
        .with_quantity("frequency", frequency, 0.0, 1000.0);
    let mut dynamic_sources = sources();
    dynamic_sources[1] = heat_source("A", vec![dynamic]);
    assert!(
        construct(
            dynamic_sources,
            axes(),
            reference()
                .with_quantity("frequency", frequency, 10.0)
                .unwrap()
        )
        .is_err()
    );
    // The shared static classifier admits legacy untyped T. This consumer
    // evolves ONLY `temperature`, so a T-dependent source cannot be frozen
    // at an independently supplied initial coordinate.
    let mut legacy = heat_claim("A", 800.0);
    legacy.value = PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: Temperature::DIMS,
        knots: vec![(200.0, 30000.0), (800.0, 90000.0)],
        dims: MassSpecificEnergyQuantityV1::DIMS,
    };
    legacy.validity = ValidityDomain::unconstrained().with("T", 200.0, 800.0);
    let mut legacy_sources = sources();
    legacy_sources[1] = heat_source("A", vec![legacy]);
    assert!(
        construct(
            legacy_sources,
            axes(),
            reference().with("T", 400.0).unwrap(),
        )
        .is_err()
    );
    // Merely having card viscosity/conductivity or caloric h never activates H.
    let mut ordinary = adapter(model(0.01, 0.03, 2000.0));
    let mut heats = [11.0; 2];
    assert!(
        !ordinary
            .heats_of_transport(
                Temperature::new(400.0),
                Pressure::new(100000.0),
                &[0.4, 0.6],
                &mut heats,
                &mut || false
            )
            .unwrap()
    );
    assert_eq!(heats, [0.0; 2]);
}

#[test]
fn g4_material_soret_cancelled_refresh_and_late_binary_refusal_keep_original_receipts() {
    let mut law = sourced();
    let original = law.last_resolved().clone();
    let mut polls = 0;
    let mut values = [0.0; 2];
    law.heats_of_transport(
        Temperature::new(500.0),
        Pressure::new(150000.0),
        &[0.7, 0.3],
        &mut values,
        &mut || {
            polls += 1;
            false
        },
    )
    .unwrap();
    assert!(polls >= 5);
    let mut cancelled = sourced();
    let mut count = 0;
    let mut output = [13.0; 2];
    assert!(matches!(
        cancelled.heats_of_transport(
            Temperature::new(500.0),
            Pressure::new(150000.0),
            &[0.7, 0.3],
            &mut output,
            &mut || {
                count += 1;
                count == polls
            }
        ),
        Err(DissipationError::Cancelled)
    ));
    assert_eq!(output, [13.0; 2]);
    assert_eq!(cancelled.last_resolved(), &original);

    let material = model(0.01, 0.03, 2000.0);
    let inner = MaterialTransport::new(
        material.clone(),
        context(),
        ObservedDiffusion {
            inner: binary(&material),
            calls: Rc::new(Cell::new(0)),
            refuse: true,
        },
        DynViscosity::new(0.0),
        || false,
    )
    .unwrap();
    let mut law = MaterialSoretDufourTransport::new(
        inner,
        ThermalConductivityConvention::ZeroDiffusiveSpeciesFlux,
        sources(),
        axes(),
        reference(),
        || false,
    )
    .unwrap();
    let original = law.last_resolved().clone();
    let mut binary = [19.0; 4];
    assert!(
        law.evaluate(
            Temperature::new(500.0),
            Pressure::new(150000.0),
            &[0.7, 0.3],
            &mut binary,
            &mut || false
        )
        .is_err()
    );
    assert_eq!(binary, [19.0; 4]);
    assert_eq!(law.last_resolved(), &original);
}
