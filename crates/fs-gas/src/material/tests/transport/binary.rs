use super::*;
use crate::material::transport::MaterialBinaryTransport;
use crate::multicomponent::NormalPrimitive;
use crate::multicomponent::dissipation::{MaxwellStefanSolveV1, PlanarTransport};
use fs_material::thermochemical_gas::energy::parcel::transfer::species::diffusion::{
    BINARY_DIFFUSIVITY_DIMS, BINARY_DIFFUSIVITY_PROPERTY, BinaryDiffusionMaterial,
};

fn reference() -> QueryPoint {
    context()
        .with("temperature", 300.0)
        .unwrap()
        .with("pressure", 80000.0)
        .unwrap()
}

fn diffusion_claim(d: f64, bounds: [f64; 2], varying: bool) -> PropertyClaim {
    let dims = BINARY_DIFFUSIVITY_DIMS;
    let value = if varying {
        PropertyValue::Surface2D {
            axes: ["temperature".into(), "pressure".into()],
            axis_dims: [Temperature::DIMS, Pressure::DIMS],
            samples: [
                [bounds[0], 50000.0],
                [bounds[1], 50000.0],
                [bounds[0], 400000.0],
                [bounds[1], 400000.0],
            ]
            .map(|p| (p, d * (p[0] / 400.0 + p[1] / 200000.0)))
            .to_vec(),
            triangles: vec![[0, 1, 2], [1, 3, 2]],
            dims,
        }
    } else {
        PropertyValue::Scalar { value: d, dims }
    };
    PropertyClaim {
        key: PropertyKey::new(BINARY_DIFFUSIVITY_PROPERTY, dims),
        value,
        validity: ValidityDomain::unconstrained()
            .with_quantity(
                "temperature",
                quantity(QuantityKind::AbsoluteTemperature),
                bounds[0],
                bounds[1],
            )
            .with_quantity(
                "pressure",
                quantity(QuantityKind::Pressure),
                50000.0,
                400000.0,
            )
            .with_quantity(
                "humidity",
                quantity(QuantityKind::Moisture(MoistureBasis::RelativeHumidity)),
                0.3,
                0.3,
            ),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: if varying {
            InterpolationPolicy::LinearInside
        } else {
            InterpolationPolicy::ConstantWithinValidity
        },
        observations: vec![],
        provenance: provenance(),
    }
}

fn diffusion_card(claims: Vec<PropertyClaim>) -> MaterialCard {
    let mut set = ClaimSet::new();
    for claim in claims {
        set.insert_claim(claim).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured binary pair".into(),
            phase: "gas".into(),
            process: "molecular diffusion fixture".into(),
            revision: 0,
        },
        set,
        vec![],
    )
    .unwrap()
}

fn pair<'a>(
    a: &str,
    b: &str,
    card: &'a MaterialCard,
    context: &'a QueryPoint,
) -> BinaryDiffusionMaterial<'a> {
    BinaryDiffusionMaterial {
        species: [id(a), id(b)],
        material: card,
        context,
        selection: MaterialPropertySelection::SingleClaimOnly,
    }
}

fn sourced(
    material: ThermochemicalGasModel,
    pairs: &[BinaryDiffusionMaterial<'_>],
) -> Result<MaterialBinaryTransport, DissipationError> {
    let source = binary(&material).source().clone();
    MaterialBinaryTransport::new(
        material,
        context(),
        source,
        pairs,
        DynViscosity::new(0.0),
        || false,
    )
}

fn original_material(unequal: bool, ternary: bool) -> ThermochemicalGasModel {
    let mut b = species("B", -2e6);
    if unequal {
        let mut card = b.regions()[0].card().clone();
        let key = required_nasa9_parameter_dims()[2].0;
        card.parameters.get_mut(key).unwrap().value = 6.0;
        b = Nasa9StandardStateModelV1::new(
            id("B"),
            MolarMass::new(0.02),
            b.convention().clone(),
            vec![Nasa9RegionV1::from_card(card).unwrap()],
        )
        .unwrap();
    }
    let mut models = vec![species("A", 1e6), b];
    let mut cards = vec![card("A", 1e-5, 0.03, 2000.0), card("B", 2e-5, 0.06, 2000.0)];
    let mut fractions = vec![0.5, 0.5];
    if ternary {
        models.push(species("C", 3e6));
        cards.push(card("C", 3e-5, 0.09, 2000.0));
        fractions = vec![0.25, 0.25, 0.5];
    }
    ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            models,
            Composition::new(CompositionBasis::MoleFraction, fractions).unwrap(),
        )
        .unwrap(),
        cards,
        None,
    )
    .unwrap()
}

fn trace(density: &[f64], pressure: f64) -> NormalPrimitive<'_> {
    NormalPrimitive {
        partial_densities: density,
        normal_velocity: 0.0,
        transverse_velocity: 0.0,
        pressure,
    }
}

#[test]
fn g1_binary_material_faces_use_actual_pair_surfaces_mass_frame_and_partial_enthalpy() {
    let material = original_material(true, false);
    let card = diffusion_card(vec![diffusion_claim(0.02, [200.0, 900.0], true)]);
    let point = reference();
    let law = sourced(material, &[pair("B", "A", &card, &point)]).unwrap();
    let gas = MixtureGas::new(
        law.source().clone(),
        CaloricDomain::new(Temperature::new(200.0), Temperature::new(1900.0)),
        || false,
    )
    .unwrap();
    let mut face = PlanarTransport::new(gas, law, || false).unwrap();
    let concentration = 100000.0 / (R * 400.0);
    let left = [concentration * 0.8 * 0.01, concentration * 0.2 * 0.02];
    let right = [concentration * 0.2 * 0.01, concentration * 0.8 * 0.02];
    let mut flux = [0.0; 5];
    let report = face
        .face_from_primitives(
            trace(&left, 100000.0),
            trace(&right, 100000.0),
            0.25,
            [0.0; 2],
            MaxwellStefanSolveV1::default(),
            &mut flux,
            || false,
        )
        .unwrap();
    let d = 0.02 * 1.5; // Actual T=400 K, p=100 kPa; not the reference query.
    let ja = -concentration * d * (0.01 * 0.02 / 0.015) * ((0.2 - 0.8) / 0.25);
    close(flux[0], ja, 1e-12);
    close(flux[1], -ja, 1e-12);
    assert!(report.mass_frame_residual_kg_m2_s.abs() < 1e-14);
    // Original species-referenced h, with each source's own Cv and molar mass.
    let ha = R * (3.0 * (400.0 - 200.0) + 400.0) / 0.01;
    let hb = R * (5.0 * (400.0 - 200.0) + 400.0) / 0.02;
    close(flux[4], ja * (ha - hb), 2e-11);
    assert!(flux[4].abs() > 100.0);
    let selected = &face.law().last_resolved()[0];
    close(selected.diffusivity_m2_s, d, 1e-15);
    close(
        selected
            .material
            .query_point()
            .iter()
            .find(|(axis, _)| axis == "temperature")
            .unwrap()
            .1,
        400.0,
        1e-13,
    );
    assert!(
        selected
            .material
            .query_point()
            .contains(&("sample_id".into(), 7.0))
    );
    let mut reverse = [0.0; 5];
    face.face_from_primitives(
        trace(&right, 100000.0),
        trace(&left, 100000.0),
        0.25,
        [0.0; 2],
        MaxwellStefanSolveV1::default(),
        &mut reverse,
        || false,
    )
    .unwrap();
    for (a, b) in flux.iter().zip(reverse) {
        close(*a, -b, 1e-12);
    }
    // Equal composition at unequal pressure still separates unequal masses.
    // The shared mass-frame pressure force is (x_A-Y_A) grad(ln p).
    let low = [
        concentration * 0.95 * 0.5 * 0.01,
        concentration * 0.95 * 0.5 * 0.02,
    ];
    let high = [
        concentration * 1.05 * 0.5 * 0.01,
        concentration * 1.05 * 0.5 * 0.02,
    ];
    face.face_from_primitives(
        trace(&low, 95000.0),
        trace(&high, 105000.0),
        0.25,
        [0.0; 2],
        MaxwellStefanSolveV1::default(),
        &mut flux,
        || false,
    )
    .unwrap();
    let pressure_force = (0.5 - 1.0 / 3.0) * (105000.0_f64 / 95000.0).ln() / 0.25;
    let pressure_ja = -concentration * d * (0.01 * 0.02 / 0.015) * pressure_force;
    close(flux[0], pressure_ja, 1e-12);
    close(flux[1], -pressure_ja, 1e-12);
    close(flux[4], pressure_ja * (ha - hb), 1e-11);
}

#[test]
fn g1_binary_cards_feed_reciprocal_thermal_diffusion_without_replacing_scalar_materials() {
    use crate::multicomponent::dissipation::SoretDufourTransport;
    use fs_thermochem::MassSpecificEnergyQuantityV1;

    let material = original_material(false, false);
    let card = diffusion_card(vec![diffusion_claim(0.02, [200.0, 900.0], false)]);
    let point = reference();
    let molecular = sourced(material, &[pair("A", "B", &card, &point)]).unwrap();
    let heat = 1e5;
    let law = SoretDufourTransport::new(
        molecular,
        "manufactured-cross-data",
        1,
        TransportDomain {
            minimum_temperature_k: 200.0,
            maximum_temperature_k: 900.0,
            minimum_pressure_pa: 50000.0,
            maximum_pressure_pa: 400000.0,
        },
        vec![
            (id("A"), MassSpecificEnergyQuantityV1::new(heat)),
            (id("B"), MassSpecificEnergyQuantityV1::new(-heat)),
        ],
    )
    .unwrap();
    let gas = MixtureGas::new(
        law.source().clone(),
        CaloricDomain::new(Temperature::new(200.0), Temperature::new(1900.0)),
        || false,
    )
    .unwrap();
    let mut face = PlanarTransport::new(gas, law, || false).unwrap();
    let left = [0.005 * 100000.0 / (R * 350.0); 2];
    let right = [0.005 * 100000.0 / (R * 450.0); 2];
    let mut flux = [0.0; 5];
    let report = face
        .face_from_primitives(
            trace(&left, 100000.0),
            trace(&right, 100000.0),
            0.25,
            [0.0; 2],
            MaxwellStefanSolveV1::default(),
            &mut flux,
            || false,
        )
        .unwrap();
    let gradient = (450.0 - 350.0) / 0.25;
    let force = 0.5 * 0.01 * heat / (R * 400.0) * (gradient / 400.0);
    let ja = -(100000.0 / (R * 400.0)) * 0.02 * 0.01 * force;
    close(flux[0], ja, 1e-12);
    close(flux[1], -ja, 1e-12);
    let reduced = -report.thermal_conductivity_w_m_k * gradient + 2.0 * heat * ja;
    close(report.reduced_heat_flux_w_m2, reduced, 1e-11);
    close(flux[4], reduced, 1e-11); // Equal calorics cancel the partial-h pair.
    assert!(report.thermal_conductivity_w_m_k > 0.03);
    assert!(report.thermal_conductivity_w_m_k < 0.06);
}

fn diffusion_mode(n: usize, diffusion: f64) -> (f64, f64) {
    let material = original_material(false, false);
    let card = diffusion_card(vec![diffusion_claim(diffusion, [200.0, 900.0], false)]);
    let point = reference();
    let law = sourced(material, &[pair("A", "B", &card, &point)]).unwrap();
    let gas = MixtureGas::new(
        law.source().clone(),
        CaloricDomain::new(Temperature::new(200.0), Temperature::new(1900.0)),
        || false,
    )
    .unwrap();
    let wave = |x: f64| (2.0 * std::f64::consts::PI * x).sin();
    let states = (0..n)
        .map(|i| {
            let a = 0.5 + 0.2 * wave((i as f64 + 0.5) / n as f64);
            gas.conserved(&[a, 1.0 - a], 0.0, R / 0.01 * 400.0).unwrap()
        })
        .collect();
    let mut flow = MixtureEuler1d::new(
        gas,
        Length::new(1.0 / n as f64),
        0.7,
        states,
        Boundaries::Periodic,
        || false,
    )
    .unwrap();
    let before = flow.totals();
    let mut operator = DissipationOperator::new(&flow, law, || false).unwrap();
    let duration = 0.2;
    let report = operator
        .advance_for(
            &mut flow,
            Time::new(duration),
            DissipationControl::default(),
            || false,
        )
        .unwrap();
    let amplitude = 0.2 * (-diffusion * (2.0 * std::f64::consts::PI).powi(2) * duration).exp();
    let mut error = 0.0;
    for i in 0..n {
        let c = flow.cell(i).unwrap();
        error += (c.partial_densities[0] - (0.5 + amplitude * wave((i as f64 + 0.5) / n as f64)))
            .abs()
            / n as f64;
        close(c.partial_densities.iter().sum(), 1.0, 1e-13);
        close(
            flow.gas().primitive(c, || false).unwrap().temperature,
            400.0,
            1e-10,
        );
    }
    for (a, b) in flow
        .totals()
        .species_masses
        .iter()
        .zip(&before.species_masses)
    {
        close(*a, *b, 1e-13);
    }
    close(flow.totals().energy, before.energy, 1e-12);
    assert!(report.absolute_energy_residual_j.abs() < 1e-6);
    (error, flow.cell(n / 4).unwrap().partial_densities[0])
}

#[test]
fn g1_g3_source_binary_diffusion_drives_conserved_field_refinement_and_data_response() {
    let coarse = diffusion_mode(16, 0.02);
    let fine = diffusion_mode(32, 0.02);
    let faster = diffusion_mode(32, 0.04);
    assert!(coarse.0 < 3e-4);
    assert!(
        fine.0 < coarse.0 * 0.3,
        "coarse={}, fine={}",
        coarse.0,
        fine.0
    );
    assert!(faster.1 < fine.1 - 0.02);
}

#[test]
fn g0_pair_axis_is_complete_symmetric_and_independent_of_assignment_order() {
    let material = original_material(false, true);
    let cards =
        [1e-5, 2e-5, 4e-5].map(|d| diffusion_card(vec![diffusion_claim(d, [200.0, 900.0], false)]));
    let point = reference();
    let pairs = [
        pair("A", "B", &cards[0], &point),
        pair("A", "C", &cards[1], &point),
        pair("B", "C", &cards[2], &point),
    ];
    let permuted = [
        pair("C", "B", &cards[2], &point),
        pair("B", "A", &cards[0], &point),
        pair("C", "A", &cards[1], &point),
    ];
    let mut a = sourced(material.clone(), &pairs).unwrap();
    let mut b = sourced(material.clone(), &permuted).unwrap();
    let mut da = [9.0; 9];
    let mut db = da;
    let ca = a
        .evaluate(
            Temperature::new(400.0),
            Pressure::new(100000.0),
            &[1.0, 0.0, 0.0],
            &mut da,
            &mut || false,
        )
        .unwrap();
    let cb = b
        .evaluate(
            Temperature::new(400.0),
            Pressure::new(100000.0),
            &[1.0, 0.0, 0.0],
            &mut db,
            &mut || false,
        )
        .unwrap();
    assert_eq!(da, [0.0, 1e-5, 2e-5, 1e-5, 0.0, 4e-5, 2e-5, 4e-5, 0.0]);
    assert_eq!(da, db);
    assert_eq!(ca, cb);
    assert_eq!(a.last_resolved(), b.last_resolved());
    assert_eq!(ca.dynamic_viscosity.value(), 1e-5);
    assert_eq!(ca.thermal_conductivity.value(), 0.03);
    assert!(sourced(material.clone(), &pairs[..2]).is_err());
    for bad in [
        pair("A", "A", &cards[0], &point),
        pair("A", "Z", &cards[0], &point),
        pairs[1].clone(),
    ] {
        assert!(sourced(material.clone(), &[bad, pairs[1].clone(), pairs[2].clone()]).is_err());
    }
    let mut wrong = diffusion_claim(1e-5, [200.0, 900.0], false);
    wrong.key = PropertyKey::new(BINARY_DIFFUSIVITY_PROPERTY, DynViscosity::DIMS);
    wrong.value = PropertyValue::Scalar {
        value: 1e-5,
        dims: DynViscosity::DIMS,
    };
    let wrong = diffusion_card(vec![wrong]);
    assert!(
        sourced(
            material,
            &[
                pair("A", "B", &wrong, &point),
                pairs[1].clone(),
                pairs[2].clone()
            ]
        )
        .is_err()
    );
}

#[test]
fn g4_pair_source_support_remains_pinned_and_all_refusals_preserve_publication() {
    let material = original_material(false, false);
    let mut replacement = diffusion_claim(0.04, [600.0, 900.0], false);
    replacement.provenance.source = "independent hot material data".into();
    let card = diffusion_card(vec![
        diffusion_claim(0.02, [200.0, 500.0], false),
        replacement,
    ]);
    let point = reference();
    let pairs = [pair("A", "B", &card, &point)];
    let mut law = sourced(material.clone(), &pairs).unwrap();
    let before = law.last_resolved().to_vec();
    for (t, p, x) in [
        (700.0, 100000.0, [0.4, 0.6]),
        (400.0, 500000.0, [0.4, 0.6]),
        (400.0, 100000.0, [0.4, 0.5]),
    ] {
        let mut output = [11.0; 4];
        assert!(
            law.evaluate(
                Temperature::new(t),
                Pressure::new(p),
                &x,
                &mut output,
                &mut || false
            )
            .is_err()
        );
        assert_eq!(output, [11.0; 4]);
        assert_eq!(law.last_resolved(), &before);
    }
    // The hot claim is selectable at a separate explicit binding, not by
    // falling through a failed refresh of the original lower-temperature one.
    let hot = point.clone().with("temperature", 700.0).unwrap();
    assert!(sourced(material.clone(), &[pair("A", "B", &card, &hot)]).is_ok());
    let mut calls = 0;
    law.evaluate(
        Temperature::new(400.0),
        Pressure::new(100000.0),
        &[0.4, 0.6],
        &mut [0.0; 4],
        &mut || {
            calls += 1;
            false
        },
    )
    .unwrap();
    for stop in [1, calls / 2, calls] {
        let mut law = sourced(material.clone(), &pairs).unwrap();
        let mut output = [11.0; 4];
        let mut count = 0;
        assert_eq!(
            law.evaluate(
                Temperature::new(400.0),
                Pressure::new(100000.0),
                &[0.4, 0.6],
                &mut output,
                &mut || {
                    count += 1;
                    count == stop
                }
            ),
            Err(DissipationError::Cancelled)
        );
        assert_eq!(output, [11.0; 4]);
        assert_eq!(law.last_resolved(), &before);
    }
}

#[test]
fn g0_binary_source_rejects_unowned_or_untyped_physical_support() {
    let material = original_material(false, false);
    let point = reference();
    let claim = diffusion_claim(0.02, [200.0, 900.0], false);
    let card = diffusion_card(vec![claim.clone()]);
    for (axis, kind, value) in [
        ("T", QuantityKind::AbsoluteTemperature, 300.0),
        ("p", QuantityKind::Pressure, 80000.0),
        (
            "X_A",
            QuantityKind::Composition(fs_qty::semantic::CompositionBasis::MoleFraction),
            0.5,
        ),
    ] {
        let ambiguous = point
            .clone()
            .with_quantity(axis, quantity(kind), value)
            .unwrap();
        assert!(sourced(material.clone(), &[pair("A", "B", &card, &ambiguous)]).is_err());
    }
    // An untyped identity can be retained in the context and query receipt,
    // but cannot stand in for physical applicability semantics in the source.
    let mut untyped = claim.clone();
    untyped.validity = untyped.validity.with("sample_id", 7.0, 7.0);
    let untyped = diffusion_card(vec![untyped]);
    assert!(sourced(material.clone(), &[pair("A", "B", &untyped, &point)]).is_err());

    let frequency = quantity(QuantityKind::Frequency(
        fs_qty::semantic::FrequencyConvention::Cyclic,
    ));
    let mut dynamic = claim;
    dynamic.validity = dynamic
        .validity
        .with_quantity("frequency", frequency, 0.0, 100.0);
    let dynamic = diffusion_card(vec![dynamic]);
    let dynamic_point = point.with_quantity("frequency", frequency, 1.0).unwrap();
    assert!(sourced(material, &[pair("A", "B", &dynamic, &dynamic_point)]).is_err());
}
