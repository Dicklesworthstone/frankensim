use super::*;
use crate::material::transport::MaterialTransport;
use crate::multicomponent::dissipation::{
    ConstantTransport, DissipationControl, DissipationError, DissipationOperator,
    TransportCoefficients, TransportDomain, TransportLaw,
};
use crate::multicomponent::reactive::{ReactingEuler1d, ReactiveError, ReactiveStepControl};
use crate::multicomponent::{Boundaries, CaloricDomain, MixtureEuler1d, MixtureGas};
use fs_qty::Length;
use fs_qty::semantic::MoistureBasis;
use fs_thermochem::transport::kinetic::ThermalConductivityQuantityV1;
use std::cell::Cell;
use std::rc::Rc;

fn card(name: &str, mu: f64, k: f64, maximum_temperature: f64) -> SpeciesTransportMaterial {
    let mut claims = ClaimSet::new();
    for (property, dims, value) in [
        ("dynamic_viscosity", DynViscosity::DIMS, mu),
        (
            "thermal_conductivity",
            fs_material::state_point::THERMAL_CONDUCTIVITY_DIMS,
            k,
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(property, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "temperature",
                        quantity(QuantityKind::AbsoluteTemperature),
                        200.0,
                        maximum_temperature,
                    )
                    .with_quantity("pressure", quantity(QuantityKind::Pressure), 100.0, 1e8)
                    .with("sample_id", 7.0, 7.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: provenance(),
            })
            .unwrap();
    }
    SpeciesTransportMaterial {
        species: id(name),
        material: MaterialCard::assemble(
            MaterialStateId {
                chemistry: name.into(),
                phase: "gas".into(),
                process: "manufactured-transport".into(),
                revision: 1,
            },
            claims,
            vec![],
        )
        .unwrap(),
        selection: MaterialPropertySelection::SingleClaimOnly,
    }
}
fn model(mu: f64, k: f64, b_maximum_temperature: f64) -> ThermochemicalGasModel {
    let (original, _, _) = fixture(0.0, 2000.0);
    ThermochemicalGasModel::new(
        original.thermodynamics().clone(),
        vec![
            card("A", mu, k, 2000.0),
            card("B", 3.0 * mu, 5.0 * k, b_maximum_temperature),
        ],
        None,
    )
    .unwrap()
}
fn context() -> QueryPoint {
    QueryPoint::new()
        .with("sample_id", 7.0)
        .unwrap()
        .with_quantity(
            "humidity",
            quantity(QuantityKind::Moisture(MoistureBasis::RelativeHumidity)),
            0.3,
        )
        .unwrap()
        .with_quantity(
            "temperature",
            quantity(QuantityKind::AbsoluteTemperature),
            999.0,
        )
        .unwrap()
        .with_quantity("pressure", quantity(QuantityKind::Pressure), 1000.0)
        .unwrap()
}
fn binary(material: &ThermochemicalGasModel) -> ConstantTransport {
    let models = material.thermodynamics().components().to_vec();
    let ids = models
        .iter()
        .map(|m| m.species().clone())
        .collect::<Vec<_>>();
    let source = IdealGasEquilibriumV1::new(
        ElementalMatrix::new(
            vec![ElementId::new("H").unwrap()],
            ids.clone(),
            vec![vec![1; ids.len()]],
        )
        .unwrap(),
        ChargeVector::new(ids.clone(), vec![0; ids.len()]).unwrap(),
        models,
    )
    .unwrap();
    let mut d = vec![1e-4; ids.len() * ids.len()];
    for i in 0..ids.len() {
        d[i * ids.len() + i] = 0.0;
    }
    ConstantTransport::new(
        source,
        "explicit-binary-diffusion",
        1,
        TransportDomain {
            minimum_temperature_k: 200.0,
            maximum_temperature_k: 2000.0,
            minimum_pressure_pa: 100.0,
            maximum_pressure_pa: 1e8,
        },
        // Intentionally unlike the material cards: none of these three scalar
        // outputs is owned by this pair-diffusion source in MaterialTransport.
        TransportCoefficients {
            dynamic_viscosity: DynViscosity::new(99.0),
            bulk_viscosity: DynViscosity::new(88.0),
            thermal_conductivity: ThermalConductivityQuantityV1::new(77.0),
        },
        d,
    )
    .unwrap()
}
fn adapter(material: ThermochemicalGasModel) -> MaterialTransport<ConstantTransport> {
    let diffusion = binary(&material);
    MaterialTransport::new(
        material,
        context(),
        diffusion,
        DynViscosity::new(0.0),
        || false,
    )
    .unwrap()
}
fn flow(
    material: &ThermochemicalGasModel,
    n: usize,
    initial: impl Fn(f64) -> (f64, f64),
) -> MixtureEuler1d {
    let g = MixtureGas::new(
        binary(material).source().clone(),
        CaloricDomain::new(Temperature::new(200.0), Temperature::new(1900.0)),
        || false,
    )
    .unwrap();
    let states = (0..n)
        .map(|i| {
            let (u, t) = initial((i as f64 + 0.5) / n as f64);
            g.conserved(&[1.0, 0.0], u, R / 0.01 * t).unwrap()
        })
        .collect();
    MixtureEuler1d::new(
        g,
        Length::new(1.0 / n as f64),
        0.7,
        states,
        Boundaries::Periodic,
        || false,
    )
    .unwrap()
}
fn cells(flow: &MixtureEuler1d) -> Vec<(Vec<f64>, f64, f64)> {
    (0..flow.cell_count())
        .map(|i| {
            let c = flow.cell(i).unwrap();
            (c.partial_densities.to_vec(), c.momentum, c.energy)
        })
        .collect()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * (1.0 + b.abs()),
        "{a:e} != {b:e}"
    );
}

#[test]
fn g0_source_cards_keep_query_context_zero_species_and_explicit_coefficient_ownership() {
    let material = model(0.01, 5.0, 300.0); // B is unsupported at the actual400K.
    let mut transport = adapter(material.clone());
    let resolved = transport
        .resolve(
            Temperature::new(400.0),
            Pressure::new(101325.0),
            &[1.0, 0.0],
            || false,
        )
        .unwrap();
    assert_eq!(resolved.components().len(), 1);
    let query = resolved.components()[0].query();
    assert_eq!(query.axes()["temperature"], 400.0);
    assert_eq!(query.axes()["pressure"], 101325.0);
    assert_eq!(query.axes()["sample_id"], 7.0);
    assert_eq!(query.axes()["humidity"], 0.3);
    assert_eq!(
        query.axis_quantities()["humidity"],
        context().axis_quantities()["humidity"]
    );
    assert_eq!(
        resolved.components()[0].card_identity(),
        material.transport_materials()[0].material.content_hash()
    );
    assert_eq!(transport.material(), &material);
    assert_eq!(transport.context(), &context());
    let mut d = [9.0; 4];
    let values = transport
        .evaluate(
            Temperature::new(400.0),
            Pressure::new(101325.0),
            &[1.0, 0.0],
            &mut d,
            &mut || false,
        )
        .unwrap();
    assert_eq!(values.dynamic_viscosity.value(), 0.01);
    assert_eq!(values.thermal_conductivity.value(), 5.0);
    assert_eq!(values.bulk_viscosity.value(), 0.0);
    assert_eq!(d, [0.0, 1e-4, 1e-4, 0.0]);
    let old = d;
    assert!(matches!(
        transport.evaluate(
            Temperature::new(400.0),
            Pressure::new(101325.0),
            &[0.5, 0.5],
            &mut d,
            &mut || false
        ),
        Err(DissipationError::Material(_))
    ));
    assert_eq!(d, old);

    let supported = model(0.01, 5.0, 2000.0);
    let mut transport = adapter(supported);
    let pure_b = transport
        .evaluate(
            Temperature::new(400.0),
            Pressure::new(101325.0),
            &[0.0, 1.0],
            &mut d,
            &mut || false,
        )
        .unwrap();
    close(pure_b.dynamic_viscosity.value(), 0.03, 1e-15);
    close(pure_b.thermal_conductivity.value(), 25.0, 1e-15);
}

#[test]
fn g1_g3_card_only_conductivity_and_viscosity_changes_reach_conserved_transport() {
    let wave = |x: f64| (2.0 * std::f64::consts::PI * x).sin();
    let mut heat_outputs = Vec::new();
    for k in [50.0, 100.0] {
        let material = model(0.01, k, 2000.0);
        let mut f = flow(&material, 24, |x| (0.0, 400.0 + 20.0 * wave(x)));
        let original = f.totals();
        let mut op = DissipationOperator::new(&f, adapter(material), || false).unwrap();
        let report = op
            .advance_for(
                &mut f,
                Time::new(0.02),
                DissipationControl::default(),
                || false,
            )
            .unwrap();
        let cv = 3.0 * R / 0.01;
        let amplitude = 20.0 * (-k / cv * (2.0 * std::f64::consts::PI).powi(2) * 0.02).exp();
        let p = f.gas().primitive(f.cell(5).unwrap(), || false).unwrap();
        let expected = 400.0 + amplitude * wave(5.5 / 24.0);
        close(p.temperature, expected, 2e-5);
        heat_outputs.push(p.temperature);
        assert_eq!(f.time().value(), 0.0);
        close(f.totals().energy, original.energy, 1e-12);
        close(report.absolute_energy_residual_j, 0.0, 1e-8);
        assert!(f.totals().species_masses[1].abs() < 1e-16);
    }
    assert!(heat_outputs[1] < heat_outputs[0] - 0.2);

    let mut shear_outputs = Vec::new();
    for mu in [0.1, 0.2] {
        let material = model(mu, 0.025, 2000.0);
        let mut f = flow(&material, 24, |x| (10.0 * wave(x), 400.0));
        let original = f.totals();
        let mut op = DissipationOperator::new(&f, adapter(material), || false).unwrap();
        op.advance_for(
            &mut f,
            Time::new(0.002),
            DissipationControl::default(),
            || false,
        )
        .unwrap();
        let p = f.gas().primitive(f.cell(5).unwrap(), || false).unwrap();
        let expected = 10.0
            * wave(5.5 / 24.0)
            * (-(4.0 / 3.0) * mu * (2.0 * std::f64::consts::PI).powi(2) * 0.002).exp();
        close(p.velocity, expected, 2e-4);
        shear_outputs.push(p.velocity);
        close(f.totals().energy, original.energy, 1e-12);
    }
    assert!(shear_outputs[1] < shear_outputs[0] - 0.08);
}

struct ObservedDiffusion {
    inner: ConstantTransport,
    calls: Rc<Cell<usize>>,
    refuse: bool,
}
impl TransportLaw for ObservedDiffusion {
    fn source(&self) -> &IdealGasEquilibriumV1 {
        self.inner.source()
    }
    fn evaluate(
        &mut self,
        t: Temperature,
        p: Pressure,
        x: &[f64],
        out: &mut [f64],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<TransportCoefficients, DissipationError> {
        self.calls.set(self.calls.get() + 1);
        let values = self.inner.evaluate(t, p, x, out, cancelled)?;
        if self.refuse {
            return Err(DissipationError::Invalid(
                "injected late binary-source refusal",
            ));
        }
        Ok(values)
    }
}

#[test]
fn g0_g4_invalid_input_does_not_query_sources_and_late_failures_preserve_outputs() {
    let material = model(0.01, 1.0, 2000.0);
    let calls = Rc::new(Cell::new(0));
    let diffusion = ObservedDiffusion {
        inner: binary(&material),
        calls: calls.clone(),
        refuse: false,
    };
    let mut transport = MaterialTransport::new(
        material.clone(),
        context(),
        diffusion,
        DynViscosity::new(0.0),
        || false,
    )
    .unwrap();
    let mut output = [7.0; 4];
    for fractions in [
        &[0.0, 0.0][..],
        &[0.8, 0.3][..],
        &[1.0][..],
        &[f64::NAN, 0.0][..],
    ] {
        assert!(matches!(
            transport.evaluate(
                Temperature::new(400.0),
                Pressure::new(101325.0),
                fractions,
                &mut output,
                &mut || false
            ),
            Err(DissipationError::Invalid(_))
        ));
    }
    assert_eq!(calls.get(), 0);
    assert_eq!(output, [7.0; 4]);
    let mut wrong_length = [7.0; 3];
    assert!(matches!(
        transport.evaluate(
            Temperature::new(400.0),
            Pressure::new(101325.0),
            &[1.0, 0.0],
            &mut wrong_length,
            &mut || false
        ),
        Err(DissipationError::Invalid(_))
    ));
    assert_eq!(calls.get(), 0);

    let mut polls = 0;
    transport
        .evaluate(
            Temperature::new(400.0),
            Pressure::new(101325.0),
            &[1.0, 0.0],
            &mut output,
            &mut || {
                polls += 1;
                false
            },
        )
        .unwrap();
    output.fill(7.0);
    let mut current = 0;
    assert_eq!(
        transport.evaluate(
            Temperature::new(400.0),
            Pressure::new(101325.0),
            &[1.0, 0.0],
            &mut output,
            &mut || {
                current += 1;
                current == polls
            }
        ),
        Err(DissipationError::Cancelled)
    );
    assert_eq!(output, [7.0; 4]);

    let diffusion = ObservedDiffusion {
        inner: binary(&material),
        calls: calls.clone(),
        refuse: true,
    };
    let mut transport = MaterialTransport::new(
        material.clone(),
        context(),
        diffusion,
        DynViscosity::new(0.0),
        || false,
    )
    .unwrap();
    assert!(matches!(
        transport.evaluate(
            Temperature::new(400.0),
            Pressure::new(101325.0),
            &[1.0, 0.0],
            &mut output,
            &mut || false
        ),
        Err(DissipationError::Invalid(_))
    ));
    assert_eq!(output, [7.0; 4]);
    let (wrong_material, _, _) = fixture(-1000.0, 2000.0);
    assert!(matches!(
        MaterialTransport::new(
            wrong_material,
            context(),
            binary(&material),
            DynViscosity::new(0.0),
            || false
        ),
        Err(DissipationError::Invalid(_))
    ));
    assert!(matches!(
        MaterialTransport::new(
            material.clone(),
            context(),
            binary(&material),
            DynViscosity::new(-1.0),
            || false
        ),
        Err(DissipationError::Invalid(_))
    ));
}

#[test]
fn g0_native_face_roundoff_is_closed_once_for_both_material_and_binary_owners() {
    struct RecordedDiffusion {
        law: ConstantTransport,
        fractions: Vec<f64>,
    }
    impl TransportLaw for RecordedDiffusion {
        fn source(&self) -> &IdealGasEquilibriumV1 {
            self.law.source()
        }
        fn evaluate(
            &mut self,
            t: Temperature,
            p: Pressure,
            x: &[f64],
            out: &mut [f64],
            cancelled: &mut dyn FnMut() -> bool,
        ) -> Result<TransportCoefficients, DissipationError> {
            self.fractions = x.to_vec();
            self.law.evaluate(t, p, x, out, cancelled)
        }
    }
    let material = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            ["A", "B", "C", "D"]
                .into_iter()
                .map(|name| species(name, 0.0))
                .collect(),
            Composition::new(CompositionBasis::MoleFraction, vec![0.25; 4]).unwrap(),
        )
        .unwrap(),
        vec![
            card("A", 0.01, 1.0, 2000.0),
            card("B", 0.02, 2.0, 2000.0),
            card("C", 0.03, 3.0, 2000.0),
            card("D", 0.04, 4.0, 300.0),
        ],
        None,
    )
    .unwrap();
    let law = RecordedDiffusion {
        law: binary(&material),
        fractions: Vec::new(),
    };
    let mut transport =
        MaterialTransport::new(material, context(), law, DynViscosity::new(0.0), || false).unwrap();
    // This is a valid native face closure: the largest component is exactly
    // one minus the other components, but canonical left-to-right addition is
    // one ulp below one. SpeciesD remains exactly absent and cannot be queried.
    let fractions: [f64; 4] = [
        0.038992402865305836,
        0.7185878433459016,
        0.24241975378879246,
        0.0,
    ];
    assert_eq!(fractions[1], 1.0 - (fractions[0] + fractions[2]));
    assert_eq!(fractions.iter().sum::<f64>(), 1.0_f64.next_down());
    let resolved = transport
        .resolve(
            Temperature::new(400.0),
            Pressure::new(101325.0),
            &fractions,
            || false,
        )
        .unwrap();
    let mut d = [9.0; 16];
    transport
        .evaluate(
            Temperature::new(400.0),
            Pressure::new(101325.0),
            &fractions,
            &mut d,
            &mut || false,
        )
        .unwrap();
    let canonical = &transport.diffusion_law().fractions;
    assert_eq!(canonical.iter().sum::<f64>(), 1.0);
    assert_eq!(canonical[3], 0.0);
    assert!(canonical[..3].iter().all(|x| *x > 0.0));
    assert_eq!(
        resolved.thermodynamics().receipt().mole_fraction_bits(),
        canonical[..3]
            .iter()
            .map(|x| x.to_bits())
            .collect::<Vec<_>>()
    );
    assert_eq!(resolved.components().len(), 3);
    for (&actual, &original) in canonical.iter().zip(&fractions) {
        assert!((actual - original).abs() <= 64.0 * f64::EPSILON * 4.0);
    }
    let previous = d;
    assert!(matches!(
        transport.evaluate(
            Temperature::new(400.0),
            Pressure::new(101325.0),
            &[0.04, 0.71, 0.24, 0.0],
            &mut d,
            &mut || false
        ),
        Err(DissipationError::Invalid(_))
    ));
    assert_eq!(d, previous);
}

#[test]
fn g4_new_reaction_product_requires_its_real_transport_card_before_atomic_publication() {
    let material = model(0.01, 1.0, 300.0); // Initially absentB has no400K support.
    let (_, mechanism, _) = fixture(0.0, 2000.0);
    let initial = flow(&material, 3, |_| (0.0, 400.0));
    let before = cells(&initial);
    let totals = initial.totals();
    let mut op = DissipationOperator::new(&initial, adapter(material), || false).unwrap();
    let mut reactive =
        ReactingEuler1d::new(initial, mechanism, controls().chemistry, || false).unwrap();
    let error = reactive
        .step_with_dissipation(
            &mut op,
            Time::new(1e-5),
            ReactiveStepControl {
                maximum_coupling_retries: 0,
                ..ReactiveStepControl::default()
            },
            DissipationControl::default(),
            || false,
        )
        .unwrap_err();
    fn material_error(error: &ReactiveError) -> bool {
        match error {
            ReactiveError::Dissipation(DissipationError::Material(_)) => true,
            ReactiveError::RetryLimit { last, .. } => material_error(last),
            _ => false,
        }
    }
    assert!(material_error(&error), "unexpected refusal: {error:?}");
    assert_eq!(cells(reactive.flow()), before);
    assert_eq!(reactive.flow().totals(), totals);
    assert_eq!(reactive.flow().time().value(), 0.0);
}
