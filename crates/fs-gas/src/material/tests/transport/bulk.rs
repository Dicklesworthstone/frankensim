//! G0/G1/G3/G4: live mixture bulk loss on the original conservative face owner.
use super::*;
use crate::material::transport::{BulkViscosityMaterial, MaterialBulkViscosityTransport};
use crate::multicomponent::NormalPrimitive;
use crate::multicomponent::dissipation::{MaxwellStefanSolveV1, PlanarTransport};
use fs_qty::Dims;
use fs_qty::semantic::{CompositionBasis as CoordinateBasis, FrequencyConvention};

fn axes() -> Vec<(SpeciesId, String)> {
    vec![(id("B"), "x_B".into()), (id("A"), "x_A".into())]
}
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
fn expected(t: f64, p: f64, x: f64, pressure_axis: bool) -> f64 {
    0.02 + 0.04 * x + if pressure_axis { 1e-7 * p } else { 1e-4 * t }
}
fn claim(pressure_axis: bool) -> PropertyClaim {
    let (axis, dims, bounds) = if pressure_axis {
        ("pressure", Pressure::DIMS, [50000.0, 800000.0])
    } else {
        ("temperature", Temperature::DIMS, [200.0, 800.0])
    };
    PropertyClaim {
        key: PropertyKey::new("mixture_bulk_viscosity", DynViscosity::DIMS),
        value: PropertyValue::Surface2D {
            axes: [axis.into(), "x_A".into()],
            axis_dims: [dims, Dims::NONE],
            samples: [
                (bounds[0], 0.0),
                (bounds[1], 0.0),
                (bounds[0], 1.0),
                (bounds[1], 1.0),
            ]
            .map(|(a, x)| ([a, x], expected(a, a, x, pressure_axis)))
            .to_vec(),
            triangles: vec![[0, 1, 2], [1, 3, 2]],
            dims: DynViscosity::DIMS,
        },
        validity: ValidityDomain::unconstrained()
            .with_quantity(
                "temperature",
                quantity(QuantityKind::AbsoluteTemperature),
                200.0,
                800.0,
            )
            .with_quantity(
                "pressure",
                quantity(QuantityKind::Pressure),
                50000.0,
                800000.0,
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
fn source(claims: Vec<PropertyClaim>) -> BulkViscosityMaterial {
    let mut set = ClaimSet::new();
    for claim in claims {
        set.insert_claim(claim).unwrap();
    }
    BulkViscosityMaterial {
        material: MaterialCard::assemble(
            MaterialStateId {
                chemistry: "A+B".into(),
                phase: "gas".into(),
                process: "manufactured-mixture-bulk-law".into(),
                revision: 0,
            },
            set,
            vec![],
        )
        .unwrap(),
        property: PropertyKey::new("mixture_bulk_viscosity", DynViscosity::DIMS),
        selection: MaterialPropertySelection::SingleClaimOnly,
    }
}
fn sourced(
    pressure_axis: bool,
) -> MaterialBulkViscosityTransport<MaterialTransport<ConstantTransport>> {
    MaterialBulkViscosityTransport::new(
        adapter(model(0.01, 0.03, 2000.0)),
        source(vec![claim(pressure_axis)]),
        axes(),
        reference(),
        || false,
    )
    .unwrap()
}

#[test]
fn g1_bulk_material_drives_actual_compression_stress_and_conservative_work() {
    for pressure_axis in [false, true] {
        let law = sourced(pressure_axis);
        let original = law.last_resolved().material().clone();
        let gas = MixtureGas::new(
            law.source().clone(),
            CaloricDomain::new(Temperature::new(200.0), Temperature::new(1900.0)),
            || false,
        )
        .unwrap();
        let mut face = PlanarTransport::new(gas, law, || false).unwrap();
        for (t, p, x) in [
            (400.0, 100000.0, 0.4),
            (600.0, 100000.0, 0.4),
            (400.0, 170000.0, 0.4),
            (400.0, 100000.0, 0.7),
            (400.0, 100000.0, 1.0),
            (400.0, 100000.0, 0.0),
        ] {
            let rho = p * 0.01 / (R * t);
            let densities = [rho * x, rho * (1.0 - x)];
            let primitive = |u, v| NormalPrimitive {
                partial_densities: &densities,
                normal_velocity: u,
                transverse_velocity: v,
                pressure: p,
            };
            for transverse_divergence in [0.0, 10.0] {
                let mut flux = [0.0; 5];
                let report = face
                    .face_from_primitives(
                        primitive(0.8, 0.1),
                        primitive(-0.2, 0.4),
                        0.1,
                        [0.0, transverse_divergence],
                        MaxwellStefanSolveV1::default(),
                        &mut flux,
                        || false,
                    )
                    .unwrap();
                let mu = report.dynamic_viscosity_pa_s;
                let zeta = expected(t, p, x, pressure_axis);
                let normal_stress = -(4.0 / 3.0 * mu + zeta) * 10.0
                    + (zeta - 2.0 / 3.0 * mu) * transverse_divergence;
                let shear_stress = mu * 3.0;
                close(flux[0], 0.0, 1e-14);
                close(flux[1], 0.0, 1e-14);
                close(flux[2], -normal_stress, 1e-12);
                close(flux[3], -shear_stress, 1e-12);
                close(flux[4], -0.3 * normal_stress - 0.25 * shear_stress, 1e-12);
                if transverse_divergence == 10.0 {
                    // Trace-free compression/shear cannot acquire bulk stress.
                    close(flux[2], 20.0 * mu, 1e-12);
                }
                let resolved = face.law().last_resolved();
                close(resolved.bulk_viscosity().value(), zeta, 1e-14);
                close(resolved.temperature().value(), t, 1e-12);
                close(resolved.pressure().value(), p, 1e-12);
                assert_eq!(
                    resolved.material().card_identity(),
                    original.card_identity()
                );
                assert_eq!(
                    resolved.material().properties()[0]
                        .answer()
                        .receipt
                        .selected,
                    original.properties()[0].answer().receipt.selected
                );
                assert_eq!(resolved.material().query().axes()["humidity"], 0.3);
            }
        }
    }
}

#[test]
fn g3_bulk_only_change_dissipates_compression_into_heat_without_changing_total_energy() {
    let mut final_kinetic = Vec::new();
    for zeta in [0.0, 0.5] {
        let material = model(0.01, 0.03, 2000.0);
        let mut flow = flow(&material, 24, |x| {
            (10.0 * (2.0 * std::f64::consts::PI * x).sin(), 400.0)
        });
        let mut c = claim(false);
        c.value = PropertyValue::Scalar {
            value: zeta,
            dims: DynViscosity::DIMS,
        };
        c.interpolation = InterpolationPolicy::ConstantWithinValidity;
        let law = MaterialBulkViscosityTransport::new(
            adapter(material),
            source(vec![c]),
            axes(),
            reference(),
            || false,
        )
        .unwrap();
        let before = flow.totals();
        let mut op = DissipationOperator::new(&flow, law, || false).unwrap();
        let report = op
            .advance_for(
                &mut flow,
                Time::new(0.002),
                DissipationControl::default(),
                || false,
            )
            .unwrap();
        let kinetic = (0..flow.cell_count())
            .map(|i| {
                let q = flow.cell(i).unwrap();
                0.5 * q.momentum * q.momentum / q.partial_densities.iter().sum::<f64>()
                    * flow.mesh().volumes()[i]
            })
            .sum::<f64>();
        final_kinetic.push(kinetic);
        close(flow.totals().energy, before.energy, 1e-12);
        close(flow.totals().mass(), before.mass(), 1e-14);
        close(report.absolute_energy_residual_j, 0.0, 1e-7);
    }
    assert!(
        final_kinetic[1] < final_kinetic[0] - 1.0,
        "bulk loss must convert additional resolved kinetic energy to thermal storage"
    );
}

#[test]
fn g0_bulk_requires_explicit_units_complete_typed_axes_and_instantaneous_support() {
    let construct = |binding, coordinates, query| {
        MaterialBulkViscosityTransport::new(
            adapter(model(0.01, 0.03, 2000.0)),
            binding,
            coordinates,
            query,
            || false,
        )
    };
    let mut wrong_units = source(vec![claim(false)]);
    wrong_units.property = PropertyKey::new("mixture_bulk_viscosity", Pressure::DIMS);
    assert!(construct(wrong_units, axes(), reference()).is_err());
    assert!(
        construct(
            source(vec![claim(false)]),
            vec![(id("A"), "x_A".into())],
            reference()
        )
        .is_err()
    );
    assert!(
        construct(
            source(vec![claim(false)]),
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
            source(vec![claim(false)]),
            axes(),
            reference()
                .with_quantity("stale_pressure", quantity(QuantityKind::Pressure), 1e5)
                .unwrap()
        )
        .is_err()
    );
    let mut negative = claim(false);
    negative.value = PropertyValue::Scalar {
        value: -0.01,
        dims: DynViscosity::DIMS,
    };
    negative.interpolation = InterpolationPolicy::ConstantWithinValidity;
    assert!(construct(source(vec![negative]), axes(), reference()).is_err());
    let frequency = quantity(QuantityKind::Frequency(FrequencyConvention::Cyclic));
    let mut dynamic = claim(false);
    dynamic.validity = dynamic
        .validity
        .with_quantity("frequency", frequency, 0.0, 1000.0);
    assert!(
        construct(
            source(vec![dynamic]),
            axes(),
            reference()
                .with_quantity("frequency", frequency, 10.0)
                .unwrap()
        )
        .is_err()
    );
}

#[test]
fn g0_g4_source_refusal_and_cancellation_preserve_outputs_and_bulk_receipt() {
    let mut law = sourced(false);
    let original = law.last_resolved().clone();
    let mut output = [7.0; 4];
    for (t, p, x) in [
        (900.0, 1e5, [0.4, 0.6]),
        (400.0, 1e6, [0.4, 0.6]),
        (400.0, 1e5, [0.4, 0.7]),
        (f64::NAN, 1e5, [0.4, 0.6]),
    ] {
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
        assert_eq!(output, [7.0; 4]);
        assert_eq!(law.last_resolved(), &original);
    }
    let mut short = [7.0; 3];
    assert!(
        law.evaluate(
            Temperature::new(500.0),
            Pressure::new(1e5),
            &[0.4, 0.6],
            &mut short,
            &mut || false
        )
        .is_err()
    );
    assert_eq!(short, [7.0; 3]);
    let mut polls = 0;
    law.evaluate(
        Temperature::new(500.0),
        Pressure::new(1e5),
        &[0.4, 0.6],
        &mut output,
        &mut || {
            polls += 1;
            false
        },
    )
    .unwrap();
    for cancel_at in 1..=polls {
        let mut law = sourced(false);
        let mut output = [7.0; 4];
        let mut count = 0;
        assert!(matches!(
            law.evaluate(
                Temperature::new(500.0),
                Pressure::new(1e5),
                &[0.4, 0.6],
                &mut output,
                &mut || {
                    count += 1;
                    count == cancel_at
                }
            ),
            Err(DissipationError::Cancelled)
        ));
        assert_eq!(output, [7.0; 4]);
        assert_eq!(law.last_resolved(), &original);
    }
}

struct CrossLaw {
    inner: ConstantTransport,
    failure: u8,
}
impl TransportLaw for CrossLaw {
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
        let result = self.inner.evaluate(t, p, x, out, cancelled)?;
        if self.failure == 1 {
            return Err(DissipationError::Invalid("late original refusal"));
        }
        if self.failure == 2 {
            out[0] = 1.0;
        }
        Ok(result)
    }
    fn heats_of_transport(
        &mut self,
        t: Temperature,
        p: Pressure,
        x: &[f64],
        out: &mut [f64],
        _: &mut dyn FnMut() -> bool,
    ) -> Result<bool, DissipationError> {
        out.copy_from_slice(&[t.value() + 100.0 * x[0], p.value() * x[1]]);
        if self.failure == 1 {
            return Err(DissipationError::Invalid("late heat refusal"));
        }
        if self.failure == 2 {
            out[0] = f64::NAN;
        }
        Ok(true)
    }
}

#[test]
fn g0_g4_bulk_delegates_all_other_coefficients_and_heats_atomically() {
    for failure in 0..3 {
        let mut law = MaterialBulkViscosityTransport::new(
            CrossLaw {
                inner: binary(&model(0.01, 0.03, 2000.0)),
                failure,
            },
            source(vec![claim(false)]),
            axes(),
            reference(),
            || false,
        )
        .unwrap();
        let original = law.last_resolved().clone();
        let mut out = [7.0; 4];
        let result = law.evaluate(
            Temperature::new(500.0),
            Pressure::new(120000.0),
            &[0.7, 0.3],
            &mut out,
            &mut || false,
        );
        let mut heats = [9.0; 2];
        let heat_result = law.heats_of_transport(
            Temperature::new(500.0),
            Pressure::new(120000.0),
            &[0.7, 0.3],
            &mut heats,
            &mut || false,
        );
        if failure == 0 {
            let c = result.unwrap();
            assert_eq!(c.dynamic_viscosity.value(), 99.0);
            assert_eq!(c.thermal_conductivity.value(), 77.0);
            close(
                c.bulk_viscosity.value(),
                expected(500.0, 120000.0, 0.7, false),
                1e-15,
            );
            assert_eq!(out, [0.0, 1e-4, 1e-4, 0.0]);
            assert!(heat_result.unwrap());
            close(heats[0], 570.0, 1e-14);
            close(heats[1], 36000.0, 1e-14);
        } else {
            assert!(result.is_err() && heat_result.is_err());
            assert_eq!(out, [7.0; 4]);
            assert_eq!(heats, [9.0; 2]);
            assert_eq!(law.last_resolved(), &original);
        }
    }
}
