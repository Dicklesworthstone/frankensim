//! Fixed-temperature field snapshots use the same scaled displacement and
//! potential as native capacitors; these solves do not advance thermal history.
use super::*;
use fs_material::field_storage::dielectric::reversible::thermal::ThermalReversibleDielectricMaterial;
use fs_qty::{
    Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

fn source() -> ThermalReversibleDielectricMaterial {
    let reference = nonlinear();
    let mut claims = reference.card().claims().clone();
    let temperature = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    let scale = PropertyKey::new("dielectric_temperature_scale", Dims([0; 6]));
    claims
        .insert_claim(PropertyClaim {
            key: scale.clone(),
            value: PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(280.0, 0.96), (300.0, 1.0), (340.0, 1.08)],
                dims: Dims([0; 6]),
            },
            validity: ValidityDomain::unconstrained().with_quantity("T", temperature, 280.0, 340.0),
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "authored nonlinear electrostatic thermal scale".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic scaled dielectric".into(),
            phase: "solid".into(),
            process: "explicit separable potential".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let point = QueryPoint::new()
        .with_quantity("E", QuantitySpec::dimensional(ELECTRIC_FIELD_DIMS), 0.0)
        .unwrap()
        .with_quantity("T", temperature, 300.0)
        .unwrap();
    let reference = ReversibleDielectricMaterial::resolve(
        &card,
        &point,
        "E",
        4.0 * E0,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let material = resolve_material_state_point(
        &card,
        &point,
        &[
            ScalarPropertyRequirement::try_new(
                REVERSIBLE_ELECTRIC_DISPLACEMENT_PROPERTY,
                ELECTRIC_DISPLACEMENT_DIMS,
                ScalarAdmissibility::NonNegative,
            )
            .unwrap(),
            ScalarPropertyRequirement::try_new(
                scale.name(),
                Dims([0; 6]),
                ScalarAdmissibility::StrictlyPositive,
            )
            .unwrap(),
        ],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    ThermalReversibleDielectricMaterial::new(&reference, &material, &scale, "T", [280.0, 340.0])
        .unwrap()
}

#[test]
fn actual_temperature_changes_picofarad_charge_tangent_and_fixed_charge_voltage() {
    let law = source();
    for temperature in [280.0, 300.0, 340.0] {
        let scale = 1.0 + 0.002 * (temperature - 300.0);
        let snapshot = law.at_temperature(temperature).unwrap();
        let problem = PlanarNonlinearElectrostatics::new(
            slab_mesh(),
            vec![DielectricMaterial::Reversible(snapshot.dielectric().clone()); 4],
            limits(),
            &mut || false,
        )
        .unwrap();
        let electrodes = problem
            .prepare_electrodes(&[&[0, 3], &[2, 5]], &mut || false)
            .unwrap();
        let fixed_voltage = electrodes
            .solve(
                &[0.0; 4],
                &[
                    ElectrodeControl::Voltage(0.0),
                    ElectrodeControl::Voltage(2.0 * V0),
                ],
                &[0.0; 6],
                controls(),
                &mut || false,
            )
            .unwrap();
        near(fixed_voltage.terminal_charge_c()[1], scale * 2.5 * Q0, Q0);
        near(fixed_voltage.energy_j(), scale * 1.75 * U0, U0);
        near(fixed_voltage.coenergy_j(), scale * 3.25 * U0, U0);
        let tangent = fixed_voltage
            .differential_capacitance(differential_controls(), &mut || false)
            .unwrap();
        let c = scale * 0.5 * EPSILON0 * AREA / GAP;
        for (actual, expected) in tangent.matrix_f().iter().zip([c, -c, -c, c]) {
            near(*actual, expected, 1e-12);
        }

        let fixed_charge = electrodes
            .solve(
                &[0.0; 4],
                &[
                    ElectrodeControl::Voltage(0.0),
                    ElectrodeControl::Charge(2.5 * Q0),
                ],
                fixed_voltage.relative_voltage_v(),
                controls(),
                &mut || false,
            )
            .unwrap();
        let z = 2.5 / scale;
        let electric = 2.0 * z - 3.0;
        let free_energy = z * z - 3.0 * z + 3.0;
        near(fixed_charge.terminal_voltage_v()[1], electric * V0, V0);
        near(fixed_charge.terminal_charge_c()[1], 2.5 * Q0, Q0);
        near(fixed_charge.energy_j(), scale * free_energy * U0, U0);
        near(
            fixed_charge.coenergy_j(),
            scale * (z * electric - free_energy) * U0,
            U0,
        );
        for material in fixed_charge.problem().materials() {
            assert_eq!(material.material().query().axes()["T"], temperature);
            assert_eq!(
                material.material().card_identity(),
                law.card().content_hash()
            );
        }
    }
}
