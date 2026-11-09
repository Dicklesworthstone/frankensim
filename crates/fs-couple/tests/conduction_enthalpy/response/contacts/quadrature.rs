//! G1/G3/G4: nonuniform original interface response shares finite thermal history.
//! The first endpoint is an exact continuum affine solution, not a fitted network.
use super::*;
use fs_em::conduction::response::{
    ConductionResponseError,
    contact::{ContactResponseControls, ContactResponseError, ContactResponseQuadrature},
};
use fs_matdb::{SurfaceSpec, SystemContext};
use fs_material::electrical_interface::{
    ELECTRICAL_CONTACT_RESISTANCE_DIMS, ELECTRICAL_CONTACT_RESISTANCE_PROPERTY,
    response::{ElectricalInterfaceCoordinate, ScalarElectricalInterfaceResponse},
};

fn interface_with_spatial_response(
    cards: &[MaterialCard; 2],
    target_mean_k: f64,
    slope: f64,
) -> InterfaceSystemCard {
    let jump = ElectricalInterfaceCoordinate::VoltageJumpMagnitude;
    let absolute = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(
                ELECTRICAL_CONTACT_RESISTANCE_PROPERTY,
                ELECTRICAL_CONTACT_RESISTANCE_DIMS,
            ),
            value: PropertyValue::Surface2D {
                axes: ["Tm".into(), "jump".into()],
                axis_dims: [Temperature::DIMS, jump.quantity().dims()],
                samples: [[295.0, 0.5], [305.0, 0.5], [295.0, 3.0], [305.0, 3.0]]
                    .into_iter()
                    .map(|p| (p, p[1] + slope * (p[0] - target_mean_k)))
                    .collect(),
                triangles: vec![[0, 1, 2], [1, 3, 2]],
                dims: ELECTRICAL_CONTACT_RESISTANCE_DIMS,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity("Tm", absolute, 295.0, 305.0)
                .with_quantity("jump", jump.quantity(), 0.5, 3.0),
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured nonuniform contact with original thermal feedback".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    InterfaceSystemCard::assemble(
        SurfaceSpec {
            material: cards[0].id().clone(),
            texture_frame: "A".into(),
        },
        SurfaceSpec {
            material: cards[1].id().clone(),
            texture_frame: "B".into(),
        },
        SystemContext {
            medium: "manufactured".into(),
            third_body: None,
            environment: "fixed pressure".into(),
            history: "memoryless interface".into(),
        },
        claims,
        vec![],
    )
    .unwrap()
}

#[test]
fn nonuniform_interface_sources_heat_once_and_continue_with_original_thermal_context() {
    let coordinate = ConductorMagnitude::ElectricField;
    let (a, ca) = material(1.0, CP[0], 0.0, 0.0, coordinate, false);
    let (b, cb) = material(1.0, CP[1], 0.0, 0.0, coordinate, false);
    let (cards, curves) = ([a, b], [ca, cb]);
    let (dt, fraction, slope) = (0.1, 0.25, 0.005);
    // At the first endpoint V_A=-x+y+2, V_B=-x+1. Both bodies have
    // J_n=1 A/m2 across the contact and jump=1+y. Its actual R''=1+y,
    // so integrated contact heat is 1.5 W. Bulk heat is 2 W and 1 W.
    // Only exterior voltage/current traces are imposed; contact nodes are free.
    let sources_w = [2.0 + fraction * 1.5, 1.0 + (1.0 - fraction) * 1.5];
    let (a, d) = (CP[0] + dt, CP[1] + dt);
    let determinant = a * d - dt * dt;
    let theta = [
        dt * (sources_w[0] * d + dt * sources_w[1]) / determinant,
        dt * (a * sources_w[1] + dt * sources_w[0]) / determinant,
    ];
    let target_mean_k = 300.0 + (theta[0] + theta[1]) / 2.0;
    let card = interface_with_spatial_response(&cards, target_mean_k, slope);
    let jump = ElectricalInterfaceCoordinate::VoltageJumpMagnitude;
    let context = point(300.0, coordinate);
    let contact_context = interface_context(300.0)
        .with_quantity("jump", jump.quantity(), 1.5)
        .unwrap();
    let source = resolve_electrical_interface_state(
        &card,
        &contact_context,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let law =
        ScalarElectricalInterfaceResponse::new(&card, source.clone(), &[("jump", jump)]).unwrap();
    let (_, bulk) = problem(&cards, coordinate).into_frozen_parts();
    let strict = PlanarConductionResponse::new_with_contact_response(
        contact_mesh(),
        bulk.bulk().to_vec(),
        vec![ConductionContact {
            first_vertices: [1, 2],
            second_vertices: [4, 7],
            source: source.clone(),
        }],
        vec![law],
        ContactResponseControls {
            absolute_resistance_tolerance_ohm_m2: 1e-12,
            relative_resistance_tolerance: 1e-12,
            voltage_tolerance_v: 1e-10,
            current_density_tolerance_a_m2: 1e-10,
        },
        limits(),
        &mut || false,
    )
    .unwrap();
    let voltage = [(0, 2.0), (3, 3.0), (5, -1.0), (6, -1.0)];
    let current = [
        ConductionCurrentTrace {
            vertices: [0, 1],
            inward_current_density_a_m2: -1.0,
        },
        ConductionCurrentTrace {
            vertices: [2, 3],
            inward_current_density_a_m2: 1.0,
        },
    ];
    assert!(matches!(
        strict.solve_candidate(
            &voltage,
            &current,
            limits(),
            field_controls(),
            controls(),
            &mut || false,
        ),
        Err(ConductionResponseError::Contact(
            ContactResponseError::NonuniformTrace { .. }
        ))
    ));
    let quadrature = ContactResponseQuadrature {
        subdivisions: 4,
        maximum_points: 32,
    };
    assert!(
        strict
            .clone()
            .with_contact_quadrature(ContactResponseQuadrature {
                maximum_points: 1,
                ..quadrature
            })
            .is_err()
    );
    let mut original = strict.with_contact_quadrature(quadrature).unwrap();
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = cells(&cards, &curves, &context, &key);
    let layout = network(&cells);
    let bindings = [binding(&cards, &card, &contact_context, fraction)];
    let drive = ConductionDrive::CurrentDensity {
        voltage: &voltage,
        currents: &current,
    };
    let mut thermal = initial(&curves, 999.0, None);
    let mut solve = config();
    solve.maximum_circuit_trials = 4096;
    for target in [dt, 2.0 * dt] {
        let before = thermal.energy_j();
        let previous_contacts = original.contacts().to_vec();
        let mut polls = 0;
        let next = layout
            .step_response_contact_candidate(
                &thermal,
                target,
                &original,
                drive,
                &bindings,
                limits(),
                field_controls(),
                controls(),
                solve,
                &[0.0; 2],
                || {
                    polls += 1;
                    false
                },
            )
            .unwrap();
        let temperatures =
            [CellId(30), CellId(10)].map(|cell| next.thermal().temperature_k(cell).unwrap());
        let step = &next.report().electrical;
        let contact = &step.contacts[0].response;
        let samples = &contact.quadrature;
        assert!(samples.len() > 2 && samples.len() <= quadrature.maximum_points);
        assert_eq!(
            samples.len(),
            step.field.evaluated_contacts[0].evaluated.len()
        );
        assert_eq!(
            samples.len(),
            step.field.evaluated_contacts[0].applied.len()
        );
        near(
            samples.iter().map(|sample| sample.area_weight_m2).sum(),
            1.0,
        );
        near(
            samples.iter().map(|sample| sample.power_w).sum(),
            contact.power_w,
        );
        near(
            samples
                .iter()
                .map(|sample| sample.area_weight_m2 * sample.current_density_a_m2)
                .sum(),
            contact.current_a,
        );
        assert!(
            samples
                .iter()
                .map(|sample| sample.resistance_ohm_m2)
                .fold(0.0, f64::max)
                - samples
                    .iter()
                    .map(|sample| sample.resistance_ohm_m2)
                    .fold(f64::INFINITY, f64::min)
                > 0.5
        );
        for ((sample, evaluated), applied) in samples
            .iter()
            .zip(&step.field.evaluated_contacts[0].evaluated)
            .zip(&step.field.evaluated_contacts[0].applied)
        {
            assert_eq!(applied.resistance_ohm_m2(), sample.resistance_ohm_m2);
            assert!(sample.area_weight_m2 > 0.0);
            near(sample.basis.iter().sum(), 1.0);
            near(
                sample.voltage_jump_v,
                sample.current_density_a_m2 * sample.resistance_ohm_m2,
            );
            near(
                sample.power_w,
                sample.area_weight_m2 * sample.voltage_jump_v * sample.current_density_a_m2,
            );
            near(
                evaluated.resistance_ohm_m2(),
                sample.voltage_jump_v.abs()
                    + slope * ((temperatures[0] + temperatures[1]) / 2.0 - target_mean_k),
            );
            near(
                evaluated.resolved().state_coordinate("jump").unwrap(),
                sample.voltage_jump_v.abs(),
            );
            near(
                evaluated.resolved().state_coordinate("Ta").unwrap(),
                temperatures[0],
            );
            near(
                evaluated.resolved().state_coordinate("Tb").unwrap(),
                temperatures[1],
            );
            near(
                evaluated.resolved().state_coordinate("Tm").unwrap(),
                (temperatures[0] + temperatures[1]) / 2.0,
            );
            assert_eq!(
                evaluated.resolved().properties()[0]
                    .answer()
                    .receipt
                    .selected,
                source.resolved().properties()[0].answer().receipt.selected,
            );
            if target == dt {
                near(sample.voltage_jump_v, 1.0 + sample.basis[1]);
                near(sample.current_density_a_m2, 1.0);
                near(sample.resistance_ohm_m2, 1.0 + sample.basis[1]);
            }
        }
        assert_eq!(step.contacts[0].cells, [CellId(30), CellId(10)]);
        near(step.contacts[0].heat_fraction_to_first, fraction);
        near(
            step.heat_j,
            dt * (step.field.cell_power_w.iter().sum::<f64>() + contact.power_w),
        );
        near(step.heat_j, step.supplied_work_j);
        near(next.thermal().energy_j() - before, step.heat_j);
        near(next.report().deposited_joule_j, step.heat_j);
        near(next.report().total_energy_residual_j, 0.0);
        if target == dt {
            for i in 0..2 {
                near(temperatures[i], 300.0 + theta[i]);
            }
            for (actual, expected) in step
                .field
                .voltage_v
                .iter()
                .zip([2.0, 1.0, 2.0, 3.0, 0.0, -1.0, -1.0, 0.0])
            {
                near(*actual, expected);
            }
            near(contact.current_a, 1.0);
            near(contact.power_w, 1.5);
            near(step.field.cell_power_w[..2].iter().sum(), 2.0);
            near(step.field.cell_power_w[2..].iter().sum(), 1.0);
            near(step.heat_j, 4.5 * dt);
            let mut bounded = solve;
            bounded.maximum_circuit_trials = next.report().circuit_trials - 1;
            assert!(matches!(
                layout.step_response_contact_candidate(
                    &thermal,
                    target,
                    &original,
                    drive,
                    &bindings,
                    limits(),
                    field_controls(),
                    controls(),
                    bounded,
                    &[0.0; 2],
                    || false,
                ),
                Err(CircuitEnthalpyError::CircuitBudget)
            ));
            let mut stop = 0;
            assert!(matches!(
                layout.step_response_contact_candidate(
                    &thermal,
                    target,
                    &original,
                    drive,
                    &bindings,
                    limits(),
                    field_controls(),
                    controls(),
                    solve,
                    &[0.0; 2],
                    || {
                        stop += 1;
                        stop >= polls - 1
                    },
                ),
                Err(CircuitEnthalpyError::Cancelled)
            ));
            assert_eq!(thermal.energy_j().to_bits(), before.to_bits());
            assert_eq!(original.contacts(), previous_contacts);
        }
        (thermal, original, _) = next.into_parts();
    }
    near(source.resolved().state_coordinate("Tm").unwrap(), 300.0);
    near(source.resolved().state_coordinate("jump").unwrap(), 1.5);
}
