//! G1/G3/G4: source phi/E and phi/J surfaces close the original latent balance.
//! Explicit manufactured data; no empirical metal or higher-dimensional fit.
use super::*;
use fs_couple::circuit_enthalpy::conduction::response::ConductionResponseEnthalpyCandidate;
use fs_em::conduction::{
    ConductionCurrentTrace,
    response::{ConductionResponseControls, PlanarConductionResponse},
};
use fs_material::conductor::{
    ResolvedConductorState,
    response::{ConductorMagnitude, ScalarConductorResponse},
};

const AXES: [(&str, ConductorMagnitude); 2] = [
    ("field", ConductorMagnitude::ElectricField),
    ("current", ConductorMagnitude::CurrentDensity),
];
fn field_point(t: f64, fraction: f64) -> QueryPoint {
    point(t, fraction)
        .with_quantity("field", AXES[0].1.quantity(), 0.)
        .unwrap()
        .with_quantity("current", AXES[1].1.quantity(), 0.)
        .unwrap()
}
fn response_controls() -> ConductionResponseControls {
    ConductionResponseControls {
        maximum_iterations: 100,
        absolute_conductivity_tolerance_s_m: 1e-12,
        relative_conductivity_tolerance: 1e-12,
        current_density_tolerance_a_m2: 1e-10,
    }
}
fn source_card(slope: f64, upper_fraction: f64, coordinate: ConductorMagnitude) -> MaterialCard {
    let basis = ElectricalPropertyBasis::Resistivity;
    let axis = if coordinate == ConductorMagnitude::ElectricField {
        "field"
    } else {
        "current"
    };
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new("declared_phase_transport", basis.dimensions()),
            value: PropertyValue::Surface2D {
                axes: [FRACTION.into(), axis.into()],
                axis_dims: [Dims([0; 6]), coordinate.quantity().dims()],
                samples: [
                    [0., 0.],
                    [upper_fraction, 0.],
                    [0., 20.],
                    [upper_fraction, 20.],
                ]
                .into_iter()
                .map(|p| (p, 1. + slope * p[0] + 0.1 * p[1]))
                .collect(),
                triangles: vec![[0, 1, 2], [1, 3, 2]],
                dims: basis.dimensions(),
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity(
                    "T",
                    QuantitySpec::dimensional(Temperature::DIMS),
                    300.,
                    350.,
                )
                .with_quantity(
                    FRACTION,
                    QuantitySpec::dimensional(Dims([0; 6])),
                    0.,
                    upper_fraction,
                )
                .with_quantity("field", AXES[0].1.quantity(), 0., 20.)
                .with_quantity("current", AXES[1].1.quantity(), 0., 20.),
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured confined phase and field conductor".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured".into(),
            phase: "solid-liquid chart".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn nonlinear(card: &MaterialCard, context: &QueryPoint) -> PlanarConductionResponse {
    let basis = ElectricalPropertyBasis::Resistivity;
    let key = PropertyKey::new("declared_phase_transport", basis.dimensions());
    let material = resolve_material_state_point(
        card,
        context,
        &[basis.requirement_for(&key).unwrap()],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let original =
        ResolvedConductorState::try_from_material_with_property(material, basis, &key).unwrap();
    let source = ScalarConductorResponse::new_with_liquid_mass_fraction_axis(
        card, original, &AXES, FRACTION,
    )
    .unwrap();
    assert!(
        PlanarConduction::new(
            mesh(),
            vec![source.original().clone(); 4],
            limits(),
            &mut || false
        )
        .is_err()
    );
    PlanarConductionResponse::new(mesh(), vec![source; 4], limits(), &mut || false).unwrap()
}
fn check_endpoint(
    before: &PlanarConductionResponse,
    after: &ConductionResponseEnthalpyCandidate<'_>,
    h: f64,
    t: f64,
    fraction: f64,
    j: f64,
    voltage: f64,
) {
    near(after.thermal().energy_j(), h);
    near(after.thermal().temperature_k(CellId(7)).unwrap(), t);
    near(
        after.report().source.cells[0].phase.liquid_mass_fraction(),
        fraction,
    );
    let step = &after.report().electrical;
    near(step.field.voltage_v[1], voltage);
    near(step.field.voltage_v[2], voltage);
    near(step.heat_j, step.supplied_work_j);
    near(step.heat_j, (step.end_s - step.start_s) * voltage * j);
    near(after.report().deposited_joule_j, step.heat_j);
    assert!(after.report().total_energy_residual_j.abs() < config().total_energy_tolerance_j);
    let direction = core::f64::consts::FRAC_1_SQRT_2;
    for triangle in 0..4 {
        near(step.field.current_density_a_m2[triangle][0], -j * direction);
        near(step.field.current_density_a_m2[triangle][1], -j * direction);
        near(step.field.cell_power_w[triangle], voltage * j / 4.);
        for source in [
            &after.problem().conductors()[triangle],
            &step.field.evaluated[triangle].conductor,
        ] {
            assert_eq!(source.liquid_mass_fraction_axis(), Some(FRACTION));
            assert_eq!(
                source.source_key(),
                before.conductors()[triangle].source_key()
            );
            near(source.material().query().axes()["T"], t);
            near(source.material().query().axes()[FRACTION], fraction);
            near(source.material().query().axes()["field"], voltage.abs());
            near(source.material().query().axes()["current"], j.abs());
            assert_eq!(
                source.source_property().answer().receipt.selected,
                before.conductors()[triangle]
                    .source_property()
                    .answer()
                    .receipt
                    .selected
            );
        }
    }
}

#[test]
fn voltage_phase_field_source_melts_then_continues_into_liquid_sensible_heat() {
    let mut fractions = Vec::new();
    for slope in [1., 3.] {
        let card = source_card(slope, 1., ConductorMagnitude::ElectricField);
        let caloric = curve(&card, 1.);
        let thermal = initial(&caloric, 10.);
        let context = field_point(305., 0.);
        let old = nonlinear(&card, &context);
        let cells = [binding(&card, &caloric, &context)];
        let layout = network(&cells);
        let voltage = [(0, 0.), (3, 0.), (1, 4.), (2, 4.)];
        let next = layout
            .step_response_candidate(
                &thermal,
                5.,
                &old,
                ConductionDrive::Voltage(&voltage),
                limits(),
                controls(),
                response_controls(),
                config(),
                &[30.],
                || false,
            )
            .unwrap();
        // Unit geometry and 1 kg: (10+100*phi)*(1.4+s*phi)=5*4^2.
        let b: f64 = 140. + 10. * slope;
        let fraction = 132. / (b + (b * b + 4. * 100. * slope * 66.).sqrt());
        let h = 20. + 100. * fraction;
        let j = 4. / (1.4 + slope * fraction);
        check_endpoint(&old, &next, h, 310., fraction, j, 4.);
        near(next.report().electrical.heat_j, h - 10.);
        // Independently choose V for h=140, T=320, phi=1 at the next BE endpoint.
        let gain = 140. - h;
        let linear = 0.1 * gain / 5.;
        let v = 0.5 * (linear + (linear * linear + 4. * gain * (1. + slope) / 5.).sqrt());
        let voltage = [(0, 0.), (3, 0.), (1, v), (2, v)];
        let liquid = layout
            .step_response_candidate(
                next.thermal(),
                10.,
                next.problem(),
                ConductionDrive::Voltage(&voltage),
                limits(),
                controls(),
                response_controls(),
                config(),
                &[gain],
                || false,
            )
            .unwrap();
        check_endpoint(
            next.problem(),
            &liquid,
            140.,
            320.,
            1.,
            v / (1. + slope + 0.1 * v),
            v,
        );
        near(liquid.report().electrical.heat_j, gain);
        near(thermal.energy_j(), 10.);
        assert_eq!(thermal.time_s(), 0.);
        assert_eq!(context.axes()[FRACTION], 0.);
        assert_eq!(old.conductors()[0].material().query().axes()["field"], 0.);
        fractions.push(fraction);
    }
    assert!(fractions[1] < fractions[0]);
}

#[test]
fn prescribed_current_density_uses_phase_current_source_and_latent_history() {
    let mut fractions = Vec::new();
    for slope in [1., 3.] {
        let card = source_card(slope, 1., ConductorMagnitude::CurrentDensity);
        let caloric = curve(&card, 1.);
        let thermal = initial(&caloric, 10.);
        let context = field_point(305., 0.);
        let old = nonlinear(&card, &context);
        let cells = [binding(&card, &caloric, &context)];
        let layout = network(&cells);
        let traces = [ConductionCurrentTrace {
            vertices: [1, 2],
            inward_current_density_a_m2: 2.,
        }];
        let drive = ConductionDrive::CurrentDensity {
            voltage: &[(0, 0.), (3, 0.)],
            currents: &traces,
        };
        let next = layout
            .step_response_candidate(
                &thermal,
                5.,
                &old,
                drive,
                limits(),
                controls(),
                response_controls(),
                config(),
                &[30.],
                || false,
            )
            .unwrap();
        // h=10+5*2^2*(1+s*phi+0.1*2)=20+100*phi on the plateau.
        let fraction = 14. / (100. - 20. * slope);
        let h = 20. + 100. * fraction;
        check_endpoint(
            &old,
            &next,
            h,
            310.,
            fraction,
            2.,
            2. * (1.2 + slope * fraction),
        );
        near(next.report().electrical.heat_j, h - 10.);
        let trial_fraction = (h - 20. + 24.) / (100. - 20. * slope);
        let (h2, t2, fraction2) = if trial_fraction < 1. {
            (20. + 100. * trial_fraction, 310., trial_fraction)
        } else {
            let h2 = h + 20. * (1.2 + slope);
            (h2, 310. + (h2 - 120.) / 2., 1.)
        };
        let continued = layout
            .step_response_candidate(
                next.thermal(),
                10.,
                next.problem(),
                drive,
                limits(),
                controls(),
                response_controls(),
                config(),
                &[h2 - h],
                || false,
            )
            .unwrap();
        check_endpoint(
            next.problem(),
            &continued,
            h2,
            t2,
            fraction2,
            2.,
            2. * (1.2 + slope * fraction2),
        );
        near(continued.report().electrical.heat_j, h2 - h);
        near(thermal.energy_j(), 10.);
        fractions.push(fraction);
    }
    assert!(fractions[1] > fractions[0]);
}

#[test]
fn phase_field_admission_domains_and_late_cancel_publish_neither_history() {
    let card = source_card(3., 1., ConductorMagnitude::ElectricField);
    let caloric = curve(&card, 1.);
    let thermal = initial(&caloric, 10.);
    let context = field_point(305., 0.);
    let old = nonlinear(&card, &context);
    let sources = old.conductors().to_vec();
    let cells = [binding(&card, &caloric, &context)];
    let layout = network(&cells);
    let voltage = [(0, 0.), (3, 0.), (1, 4.), (2, 4.)];
    let drive = ConductionDrive::Voltage(&voltage);
    let mut polls = 0;
    let good = layout
        .step_response_candidate(
            &thermal,
            5.,
            &old,
            drive,
            limits(),
            controls(),
            response_controls(),
            config(),
            &[30.],
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut calls = 0;
    assert!(matches!(
        layout.step_response_candidate(
            &thermal,
            5.,
            &old,
            drive,
            limits(),
            controls(),
            response_controls(),
            config(),
            &[30.],
            || {
                calls += 1;
                calls >= polls
            }
        ),
        Err(CircuitEnthalpyError::Cancelled)
    ));
    for axis in [None, Some("another_phase")] {
        let mut cell = binding(&card, &caloric, &context);
        cell.point.liquid_fraction_axis = axis;
        assert!(
            network(&[cell])
                .step_response_candidate(
                    &thermal,
                    5.,
                    &old,
                    drive,
                    limits(),
                    controls(),
                    response_controls(),
                    config(),
                    &[30.],
                    || false
                )
                .is_err()
        );
    }
    assert!(
        layout
            .step_response_candidate(
                good.thermal(),
                10.,
                &old,
                drive,
                limits(),
                controls(),
                response_controls(),
                config(),
                &[30.],
                || false
            )
            .is_err()
    );
    let limited = source_card(3., 0.2, ConductorMagnitude::ElectricField);
    let limited_curve = curve(&limited, 1.);
    let limited_thermal = initial(&limited_curve, 10.);
    let limited_old = nonlinear(&limited, &context);
    let limited_cell = binding(&limited, &limited_curve, &context);
    assert!(
        network(&[limited_cell])
            .step_response_candidate(
                &limited_thermal,
                5.,
                &limited_old,
                drive,
                limits(),
                controls(),
                response_controls(),
                config(),
                &[10.],
                || false
            )
            .is_err()
    );
    let retry = layout
        .step_response_candidate(
            &thermal,
            5.,
            &old,
            drive,
            limits(),
            controls(),
            response_controls(),
            config(),
            &[30.],
            || false,
        )
        .unwrap();
    near(retry.thermal().energy_j(), good.thermal().energy_j());
    assert_eq!(retry.problem().conductors(), good.problem().conductors());
    assert_eq!(old.conductors(), sources.as_slice());
    assert_eq!(thermal.time_s(), 0.);
    near(thermal.energy_j(), 10.);
    near(limited_thermal.energy_j(), 10.);
}

#[test]
fn total_current_electrodes_close_field_phase_and_latent_heat_for_both_signs() {
    let mut fractions = Vec::new();
    for slope in [1., 3.] {
        let card = source_card(slope, 1., ConductorMagnitude::ElectricField);
        let caloric = curve(&card, 1.);
        let thermal = initial(&caloric, 10.);
        let context = field_point(305., 0.);
        let old = nonlinear(&card, &context);
        let cells = [binding(&card, &caloric, &context)];
        let layout = network(&cells);
        for current in [2.0_f64, -2.0] {
            let drive = ConductionDrive::Current {
                positive: &[1, 2],
                negative: &[0, 3],
                current_a: current,
            };
            // rho=1+s*phi+0.1*|E| with |J|=2 implies
            // rho=(1+s*phi)/0.8. Original latent balance is
            // 20+100*phi=10+5*2²*(1+s*phi)/0.8.
            let fraction = 15. / (100. - 25. * slope);
            let h = 20. + 100. * fraction;
            let voltage = current * (1. + slope * fraction) / 0.8;
            let next = layout
                .step_response_candidate(
                    &thermal,
                    5.,
                    &old,
                    drive,
                    limits(),
                    controls(),
                    response_controls(),
                    config(),
                    &[h - 10.],
                    || false,
                )
                .unwrap();
            check_endpoint(&old, &next, h, 310., fraction, current, voltage);
            near(next.report().electrical.heat_j, h - 10.);
            let field = &next.report().electrical.field;
            near(
                [1, 2]
                    .into_iter()
                    .map(|v| field.nodal_boundary_current_a[v].unwrap())
                    .sum(),
                current,
            );
            near(
                [0, 3]
                    .into_iter()
                    .map(|v| field.nodal_boundary_current_a[v].unwrap())
                    .sum(),
                -current,
            );
            assert!(field.nodal_applied_current_a.iter().all(Option::is_none));
            near(
                next.report().electrical.supplied_work_j,
                5. * current * voltage,
            );
            // Continue the same original h history. The high-resistivity source
            // reaches h=180 (T=340); the other remains on the latent plateau.
            let trial_fraction = (h - 20. + 25.) / (100. - 25. * slope);
            let (h2, t2, fraction2) = if trial_fraction < 1. {
                (20. + 100. * trial_fraction, 310., trial_fraction)
            } else {
                let h2 = h + 25. * (1. + slope);
                (h2, 310. + (h2 - 120.) / 2., 1.)
            };
            let voltage2 = current * (1. + slope * fraction2) / 0.8;
            let continued = layout
                .step_response_candidate(
                    next.thermal(),
                    10.,
                    next.problem(),
                    drive,
                    limits(),
                    controls(),
                    response_controls(),
                    config(),
                    &[h2 - h],
                    || false,
                )
                .unwrap();
            check_endpoint(
                next.problem(),
                &continued,
                h2,
                t2,
                fraction2,
                current,
                voltage2,
            );
            near(continued.report().electrical.heat_j, h2 - h);
            near(
                continued.report().electrical.supplied_work_j,
                5. * current * voltage2,
            );
            if current > 0. {
                fractions.push(fraction);
            }
        }
        near(thermal.energy_j(), 10.);
        assert_eq!(thermal.time_s(), 0.);
        assert_eq!(context.axes()[FRACTION], 0.);
        assert_eq!(old.conductors()[0].material().query().axes()["field"], 0.);
    }
    assert!(fractions[1] > fractions[0]);
}

#[test]
fn total_current_nested_budget_refuses_atomically_and_zero_current_preserves_phase() {
    let card = source_card(3., 1., ConductorMagnitude::ElectricField);
    let caloric = curve(&card, 1.);
    let thermal = initial(&caloric, 10.);
    let context = field_point(305., 0.);
    let old = nonlinear(&card, &context);
    let sources = old.conductors().to_vec();
    let cells = [binding(&card, &caloric, &context)];
    let layout = network(&cells);
    let mut short = config();
    short.maximum_circuit_trials = 1;
    assert!(matches!(
        layout.step_response_candidate(
            &thermal,
            5.,
            &old,
            ConductionDrive::Current {
                positive: &[1, 2],
                negative: &[0, 3],
                current_a: 2.
            },
            limits(),
            controls(),
            response_controls(),
            short,
            &[0.],
            || false
        ),
        Err(CircuitEnthalpyError::CircuitBudget)
    ));
    let zero = layout
        .step_response_candidate(
            &thermal,
            5.,
            &old,
            ConductionDrive::Current {
                positive: &[1, 2],
                negative: &[0, 3],
                current_a: 0.,
            },
            limits(),
            controls(),
            response_controls(),
            config(),
            &[0.],
            || false,
        )
        .unwrap();
    check_endpoint(&old, &zero, 10., 305., 0., 0., 0.);
    assert_eq!(zero.thermal().time_s(), 5.);
    assert_eq!(old.conductors(), sources.as_slice());
    assert_eq!(thermal.time_s(), 0.);
    near(thermal.energy_j(), 10.);
}
