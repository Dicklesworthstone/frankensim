//! Original nonlinear bulk and surface conduction share one native field.
use super::super::contacts::{binding, interface, mesh};
use super::*;
use fs_em::conduction::ConductionContact;
use fs_matdb::{InterfaceSystemCard, SurfaceSpec, SystemContext};
use fs_material::electrical_interface::{
    ELECTRICAL_CONTACT_RESISTANCE_DIMS, ELECTRICAL_CONTACT_RESISTANCE_PROPERTY,
    resolve_electrical_interface_state,
};

fn contact(
    a: &ScalarConductorResponse,
    b: &ScalarConductorResponse,
    resistance: f64,
) -> ConductionContact {
    contact_at(a, b, resistance, &QueryPoint::new())
}

fn contact_at(
    a: &ScalarConductorResponse,
    b: &ScalarConductorResponse,
    resistance: f64,
    point: &QueryPoint,
) -> ConductionContact {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(
                ELECTRICAL_CONTACT_RESISTANCE_PROPERTY,
                ELECTRICAL_CONTACT_RESISTANCE_DIMS,
            ),
            value: PropertyValue::Scalar {
                value: resistance,
                dims: ELECTRICAL_CONTACT_RESISTANCE_DIMS,
            },
            validity: ValidityDomain::unconstrained(),
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured constant contact".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = InterfaceSystemCard::assemble(
        SurfaceSpec {
            material: a.original().material().material().clone(),
            texture_frame: "left".into(),
        },
        SurfaceSpec {
            material: b.original().material().material().clone(),
            texture_frame: "right".into(),
        },
        SystemContext {
            medium: "manufactured".into(),
            third_body: None,
            environment: "fixed".into(),
            history: "fixed".into(),
        },
        claims,
        vec![],
    )
    .unwrap();
    ConductionContact {
        first_vertices: [1, 2],
        second_vertices: [4, 7],
        source: resolve_electrical_interface_state(
            &card,
            point,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap(),
    }
}

#[test]
fn contact_refresh_cannot_drop_or_rename_original_untyped_context() {
    let source = response(
        ElectricalPropertyBasis::Resistivity,
        ConductorMagnitude::CurrentDensity,
        1.0,
        0.1,
        10.0,
    );
    let context = QueryPoint::new().with("original_context", 1.0).unwrap();
    let original = contact_at(&source, &source, 0.5, &context);
    let model = PlanarConductionResponse::new_with_phase_contacts(
        mesh(1.0),
        vec![source.clone(); 4],
        vec![original.clone()],
        vec![[None, None]],
        limits(),
        &mut || false,
    )
    .unwrap();
    for point in [
        QueryPoint::new(),
        QueryPoint::new().with("substituted_context", 1.0).unwrap(),
    ] {
        let proposed = contact_at(&source, &source, 0.5, &point);
        // Typed quantity maps alone do not cover original untyped context.
        assert_eq!(
            proposed.source.resolved().axis_quantities(),
            original.source.resolved().axis_quantities()
        );
        assert!(
            model
                .solve_current_with_materials_and_contacts_candidate(
                    model.conductors().to_vec(),
                    vec![proposed],
                    &[5, 6],
                    &[0, 3],
                    1.0,
                    limits(),
                    controls(),
                    policy(),
                    &mut || false,
                )
                .is_err()
        );
    }
    assert_eq!(model.contacts(), &[original]);
}

#[test]
fn nonlinear_series_field_retains_surface_current_heat_and_source_identity() {
    let a = response(
        ElectricalPropertyBasis::Conductivity,
        ConductorMagnitude::ElectricField,
        1.0,
        0.2,
        10.0,
    );
    let b = response(
        ElectricalPropertyBasis::Conductivity,
        ConductorMagnitude::ElectricField,
        3.0,
        0.4,
        10.0,
    );
    let mut currents = Vec::new();
    for rc in [0.2, 1.3] {
        let model = PlanarConductionResponse::new_with_phase_contacts(
            mesh(1.0),
            vec![a.clone(), a.clone(), b.clone(), b.clone()],
            vec![contact(&a, &b, rc)],
            vec![[None, None]],
            limits(),
            &mut || false,
        )
        .unwrap();
        let fields = |j: f64| {
            [
                2.0 * j / (1.0 + (1.0 + 0.8 * j).sqrt()),
                2.0 * j / (3.0 + (9.0 + 1.6 * j).sqrt()),
            ]
        };
        let (mut lo, mut hi) = (0.0, 10.0);
        for _ in 0..100 {
            let j = 0.5 * (lo + hi);
            if fields(j).iter().sum::<f64>() + rc * j < 2.0 {
                lo = j;
            } else {
                hi = j;
            }
        }
        let j = 0.5 * (lo + hi);
        let e = fields(j);
        currents.push(j);
        for sign in [1.0, -1.0] {
            let voltage = [(0, 0.0), (3, 0.0), (5, sign * 2.0), (6, sign * 2.0)];
            let solved = model
                .solve_candidate(&voltage, &[], limits(), controls(), policy(), &mut || false)
                .unwrap();
            let report = solved.report();
            assert!(report.field_solves > 1);
            for cell in 0..4 {
                near(report.current_density_a_m2[cell][0], -sign * j);
                near(report.electric_field_v_m[cell][0], -sign * e[cell / 2]);
                near(report.cell_power_w[cell], 0.5 * j * e[cell / 2]);
            }
            near(report.contact_response[0].current_a, -sign * j);
            near(report.contact_response[0].voltage_jump_v[0], -sign * j * rc);
            near(report.contact_response[0].power_w, j * j * rc);
            near(report.power_w, 2.0 * j);
            near(
                report.cell_power_w.iter().sum::<f64>() + report.contact_response[0].power_w,
                report.power_w,
            );
            let current = model
                .solve_current_candidate(
                    &[5, 6],
                    &[0, 3],
                    sign * j,
                    limits(),
                    controls(),
                    policy(),
                    &mut || false,
                )
                .unwrap();
            near(current.report().voltage_v[5], sign * 2.0);
            let (native, sources) = solved.into_parts().0.into_frozen_parts();
            let continued =
                PlanarConductionResponse::from_frozen_parts(native, sources, &mut || false)
                    .unwrap();
            assert_eq!(continued.contacts(), model.contacts());
            assert!(
                continued
                    .solve_with_materials_and_contacts_candidate(
                        continued.conductors().to_vec(),
                        vec![contact(&a, &b, rc + 1.0)],
                        &voltage,
                        &[],
                        limits(),
                        controls(),
                        policy(),
                        &mut || false,
                    )
                    .is_err()
            );
        }
    }
    assert!(currents[0] > currents[1]);
}

#[test]
fn nonlinear_bulk_and_source_contacts_share_current_voltage_jump_and_heat() {
    let r = [1.0, 2.0];
    let slope = [0.05, 0.1];
    let depth = 0.4;
    let mut currents = Vec::new();
    for contact_scale in [1.0, 3.0] {
        let a = response(
            ElectricalPropertyBasis::Resistivity,
            ConductorMagnitude::CurrentDensity,
            r[0],
            slope[0],
            10.0,
        );
        let b = response(
            ElectricalPropertyBasis::Resistivity,
            ConductorMagnitude::CurrentDensity,
            r[1],
            slope[1],
            10.0,
        );
        let contact = binding(interface(a.original(), b.original(), 1e5, contact_scale));
        let model = PlanarConductionResponse::new_with_phase_contacts(
            mesh(depth),
            vec![a.clone(), a, b.clone(), b],
            vec![contact.clone()],
            vec![[None, None]],
            limits(),
            &mut || false,
        )
        .unwrap();
        let original = model.conductors().to_vec();
        let resistance = r[0] + r[1] + contact.source.resistance_ohm_m2();
        // Independent series equation: |V| = (rho1+rho2+R'')j+(b1+b2)j^2.
        let j =
            8.0 / (resistance + (resistance * resistance + 16.0 * (slope[0] + slope[1])).sqrt());
        for voltage in [4.0, -4.0] {
            let solved = model
                .solve_candidate(
                    &[(0, 0.0), (3, 0.0), (5, voltage), (6, voltage)],
                    &[],
                    limits(),
                    controls(),
                    policy(),
                    &mut || false,
                )
                .unwrap();
            let report = solved.report();
            assert!(report.field_solves > 1);
            assert_eq!(solved.problem().contacts(), &[contact.clone()]);
            for triangle in 0..4 {
                let region = triangle / 2;
                near(
                    report.current_density_a_m2[triangle][0],
                    -voltage.signum() * j,
                );
                near(
                    report.cell_power_w[triangle],
                    0.5 * depth * j * j * (r[region] + slope[region] * j),
                );
                near(
                    report.evaluated[triangle]
                        .conductor
                        .material()
                        .query()
                        .axes()["operating"],
                    j,
                );
            }
            let surface = &report.contact_response[0];
            near(surface.current_a, -voltage.signum() * depth * j);
            near(
                surface.power_w,
                depth * j * j * contact.source.resistance_ohm_m2(),
            );
            near(
                surface.voltage_jump_v[0],
                -voltage.signum() * j * contact.source.resistance_ohm_m2(),
            );
            near(surface.voltage_jump_v[1], surface.voltage_jump_v[0]);
            near(report.power_w, voltage.abs() * depth * j);
            near(
                report.cell_power_w.iter().sum::<f64>() + surface.power_w,
                report.power_w,
            );
            let driven = model
                .solve_current_candidate(
                    &[5, 6],
                    &[0, 3],
                    voltage.signum() * depth * j,
                    limits(),
                    controls(),
                    policy(),
                    &mut || false,
                )
                .unwrap();
            near(driven.report().voltage_v[5], voltage);
            near(driven.report().power_w, report.power_w);
            assert_eq!(model.conductors(), original);
        }
        currents.push(j);
        let (native, sources) = model.into_frozen_parts();
        let restored =
            PlanarConductionResponse::from_frozen_parts(native, sources, &mut || false).unwrap();
        assert_eq!(restored.contacts(), &[contact]);
    }
    assert!(currents[0] > currents[1]);
}

#[test]
fn contact_context_refresh_keeps_original_source_and_refuses_missing_or_replaced_owners() {
    let a = response(
        ElectricalPropertyBasis::Resistivity,
        ConductorMagnitude::CurrentDensity,
        1.0,
        0.1,
        10.0,
    );
    let b = response(
        ElectricalPropertyBasis::Resistivity,
        ConductorMagnitude::CurrentDensity,
        2.0,
        0.1,
        10.0,
    );
    let original = binding(interface(a.original(), b.original(), 1e5, 1.0));
    let refreshed = binding(interface(a.original(), b.original(), 2e5, 1.0));
    let foreign = binding(interface(a.original(), b.original(), 2e5, 2.0));
    let model = PlanarConductionResponse::new_with_phase_contacts(
        mesh(1.0),
        vec![a.clone(), a, b.clone(), b],
        vec![original.clone()],
        vec![[None, None]],
        limits(),
        &mut || false,
    )
    .unwrap();
    let solve = |contacts: &[fs_em::conduction::ConductionContact]| {
        model.solve_coupled_with_materials_and_contacts_candidate(
            model.conductors().to_vec(),
            contacts.to_vec(),
            limits(),
            controls(),
            policy(),
            &mut || false,
            |native, cancelled| {
                let field =
                    native.solve_current(&[5, 6], &[0, 3], 2.0, controls(), &mut || cancelled())?;
                Ok::<_, ConductionResponseError>((field, ()))
            },
        )
    };
    let (next, ()) = solve(&[refreshed.clone()]).unwrap();
    near(next.report().voltage_v[5], 2.0 * (1.2 + 2.2 + 0.25));
    assert_eq!(next.problem().contacts(), &[refreshed]);
    assert!(solve(&[]).is_err());
    assert!(solve(&[foreign]).is_err());
    assert_eq!(model.contacts(), &[original]);
}
