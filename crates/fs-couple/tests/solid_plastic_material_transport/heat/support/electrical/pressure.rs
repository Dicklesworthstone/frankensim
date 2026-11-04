//! G1/G3: native support compression closes the electrical/thermal/plastic loop.
#[path = "pressure/dynamics.rs"]
mod dynamics;
use super::*;
use fs_material::electrical_interface::{
    ELECTRICAL_CONTACT_RESISTANCE_DIMS as R_DIMS, ELECTRICAL_CONTACT_RESISTANCE_PROPERTY as R_NAME,
};

fn with_current_support<R>(
    resistance: f64,
    electrical_pressure_limit: f64,
    charge: f64,
    fraction: f64,
    run: impl FnOnce(
        &Cx<'_>,
        &CancelGate,
        &DepositionSimulation<'_>,
        &ThermoplasticSolid<'_>,
        &ThermalSupportContact<'_>,
        ElectricalSupportContact,
        &PreparedHeatNetwork<'_, '_>,
        &MaterialCircuit,
    ) -> R,
) -> R {
    powered(
        5.0,
        8e6,
        false,
        charge,
        |cx, gate, t, s, support, prepared, original| {
            let mut claim = scalar(R_NAME, R_DIMS, resistance);
            claim.value = PropertyValue::Curve {
                abscissa: "pressure".into(),
                abscissa_dims: Pressure::DIMS,
                knots: vec![
                    (0.0, resistance),
                    (electrical_pressure_limit, 0.25 * resistance),
                ],
                dims: R_DIMS,
            };
            claim.validity = ValidityDomain::unconstrained()
                .with("pressure", 0.0, electrical_pressure_limit)
                .with("Ta", 280.0, 320.0)
                .with("Tb", 280.0, 320.0);
            claim.interpolation = InterpolationPolicy::LinearInside;
            let old = support.point.interface;
            let mut claims = old.claims().clone();
            let pin = claims.insert_claim(claim).unwrap();
            let card = InterfaceSystemCard::assemble(
                old.surface_a().clone(),
                old.surface_b().clone(),
                old.context().clone(),
                claims,
                vec![],
            )
            .unwrap();
            let support = ThermalSupportContact {
                link: support.link,
                face_nodes: support.face_nodes,
                point: ThermalInterfacePoint {
                    interface: &card,
                    ..support.point
                },
                caloric: support.caloric,
                claim: support.claim,
                pressure_axis: support.pressure_axis,
            };
            let electrical = ElectricalSupportContact {
                link: support.link,
                branch: 1,
                claim: pin,
                heat_fraction_to_first: fraction,
            };
            // The context still says 7 MPa. This must resolve from actual accepted
            // stress (zero), before the circuit even exists with this contact.
            let initial = s
                .resolve_electrical_supports(
                    t,
                    std::slice::from_ref(&support),
                    &[electrical],
                    &displacements(0.0, 0.0),
                    cx,
                )
                .unwrap();
            close(initial[0].resistance_ohm(), resistance / LINK_AREA, 1e-10);
            let circuit = original.rebind_contacts_candidate(initial, &[]).unwrap();
            run(cx, gate, t, s, &support, electrical, prepared, &circuit)
        },
    )
}

// Independent scalar equations, not a second invocation of any coupled owner.
// RC midpoint fixes current; two finite heat capacities receive bulk/contact
// Joule heat once; an analytic scalar J2 return supplies physical plastic heat.
fn current_oracle(
    resistance: f64,
    compression: f64,
    charge: f64,
    fraction: f64,
    plastic: bool,
) -> [f64; 5] {
    let capacity = MASS * 5.0;
    let pressure = 2e9 * compression;
    let contact_r = resistance * (1.0 - 0.75 * pressure / 8e6) / LINK_AREA;
    let thermal_r = 0.5 - 0.4 * pressure / 8e6;
    let bulk_g = LINK_AREA / (L / (2.0 * 0.01));
    let g = 1.0 / (1.0 / bulk_g + thermal_r / LINK_AREA);
    let radiation = 5.670374419e-8 * L * L / (2.0 / 0.7 - 1.0);
    let exchange = |a: f64, b: f64| g * (a - b) + radiation * (a.powi(4) - b.powi(4));
    let electrical = |a: f64| {
        let bulk_r = 5.0 * (1.0 + 0.01 * (a - 300.0)) * L / CURRENT_AREA;
        let r = bulk_r + contact_r;
        let z = 1.0 / (2.0 * r * CAP);
        let q = charge * (1.0 - z) / (1.0 + z);
        let current = (charge + q) / (2.0 * CAP * r);
        (q, current * current * bulk_r, current * current * contact_r)
    };
    let lower = |a| {
        let contact_heat = electrical(a).2;
        let (mut lo, mut hi) = (280.0, 320.0);
        for _ in 0..80 {
            let b = 0.5 * (lo + hi);
            if capacity * (b - 300.0) < (1.0 - fraction) * contact_heat + exchange(a, b) {
                lo = b;
            } else {
                hi = b;
            }
        }
        0.5 * (lo + hi)
    };
    let residual = |a| {
        let (_, bulk_heat, contact_heat) = electrical(a);
        let plastic_heat = if plastic {
            let (_, alpha, fy) = apex_oracle(a);
            VOLUME * fy * alpha
        } else {
            0.0
        };
        capacity * (a - 300.0) + exchange(a, lower(a))
            - plastic_heat
            - bulk_heat
            - fraction * contact_heat
    };
    let (mut lo, mut hi) = (300.0, 319.0);
    assert!(residual(lo) < 0.0 && residual(hi) > 0.0);
    for _ in 0..80 {
        let a = 0.5 * (lo + hi);
        if residual(a) < 0.0 {
            lo = a;
        } else {
            hi = a;
        }
    }
    let a = 0.5 * (lo + hi);
    [
        a,
        lower(a),
        electrical(a).0,
        electrical(a).2,
        if plastic {
            2.0 * L * apex_oracle(a).0
        } else {
            0.0
        },
    ]
}

#[test]
fn solved_pressure_and_material_only_interface_swaps_reach_current_heat_and_plasticity() {
    let mut outcomes = Vec::new();
    for (resistance, compression, charge) in [
        (0.2, 0.001, CHARGE),
        (0.2, 0.003, CHARGE),
        (0.4, 0.001, CHARGE),
        (0.2, 0.001, -CHARGE),
    ] {
        with_current_support(
            resistance,
            8e6,
            charge,
            0.3,
            |cx, _, t, s, support, electrical, p, c| {
                let next = s
                    .advance_electrical_with_pressure_contacts_candidate(
                        t,
                        c,
                        p,
                        std::slice::from_ref(support),
                        &[electrical],
                        &[],
                        1.0,
                        &displacements(compression, SHEAR),
                        &applied(),
                        IncrementSettings::default(),
                        powered_controls(),
                        &[0.0; 2],
                        cx,
                    )
                    .unwrap();
                let expected = current_oracle(resistance, compression, charge, 0.3, true);
                let r = next.report();
                let loss = &r.circuit_thermal.electrical.contacts[0];
                let actual = [
                    next.thermal().temperature_k(CELLS[0]).unwrap(),
                    next.thermal().temperature_k(CELLS[1]).unwrap(),
                    next.circuit().charges_c()[0],
                    loss.heat_j,
                    next.solid().solution().displacement_m[3][0],
                ];
                for i in 0..2 {
                    close(actual[i], expected[i], 3e-6);
                }
                close(actual[2], expected[2], 1e-10);
                close(actual[3], expected[3], 1e-10);
                close(actual[4], expected[4], 3e-11);
                close(
                    loss.contact.resistance_ohm(),
                    resistance * (1.0 - 0.75 * 2e9 * compression / 8e6) / LINK_AREA,
                    1e-8,
                );
                let property = loss.contact.source().resolved().property(R_NAME).unwrap();
                let point = &property.answer().receipt.query_point;
                assert!(point.contains(&("pressure".into(), r.support_contacts[0].compression_pa)));
                assert!(point.contains(&("Ta".into(), actual[0])));
                assert!(point.contains(&("Tb".into(), actual[1])));
                assert_eq!(property.answer().receipt.selected, electrical.claim);
                assert_eq!(r.circuit_thermal.source.destinations[0].cells, CELLS);
                close(
                    r.circuit_thermal.source.destinations[0].heat_fraction_to_first,
                    0.3,
                    0.0,
                );
                let electrical_loss = (charge * charge - actual[2] * actual[2]) / (2.0 * CAP);
                close(r.circuit_thermal.deposited_joule_j, electrical_loss, 1e-10);
                close(
                    r.circuit_thermal.thermal.stored_energy_change_j,
                    electrical_loss + r.deposited_mechanical_heat_j,
                    1e-8,
                );
                assert!(r.combined_residual_j.abs() <= 1e-8 && r.thermal_residual_j.abs() <= 1e-8);
                assert_eq!(t.time_s(), 0.0);
                assert_eq!(s.time_s(), 0.0);
                assert_eq!(c.charges_c(), &[charge]);
                outcomes.push(actual);
            },
        );
    }
    assert!((outcomes[1][2] - outcomes[0][2]).abs() > 1e-6);
    assert!((outcomes[2][2] - outcomes[0][2]).abs() > 1e-6);
    assert!((outcomes[2][4] - outcomes[0][4]).abs() > 1e-10);
    for i in [0, 1, 3, 4] {
        close(outcomes[0][i], outcomes[3][i], 1e-12);
    }
    close(outcomes[0][2], -outcomes[3][2], 1e-14);
}

#[test]
fn zero_compression_retains_sourced_finite_resistance_and_both_heat_destinations() {
    for fraction in [0.0, 1.0] {
        with_current_support(
            0.2,
            8e6,
            CHARGE,
            fraction,
            |cx, _, t, s, support, electrical, p, c| {
                let next = s
                    .advance_electrical_with_pressure_contacts_candidate(
                        t,
                        c,
                        p,
                        std::slice::from_ref(support),
                        &[electrical],
                        &[],
                        1.0,
                        &displacements(0.0, 0.0),
                        &[],
                        IncrementSettings::default(),
                        powered_controls(),
                        &[0.0; 2],
                        cx,
                    )
                    .unwrap();
                let expected = current_oracle(0.2, 0.0, CHARGE, fraction, false);
                for i in 0..2 {
                    close(
                        next.thermal().temperature_k(CELLS[i]).unwrap(),
                        expected[i],
                        1e-7,
                    );
                }
                close(next.circuit().charges_c()[0], expected[2], 1e-10);
                assert_eq!(next.report().support_contacts[0].compression_pa, 0.0);
                assert_eq!(next.report().physical_plastic_heat_j, 0.0);
                assert!(next.report().circuit_thermal.electrical.contacts[0].heat_j > 0.0);
                assert!(next.report().combined_residual_j.abs() <= 1e-8);
            },
        );
    }
}

#[test]
fn electrical_pressure_support_failure_and_bad_bindings_preserve_all_accepted_histories() {
    with_current_support(
        0.2,
        1e6,
        CHARGE,
        0.3,
        |cx, _, t, s, support, electrical, p, c| {
            let old_cells = t.cells().clone();
            let old_solid = s.solution().clone();
            let old_circuit = c.clone();
            // The thermal support admits 4 MPa; ONLY its electrical claim stops at 1 MPa.
            assert!(
                s.advance_electrical_with_pressure_contacts_candidate(
                    t,
                    c,
                    p,
                    std::slice::from_ref(support),
                    &[electrical],
                    &[],
                    1.0,
                    &displacements(0.002, SHEAR),
                    &applied(),
                    IncrementSettings::default(),
                    powered_controls(),
                    &[0.0; 2],
                    cx
                )
                .is_err()
            );
            for bad in [
                ElectricalSupportContact {
                    heat_fraction_to_first: 1.1,
                    ..electrical
                },
                ElectricalSupportContact {
                    branch: 0,
                    ..electrical
                },
                ElectricalSupportContact {
                    link: 92,
                    ..electrical
                },
                ElectricalSupportContact {
                    claim: support.claim,
                    ..electrical
                },
            ] {
                assert!(
                    s.advance_electrical_with_pressure_contacts_candidate(
                        t,
                        c,
                        p,
                        std::slice::from_ref(support),
                        &[bad],
                        &[],
                        1.0,
                        &displacements(0.0, 0.0),
                        &[],
                        IncrementSettings::default(),
                        powered_controls(),
                        &[0.0; 2],
                        cx
                    )
                    .is_err()
                );
            }
            assert_eq!(t.cells(), &old_cells);
            assert_eq!(s.solution(), &old_solid);
            assert_eq!(c, &old_circuit);
        },
    );
    with_current_support(
        0.2,
        8e6,
        CHARGE,
        0.3,
        |cx, gate, t, s, support, electrical, p, c| {
            let saved = t.cells().clone();
            let old_solid = s.solution().clone();
            let old_circuit = c.clone();
            for (compression, bindings) in [
                (-0.001, vec![electrical]),
                (0.001, vec![electrical, electrical]),
            ] {
                assert!(
                    s.advance_electrical_with_pressure_contacts_candidate(
                        t,
                        c,
                        p,
                        std::slice::from_ref(support),
                        &bindings,
                        &[],
                        1.0,
                        &displacements(compression, SHEAR),
                        &applied(),
                        IncrementSettings::default(),
                        powered_controls(),
                        &[0.0; 2],
                        cx,
                    )
                    .is_err()
                );
            }
            assert_eq!(t.cells(), &saved);
            assert_eq!(s.solution(), &old_solid);
            assert_eq!(c, &old_circuit);
            gate.request();
            assert!(matches!(
                s.advance_electrical_with_pressure_contacts_candidate(
                    t,
                    c,
                    p,
                    std::slice::from_ref(support),
                    &[electrical],
                    &[],
                    1.0,
                    &displacements(0.001, SHEAR),
                    &applied(),
                    IncrementSettings::default(),
                    powered_controls(),
                    &[0.0; 2],
                    cx
                ),
                Err(ElectricalSolidError::Cancelled)
            ));
            assert_eq!(t.time_s(), 0.0);
            assert_eq!(s.time_s(), 0.0);
            assert_eq!(c.charges_c(), &[CHARGE]);
        },
    );
}
