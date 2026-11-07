//! Authored radial specimens and native gas inventories; independent balances.
use super::*;
use crate::atmosphere::{
    atmosphere as gas, carrier_mass, controls as finite_controls, sorption as sorption_law,
};
use fs_couple::moisture::atmosphere::{
    AtmosphereMoistureError, DryBasisSorption,
    pair::{AtmosphereMoistureFace, MoistureAtmospherePort, PairedAtmosphereMoistureSolve},
};
use fs_couple::moisture::radial::atmosphere::RadialMoistureParticipant;

#[path = "atmosphere/sdirk.rs"]
mod sdirk;

fn paired_controls() -> PairedAtmosphereMoistureSolve {
    PairedAtmosphereMoistureSolve {
        atmospheres: [finite_controls(); 2],
        maximum_evaluations: 10_000,
        mass_tolerance_kg: 5e-12,
    }
}

fn port(atmosphere: MoistureAtmospherePort, slope: f64) -> Option<AtmosphereMoistureFace> {
    Some(AtmosphereMoistureFace {
        atmosphere,
        sorption: sorption_law(slope),
    })
}

fn participant<'a>(
    id: u64,
    moisture: &'a MoistureRadial,
    faces: &'a [Option<DryBasisSorption>; 2],
) -> RadialMoistureParticipant<'a> {
    RadialMoistureParticipant {
        id,
        moisture,
        faces,
        transport: config(),
        nonlinear: nonlinear(),
    }
}

fn local_water_change(before: &MoistureRadial, after: &MoistureRadial) -> f64 {
    before
        .moisture()
        .iter()
        .zip(after.moisture())
        .zip(before.dry_masses_kg())
        .map(|((&old, &new), &mass)| mass * (new - old))
        .sum()
}

#[test]
fn finite_radial_vapor_matches_analytic_uptake_desorption_and_continuation() {
    with_cx(|cx| {
        let radius = 0.05;
        let diffusion = 0.001;
        let dt = 0.2;
        let slope = 0.25;
        for finite_surface in [false, true] {
            let h = 0.004;
            let law = if finite_surface {
                crate::surface::law([h; 2])
            } else {
                sorption_law(slope)
            };
            let speed = if finite_surface {
                1.0 / (radius / diffusion + 1.0 / h)
            } else {
                diffusion / radius
            };
            let alpha = dt * 3.0 * speed / radius;
            let faces = [None, Some(law)];
            for initial_rh in [0.8, 0.2] {
                let mut air = gas(initial_rh, 1.0);
                let mut field = specimen(
                    RadialMoistureSymmetry::Sphere,
                    vec![0.0, radius],
                    vec![0.1],
                    [diffusion; 2],
                    cx,
                );
                let carrier = carrier_mass(&air);
                let identity = field.identity();
                for interval in 0..2 {
                    let rh = air.relative_humidity().unwrap();
                    let capacity = air.vapor_mass_kg().unwrap() / rh;
                    let dry = field.dry_mass_kg();
                    let initial = field.moisture()[0];
                    // Eliminate M' from B(M'-M)=B*alpha*(s*RH'-M')
                    // and C(RH'-RH)=-B(M'-M), without a diffusion surrogate.
                    let weight = dry * alpha / (1.0 + alpha);
                    let exchange =
                        weight * (slope * rh - initial) / (1.0 + weight * slope / capacity);
                    let expected = initial + exchange / dry;
                    let expected_rh = rh - exchange / capacity;
                    let before_mass = air.gas().mass_kg() + field.wet_mass_kg();
                    let candidate = air
                        .advance_radial(
                            &field,
                            dt,
                            &faces,
                            config(),
                            nonlinear(),
                            finite_controls(),
                            cx,
                        )
                        .unwrap();
                    near(candidate.moisture().moisture()[0], expected, 3e-12);
                    near(
                        candidate.atmosphere().relative_humidity().unwrap(),
                        expected_rh,
                        3e-11,
                    );
                    near(candidate.step().water_to_solid_kg, exchange, 3e-12);
                    near(
                        candidate.step().transport.boundary_water_kg[1],
                        exchange,
                        3e-12,
                    );
                    assert_eq!(candidate.step().transport.boundary_water_kg[0], 0.0);
                    near(carrier_mass(candidate.atmosphere()), carrier, 3e-13);
                    near(
                        candidate.atmosphere().gas().mass_kg() + candidate.moisture().wet_mass_kg(),
                        before_mass,
                        3e-12,
                    );
                    near(
                        candidate.step().water_balance_residual_kg,
                        0.0,
                        finite_controls().mass_tolerance_kg,
                    );
                    let source = candidate.step().sorption[1].as_ref().unwrap();
                    near(
                        source.query().axes()["relative_humidity"],
                        expected_rh,
                        3e-11,
                    );
                    assert_eq!(
                        candidate.step().transport.surface_transfer[1].is_some(),
                        finite_surface
                    );
                    assert_eq!(candidate.moisture().identity(), identity);
                    assert_eq!(
                        candidate.step().transport.interval_s,
                        [interval as f64 * dt, (interval + 1) as f64 * dt]
                    );
                    let unlimited_change = alpha / (1.0 + alpha) * (slope * rh - initial);
                    assert!(
                        (candidate.moisture().moisture()[0] - initial).abs()
                            < unlimited_change.abs()
                    );
                    let (next_air, next_field, _) = candidate.into_parts();
                    air = next_air;
                    field = next_field;
                }
            }
        }
    });
}

#[test]
fn paired_radial_chambers_retain_through_flow_when_solid_water_is_unchanged() {
    with_cx(|cx| {
        let a: f64 = 0.03;
        let b: f64 = 0.05;
        let center = 0.5 * (a + b);
        let diffusion = 1e-5;
        let dt = 0.5;
        let slope = [0.22, 0.28];
        let air = [gas(0.85, 0.6), gas(0.15, 1.1)];
        let rh = air.each_ref().map(|g| g.relative_humidity().unwrap());
        let capacity: [f64; 2] = std::array::from_fn(|i| air[i].vapor_mass_kg().unwrap() / rh[i]);
        for symmetry in [
            RadialMoistureSymmetry::Sphere,
            RadialMoistureSymmetry::Cylinder { length_m: 0.3 },
        ] {
            let (factor, resistance) = match symmetry {
                RadialMoistureSymmetry::Sphere => (
                    4.0 * core::f64::consts::PI,
                    [1.0 / a - 1.0 / center, 1.0 / center - 1.0 / b],
                ),
                RadialMoistureSymmetry::Cylinder { length_m } => (
                    2.0 * core::f64::consts::PI * length_m,
                    [(center / a).ln(), (b / center).ln()],
                ),
            };
            // Native curved half-shell resistances give q_i=lambda_i*(s_i*r_i-M).
            // Eliminate each finite gas to obtain q_i=w_i*(s_i*r_i0-M).
            let lambda = resistance.map(|r| dt * 500.0 * factor * diffusion / r);
            let weight: [f64; 2] =
                std::array::from_fn(|i| lambda[i] / (1.0 + lambda[i] * slope[i] / capacity[i]));
            let initial = (weight[0] * slope[0] * rh[0] + weight[1] * slope[1] * rh[1])
                / (weight[0] + weight[1]);
            let transfer: [f64; 2] =
                std::array::from_fn(|i| weight[i] * (slope[i] * rh[i] - initial));
            let expected_rh: [f64; 2] = std::array::from_fn(|i| rh[i] - transfer[i] / capacity[i]);
            let field = specimen(symmetry, vec![a, b], vec![initial], [diffusion; 2], cx);
            let faces = [
                port(MoistureAtmospherePort::First, slope[0]),
                port(MoistureAtmospherePort::Second, slope[1]),
            ];
            let candidate = air[0]
                .advance_paired_radial(
                    &air[1],
                    &field,
                    dt,
                    &faces,
                    config(),
                    nonlinear(),
                    paired_controls(),
                    cx,
                )
                .unwrap();
            near(candidate.moisture().moisture()[0], initial, 3e-12);
            near(candidate.step().water_change_kg, 0.0, 3e-12);
            assert!(candidate.step().atmospheres[0].water_to_solids_kg > 1e-6);
            assert!(candidate.step().atmospheres[1].water_to_solids_kg < -1e-6);
            for i in 0..2 {
                near(
                    candidate.atmospheres()[i].relative_humidity().unwrap(),
                    expected_rh[i],
                    3e-11,
                );
                near(
                    candidate.step().atmospheres[i].water_to_solids_kg,
                    transfer[i],
                    3e-12,
                );
                near(
                    candidate.step().transport.boundary_water_kg[i],
                    transfer[i],
                    3e-12,
                );
                near(
                    carrier_mass(&candidate.atmospheres()[i]),
                    carrier_mass(&air[i]),
                    3e-13,
                );
            }
            near(
                candidate.step().water_balance_residual_kg,
                0.0,
                paired_controls().mass_tolerance_kg,
            );
            near(
                candidate.step().total_mass_residual_kg,
                0.0,
                paired_controls().mass_tolerance_kg,
            );
            // Swap owner order and labels together, preserving physical inner/outer contact.
            let reversed = air[1]
                .advance_paired_radial(
                    &air[0],
                    &field,
                    dt,
                    &[
                        port(MoistureAtmospherePort::Second, slope[0]),
                        port(MoistureAtmospherePort::First, slope[1]),
                    ],
                    config(),
                    nonlinear(),
                    paired_controls(),
                    cx,
                )
                .unwrap();
            near(reversed.moisture().moisture()[0], initial, 3e-12);
            near(
                reversed.atmospheres()[1].relative_humidity().unwrap(),
                expected_rh[0],
                3e-11,
            );
            near(
                reversed.atmospheres()[0].relative_humidity().unwrap(),
                expected_rh[1],
                3e-11,
            );
        }
    });
}

#[test]
fn shared_radial_chamber_conserves_actual_shell_inventories_and_order() {
    with_cx(|cx| {
        let dry = specimen(
            RadialMoistureSymmetry::Sphere,
            vec![0.0, 0.003, 0.009, 0.02],
            vec![0.01, 0.03, 0.06],
            [1e-5; 2],
            cx,
        );
        let wet = specimen(
            RadialMoistureSymmetry::Cylinder { length_m: 0.04 },
            vec![0.0, 0.005, 0.015, 0.025],
            vec![0.24, 0.22, 0.20],
            [2e-5; 2],
            cx,
        );
        let air = gas(0.5, 0.5);
        let faces = [None, Some(sorption_law(0.25))];
        let candidate = air
            .advance_shared_radial(
                &[participant(9, &dry, &faces), participant(2, &wet, &faces)],
                0.1,
                6,
                finite_controls(),
                cx,
            )
            .unwrap();
        let reordered = air
            .advance_shared_radial(
                &[participant(2, &wet, &faces), participant(9, &dry, &faces)],
                0.1,
                6,
                finite_controls(),
                cx,
            )
            .unwrap();
        assert_eq!(candidate.atmosphere(), reordered.atmosphere());
        assert_eq!(
            candidate.solids().iter().map(|s| s.id).collect::<Vec<_>>(),
            vec![2, 9]
        );
        let mut total = 0.0;
        for ((before, solid), replay) in [&wet, &dry]
            .into_iter()
            .zip(candidate.solids())
            .zip(reordered.solids())
        {
            let change = local_water_change(before, &solid.moisture);
            near(
                change,
                solid.transport.boundary_water_kg.iter().sum(),
                3e-12,
            );
            assert_eq!(solid.moisture.moisture(), replay.moisture.moisture());
            assert_eq!(solid.moisture.time_s(), 0.1);
            total += change;
            assert_eq!(change.is_sign_positive(), solid.id == 9);
        }
        near(candidate.step().water_to_solids_kg, total, 1e-14);
        near(
            candidate.atmosphere().vapor_mass_kg().unwrap() - air.vapor_mass_kg().unwrap() + total,
            0.0,
            finite_controls().mass_tolerance_kg,
        );
        near(
            candidate.step().total_mass_residual_kg,
            0.0,
            finite_controls().mass_tolerance_kg,
        );
        let mixed_clocks = [
            participant(9, &dry, &faces),
            participant(2, &candidate.solids()[0].moisture, &faces),
        ];
        assert!(matches!(
            air.advance_shared_radial(&mixed_clocks, 0.1, 6, finite_controls(), cx),
            Err(AtmosphereMoistureError::Invalid(_))
        ));
        assert!(matches!(
            air.advance_shared_radial(
                &[participant(9, &dry, &faces), participant(2, &wet, &faces)],
                0.1,
                5,
                finite_controls(),
                cx
            ),
            Err(AtmosphereMoistureError::Invalid(_))
        ));
    });
}

#[test]
fn sealed_radial_field_redistributes_water_without_debiting_the_gas() {
    with_cx(|cx| {
        let field = specimen(
            RadialMoistureSymmetry::Sphere,
            vec![0.0, 0.01, 0.03],
            vec![0.2, 0.05],
            [1e-4; 2],
            cx,
        );
        let air = gas(0.7, 0.3);
        let candidate = air
            .advance_radial(
                &field,
                0.1,
                &[None, None],
                config(),
                nonlinear(),
                finite_controls(),
                cx,
            )
            .unwrap();
        assert_eq!(candidate.atmosphere(), &air);
        assert_eq!(candidate.step().iterations, 0);
        assert!(candidate.moisture().moisture()[0] < field.moisture()[0]);
        assert!(candidate.moisture().moisture()[1] > field.moisture()[1]);
        near(
            candidate.moisture().water_mass_kg(),
            field.water_mass_kg(),
            1e-13,
        );
        near(
            candidate.step().water_to_solid_kg,
            0.0,
            finite_controls().mass_tolerance_kg,
        );
    });
}

#[test]
fn finite_radial_atmosphere_wets_a_degenerate_dry_front_without_creating_water() {
    with_cx(|cx| {
        let field = specimen(
            RadialMoistureSymmetry::Sphere,
            vec![0.0, 0.01, 0.03],
            vec![0.0; 2],
            [0.0, 0.03],
            cx,
        );
        let air = gas(0.8, 1.0);
        let candidate = air
            .advance_radial(
                &field,
                0.01,
                &[None, Some(sorption_law(0.25))],
                config(),
                nonlinear(),
                finite_controls(),
                cx,
            )
            .unwrap();
        assert!(candidate.moisture().moisture()[1] > 0.0);
        assert!(candidate.moisture().moisture()[0] <= candidate.moisture().moisture()[1]);
        assert!(candidate.atmosphere().vapor_mass_kg().unwrap() < air.vapor_mass_kg().unwrap());
        assert!(candidate.step().transport.kirchhoff_source.is_some());
        near(
            candidate.atmosphere().vapor_mass_kg().unwrap() + candidate.moisture().water_mass_kg(),
            air.vapor_mass_kg().unwrap(),
            3e-12,
        );
        assert_eq!(field.moisture(), &[0.0; 2]);
        assert_eq!(field.time_s(), 0.0);
    });
}

#[test]
fn radial_source_center_and_nested_budget_refusals_are_atomic() {
    with_cx(|cx| {
        let field = specimen(
            RadialMoistureSymmetry::Sphere,
            vec![0.0, 0.05],
            vec![0.1],
            [0.001; 2],
            cx,
        );
        let air = gas(0.5, 1.0);
        let saved = air.clone();
        assert!(matches!(
            air.advance_radial(
                &field,
                0.2,
                &[Some(sorption_law(0.25)), None],
                config(),
                nonlinear(),
                finite_controls(),
                cx
            ),
            Err(AtmosphereMoistureError::Invalid(_))
        ));
        // Initial RH and the lower probe are supported; the upper probe asks
        // the radial source for M_eq=0.4, beyond its explicit M<=0.3 domain.
        assert!(matches!(
            air.advance_radial(
                &field,
                0.2,
                &[None, Some(sorption_law(0.4))],
                config(),
                nonlinear(),
                finite_controls(),
                cx
            ),
            Err(AtmosphereMoistureError::Moisture(_))
        ));
        assert_eq!(air, saved);
        assert_eq!(field.moisture(), &[0.1]);
        assert_eq!(field.time_s(), 0.0);
        let shell = specimen(
            RadialMoistureSymmetry::Sphere,
            vec![0.03, 0.05],
            vec![0.1],
            [0.001; 2],
            cx,
        );
        let second = gas(0.2, 0.7);
        let saved_second = second.clone();
        let mut limited = paired_controls();
        limited.maximum_evaluations = 1;
        assert!(matches!(
            air.advance_paired_radial(
                &second,
                &shell,
                0.2,
                &[
                    port(MoistureAtmospherePort::First, 0.25),
                    port(MoistureAtmospherePort::Second, 0.25)
                ],
                config(),
                nonlinear(),
                limited,
                cx
            ),
            Err(AtmosphereMoistureError::EvaluationBudget { .. })
        ));
        assert_eq!(air, saved);
        assert_eq!(second, saved_second);
        assert_eq!(shell.moisture(), &[0.1]);
        assert_eq!(shell.time_s(), 0.0);
    });
}

#[test]
fn finite_radial_vapor_refuses_supersaturation_without_condensation_physics() {
    with_cx(|cx| {
        let field = specimen(
            RadialMoistureSymmetry::Sphere,
            vec![0.0, 0.1],
            vec![0.28],
            [0.01; 2],
            cx,
        );
        let air = gas(0.95, 0.001);
        let saved = air.clone();
        assert!(matches!(
            air.advance_radial(
                &field,
                0.5,
                &[None, Some(sorption_law(0.25))],
                config(),
                nonlinear(),
                finite_controls(),
                cx
            ),
            Err(AtmosphereMoistureError::Unbracketed { .. })
        ));
        assert_eq!(air, saved);
        assert_eq!(field.moisture(), &[0.28]);
        assert_eq!(field.time_s(), 0.0);
    });
}
