//! Independent simplex balances with signed heat, live contact and radiation.

use super::*;
use fs_conduction::transient::enthalpy::heterogeneous::EnthalpyRadiation;
use fs_math::dd::Dd;

const HEAT: [f64; 8] = [0.3, -0.2, 0.1, -0.4, -0.15, 0.25, -0.35, 0.65];

fn history_with_heat(exterior: [f64; 8]) -> Vec<f64> {
    manufactured_history(exterior)
        .into_iter()
        .enumerate()
        .map(|(v, h)| h - HEAT[v] / if v < 4 { 1. / 24. } else { 1. / 12. })
        .collect()
}

fn storage<'m, 'c>(
    cx: &Cx<'_>,
    mesh: &'m ConductionMesh,
    cold: &'c EquilibriumEnthalpyPhaseCurve,
    hot: &'c EquilibriumEnthalpyPhaseCurve,
) -> HeterogeneousEnthalpyBackwardEuler<'m, 'c> {
    HeterogeneousEnthalpyBackwardEuler::new(
        cx,
        mesh,
        &[
            ReferenceEnthalpyMaterial {
                curve: cold,
                reference_density_kg_m3: 1.,
            },
            ReferenceEnthalpyMaterial {
                curve: hot,
                reference_density_kg_m3: 2.,
            },
        ],
        &[0, 1],
        EnthalpyBudget {
            max_vertices: 8,
            max_elements: 2,
        },
    )
    .unwrap()
}

fn assert_endpoint(
    next: &fs_conduction::transient::enthalpy::nodal_heat::NodalHeatStep,
    offsets: [f64; 2],
) {
    for v in 0..8 {
        close(
            next.specific_enthalpy_j_kg[v],
            TARGET_H[v] + offsets[v / 4],
            3e-6,
        );
        close(next.endpoint.temperature[v], TARGET_T[v], 3e-7);
    }
    assert!(next.endpoint.residual_norm_j <= coupled_newton().residual_tolerance_j);
    assert!(next.endpoint.maximum_free_residual_j <= coupled_newton().step.energy_tolerance_j);
    close(next.endpoint.energy.residual_j.to_f64(), 0., 1e-8);
    close(
        next.endpoint.energy.nodal_heat_j.to_f64(),
        HEAT.iter().sum(),
        1e-14,
    );
}

#[test]
fn g1_signed_heat_live_contact_and_ambient_recover_original_nonuniform_latent_endpoint() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh();
        let boundary = ambient_boundary(&mesh);
        let card = card();
        let contacts =
            TemperatureContacts::new(cx, &mesh, &boundary, vec![], vec![law(&mesh, &card)], 1)
                .unwrap();
        let material = ConductivityModel::isotropic_declared(2.).unwrap();
        let source = ScalarField::Uniform(0.);
        let problem = ConductionProblem {
            mesh: &mesh,
            boundary: &boundary,
            material: &material,
            element_materials: None,
            source: &source,
        };
        let heat = HEAT.map(Dd::from_f64);
        let convection = convection_load();
        for ambient in [300., 340.] {
            let radiation_load = ambient_load(ambient);
            let exterior = std::array::from_fn(|v| radiation_load[v] + convection[v]);
            let history = history_with_heat(exterior);
            let patches =
                [
                    AmbientRadiationPatch::new("exterior", emissivity("exterior", 0.2), ambient)
                        .unwrap(),
                ];
            let radiation = PointwiseAmbientRadiation::new(cx, &mesh, &boundary, &patches).unwrap();
            for offsets in [[0., 0.], [1200., -800.]] {
                let cold = shifted_chart(false, offsets[0]);
                let hot = shifted_chart(true, offsets[1]);
                let engine = storage(cx, &mesh, &cold, &hot);
                let old: Vec<_> = history
                    .iter()
                    .enumerate()
                    .map(|(v, h)| h + offsets[v / 4])
                    .collect();
                let before = old.clone();
                let next = engine
                    .advance_with_nodal_heat_and_temperature_contacts(
                        cx,
                        problem,
                        &contacts,
                        &old,
                        DT,
                        coupled_newton(),
                        &heat,
                        EnthalpyRadiation::PointwiseAmbient(&radiation),
                    )
                    .unwrap();
                assert_endpoint(&next, offsets);
                let energy = &next.endpoint.energy;
                close(
                    energy.robin_out_j.to_f64(),
                    DT * convection.iter().sum::<f64>(),
                    1e-8,
                );
                close(
                    energy.radiation_out_j.to_f64(),
                    DT * radiation_load.iter().sum::<f64>(),
                    1e-8,
                );
                close(
                    energy.stored_energy_change_j.to_f64(),
                    HEAT.iter().sum::<f64>() - DT * exterior.iter().sum::<f64>(),
                    2e-8,
                );
                assert!(next.endpoint.pointwise_radiation.is_some());
                assert!(next.endpoint.pointwise_enclosure_radiation.is_none());
                let actual = next.endpoint.temperature_contacts.as_ref().unwrap();
                let fresh = contacts
                    .evaluate_enthalpy(cx, &boundary, &engine, &next.specific_enthalpy_j_kg)
                    .unwrap();
                assert_eq!(actual.states(), fresh.states());
                assert_eq!(
                    actual.states()[0].resistance.card_identity(),
                    card.content_hash()
                );
                let old_contacts = contacts
                    .evaluate_enthalpy(cx, &boundary, &engine, &old)
                    .unwrap();
                assert!(
                    (actual.states()[0].interface_temperature_k
                        - old_contacts.states()[0].interface_temperature_k)
                        .abs()
                        > 1e-3
                );
                assert_eq!(old, before);
            }
        }
    });
}

#[test]
fn g1_signed_heat_and_live_contact_retain_original_reflected_enclosure_exchange() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh();
        let boundary = enclosure_boundary(&mesh);
        let card = card();
        let contacts =
            TemperatureContacts::new(cx, &mesh, &boundary, vec![], vec![law(&mesh, &card)], 1)
                .unwrap();
        let cold = chart(false);
        let hot = chart(true);
        let engine = storage(cx, &mesh, &cold, &hot);
        let material = ConductivityModel::isotropic_declared(2.).unwrap();
        let source = ScalarField::Uniform(0.);
        let problem = ConductionProblem {
            mesh: &mesh,
            boundary: &boundary,
            material: &material,
            element_materials: None,
            source: &source,
        };
        let epsilon = [0.2, 0.7];
        let enclosure = enclosure(&mesh, epsilon);
        let radiation = PointwiseEnclosureRadiation::new(
            cx,
            &mesh,
            &boundary,
            &enclosure,
            PointwiseEnclosureTolerance::default(),
        )
        .unwrap();
        let old = history_with_heat(enclosure_load(epsilon));
        let before = old.clone();
        let next = engine
            .advance_with_nodal_heat_and_temperature_contacts(
                cx,
                problem,
                &contacts,
                &old,
                DT,
                coupled_newton(),
                &HEAT.map(Dd::from_f64),
                EnthalpyRadiation::PointwiseEnclosure(&radiation),
            )
            .unwrap();
        assert_endpoint(&next, [0., 0.]);
        close(next.endpoint.energy.radiation_out_j.to_f64(), 0., 2e-8);
        close(
            next.endpoint.energy.stored_energy_change_j.to_f64(),
            HEAT.iter().sum(),
            2e-8,
        );
        assert!(next.endpoint.pointwise_enclosure_radiation.is_some());
        assert!(next.endpoint.pointwise_radiation.is_none());
        assert!(next.endpoint.temperature_contacts.is_some());
        assert_eq!(old, before);
    });
}

#[test]
fn g3_actual_heat_recheck_exposes_equal_and_opposite_local_defects() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let card = card();
        let contacts =
            TemperatureContacts::new(cx, &mesh, &boundary, vec![], vec![law(&mesh, &card)], 1)
                .unwrap();
        let cold = chart(false);
        let hot = chart(true);
        let engine = storage(cx, &mesh, &cold, &hot);
        let material = ConductivityModel::isotropic_declared(2.).unwrap();
        let source = ScalarField::Uniform(0.);
        let problem = ConductionProblem {
            mesh: &mesh,
            boundary: &boundary,
            material: &material,
            element_materials: None,
            source: &source,
        };
        let old = history_with_heat([0.; 8]);
        let heat = HEAT.map(Dd::from_f64);
        let next = engine
            .advance_with_nodal_heat_and_temperature_contacts(
                cx,
                problem,
                &contacts,
                &old,
                DT,
                coupled_newton(),
                &heat,
                EnthalpyRadiation::None,
            )
            .unwrap();
        assert_endpoint(&next, [0., 0.]);
        let mut actual_heat = heat;
        actual_heat[0] = actual_heat[0] + Dd::from_f64(0.5);
        actual_heat[4] = actual_heat[4] - Dd::from_f64(0.5);
        let measured = engine
            .evaluate_endpoint_with_nodal_heat_and_temperature_contacts(
                cx,
                problem,
                &contacts,
                &old,
                &next.specific_enthalpy_j_kg,
                DT,
                coupled_newton(),
                &actual_heat,
                EnthalpyRadiation::None,
            )
            .unwrap();
        close(measured.energy.residual_j.to_f64(), 0., 1e-8);
        close(measured.residual_j[0], -0.5, 1e-8);
        close(measured.residual_j[4], 0.5, 1e-8);
        assert!(measured.maximum_free_residual_j > 0.49);
        assert!(measured.residual_norm_j > 0.70);
        assert_eq!(
            measured.temperature_contacts.as_ref().unwrap().states(),
            next.endpoint
                .temperature_contacts
                .as_ref()
                .unwrap()
                .states()
        );
    });
}

#[test]
fn g4_signed_heat_live_contact_keeps_work_limits_and_history_atomic() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let card = card();
        let contacts =
            TemperatureContacts::new(cx, &mesh, &boundary, vec![], vec![law(&mesh, &card)], 1)
                .unwrap();
        let cold = chart(false);
        let hot = chart(true);
        let engine = storage(cx, &mesh, &cold, &hot);
        let material = ConductivityModel::isotropic_declared(2.).unwrap();
        let source = ScalarField::Uniform(0.);
        let problem = ConductionProblem {
            mesh: &mesh,
            boundary: &boundary,
            material: &material,
            element_materials: None,
            source: &source,
        };
        let old = history_with_heat([0.; 8]);
        let before = old.clone();
        let heat = HEAT.map(Dd::from_f64);
        let limited =
            fs_conduction::interface::temperature::enthalpy::ContactEnthalpyNewtonConfig {
                max_residual_evaluations: 1,
                ..coupled_newton()
            };
        assert!(
            engine
                .advance_with_nodal_heat_and_temperature_contacts(
                    cx,
                    problem,
                    &contacts,
                    &old,
                    DT,
                    limited,
                    &heat,
                    EnthalpyRadiation::None,
                )
                .is_err()
        );
        assert!(
            engine
                .advance_with_nodal_heat_and_temperature_contacts(
                    cx,
                    problem,
                    &contacts,
                    &old,
                    DT,
                    coupled_newton(),
                    &heat[..7],
                    EnthalpyRadiation::None,
                )
                .is_err()
        );
        let mut invalid = heat;
        invalid[2] = Dd::from_f64(f64::NAN);
        assert!(
            engine
                .advance_with_nodal_heat_and_temperature_contacts(
                    cx,
                    problem,
                    &contacts,
                    &old,
                    DT,
                    coupled_newton(),
                    &invalid,
                    EnthalpyRadiation::None,
                )
                .is_err()
        );
        assert!(
            engine
                .evaluate_endpoint_with_nodal_heat_and_temperature_contacts(
                    cx,
                    problem,
                    &contacts,
                    &old,
                    &TARGET_H[..7],
                    DT,
                    coupled_newton(),
                    &heat,
                    EnthalpyRadiation::None,
                )
                .is_err()
        );
        let gate = CancelGate::new_clock_free();
        gate.request();
        context(&gate, |cancelled| {
            assert!(matches!(
                engine.advance_with_nodal_heat_and_temperature_contacts(
                    cancelled,
                    problem,
                    &contacts,
                    &old,
                    DT,
                    coupled_newton(),
                    &heat,
                    EnthalpyRadiation::None,
                ),
                Err(EnthalpyError::Conduction(ConductionError::Cancelled { .. }))
            ));
        });
        assert_eq!(old, before);
    });
}
