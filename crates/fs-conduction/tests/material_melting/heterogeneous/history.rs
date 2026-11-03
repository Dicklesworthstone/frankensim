//! Scheduled original-source melting uses the existing physical history owners.
use super::*;
use fs_conduction::bc::history::ThermalBoundaryHistory;
use fs_conduction::bc::history::march::HeatSourceHistory;
use fs_conduction::bc::history::march::enthalpy::radiation::RadiativeEnthalpyHistoryRun;
use fs_conduction::bc::history::march::enthalpy::{
    EnthalpyHistoryBudget, EnthalpyHistoryRun, PhaseEnthalpyHistoryProblem,
};
use fs_conduction::bc::history::march::radiation::RadiationRegionHistory;
use fs_conduction::{
    AmbientRadiationConfig, EMISSIVITY_DIMS, STEFAN_BOLTZMANN_W_M2_K4 as SIGMA,
    SURFACE_EMISSIVITY_PROPERTY, SurfaceEmissivity, TEMPERATURE_DIMS, VOLUMETRIC_SOURCE_DIMS,
};
use fs_scenario::scenario::Environment;
use fs_scenario::signal::{Interp, TimeSignal};

fn boundary<'m>(cx: &fs_exec::Cx<'_>, mesh: &'m ConductionMesh) -> ThermalBoundaryHistory<'m> {
    ThermalBoundaryHistory::new(
        cx,
        ThermalBoundaryBuilder::new(mesh)
            .remainder("surface", ThermalBc::adiabatic())
            .unwrap(),
        &Environment::earth_lab(),
        vec![],
        0.0,
    )
    .unwrap()
}
fn held(times: &[f64], values: &[f64], dims: Dims) -> TimeSignal {
    TimeSignal::Table {
        times: times.to_vec(),
        values: values.to_vec(),
        dims,
        interp: Interp::Hold,
    }
}
fn limits(mesh: &ConductionMesh, steps: usize) -> EnthalpyHistoryBudget {
    EnthalpyHistoryBudget {
        spatial: EnthalpyBudget {
            max_vertices: mesh.vertex_count(),
            max_elements: mesh.element_count(),
        },
        max_steps: steps,
    }
}

#[test]
fn held_heating_and_cooling_cross_original_fusion_in_be_and_sdirk_histories() {
    let mesh = mesh();
    let materials = [model(CoexistenceConductivity::Parallel)];
    with_cx(|cx| {
        let boundary = boundary(cx, &mesh);
        let problem = || PhaseEnthalpyHistoryProblem {
            boundary: &boundary,
            interfaces: None,
            source: HeatSourceHistory::Uniform(held(
                &[0.0, 0.5, 1.5, 2.5, 3.0],
                &[400.0, 400.0, -400.0, -400.0, 0.0],
                VOLUMETRIC_SOURCE_DIMS,
            )),
        };
        let mut be = EnthalpyHistoryRun::from_melting_materials(
            cx,
            problem(),
            &materials,
            &[0],
            &[50.0; 4],
            limits(&mesh, 4),
        )
        .unwrap();
        let mut rk = EnthalpyHistoryRun::from_melting_materials_with_nodal_reference_density(
            cx,
            problem(),
            &materials,
            &[0],
            &[2.0; 4],
            &[50.0; 4],
            limits(&mesh, 4),
        )
        .unwrap()
        .into_sdirk2();
        let mut old = 50.0;
        let mut total = [0.0; 2];
        for (end, expected, q) in [
            (0.5, 150.0, 400.0),
            (1.5, 350.0, 400.0),
            (2.5, 150.0, -400.0),
            (3.0, 50.0, -400.0),
        ] {
            let b = be.advance_event_step(cx, 3.0, 3.0, config()).unwrap();
            let r = rk.advance_event_step(cx, 3.0, 3.0, config()).unwrap();
            assert_eq!(b.end_time_s, end);
            assert_eq!(r.end_time_s, end);
            let joules = (end - b.start_time_s) * q / 6.0;
            near(
                b.report.stored_energy_change_j,
                (expected - old) / 3.0,
                2e-9,
            );
            near(
                r.report.stored_energy_change_j,
                (expected - old) / 3.0,
                2e-9,
            );
            near(r.report.source_j, joules, 2e-9);
            near(b.report.source_w * (end - b.start_time_s), joules, 2e-9);
            total[0] += b.report.stored_energy_change_j;
            total[1] += r.report.stored_energy_change_j;
            for run in [&be, rk.conduction()] {
                for i in 0..4 {
                    near(run.specific_enthalpy_j_kg()[i], expected, 2e-8);
                    let phase = reference(expected, CoexistenceConductivity::Parallel);
                    near(run.temperature()[i], phase.0, 1e-8);
                    near(run.liquid_mass_fraction()[i], phase.1, 1e-10);
                }
                assert!(std::ptr::eq(
                    run.storage().melting_transport_source(0).unwrap(),
                    &materials[0]
                ));
            }
            old = expected;
        }
        for heat in total {
            near(heat, 0.0, 2e-9);
        }
        assert_eq!(be.accepted_steps(), 4);
        assert_eq!(rk.conduction().accepted_steps(), 4);
        assert!(
            EnthalpyHistoryRun::from_melting_materials(
                cx,
                problem(),
                &materials,
                &[0],
                &[800.0; 4],
                limits(&mesh, 4)
            )
            .is_err()
        );
    });
}

fn scaled_balance(h: &[f64], old: &[f64], dt: f64, scale: f64) -> [f64; 8] {
    let mut residual = balance(h, old, dt);
    for i in 0..8 {
        let sum = LOAD[4 * (i / 4)..4 * (i / 4) + 4].iter().sum::<f64>();
        residual[i] += dt * (1.0 - scale) * (sum + LOAD[i]) / 120.0;
    }
    residual
}

#[test]
fn regional_contact_history_preserves_source_laws_at_events_and_after_failed_future_heat() {
    let mesh = coupon();
    let materials = [model(CoexistenceConductivity::Parallel), second_source()];
    with_cx(|cx| {
        let boundary = ThermalBoundaryHistory::new(
            cx,
            ThermalBoundaryBuilder::new(&mesh).adiabatic_remainder(),
            &Environment::earth_lab(),
            vec![],
            0.0,
        )
        .unwrap();
        let initial = boundary.at(cx, 0.0).unwrap();
        let contact = contact(&mesh, &initial);
        let problem = || PhaseEnthalpyHistoryProblem {
            boundary: &boundary,
            interfaces: Some(&contact),
            source: HeatSourceHistory::Scaled {
                field: ScalarField::nodal("heater footprint", 8, LOAD.to_vec()).unwrap(),
                scale: held(
                    &[0.0, 0.01, 0.025, 0.1],
                    &[1.0, -0.5, 2.0, 1e12],
                    Dims::NONE,
                ),
            },
        };
        let create = || {
            EnthalpyHistoryRun::from_melting_materials(
                cx,
                problem(),
                &materials,
                &[0, 1],
                &INITIAL,
                limits(&mesh, 8),
            )
            .unwrap()
        };
        let mut be = create();
        let mut rk = create().into_sdirk2();
        for (end, scale) in [(0.01, 1.0), (0.025, -0.5), (0.04, 2.0)] {
            let old = be.specific_enthalpy_j_kg().to_vec();
            let old_rk = rk.conduction().specific_enthalpy_j_kg().to_vec();
            let b = be.advance_event_step(cx, 0.04, 0.04, policy()).unwrap();
            let r = rk.advance_event_step(cx, 0.04, 0.04, policy()).unwrap();
            assert_eq!(b.end_time_s, end);
            assert_eq!(r.end_time_s, end);
            let dt = end - b.start_time_s;
            for residual in scaled_balance(&b.report.specific_enthalpy_j_kg, &old, dt, scale) {
                near(residual, 0.0, 2e-10);
            }
            let h1 = &r.report.first_stage.specific_enthalpy_j_kg;
            let h2 = &r.report.second_stage.specific_enthalpy_j_kg;
            for residual in scaled_balance(h1, &old_rk, G * dt, scale) {
                near(residual, 0.0, 2e-10);
            }
            let second = scaled_balance(h2, &old_rk, G * dt, scale);
            let carry = scaled_balance(h1, h1, (1.0 - G) * dt, scale);
            for i in 0..8 {
                near(second[i] + carry[i], 0.0, 2e-10);
            }
            let source_j = dt * scale * LOAD.iter().sum::<f64>() / 24.0;
            near(b.report.stored_energy_change_j, source_j, 2e-10);
            near(r.report.stored_energy_change_j, source_j, 2e-10);
            assert!(
                be.storage()
                    .melting_transport_heat_flux(cx, be.specific_enthalpy_j_kg())
                    .unwrap()
                    .iter()
                    .flatten()
                    .any(|q| q.abs() > 1.0)
            );
        }
        let saved = rk.conduction().specific_enthalpy_j_kg().to_vec();
        assert!(rk.advance_to(cx, 0.2, policy()).is_err());
        with_cancelled_cx(|cancelled| assert!(rk.advance_to(cancelled, 0.05, policy()).is_err()));
        assert_eq!(rk.conduction().specific_enthalpy_j_kg(), saved);
        assert_eq!(rk.conduction().time_s(), 0.04);
        assert_eq!(rk.conduction().accepted_steps(), 3);
        let resumed = rk.advance_to(cx, 0.05, policy()).unwrap();
        let mut replay = create().into_sdirk2();
        for _ in 0..3 {
            replay.advance_event_step(cx, 0.04, 0.04, policy()).unwrap();
        }
        assert_eq!(
            resumed.second_stage.specific_enthalpy_j_kg,
            replay
                .advance_to(cx, 0.05, policy())
                .unwrap()
                .second_stage
                .specific_enthalpy_j_kg
        );
    });
}

fn coating() -> SurfaceEmissivity {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS),
            value: PropertyValue::Scalar {
                value: 0.5,
                dims: EMISSIVITY_DIMS,
            },
            validity: ValidityDomain::unconstrained().with("T", 300.0, 450.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic fixed surface coating".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "gray coating".into(),
            phase: "solid".into(),
            process: "fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    SurfaceEmissivity::from_card("surface", &card, 350.0, SelectionPolicy::SingleClaimOnly).unwrap()
}

#[test]
fn radiative_melting_history_counts_source_and_reservoir_heat_once_and_resumes_atomically() {
    let mesh = mesh();
    let materials = [model(CoexistenceConductivity::Series)];
    with_cx(|cx| {
        let boundary = boundary(cx, &mesh);
        let base = EnthalpyHistoryRun::from_melting_materials(
            cx,
            PhaseEnthalpyHistoryProblem {
                boundary: &boundary,
                interfaces: None,
                source: HeatSourceHistory::Uniform(held(
                    &[0.0, 0.001],
                    &[4.0, -6.0],
                    VOLUMETRIC_SOURCE_DIMS,
                )),
            },
            &materials,
            &[0],
            &[200.0; 4],
            limits(&mesh, 4),
        )
        .unwrap();
        let mut run = RadiativeEnthalpyHistoryRun::new(
            cx,
            base,
            vec![RadiationRegionHistory {
                region: "surface".into(),
                emissivity: coating(),
                ambient_temperature: held(
                    &[0.0, 0.001, 0.002],
                    &[300.0, 400.0, -1.0],
                    TEMPERATURE_DIMS,
                ),
            }],
            AmbientRadiationConfig::default(),
        )
        .unwrap();
        let mut expected = [200.0; 4];
        for (to, ambient, q, dt) in [
            (0.001, 300.0_f64, 4.0, 0.001),
            (0.0015, 400.0, -6.0, 0.0005),
        ] {
            if to > 0.001 {
                let saved = run.conduction().specific_enthalpy_j_kg().to_vec();
                assert!(run.advance_to(cx, 0.0025, config()).is_err());
                with_cancelled_cx(|cancelled| {
                    assert!(run.advance_to(cancelled, to, config()).is_err())
                });
                assert_eq!(run.conduction().specific_enthalpy_j_kg(), saved);
                assert_eq!(run.conduction().time_s(), 0.001);
                assert_eq!(run.conduction().accepted_steps(), 1);
            }
            let step = run.advance_event_step(cx, to, to, config()).unwrap();
            let flux = 0.5 * SIGMA * (350.0_f64.powi(4) - ambient.powi(4));
            for i in 0..4 {
                let area: f64 = mesh
                    .boundary()
                    .iter()
                    .filter(|face| face.vertices.contains(&(i as u32)))
                    .map(|face| face.area / 3.0)
                    .sum();
                expected[i] += dt * (q / 24.0 - flux * area) / MASS;
                near(
                    run.conduction().specific_enthalpy_j_kg()[i],
                    expected[i],
                    2e-8,
                );
                near(run.conduction().temperature()[i], 350.0, 0.0);
            }
            let area: f64 = mesh.boundary().iter().map(|face| face.area).sum();
            near(
                step.report.radiation.nonlinear_radiation_out_w,
                flux * area,
                2e-9,
            );
            near(
                step.report.conduction.stored_energy_change_j,
                dt * (q / 6.0 - flux * area),
                2e-9,
            );
            near(
                step.report.physical_energy_residual_j,
                0.0,
                config().energy_tolerance_j,
            );
            assert!(std::ptr::eq(
                run.conduction()
                    .storage()
                    .melting_transport_source(0)
                    .unwrap(),
                &materials[0]
            ));
        }
    });
}
