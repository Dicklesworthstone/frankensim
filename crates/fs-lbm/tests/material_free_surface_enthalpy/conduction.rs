//! Source conductivity reaches the original finite-volume heat and VOF owners.

use super::*;
use fs_lbm::material::free_surface::enthalpy::MaterialFreeSurfaceEnthalpy3;
use fs_lbm::material::free_surface::enthalpy::conduction::{
    FreeSurfaceConductionControls, FreeSurfaceConductionFace,
};

#[path = "conduction/body.rs"]
mod body;

#[path = "conduction/environment.rs"]
mod environment;

const PAIR: [[usize; 3]; 2] = [[3, 3, 3], [4, 3, 3]];

fn conducting_card(scale: f64) -> MaterialCard {
    let base = card(0.02, 1000.0);
    let mut claims = ClaimSet::new();
    for (_, claim) in base.claims().claims_ordered() {
        let mut claim = claim.clone();
        if claim.key.name() == "thermal_conductivity" {
            // k(T) = scale * (T - 200) over the original source's support.
            claim.value = PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(300.0, scale * 100.0), (500.0, scale * 300.0)],
                dims: Dims([1, 1, -3, -1, 0, 0]),
            };
        }
        claims.insert_claim(claim).unwrap();
    }
    MaterialCard::assemble(base.id().clone(), claims, base.models().to_vec()).unwrap()
}

fn face() -> FreeSurfaceConductionFace {
    FreeSurfaceConductionFace {
        cells: PAIR,
        area_m2: 1e-6,
        distances_m: [0.0005; 2],
    }
}

fn controls() -> FreeSurfaceConductionControls {
    FreeSurfaceConductionControls {
        maximum_faces: 8,
        maximum_substeps: 1000,
        maximum_dt_s: 0.00025,
        energy_tolerance_j: 1e-9,
    }
}

fn with_liquid<T>(
    card: &MaterialCard,
    moving: bool,
    f: impl FnOnce(&mut MaterialFreeSurfaceEnthalpy3<'_, '_>, &CancelGate) -> T,
) -> T {
    let law = LiquidTemperatureLaw::try_new(
        card,
        &point(400.0),
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let source = resolve_liquid_surface_state(
        card,
        &point(400.0),
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Neglected,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let chart = caloric(card);
    let gate = CancelGate::new_clock_free();
    let mut cells = initial(true);
    for cell in &mut cells {
        cell.velocity_m_s = if moving {
            [0.0, 0.003 * (cell.coordinates[0] as f64 - 3.5), 0.0]
        } else {
            [0.0; 3]
        };
    }
    with_cx(&gate, |cx| {
        let flow = plan_free_surface_from_material(&source, 0.001, 0.001, 0.1)
            .unwrap()
            .free_surface3_from_cells(domain(), &cells, cx)
            .unwrap();
        let mut thermal = flow
            .with_enthalpy(
                &law,
                &chart,
                FreeSurfaceEnthalpyControls::default(),
                |p| {
                    if p == PAIR[0] {
                        100000.0
                    } else if p == PAIR[1] {
                        300000.0
                    } else {
                        200000.0
                    }
                },
                &gate,
            )
            .unwrap();
        f(&mut thermal, &gate)
    })
}

fn temperature(state: &fs_material::liquid::ResolvedLiquidState) -> f64 {
    state
        .material()
        .query_point()
        .iter()
        .find(|(axis, _)| axis == "T")
        .unwrap()
        .1
}

#[test]
fn g1_two_cell_nonlinear_source_conduction_matches_independent_heat_balance() {
    let card = conducting_card(1.0);
    with_liquid(&card, false, |thermal, gate| {
        let energy = thermal.total_enthalpy_j().unwrap();
        let report = thermal
            .step_with_conduction(&[face()], &[], controls(), gate)
            .unwrap();

        // Independent two-cell explicit reference: each cell contains 1 mg,
        // cp = 2000 J/(kg K); heat is equal and opposite at the shared face.
        // Do not reuse the production source evaluation or conductance helper.
        let mut t = [350.0, 450.0];
        let mut applied_g = 0.0;
        for _ in 0..4 {
            let k = t.map(|v| v - 200.0);
            applied_g = 1e-6 / (0.0005 / k[0] + 0.0005 / k[1]);
            let heat = 0.00025 * applied_g * (t[1] - t[0]);
            t[0] += heat / 0.002;
            t[1] -= heat / 0.002;
        }
        assert_eq!(report.thermal.substeps, 4);
        assert_eq!(report.thermal.external_heat_j, 0.0);
        assert_eq!(report.flow.external_heat_j, 0.0);
        near(report.last_applied[0].conductance_w_k, applied_g, 2e-13);
        let endpoint = &report.thermal_endpoint[0];
        for i in 0..2 {
            let source = &endpoint.materials[i];
            assert_eq!(source.material().card_identity(), card.content_hash());
            near(temperature(source), t[i], 2e-13);
            near(source.thermal_conductivity_w_m_k(), t[i] - 200.0, 2e-13);
            near(
                thermal.specific_enthalpy_at(PAIR[i]).unwrap(),
                (t[i] - 300.0) * 2000.0,
                2e-12,
            );
        }
        let endpoint_g = 1e-6 / (0.0005 / (t[0] - 200.0) + 0.0005 / (t[1] - 200.0));
        near(endpoint.conductance_w_k, endpoint_g, 2e-13);
        assert_ne!(
            endpoint.conductance_w_k,
            report.last_applied[0].conductance_w_k
        );
        near(thermal.total_enthalpy_j().unwrap(), energy, 2e-12);
        near(report.flow.mass_before_kg, report.flow.mass_after_kg, 2e-13);
        assert!(report.thermal.residual_j.abs() < 1e-9);
        assert!(report.flow.balance_defect_j.abs() < 1e-9);
    });
}

#[test]
fn g3_conductivity_only_card_swap_changes_heat_viscosity_and_native_momentum() {
    let run = |scale| {
        let card = conducting_card(scale);
        with_liquid(&card, true, |thermal, gate| {
            let initial_energy = thermal.total_enthalpy_j().unwrap();
            for _ in 0..3 {
                thermal
                    .step_with_conduction(&[face()], &[], controls(), gate)
                    .unwrap();
            }
            near(thermal.total_enthalpy_j().unwrap(), initial_energy, 2e-11);
            (
                thermal.flow().lattice().grid().state_bits(),
                thermal
                    .flow()
                    .sample_at_si(3, 3, 3)
                    .unwrap()
                    .unwrap()
                    .velocity_m_s,
                temperature(thermal.material_at(PAIR[0]).unwrap()),
                thermal
                    .material_at(PAIR[0])
                    .unwrap()
                    .dynamic_viscosity_pa_s(),
                thermal.total_enthalpy_j().unwrap(),
            )
        })
    };
    let slow = run(1.0);
    let fast = run(3.0);
    assert_ne!(
        slow.0, fast.0,
        "native populations must change, not just receipts"
    );
    assert!((slow.1[1] - fast.1[1]).abs() > 1e-10);
    assert!(fast.2 > slow.2, "higher k heats the cold cell faster");
    assert!(fast.3 < slow.3, "the original same-card mu(T) must refresh");
    near(slow.4, fast.4, 2e-11);
}

#[test]
fn g4_thermal_work_source_geometry_and_cancellation_failures_are_atomic_and_retryable() {
    let card = conducting_card(1.0);
    let clean = with_liquid(&card, false, |thermal, gate| {
        thermal
            .step_with_conduction(&[face()], &[], controls(), gate)
            .unwrap();
        (
            thermal.flow().lattice().grid().state_bits(),
            thermal.specific_enthalpy_at(PAIR[0]),
            thermal.total_enthalpy_j().unwrap(),
        )
    });
    with_liquid(&card, false, |thermal, gate| {
        let bits = thermal.flow().lattice().grid().state_bits();
        let energy = thermal.total_enthalpy_j().unwrap();
        let material = thermal.material_at(PAIR[0]).unwrap().clone();
        let unchanged = |thermal: &MaterialFreeSurfaceEnthalpy3<'_, '_>| {
            assert_eq!(thermal.flow().lattice().grid().state_bits(), bits);
            assert_eq!(thermal.total_enthalpy_j().unwrap(), energy);
            assert_eq!(thermal.specific_enthalpy_at(PAIR[0]), Some(100000.0));
            assert_eq!(thermal.material_at(PAIR[0]).unwrap(), &material);
        };

        let reverse = FreeSurfaceConductionFace {
            cells: [PAIR[1], PAIR[0]],
            ..face()
        };
        assert!(
            thermal
                .step_with_conduction(&[face(), reverse], &[], controls(), gate)
                .is_err()
        );
        unchanged(thermal);
        let oversized = FreeSurfaceConductionFace {
            area_m2: 1.1e-6,
            ..face()
        };
        assert!(
            thermal
                .step_with_conduction(&[oversized], &[], controls(), gate)
                .is_err()
        );
        unchanged(thermal);

        // One real thermal substep happens before the original owner refuses
        // the remaining work. Neither accepted heat nor flow may be published.
        let small_budget = FreeSurfaceConductionControls {
            maximum_substeps: 1,
            ..controls()
        };
        assert!(matches!(
            thermal.step_with_conduction(&[face()], &[], small_budget, gate),
            Err(FreeSurfaceEnthalpyError::Thermal(_))
        ));
        unchanged(thermal);

        // An isolated heated cell is covered by source admission too. T=550 K
        // is supported by the caloric chart, but not the original source card.
        let heat = [FreeSurfaceCellHeat {
            coordinates: [6, 6, 6],
            energy_j: 0.3,
        }];
        assert!(matches!(
            thermal.step_with_conduction(&[face()], &heat, controls(), gate),
            Err(FreeSurfaceEnthalpyError::Material(_))
        ));
        unchanged(thermal);
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        assert!(
            thermal
                .step_with_conduction(&[face()], &[], controls(), &cancelled)
                .is_err()
        );
        unchanged(thermal);

        thermal
            .step_with_conduction(&[face()], &[], controls(), gate)
            .unwrap();
        assert_eq!(thermal.flow().lattice().grid().state_bits(), clean.0);
        assert_eq!(thermal.specific_enthalpy_at(PAIR[0]), clean.1);
        assert_eq!(thermal.total_enthalpy_j().unwrap(), clean.2);
    });
}

#[test]
fn g1_external_cell_heat_is_counted_once_beside_internal_conduction() {
    let card = conducting_card(1.0);
    with_liquid(&card, false, |thermal, gate| {
        let energy = thermal.total_enthalpy_j().unwrap();
        let heat = [FreeSurfaceCellHeat {
            coordinates: [6, 6, 6],
            energy_j: 0.02,
        }];
        let report = thermal
            .step_with_conduction(&[face()], &heat, controls(), gate)
            .unwrap();
        near(report.thermal.external_heat_j, 0.02, 1e-14);
        near(report.flow.external_heat_j, 0.02, 1e-14);
        near(thermal.total_enthalpy_j().unwrap() - energy, 0.02, 1e-9);
        near(
            temperature(thermal.material_at([6, 6, 6]).unwrap()),
            410.0,
            2e-13,
        );
        assert!(report.flow.balance_defect_j.abs() < 1e-9);
    });
}

#[test]
fn g4_late_native_mass_route_refusal_discards_the_completed_conduction_stage() {
    let card = conducting_card(1.0);
    let law = LiquidTemperatureLaw::try_new(
        &card,
        &point(400.0),
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let source = resolve_liquid_surface_state(
        &card,
        &point(400.0),
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Neglected,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let chart = caloric(&card);
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let flow = plan_free_surface_from_material(&source, 0.001, 0.001, 0.1)
            .unwrap()
            .free_surface3_from_cells(domain(), &initial(false), cx)
            .unwrap();
        let mut thermal = flow
            .with_enthalpy(
                &law,
                &chart,
                FreeSurfaceEnthalpyControls {
                    max_mass_routes: 1,
                    ..FreeSurfaceEnthalpyControls::default()
                },
                |p| if p == [2, 3, 3] { 100000.0 } else { 300000.0 },
                &gate,
            )
            .unwrap();
        let bits = thermal.flow().lattice().grid().state_bits();
        let energy = thermal.total_enthalpy_j().unwrap();
        let material = thermal.material_at([2, 3, 3]).unwrap().clone();
        let face = FreeSurfaceConductionFace {
            cells: [[2, 3, 3], [3, 3, 3]],
            ..face()
        };
        let controls = FreeSurfaceConductionControls {
            maximum_faces: 1,
            ..controls()
        };
        assert!(matches!(
            thermal.step_with_conduction(&[face], &[], controls, &gate),
            Err(FreeSurfaceEnthalpyError::Flow(_))
        ));
        assert_eq!(thermal.flow().lattice().grid().state_bits(), bits);
        assert_eq!(thermal.total_enthalpy_j().unwrap(), energy);
        assert_eq!(thermal.specific_enthalpy_at([2, 3, 3]), Some(100000.0));
        assert_eq!(thermal.material_at([2, 3, 3]).unwrap(), &material);
        assert_eq!(thermal.flow().lattice().conversions().to_fluid, 0);
    });
}
