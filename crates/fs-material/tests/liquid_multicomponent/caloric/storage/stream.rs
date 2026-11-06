//! Original ternary mixing energy and entropy on native BE and SDIRK stages.

use super::*;
use fs_thermal::enthalpy::advection::species::InletSpecies;
use fs_thermal::enthalpy::advection::{Advected, MassBoundary};

const COMPOSITIONS: [[f64; 3]; 4] = [
    [0.25, 0.25, 0.5],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];
fn coefficients(y: [f64; 3]) -> [f64; 3] {
    let n = core::array::from_fn::<_, 3, _>(|i| y[i] / MASSES[i]);
    let total = n.iter().sum::<f64>();
    let mixing = R * (120.0 * n[0] * n[1] - 60.0 * n[0] * n[2] + 80.0 * n[1] * n[2]) / total;
    [
        mixing,
        80.0 * n[0] + 120.0 * n[1] + 200.0 * n[2],
        0.1 * n[0],
    ]
}
fn expected_h(y: [f64; 3], t: f64) -> f64 {
    let [c, a, b] = coefficients(y);
    let d = t - 300.0;
    c + d * (a + b * d)
}
fn expected_t(y: [f64; 3], h: f64) -> f64 {
    let [c, a, b] = coefficients(y);
    300.0 + 2.0 * (h - c) / (a + (a * a + 4.0 * b * (h - c)).sqrt())
}
fn expected_s(y: [f64; 3], t: f64) -> f64 {
    let n = core::array::from_fn::<_, 3, _>(|i| y[i] / MASSES[i]);
    let total = n.iter().sum::<f64>();
    let ideal = -R
        * n.iter()
            .filter(|&&v| v > 0.0)
            .map(|v| v * (v / total).ln())
            .sum::<f64>();
    let excess = -R * (0.3 * n[0] * n[1] - 0.1 * n[0] * n[2] + 0.6 * n[1] * n[2]) / total;
    let sensible =
        (20.0 * n[0] + 120.0 * n[1] + 200.0 * n[2]) * (t / 300.0).ln() + 0.2 * n[0] * (t - 300.0);
    sensible + ideal + excess
}
fn storages(
    law: &MulticomponentLiquidCaloricLaw<'_, '_>,
) -> Vec<MulticomponentLiquidThermalStorage> {
    COMPOSITIONS
        .iter()
        .enumerate()
        .map(|(i, y)| {
            let mut control = thermal_controls();
            control.temperature_bracket[0] = Temperature::new([310.0, 300.0, 320.0, 330.0][i]);
            law.homogeneous_stream_storage_from_mass_fractions(y, control, || false)
                .unwrap()
        })
        .collect()
}
fn flow(inlet_h: [f64; 3]) -> Advected<()> {
    let rates = COMPOSITIONS[0];
    let mut boundaries: Vec<_> = (0..3)
        .map(|i| MassBoundary::Inlet {
            id: i as u64,
            cell: CellId(1),
            kg_s: rates[i],
            material_id: MaterialId(i as u64 + 2),
            specific_h_j_kg: inlet_h[i],
        })
        .collect();
    boundaries.push(MassBoundary::Outlet {
        id: 3,
        cell: CellId(1),
        kg_s: 1.0,
    });
    Advected::new((), vec![], 0)
        .unwrap()
        .with_boundaries(boundaries, 4)
        .unwrap()
        .with_stationary_species(
            vec!["A".into(), "B".into(), "C".into()],
            (0..3)
                .map(|i| InletSpecies {
                    boundary_id: i as u64,
                    mass_fractions: COMPOSITIONS[i + 1].to_vec(),
                })
                .collect(),
            1,
            3,
        )
        .unwrap()
}

#[test]
fn common_references_keep_original_nonideal_energy_entropy_and_exact_species() {
    let card = thermal_card(true);
    let liquid = liquid(&card);
    let law =
        MulticomponentLiquidCaloricLaw::try_new(&liquid, variable_sources(300.0), "p").unwrap();
    let storage = storages(&law);
    let reference = storage[0]
        .advected_composition()
        .unwrap()
        .enthalpy_reference;
    let entropy_reference = storage[0]
        .advected_entropy_anchor()
        .unwrap()
        .entropy_reference;
    for (i, item) in storage.iter().enumerate() {
        let composition = item.advected_composition().unwrap();
        assert_eq!(composition.species, ["A", "B", "C"]);
        assert_eq!(composition.mass_fractions, COMPOSITIONS[i]);
        assert_eq!(composition.enthalpy_reference, reference);
        let anchor = item.advected_entropy_anchor().unwrap();
        assert_eq!(anchor.entropy_reference, entropy_reference);
        let t0 = item.temperature_domain()[0].value();
        close(
            anchor.specific_h_j_kg,
            expected_h(COMPOSITIONS[i], t0),
            2e-8,
        );
        close(
            anchor.reduced_specific_entropy_j_kg_k,
            expected_s(COMPOSITIONS[i], t0),
            2e-9,
        );
        for t in [350.0, 400.0, 550.0] {
            close(h(item, t), expected_h(COMPOSITIONS[i], t), 2e-8);
            let entropy = anchor.reduced_specific_entropy_j_kg_k
                + item
                    .specific_entropy_change_j_kg_k(anchor.specific_h_j_kg, h(item, t))
                    .unwrap()
                    .unwrap();
            close(entropy, expected_s(COMPOSITIONS[i], t), 2e-8);
        }
        assert!(item.caloric_identity().is_some());
        if i > 0 {
            assert_ne!(item.caloric_identity(), storage[0].caloric_identity());
        }
    }
    let repeat = storages(&law);
    assert_eq!(storage, repeat);
    let mut changed = variable_sources(300.0);
    changed[0].1 = cp_source(0, [80.0, 100.0, 150.0], 100_000.0, 300.0, 1);
    let foreign = MulticomponentLiquidCaloricLaw::try_new(&liquid, changed, "p").unwrap();
    let foreign_storage = storages(&foreign);
    assert_ne!(
        foreign_storage[0]
            .advected_composition()
            .unwrap()
            .enthalpy_reference,
        reference
    );
    assert_ne!(
        foreign_storage[0]
            .advected_entropy_anchor()
            .unwrap()
            .entropy_reference,
        entropy_reference
    );
    for invalid in [
        &[0.25, 0.25, 0.4][..],
        &[1.0, 0.0][..],
        &[f64::NAN, 0.0, 1.0][..],
    ] {
        assert!(
            law.homogeneous_stream_storage_from_mass_fractions(invalid, thermal_controls(), || {
                false
            })
            .is_err()
        );
    }
    assert!(matches!(
        law.homogeneous_stream_storage_from_mass_fractions(
            &COMPOSITIONS[0],
            thermal_controls(),
            || true
        ),
        Err(MulticomponentLiquidCaloricError::Cancelled)
    ));
}

#[test]
fn native_ternary_mixer_closes_original_energy_and_entropy_on_be_and_sdirk_stages() {
    let card = thermal_card(true);
    let liquid = liquid(&card);
    let law =
        MulticomponentLiquidCaloricLaw::try_new(&liquid, variable_sources(300.0), "p").unwrap();
    let storage = storages(&law);
    let adapters: Vec<_> = storage
        .iter()
        .map(|s| MulticomponentLiquidConductionLaw::try_new(&card, s).unwrap())
        .collect();
    let calorics: Vec<_> = storage
        .iter()
        .enumerate()
        .map(|(i, s)| (MaterialId(i as u64 + 1), s as &dyn ThermalLaw))
        .collect();
    let conductivities: Vec<_> = adapters
        .iter()
        .enumerate()
        .map(|(i, k)| (MaterialId(i as u64 + 1), k as &dyn ConductivityLaw))
        .collect();
    let implicit_laws: Vec<_> = adapters
        .iter()
        .enumerate()
        .map(|(i, k)| (MaterialId(i as u64 + 1), k as &dyn ImplicitConductionLaw))
        .collect();
    for (initial_t, inlet_t) in [(380.0, [350.0, 450.0, 500.0]), (400.0, [400.0; 3])] {
        let inlet_h = core::array::from_fn(|i| h(&storage[i + 1], inlet_t[i]));
        let flow = flow(inlet_h);
        let timeline = DepositionSimulation::new(
            vec![birth(1, 2.0, h(&storage[0], initial_t))],
            &calorics,
            DepositionBudget {
                max_cells: storage.len(),
                max_links: 1,
                max_events: 1,
                max_substeps: 8,
            },
        )
        .unwrap()
        .advance_to(
            0.0,
            StepConfig {
                max_dt_s: 1.0,
                energy_tolerance_j: 1e-6,
            },
            || Ok(()),
        )
        .unwrap()
        .0;
        let conduction =
            ConductionObserver::new(&timeline, &conductivities, &[], &[], || Ok(())).unwrap();
        let plain = ImplicitConductionChain::new(&timeline, &conduction, &implicit_laws, || Ok(()))
            .unwrap()
            .with_implicit_advection(&timeline, &flow, || Ok(()))
            .unwrap();
        let audited =
            ImplicitConductionChain::new(&timeline, &conduction, &implicit_laws, || Ok(()))
                .unwrap()
                .with_implicit_advection(&timeline, &flow, || Ok(()))
                .unwrap()
                .with_implicit_advection_entropy(&timeline, || Ok(()))
                .unwrap();
        let incoming_h: f64 = (0..3).map(|i| COMPOSITIONS[0][i] * inlet_h[i]).sum();
        let incoming_s: f64 = (0..3)
            .map(|i| COMPOSITIONS[0][i] * expected_s(COMPOSITIONS[i + 1], inlet_t[i]))
            .sum();
        let initial_h = expected_h(COMPOSITIONS[0], initial_t);
        let original = timeline.cells().clone();
        for second_order in [false, true] {
            let dt = 20.0;
            let gamma = 1.0 - core::f64::consts::FRAC_1_SQRT_2;
            let stage_h = (2.0 * initial_h + gamma * dt * incoming_h) / (2.0 + gamma * dt);
            let final_h = if second_order {
                (2.0 * initial_h
                    + (1.0 - gamma) * dt * (incoming_h - stage_h)
                    + gamma * dt * incoming_h)
                    / (2.0 + gamma * dt)
            } else {
                (2.0 * initial_h + dt * incoming_h) / (2.0 + dt)
            };
            let controls = ImplicitConductionControls {
                maximum_dt_s: dt,
                maximum_halvings: 0,
                cell_tolerance_j: 1e-5,
                total_tolerance_j: 2e-5,
                ..ImplicitConductionControls::default()
            };
            let advance = |system: &ImplicitConductionChain<'_>| {
                let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
                if second_order {
                    system.advance_sdirk2_to(&timeline, dt, controls, &mut workspace, || Ok(()))
                } else {
                    system.advance_to(&timeline, dt, controls, &mut workspace, || Ok(()))
                }
            };
            let (next, report) = advance(&audited).unwrap();
            let (unclassified, unclassified_report) = advance(&plain).unwrap();
            assert_eq!(next.cells(), unclassified.cells());
            assert_eq!(report.energy, unclassified_report.energy);
            close(next.cells()[&CellId(1)].specific_h_j_kg, final_h, 2e-6);
            let final_t = expected_t(COMPOSITIONS[0], final_h);
            close(next.temperature_k(CellId(1)).unwrap(), final_t, 2e-8);
            if inlet_t == [400.0; 3] {
                assert!(final_t < 399.5);
            }
            let stages = if second_order {
                [((1.0 - gamma) * dt, stage_h), (gamma * dt, final_h)]
            } else {
                [(dt, final_h), (0.0, final_h)]
            };
            let (mut boundary_s, mut mixing_s, mut carried_heat_s) = (0.0, 0.0, 0.0);
            for (weight, energy) in stages {
                let t = expected_t(COMPOSITIONS[0], energy);
                let outgoing_s = expected_s(COMPOSITIONS[0], t);
                let heat_s = (incoming_h - energy) / t;
                boundary_s += weight * (incoming_s - outgoing_s);
                mixing_s += weight * (heat_s - (incoming_s - outgoing_s));
                carried_heat_s += weight * heat_s;
            }
            let entropy = report.entropy;
            let flow_entropy = entropy.advection.unwrap();
            close(flow_entropy.boundary_entropy_input_j_k, boundary_s, 2e-7);
            close(flow_entropy.upwind_mixing_j_k, mixing_s, 2e-7);
            close(flow_entropy.balance_residual_j_k, 0.0, 2e-9);
            let storage_s = 2.0
                * (expected_s(COMPOSITIONS[0], final_t) - expected_s(COMPOSITIONS[0], initial_t));
            close(entropy.thermal_storage_change_j_k, storage_s, 2e-7);
            close(
                entropy.temporal_defect_j_k,
                storage_s - carried_heat_s,
                2e-7,
            );
            assert!(flow_entropy.upwind_mixing_j_k > 0.0);
            assert_eq!(entropy.conduction_generation_j_k, 0.0);
            assert!(!entropy.has_unclassified_sources);
            assert_eq!(report.energy.substeps, 1);
            assert_eq!(report.rejected_intervals, 0);
        }
        assert_eq!(timeline.cells(), &original);
        assert_eq!(timeline.time_s(), 0.0);
        assert_eq!(flow.audit().substeps, 0);
    }
}
