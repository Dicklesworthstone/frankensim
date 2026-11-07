//! Actual homogeneous-to-separated mixing on the original implicit flow clock.

use super::*;
use fs_material::liquid::solution::multicomponent::caloric::equilibrium::storage::conduction::MulticomponentLiquidEquilibriumConductionLaw;
use fs_material::liquid::solution::multicomponent::caloric::storage::{
    MulticomponentLiquidThermalStorage, conduction::MulticomponentLiquidConductionLaw,
};
use fs_material::thermal::conduction::ConductivityCurve;
use fs_thermal::enthalpy::advection::species::InletSpecies;
use fs_thermal::enthalpy::advection::{Advected, MassBoundary};
use fs_thermal::enthalpy::conduction::implicit::{
    ImplicitConductionChain, ImplicitConductionControls, ImplicitConductionLaw,
    ImplicitConductionWorkspace,
};
use fs_thermal::enthalpy::conduction::{ConductionObserver, ConductivityLaw};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionSimulation, Event,
    MaterialId, StepConfig, ThermalLaw, TimedEvent,
};

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.17e}, expected={expected:.17e}, tolerance={tolerance:.3e}"
    );
}

fn conductivity() -> ConductivityCurve {
    ConductivityCurve::new(vec![[380.0, 5.0], [450.0, 12.0], [490.0, 14.0]]).unwrap()
}

fn temperature_oracle(h: f64) -> f64 {
    let (mut lo, mut hi) = (380.0, 490.0);
    for _ in 0..80 {
        let t = 0.5 * (lo + hi);
        if enthalpy_oracle(t) < h {
            lo = t;
        } else {
            hi = t;
        }
    }
    0.5 * (lo + hi)
}

fn fractions() -> [f64; 3] {
    let masses = core::array::from_fn::<_, 3, _>(|i| AMOUNTS[i] * MASSES[i]);
    let total = masses.iter().sum::<f64>();
    let a = masses[0] / total;
    let b = masses[1] / total;
    [a, b, 1.0 - (a + b)]
}

fn pure_h(index: usize, t: f64) -> f64 {
    let factor = if t <= 450.0 {
        let d = t - 300.0;
        0.6 * d + 0.002 * d * d
    } else {
        let d = t - 450.0;
        135.0 + 1.2 * d + 0.003 * d * d
    };
    MOLAR_CP_AT_400[index] / MASSES[index] * factor
}

fn pure_s(index: usize, t: f64) -> f64 {
    let factor = if t <= 450.0 {
        -0.6 * (t / 300.0).ln() + 0.004 * (t - 300.0)
    } else {
        -0.6 * 1.5_f64.ln() + 0.6 - 1.5 * (t / 450.0).ln() + 0.006 * (t - 450.0)
    };
    MOLAR_CP_AT_400[index] / MASSES[index] * factor
}

fn inlets(
    caloric: &MulticomponentLiquidCaloricLaw<'_, '_>,
) -> Vec<MulticomponentLiquidThermalStorage> {
    (0..3)
        .map(|i| {
            let mut c = controls().thermal;
            c.temperature_bracket[0] = Temperature::new(310.0 + 10.0 * i as f64);
            caloric
                .homogeneous_stream_storage_from_mass_fractions(
                    &core::array::from_fn::<_, 3, _>(|j| if i == j { 1.0 } else { 0.0 }),
                    c,
                    || false,
                )
                .unwrap()
        })
        .collect()
}

#[test]
fn g1_equilibrium_streams_share_original_pure_references_and_actual_phase_entropy() {
    let card = card(1.0);
    let liquid = law(&card, enthalpic_model());
    let caloric = MulticomponentLiquidCaloricLaw::try_new(&liquid, sources(), "p").unwrap();
    let initial = initial(&caloric, 400.0, AMOUNTS);
    let mass = initial.equilibrium().total_mass_kg();
    let original = initial.clone();
    let storage = caloric
        .equilibrium_stream_storage_from_mass_fractions(&initial, &fractions(), controls(), || {
            false
        })
        .unwrap();
    let composition = storage.advected_composition().unwrap();
    let anchor = storage.advected_entropy_anchor().unwrap();
    assert_eq!(composition.mass_fractions, fractions());
    assert_eq!(composition.species, ["A", "B", "C"]);
    assert!(storage.caloric_identity().is_some());
    for inlet in inlets(&caloric) {
        assert_eq!(
            inlet.advected_composition().unwrap().enthalpy_reference,
            composition.enthalpy_reference
        );
        assert_eq!(
            inlet.advected_entropy_anchor().unwrap().entropy_reference,
            anchor.entropy_reference
        );
        assert_ne!(inlet.caloric_identity(), storage.caloric_identity());
    }
    close(anchor.specific_h_j_kg, enthalpy_oracle(380.0) / mass, 2e-5);
    close(
        anchor.reduced_specific_entropy_j_kg_k,
        entropy_oracle(380.0) / mass,
        2e-7,
    );
    for t in [380.0, 400.0, 450.0, 470.0, 490.0] {
        let h = storage
            .specific_enthalpy_j_kg(Temperature::new(t), || false)
            .unwrap();
        close(h, enthalpy_oracle(t) / mass, 2e-5);
        close(storage.temperature_k(h).unwrap(), t, 2e-7);
        let s = anchor.reduced_specific_entropy_j_kg_k
            + storage
                .specific_entropy_change_j_kg_k(anchor.specific_h_j_kg, h)
                .unwrap()
                .unwrap();
        close(s, entropy_oracle(t) / mass, 2e-7);
        let state = storage
            .caloric_state_at_specific_enthalpy(h, || false)
            .unwrap();
        close(
            state.reduced_entropy_j_per_k().unwrap(),
            entropy_oracle(t),
            2e-7,
        );
        assert_eq!(state.candidate().phases().len(), 3);
    }
    let homogeneous = caloric
        .model()
        .homogeneous(Temperature::new(380.0), &AMOUNTS, || false)
        .unwrap();
    assert!(
        (homogeneous.reduced_entropy_j_per_k().unwrap() / mass
            - anchor.reduced_specific_entropy_j_kg_k)
            .abs()
            > 1.0
    );
    for invalid in [
        &[0.2, 0.3, 0.5][..],
        &[0.0, 0.0, 1.0],
        &[f64::NAN, 0.0, 1.0],
        &[1.0, 0.0],
    ] {
        assert!(
            caloric
                .equilibrium_stream_storage_from_mass_fractions(
                    &initial,
                    invalid,
                    controls(),
                    || false
                )
                .is_err()
        );
    }
    let foreign_sources = (0..3)
        .map(|i| (species(i), cp_source(i, 100_000.0, 310.0)))
        .collect();
    let foreign = MulticomponentLiquidCaloricLaw::try_new(&liquid, foreign_sources, "p").unwrap();
    assert!(
        foreign
            .equilibrium_stream_storage_from_mass_fractions(
                &initial,
                &fractions(),
                controls(),
                || false
            )
            .is_err()
    );
    assert!(matches!(
        caloric.equilibrium_stream_storage_from_mass_fractions(
            &initial,
            &fractions(),
            controls(),
            || true
        ),
        Err(MulticomponentLiquidCaloricError::Cancelled)
    ));
    assert_eq!(initial, original);
}

fn flow(inlet_h: [f64; 3]) -> Advected<()> {
    let mut boundaries: Vec<_> = (0..3)
        .map(|i| MassBoundary::Inlet {
            id: i as u64,
            cell: CellId(1),
            kg_s: fractions()[i],
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
                    mass_fractions: (0..3).map(|j| if i == j { 1.0 } else { 0.0 }).collect(),
                })
                .collect(),
            1,
            3,
        )
        .unwrap()
}

#[test]
fn g1_native_be_and_sdirk_mix_homogeneous_inlets_into_changing_equilibrium_phases() {
    let card = card(1.0);
    let liquid = law(&card, enthalpic_model());
    let caloric = MulticomponentLiquidCaloricLaw::try_new(&liquid, sources(), "p").unwrap();
    let initial = initial(&caloric, 400.0, AMOUNTS);
    let mass = initial.equilibrium().total_mass_kg();
    let receiver = caloric
        .equilibrium_stream_storage_from_mass_fractions(&initial, &fractions(), controls(), || {
            false
        })
        .unwrap();
    let inlets = inlets(&caloric);
    let receiver_law =
        MulticomponentLiquidEquilibriumConductionLaw::try_new(&receiver, conductivity()).unwrap();
    let inlet_laws: Vec<_> = inlets
        .iter()
        .map(|s| MulticomponentLiquidConductionLaw::try_new(&card, s).unwrap())
        .collect();
    let mut calorics: Vec<(MaterialId, &dyn ThermalLaw)> = vec![(MaterialId(1), &receiver)];
    let mut conductivities: Vec<(MaterialId, &dyn ConductivityLaw)> =
        vec![(MaterialId(1), &receiver_law)];
    let mut implicit_laws: Vec<(MaterialId, &dyn ImplicitConductionLaw)> =
        vec![(MaterialId(1), &receiver_law)];
    for i in 0..3 {
        calorics.push((MaterialId(i as u64 + 2), &inlets[i]));
        conductivities.push((MaterialId(i as u64 + 2), &inlet_laws[i]));
        implicit_laws.push((MaterialId(i as u64 + 2), &inlet_laws[i]));
    }
    for inlet_t in [[430.0, 470.0, 485.0], [400.0; 3]] {
        let initial_h = receiver
            .specific_enthalpy_j_kg(Temperature::new(400.0), || false)
            .unwrap();
        let inlet_h = core::array::from_fn(|i| {
            inlets[i]
                .specific_enthalpy_j_kg(Temperature::new(inlet_t[i]), || false)
                .unwrap()
        });
        let flow = flow(inlet_h);
        let timeline = DepositionSimulation::new(
            vec![TimedEvent {
                time_s: 0.0,
                event: Event::Birth(DepositedCell {
                    id: CellId(1),
                    segment_id: 1,
                    body_id: 1,
                    material_id: MaterialId(1),
                    mass_kg: 2.0,
                    specific_h_j_kg: initial_h,
                    geometry: CellGeometry {
                        start_m: [0.0; 3],
                        end_m: [1.0, 0.0, 0.0],
                        volume_m3: 0.002,
                        surface_area_m2: 1.0,
                    },
                    boundary: Boundary {
                        ambient_k: 300.0,
                        htc_w_m2_k: 0.0,
                        reservoir: None,
                        source_w: 0.0,
                    },
                }),
            }],
            &calorics,
            DepositionBudget {
                max_cells: 4,
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
        let incoming_h: f64 = (0..3).map(|i| fractions()[i] * pure_h(i, inlet_t[i])).sum();
        let incoming_s: f64 = (0..3).map(|i| fractions()[i] * pure_s(i, inlet_t[i])).sum();
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
            let c = ImplicitConductionControls {
                maximum_dt_s: dt,
                maximum_halvings: 0,
                cell_tolerance_j: 1e-5,
                total_tolerance_j: 2e-5,
                ..Default::default()
            };
            let advance = |system: &ImplicitConductionChain<'_>| {
                let mut work = ImplicitConductionWorkspace::new(1).unwrap();
                if second_order {
                    system.advance_sdirk2_to(&timeline, dt, c, &mut work, || Ok(()))
                } else {
                    system.advance_to(&timeline, dt, c, &mut work, || Ok(()))
                }
            };
            let (next, report) = advance(&audited).unwrap();
            let (plain_next, plain_report) = advance(&plain).unwrap();
            assert_eq!(next.cells(), plain_next.cells());
            assert_eq!(report.energy, plain_report.energy);
            close(next.cells()[&CellId(1)].specific_h_j_kg, final_h, 2e-5);
            let final_t = temperature_oracle(final_h * mass);
            close(next.temperature_k(CellId(1)).unwrap(), final_t, 5e-7);
            if inlet_t == [400.0; 3] {
                assert!(final_t < 399.5);
            } else {
                assert!(final_t > 450.0);
            }
            let stages = if second_order {
                [((1.0 - gamma) * dt, stage_h), (gamma * dt, final_h)]
            } else {
                [(dt, final_h), (0.0, final_h)]
            };
            let (mut boundary_s, mut mixing_s, mut carried_s) = (0.0, 0.0, 0.0);
            for (weight, h) in stages {
                let t = temperature_oracle(h * mass);
                let outgoing_s = entropy_oracle(t) / mass;
                let heat_s = (incoming_h - h) / t;
                boundary_s += weight * (incoming_s - outgoing_s);
                mixing_s += weight * (heat_s - (incoming_s - outgoing_s));
                carried_s += weight * heat_s;
            }
            let entropy = report.entropy;
            let transported = entropy.advection.unwrap();
            close(transported.boundary_entropy_input_j_k, boundary_s, 2e-5);
            close(transported.upwind_mixing_j_k, mixing_s, 2e-5);
            close(transported.balance_residual_j_k, 0.0, 2e-8);
            let stored_s = 2.0 * (entropy_oracle(final_t) - entropy_oracle(400.0)) / mass;
            close(entropy.thermal_storage_change_j_k, stored_s, 2e-6);
            close(entropy.temporal_defect_j_k, stored_s - carried_s, 2e-6);
            assert_eq!(report.energy.substeps, 1);
            assert_eq!(report.rejected_intervals, 0);
            assert!(!entropy.has_unclassified_sources);
        }
        assert_eq!(timeline.time_s(), 0.0);
        assert_eq!(flow.audit().substeps, 0);
    }
}
