//! Original nonideal calorics carry heat of mixing through native flow graphs.

use super::*;
use fs_material::liquid::solution::caloric::storage::BinaryLiquidThermalStorage;
use fs_material::liquid::solution::caloric::storage::conduction::BinaryLiquidConductionLaw;
use fs_material::thermal::conduction::ConductivityCurve;
use fs_thermal::enthalpy::advection::species::InletSpecies;
use fs_thermal::enthalpy::advection::{Advected, MassBoundary, MassFlow};
use fs_thermal::enthalpy::conduction::implicit::{
    ImplicitConductionChain, ImplicitConductionControls, ImplicitConductionGraph,
    ImplicitConductionLaw, ImplicitConductionWorkspace,
};
use fs_thermal::enthalpy::conduction::{ConductionObserver, ConductiveLink, ConductivityLaw};

#[path = "stream/entropy.rs"]
mod entropy;

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.17e}, expected={expected:.17e}, tolerance={tolerance:.3e}"
    );
}

fn variable_sources(datum: f64, revision: u32) -> [(SpeciesId, IsobaricHeatCapacityLaw); 2] {
    [
        (
            SpeciesId::new("A").unwrap(),
            cp_source("A", 0.018, [89.0, 100.0, 120.0], 100_000.0, datum, revision),
        ),
        (
            SpeciesId::new("B").unwrap(),
            cp_source("B", 0.10, [138.0, 160.0, 200.0], 100_000.0, datum, revision),
        ),
    ]
}

// Independent homogeneous calorics, per kg: Cp_A=60+.1T, Cp_B=80+.2T,
// pure datum 300 K and H_excess=x_A*x_B*(R*100+C*(T-400)) per mole.
// Both signed C fixtures remain homogeneous over the complete source domain.
fn coefficients(y: [f64; 2], excess_cp: f64) -> [f64; 3] {
    let n = [y[0] / 0.018, y[1] / 0.10];
    let mixed = n[0] * n[1] / (n[0] + n[1]);
    [
        mixed * (R * 100.0 - 100.0 * excess_cp),
        90.0 * n[0] + 140.0 * n[1] + excess_cp * mixed,
        0.05 * n[0] + 0.10 * n[1],
    ]
}
fn expected_h(y: [f64; 2], excess_cp: f64, t: f64) -> f64 {
    let [d, p, q] = coefficients(y, excess_cp);
    let dt = t - 300.0;
    d + dt * (p + q * dt)
}
fn expected_t(y: [f64; 2], excess_cp: f64, h: f64) -> f64 {
    let [d, p, q] = coefficients(y, excess_cp);
    300.0 + 2.0 * (h - d) / (p + (p * p + 4.0 * q * (h - d)).sqrt())
}
fn h(law: &BinaryLiquidThermalStorage, t: f64) -> f64 {
    law.specific_enthalpy_j_kg(Temperature::new(t), || false)
        .unwrap()
}
fn implicit_controls(dt: f64) -> ImplicitConductionControls {
    ImplicitConductionControls {
        maximum_dt_s: dt,
        maximum_halvings: 0,
        cell_tolerance_j: 1e-5,
        total_tolerance_j: 2e-5,
        ..Default::default()
    }
}
fn conduction_laws(storage: &[BinaryLiquidThermalStorage]) -> Vec<BinaryLiquidConductionLaw<'_>> {
    storage
        .iter()
        .map(|law| {
            BinaryLiquidConductionLaw::try_new(
                law,
                ConductivityCurve::new(vec![[290.0, 2000.0], [600.0, 2000.0]]).unwrap(),
            )
            .unwrap()
        })
        .collect()
}
fn bindings<'a>(
    laws: &'a [BinaryLiquidConductionLaw<'_>],
) -> Vec<(MaterialId, &'a dyn ImplicitConductionLaw)> {
    laws.iter()
        .enumerate()
        .map(|(i, law)| (MaterialId(i as u64), law as &dyn ImplicitConductionLaw))
        .collect()
}
fn prepare<'a>(
    laws: &'a [BinaryLiquidConductionLaw<'_>],
    temperatures: &[f64],
    masses: &[f64],
    conducting: bool,
) -> (DepositionSimulation<'a>, ConductionObserver<'a>) {
    let mut events: Vec<_> = temperatures
        .iter()
        .enumerate()
        .map(|(i, &t)| {
            let mut event = birth(i as u64, masses[i], h(laws[i].caloric(), t), 0.0);
            if let Event::Birth(cell) = &mut event.event {
                cell.material_id = MaterialId(i as u64);
            }
            event
        })
        .collect();
    let links = if conducting {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: 0,
                a: CellId(0),
                b: CellId(1),
                area_m2: 0.1,
                conductance_w_k: 0.0,
            }),
        });
        vec![ConductiveLink {
            link_id: 0,
            half_distances_m: [0.1; 2],
            contact_resistance_m2_k_w: 0.0,
        }]
    } else {
        vec![]
    };
    let calorics: Vec<_> = laws
        .iter()
        .enumerate()
        .map(|(i, law)| (MaterialId(i as u64), law.caloric() as &dyn ThermalLaw))
        .collect();
    let conductivities: Vec<_> = laws
        .iter()
        .enumerate()
        .map(|(i, law)| (MaterialId(i as u64), law as &dyn ConductivityLaw))
        .collect();
    let simulation = DepositionSimulation::new(
        events,
        &calorics,
        DepositionBudget {
            max_cells: laws.len(),
            max_links: 4,
            ..budget()
        },
    )
    .unwrap();
    let thermal = simulation
        .advance_to(0.0, config(1.0), || Ok(()))
        .unwrap()
        .0;
    let conduction =
        ConductionObserver::new(&thermal, &conductivities, &links, &[], || Ok(())).unwrap();
    (thermal, conduction)
}
fn mixer_flow(inlet_h: [f64; 2]) -> Advected<()> {
    Advected::new((), vec![], 0)
        .unwrap()
        .with_boundaries(
            vec![
                MassBoundary::Inlet {
                    id: 10,
                    cell: CellId(0),
                    kg_s: 0.25,
                    material_id: MaterialId(1),
                    specific_h_j_kg: inlet_h[0],
                },
                MassBoundary::Inlet {
                    id: 11,
                    cell: CellId(0),
                    kg_s: 0.75,
                    material_id: MaterialId(2),
                    specific_h_j_kg: inlet_h[1],
                },
                MassBoundary::Outlet {
                    id: 12,
                    cell: CellId(0),
                    kg_s: 1.0,
                },
            ],
            3,
        )
        .unwrap()
        .with_stationary_species(
            vec!["A".into(), "B".into()],
            vec![
                InletSpecies {
                    boundary_id: 10,
                    mass_fractions: vec![1.0, 0.0],
                },
                InletSpecies {
                    boundary_id: 11,
                    mass_fractions: vec![0.0, 1.0],
                },
            ],
            1,
            2,
        )
        .unwrap()
}

#[test]
fn native_open_mixer_preserves_heat_of_mixing_and_variable_signed_specific_heat() {
    for excess_cp in [-40.0, 40.0] {
        let card = card(CompositionBasis::MoleFraction, 1.0);
        let liquid = liquid_with_excess_cp(&card, 100.0, excess_cp);
        let caloric =
            BinaryLiquidCaloricLaw::try_new(&liquid, variable_sources(300.0, 0), "pressure")
                .unwrap();
        let storage: Vec<_> = [[0.25, 0.75], [1.0, 0.0], [0.0, 1.0]]
            .map(|fractions| {
                caloric
                    .stream_storage_from_mass_fractions(fractions, controls(), || false)
                    .unwrap()
            })
            .into();
        let laws = conduction_laws(&storage);
        for (initial_t, inlet_t) in [(380.0, [350.0, 500.0]), (400.0, [400.0; 2])] {
            let inlet_h = [h(&storage[1], inlet_t[0]), h(&storage[2], inlet_t[1])];
            let flow = mixer_flow(inlet_h);
            let (thermal, conduction) = prepare(&laws, &[initial_t], &[2.0], false);
            let system =
                ImplicitConductionChain::new(&thermal, &conduction, &bindings(&laws), || Ok(()))
                    .unwrap()
                    .with_implicit_advection(&thermal, &flow, || Ok(()))
                    .unwrap();
            let dt = 20.0; // Ten complete cell residence times in one accepted step.
            let h_initial = expected_h([0.25, 0.75], excess_cp, initial_t);
            let incoming = 0.25 * expected_h([1.0, 0.0], excess_cp, inlet_t[0])
                + 0.75 * expected_h([0.0, 1.0], excess_cp, inlet_t[1]);
            let expected = (2.0 * h_initial + dt * incoming) / (2.0 + dt);
            let expected_temperature = expected_t([0.25, 0.75], excess_cp, expected);
            let (next, report) = system
                .advance_to(
                    &thermal,
                    dt,
                    implicit_controls(dt),
                    &mut ImplicitConductionWorkspace::new(1).unwrap(),
                    || Ok(()),
                )
                .unwrap();
            near(next.cells()[&CellId(0)].specific_h_j_kg, expected, 2e-6);
            near(
                next.temperature_k(CellId(0)).unwrap(),
                expected_temperature,
                2e-8,
            );
            near(
                report.energy.external_heat_j,
                2.0 * (expected - h_initial),
                2e-5,
            );
            near(report.energy.residual_j, 0.0, 2e-5);
            assert_eq!(report.energy.substeps, 1);
            assert_eq!(report.rejected_intervals, 0);
            assert_eq!(next.mass_kg(), thermal.mass_kg());
            assert_eq!(thermal.time_s(), 0.0);
            assert_eq!(
                flow.stationary_species().unwrap().cells()[&CellId(0)],
                [0.25, 0.75]
            );
            assert!(report.entropy.advection.is_none());
            assert!(report.entropy.has_unclassified_sources);
            if inlet_t == [400.0; 2] {
                // Equal feed/cell temperatures still cool because the positive
                // original mixing enthalpy consumes sensible heat.
                assert!(expected_temperature < 399.0);
            }
            assert!(
                system
                    .with_implicit_advection_entropy(&thermal, || Ok(()))
                    .is_ok()
            );
        }
    }
}

#[test]
fn native_recirculation_couples_composition_specific_cp_with_conduction_conservatively() {
    let excess_cp = -40.0;
    let card = card(CompositionBasis::MoleFraction, 1.0);
    let liquid = liquid_with_excess_cp(&card, 100.0, excess_cp);
    let caloric =
        BinaryLiquidCaloricLaw::try_new(&liquid, variable_sources(300.0, 0), "pressure").unwrap();
    let a = caloric
        .stream_storage_from_mass_fractions([1.0, 0.0], controls(), || false)
        .unwrap();
    let b = caloric
        .stream_storage_from_mass_fractions([0.0, 1.0], controls(), || false)
        .unwrap();
    let inlet = [h(&a, 460.0), h(&b, 330.0)];
    let flow = Advected::new(
        (),
        vec![
            MassFlow {
                id: 10,
                a: CellId(0),
                b: CellId(1),
                kg_s: 2.0,
            },
            MassFlow {
                id: 11,
                a: CellId(1),
                b: CellId(0),
                kg_s: 1.0,
            },
        ],
        2,
    )
    .unwrap()
    .with_boundaries(
        vec![
            MassBoundary::Inlet {
                id: 20,
                cell: CellId(0),
                kg_s: 1.0,
                material_id: MaterialId(2),
                specific_h_j_kg: inlet[0],
            },
            MassBoundary::Inlet {
                id: 21,
                cell: CellId(1),
                kg_s: 1.0,
                material_id: MaterialId(3),
                specific_h_j_kg: inlet[1],
            },
            MassBoundary::Outlet {
                id: 22,
                cell: CellId(1),
                kg_s: 2.0,
            },
        ],
        5,
    )
    .unwrap()
    .with_stationary_species(
        vec!["A".into(), "B".into()],
        vec![
            InletSpecies {
                boundary_id: 20,
                mass_fractions: vec![1.0, 0.0],
            },
            InletSpecies {
                boundary_id: 21,
                mass_fractions: vec![0.0, 1.0],
            },
        ],
        2,
        2,
    )
    .unwrap();
    let species = flow.stationary_species().unwrap();
    let fractions: [[f64; 2]; 2] =
        [0, 1].map(|i| species.cells()[&CellId(i)].as_slice().try_into().unwrap());
    near(fractions[0][0], 0.75, 1e-14);
    near(fractions[1][0], 0.5, 1e-14);
    assert!(species.maximum_residual_kg_s() < 1e-14);
    let storage = vec![
        caloric
            .stream_storage_from_mass_fractions(fractions[0], controls(), || false)
            .unwrap(),
        caloric
            .stream_storage_from_mass_fractions(fractions[1], controls(), || false)
            .unwrap(),
        a,
        b,
    ];
    let laws = conduction_laws(&storage);
    let masses = [2.0, 3.0];
    let initial_t = [360.0, 420.0];
    let old = [0, 1].map(|i| expected_h(fractions[i], excess_cp, initial_t[i]));
    let (thermal, conduction) = prepare(&laws, &initial_t, &masses, true);
    let system = ImplicitConductionGraph::new(
        &thermal,
        &conduction,
        &bindings(&laws),
        Default::default(),
        || Ok(()),
    )
    .unwrap()
    .with_implicit_advection(&thermal, &flow, || Ok(()))
    .unwrap();
    let dt = 4.0;
    // Independent two-cell BE balance at a trial conductive power into cell 0.
    // Mass fluxes carry common-datum h; G=area/(2*distance/k)=1000 W/K.
    let solve_at_power = |power| {
        let a = masses[0] + 2.0 * dt;
        let b = -dt;
        let c = -2.0 * dt;
        let d = masses[1] + 3.0 * dt;
        let rhs = [
            masses[0] * old[0] + dt * (inlet[0] + power),
            masses[1] * old[1] + dt * (inlet[1] - power),
        ];
        [
            (d * rhs[0] - b * rhs[1]) / (a * d - b * c),
            (a * rhs[1] - c * rhs[0]) / (a * d - b * c),
        ]
    };
    let residual = |power| {
        let energy = solve_at_power(power);
        power
            - 1000.0
                * (expected_t(fractions[1], excess_cp, energy[1])
                    - expected_t(fractions[0], excess_cp, energy[0]))
    };
    let mut bracket = [-310_000.0, 310_000.0];
    assert!(residual(bracket[0]) < 0.0 && residual(bracket[1]) > 0.0);
    for _ in 0..80 {
        let midpoint = 0.5 * (bracket[0] + bracket[1]);
        if residual(midpoint) < 0.0 {
            bracket[0] = midpoint;
        } else {
            bracket[1] = midpoint;
        }
    }
    let expected = solve_at_power(0.5 * (bracket[0] + bracket[1]));
    let (next, report) = system
        .advance_to(
            &thermal,
            dt,
            implicit_controls(dt),
            &mut ImplicitConductionWorkspace::new(2).unwrap(),
            || Ok(()),
        )
        .unwrap();
    for i in 0..2 {
        near(
            next.cells()[&CellId(i as u64)].specific_h_j_kg,
            expected[i],
            3e-6,
        );
        near(
            next.temperature_k(CellId(i as u64)).unwrap(),
            expected_t(fractions[i], excess_cp, expected[i]),
            2e-8,
        );
        assert_eq!(
            storage[i].advected_composition().unwrap().mass_fractions,
            fractions[i]
        );
    }
    let external = dt * (inlet[0] + inlet[1] - 2.0 * expected[1]);
    near(report.energy.external_heat_j, external, 2e-5);
    near(report.energy.stored_energy_change_j, external, 2e-5);
    near(report.energy.residual_j, 0.0, 2e-5);
    assert_eq!(report.energy.substeps, 1);
    assert_eq!(report.rejected_intervals, 0);
    assert_eq!(next.mass_kg(), thermal.mass_kg());
    assert_eq!(thermal.time_s(), 0.0);
    assert!(report.entropy.advection.is_none());
    assert!(report.entropy.has_unclassified_sources);

    let audited = system
        .with_implicit_advection_entropy(&thermal, || Ok(()))
        .unwrap();
    let (classified, classified_report) = audited
        .advance_to(
            &thermal,
            dt,
            implicit_controls(dt),
            &mut ImplicitConductionWorkspace::new(2).unwrap(),
            || Ok(()),
        )
        .unwrap();
    assert_eq!(classified.cells(), next.cells());
    assert_eq!(classified_report.energy, report.energy);
    let temperatures = [0, 1].map(|i| expected_t(fractions[i], excess_cp, expected[i]));
    let s = [0, 1].map(|i| entropy::expected_s(fractions[i], excess_cp, temperatures[i]));
    let inlet_s = [
        entropy::expected_s([1.0, 0.0], excess_cp, 460.0),
        entropy::expected_s([0.0, 1.0], excess_cp, 330.0),
    ];
    let expected_boundary = dt * (inlet_s[0] + inlet_s[1] - 2.0 * s[1]);
    let expected_mixing = dt
        * ((inlet[0] - expected[0]) / temperatures[0] - (inlet_s[0] - s[0])
            + (inlet[1] - expected[1]) / temperatures[1]
            - (inlet_s[1] - s[1])
            + 2.0 * ((expected[0] - expected[1]) / temperatures[1] - (s[0] - s[1]))
            + (expected[1] - expected[0]) / temperatures[0]
            - (s[1] - s[0]));
    let expected_storage: f64 = (0..2)
        .map(|i| masses[i] * (s[i] - entropy::expected_s(fractions[i], excess_cp, initial_t[i])))
        .sum();
    let generation = dt
        * 1000.0
        * (temperatures[1] - temperatures[0])
        * (1.0 / temperatures[0] - 1.0 / temperatures[1]);
    let account = classified_report.entropy;
    let transported = account.advection.unwrap();
    near(
        transported.boundary_entropy_input_j_k,
        expected_boundary,
        2e-7,
    );
    near(transported.upwind_mixing_j_k, expected_mixing, 2e-7);
    near(transported.balance_residual_j_k, 0.0, 2e-10);
    near(account.thermal_storage_change_j_k, expected_storage, 2e-7);
    near(account.conduction_generation_j_k, generation, 2e-8);
    near(
        account.temporal_defect_j_k,
        expected_storage - expected_boundary - expected_mixing - generation,
        3e-7,
    );
    assert!(!account.has_unclassified_sources);
    assert_eq!(account.unclassified_source_at_cell_j_k, 0.0);
}

#[test]
fn common_stream_reference_retains_sources_and_exact_fractions_without_changing_closed_storage() {
    let card = card(CompositionBasis::MoleFraction, 1.0);
    let liquid = liquid_with_excess_cp(&card, 100.0, -40.0);
    let caloric =
        BinaryLiquidCaloricLaw::try_new(&liquid, variable_sources(300.0, 0), "pressure").unwrap();
    let y = [0.1, 0.9];
    let storage = caloric
        .stream_storage_from_mass_fractions(y, controls(), || false)
        .unwrap();
    let descriptor = storage.advected_composition().unwrap();
    let entropy_anchor = storage.advected_entropy_anchor().unwrap();
    near(
        entropy_anchor.reduced_specific_entropy_j_kg_k,
        entropy::expected_s(y, -40.0, 290.0),
        2e-10,
    );
    assert_eq!(
        entropy_anchor.specific_h_j_kg,
        storage.specific_enthalpy_domain_j_kg()[0]
    );
    assert_eq!(descriptor.species, ["A", "B"]);
    assert_eq!(descriptor.mass_fractions, y);
    assert_eq!(storage.sources(), caloric.sources());
    let pure = caloric
        .stream_storage_from_mass_fractions([1.0, 0.0], controls(), || false)
        .unwrap();
    assert_eq!(
        descriptor.enthalpy_reference,
        pure.advected_composition().unwrap().enthalpy_reference
    );
    assert_ne!(storage.caloric_identity(), pure.caloric_identity());
    assert_eq!(
        entropy_anchor.entropy_reference,
        pure.advected_entropy_anchor().unwrap().entropy_reference
    );
    let mut varied = controls();
    varied.maximum_iterations += 1;
    varied.energy_tolerance_j *= 2.0;
    let other = caloric
        .stream_storage_from_mass_fractions(y, varied, || false)
        .unwrap();
    assert_eq!(storage.caloric_identity(), other.caloric_identity());
    varied.temperature_bracket = [Temperature::new(300.0), Temperature::new(550.0)];
    let narrow = caloric
        .stream_storage_from_mass_fractions(y, varied, || false)
        .unwrap();
    assert_eq!(
        descriptor.enthalpy_reference,
        narrow.advected_composition().unwrap().enthalpy_reference
    );
    assert_ne!(storage.caloric_identity(), narrow.caloric_identity());
    let narrow_anchor = narrow.advected_entropy_anchor().unwrap();
    assert_eq!(
        entropy_anchor.entropy_reference,
        narrow_anchor.entropy_reference
    );
    near(
        narrow_anchor.reduced_specific_entropy_j_kg_k,
        entropy::expected_s(y, -40.0, 300.0),
        2e-10,
    );
    for t in [290.0, 300.0, 350.0, 400.0, 500.0, 600.0] {
        let enthalpy = h(&storage, t);
        near(enthalpy, expected_h(y, -40.0, t), 2e-8);
        near(storage.temperature_k(enthalpy).unwrap(), t, 2e-8);
        let [_, p, q] = coefficients(y, -40.0);
        near(
            storage
                .specific_heat_cp_j_kg_k(Temperature::new(t), || false)
                .unwrap(),
            p + 2.0 * q * (t - 300.0),
            2e-10,
        );
    }
    let closed = caloric
        .thermal_storage([y[0] / 0.018, y[1] / 0.10], controls(), || false)
        .unwrap();
    assert_eq!(closed.specific_enthalpy_domain_j_kg()[0], 0.0);
    assert!(closed.advected_composition().is_none());
    assert!(closed.advected_entropy_anchor().is_none());
    assert_eq!(closed.caloric_identity(), None);
    near(
        h(&storage, 500.0) - h(&storage, 350.0),
        h(&closed, 500.0) - h(&closed, 350.0),
        2e-8,
    );
    near(
        storage
            .specific_entropy_change_j_kg_k(h(&storage, 350.0), h(&storage, 500.0))
            .unwrap()
            .unwrap(),
        closed
            .specific_entropy_change_j_kg_k(h(&closed, 350.0), h(&closed, 500.0))
            .unwrap()
            .unwrap(),
        2e-9,
    );
    let [lower, upper] = storage.specific_enthalpy_domain_j_kg();
    for invalid in [lower.next_down(), upper.next_up(), f64::NAN] {
        assert_eq!(
            storage.temperature_k(invalid),
            Err(DepositionError::MaterialDomain)
        );
    }
    for invalid in [[0.2, 0.7], [-0.1, 1.1], [f64::NAN, 1.0], [0.0, 0.0]] {
        assert!(
            caloric
                .stream_storage_from_mass_fractions(invalid, controls(), || false)
                .is_err()
        );
    }
    assert!(matches!(
        caloric.stream_storage_from_mass_fractions(y, controls(), || true),
        Err(BinaryLiquidCaloricError::Cancelled)
    ));
}

#[test]
fn actual_flow_binding_refuses_changed_source_datum_or_composition_atomically() {
    let card = card(CompositionBasis::MoleFraction, 1.0);
    let liquid = liquid_with_excess_cp(&card, 100.0, 40.0);
    let original =
        BinaryLiquidCaloricLaw::try_new(&liquid, variable_sources(300.0, 0), "pressure").unwrap();
    for (datum, revision, wrong_composition) in
        [(300.0, 1, false), (350.0, 0, false), (300.0, 0, true)]
    {
        let other =
            BinaryLiquidCaloricLaw::try_new(&liquid, variable_sources(datum, revision), "pressure")
                .unwrap();
        let y = if wrong_composition {
            [0.5, 0.5]
        } else {
            [0.25, 0.75]
        };
        let storage = vec![
            original
                .stream_storage_from_mass_fractions(y, controls(), || false)
                .unwrap(),
            other
                .stream_storage_from_mass_fractions([1.0, 0.0], controls(), || false)
                .unwrap(),
            original
                .stream_storage_from_mass_fractions([0.0, 1.0], controls(), || false)
                .unwrap(),
        ];
        let flow = mixer_flow([h(&storage[1], 400.0), h(&storage[2], 400.0)]);
        if !wrong_composition {
            assert_ne!(
                storage[0]
                    .advected_entropy_anchor()
                    .unwrap()
                    .entropy_reference,
                storage[1]
                    .advected_entropy_anchor()
                    .unwrap()
                    .entropy_reference
            );
        }
        let laws = conduction_laws(&storage);
        let (thermal, conduction) = prepare(&laws, &[400.0], &[2.0], false);
        let before = thermal.cells().clone();
        let system =
            ImplicitConductionChain::new(&thermal, &conduction, &bindings(&laws), || Ok(()))
                .unwrap();
        assert!(
            system
                .with_implicit_advection(&thermal, &flow, || Ok(()))
                .is_err()
        );
        assert_eq!(thermal.cells(), &before);
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(flow.audit().substeps, 0);
    }
    let changed_pair = liquid_with_excess_cp(&card, 100.0, -40.0);
    let changed =
        BinaryLiquidCaloricLaw::try_new(&changed_pair, variable_sources(300.0, 0), "pressure")
            .unwrap();
    let a = original
        .stream_storage_from_mass_fractions([0.25, 0.75], controls(), || false)
        .unwrap();
    let b = changed
        .stream_storage_from_mass_fractions([0.25, 0.75], controls(), || false)
        .unwrap();
    assert_ne!(
        a.advected_composition().unwrap().enthalpy_reference,
        b.advected_composition().unwrap().enthalpy_reference
    );
}
