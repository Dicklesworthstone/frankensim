//! G1/G3/G4: source-driven flow, mixed chemistry and one spatial material wall.
use super::*;
use fs_couple::solid_enthalpy::gas::network::{exchange::*, *};
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::{
    GasNozzle, network::snapshot::SnapshotNozzleLink,
};
use fs_material::thermochemical_gas::energy::parcel::transfer::snapshot::GasSnapshotSolve;
use fs_qty::QuantitySpec;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_thermochem::{Composition, CompositionBasis};

const R: f64 = fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K;
const MOLAR_MASS: f64 = 0.028;
const RS: f64 = R / MOLAR_MASS;
const CV: f64 = 2.5 * RS;

fn controls() -> OpenSolidGasSolve {
    OpenSolidGasSolve {
        flow: GasSnapshotSolve {
            maximum_withdrawn_fraction: 0.2,
            mass_tolerance_kg: 1e-12,
            energy_tolerance_j: 1e-7,
        },
        thermal: SolidGasNetworkSolve {
            maximum_gases: 2,
            coupled: solve(),
        },
        total_energy_tolerance_j: 5e-7,
    }
}
fn link(area: f64) -> SnapshotNozzleLink {
    SnapshotNozzleLink {
        id: 77,
        parcels: [20, 10],
        nozzle: GasNozzle {
            area_m2: area,
            discharge_coefficient: 0.8,
        },
        nozzle_solve: solve().gas,
    }
}
fn gas_at(
    law: &GasReaction,
    b_fraction: f64,
    temperature: f64,
    pressure: f64,
) -> IsochoricGasParcel {
    let model = law
        .catalog()
        .with_composition(
            Composition::new(
                CompositionBasis::MoleFraction,
                vec![1.0 - b_fraction, b_fraction],
            )
            .unwrap(),
        )
        .unwrap();
    let quantity = |kind| QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
    let point = fs_matdb::QueryPoint::new()
        .with_quantity(
            "temperature",
            quantity(QuantityKind::AbsoluteTemperature),
            temperature,
        )
        .unwrap()
        .with_quantity("pressure", quantity(QuantityKind::Pressure), pressure)
        .unwrap()
        .with("catalyst-lot", 7.0)
        .unwrap();
    IsochoricGasParcel::new(model, point, VOLUME).unwrap()
}
fn b_fraction(gas: &IsochoricGasParcel) -> f64 {
    // These fixtures retain both positive species in the original A,B order.
    gas.model().thermodynamics().mole_composition().fractions()[1]
}

// Exact choked ideal-gas mass flux for constant gamma=1.4, in SI units.
fn mass_flow(p: f64, t: f64, area: f64) -> f64 {
    0.8 * area * p * (1.4 / (RS * t)).sqrt() * (2.0_f64 / 2.4).powf(2.4 / 0.8)
}

// Independent four-store BE elimination after the conservative transfer.
// Temperatures: two wall cells, donor, receiver. Cp/Cv and film k are source
// fixture constants; capacities are recomputed from the CHANGED gas masses.
fn thermal_reference(t: [f64; 4], mass: [f64; 2], b: [f64; 2], wall_k: f64) -> [f64; 4] {
    let capacity = mass.map(|m| m * CV);
    let gas_k = b.map(|x| 0.03 * (1.0 + 3.0 * x));
    let g = [
        0.1 / (0.001 / gas_k[0] + 0.02 / wall_k),
        0.2 / (0.001 / gas_k[1] + 0.02 / wall_k),
    ];
    let a = [
        g[0] * capacity[0] / (capacity[0] + g[0]),
        g[1] * capacity[1] / (capacity[1] + g[1]),
    ];
    let d = [10.5 + a[0], 10.5 + a[1]];
    let rhs = [10.0 * t[0] + a[0] * t[2], 10.0 * t[1] + a[1] * t[3]];
    let det = d[0] * d[1] - 0.25;
    let wall = [
        (rhs[0] * d[1] + 0.5 * rhs[1]) / det,
        (rhs[1] * d[0] + 0.5 * rhs[0]) / det,
    ];
    [
        wall[0],
        wall[1],
        (capacity[0] * t[2] + g[0] * wall[0]) / (capacity[0] + g[0]),
        (capacity[1] * t[3] + g[1] * wall[1]) / (capacity[1] + g[1]),
    ]
}

#[test]
fn material_wall_changes_continued_nozzle_flow_species_heat_and_native_deformation() {
    with_cx(|cx, _| {
        let law = reaction(1.0, 0.0, 4.0);
        let mut observed = Vec::new();
        for wall_k in [0.1, 10.0] {
            let (card, caloric) = wall(wall_k, false);
            let mut thermal = initial(caloric.curve());
            let context = point(300.0);
            let model = solid(&thermal, &card, caloric.curve(), &context, cx);
            let mut donor = gas_at(&law, 0.8, 400.0, 300000.0);
            let mut receiver = gas_at(&law, 0.1, 300.0, 100000.0);
            let initial_energy = thermal.energy_j()
                + donor.mass_kg() * donor.state().specific_internal_energy_j_kg()
                + receiver.mass_kg() * receiver.state().specific_internal_energy_j_kg();
            let mut expected = [300.0, 300.0, 400.0, 300.0];
            let mut masses = [
                300000.0 * VOLUME / (RS * 400.0),
                100000.0 * VOLUME / (RS * 300.0),
            ];
            let mut composition = [0.8, 0.1];
            let patches = contacts();
            let mut last_flow = 0.0;
            for interval in 1..=3 {
                let inputs = [
                    SolidGasCompartment {
                        id: 20,
                        gas: &donor,
                        external_heat_j: 0.0,
                        contacts: &patches[..1],
                        reaction: None,
                    },
                    SolidGasCompartment {
                        id: 10,
                        gas: &receiver,
                        external_heat_j: 0.0,
                        contacts: &patches[1..],
                        reaction: None,
                    },
                ];
                let candidate = model
                    .step_with_open_gas_network(
                        &thermal,
                        &inputs,
                        &[link(1e-5)],
                        interval as f64,
                        None,
                        controls(),
                        &[0.0; 4],
                        cx,
                    )
                    .unwrap();
                let dm = mass_flow(masses[0] * RS * expected[2] / VOLUME, expected[2], 1e-5);
                let transferred_mass = [masses[0] - dm, masses[1] + dm];
                let after_flow = [
                    expected[0],
                    expected[1],
                    expected[2] * (masses[0] - 1.4 * dm) / transferred_mass[0],
                    (masses[1] * expected[3] + 1.4 * dm * expected[2]) / transferred_mass[1],
                ];
                composition[1] =
                    (masses[1] * composition[1] + dm * composition[0]) / transferred_mass[1];
                masses = transferred_mass;
                expected = thermal_reference(after_flow, masses, composition, wall_k);
                let report = candidate.report();
                assert!(report.flow.links[0].flow.choked);
                near(report.flow.links[0].transferred_mass_kg, dm, 1e-11);
                last_flow = report.flow.links[0].flow.mass_flow_kg_s;
                for (index, id) in [20, 10].into_iter().enumerate() {
                    let gas = candidate.gas(id).unwrap();
                    near(gas.mass_kg(), masses[index], 1e-12);
                    near(gas.state().state().temperature, expected[index + 2], 3e-7);
                    near(b_fraction(gas), composition[index], 1e-11);
                    near(
                        gas.state().state().thermal_conductivity,
                        0.03 * (1.0 + 3.0 * composition[index]),
                        1e-12,
                    );
                }
                near(
                    candidate.thermal().temperature_k(CellId(30)).unwrap(),
                    expected[0],
                    3e-7,
                );
                near(
                    candidate.thermal().temperature_k(CellId(10)).unwrap(),
                    expected[1],
                    3e-7,
                );
                assert!(report.total_mass_residual_kg.abs() <= controls().flow.mass_tolerance_kg);
                assert!(
                    report.total_energy_residual_j.abs() <= controls().total_energy_tolerance_j
                );
                let loads = model
                    .geometry()
                    .assemble(&[], VolumeLoad3::None, &[], cx)
                    .unwrap();
                let free = candidate
                    .state()
                    .solve(&loads, &FREE, config(), cx)
                    .unwrap();
                near(
                    free.solution.displacement.displacement_m()[3][2],
                    1e-5 * (expected[0] - 300.0),
                    1e-10,
                );
                if interval == 1 {
                    let reversed = [inputs[1], inputs[0]];
                    let renamed = SnapshotNozzleLink {
                        id: 2,
                        parcels: [10, 20],
                        ..link(1e-5)
                    };
                    let replay = model
                        .step_with_open_gas_network(
                            &thermal,
                            &reversed,
                            &[renamed],
                            1.0,
                            None,
                            controls(),
                            &[0.0; 4],
                            cx,
                        )
                        .unwrap();
                    assert_eq!(candidate.gas(20), replay.gas(20));
                    assert_eq!(candidate.gas(10), replay.gas(10));
                    assert_eq!(candidate.thermal().cells(), replay.thermal().cells());
                }
                let (accepted, mut gases, _, _) = candidate.into_parts();
                donor = gases.pop().unwrap().gas;
                receiver = gases.pop().unwrap().gas;
                thermal = accepted;
                near(
                    thermal.energy_j()
                        + donor.mass_kg() * donor.state().specific_internal_energy_j_kg()
                        + receiver.mass_kg() * receiver.state().specific_internal_energy_j_kg()
                        - initial_energy,
                    0.0,
                    2e-6,
                );
            }
            observed.push((last_flow, expected[0]));
        }
        assert!((observed[0].0 - observed[1].0).abs() > 1e-7);
        assert!(observed[1].1 - observed[0].1 > 10.0);
    });
}

#[test]
fn mixed_inlet_inventory_feeds_native_reaction_and_complete_heat_network_once() {
    with_cx(|cx, _| {
        let law = super::coupled::implicit(reaction(1.0, 20000.0, 4.0));
        let donor = gas_at(&law, 0.8, 400.0, 300000.0);
        let receiver = gas_at(&law, 0.1, 300.0, 100000.0);
        let (card, caloric) = wall(1.0, false);
        let thermal = initial(caloric.curve());
        let context = point(300.0);
        let model = solid(&thermal, &card, caloric.curve(), &context, cx);
        let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
        let cells = [30, 10].map(|id| MaterialCell {
            cell: CellId(id),
            material_id: MaterialId(0),
            curve: caloric.curve(),
            point: ThermalMaterialPoint {
                material: &card,
                context: &context,
                temperature_axis: "T",
                liquid_fraction_axis: None,
            },
            conductivity: Some((&key, MaterialPropertySelection::SingleClaimOnly)),
        });
        let links = [ConductionLink {
            link: 7,
            distance_a_m: 0.3,
            distance_b_m: 0.7,
        }];
        let bulk = MaterialNetwork {
            cells: &cells,
            resistors: &[],
            links: &links,
        };
        let electrical =
            fs_couple::circuit_enthalpy::contacts::material::ElectricalContactNetwork {
                bulk: &bulk,
                contacts: &[],
                electrical_bases: &[],
            };
        let network = fs_couple::circuit_enthalpy::contacts::heat::MaterialHeatNetwork {
            electrical: &electrical,
            thermal_contacts: &[],
            radiation: None,
        };
        let prepared = network.prepare(cx, &thermal).unwrap();
        let patches = contacts();
        let inputs = [
            SolidGasCompartment {
                id: 20,
                gas: &donor,
                external_heat_j: 0.0,
                contacts: &patches[..1],
                reaction: None,
            },
            SolidGasCompartment {
                id: 10,
                gas: &receiver,
                external_heat_j: 0.0,
                contacts: &patches[1..],
                reaction: Some((&law, chemical())),
            },
        ];
        let candidate = model
            .step_with_open_gas_heat_network(
                &thermal,
                &inputs,
                &[link(1e-5)],
                1.0,
                &prepared,
                controls(),
                &[0.0; 4],
                cx,
            )
            .unwrap();
        let dm = mass_flow(300000.0, 400.0, 1e-5);
        let mass = [donor.mass_kg() - dm, receiver.mass_kg() + dm];
        let initial_b = (receiver.mass_kg() * 0.1 + dm * 0.8) / MOLAR_MASS;
        let amount = mass[1] / MOLAR_MASS;
        let final_t = candidate.gas(10).unwrap().state().state().temperature;
        // Eliminate the native BE species equation at the ACTUAL final T and
        // changed mole inventory. This formula does not query production kinetics.
        let q = R * final_t * final_t / (100000.0 * 300.0);
        let reverse = (-20000.0 / (R * final_t)).exp();
        let xi = q * (amount - (1.0 + reverse) * initial_b) / (1.0 + q * (1.0 + reverse));
        let chemistry = candidate.report().thermal.gases[0]
            .endpoint
            .as_ref()
            .unwrap();
        near(chemistry.reaction.extent_mol, xi, 1e-9);
        near(chemistry.reaction.species[1].initial_mol, initial_b, 1e-10);
        near(
            chemistry.reaction.species[1].final_mol,
            initial_b + xi,
            1e-9,
        );
        assert!(candidate.report().thermal.gases[1].endpoint.is_none());
        let before = thermal.energy_j()
            + donor.mass_kg() * donor.state().specific_internal_energy_j_kg()
            + receiver.mass_kg() * receiver.state().specific_internal_energy_j_kg();
        let after = candidate.thermal().energy_j()
            + candidate
                .gases()
                .iter()
                .map(|g| g.gas.mass_kg() * g.gas.state().specific_internal_energy_j_kg())
                .sum::<f64>();
        near(after - before, 0.0, 1e-6);
        near(candidate.thermal().time_s(), 1.0, 0.0);
        assert!(candidate.report().thermal.transport.is_some());
    });
}

#[test]
fn adiabatic_reservoir_closed_nozzle_parity_and_late_refusal_preserve_originals() {
    with_cx(|cx, gate| {
        let law = reaction(1.0, 0.0, 4.0);
        let donor = gas_at(&law, 0.8, 400.0, 300000.0);
        let receiver = gas_at(&law, 0.1, 300.0, 100000.0);
        let saved = [donor.clone(), receiver.clone()];
        let (card, caloric) = wall(1.0, false);
        let thermal = initial(caloric.curve());
        let context = point(300.0);
        let model = solid(&thermal, &card, caloric.curve(), &context, cx);
        let patches = contacts();
        let inputs = [
            SolidGasCompartment {
                id: 20,
                gas: &donor,
                external_heat_j: 0.0,
                contacts: &[],
                reaction: None,
            },
            SolidGasCompartment {
                id: 10,
                gas: &receiver,
                external_heat_j: 0.0,
                contacts: &patches,
                reaction: None,
            },
        ];
        let closed = model
            .step_with_open_gas_network(
                &thermal,
                &inputs,
                &[link(0.0)],
                1.0,
                None,
                controls(),
                &[0.0; 4],
                cx,
            )
            .unwrap();
        let original = model
            .step_with_isochoric_gases(
                &thermal,
                &inputs,
                1.0,
                None,
                controls().thermal,
                &[0.0; 4],
                cx,
            )
            .unwrap();
        assert_eq!(closed.gas(20), Some(&donor));
        assert_eq!(closed.gas(10), original.gas(10));
        assert_eq!(closed.thermal().cells(), original.thermal().cells());
        let accepted = model
            .step_with_open_gas_network(
                &thermal,
                &inputs,
                &[link(1e-5)],
                1.0,
                None,
                controls(),
                &[0.0; 4],
                cx,
            )
            .unwrap();
        assert!(accepted.report().flow.links[0].transferred_mass_kg > 0.0);
        assert!(accepted.report().thermal.gases[1].contacts.is_empty());
        let mut limited = controls();
        limited.thermal.coupled.maximum_trials = 1;
        assert!(matches!(
            model.step_with_open_gas_network(
                &thermal,
                &inputs,
                &[link(1e-5)],
                1.0,
                None,
                limited,
                &[0.0; 4],
                cx
            ),
            Err(OpenSolidGasError::Thermal(SolidGasError::TrialBudget))
        ));
        let mut withdrawal = controls();
        withdrawal.flow.maximum_withdrawn_fraction = 1e-5;
        assert!(matches!(
            model.step_with_open_gas_network(
                &thermal,
                &inputs,
                &[link(1e-5)],
                1.0,
                None,
                withdrawal,
                &[0.0; 4],
                cx
            ),
            Err(OpenSolidGasError::Flow(_))
        ));
        let retry = model
            .step_with_open_gas_network(
                &thermal,
                &inputs,
                &[link(1e-5)],
                1.0,
                None,
                controls(),
                &[0.0; 4],
                cx,
            )
            .unwrap();
        assert_eq!(accepted.gas(20), retry.gas(20));
        assert_eq!(accepted.gas(10), retry.gas(10));
        gate.request();
        assert!(matches!(
            model.step_with_open_gas_network(
                &thermal,
                &inputs,
                &[link(1e-5)],
                1.0,
                None,
                controls(),
                &[0.0; 4],
                cx
            ),
            Err(OpenSolidGasError::Cancelled)
        ));
        assert_eq!(donor, saved[0]);
        assert_eq!(receiver, saved[1]);
        near(thermal.time_s(), 0.0, 0.0);
    });
}
