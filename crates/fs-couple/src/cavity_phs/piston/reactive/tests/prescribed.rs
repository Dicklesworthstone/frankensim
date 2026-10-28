//! Driven open chemistry against independent inventory, caloric and rate equations.
//! Manufactured isomers use the original certified stoichiometry and source laws.
use super::*;
use crate::prescribed_gas::reactive::{self, DrivenGasReaction, ReactivePrescribedGasSolve};
use crate::prescribed_gas::{self, PrescribedGasDrive, PrescribedGasSolve};
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::GasNozzle;
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::network::snapshot::SnapshotNozzleLink;
use fs_material::thermochemical_gas::energy::parcel::transfer::snapshot::GasSnapshotSolve;
use std::collections::BTreeMap;

fn controls() -> ReactivePrescribedGasSolve {
    ReactivePrescribedGasSolve {
        network: PrescribedGasSolve {
            exchange: GasSnapshotSolve {
                maximum_withdrawn_fraction: 0.1,
                mass_tolerance_kg: 1e-14,
                energy_tolerance_j: 1e-8,
            },
            total_energy_tolerance_j: 1e-7,
        },
        species_tolerance_mol: 1e-12,
    }
}
fn prescribed(volume: f64, heat: f64) -> PrescribedGasDrive {
    PrescribedGasDrive {
        final_volume_m3: volume,
        heat_j: heat,
        exterior_pressure_pa: 80000.0,
    }
}
fn at(law: &GasReaction, fractions: [f64; 3], pressure: f64) -> ClosedGasParcel {
    ClosedGasParcel::new(
        law.catalog()
            .with_composition(
                Composition::new(CompositionBasis::MoleFraction, fractions.to_vec()).unwrap(),
            )
            .unwrap(),
        QueryPoint::new()
            .with_quantity(
                "temperature",
                quantity(QuantityKind::AbsoluteTemperature),
                600.0,
            )
            .unwrap()
            .with_quantity("pressure", quantity(QuantityKind::Pressure), pressure)
            .unwrap()
            .with("specimen-lot", 42.0)
            .unwrap(),
        0.001,
    )
    .unwrap()
}
fn endpoints(gases: &BTreeMap<u64, ClosedGasParcel>) -> BTreeMap<u64, GasEnergySolve> {
    gases
        .keys()
        .map(|&id| (id, chemical_controls().gas))
        .collect()
}
fn binding(law: &GasReaction, chamber: u64) -> DrivenGasReaction<'_> {
    DrivenGasReaction {
        chamber,
        reaction: law,
        controls: chemical_controls(),
    }
}

#[test]
fn formation_energy_drives_expansion_without_a_second_heat_supply() {
    let law = reaction(4000.0, 4.0, 1000.0);
    let initial = gas(&law, [0.75, 0.0, 0.25]);
    let original = initial.clone();
    let gases = BTreeMap::from([(1, initial.clone())]);
    let duration = 5e-5;
    let final_volume = 0.00108;
    let supplied_heat = 3.0;
    let result = reactive::step_candidate(
        &gases,
        &[],
        duration,
        &endpoints(&gases),
        &BTreeMap::from([(1, prescribed(final_volume, supplied_heat))]),
        &[binding(&law, 1)],
        controls(),
        || false,
    )
    .unwrap();
    let n0 = law.inventory(&initial).unwrap();
    let extent = duration * initial.volume_m3() * 1000.0 * 0.75;
    let amounts = [n0[0] - extent, extent, n0[2]];
    let capacity = R * (2.5 * amounts[0] + 3.0 * amounts[1] + 1.5 * amounts[2]);
    let u0 = R * 600.0 * (2.5 * n0[0] + 1.5 * n0[2]);
    let reacted_temperature = (u0 + 4000.0 * extent) / capacity;
    let midpoint_volume = 0.5 * (initial.volume_m3() + final_volume);
    let dv = final_volume - initial.volume_m3();
    let nr = amounts.iter().sum::<f64>() * R;
    let midpoint_temperature = (capacity * reacted_temperature + 0.5 * supplied_heat)
        / (capacity + 0.5 * nr * dv / midpoint_volume);
    let expected_temperature = 2.0 * midpoint_temperature - reacted_temperature;
    let expected_work = nr * midpoint_temperature * dv / midpoint_volume;
    near(result.reactions()[&1].extent_mol, extent, 1e-17);
    near(
        result.parcels()[&1].state().state().temperature,
        expected_temperature,
        2e-7,
    );
    near(result.step().gross_gas_work_j, expected_work, 1e-8);
    near(
        result.step().net_drive_work_j,
        expected_work - 80000.0 * dv,
        1e-8,
    );
    near(result.step().supplied_heat_j, supplied_heat, 0.0);
    near(
        result.step().gas_energy_change_j + expected_work,
        supplied_heat,
        1e-8,
    );
    assert!(expected_work > 0.0);
    assert!(reacted_temperature > 600.0);
    assert!(result.maximum_species_residual_mol() <= controls().species_tolerance_mol);
    assert_eq!(gases[&1], original);
    assert_eq!(
        result.parcels()[&1].point().axes().get("specimen-lot"),
        Some(&42.0)
    );
}

fn temperature_rate(law: GasReaction, scale: f64) -> GasReaction {
    use crate::cavity_phs::piston::reactive::material_rate::{
        ACTIVITY_RATE_SCALE_DIMS, ACTIVITY_RATE_SCALE_PROPERTY, MaterialRateScale,
    };
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(ACTIVITY_RATE_SCALE_PROPERTY, ACTIVITY_RATE_SCALE_DIMS),
            value: PropertyValue::Curve {
                abscissa: "temperature".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(200.0, scale / 3.0), (1800.0, 3.0 * scale)],
                dims: ACTIVITY_RATE_SCALE_DIMS,
            },
            validity: ValidityDomain::unconstrained().with_quantity(
                "temperature",
                quantity(QuantityKind::AbsoluteTemperature),
                200.0,
                1800.0,
            ),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: provenance(),
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured-driven-isomer-rate".into(),
            phase: "gas".into(),
            process: "linear-temperature-scale".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    law.with_material_rate(
        MaterialRateScale::new(
            card,
            MaterialPropertySelection::SingleClaimOnly,
            Pressure::new(100000.0),
        )
        .unwrap(),
    )
    .unwrap()
}

#[test]
fn simultaneous_inlet_outlet_reacts_actual_incoming_species_and_closes_total_enthalpy() {
    let law = temperature_rate(reaction(4000.0, 3.5, 1.0), 1000.0);
    let gases = BTreeMap::from([
        (1, at(&law, [1.0, 0.0, 0.0], 180000.0)),
        (2, at(&law, [0.0, 0.0, 1.0], 100000.0)),
        (3, at(&law, [0.0, 0.0, 1.0], 50000.0)),
    ]);
    let original = gases.clone();
    let links = [
        SnapshotNozzleLink {
            id: 10,
            parcels: [1, 2],
            nozzle: GasNozzle {
                area_m2: 2e-7,
                discharge_coefficient: 0.8,
            },
            nozzle_solve: chemical_controls().gas,
        },
        SnapshotNozzleLink {
            id: 20,
            parcels: [2, 3],
            nozzle: GasNozzle {
                area_m2: 1e-7,
                discharge_coefficient: 0.8,
            },
            nozzle_solve: chemical_controls().gas,
        },
    ];
    let drives = BTreeMap::from([
        (1, prescribed(0.001, 0.0)),
        (2, prescribed(0.000999, 0.0)),
        (3, prescribed(0.001, 0.0)),
    ]);
    let duration = 5e-5;
    let run = |links: &[SnapshotNozzleLink]| {
        reactive::step_candidate(
            &gases,
            links,
            duration,
            &endpoints(&gases),
            &drives,
            &[binding(&law, 2)],
            controls(),
            || false,
        )
        .unwrap()
    };
    let result = run(&links);
    let flow = &result.step().flow;
    let inlet = flow
        .links
        .iter()
        .find(|r| r.id == 10)
        .unwrap()
        .transferred_mass_kg;
    let outlet = flow
        .links
        .iter()
        .find(|r| r.id == 20)
        .unwrap()
        .transferred_mass_kg;
    let a = inlet / 0.028;
    let inert = gases[&2].mass_kg() / 0.004 - outlet / 0.004;
    let incoming_enthalpy = a * R * 3.5 * 600.0;
    let outgoing_enthalpy = outlet / 0.004 * R * 2.5 * 600.0;
    let internal =
        gases[&2].mass_kg() / 0.004 * R * 1.5 * 600.0 + incoming_enthalpy - outgoing_enthalpy;
    let mixed_temperature = internal / (R * (2.5 * a + 1.5 * inert));
    let mixed_pressure_a = a * R * mixed_temperature / 0.001;
    let expected_extent =
        duration * 0.001 * 1000.0 * (mixed_temperature / 600.0) * (mixed_pressure_a / 100000.0);
    let report = &result.reactions()[&2];
    near(report.species[0].initial_mol, a, 2e-15);
    near(report.species[2].initial_mol, inert, 2e-15);
    near(report.extent_mol, expected_extent, 2e-16);
    assert!(report.extent_mol > 0.0);
    let rate_point = report.rate_material.as_ref().unwrap().query();
    near(rate_point.axes()["temperature"], mixed_temperature, 1e-7);
    near(
        result.step().parcels[&2].net_enthalpy_into_gas_j,
        incoming_enthalpy - outgoing_enthalpy,
        1e-8,
    );
    assert!(result.step().total_energy_residual_j.abs() < 1e-7);
    assert!(result.maximum_species_residual_mol() < 1e-12);
    assert_eq!(result, run(&[links[1], links[0]]));
    assert_eq!(gases, original);
}

#[test]
fn existing_implicit_policy_survives_driven_compression() {
    use crate::cavity_phs::piston::reactive::implicit::ImplicitReactionSolve;
    let law = reaction(0.0, 3.5, 1000.0)
        .with_implicit(ImplicitReactionSolve {
            extent_bracket_mol: [-0.01, 0.01],
            extent_tolerance_mol: 1e-13,
            maximum_iterations: 128,
        })
        .unwrap();
    let gases = BTreeMap::from([(1, gas(&law, [0.75, 0.0, 0.25]))]);
    let duration = 5e-5;
    let result = reactive::step_candidate(
        &gases,
        &[],
        duration,
        &endpoints(&gases),
        &BTreeMap::from([(1, prescribed(0.0009, 0.0))]),
        &[binding(&law, 1)],
        controls(),
        || false,
    )
    .unwrap();
    let n0 = law.inventory(&gases[&1]).unwrap();
    let z = duration * 1000.0 * R * 600.0 / 100000.0;
    let extent = z * (n0[0] - n0[1]) / (1.0 + 2.0 * z);
    near(result.reactions()[&1].extent_mol, extent, 2e-13);
    assert!(result.reactions()[&1].implicit.is_some());
    assert!(result.parcels()[&1].state().state().temperature > 600.0);
    assert!(result.step().gross_gas_work_j < 0.0);
    assert!(result.step().total_energy_residual_j.abs() < 1e-7);
}

#[test]
fn empty_chemistry_preserves_original_driven_solution_exactly() {
    let law = reaction(4000.0, 4.0, 1000.0);
    let gases = BTreeMap::from([(1, gas(&law, [0.75, 0.0, 0.25]))]);
    let drives = BTreeMap::from([(1, prescribed(0.00095, 5.0))]);
    let base = prescribed_gas::step_candidate(
        &gases,
        &[],
        5e-5,
        &endpoints(&gases),
        &drives,
        controls().network,
        || false,
    )
    .unwrap();
    let result = reactive::step_candidate(
        &gases,
        &[],
        5e-5,
        &endpoints(&gases),
        &drives,
        &[],
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(result.parcels(), base.parcels());
    assert_eq!(result.step(), base.step());
    assert!(result.reactions().is_empty());
}

#[test]
fn full_step_refinement_converges_to_independent_adiabatic_reaction_solution() {
    let law = reaction(0.0, 3.5, 1000.0);
    let initial = gas(&law, [0.75, 0.0, 0.25]);
    let n0 = law.inventory(&initial).unwrap();
    let total = n0.iter().sum::<f64>();
    let alpha = total / (2.5 * n0[0] + 1.5 * n0[2]);
    let ratio = 0.7_f64;
    let duration = 0.002;
    let integral_temperature =
        600.0 * duration * (ratio.powf(1.0 - alpha) - 1.0) / ((ratio - 1.0) * (1.0 - alpha));
    let expected_b = 0.5
        * (n0[0] + n0[1]
            - (n0[0] - n0[1]) * (-2.0 * 1000.0 * R * integral_temperature / 100000.0).exp());
    let expected_t = 600.0 * ratio.powf(-alpha);
    let mut errors = Vec::new();
    for count in [8, 16, 32] {
        let mut gases = BTreeMap::from([(1, initial.clone())]);
        let mut work = 0.0;
        for i in 1..=count {
            let volume = initial.volume_m3() * (1.0 + (ratio - 1.0) * i as f64 / count as f64);
            let result = reactive::step_candidate(
                &gases,
                &[],
                duration / count as f64,
                &endpoints(&gases),
                &BTreeMap::from([(1, prescribed(volume, 0.0))]),
                &[binding(&law, 1)],
                controls(),
                || false,
            )
            .unwrap();
            work += result.step().gross_gas_work_j;
            gases = result.into_parts().0;
        }
        let final_gas = &gases[&1];
        let amounts = law.inventory(final_gas).unwrap();
        errors.push((amounts[1] - expected_b).abs());
        near(
            initial.mass_kg()
                * (final_gas.state().specific_internal_energy_j_kg()
                    - initial.state().specific_internal_energy_j_kg())
                + work,
            0.0,
            1e-7,
        );
        near(final_gas.state().state().temperature, expected_t, 0.06);
    }
    assert!(
        errors[0] / errors[1] > 1.7 && errors[1] / errors[2] > 1.7,
        "{errors:?}"
    );
    assert!(
        errors[0] / errors[1] < 2.3 && errors[1] / errors[2] < 2.3,
        "{errors:?}"
    );
    assert!(errors[2] < 2e-6, "{errors:?}");
}

#[test]
fn late_source_failure_duplicate_bindings_and_cancellation_are_atomic() {
    let law = reaction(4000.0, 4.0, 1000.0);
    let gases = BTreeMap::from([
        (1, gas(&law, [0.75, 0.0, 0.25])),
        (2, gas(&law, [0.75, 0.0, 0.25])),
    ]);
    let original = gases.clone();
    let drives = BTreeMap::from([(1, prescribed(0.00099, 0.0)), (2, prescribed(0.00099, 0.0))]);
    let bindings = [binding(&law, 1), binding(&law, 2)];
    let run = |drives: &BTreeMap<u64, PrescribedGasDrive>, bindings: &[DrivenGasReaction<'_>]| {
        reactive::step_candidate(
            &gases,
            &[],
            5e-5,
            &endpoints(&gases),
            drives,
            bindings,
            controls(),
            || false,
        )
    };
    let good = run(&drives, &bindings).unwrap();
    let mut bad = drives.clone();
    bad.get_mut(&2).unwrap().heat_j = 1e12;
    assert!(run(&bad, &bindings).is_err());
    assert!(run(&drives, &[binding(&law, 1), binding(&law, 1)]).is_err());
    assert!(run(&drives, &[binding(&law, 3)]).is_err());
    let mut polls = 0;
    reactive::step_candidate(
        &gases,
        &[],
        5e-5,
        &endpoints(&gases),
        &drives,
        &bindings,
        controls(),
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    let mut current = 0;
    assert_eq!(
        reactive::step_candidate(
            &gases,
            &[],
            5e-5,
            &endpoints(&gases),
            &drives,
            &bindings,
            controls(),
            || {
                current += 1;
                current >= polls - 2
            }
        )
        .unwrap_err(),
        prescribed_gas::PrescribedGasError::Cancelled
    );
    assert_eq!(gases, original);
    assert_eq!(run(&drives, &bindings).unwrap(), good);
}
