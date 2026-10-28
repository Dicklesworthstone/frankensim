//! Public gas owners against independent ideal-gas and cycle equations.
mod source;
use super::*;
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::GasNozzle;
use source::gas;

fn energy_solve() -> GasEnergySolve {
    GasEnergySolve {
        minimum_temperature_k: 200.0,
        maximum_temperature_k: 1800.0,
        energy_tolerance_j_kg: 1e-7,
        maximum_iterations: 128,
    }
}
fn solve() -> PrescribedGasSolve {
    PrescribedGasSolve {
        exchange: GasSnapshotSolve {
            maximum_withdrawn_fraction: 0.1,
            mass_tolerance_kg: 1e-14,
            energy_tolerance_j: 1e-8,
        },
        total_energy_tolerance_j: 1e-7,
    }
}
fn drive(volume: f64, heat: f64) -> PrescribedGasDrive {
    PrescribedGasDrive {
        final_volume_m3: volume,
        heat_j: heat,
        exterior_pressure_pa: 100000.0,
    }
}
fn near(a: f64, b: f64, absolute: f64) {
    assert!(
        (a - b).abs() <= absolute,
        "{a} != {b}, tolerance {absolute}"
    );
}
fn single(gas: &ClosedGasParcel, drive: PrescribedGasDrive) -> PrescribedGasCandidate {
    step_candidate(
        &BTreeMap::from([(1, gas.clone())]),
        &[],
        0.001,
        &BTreeMap::from([(1, energy_solve())]),
        &BTreeMap::from([(1, drive)]),
        solve(),
        || false,
    )
    .unwrap()
}

#[test]
fn prescribed_compression_and_heat_match_independent_constant_capacity_solution() {
    let initial = gas(0.75, 0.0, 400.0, 100000.0, 0.001);
    let r = fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / (0.75 * 0.028 + 0.25 * 0.004);
    let cv = 2.25 * r;
    for (volume, heat) in [(0.0006, 0.0), (0.0014, 20.0), (0.001, -10.0)] {
        let report = single(&initial, drive(volume, heat));
        let vm = 0.5 * (0.001 + volume);
        let ratio = (volume - 0.001) / vm;
        let tm = (cv * 400.0 + heat / (2.0 * initial.mass_kg())) / (cv + 0.5 * ratio * r);
        let t1 = 2.0 * tm - 400.0;
        let pressure = initial.mass_kg() * r * tm / vm;
        let work = pressure * (volume - 0.001);
        near(report.parcels()[&1].state().state().temperature, t1, 1e-8);
        near(
            report.step().parcels[&1].midpoint.state().pressure,
            pressure,
            1e-6,
        );
        near(report.step().gross_gas_work_j, work, 1e-8);
        near(initial.mass_kg() * cv * (t1 - 400.0) + work, heat, 1e-8);
        near(
            report.step().exterior_pressure_work_j,
            100000.0 * (volume - 0.001),
            1e-12,
        );
        near(
            report.step().net_drive_work_j,
            work - 100000.0 * (volume - 0.001),
            1e-8,
        );
        assert_eq!(report.parcels()[&1].mass_kg(), initial.mass_kg());
        assert_eq!(
            report.parcels()[&1].point().axes().get("specimen-lot"),
            Some(&42.0)
        );
        if volume < 0.001 {
            assert!(work < 0.0);
            assert!(t1 > 400.0);
        }
    }
}

fn march(initial: &ClosedGasParcel, final_volume: f64, steps: usize) -> (ClosedGasParcel, f64) {
    let mut gas = initial.clone();
    let mut work = 0.0;
    for index in 1..=steps {
        let volume = initial.volume_m3()
            + (final_volume - initial.volume_m3()) * (index as f64 / steps as f64);
        let candidate = single(&gas, drive(volume, 0.0));
        work += candidate.step().gross_gas_work_j;
        gas = candidate.into_parts().0.remove(&1).unwrap();
    }
    (gas, work)
}

#[test]
fn air_standard_otto_cycle_converges_to_analytic_efficiency_and_closes_work() {
    let initial = gas(0.75, 0.0, 300.0, 100000.0, 0.001);
    let cv = initial.state().state().specific_gas_constant * 2.25;
    let ratio = 3.0_f64;
    let expected = 1.0 - ratio.powf(-1.0 / 2.25);
    let mut errors = Vec::new();
    for count in [16, 32, 64] {
        let (compressed, compression_work) = march(&initial, initial.volume_m3() / ratio, count);
        let q_in = initial.mass_kg() * cv * (850.0 - compressed.state().state().temperature);
        let heated = single(&compressed, drive(compressed.volume_m3(), q_in))
            .into_parts()
            .0
            .remove(&1)
            .unwrap();
        let (expanded, expansion_work) = march(&heated, initial.volume_m3(), count);
        let q_out = initial.mass_kg() * cv * (300.0 - expanded.state().state().temperature);
        let cooled = single(&expanded, drive(initial.volume_m3(), q_out));
        let work = compression_work + expansion_work;
        assert!(compression_work < 0.0 && expansion_work > 0.0 && work > 0.0);
        near(work, q_in + q_out, 2e-7);
        near(
            cooled.parcels()[&1].state().state().temperature,
            300.0,
            1e-8,
        );
        near(
            cooled.parcels()[&1].state().state().pressure,
            100000.0,
            1e-5,
        );
        errors.push((work / q_in - expected).abs());
    }
    assert!(
        errors[0] / errors[1] > 3.7 && errors[1] / errors[2] > 3.7,
        "{errors:?}"
    );
    assert!(errors[2] < 3e-5, "{errors:?}");
}

#[test]
fn nonlinear_caloric_compression_converges_to_integrated_entropy_oracle() {
    let initial = gas(0.75, 0.002, 400.0, 100000.0, 0.001);
    // cp/R = 3.25 + 0.0015*T. Independently integrate cv/T at fixed composition.
    let (mut lo, mut hi) = (400.0_f64, 800.0_f64);
    for _ in 0..60 {
        let mid = 0.5 * (lo + hi);
        let entropy = 2.25 * (mid / 400.0).ln() + 0.0015 * (mid - 400.0) + 0.5_f64.ln();
        if entropy > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    let oracle = 0.5 * (lo + hi);
    let (coarse, _) = march(&initial, 0.0005, 16);
    let (fine, work) = march(&initial, 0.0005, 32);
    let a = (coarse.state().state().temperature - oracle).abs();
    let b = (fine.state().state().temperature - oracle).abs();
    assert!(a / b > 3.8 && b < 0.03, "{a}, {b}");
    let t = fine.state().state().temperature;
    let r = initial.state().state().specific_gas_constant;
    near(
        initial.mass_kg() * r * (2.25 * (t - 400.0) + 0.00075 * (t * t - 400.0 * 400.0)) + work,
        0.0,
        2e-7,
    );
    let (constant, _) = march(&gas(0.75, 0.0, 400.0, 100000.0, 0.001), 0.0005, 32);
    assert!(constant.state().state().temperature - t > 20.0);
}

fn open_network() -> (
    BTreeMap<u64, ClosedGasParcel>,
    Vec<SnapshotNozzleLink>,
    BTreeMap<u64, GasEnergySolve>,
    BTreeMap<u64, PrescribedGasDrive>,
) {
    let gases = BTreeMap::from([
        (1, gas(0.75, 0.0, 400.0, 180000.0, 0.002)),
        (2, gas(0.75, 0.0, 400.0, 100000.0, 0.001)),
        (3, gas(0.75, 0.0, 400.0, 60000.0, 0.002)),
    ]);
    let links = [[1, 2], [2, 3]]
        .into_iter()
        .enumerate()
        .map(|(id, parcels)| SnapshotNozzleLink {
            id: id as u64,
            parcels,
            nozzle: GasNozzle {
                area_m2: 1e-5,
                discharge_coefficient: 0.8,
            },
            nozzle_solve: energy_solve(),
        })
        .collect();
    let endpoints = gases.keys().map(|&id| (id, energy_solve())).collect();
    let mut drives: BTreeMap<_, _> = gases
        .iter()
        .map(|(&id, g)| (id, drive(g.volume_m3(), 0.0)))
        .collect();
    drives.insert(2, drive(0.00099, 0.3));
    (gases, links, endpoints, drives)
}

#[test]
fn simultaneous_inlet_outlet_and_motion_pay_upstream_enthalpy_once() {
    let (gases, links, endpoints, drives) = open_network();
    let actual =
        step_candidate(&gases, &links, 1e-5, &endpoints, &drives, solve(), || false).unwrap();
    let dm_in = actual.step().flow.links[0].transferred_mass_kg;
    let dm_out = actual.step().flow.links[1].transferred_mass_kg;
    assert!(dm_in > 0.0 && dm_out > 0.0);
    let old = &gases[&2];
    let (r, cv) = (
        old.state().state().specific_gas_constant,
        old.state().state().specific_gas_constant * 2.25,
    );
    let mass = old.mass_kg() + dm_in - dm_out;
    let inflow_enthalpy = (dm_in - dm_out) * (cv + r) * 400.0;
    let t_flow = (old.mass_kg() * cv * 400.0 + inflow_enthalpy) / (mass * cv);
    let volume = drives[&2].final_volume_m3;
    let ratio = (volume - old.volume_m3()) / (0.5 * (volume + old.volume_m3()));
    let tm = (cv * t_flow + 0.3 / (2.0 * mass)) / (cv + 0.5 * ratio * r);
    near(actual.parcels()[&2].mass_kg(), mass, 1e-16);
    near(
        actual.parcels()[&2].state().state().temperature,
        2.0 * tm - t_flow,
        1e-8,
    );
    near(
        actual.step().parcels[&2].net_enthalpy_into_gas_j,
        inflow_enthalpy,
        1e-10,
    );
    assert!(actual.parcels()[&1].mass_kg() < gases[&1].mass_kg());
    assert!(actual.parcels()[&1].state().state().temperature < 400.0);
    assert!(actual.parcels()[&3].state().state().temperature > 400.0);
    let delta_u: f64 = actual
        .parcels()
        .iter()
        .map(|(id, g)| {
            cv * (g.mass_kg() * g.state().state().temperature - gases[id].mass_kg() * 400.0)
        })
        .sum();
    near(delta_u + actual.step().gross_gas_work_j, 0.3, 2e-8);
    near(
        actual.step().net_drive_work_j + actual.step().exterior_pressure_work_j,
        actual.step().gross_gas_work_j,
        1e-12,
    );
    // Exterior effort partitions work, without changing any thermodynamic state.
    let mut vacuum = drives.clone();
    for d in vacuum.values_mut() {
        d.exterior_pressure_pa = 0.0;
    }
    let alternate =
        step_candidate(&gases, &links, 1e-5, &endpoints, &vacuum, solve(), || false).unwrap();
    assert_eq!(alternate.parcels(), actual.parcels());
    near(
        alternate.step().net_drive_work_j,
        actual.step().gross_gas_work_j,
        1e-12,
    );
}

fn species_inventory(gases: &BTreeMap<u64, ClosedGasParcel>) -> BTreeMap<String, f64> {
    let mut totals = BTreeMap::new();
    for gas in gases.values() {
        let thermo = gas.model().thermodynamics();
        let molar_mass: f64 = thermo
            .components()
            .iter()
            .zip(thermo.mole_composition().fractions())
            .map(|(s, x)| s.molar_mass().value() * x)
            .sum();
        for (s, x) in thermo
            .components()
            .iter()
            .zip(thermo.mole_composition().fractions())
        {
            *totals.entry(s.species().as_str().to_owned()).or_insert(0.0) +=
                gas.mass_kg() * x * s.molar_mass().value() / molar_mass;
        }
    }
    totals
}

#[test]
fn mixing_preserves_each_species_and_link_order_does_not_select_physics() {
    let (mut gases, links, endpoints, drives) = open_network();
    gases.insert(1, gas(0.25, 0.0, 400.0, 180000.0, 0.002));
    let before = species_inventory(&gases);
    let actual =
        step_candidate(&gases, &links, 1e-5, &endpoints, &drives, solve(), || false).unwrap();
    let after = species_inventory(actual.parcels());
    for (id, mass) in before {
        near(after[&id], mass, 1e-14);
    }
    assert_ne!(
        actual.parcels()[&2].state().state().specific_gas_constant,
        gases[&2].state().state().specific_gas_constant
    );
    // Outlet carries the accepted chamber mixture, never the just-arrived material.
    assert_eq!(
        actual.step().flow.links[1].flow.upstream_enthalpy_j_kg,
        gases[&2].state().specific_enthalpy_j_kg()
    );
    let reversed: Vec<_> = links
        .iter()
        .rev()
        .enumerate()
        .map(|(id, l)| SnapshotNozzleLink {
            id: 100 + id as u64,
            parcels: [l.parcels[1], l.parcels[0]],
            ..*l
        })
        .collect();
    let alternate = step_candidate(
        &gases,
        &reversed,
        1e-5,
        &endpoints,
        &drives,
        solve(),
        || false,
    )
    .unwrap();
    assert_eq!(actual.parcels(), alternate.parcels());
    assert_eq!(actual.step().parcels, alternate.step().parcels);
}

#[test]
fn idle_and_held_volume_paths_preserve_original_owner_behavior() {
    let original = gas(0.75, 0.001, 400.0, 100000.0, 0.001);
    assert_eq!(single(&original, drive(0.001, 0.0)).parcels()[&1], original);
    let mut expected = original.clone();
    expected.add_heat(12.0, energy_solve(), || false).unwrap();
    let actual = single(&original, drive(0.001, 12.0));
    near(
        actual.parcels()[&1].state().state().temperature,
        expected.state().state().temperature,
        1e-8,
    );
    assert_eq!(actual.step().gross_gas_work_j, 0.0);
}

#[test]
fn late_caloric_failure_pressure_crossing_and_cancellation_publish_no_parcels() {
    let (gases, links, endpoints, drives) = open_network();
    let saved = gases.clone();
    let mut bad = drives.clone();
    bad.get_mut(&3).unwrap().heat_j = 1e12;
    assert!(matches!(
        step_candidate(&gases, &links, 1e-5, &endpoints, &bad, solve(), || false),
        Err(Error::Parcel { id: 3, .. })
    ));
    bad = drives.clone();
    bad.remove(&3);
    assert!(matches!(
        step_candidate(&gases, &links, 1e-5, &endpoints, &bad, solve(), || false),
        Err(Error::Invalid(_))
    ));
    bad = drives.clone();
    bad.get_mut(&2).unwrap().exterior_pressure_pa = -1.0;
    assert!(step_candidate(&gases, &links, 1e-5, &endpoints, &bad, solve(), || false).is_err());
    assert!(step_candidate(&gases, &links, 1e3, &endpoints, &drives, solve(), || false).is_err());
    let mut total = 0;
    step_candidate(&gases, &links, 1e-5, &endpoints, &drives, solve(), || {
        total += 1;
        false
    })
    .unwrap();
    for stop in [1, total / 2, total] {
        let mut count = 0;
        assert_eq!(
            step_candidate(&gases, &links, 1e-5, &endpoints, &drives, solve(), || {
                count += 1;
                count == stop
            }),
            Err(Error::Cancelled)
        );
    }
    assert_eq!(gases, saved);
}
