use super::*;
use fs_material::thermochemical_gas::energy::parcel::transfer::species::network::{
    SpeciesNetworkAdvanceControls, transport::GasTransportSolve,
};

fn physical() -> GasTransportSolve {
    GasTransportSolve {
        exchange: controls(),
        entropy_tolerance_j_k: 1e-9,
    }
}

fn temporal(scale: f64) -> SpeciesNetworkAdvanceControls {
    SpeciesNetworkAdvanceControls {
        maximum_step_s: 4.0,
        minimum_step_s: 1e-7,
        maximum_evaluations: 20_000,
        absolute_species_tolerance_kg: scale * 1e-9,
        relative_species_tolerance: scale * 1e-4,
        temperature_tolerance_k: scale * 0.01,
    }
}

fn heated_flow(domain: &VoxelDomain, species: [f64; 2]) -> FvFlow {
    let mut field = flow(domain, &species);
    field
        .heat_eddy_viscosity_m2_s
        .copy_from_slice(&[0.0001, 0.0002]);
    field.turbulent_prandtl = 0.9;
    field
}

fn gas_constant(water: f64) -> f64 {
    R * ((1.0 - water) / MOLAR_MASS[0] + water / MOLAR_MASS[1])
}

fn caloric(t: f64) -> f64 {
    2.0 * t + 0.0005 * t * t
}

fn temperature(target: f64) -> f64 {
    2.0 * target / (2.0 + (4.0 + 0.002 * target).sqrt())
}

fn internal_energy(parcel: &ClosedGasParcel) -> f64 {
    parcel.mass_kg()
        * gas_constant(fractions(parcel)[1])
        * caloric(parcel.state().state().temperature)
}

fn entropy(parcels: &BTreeMap<u64, ClosedGasParcel>) -> f64 {
    parcels
        .values()
        .map(|p| {
            let state = p.state().state();
            let moles =
                core::array::from_fn::<_, 2, _>(|k| p.mass_kg() * fractions(p)[k] / MOLAR_MASS[k]);
            let total = moles.iter().sum::<f64>();
            moles
                .iter()
                .map(|n| {
                    R * n
                        * (3.0 * state.temperature.ln() + 0.001 * state.temperature
                            - (n / total * state.pressure / 100_000.0).ln())
                })
                .sum::<f64>()
        })
        .sum()
}

#[test]
fn heat_only_snapshot_integrates_original_variable_cp_and_preserves_sources() {
    let domain = VoxelDomain::graded([vec![0.05, 0.15], vec![0.1], vec![0.1]]).unwrap();
    let field = heated_flow(&domain, [0.0; 2]);
    for reference in [0.0, 1e18] {
        let parcels = BTreeMap::from([
            (0, parcel(&domain, 0, 0.25, 500.0, reference)),
            (1, parcel(&domain, 1, 0.25, 350.0, reference)),
        ]);
        let result = field
            .diffuse_turbulent_gas_transport(
                &domain,
                &carrier(&domain, &parcels),
                &parcels,
                0.7,
                0.2,
                &endpoints(&parcels),
                controls(),
                &CancelGate::new_clock_free(),
            )
            .unwrap();
        assert!(result.exchange.step().fluxes.is_empty());
        assert_eq!(result.exchange.step().heat.len(), 1);
        let conductivity: [f64; 2] = core::array::from_fn(|i| {
            let p = &parcels[&(i as u64)];
            p.mass_kg() / p.volume_m3()
                * gas_constant(fractions(p)[1])
                * (3.0 + 0.001 * p.state().state().temperature)
                * field.heat_eddy_viscosity_m2_s[i]
                / field.turbulent_prandtl
        });
        let heat = 0.2 * 0.01 / (0.025 / conductivity[0] + 0.075 / conductivity[1]) * 150.0;
        let edge = &result.exchange.step().heat[0];
        assert_eq!(edge.parcels, [0, 1]);
        near(edge.heat_left_to_right_j, heat, 2e-14);
        let mut entropy_change = 0.0;
        for (id, sign) in [(0, -1.0), (1, 1.0)] {
            let old = &parcels[&id];
            let new = &result.exchange.parcels()[&id];
            let old_t = old.state().state().temperature;
            let capacity = old.mass_kg() * gas_constant(fractions(old)[1]);
            let expected = temperature(caloric(old_t) + sign * heat / capacity);
            near(new.state().state().temperature, expected, 2e-7);
            near(
                old.stored_energy_change_to(new, || false).unwrap(),
                sign * heat,
                1e-7,
            );
            near(
                new.state().state().pressure,
                100_000.0 * expected / old_t,
                1e-4,
            );
            assert_eq!(new.mass_kg(), old.mass_kg());
            assert_eq!(new.volume_m3(), old.volume_m3());
            assert_eq!(new.model(), old.model());
            assert_eq!(new.model().source_packs(), old.model().source_packs());
            assert_eq!(
                new.point().axes().get("cell-lot"),
                old.point().axes().get("cell-lot")
            );
            entropy_change +=
                capacity * (2.0 * (expected / old_t).ln() + 0.001 * (expected - old_t));
        }
        assert!(entropy_change > 0.0);
        assert!(
            result.exchange.step().total_energy_residual_j.abs() <= controls().energy_tolerance_j
        );
    }
}

#[test]
fn adaptive_heat_retries_source_bracket_overshoot_and_completes_the_window() {
    let domain = VoxelDomain::graded([vec![0.05, 0.15], vec![0.1], vec![0.1]]).unwrap();
    let parcels = BTreeMap::from([
        (0, parcel(&domain, 0, 0.25, 500.0, 0.0)),
        (1, parcel(&domain, 1, 0.25, 350.0, 0.0)),
    ]);
    let accepted = carrier(&domain, &parcels);
    let saved = (parcels.clone(), accepted.clone());
    let field = heated_flow(&domain, [0.0; 2]);
    let solvers = endpoints(&parcels);
    let gate = CancelGate::new_clock_free();
    // The frozen law sends 222.54 J in 100 s. Original hot-cell Cv permits
    // only 70.5 J before the source's 200 K lower bracket is reached.
    let refused = field
        .diffuse_turbulent_gas_transport(
            &domain,
            &accepted,
            &parcels,
            0.7,
            100.0,
            &solvers,
            controls(),
            &gate,
        )
        .unwrap_err();
    assert!(
        matches!(refused, fs_lbm::conjugate::ChtError::InvalidInput { reason, .. }
        if reason.contains("OutsideBracket"))
    );
    let settings = SpeciesNetworkAdvanceControls {
        maximum_step_s: 100.0,
        maximum_evaluations: 1000,
        ..temporal(2.0)
    };
    let result = field
        .advance_turbulent_gas_transport(
            &domain,
            &accepted,
            &parcels,
            0.7,
            100.0,
            &solvers,
            physical(),
            settings,
            &gate,
        )
        .unwrap();
    let report = result.exchange.report();
    assert_eq!(report.duration_s, 100.0);
    assert!(report.rejected_intervals >= 2);
    assert!(report.accepted_intervals > 1);
    assert!(report.transfers.is_empty());
    assert_eq!(report.heat.len(), 1);
    let mut energy_change = 0.0;
    for id in [0, 1] {
        let old = &parcels[&id];
        let new = &result.exchange.parcels()[&id];
        assert!((350.0..=500.0).contains(&new.state().state().temperature));
        assert_eq!(new.mass_kg(), old.mass_kg());
        assert_eq!(new.model(), old.model());
        let change = old.stored_energy_change_to(new, || false).unwrap();
        let sign = if id == 0 { -1.0 } else { 1.0 };
        near(
            change,
            sign * report.heat[0].heat_left_to_right_j,
            controls().energy_tolerance_j,
        );
        energy_change += change;
    }
    near(energy_change, 0.0, controls().energy_tolerance_j);
    assert!(report.endpoint_energy_change_j.abs() <= controls().energy_tolerance_j);
    assert!(report.total_energy_residual_j.abs() <= controls().energy_tolerance_j);
    assert!(report.entropy_change_j_k > 0.0);
    near(
        report.entropy_change_j_k,
        entropy(result.exchange.parcels()) - entropy(&parcels),
        2e-10,
    );
    assert_eq!((parcels, accepted), saved);
}

// Independent classical RK4 of the authored ideal-gas laws. Coordinates are
// left water mass, each cell's internal energy, and integrated additional heat.
// No production transport, caloric inversion or adaptive controller is called.
fn rk4(parcels: &BTreeMap<u64, ClosedGasParcel>, duration: f64) -> ([f64; 2], f64) {
    let mass = [parcels[&0].mass_kg(), parcels[&1].mass_kg()];
    let rho = [mass[0] / 0.0005, mass[1] / 0.0015];
    let water = mass[0] * fractions(&parcels[&0])[1] + mass[1] * fractions(&parcels[&1])[1];
    let conductance = 0.01 / (0.025 / (rho[0] * 0.0001 / 0.7) + 0.075 / (rho[1] * 0.0002 / 0.7));
    let unpack = |s: [f64; 4]| {
        let y = [s[0] / mass[0], (water - s[0]) / mass[1]];
        let t = core::array::from_fn::<_, 2, _>(|i| {
            temperature(s[i + 1] / (mass[i] * gas_constant(y[i])))
        });
        (y, t)
    };
    let rhs = |s: [f64; 4]| {
        let (y, t) = unpack(s);
        let water_rate = conductance * (y[0] - y[1]);
        let rates = [-water_rate, water_rate];
        let enthalpy = (0..2)
            .map(|k| {
                let donor_t = t[usize::from(rates[k] < 0.0)];
                rates[k] * R / MOLAR_MASS[k] * (3.0 * donor_t + 0.0005 * donor_t * donor_t)
            })
            .sum::<f64>();
        let k: [f64; 2] = core::array::from_fn(|i| {
            rho[i] * gas_constant(y[i]) * (3.0 + 0.001 * t[i]) * [0.0001, 0.0002][i] / 0.9
        });
        let heat = 0.01 / (0.025 / k[0] + 0.075 / k[1]) * (t[0] - t[1]);
        [-water_rate, -enthalpy - heat, enthalpy + heat, heat]
    };
    let mut state = [
        mass[0] * fractions(&parcels[&0])[1],
        internal_energy(&parcels[&0]),
        internal_energy(&parcels[&1]),
        0.0,
    ];
    let dt = duration / 1024.0;
    for _ in 0..1024 {
        let a = rhs(state);
        let b = rhs(core::array::from_fn(|i| state[i] + 0.5 * dt * a[i]));
        let c = rhs(core::array::from_fn(|i| state[i] + 0.5 * dt * b[i]));
        let d = rhs(core::array::from_fn(|i| state[i] + dt * c[i]));
        state =
            core::array::from_fn(|i| state[i] + dt / 6.0 * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]));
    }
    (unpack(state).1, state[3])
}

#[test]
fn combined_adaptive_transport_converges_and_ledgers_only_accepted_heat() {
    let domain = VoxelDomain::graded([vec![0.05, 0.15], vec![0.1], vec![0.1]]).unwrap();
    let parcels = pair(&domain, 0.0);
    let accepted = carrier(&domain, &parcels);
    let saved = (parcels.clone(), accepted.clone());
    let field = heated_flow(&domain, [0.0001, 0.0002]);
    let solve = |scale| {
        field
            .advance_turbulent_gas_transport(
                &domain,
                &accepted,
                &parcels,
                0.7,
                4.0,
                &endpoints(&parcels),
                physical(),
                temporal(scale),
                &CancelGate::new_clock_free(),
            )
            .unwrap()
    };
    let coarse = solve(1.0);
    let fine = solve(0.01);
    let (expected_t, expected_heat) = rk4(&parcels, 4.0);
    let error = |state: &BTreeMap<u64, ClosedGasParcel>| {
        (0..2)
            .map(|i| (state[&(i as u64)].state().state().temperature - expected_t[i]).abs())
            .fold(0.0_f64, f64::max)
    };
    assert!(error(fine.exchange.parcels()) < 0.25 * error(coarse.exchange.parcels()));
    assert!(error(fine.exchange.parcels()) < 0.05);
    let report = fine.exchange.report();
    assert_eq!(report.duration_s, 4.0);
    assert_eq!(report.accepted_substeps, 2 * report.accepted_intervals);
    assert!(report.rejected_intervals > 0);
    assert!(report.evaluations > report.accepted_substeps);
    assert!(
        report.maximum_species_error_ratio <= 1.0 && report.maximum_temperature_error_ratio <= 1.0
    );
    assert_eq!(report.heat.len(), 1);
    assert_eq!(report.heat[0].parcels, [0, 1]);
    let heat = report.heat[0].heat_left_to_right_j;
    near(heat, expected_heat, 0.01);
    let enthalpy = report
        .transfers
        .iter()
        .map(|row| {
            if row.donor == 0 {
                row.carried_enthalpy_j
            } else {
                -row.carried_enthalpy_j
            }
        })
        .sum::<f64>();
    let endpoint = fine.exchange.parcels();
    near(
        internal_energy(&endpoint[&0]) - internal_energy(&parcels[&0]) + enthalpy + heat,
        0.0,
        2e-7,
    );
    near(
        endpoint.values().map(internal_energy).sum::<f64>()
            - parcels.values().map(internal_energy).sum::<f64>(),
        0.0,
        2e-7,
    );
    near(
        report.entropy_change_j_k,
        entropy(endpoint) - entropy(&parcels),
        2e-10,
    );
    assert!(report.entropy_change_j_k > 0.0);
    assert!(report.endpoint_energy_change_j.abs() <= physical().exchange.energy_tolerance_j);
    assert!(report.total_energy_residual_j.abs() <= physical().exchange.energy_tolerance_j);
    assert!(fine.maximum_carrier_mass_defect_kg <= controls().mass_tolerance_kg);
    assert_eq!((parcels, accepted), saved);
}

#[test]
fn prandtl_and_schmidt_remain_independent_with_atomic_heat_and_work_refusals() {
    let domain = VoxelDomain::graded([vec![0.05, 0.15], vec![0.1], vec![0.1]]).unwrap();
    let parcels = pair(&domain, 0.0);
    let accepted = carrier(&domain, &parcels);
    let saved = (parcels.clone(), accepted.clone());
    let field = heated_flow(&domain, [0.0001, 0.0002]);
    let gate = CancelGate::new_clock_free();
    let solvers = endpoints(&parcels);
    let snapshot = |f: &FvFlow, sc, dt, gate: &CancelGate| {
        f.diffuse_turbulent_gas_transport(
            &domain,
            &accepted,
            &parcels,
            sc,
            dt,
            &solvers,
            controls(),
            gate,
        )
    };
    let baseline = snapshot(&field, 0.7, 0.2, &gate).unwrap();
    let mut changed = field.clone();
    changed.turbulent_prandtl *= 2.0;
    let prandtl = snapshot(&changed, 0.7, 0.2, &gate).unwrap();
    assert_eq!(
        prandtl.exchange.step().fluxes,
        baseline.exchange.step().fluxes
    );
    near(
        prandtl.exchange.step().heat[0].heat_left_to_right_j * 2.0,
        baseline.exchange.step().heat[0].heat_left_to_right_j,
        1e-14,
    );
    let schmidt = snapshot(&field, 1.4, 0.2, &gate).unwrap();
    assert_eq!(schmidt.exchange.step().heat, baseline.exchange.step().heat);
    for (a, b) in schmidt
        .exchange
        .step()
        .fluxes
        .iter()
        .zip(&baseline.exchange.step().fluxes)
    {
        assert_eq!(a.transfer.species, b.transfer.species);
        near(a.transfer.mass_kg * 2.0, b.transfer.mass_kg, 1e-18);
        assert_eq!(a.source_enthalpy_j_kg, b.source_enthalpy_j_kg);
    }
    changed.heat_eddy_viscosity_m2_s.fill(0.0);
    let species_only = field
        .diffuse_turbulent_gas_species(
            &domain,
            &accepted,
            &parcels,
            0.7,
            0.2,
            &solvers,
            controls(),
            &gate,
        )
        .unwrap();
    assert_eq!(snapshot(&changed, 0.7, 0.2, &gate).unwrap(), species_only);
    let isothermal = BTreeMap::from([
        (0, parcel(&domain, 0, 0.1, 400.0, 0.0)),
        (1, parcel(&domain, 1, 0.6, 400.0, 0.0)),
    ]);
    let mixed = field
        .diffuse_turbulent_gas_transport(
            &domain,
            &carrier(&domain, &isothermal),
            &isothermal,
            0.7,
            0.2,
            &endpoints(&isothermal),
            controls(),
            &gate,
        )
        .unwrap();
    assert!(mixed.exchange.step().heat.is_empty());
    assert!(!mixed.exchange.step().fluxes.is_empty());
    for pr in [0.0, f64::NAN] {
        changed.turbulent_prandtl = pr;
        assert!(snapshot(&changed, 0.7, 0.2, &gate).is_err());
    }
    changed = field.clone();
    changed.heat_eddy_viscosity_m2_s[0] = -0.01;
    assert!(snapshot(&changed, 0.7, 0.2, &gate).is_err());
    assert!(snapshot(&heated_flow(&domain, [0.0; 2]), 0.7, 1e6, &gate).is_err());
    let mut limited = temporal(0.01);
    limited.maximum_evaluations = 3;
    assert!(
        field
            .advance_turbulent_gas_transport(
                &domain,
                &accepted,
                &parcels,
                0.7,
                4.0,
                &solvers,
                physical(),
                limited,
                &gate,
            )
            .is_err()
    );
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(snapshot(&field, 0.7, 0.2, &cancelled).is_err());
    assert!(
        field
            .advance_turbulent_gas_transport(
                &domain,
                &accepted,
                &parcels,
                0.7,
                4.0,
                &solvers,
                physical(),
                temporal(0.01),
                &cancelled,
            )
            .is_err()
    );
    assert_eq!((parcels, accepted), saved);
}
