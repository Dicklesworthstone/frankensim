use super::*;
use fs_material::thermochemical_gas::energy::parcel::transfer::species::network::{
    SpeciesNetworkAdvanceControls, transport::GasTransportSolve,
};

// Independent integral and inverse of this fixture's Cv/R = 2 + 0.001 T.
fn caloric(t: f64) -> f64 {
    2.0 * t + 0.0005 * t * t
}

fn temperature(u: f64) -> f64 {
    2.0 * u / (2.0 + (4.0 + 0.002 * u).sqrt())
}

fn capacity(parcel: &ClosedGasParcel) -> f64 {
    parcel.mass_kg()
        * fractions(parcel)
            .into_iter()
            .zip(MOLAR_MASS)
            .map(|(y, molar_mass)| y * R / molar_mass)
            .sum::<f64>()
}

fn conductance(domain: &VoxelDomain, k: [f64; 2]) -> f64 {
    domain.face_area(0, 0) / (0.5 * domain.widths(0)[0] / k[0] + 0.5 * domain.widths(1)[0] / k[1])
}

#[test]
fn zero_sgs_retains_molecular_heat_variable_cp_and_original_source_references() {
    let domain = VoxelDomain::graded([vec![0.05, 0.15], vec![0.1], vec![0.1]]).unwrap();
    let field = flow(&domain, &[0.0; 2]);
    let gate = CancelGate::new_clock_free();
    let mut unshifted_temperature = [0.0; 2];
    for reference in [0.0, 1e18] {
        let parcels = pair(&domain, reference);
        let accepted = carrier(&domain, &parcels);
        let solvers = endpoints(&parcels);
        let old_component = field
            .diffuse_turbulent_gas_transport(
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
        assert_eq!(old_component.exchange.parcels(), &parcels);
        assert!(old_component.exchange.step().heat.is_empty());
        let result = field
            .diffuse_turbulent_gas_transport_with_molecular_heat(
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
        // There is a composition gradient, but this API adds no molecular
        // species model when the original momentum eddy field is zero.
        assert!(result.exchange.step().fluxes.is_empty());
        assert_eq!(result.exchange.step().heat.len(), 1);
        let k = core::array::from_fn(|i| parcels[&(i as u64)].state().state().thermal_conductivity);
        let heat = 0.2 * conductance(&domain, k) * 150.0;
        let edge = &result.exchange.step().heat[0];
        assert_eq!(edge.parcels, [0, 1]);
        near(edge.heat_left_to_right_j, heat, 2e-14);
        let mut energy_change = 0.0;
        for (id, sign) in [(0, -1.0), (1, 1.0)] {
            let old = &parcels[&id];
            let new = &result.exchange.parcels()[&id];
            let initial_t = old.state().state().temperature;
            let expected = temperature(caloric(initial_t) + sign * heat / capacity(old));
            near(new.state().state().temperature, expected, 2e-7);
            let change = old.stored_energy_change_to(new, || false).unwrap();
            near(change, sign * heat, controls().energy_tolerance_j);
            energy_change += change;
            assert_eq!(new.model(), old.model());
            assert_eq!(new.model().source_packs(), old.model().source_packs());
            assert_eq!(new.mass_kg(), old.mass_kg());
            assert_eq!(new.volume_m3(), old.volume_m3());
            assert_eq!(new.point().axes().get("cell-lot"), Some(&(id as f64)));
            near(
                new.state().state().pressure,
                100_000.0 * expected / initial_t,
                1e-4,
            );
            if reference == 0.0 {
                unshifted_temperature[id as usize] = new.state().state().temperature;
            } else {
                near(
                    new.state().state().temperature,
                    unshifted_temperature[id as usize],
                    2e-7,
                );
            }
        }
        near(energy_change, 0.0, controls().energy_tolerance_j);
        assert!(
            result.exchange.step().total_energy_residual_j.abs() <= controls().energy_tolerance_j
        );
    }
    let uniform = BTreeMap::from([
        (0, parcel(&domain, 0, 0.1, 400.0, 0.0)),
        (1, parcel(&domain, 1, 0.6, 400.0, 0.0)),
    ]);
    let identity = field
        .diffuse_turbulent_gas_transport_with_molecular_heat(
            &domain,
            &carrier(&domain, &uniform),
            &uniform,
            0.7,
            0.2,
            &endpoints(&uniform),
            controls(),
            &gate,
        )
        .unwrap();
    assert_eq!(identity.exchange.parcels(), &uniform);
    assert!(identity.exchange.step().heat.is_empty());
}

#[test]
fn source_and_eddy_conductivities_combine_before_the_half_cell_resistances() {
    let domain = VoxelDomain::graded([vec![0.05, 0.15], vec![0.1], vec![0.1]]).unwrap();
    let parcels = pair(&domain, 0.0);
    let accepted = carrier(&domain, &parcels);
    let mut field = flow(&domain, &[0.0; 2]);
    // A zero-SGS cell closes the old turbulent-only face. The original
    // molecular source still transmits heat through that same half cell.
    field
        .heat_eddy_viscosity_m2_s
        .copy_from_slice(&[0.0, 0.0003]);
    field.turbulent_prandtl = 0.9;
    let solvers = endpoints(&parcels);
    let gate = CancelGate::new_clock_free();
    let old_component = field
        .diffuse_turbulent_gas_transport(
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
    assert!(old_component.exchange.step().heat.is_empty());
    assert_eq!(old_component.exchange.parcels(), &parcels);
    let molecular: [f64; 2] =
        core::array::from_fn(|i| parcels[&(i as u64)].state().state().thermal_conductivity);
    let total = core::array::from_fn(|i| {
        let state = parcels[&(i as u64)].state().state();
        molecular[i]
            + state.density * state.specific_heat_cp * field.heat_eddy_viscosity_m2_s[i]
                / field.turbulent_prandtl
    });
    let expected = 0.2 * conductance(&domain, total) * 150.0;
    // Adding separate face conductances would add zero turbulent conductance
    // to this molecular value. It misses the reduced right-cell resistance.
    let incorrectly_separate = 0.2 * conductance(&domain, molecular) * 150.0;
    assert!(expected > 1.5 * incorrectly_separate);
    let result = field
        .diffuse_turbulent_gas_transport_with_molecular_heat(
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
    assert!(result.exchange.step().fluxes.is_empty());
    assert_eq!(result.exchange.step().heat.len(), 1);
    near(
        result.exchange.step().heat[0].heat_left_to_right_j,
        expected,
        2e-14,
    );
    for (id, sign) in [(0, -1.0), (1, 1.0)] {
        near(
            parcels[&id]
                .stored_energy_change_to(&result.exchange.parcels()[&id], || false)
                .unwrap(),
            sign * expected,
            controls().energy_tolerance_j,
        );
    }
}

// Scalar RK4 of heat exchange with the authored pressure-linear source law.
// For fixed composition all species share its pressure multiplier, so Wilke
// viscosity ratios and the WMS mixture factor stay constant. No production
// caloric solve, transport evaluator or adaptive controller is used here.
fn reference_heat(
    domain: &VoxelDomain,
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    duration: f64,
    refresh: bool,
) -> ([f64; 2], f64) {
    let original: [_; 2] = core::array::from_fn(|i| parcels[&(i as u64)].state().state());
    let capacities: [_; 2] = core::array::from_fn(|i| capacity(&parcels[&(i as u64)]));
    let curve = |p: f64| 0.5 + (p - 1000.0) / 249_000.0;
    let unpack = |heat: f64| {
        core::array::from_fn::<_, 2, _>(|i| {
            temperature(caloric(original[i].temperature) + [-1.0, 1.0][i] * heat / capacities[i])
        })
    };
    let rhs = |heat: f64| {
        let t = unpack(heat);
        let k = core::array::from_fn(|i| {
            let factor = if refresh {
                curve(original[i].pressure * t[i] / original[i].temperature)
                    / curve(original[i].pressure)
            } else {
                1.0
            };
            original[i].thermal_conductivity * factor
        });
        conductance(domain, k) * (t[0] - t[1])
    };
    let dt = duration / 1024.0;
    let mut heat = 0.0;
    for _ in 0..1024 {
        let a = rhs(heat);
        let b = rhs(heat + 0.5 * dt * a);
        let c = rhs(heat + 0.5 * dt * b);
        let d = rhs(heat + dt * c);
        heat += dt / 6.0 * (a + 2.0 * b + 2.0 * c + d);
    }
    (unpack(heat), heat)
}

#[test]
fn adaptive_molecular_heat_refreshes_source_conductivity_and_converges_with_fixed_budgets() {
    let domain = VoxelDomain::graded([vec![0.05, 0.15], vec![0.1], vec![0.1]]).unwrap();
    // Unequal original pressures make source-conductivity refresh observable.
    // Both sources remain within their authored pressure/temperature support.
    let parcels = [(0, 500.0, 60_000.0), (1, 350.0, 180_000.0)]
        .into_iter()
        .map(|(id, t, p)| {
            (
                id,
                ClosedGasParcel::new(
                    source::moist_model(0.25, 0.0),
                    source::point(t, p),
                    domain.volume(id as usize),
                )
                .unwrap(),
            )
        })
        .collect::<BTreeMap<_, _>>();
    let accepted = carrier(&domain, &parcels);
    let saved = (parcels.clone(), accepted.clone());
    let field = flow(&domain, &[0.0; 2]);
    let physical = GasTransportSolve {
        exchange: controls(),
        entropy_tolerance_j_k: 1e-9,
    };
    let solve = |tolerance| {
        field
            .advance_turbulent_gas_transport_with_molecular_heat(
                &domain,
                &accepted,
                &parcels,
                0.7,
                40.0,
                &endpoints(&parcels),
                physical,
                SpeciesNetworkAdvanceControls {
                    maximum_step_s: 40.0,
                    minimum_step_s: 1e-7,
                    maximum_evaluations: 20_000,
                    absolute_species_tolerance_kg: 1e-12,
                    relative_species_tolerance: 0.0,
                    temperature_tolerance_k: tolerance,
                },
                &CancelGate::new_clock_free(),
            )
            .unwrap()
    };
    let coarse = solve(0.02);
    let fine = solve(0.0002);
    let (expected, expected_heat) = reference_heat(&domain, &parcels, 40.0, true);
    let (frozen, _) = reference_heat(&domain, &parcels, 40.0, false);
    assert!((expected[0] - frozen[0]).abs() > 0.1);
    let error = |values: &BTreeMap<u64, ClosedGasParcel>| {
        (0..2)
            .map(|i| (values[&(i as u64)].state().state().temperature - expected[i]).abs())
            .fold(0.0_f64, f64::max)
    };
    assert!(error(fine.exchange.parcels()) < 0.2 * error(coarse.exchange.parcels()));
    assert!(error(fine.exchange.parcels()) < 0.05);
    let report = fine.exchange.report();
    assert_eq!(report.duration_s, 40.0);
    assert!(report.rejected_intervals > 0);
    assert_eq!(report.accepted_substeps, 2 * report.accepted_intervals);
    assert!(report.maximum_temperature_error_ratio <= 1.0);
    assert!(report.transfers.is_empty());
    assert_eq!(report.heat.len(), 1);
    assert_eq!(report.heat[0].parcels, [0, 1]);
    near(report.heat[0].heat_left_to_right_j, expected_heat, 0.01);
    let mut entropy_change = 0.0;
    for (id, sign) in [(0, -1.0), (1, 1.0)] {
        let old = &parcels[&id];
        let new = &fine.exchange.parcels()[&id];
        near(
            old.stored_energy_change_to(new, || false).unwrap(),
            sign * report.heat[0].heat_left_to_right_j,
            controls().energy_tolerance_j,
        );
        assert_eq!(new.model(), old.model());
        assert_eq!(new.mass_kg(), old.mass_kg());
        assert_eq!(new.volume_m3(), old.volume_m3());
        let before = old.state().state().temperature;
        let after = new.state().state().temperature;
        entropy_change += capacity(old) * (2.0 * (after / before).ln() + 0.001 * (after - before));
    }
    assert!(entropy_change > 0.0);
    near(report.entropy_change_j_k, entropy_change, 2e-10);
    assert!(report.endpoint_energy_change_j.abs() <= controls().energy_tolerance_j);
    assert!(report.total_energy_residual_j.abs() <= controls().energy_tolerance_j);
    assert_eq!((parcels, accepted), saved);
}
