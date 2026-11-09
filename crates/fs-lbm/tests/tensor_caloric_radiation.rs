#![cfg(feature = "enthalpy-coupling")]
//! G0/G1/G4: original enthalpy storage, physical radiation and reservoir heat.

use fs_blake3::ContentHash;
use fs_exec::CancelGate;
use fs_lbm::conjugate::{
    CellSink, ChtError, ConductivityTensor, EnergyConfig, RadiationConfig, STEFAN_BOLTZMANN,
    SolidMaterial, TensorCaloricLaw, TensorCaloricStorage, TensorNonlinearConfig, ThermalFace,
    ThermalSetup, Voxel, VoxelDomain, advance_caloric_tensor_conduction_radiating,
};
use fs_material::phase::{
    EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, HeatCapacityKnot, SolidLiquidPhase,
};

const EMISSIVITY: f64 = 0.7;
const MASS: f64 = 2.0;
const CP: f64 = 2.0;

fn controls() -> TensorNonlinearConfig {
    TensorNonlinearConfig {
        linear: EnergyConfig {
            tolerance: 1e-10,
            ..EnergyConfig::default()
        },
        max_newton_iterations: 256,
        max_backtracks: 256,
    }
}

fn close(actual: f64, expected: f64) {
    let tolerance = 2e-6 + 2e-9 * expected.abs();
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, tolerance {tolerance:e}",
    );
}

fn conductivity() -> ConductivityTensor {
    ConductivityTensor::new([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]).unwrap()
}

fn constant_cp(datum: f64) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_from_heat_capacity_exact(
        ContentHash([0xD1; 32]),
        SolidLiquidPhase::Solid,
        datum - 200.0,
        &[200.0, 500.0].map(|temperature_k| HeatCapacityKnot {
            temperature_k,
            specific_heat_capacity_j_kg_k: CP,
            bulk_density_kg_m3: MASS,
        }),
    )
    .unwrap()
}

fn melting() -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([0xD2; 32]),
        [
            (-200.0, 200.0, 0.0),
            (0.0, 300.0, 0.0),
            (100.0, 300.0, 1.0),
            (500.0, 500.0, 1.0),
        ]
        .map(
            |(h, temperature_k, liquid_mass_fraction)| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: h,
                temperature_k,
                liquid_mass_fraction,
                bulk_density_kg_m3: MASS,
            },
        )
        .to_vec(),
    )
    .unwrap()
}

fn exposed_pair(
    ambient_k: f64,
) -> (
    VoxelDomain,
    [SolidMaterial; 1],
    ThermalSetup,
    RadiationConfig,
) {
    // Only the x+ face of cell zero radiates: native radiation emits from
    // solid faces bordering fluid. All rays escape to one temperature, F=1,
    // so this oracle has no Monte Carlo view-factor uncertainty.
    let domain = VoxelDomain::from_fn(2, 1, 1, 1.0, |p| {
        if p[0] < 1.0 {
            Voxel::Solid(0)
        } else {
            Voxel::Fluid
        }
    })
    .unwrap();
    let solids = [SolidMaterial::new("finite radiating solid", 1.0)];
    let mut setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    setup.fixed_temperature = vec![(1, 300.0)];
    let mut radiation = RadiationConfig::new(vec![EMISSIVITY], [Some(ambient_k); 6]);
    radiation.surface_exchange = false;
    radiation.rays_per_face = 32;
    radiation.patch_size = 1;
    radiation.max_iterations = 64;
    radiation.tolerance = 1e-10;
    (domain, solids, setup, radiation)
}

fn radiant_power(temperature_k: f64, ambient_k: f64) -> f64 {
    EMISSIVITY * STEFAN_BOLTZMANN * (temperature_k.powi(4) - ambient_k.powi(4))
}

fn root(mut low: f64, mut high: f64, residual: impl Fn(f64) -> f64) -> f64 {
    assert!(residual(low) <= 0.0 && residual(high) >= 0.0);
    for _ in 0..100 {
        let middle = low + 0.5 * (high - low);
        if residual(middle) > 0.0 {
            high = middle;
        } else {
            low = middle;
        }
    }
    low + 0.5 * (high - low)
}

#[test]
fn backward_euler_matches_quartic_balance_with_authored_heat_and_energy_datums() {
    let (domain, solids, mut setup, mut radiation) = exposed_pair(300.0);
    setup.power_w = vec![17.0, 0.0];
    setup.cell_sinks.push(CellSink {
        cell: 0,
        conductance_w_k: 0.3,
        temperature_k: 310.0,
    });
    let saved_setup = setup.clone();
    let k = conductivity();
    let dt = 0.25;
    // Unit half cells in series give exactly 1 W/K to the fluid thermostat.
    let outflow = |t: f64| (t - 300.0) + 0.3 * (t - 310.0) + radiant_power(t, 300.0);
    let expected = root(300.0, 400.0, |t| {
        MASS * CP * (t - 400.0) + dt * (outflow(t) - 17.0)
    });

    for surface_exchange in [false, true] {
        radiation.surface_exchange = surface_exchange;
        for datum in [0.0, 1e6] {
            let curve = constant_cp(datum);
            let laws: [&dyn TensorCaloricLaw; 2] = [&curve, &curve];
            let previous = [datum + 200.0, datum];
            let (step, report) = advance_caloric_tensor_conduction_radiating(
                &domain,
                &solids,
                &setup,
                &[&k, &k],
                TensorCaloricStorage {
                    laws: &laws,
                    reference_mass_kg: &[MASS; 2],
                    previous_specific_energy_j_kg: &previous,
                },
                &[(1, datum)],
                5.0,
                dt,
                &controls(),
                &radiation,
                &CancelGate::new(),
            )
            .unwrap();
            let temperature = &step.endpoint.solution.solution.temperature;
            close(temperature[0], expected);
            assert_eq!(temperature[1], 300.0);
            close(
                step.endpoint.specific_energy_j_kg[0] - datum,
                CP * (expected - 300.0),
            );
            close(
                step.record.stored_energy_change_j,
                MASS * CP * (expected - 400.0),
            );
            close(step.record.source_j, 17.0 * dt);
            close(step.record.boundary_outflow_j, dt * outflow(expected));
            close(step.record.radiated_w, radiant_power(expected, 300.0));
            close(step.record.closure_j, 0.0);
            close(step.record.time_s, 5.0 + dt);
            assert_eq!(report.exposed_faces, 1);
            assert!(report.iterations > 1 && report.iterations <= radiation.max_iterations);
            assert!(report.relative_power_defect <= radiation.tolerance);
            assert!(report.temperature_change <= radiation.tolerance);
            assert_eq!(previous, [datum + 200.0, datum]);
        }
    }
    assert_eq!(setup, saved_setup);
}

#[test]
fn radiation_changes_latent_energy_at_exactly_constant_temperature() {
    let curve = melting();
    let laws: [&dyn TensorCaloricLaw; 2] = [&curve, &curve];
    let k = conductivity();
    let dt = 0.25;
    for (ambient, start_h) in [(250.0, 80.0), (350.0, 20.0)] {
        let (domain, solids, setup, radiation) = exposed_pair(ambient);
        let previous = [start_h, 0.0];
        let power = radiant_power(300.0, ambient);
        let expected_h = start_h - dt * power / MASS;
        assert!((0.0..100.0).contains(&expected_h));
        let (step, report) = advance_caloric_tensor_conduction_radiating(
            &domain,
            &solids,
            &setup,
            &[&k, &k],
            TensorCaloricStorage {
                laws: &laws,
                reference_mass_kg: &[MASS; 2],
                previous_specific_energy_j_kg: &previous,
            },
            &[(1, 0.0)],
            0.0,
            dt,
            &controls(),
            &radiation,
            &CancelGate::new(),
        )
        .unwrap();
        assert_eq!(step.endpoint.solution.solution.temperature, vec![300.0; 2]);
        let h = step.endpoint.specific_energy_j_kg[0];
        close(h, expected_h);
        let phase = curve.values_at_specific_enthalpy(h).unwrap();
        close(phase.liquid_mass_fraction, expected_h / 100.0);
        assert_eq!(curve.temperature_slope_kg_k_j(h).unwrap(), 0.0);
        close(step.record.stored_energy_change_j, -dt * power);
        close(step.record.boundary_outflow_j, dt * power);
        close(step.record.source_j, 0.0);
        close(step.record.closure_j, 0.0);
        close(report.radiated_w, power);
        for faces in &step.endpoint.solution.solution.outward_heat_w {
            for &heat in faces {
                close(heat, 0.0);
            }
        }
        assert_eq!(previous, [start_h, 0.0]);
    }
}

#[test]
fn closed_enclosure_exchanges_radiation_without_exporting_or_creating_heat() {
    let domain = VoxelDomain::from_fn(3, 1, 1, 1.0, |p| {
        if (1.0..2.0).contains(&p[0]) {
            Voxel::Fluid
        } else {
            Voxel::Solid(0)
        }
    })
    .unwrap();
    let solids = [SolidMaterial::new("closed radiating plates", 1.0)];
    let setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    let curve = constant_cp(0.0);
    let laws: [&dyn TensorCaloricLaw; 3] = [&curve, &curve, &curve];
    let previous = [200.0, 100.0, 0.0];
    let k = conductivity();
    let mut radiation = RadiationConfig::new(vec![0.8], [None; 6]);
    radiation.patch_size = 1;
    radiation.rays_per_face = 64;
    radiation.max_iterations = 128;
    radiation.tolerance = 1e-10;
    let solve = |radiation: &RadiationConfig| {
        advance_caloric_tensor_conduction_radiating(
            &domain,
            &solids,
            &setup,
            &[&k, &k, &k],
            TensorCaloricStorage {
                laws: &laws,
                reference_mass_kg: &[MASS; 3],
                previous_specific_energy_j_kg: &previous,
            },
            &[],
            0.0,
            0.1,
            &controls(),
            radiation,
            &CancelGate::new(),
        )
        .unwrap()
    };
    let (radiating, report) = solve(&radiation);
    let mut dark = radiation.clone();
    dark.emissivity[0] = 0.0;
    let (conducting, _) = solve(&dark);
    let temperature = &radiating.endpoint.solution.solution.temperature;
    // Reciprocal equal-area exchange transfers more heat than conduction alone.
    assert!(temperature[0] < conducting.endpoint.solution.solution.temperature[0]);
    assert!(temperature[2] > conducting.endpoint.solution.solution.temperature[2]);
    assert!(temperature[0] > temperature[1] && temperature[1] > temperature[2]);
    close(temperature[1], 350.0);
    close(temperature[0] + temperature[2], 700.0);
    close(
        MASS * radiating.endpoint.specific_energy_j_kg.iter().sum::<f64>(),
        MASS * previous.iter().sum::<f64>(),
    );
    close(radiating.record.stored_energy_change_j, 0.0);
    close(radiating.record.source_j, 0.0);
    close(radiating.record.boundary_outflow_j, 0.0);
    close(radiating.record.radiated_w, 0.0);
    close(radiating.record.closure_j, 0.0);
    assert_eq!(report.exposed_faces, 2);
    assert_eq!(report.patches, 2);
    assert_eq!(previous, [200.0, 100.0, 0.0]);
}

#[test]
fn fixed_temperature_keeps_explicit_phase_energy_and_thermostat_heat() {
    let (domain, solids, mut setup, radiation) = exposed_pair(250.0);
    setup.fixed_temperature.push((0, 300.0));
    let curve = melting();
    let laws: [&dyn TensorCaloricLaw; 2] = [&curve, &curve];
    let previous = [20.0, 0.0];
    let k = conductivity();
    let dt = 0.25;
    let power = radiant_power(300.0, 250.0);
    for fixed_h in [40.0, 80.0] {
        let (step, report) = advance_caloric_tensor_conduction_radiating(
            &domain,
            &solids,
            &setup,
            &[&k, &k],
            TensorCaloricStorage {
                laws: &laws,
                reference_mass_kg: &[MASS; 2],
                previous_specific_energy_j_kg: &previous,
            },
            &[(0, fixed_h), (1, 0.0)],
            0.0,
            dt,
            &controls(),
            &radiation,
            &CancelGate::new(),
        )
        .unwrap();
        let stored = MASS * (fixed_h - previous[0]);
        assert_eq!(step.endpoint.specific_energy_j_kg, vec![fixed_h, 0.0]);
        assert_eq!(step.endpoint.solution.solution.temperature, vec![300.0; 2]);
        close(step.record.stored_energy_change_j, stored);
        close(step.record.source_j, 0.0);
        close(step.record.boundary_outflow_j, -stored);
        close(step.record.radiated_w, power);
        close(step.record.closure_j, 0.0);
        // Thermostat input is storage plus radiation, while the net external
        // outflow already subtracts that same thermostat once.
        let balance = step.endpoint.solution.solution.report.balance;
        close(balance.sink_outflow_w, power);
        close(
            balance.fixed_cell_injection_w * dt + stored,
            power * dt + stored,
        );
        close(report.net_outflow_correction_w, 0.0);
    }
    assert_eq!(previous, [20.0, 0.0]);
}

#[test]
fn exhausted_radiation_and_native_work_or_cancellation_return_no_interval() {
    let (domain, solids, setup, radiation) = exposed_pair(300.0);
    let curve = constant_cp(0.0);
    let laws: [&dyn TensorCaloricLaw; 2] = [&curve, &curve];
    let previous = [200.0, 0.0];
    let saved = setup.clone();
    let k = conductivity();
    let solve =
        |control: &TensorNonlinearConfig, radiation: &RadiationConfig, gate: &CancelGate| {
            advance_caloric_tensor_conduction_radiating(
                &domain,
                &solids,
                &setup,
                &[&k, &k],
                TensorCaloricStorage {
                    laws: &laws,
                    reference_mass_kg: &[MASS; 2],
                    previous_specific_energy_j_kg: &previous,
                },
                &[(1, 0.0)],
                0.0,
                0.25,
                control,
                radiation,
                gate,
            )
        };
    let mut one_radiative = radiation.clone();
    one_radiative.max_iterations = 1;
    assert!(matches!(
        solve(&controls(), &one_radiative, &CancelGate::new()),
        Err(ChtError::SolverNotConverged { .. })
    ));
    let mut one_newton = controls();
    one_newton.max_newton_iterations = 1;
    assert!(solve(&one_newton, &radiation, &CancelGate::new()).is_err());
    let cancelled = CancelGate::new();
    cancelled.request();
    assert_eq!(
        solve(&controls(), &radiation, &cancelled).unwrap_err(),
        ChtError::Cancelled,
    );
    assert_eq!(previous, [200.0, 0.0]);
    assert_eq!(setup, saved);
}
