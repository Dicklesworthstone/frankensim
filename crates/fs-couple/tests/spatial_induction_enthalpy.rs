//! G1/G3/G4: actual 3D A–phi heat and magnetic entropy share finite enthalpy.
//! Synthetic source cards qualify the equations, not an empirical device.
#[path = "../../fs-em/tests/winding_ohmic/common.rs"]
mod source;

use fs_couple::circuit_enthalpy::{
    CircuitEnthalpyError, CircuitEnthalpySolve,
    material::MaterialCell,
    spatial_induction::{
        SpatialMagneticRegion, SpatialMagneticThermalLaw, SpatialMagneticThermalNetwork,
        ThermalMagneticMaterial,
    },
};
use fs_couple::material_circuit::MaterialCircuitSolve;
use fs_em::spatial::{
    SpatialMagneticMaterial,
    vector_potential::{
        SpatialVectorMagnetostatics, VectorMagneticBoundary, VectorMagneticControls,
        VectorPotentialFace,
        induction::{SpatialInduction, SpatialInductionControls, SpatialInductionLimits},
    },
};
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_matdb::{MaterialCard, PropertyKey, SelectionPolicy};
use fs_material::{
    conductor::{ElectricalPropertyBasis, ResolvedConductorState},
    phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase},
    state_point::SPECIFIC_HEAT_CAPACITY_DIMS,
    thermal::source::ThermalMaterialPoint,
};
use fs_qty::Density;
use fs_solver::{CheckedCgConfig, NewtonKrylovConfig};
use fs_thermal::enthalpy::*;

fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.16e} != {b:.16e}; error={:.6e}",
        (a - b).abs()
    );
}
fn assembly() -> WeightedAssemblyLimits {
    WeightedAssemblyLimits {
        max_cells: 100,
        max_dofs: 1000,
        max_triplets: 10000,
    }
}
fn thermal_controls() -> CircuitEnthalpySolve {
    CircuitEnthalpySolve {
        maximum_cells: 8,
        maximum_circuit_trials: 1000,
        maximum_newton_iterations: 24,
        difference_step_j: 1e-3,
        cell_energy_tolerance_j: 1e-8,
        thermal_energy_tolerance_j: 1e-8,
        total_energy_tolerance_j: 2e-8,
        circuit: MaterialCircuitSolve {
            energy_tolerance_j: 1e-9,
            power_tolerance_w: 1e-9,
        },
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-10,
            relative_tolerance: 1e-14,
            forcing_minimum: 1e-4,
            linear_restart: 8,
            ..NewtonKrylovConfig::default()
        },
    }
}
fn field_controls() -> SpatialInductionControls {
    SpatialInductionControls {
        relative_tolerance: 1e-12,
        absolute_residual_a: 1e-12,
        restart: 32,
        maximum_cycles: 100,
        maximum_krylov_scalars: 100000,
        current_tolerance_a: 1e-9,
        flux_tolerance_wb: 1e-10,
        energy_tolerance_j: 1e-9,
        maximum_cell_magnetic_peclet: 100.0,
    }
}
fn curve(card: &MaterialCard) -> ResolvedHeatCapacityCurve {
    ResolvedHeatCapacityCurve::try_new(
        card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
            density: &PropertyKey::new("density", Density::DIMS),
            temperature_axis: "T",
            lower: &source::absolute_point(200.0),
            upper: &source::absolute_point(2000.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-6,
    )
    .unwrap()
}
fn thermal(curve: &ResolvedHeatCapacityCurve, source_w: f64) -> DepositionSimulation<'_> {
    DepositionSimulation::new(
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(7),
                segment_id: 0,
                body_id: 0,
                material_id: MaterialId(9),
                mass_kg: 1.0,
                specific_h_j_kg: 400.0,
                geometry: CellGeometry {
                    start_m: [0.0; 3],
                    end_m: [1.0, 0.0, 0.0],
                    volume_m3: 1.0,
                    surface_area_m2: 1.0,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w,
                },
            }),
        }],
        &[(MaterialId(9), curve.curve() as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 1,
            max_links: 1,
            max_events: 1,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 0.1,
            energy_tolerance_j: 1e-9,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}
fn problem(card: &MaterialCard, basis: ElectricalPropertyBasis) -> SpatialVectorMagnetostatics {
    let p = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [1.0, 0.0, 1.0],
        [0.0, 1.0, 1.0],
        [1.0, 1.0, 1.0],
    ];
    let t = vec![
        [0, 1, 3, 7],
        [0, 3, 2, 7],
        [0, 2, 6, 7],
        [0, 6, 4, 7],
        [0, 4, 5, 7],
        [0, 5, 1, 7],
    ];
    let m = source::bundle_at(card, &source::absolute_point(300.0), basis);
    SpatialVectorMagnetostatics::new(
        p,
        t,
        vec![SpatialMagneticMaterial::Linear(m); 6],
        assembly(),
        &mut || false,
    )
    .unwrap()
}
fn initial(
    problem: &SpatialVectorMagnetostatics,
    basis: ElectricalPropertyBasis,
    amplitude: f64,
    mu: f64,
) -> SpatialInduction<'_> {
    let boundary = VectorMagneticBoundary {
        potential: problem
            .exterior_faces()
            .into_iter()
            .map(|face| VectorPotentialFace {
                face,
                edge_integrals_wb: [0.0; 3],
            })
            .collect(),
    };
    let solve = CheckedCgConfig {
        tolerance: 1e-12,
        max_iterations: 1000,
        max_corrections: 2,
    };
    let initial = problem
        .solve(
            &vec![[8.0 * amplitude / mu; 3]; 6],
            &boundary,
            VectorMagneticControls {
                field: solve,
                gauge: solve,
                maximum_gauge_product_terms: 100000,
                current_tolerance_a: 1e-9,
                flux_tolerance_wb: 1e-10,
                energy_tolerance_j: 1e-9,
            },
            &mut || false,
        )
        .unwrap();
    let sources = problem
        .materials()
        .iter()
        .map(|m| {
            let SpatialMagneticMaterial::Linear(m) = m else {
                unreachable!()
            };
            Some(ResolvedConductorState::try_from_material_with_basis(m.clone(), basis).unwrap())
        })
        .collect();
    let state = SpatialInduction::new(
        &initial,
        sources,
        SpatialInductionLimits {
            maximum_unknowns: 1000,
            maximum_conducting_cells: 100,
            maximum_assembly_terms: 10000,
        },
        1e-10,
        &mut || false,
    )
    .unwrap();
    near(state.potential_wb()[diagonal(problem)], amplitude, 1e-10);
    state
}
fn diagonal(problem: &SpatialVectorMagnetostatics) -> usize {
    problem.mesh().edges.binary_search(&[0, 7]).unwrap()
}
fn regions() -> Vec<SpatialMagneticRegion> {
    (0..6)
        .map(|tetrahedron| SpatialMagneticRegion::Cell {
            tetrahedron,
            cell: CellId(7),
        })
        .collect()
}
fn laws(
    card: &MaterialCard,
    problem: &SpatialVectorMagnetostatics,
    maximum: f64,
) -> Vec<SpatialMagneticThermalLaw> {
    problem
        .materials()
        .iter()
        .enumerate()
        .map(|(tetrahedron, m)| {
            let SpatialMagneticMaterial::Linear(m) = m else {
                unreachable!()
            };
            SpatialMagneticThermalLaw {
                tetrahedron,
                law: ThermalMagneticMaterial::new(card, m, "T", [280.0, 350.0]).unwrap(),
                maximum_flux_density_t: maximum,
            }
        })
        .collect()
}
fn bracket(mut lo: f64, mut hi: f64, residual: impl Fn(f64) -> f64) -> f64 {
    assert!(residual(lo) <= 0.0 && residual(hi) >= 0.0);
    for _ in 0..100 {
        let mid = f64::midpoint(lo, hi);
        if residual(mid) > 0.0 {
            hi = mid
        } else {
            lo = mid
        }
    }
    f64::midpoint(lo, hi)
}
fn sigma(t: f64, rho0: f64, alpha: f64, basis: ElectricalPropertyBasis) -> f64 {
    match basis {
        ElectricalPropertyBasis::Resistivity => 1.0 / (rho0 * (1.0 + alpha * (t - 300.0))),
        ElectricalPropertyBasis::Conductivity => {
            // The authored conductivity card linearly interpolates its actual
            // tabulated conductivity values over [300,2000], not resistivity.
            assert!(t >= 300.0);
            (1.0 / rho0)
                + (t - 300.0) * ((1.0 / (rho0 * (1.0 + 1700.0 * alpha))) - (1.0 / rho0)) / 1700.0
        }
    }
}

#[test]
fn actual_three_dimensional_heat_and_source_temperature_match_independent_scalar_elimination() {
    let dt = 0.1;
    let a0 = 3.0;
    let rho = 0.2;
    let alpha = 0.005;
    let mu0 = 1.0;
    for basis in [
        ElectricalPropertyBasis::Resistivity,
        ElectricalPropertyBasis::Conductivity,
    ] {
        let mut temperatures = Vec::new();
        for slope in [-0.0004, 0.0, 0.0004] {
            let card = source::absolute_affine_magnetic_card(rho, alpha, mu0, slope, basis);
            let curve = curve(&card);
            let thermal = thermal(&curve, 0.0);
            let problem = problem(&card, basis);
            let state = initial(&problem, basis, a0, mu0);
            let context = source::absolute_point(300.0);
            let cells = [MaterialCell {
                cell: CellId(7),
                material_id: MaterialId(9),
                curve: curve.curve(),
                point: ThermalMaterialPoint {
                    material: &card,
                    context: &context,
                    temperature_axis: "T",
                    liquid_fraction_axis: None,
                },
                conductivity: None,
            }];
            let regions = regions();
            let bindings = laws(&card, &problem, 6.1);
            let network = SpatialMagneticThermalNetwork {
                cells: &cells,
                links: &[],
                magnetic_regions: &regions,
                volume_tolerance_m3: 1e-12,
            };
            let candidate = network
                .step_candidate(
                    &thermal,
                    dt,
                    &state,
                    &[[0.0; 3]; 6],
                    &[[0.0; 3]; 8],
                    &bindings,
                    thermal_controls(),
                    field_controls(),
                    assembly(),
                    &[0.0],
                    || false,
                )
                .unwrap();
            // Independent exact Whitney mode: K=4/mu, charge-eliminated mass
            // M=3*sigma/40. Initial J is ramped to zero, so j_mid=K0*a0/2.
            let amplitude = |t: f64| {
                let m = 3.0 * sigma(t, rho, alpha, basis) / 40.0;
                let mu = mu0 + slope * (t - 300.0);
                a0 * m / (m + 2.0 * dt / mu)
            };
            let internal = |t: f64, a: f64| {
                let mu = mu0 + slope * (t - 300.0);
                2.0 * a * a / mu * (1.0 + t * slope / mu)
            };
            let expected = bracket(300.0, 310.0, |t| {
                let a = amplitude(t);
                4.0 * (t - 300.0) + internal(t, a) - internal(300.0, a0) - 2.0 * a0 / mu0 * (a - a0)
            });
            let t = candidate.thermal().temperature_k(CellId(7)).unwrap();
            let a = candidate.field().potential_wb()[diagonal(&problem)];
            near(t, expected, 2e-8);
            near(a, amplitude(expected), 2e-9);
            let report = candidate.report();
            let native = &report.coupled.electrical;
            let m = 3.0 * sigma(t, rho, alpha, basis) / 40.0;
            near(native.heat_j, m * (a - a0).powi(2) / dt, 2e-8);
            near(native.impressed_work_j, 2.0 * a0 / mu0 * (a - a0), 2e-8);
            near(native.mechanical_work_j, 0.0, 1e-12);
            near(report.total_energy_residual_j, 0.0, 2e-8);
            near(report.external_heat_j, 0.0, 2e-8);
            near(report.coupled.deposited_joule_j, native.heat_j, 2e-8);
            near(
                4.0 * (t - 300.0),
                native.heat_j + report.coupled.source.reversible_heat_j,
                2e-8,
            );
            for sample in &report.coupled.source.fields {
                near(sample.temperature_k[1], t, 1e-12);
                near(sample.flux_density_t[1], 2.0 * a.abs(), 2e-9);
            }
            for conductor in candidate.field().conductors().iter().flatten() {
                near(
                    conductor.resistivity_ohm_m(),
                    sigma(t, rho, alpha, basis).recip(),
                    2e-12,
                );
            }
            assert!(std::ptr::eq(
                candidate.field().problem().mesh(),
                state.problem().mesh()
            ));
            assert!(std::ptr::eq(
                candidate.field().problem().positions_m(),
                state.problem().positions_m()
            ));
            assert_eq!(state.time_s(), 0.0);
            near(thermal.temperature_k(CellId(7)).unwrap(), 300.0, 1e-12);
            let continued = network
                .step_candidate(
                    candidate.thermal(),
                    2.0 * dt,
                    candidate.field(),
                    &[[0.0; 3]; 6],
                    &[[0.0; 3]; 8],
                    &bindings,
                    thermal_controls(),
                    field_controls(),
                    assembly(),
                    &[0.0],
                    || false,
                )
                .unwrap();
            near(continued.report().total_energy_residual_j, 0.0, 2e-8);
            assert_eq!(continued.field().time_s(), 2.0 * dt);
            temperatures.push(t);
        }
        assert!(temperatures[0] < temperatures[1] && temperatures[1] < temperatures[2]);
    }
}

#[test]
fn prescribed_motion_heat_and_external_lorentz_work_close_with_real_magnetic_entropy() {
    let basis = ElectricalPropertyBasis::Resistivity;
    let slope = 0.0004;
    let card = source::absolute_affine_magnetic_card(0.2, 0.005, 1.0, slope, basis);
    let curve = curve(&card);
    let thermal = thermal(&curve, 0.0);
    let problem = problem(&card, basis);
    let a0 = 1.0;
    let state = initial(&problem, basis, a0, 1.0);
    let context = source::absolute_point(300.0);
    let cells = [MaterialCell {
        cell: CellId(7),
        material_id: MaterialId(9),
        curve: curve.curve(),
        point: ThermalMaterialPoint {
            material: &card,
            context: &context,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        conductivity: None,
    }];
    let regions = regions();
    let bindings = laws(&card, &problem, 2.1);
    let network = SpatialMagneticThermalNetwork {
        cells: &cells,
        links: &[],
        magnetic_regions: &regions,
        volume_tolerance_m3: 1e-12,
    };
    let dt = 0.05;
    let alpha = 0.6;
    let velocity: Vec<_> = problem
        .positions_m()
        .iter()
        .map(|p| p.map(|x| alpha * x))
        .collect();
    let candidate = network
        .step_candidate(
            &thermal,
            dt,
            &state,
            &[[8.0; 3]; 6],
            &velocity,
            &bindings,
            thermal_controls(),
            field_controls(),
            assembly(),
            &[0.0],
            || false,
        )
        .unwrap();
    let amplitude = |t: f64| {
        let m = 3.0 * sigma(t, 0.2, 0.005, basis) / 40.0;
        let k = 4.0 / (1.0 + slope * (t - 300.0));
        (m / dt + 2.0 - alpha * m) * a0 / (m / dt + 0.5 * k + alpha * m)
    };
    let work = |t: f64| {
        let a = amplitude(t);
        let m = 3.0 * sigma(t, 0.2, 0.005, basis) / 40.0;
        let mid = f64::midpoint(a, a0);
        let rate = (a - a0) / dt + 2.0 * alpha * mid;
        -dt * 2.0 * alpha * mid * m * rate
    };
    let internal = |t: f64, a: f64| {
        let mu = 1.0 + slope * (t - 300.0);
        2.0 * a * a / mu * (1.0 + t * slope / mu)
    };
    let expected = bracket(300.0, 305.0, |t| {
        4.0 * (t - 300.0) + internal(t, amplitude(t)) - internal(300.0, a0) + work(t)
            - 4.0 * a0 * (amplitude(t) - a0)
    });
    let t = candidate.thermal().temperature_k(CellId(7)).unwrap();
    let a = amplitude(t);
    near(t, expected, 2e-8);
    near(
        candidate.field().potential_wb()[diagonal(&problem)],
        a,
        2e-9,
    );
    let native = &candidate.report().coupled.electrical;
    near(native.mechanical_work_j, work(t), 2e-9);
    let m = 3.0 * sigma(t, 0.2, 0.005, basis) / 40.0;
    near(
        native.heat_j,
        dt * m * ((a - a0) / dt + alpha * (a + a0)).powi(2),
        2e-9,
    );
    assert!(native.heat_j > 0.0 && native.mechanical_work_j.abs() > 1e-4);
    near(candidate.report().total_energy_residual_j, 0.0, 2e-8);
}

#[test]
fn joint_spatial_field_and_temperature_refine_to_independent_continuous_equations() {
    let basis = ElectricalPropertyBasis::Resistivity;
    let rho = 0.2;
    let alpha = 0.009;
    let slope = 0.0004;
    let total = 0.1;
    // A substantial resolved thermal transient exposes the endpoint material
    // feedback error at practical refinement resolutions.
    let external = 1400.0;
    let card = source::absolute_affine_magnetic_card(rho, alpha, 1.0, slope, basis);
    let curve = curve(&card);
    let problem = problem(&card, basis);
    let context = source::absolute_point(300.0);
    let cells = [MaterialCell {
        cell: CellId(7),
        material_id: MaterialId(9),
        curve: curve.curve(),
        point: ThermalMaterialPoint {
            material: &card,
            context: &context,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        conductivity: None,
    }];
    let regions = regions();
    let bindings = laws(&card, &problem, 6.1);
    let network = SpatialMagneticThermalNetwork {
        cells: &cells,
        links: &[],
        magnetic_regions: &regions,
        volume_tolerance_m3: 1e-12,
    };
    // Independently derived continuous one-mode induction plus caloric/field
    // entropy equations. No production matrices or thermal residual are used.
    // S_field=2*a^2*mu'/mu^2, C_total=4+T*dS_field/dT.
    let rhs = |time: f64, y: [f64; 2]| {
        let [a, t] = y;
        let mu = 1.0 + slope * (t - 300.0);
        let m = 3.0 * sigma(t, rho, alpha, basis) / 40.0;
        let da = (12.0 * (1.0 - time / total) - 4.0 * a / mu) / m;
        let capacity = 4.0 - 4.0 * t * a * a * slope * slope / (mu * mu * mu);
        let dt = (external + m * da * da - 4.0 * t * a * slope / (mu * mu) * da) / capacity;
        [da, dt]
    };
    let rk4 = |steps: usize| {
        let h = total / steps as f64;
        let mut y = [3.0, 300.0];
        for k in 0..steps {
            let t = k as f64 * h;
            let k1 = rhs(t, y);
            let k2 = rhs(t + 0.5 * h, std::array::from_fn(|i| y[i] + 0.5 * h * k1[i]));
            let k3 = rhs(t + 0.5 * h, std::array::from_fn(|i| y[i] + 0.5 * h * k2[i]));
            let k4 = rhs(t + h, std::array::from_fn(|i| y[i] + h * k3[i]));
            y = std::array::from_fn(|i| {
                y[i] + h * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0
            });
        }
        y
    };
    let exact = rk4(20000);
    let check = rk4(10000);
    for i in 0..2 {
        near(exact[i], check[i], 1e-9);
    }
    let mut errors = Vec::new();
    for steps in [40, 80, 160] {
        let mut thermal = self::thermal(&curve, external);
        let mut state = initial(&problem, basis, 3.0, 1.0);
        for i in 1..=steps {
            let target = i as f64 * total / steps as f64;
            let current = [[24.0 * (1.0 - i as f64 / steps as f64); 3]; 6];
            let result = network
                .step_candidate(
                    &thermal,
                    target,
                    &state,
                    &current,
                    &[[0.0; 3]; 8],
                    &bindings,
                    thermal_controls(),
                    field_controls(),
                    assembly(),
                    &[0.0],
                    || false,
                )
                .unwrap();
            let (next_thermal, next_state, report) = result.into_parts();
            near(report.total_energy_residual_j, 0.0, 2e-8);
            thermal = next_thermal;
            state = next_state;
        }
        let result = [
            state.potential_wb()[diagonal(&problem)],
            thermal.temperature_k(CellId(7)).unwrap(),
        ];
        errors.push([(result[0] - exact[0]).abs(), (result[1] - exact[1]).abs()]);
    }
    for pair in errors.windows(2) {
        for coordinate in 0..2 {
            let order = (pair[0][coordinate] / pair[1][coordinate]).log2();
            assert!(
                (0.8..1.2).contains(&order),
                "joint coordinate {coordinate}: {errors:?}; order={order}"
            );
        }
    }
}

#[test]
fn distinct_tetrahedral_receivers_follow_local_magnetic_entropy_and_change_the_same_field() {
    let basis = ElectricalPropertyBasis::Resistivity;
    let slope = 0.001;
    let card = source::absolute_affine_magnetic_card(0.2, 0.0, 1.0, slope, basis);
    let curve = curve(&card);
    let problem = problem(&card, basis);
    let a0 = 3.0;
    let state = initial(&problem, basis, a0, 1.0);
    let context = source::absolute_point(300.0);
    // Receiver order deliberately differs from field order and CellId order.
    let ids = [
        CellId(80),
        CellId(3),
        CellId(50),
        CellId(7),
        CellId(20),
        CellId(9),
    ];
    let volumes = (0..6)
        .map(|i| problem.cell_volume_m3(i).unwrap())
        .collect::<Vec<_>>();
    let powers = [0.0, 20.0, 40.0, 60.0, 80.0, 100.0];
    let events = (0..6)
        .map(|i| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: ids[i],
                segment_id: i as u64,
                body_id: 0,
                material_id: MaterialId(9),
                mass_kg: volumes[i],
                specific_h_j_kg: 400.0,
                geometry: CellGeometry {
                    start_m: [0.0; 3],
                    end_m: [1.0, 0.0, 0.0],
                    volume_m3: volumes[i],
                    surface_area_m2: 1.0,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: powers[i],
                },
            }),
        })
        .collect();
    let thermal = DepositionSimulation::new(
        events,
        &[(MaterialId(9), curve.curve() as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 6,
            max_links: 1,
            max_events: 6,
            max_substeps: 100,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 0.1,
            energy_tolerance_j: 1e-9,
        },
        || Ok(()),
    )
    .unwrap()
    .0;
    let mut cells = ids
        .iter()
        .map(|&cell| MaterialCell {
            cell,
            material_id: MaterialId(9),
            curve: curve.curve(),
            point: ThermalMaterialPoint {
                material: &card,
                context: &context,
                temperature_axis: "T",
                liquid_fraction_axis: None,
            },
            conductivity: None,
        })
        .collect::<Vec<_>>();
    cells.rotate_left(2);
    let mut regions = ids
        .iter()
        .enumerate()
        .map(|(tetrahedron, &cell)| SpatialMagneticRegion::Cell { tetrahedron, cell })
        .collect::<Vec<_>>();
    regions.reverse();
    let bindings = laws(&card, &problem, 6.1);
    let network = SpatialMagneticThermalNetwork {
        cells: &cells,
        links: &[],
        magnetic_regions: &regions,
        volume_tolerance_m3: 1e-12,
    };
    let dt = 0.1;
    let result = network
        .step_candidate(
            &thermal,
            dt,
            &state,
            &[[0.0; 3]; 6],
            &[[0.0; 3]; 8],
            &bindings,
            thermal_controls(),
            field_controls(),
            assembly(),
            &[0.0; 6],
            || false,
        )
        .unwrap();
    // Uniform sigma leaves the exact charge-eliminated M=3*sigma/40.
    // Each tet has V=1/6 and |B|=2|a|, so K=sum_i 2/(3*mu_i).
    // Eliminate each distinct local caloric equation before the one field root.
    let m = 3.0 * 5.0 / 40.0;
    let temperatures = |a: f64| {
        powers.map(|power| {
            bracket(280.0, 350.0, |t| {
                let mu = 1.0 + slope * (t - 300.0);
                let parameter = a0 * a * (1.0 / mu - 1.0) / 3.0;
                let delta_ts = slope / 3.0 * (t * a * a / (mu * mu) - 300.0 * a0 * a0);
                let reversible = -parameter - delta_ts;
                (4.0 / 6.0) * (t - 300.0)
                    - m * (a - a0).powi(2) / (6.0 * dt)
                    - reversible
                    - dt * power
            })
        })
    };
    let expected_a = bracket(0.0, a0, |a| {
        let reluctivities: f64 = temperatures(a)
            .iter()
            .map(|t| (1.0 + slope * (t - 300.0)).recip())
            .sum();
        a * (m + dt * reluctivities / 3.0) - a0 * m
    });
    let expected_t = temperatures(expected_a);
    near(
        result.field().potential_wb()[diagonal(&problem)],
        expected_a,
        2e-9,
    );
    // Local temperature-dependent mu affects the actual retained field mode.
    assert!((expected_a - a0 * m / (m + 2.0 * dt)).abs() > 1e-3);
    for i in 0..6 {
        let t = result.thermal().temperature_k(ids[i]).unwrap();
        near(t, expected_t[i], 3e-8);
        let receipt = result
            .report()
            .coupled
            .source
            .fields
            .iter()
            .find(|s| s.tetrahedron == i)
            .unwrap();
        assert_eq!(receipt.cell, ids[i]);
        near(receipt.temperature_k[1], expected_t[i], 3e-8);
        near(
            result.report().coupled.electrical.cells[i].heat_j,
            m * (expected_a - a0).powi(2) / (6.0 * dt),
            2e-9,
        );
    }
    assert!(expected_t.windows(2).all(|t| t[0] < t[1]));
    near(
        result.report().external_heat_j,
        dt * powers.iter().sum::<f64>(),
        2e-8,
    );
    near(result.report().total_energy_residual_j, 0.0, 2e-8);
}

#[test]
fn source_geometry_capacity_budgets_and_late_cancellation_preserve_both_histories() {
    let basis = ElectricalPropertyBasis::Resistivity;
    let card = source::absolute_affine_magnetic_card(0.2, 0.005, 1.0, 0.0004, basis);
    let curve = curve(&card);
    let thermal = thermal(&curve, 0.0);
    let problem = problem(&card, basis);
    let state = initial(&problem, basis, 3.0, 1.0);
    let context = source::absolute_point(300.0);
    let cells = [MaterialCell {
        cell: CellId(7),
        material_id: MaterialId(9),
        curve: curve.curve(),
        point: ThermalMaterialPoint {
            material: &card,
            context: &context,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        conductivity: None,
    }];
    let regions = regions();
    let bindings = laws(&card, &problem, 6.1);
    let network = SpatialMagneticThermalNetwork {
        cells: &cells,
        links: &[],
        magnetic_regions: &regions,
        volume_tolerance_m3: 1e-12,
    };
    let run = |laws: &[SpatialMagneticThermalLaw], controls| {
        network.step_candidate(
            &thermal,
            0.1,
            &state,
            &[[0.0; 3]; 6],
            &[[0.0; 3]; 8],
            laws,
            controls,
            field_controls(),
            assembly(),
            &[0.0],
            || false,
        )
    };
    assert!(run(&bindings[..5], thermal_controls()).is_err());
    let mut duplicate = bindings.clone();
    duplicate[5] = duplicate[4].clone();
    assert!(run(&duplicate, thermal_controls()).is_err());
    let mut too_small = bindings.clone();
    too_small[0].maximum_flux_density_t = 5.9;
    assert!(run(&too_small, thermal_controls()).is_err());
    let mut bad_capacity = bindings.clone();
    for law in &mut bad_capacity {
        law.maximum_flux_density_t = 10000.0;
    }
    assert!(run(&bad_capacity, thermal_controls()).is_err());
    let other = source::absolute_affine_magnetic_card(0.2, 0.005, 1.0, 0.0003, basis);
    let other_problem = self::problem(&other, basis);
    assert!(run(&laws(&other, &other_problem, 6.1), thermal_controls()).is_err());
    let restricted = problem
        .materials()
        .iter()
        .enumerate()
        .map(|(tetrahedron, m)| {
            let SpatialMagneticMaterial::Linear(m) = m else {
                unreachable!()
            };
            SpatialMagneticThermalLaw {
                tetrahedron,
                law: ThermalMagneticMaterial::new(&card, m, "T", [300.0, 300.001]).unwrap(),
                maximum_flux_density_t: 6.1,
            }
        })
        .collect::<Vec<_>>();
    // The demanded heating leaves this real source interval; no clamping or
    // source fallback may turn it into an accepted endpoint.
    assert!(run(&restricted, thermal_controls()).is_err());
    let mut budget = thermal_controls();
    budget.maximum_circuit_trials = 1;
    assert!(run(&bindings, budget).is_err());
    let short = SpatialMagneticThermalNetwork {
        cells: &cells,
        links: &[],
        magnetic_regions: &regions[..5],
        volume_tolerance_m3: 1e-12,
    };
    assert!(
        short
            .step_candidate(
                &thermal,
                0.1,
                &state,
                &[[0.0; 3]; 6],
                &[[0.0; 3]; 8],
                &bindings,
                thermal_controls(),
                field_controls(),
                assembly(),
                &[0.0],
                || false
            )
            .is_err()
    );
    let altered_thermal = |volume: f64, boundary_event: bool| {
        let mut body = thermal.cells()[&CellId(7)].clone();
        body.mass_kg = volume;
        body.geometry.volume_m3 = volume;
        let boundary = body.boundary;
        let mut schedule = vec![TimedEvent {
            time_s: 0.0,
            event: Event::Birth(body),
        }];
        if boundary_event {
            schedule.push(TimedEvent {
                time_s: 0.05,
                event: Event::Boundary(
                    CellId(7),
                    Boundary {
                        source_w: 1.0,
                        ..boundary
                    },
                ),
            });
        }
        DepositionSimulation::new(
            schedule,
            &[(MaterialId(9), curve.curve() as &dyn ThermalLaw)],
            DepositionBudget {
                max_cells: 1,
                max_links: 1,
                max_events: 2,
                max_substeps: 10,
            },
        )
        .unwrap()
        .advance_to(
            0.0,
            StepConfig {
                max_dt_s: 0.1,
                energy_tolerance_j: 1e-9,
            },
            || Ok(()),
        )
        .unwrap()
        .0
    };
    for invalid_thermal in [altered_thermal(1.01, false), altered_thermal(1.0, true)] {
        assert!(
            network
                .step_candidate(
                    &invalid_thermal,
                    0.1,
                    &state,
                    &[[0.0; 3]; 6],
                    &[[0.0; 3]; 8],
                    &bindings,
                    thermal_controls(),
                    field_controls(),
                    assembly(),
                    &[0.0],
                    || false
                )
                .is_err()
        );
        assert_eq!(invalid_thermal.time_s(), 0.0);
    }
    let mut calls = 0;
    let accepted = network
        .step_candidate(
            &thermal,
            0.1,
            &state,
            &[[0.0; 3]; 6],
            &[[0.0; 3]; 8],
            &bindings,
            thermal_controls(),
            field_controls(),
            assembly(),
            &[0.0],
            || {
                calls += 1;
                false
            },
        )
        .unwrap();
    let mut count = 0;
    assert!(matches!(
        network.step_candidate(
            &thermal,
            0.1,
            &state,
            &[[0.0; 3]; 6],
            &[[0.0; 3]; 8],
            &bindings,
            thermal_controls(),
            field_controls(),
            assembly(),
            &[0.0],
            || {
                count += 1;
                count == calls - 1
            }
        ),
        Err(CircuitEnthalpyError::Cancelled)
    ));
    assert!(
        network
            .step_candidate(
                &thermal,
                0.2,
                accepted.field(),
                &[[0.0; 3]; 6],
                &[[0.0; 3]; 8],
                &bindings,
                thermal_controls(),
                field_controls(),
                assembly(),
                &[0.0],
                || false
            )
            .is_err()
    );
    let repeated = run(&bindings, thermal_controls()).unwrap();
    assert_eq!(
        accepted.field().potential_wb(),
        repeated.field().potential_wb()
    );
    assert_eq!(accepted.thermal().cells(), repeated.thermal().cells());
    assert_eq!(state.time_s(), 0.0);
    near(thermal.temperature_k(CellId(7)).unwrap(), 300.0, 1e-12);
}
