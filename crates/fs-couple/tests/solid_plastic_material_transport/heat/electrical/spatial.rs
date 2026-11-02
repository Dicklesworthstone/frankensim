//! G1/G3/G4: deformed native P1 current, original finite heat and loaded J2 history.
//! All sources are manufactured. This file's full-boundary material-carried
//! Dirichlet traces supply distributed current; the terminals child instead
//! exercises insulated equipotential electrodes. The heterogeneous loaded case
//! is an independent P1 weak-form benchmark on this mesh, not an exact continuum
//! transmission solution.
#[path = "spatial/dynamics.rs"]
mod dynamics;
#[path = "spatial/terminals.rs"]
mod terminals;
#[path = "spatial/contacts.rs"]
mod contacts;
use super::*;
use fs_couple::solid_enthalpy::plastic::electrical::spatial::*;
use fs_em::{conduction::ConductionControls, spatial::conduction::*};
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_material::conductor::resolve_conductor_state_with_basis;
use fs_solver::CheckedCgConfig;

const DRIVE: f64 = 1.0;
const RHO_SLOPE: f64 = 0.01;

fn spatial_controls() -> SpatialElectricalSolidSolve {
    SpatialElectricalSolidSolve {
        coupled: controls(),
        field: ConductionControls {
            linear: CheckedCgConfig {
                tolerance: 1e-12,
                max_iterations: 200,
                max_corrections: 2,
            },
            current_tolerance_a: 1e-10,
            power_tolerance_w: 1e-10,
        },
        assembly: WeightedAssemblyLimits {
            max_cells: 8,
            max_dofs: 64,
            max_triplets: 288,
        },
        reference_volume_tolerance_m3: 1e-18,
    }
}
fn electrical_boundary(problem: &SpatialConduction) -> SpatialConductionBoundary {
    SpatialConductionBoundary {
        voltage: problem
            .exterior_faces()
            .map(|face| ConductionVoltageFace {
                face,
                voltage_v: problem.mesh().faces[face].map(|v| DRIVE * NODES[v as usize][0] / L),
            })
            .collect(),
        current: Vec::new(),
    }
}
fn with_spatial<R>(
    rho: f64,
    short_domain: bool,
    run: impl FnOnce(
        &Cx<'_>,
        &CancelGate,
        &DepositionSimulation<'_>,
        &ThermoplasticSolid<'_>,
        &SpatialConduction,
        &PreparedHeatNetwork<'_, '_>,
    ) -> R,
) -> R {
    let mut property = scalar("electrical_resistivity", ELECTRICAL_RESISTIVITY_DIMS, rho);
    property.value = PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: Temperature::DIMS,
        knots: vec![(280.0, 0.8 * rho), (320.0, 1.2 * rho)],
        dims: ELECTRICAL_RESISTIVITY_DIMS,
    };
    property.interpolation = InterpolationPolicy::LinearInside;
    property.validity =
        ValidityDomain::unconstrained().with("T", 280.0, if short_domain { 300.05 } else { 320.0 });
    with_heat_fixture(
        fixture_with_property(EMISSIVITY, Some(property)),
        THERMAL_CONTACT,
        |cx, gate, f, thermal, solid, network| {
            let prepared = network.prepare(cx, thermal).unwrap();
            let conductors = network
                .electrical
                .bulk
                .cells
                .iter()
                .enumerate()
                .map(|(i, binding)| {
                    let state = binding
                        .curve
                        .state_at_specific_enthalpy(thermal.cells()[&binding.cell].specific_h_j_kg)
                        .unwrap();
                    let point = binding.point.at(state).unwrap();
                    resolve_conductor_state_with_basis(
                        &f.cards[i],
                        &point,
                        MaterialPropertySelection::SingleClaimOnly,
                        ElectricalPropertyBasis::Resistivity,
                    )
                    .unwrap()
                })
                .collect();
            let problem = SpatialConduction::new(
                NODES.to_vec(),
                TETS.map(|tet| tet.map(|v| v as u32)).to_vec(),
                conductors,
                spatial_controls().assembly,
                &mut || false,
            )
            .unwrap();
            run(cx, gate, thermal, solid, &problem, &prepared)
        },
    )
}

fn affine(f: [[f64; 3]; 3]) -> Vec<(usize, f64)> {
    NODES
        .iter()
        .enumerate()
        .flat_map(|(vertex, x)| {
            (0..3).map(move |axis| {
                let displacement: f64 = (0..3).map(|j| f[axis][j] * x[j]).sum::<f64>() - x[axis];
                (3 * vertex + axis, displacement)
            })
        })
        .collect()
}

fn uniform_temperature(old_t: f64, rho: f64, geometry_factor: f64) -> f64 {
    // BE finite capacity: x*(1+a*(Told-300)+a*x)=dt*P0/C, dt=1.
    let c = (DRIVE / L).powi(2) * geometry_factor / (1000.0 * 5.0 * rho);
    let b = 1.0 + RHO_SLOPE * (old_t - 300.0);
    old_t + 2.0 * c / (b + (b * b + 4.0 * RHO_SLOPE * c).sqrt())
}

fn check_field(
    problem: &SpatialConduction,
    field: &SpatialConductionField,
    temperatures: [f64; 2],
    rho: f64,
    expected_e: [[f64; 3]; 2],
    jacobian: [f64; 2],
) {
    check_field_with_current_tolerance(
        problem,
        field,
        temperatures,
        rho,
        expected_e,
        jacobian,
        2e-9,
    )
}

fn check_field_with_current_tolerance(
    problem: &SpatialConduction,
    field: &SpatialConductionField,
    temperatures: [f64; 2],
    rho: f64,
    expected_e: [[f64; 3]; 2],
    jacobian: [f64; 2],
    current_tolerance_a_m2: f64,
) {
    let mut expected_power = 0.0;
    for i in 0..2 {
        let sigma = 1.0 / (rho * (1.0 + RHO_SLOPE * (temperatures[i] - 300.0)));
        let power = VOLUME * jacobian[i] * sigma * expected_e[i].iter().map(|x| x * x).sum::<f64>();
        for axis in 0..3 {
            close(field.electric_field_v_m[i][axis], expected_e[i][axis], 2e-9);
            close(
                field.current_density_a_m2[i][axis],
                sigma * expected_e[i][axis],
                current_tolerance_a_m2,
            );
        }
        close(
            problem.cell_volume_m3(i).unwrap(),
            VOLUME * jacobian[i],
            1e-18,
        );
        close(field.cell_power_w[i], power, 5e-11);
        close(
            problem.conductors()[i].material().query().axes()["T"],
            temperatures[i],
            1e-10,
        );
        expected_power += power;
    }
    close(field.power_w, expected_power, 1e-10);
    close(field.supplied_power_w, expected_power, 1e-10);
    close(
        field.nodal_boundary_current_a.iter().flatten().sum(),
        0.0,
        1e-10,
    );
    close(field.current_residual_a, 0.0, 1e-10);
    close(field.power_residual_w, 0.0, 1e-10);
}

#[test]
fn affine_shear_and_volume_change_use_actual_field_quadrature_with_reference_mass() {
    let (sx, sy, sz, gamma) = (1.0004, 0.9998, 1.0002, 0.002);
    let deformation = [[sx, gamma, 0.0], [0.0, sy, 0.0], [0.0, 0.0, sz]];
    let prescribed = affine(deformation);
    let j = sx * sy * sz;
    let geometry_factor = j / (sx * sx) * (1.0 + gamma * gamma / (sy * sy));
    let e = [-DRIVE / (L * sx), DRIVE * gamma / (L * sx * sy), 0.0];
    let mut temperatures = Vec::new();
    for rho in [5.0, 20.0] {
        with_spatial(rho, false, |cx, _, thermal, solid, problem, prepared| {
            let boundary = electrical_boundary(problem);
            let saved_cells = thermal.cells().clone();
            let saved_solid = solid.solution().clone();
            let next = solid
                .advance_spatial_electrical_candidate(
                    thermal,
                    problem,
                    prepared,
                    &boundary,
                    1.0,
                    &prescribed,
                    &[],
                    IncrementSettings::default(),
                    spatial_controls(),
                    &[0.0; 2],
                    cx,
                )
                .unwrap();
            let t = uniform_temperature(300.0, rho, geometry_factor);
            for id in CELLS {
                close(next.thermal().temperature_k(id).unwrap(), t, 3e-7);
                assert_eq!(
                    next.thermal().cells()[&id].mass_kg,
                    thermal.cells()[&id].mass_kg
                );
                assert_eq!(
                    next.thermal().cells()[&id].geometry,
                    thermal.cells()[&id].geometry
                );
            }
            // P=integral sigma |F^-T grad_X phi|² J dV0, not a wire L/A fit.
            let report = next.report();
            check_field(
                next.problem(),
                &report.circuit_thermal.electrical.field,
                CELLS.map(|id| next.thermal().temperature_k(id).unwrap()),
                rho,
                [e; 2],
                [j; 2],
            );
            close(report.physical_plastic_heat_j, 0.0, 1e-13);
            close(report.reversible_heat_j, 0.0, 1e-13);
            assert!(report.numerical_loss_j > 0.0);
            close(
                report.circuit_thermal.thermal.stored_energy_change_j,
                report.circuit_thermal.deposited_joule_j,
                1e-9,
            );
            close(report.combined_residual_j, 0.0, 1e-8);
            for (node, x) in NODES.iter().enumerate() {
                for axis in 0..3 {
                    let expected: f64 = (0..3).map(|k| deformation[axis][k] * x[k]).sum();
                    close(next.problem().positions_m()[node][axis], expected, 1e-14);
                }
            }
            let continued = next
                .solid()
                .advance_spatial_electrical_candidate(
                    next.thermal(),
                    next.problem(),
                    prepared,
                    &boundary,
                    2.0,
                    &prescribed,
                    &[],
                    IncrementSettings::default(),
                    spatial_controls(),
                    &[0.0; 2],
                    cx,
                )
                .unwrap();
            let t2 = uniform_temperature(t, rho, geometry_factor);
            check_field(
                continued.problem(),
                &continued.report().circuit_thermal.electrical.field,
                CELLS.map(|id| continued.thermal().temperature_k(id).unwrap()),
                rho,
                [e; 2],
                [j; 2],
            );
            for id in CELLS {
                close(continued.thermal().temperature_k(id).unwrap(), t2, 4e-7);
            }
            assert_eq!(
                continued.problem().positions_m(),
                next.problem().positions_m()
            );
            for i in 0..2 {
                assert_eq!(
                    continued.problem().conductors()[i]
                        .source_property()
                        .answer()
                        .receipt
                        .selected,
                    problem.conductors()[i]
                        .source_property()
                        .answer()
                        .receipt
                        .selected
                );
            }
            assert_eq!(thermal.cells(), &saved_cells);
            assert_eq!(solid.solution(), &saved_solid);
            assert_eq!(problem.positions_m(), &NODES);
            temperatures.push(t);
        });
    }
    assert!(temperatures[0] > temperatures[1] + 0.2);
}

// Independently eliminate the native J2 free-apex equation, P1 field gradient,
// source rho(T), and conservative two-cell contact/radiation exchange.
// Both tetrahedra have all voltage nodes prescribed, so cell E is exactly the
// affine pullback for this discrete mesh even when the two temperatures differ.
fn loaded_reference(rho: f64) -> [f64; 4] {
    let capacity = MASS * 5.0;
    let bulk = LINK_AREA / (L / (2.0 * 0.01));
    let conductance = 1.0 / (1.0 / bulk + THERMAL_CONTACT / LINK_AREA);
    let radiation = 5.670374419e-8 * L * L / (2.0 / EMISSIVITY - 1.0);
    let exchange =
        |ta: f64, tb: f64| conductance * (ta - tb) + radiation * (ta.powi(4) - tb.powi(4));
    let power = |temperature: f64, extra_shear: f64| {
        VOLUME * (DRIVE / L).powi(2) * (1.0 + 4.0 * SHEAR * SHEAR + extra_shear * extra_shear)
            / (rho * (1.0 + RHO_SLOPE * (temperature - 300.0)))
    };
    let lower = |ta: f64| {
        let (mut lo, mut hi) = (300.0, 320.0);
        for _ in 0..80 {
            let tb = 0.5 * (lo + hi);
            if capacity * (tb - 300.0) - power(tb, 0.0) - exchange(ta, tb) < 0.0 {
                lo = tb;
            } else {
                hi = tb;
            }
        }
        0.5 * (lo + hi)
    };
    let residual = |ta: f64| {
        let (strain, alpha, fy) = apex_oracle(ta);
        capacity * (ta - 300.0) - VOLUME * fy * alpha - power(ta, 2.0 * strain)
            + exchange(ta, lower(ta))
    };
    let (mut lo, mut hi) = (300.0, 319.0);
    assert!(residual(lo) < 0.0 && residual(hi) > 0.0);
    for _ in 0..80 {
        let ta = 0.5 * (lo + hi);
        if residual(ta) < 0.0 {
            lo = ta;
        } else {
            hi = ta;
        }
    }
    let ta = 0.5 * (lo + hi);
    let (strain, alpha, _) = apex_oracle(ta);
    [ta, lower(ta), 2.0 * L * strain, alpha]
}

#[test]
fn spatial_current_heating_changes_actual_free_apex_through_original_yield_data() {
    let mut endpoints = Vec::new();
    for rho in [1.0, 4.0] {
        with_spatial(rho, false, |cx, _, thermal, solid, problem, prepared| {
            let expected = loaded_reference(rho);
            let boundary = electrical_boundary(problem);
            let next = solid
                .advance_spatial_electrical_candidate(
                    thermal,
                    problem,
                    prepared,
                    &boundary,
                    1.0,
                    &supports(SHEAR),
                    &applied(),
                    IncrementSettings::default(),
                    spatial_controls(),
                    &[0.0; 2],
                    cx,
                )
                .unwrap();
            let actual = [
                next.thermal().temperature_k(CELLS[0]).unwrap(),
                next.thermal().temperature_k(CELLS[1]).unwrap(),
                next.solid().solution().displacement_m[3][0],
                next.solid().solution().element_states[0].alpha,
            ];
            for i in 0..2 {
                close(actual[i], expected[i], 3e-6);
            }
            close(actual[2], expected[2], 3e-11);
            close(actual[3], expected[3], 2e-9);
            let scale = DRIVE / L;
            let e = [
                [-scale, 2.0 * SHEAR * scale, actual[2] / L * scale],
                [-scale, 2.0 * SHEAR * scale, 0.0],
            ];
            check_field(
                next.problem(),
                &next.report().circuit_thermal.electrical.field,
                [actual[0], actual[1]],
                rho,
                e,
                [1.0; 2],
            );
            let report = next.report();
            assert!(report.physical_plastic_heat_j > 0.0);
            assert!(report.numerical_loss_j > 0.0);
            close(
                report.circuit_thermal.thermal.stored_energy_change_j,
                report.deposited_mechanical_heat_j + report.circuit_thermal.deposited_joule_j,
                1e-8,
            );
            close(report.combined_residual_j, 0.0, 1e-8);
            close(report.external_heat_j, 0.0, 1e-11);
            assert_eq!(next.solid().time_s(), next.thermal().time_s());
            endpoints.push(actual);
        });
    }
    assert!(endpoints[0][0] > endpoints[1][0] + 0.5);
    assert!(endpoints[0][2] > endpoints[1][2] + 1e-8);
    assert!(endpoints[0][3] > endpoints[1][3] + 1e-7);
}

#[test]
fn late_spatial_source_or_budget_failure_preserves_heat_configuration_and_material_history() {
    for short_domain in [false, true] {
        with_spatial(
            1.0,
            short_domain,
            |cx, _, thermal, solid, problem, prepared| {
                let original_cells = thermal.cells().clone();
                let original_solid = solid.solution().clone();
                let original_materials = solid.materials().to_vec();
                let original_conductors = problem.conductors().to_vec();
                let boundary = electrical_boundary(problem);
                if !short_domain {
                    let good = solid
                        .advance_spatial_electrical_candidate(
                            thermal,
                            problem,
                            prepared,
                            &boundary,
                            1.0,
                            &supports(SHEAR),
                            &applied(),
                            IncrementSettings::default(),
                            spatial_controls(),
                            &[0.0; 2],
                            cx,
                        )
                        .unwrap();
                    let mut controls = spatial_controls();
                    controls.coupled.maximum_mechanical_trials =
                        good.report().mechanical_trials - 1;
                    assert!(matches!(
                        solid.advance_spatial_electrical_candidate(
                            thermal,
                            problem,
                            prepared,
                            &boundary,
                            1.0,
                            &supports(SHEAR),
                            &applied(),
                            IncrementSettings::default(),
                            controls,
                            &[0.0; 2],
                            cx
                        ),
                        Err(ElectricalSolidError::Coupled(CircuitEnthalpyError::Source(
                            ElectricalSolidSourceError::MechanicalBudget
                        )))
                    ));
                    // A warmed, moved solid must not reuse its old current geometry.
                    assert!(
                        good.solid()
                            .advance_spatial_electrical_candidate(
                                good.thermal(),
                                problem,
                                prepared,
                                &boundary,
                                2.0,
                                &supports(SHEAR),
                                &applied(),
                                IncrementSettings::default(),
                                spatial_controls(),
                                &[0.0; 2],
                                cx
                            )
                            .is_err()
                    );
                } else {
                    assert!(
                        solid
                            .advance_spatial_electrical_candidate(
                                thermal,
                                problem,
                                prepared,
                                &boundary,
                                1.0,
                                &supports(SHEAR),
                                &applied(),
                                IncrementSettings::default(),
                                spatial_controls(),
                                &[0.0; 2],
                                cx
                            )
                            .is_err()
                    );
                }
                assert_eq!(thermal.cells(), &original_cells);
                assert_eq!(solid.solution(), &original_solid);
                assert_eq!(solid.materials(), &original_materials);
                assert_eq!(problem.conductors(), &original_conductors);
                assert_eq!(problem.positions_m(), &NODES);
            },
        );
    }
}
