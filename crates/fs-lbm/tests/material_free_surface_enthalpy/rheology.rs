//! T/shear data changes real native momentum, molecular heat and thermal links.

use super::*;
use fs_lbm::d3q19::{E3, Q3, equilibrium3};
use fs_lbm::material::free_surface::enthalpy::MaterialFreeSurfaceEnthalpy3;
use fs_lbm::material::free_surface::enthalpy::rheology::RheologySolve;
use fs_material::liquid::surface::ResolvedLiquidSurfaceState;
use fs_qty::Frequency;

const CELL: [usize; 3] = [3; 3];

fn solve() -> RheologySolve {
    RheologySolve {
        maximum_iterations: 128,
        relaxation: 0.5,
        absolute_relaxation_tolerance: 1e-13,
        relative_relaxation_tolerance: 1e-13,
    }
}

fn query_point() -> QueryPoint {
    point(400.0)
        .with_quantity("rate", QuantitySpec::dimensional(Frequency::DIMS), 0.0)
        .unwrap()
}

fn rheological_card(
    temperature_slope: f64,
    shear_slope: f64,
    conductivity_slope: f64,
    max_rate: f64,
) -> MaterialCard {
    let original = card(0.02, 1000.0);
    let mut claims = ClaimSet::new();
    for (_, original) in original.claims().claims_ordered() {
        let mut claim = original.clone();
        let viscosity = claim.key.name() == "dynamic_viscosity";
        if viscosity || claim.key.name() == "thermal_conductivity" {
            let dims = claim.key.dims();
            claim.value = PropertyValue::Surface2D {
                axes: ["T".into(), "rate".into()],
                axis_dims: [Temperature::DIMS, Frequency::DIMS],
                samples: [
                    [300.0, 0.0],
                    [500.0, 0.0],
                    [300.0, max_rate],
                    [500.0, max_rate],
                ]
                .into_iter()
                .map(|p| {
                    (
                        p,
                        if viscosity {
                            0.1 + temperature_slope * (p[0] - 400.0) + shear_slope * p[1]
                        } else {
                            100.0 + 0.5 * (p[0] - 300.0) + conductivity_slope * p[1]
                        },
                    )
                })
                .collect(),
                triangles: vec![[0, 1, 2], [1, 3, 2]],
                dims,
            };
            claim.validity = claim.validity.with_quantity(
                "rate",
                QuantitySpec::dimensional(Frequency::DIMS),
                0.0,
                max_rate,
            );
        }
        claims.insert_claim(claim).unwrap();
    }
    MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap()
}

#[allow(clippy::too_many_arguments)]
fn with_liquid<R>(
    card: &MaterialCard,
    evolving: bool,
    full: bool,
    gravity: [f64; 3],
    solve: RheologySolve,
    initial_h: impl Fn([usize; 3]) -> f64,
    exercise: impl FnOnce(&mut MaterialFreeSurfaceEnthalpy3<'_, '_>, &CancelGate) -> R,
) -> R {
    let law = LiquidTemperatureLaw::try_new(
        card,
        &query_point(),
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let source = if evolving {
        ResolvedLiquidSurfaceState::try_for_shear_rate_evolution(
            law.reference_state().material().clone(),
            LiquidViscosityBasis::Dynamic,
            "T",
            "rate",
        )
    } else {
        ResolvedLiquidSurfaceState::try_from_material(
            law.reference_state().material().clone(),
            LiquidViscosityBasis::Dynamic,
            LiquidCapillarity::Neglected,
        )
    }
    .unwrap();
    let chart = caloric(card);
    let gate = CancelGate::new_clock_free();
    let mut cells = initial(full);
    for cell in &mut cells {
        cell.velocity_m_s = [0.003 * (cell.coordinates[2] as f64 - 3.5), 0.0, 0.0];
    }
    with_cx(&gate, |cx| {
        let scaling = plan_free_surface_from_material(&source, 0.001, 0.001, 0.1).unwrap();
        let domain = MaterialFreeSurfaceDomain3 {
            gravity_m_s2: gravity,
            ..domain()
        };
        let mut liquid = if evolving {
            // An admitted rate sample cannot leak into a frozen solver.
            assert!(scaling.free_surface3(domain, |_, _, _| true).is_err());
            assert!(
                scaling
                    .free_surface3_from_cells(domain, &cells, cx)
                    .is_err()
            );
            scaling.free_surface3_from_cells_with_rheology(
                domain,
                &cells,
                &law,
                &chart,
                FreeSurfaceEnthalpyControls::default(),
                solve,
                initial_h,
                cx,
                &gate,
            )
        } else {
            scaling
                .free_surface3_from_cells(domain, &cells, cx)
                .unwrap()
                .with_enthalpy(
                    &law,
                    &chart,
                    FreeSurfaceEnthalpyControls::default(),
                    initial_h,
                    &gate,
                )
        }
        .unwrap();
        exercise(&mut liquid, &gate)
    })
}

fn coordinate(state: &fs_material::liquid::ResolvedLiquidState, axis: &str) -> f64 {
    state
        .material()
        .query_point()
        .iter()
        .find(|(name, _)| name == axis)
        .unwrap()
        .1
}

// Independently reduce the physical nonequilibrium tensor from diagnostic
// populations. The full 8^3 fixture has all eight Morton-ordered tiles active.
fn stress_numerator(
    liquid: &MaterialFreeSurfaceEnthalpy3<'_, '_>,
    cell: [usize; 3],
    gravity: [f64; 3],
) -> f64 {
    let lattice = liquid.flow().lattice();
    let bits = lattice.grid().state_bits();
    let slot = (cell[0] / 4) | ((cell[1] / 4) << 1) | ((cell[2] / 4) << 2);
    let lane = cell[0] % 4 + 4 * (cell[1] % 4 + 4 * (cell[2] % 4));
    let state = lattice.sample_at(cell[0], cell[1], cell[2]).unwrap();
    let equilibrium = equilibrium3(state.density, state.velocity);
    let force = gravity.map(|g| g * 0.001);
    let mut tensor = [[0.0; 3]; 3];
    for q in 0..Q3 {
        let e = [f64::from(E3[q].0), f64::from(E3[q].1), f64::from(E3[q].2)];
        let neq = f64::from_bits(bits[(slot * Q3 + q) * 64 + lane]) - equilibrium[q];
        for a in 0..3 {
            for b in 0..3 {
                tensor[a][b] += e[a] * e[b] * neq;
            }
        }
    }
    for a in 0..3 {
        for b in 0..3 {
            tensor[a][b] += 0.5 * (state.velocity[a] * force[b] + state.velocity[b] * force[a]);
        }
    }
    let squared = tensor.iter().flatten().map(|x| x * x).sum::<f64>();
    3.0 * (0.5 * squared).sqrt() / state.density
}

fn independent_rate(
    numerator: f64,
    temperature: f64,
    temperature_slope: f64,
    shear_slope: f64,
) -> f64 {
    // gamma*(a + 3*b*gamma) = numerator/dt, because tau=0.5+3*mu here.
    let a = 0.5 + 3.0 * (0.1 + temperature_slope * (temperature - 400.0));
    let rhs = numerator / 0.001;
    2.0 * rhs / (a + (a * a + 12.0 * shear_slope * rhs).sqrt())
}

#[test]
fn g1_native_stress_closes_the_constitutive_root_and_same_molecular_heat() {
    let source = rheological_card(0.0002, 0.003, 0.0, 1000.0);
    for gravity in [[0.0; 3], [0.03, -0.02, 0.01]] {
        with_liquid(
            &source,
            true,
            true,
            gravity,
            solve(),
            |_| 200000.0,
            |liquid, gate| {
                liquid.step(gate).unwrap();
                let material = liquid.material_at(CELL).unwrap();
                let sample = *liquid.rheology_at(CELL).unwrap();
                let temperature = coordinate(material, "T");
                let expected_rate = independent_rate(
                    stress_numerator(liquid, CELL, gravity),
                    temperature,
                    0.0002,
                    0.003,
                );
                near(sample.shear_rate_per_s, expected_rate, 2e-10);
                near(coordinate(material, "rate"), expected_rate, 2e-10);
                near(
                    sample.collision_tau,
                    0.5 + 3.0 * material.dynamic_viscosity_pa_s(),
                    2e-14,
                );
                let mu = material.dynamic_viscosity_pa_s();
                let energy = liquid.total_enthalpy_j().unwrap();
                let report = liquid.step_with_viscous_heating(&[], gate).unwrap();
                let applied = report
                    .applied
                    .iter()
                    .find(|s| s.coordinates == CELL)
                    .unwrap();
                near(applied.material.dynamic_viscosity_pa_s(), mu, 2e-12);
                near(applied.shear_rate_per_s, expected_rate, 2e-10);
                near(
                    applied.heat_j,
                    mu * expected_rate.powi(2) * applied.liquid_volume_m3 * 0.001,
                    4e-10,
                );
                assert!(report.flow.viscous_heat_j > 0.0);
                assert_eq!(report.flow.external_heat_j, 0.0);
                assert!(
                    (liquid.total_enthalpy_j().unwrap() - energy - report.flow.viscous_heat_j)
                        .abs()
                        < 1e-8
                );
                // Endpoint rate must use the NEW native populations, not the old
                // flow hidden behind the still-unpublished candidate.
                let material = liquid.material_at(CELL).unwrap();
                near(
                    liquid.rheology_at(CELL).unwrap().shear_rate_per_s,
                    independent_rate(
                        stress_numerator(liquid, CELL, gravity),
                        coordinate(material, "T"),
                        0.0002,
                        0.003,
                    ),
                    2e-10,
                );
            },
        );
    }
}

#[test]
fn g3_data_only_temperature_and_shear_swap_changes_pour_and_finite_heat() {
    let mut results = Vec::new();
    for (temperature_slope, shear_slope) in
        [(0.0, 0.0), (0.0, 0.003), (0.0003, 0.003), (0.0, -0.00005)]
    {
        let source = rheological_card(temperature_slope, shear_slope, 0.0, 1000.0);
        results.push(with_liquid(
            &source,
            true,
            false,
            [0.0; 3],
            solve(),
            |_| 100000.0,
            |liquid, gate| {
                liquid.step(gate).unwrap();
                let mass = liquid.flow().ledger_mass_kg().unwrap();
                let report = liquid
                    .step_with_viscous_heating(
                        &[FreeSurfaceCellHeat {
                            coordinates: CELL,
                            energy_j: 0.01,
                        }],
                        gate,
                    )
                    .unwrap();
                near(report.flow.mass_after_kg, mass, 1e-12);
                assert!(
                    (report.flow.enthalpy_after_j
                        - report.flow.enthalpy_before_j
                        - 0.01
                        - report.flow.viscous_heat_j)
                        .abs()
                        < 1e-8
                );
                (
                    liquid.flow().lattice().grid().state_bits(),
                    liquid.specific_enthalpy_at(CELL).unwrap(),
                    report.flow.viscous_heat_j,
                )
            },
        ));
    }
    for pair in results.windows(2) {
        assert_ne!(pair[0].0, pair[1].0);
        assert_ne!(pair[0].1.to_bits(), pair[1].1.to_bits());
        assert_ne!(pair[0].2.to_bits(), pair[1].2.to_bits());
    }
}

#[test]
fn g5_temperature_only_newtonian_data_preserves_existing_native_and_heat_bits() {
    let source = card(0.02, 1000.0);
    let run = |evolving| {
        with_liquid(
            &source,
            evolving,
            true,
            [0.0; 3],
            solve(),
            |_| 200000.0,
            |liquid, gate| {
                liquid.step(gate).unwrap();
                let report = liquid
                    .step_with_viscous_heating(
                        &[FreeSurfaceCellHeat {
                            coordinates: CELL,
                            energy_j: 0.01,
                        }],
                        gate,
                    )
                    .unwrap();
                (
                    liquid.flow().lattice().grid().state_bits(),
                    liquid.total_enthalpy_j().unwrap().to_bits(),
                    report.flow.viscous_heat_j.to_bits(),
                )
            },
        )
    };
    assert_eq!(run(false), run(true));
}

#[test]
fn g4_endpoint_rate_refusal_and_nonconvergence_preserve_entire_liquid() {
    for (maximum_rate, iterations) in [(0.05, 128), (1000.0, 1)] {
        let source = rheological_card(0.0, 0.003, 0.0, maximum_rate);
        with_liquid(
            &source,
            true,
            true,
            [0.0; 3],
            RheologySolve {
                maximum_iterations: iterations,
                ..solve()
            },
            |_| 200000.0,
            |liquid, gate| {
                let bits = liquid.flow().lattice().grid().state_bits();
                let h = liquid.total_enthalpy_j().unwrap();
                let material = liquid.material_at(CELL).unwrap().clone();
                let rheology = *liquid.rheology_at(CELL).unwrap();
                let failure = liquid
                    .step_with_heat(
                        &[FreeSurfaceCellHeat {
                            coordinates: CELL,
                            energy_j: 0.001,
                        }],
                        gate,
                    )
                    .unwrap_err();
                if iterations == 1 {
                    assert!(matches!(
                        failure,
                        FreeSurfaceEnthalpyError::RheologyConvergence { .. }
                    ));
                } else {
                    assert!(matches!(failure, FreeSurfaceEnthalpyError::Material(_)));
                }
                assert_eq!(liquid.flow().lattice().grid().state_bits(), bits);
                assert_eq!(liquid.total_enthalpy_j().unwrap(), h);
                assert_eq!(liquid.material_at(CELL), Some(&material));
                assert_eq!(liquid.rheology_at(CELL), Some(&rheology));
            },
        );
    }
    let source = rheological_card(0.0, 0.003, 0.0, 1000.0);
    with_liquid(
        &source,
        true,
        true,
        [0.0; 3],
        solve(),
        |_| 200000.0,
        |liquid, gate| {
            let bits = liquid.flow().lattice().grid().state_bits();
            let h = liquid.total_enthalpy_j().unwrap();
            let cancelled = CancelGate::new_clock_free();
            cancelled.request();
            assert!(liquid.step_with_viscous_heating(&[], &cancelled).is_err());
            assert_eq!(liquid.flow().lattice().grid().state_bits(), bits);
            assert_eq!(liquid.total_enthalpy_j().unwrap(), h);
            liquid.step_with_viscous_heating(&[], gate).unwrap();
        },
    );
}

#[cfg(feature = "enthalpy-coupling")]
#[test]
fn g1_thermal_substeps_share_temperature_shear_state_and_energy_account() {
    use fs_lbm::material::free_surface::enthalpy::conduction::{
        FreeSurfaceConductionControls, FreeSurfaceConductionFace,
    };
    let pair = [CELL, [4, 3, 3]];
    let mut endpoints = Vec::new();
    for conductivity_slope in [0.0, 10.0] {
        let source = rheological_card(0.0002, 0.003, conductivity_slope, 1000.0);
        endpoints.push(with_liquid(
            &source,
            true,
            true,
            [0.0; 3],
            solve(),
            |p| {
                if p == pair[0] {
                    100000.0
                } else if p == pair[1] {
                    300000.0
                } else {
                    200000.0
                }
            },
            |liquid, gate| {
                liquid.step(gate).unwrap();
                let numerator = pair.map(|p| stress_numerator(liquid, p, [0.0; 3]));
                let report = liquid
                    .step_with_conduction_and_viscous_heating(
                        &[FreeSurfaceConductionFace {
                            cells: pair,
                            area_m2: 1e-6,
                            distances_m: [0.0005; 2],
                        }],
                        &[],
                        FreeSurfaceConductionControls {
                            maximum_faces: 1,
                            maximum_substeps: 100,
                            maximum_dt_s: 0.00025,
                            energy_tolerance_j: 1e-8,
                        },
                        gate,
                    )
                    .unwrap();
                for samples in [
                    &report.conduction.last_applied,
                    &report.conduction.thermal_endpoint,
                ] {
                    let sample = &samples[0];
                    let mut k = [0.0; 2];
                    for i in 0..2 {
                        let material = &sample.materials[i];
                        let temperature = coordinate(material, "T");
                        let rate = coordinate(material, "rate");
                        near(
                            rate,
                            independent_rate(numerator[i], temperature, 0.0002, 0.003),
                            2e-10,
                        );
                        k[i] = 100.0 + 0.5 * (temperature - 300.0) + conductivity_slope * rate;
                        near(material.thermal_conductivity_w_m_k(), k[i], 2e-13);
                    }
                    near(
                        sample.conductance_w_k,
                        1e-6 / (0.0005 / k[0] + 0.0005 / k[1]),
                        2e-13,
                    );
                }
                assert!(report.conduction.thermal.substeps >= 4);
                assert_eq!(report.conduction.flow.external_heat_j, 0.0);
                near(
                    report.conduction.flow.viscous_heat_j,
                    report.applied.iter().map(|s| s.heat_j).sum(),
                    1e-13,
                );
                assert!(report.conduction.flow.balance_defect_j.abs() < 1e-8);
                liquid.specific_enthalpy_at(CELL).unwrap()
            },
        ));
    }
    assert_ne!(endpoints[0].to_bits(), endpoints[1].to_bits());
}
