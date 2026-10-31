//! Complete finite coolant closes the original three-dimensional A-phi field.
//! Manufactured source data check equations; they are not measured device data.
use super::*;
#[path = "../../../../../../fs-em/tests/winding_ohmic/common.rs"]
mod source;

use fs_couple::circuit_enthalpy::magnetic::caloric::MagneticConstitutiveThermalLaw;
use fs_couple::circuit_enthalpy::spatial_induction::{
    SpatialInductionHeatError, SpatialMagneticRegion, SpatialMagneticThermalLaw,
    ThermalMagneticMaterial, heat_assembly_candidate, nonlinear_heat_assembly_candidate,
    nonlinear_thermal_heat_assembly_candidate,
};
use fs_em::spatial::{
    SpatialMagneticMaterial,
    nonlinear::SpatialNonlinearMaterial,
    vector_potential::{
        SpatialVectorMagnetostatics, VectorMagneticBoundary, VectorMagneticControls,
        VectorPotentialFace,
        induction::{
            SpatialInduction, SpatialInductionControls, SpatialInductionLimits,
            nonlinear::{SpatialNonlinearInduction, SpatialNonlinearInductionControls},
        },
        nonlinear::{NonlinearVectorControls, SpatialNonlinearVectorMagnetostatics},
    },
};
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_material::conductor::ResolvedConductorState;
use fs_solver::CheckedCgConfig;

const AMPLITUDE: f64 = 3.0;
const SOLID_CAPACITY: f64 = 4.0;

fn assembly_limits() -> WeightedAssemblyLimits {
    WeightedAssemblyLimits {
        max_cells: 100,
        max_dofs: 1000,
        max_triplets: 10000,
    }
}
fn geometry() -> (Vec<[f64; 3]>, Vec<[u32; 4]>) {
    (
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [0.0, 1.0, 1.0],
            [1.0, 1.0, 1.0],
        ],
        vec![
            [0, 1, 3, 7],
            [0, 3, 2, 7],
            [0, 2, 6, 7],
            [0, 6, 4, 7],
            [0, 4, 5, 7],
            [0, 5, 1, 7],
        ],
    )
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
fn vector_controls() -> VectorMagneticControls {
    let solve = CheckedCgConfig {
        tolerance: 1e-12,
        max_iterations: 1000,
        max_corrections: 2,
    };
    VectorMagneticControls {
        field: solve,
        gauge: solve,
        maximum_gauge_product_terms: 100000,
        current_tolerance_a: 1e-9,
        flux_tolerance_wb: 1e-10,
        energy_tolerance_j: 1e-9,
    }
}
fn induction_limits() -> SpatialInductionLimits {
    SpatialInductionLimits {
        maximum_unknowns: 1000,
        maximum_conducting_cells: 100,
        maximum_assembly_terms: 10000,
    }
}
fn solve_controls() -> CircuitEnthalpySolve {
    let mut solve = controls();
    solve.maximum_cells = 2;
    solve.maximum_circuit_trials = 1000;
    solve.newton.linear_restart = 2;
    solve
}
fn solid_curve(card: &MaterialCard) -> ResolvedHeatCapacityCurve {
    ResolvedHeatCapacityCurve::try_new(
        card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new(
                SPECIFIC_HEAT_CAPACITY_PROPERTY,
                SPECIFIC_HEAT_CAPACITY_DIMS,
            ),
            density: &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
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
fn initial_thermal<'a>(
    solid: &'a ResolvedHeatCapacityCurve,
    fluid: &'a ResolvedHeatCapacityCurve,
    cp: f64,
) -> DepositionSimulation<'a> {
    let mut events = [
        (BODY, 0, 1.0, 1.0, SOLID_CAPACITY * 100.0),
        (FLUID, 1, MASS, 1e-6, cp * 10.0),
    ]
    .map(
        |(id, material, mass_kg, volume_m3, specific_h_j_kg)| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id,
                segment_id: material,
                body_id: material,
                material_id: MaterialId(material),
                mass_kg,
                specific_h_j_kg,
                geometry: CellGeometry {
                    start_m: [0.0, material as f64, 0.0],
                    end_m: [1.0, material as f64, 0.0],
                    volume_m3,
                    surface_area_m2: SURFACE,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        },
    )
    .to_vec();
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: LINK,
            a: BODY,
            b: FLUID,
            area_m2: SURFACE,
            conductance_w_k: 0.0,
        }),
    });
    DepositionSimulation::new(
        events,
        &[
            (MaterialId(0), solid.curve() as &dyn ThermalLaw),
            (MaterialId(1), fluid.curve() as &dyn ThermalLaw),
        ],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
            max_substeps: 100,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: DT,
            energy_tolerance_j: 1e-9,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}
fn initial_linear(problem: &SpatialVectorMagnetostatics) -> SpatialInduction<'_> {
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
    let solved = problem
        .solve(
            &[[8.0 * AMPLITUDE; 3]; 6],
            &boundary,
            vector_controls(),
            &mut || false,
        )
        .unwrap();
    let conductors = problem
        .materials()
        .iter()
        .map(|m| {
            let SpatialMagneticMaterial::Linear(m) = m else {
                unreachable!()
            };
            Some(
                ResolvedConductorState::try_from_material_with_basis(
                    m.clone(),
                    ElectricalPropertyBasis::Resistivity,
                )
                .unwrap(),
            )
        })
        .collect();
    SpatialInduction::new(&solved, conductors, induction_limits(), 1e-10, &mut || {
        false
    })
    .unwrap()
}

#[derive(Clone, Copy)]
struct Reference {
    amplitude: f64,
    solid: f64,
    fluid: f64,
    joule: f64,
    reversible: f64,
    coolant_power: f64,
    impressed: f64,
    numerical: f64,
    internal: f64,
    parameter: f64,
}
impl Reference {
    fn state(self) -> [f64; 3] {
        [self.amplitude, self.solid, self.fluid]
    }
}

// Exact unit-cube Whitney integration gives K=4/mu and the charge-eliminated
// conductor mass M=3*sigma/40. The linear owner uses endpoint-averaged magnetic
// effort; the nonlinear owner's linear limit uses BE. No production matrices,
// thermal residual or convection implementation enter this scalar reference.
fn reference(
    old: [f64; 3],
    viscosity: f64,
    cp: f64,
    resistivity: f64,
    first: bool,
    nonlinear: bool,
    slope: f64,
) -> Reference {
    let [a0, t0, f0] = old;
    let mu0 = 1.0 + slope * (t0 - 300.0);
    let current = if first && !nonlinear {
        2.0 * AMPLITUDE
    } else {
        0.0
    };
    let evaluate = |t: f64| {
        let mu = 1.0 + slope * (t - 300.0);
        let m = 3.0 / (resistivity * (1.0 + 0.005 * (t - 300.0))) / 40.0;
        let a = if nonlinear {
            m * a0 / (m + 4.0 * DT / mu)
        } else {
            (m * a0 + DT * current - 2.0 * DT * a0 / mu0) / (m + 2.0 * DT / mu)
        };
        let joule = m * (a - a0).powi(2) / DT;
        let internal = 2.0 * a * a / mu * (1.0 + t * slope / mu);
        let internal0 = 2.0 * a0 * a0 / mu0 * (1.0 + t0 * slope / mu0);
        let parameter = 2.0 * a0 * (if nonlinear { a0 } else { a }) * (1.0 / mu - 1.0 / mu0);
        let reversible =
            -parameter - 2.0 * slope * (t * a * a / (mu * mu) - t0 * a0 * a0 / (mu0 * mu0));
        let numerical = if nonlinear {
            2.0 * (a - a0).powi(2) / mu
        } else {
            0.0
        };
        let impressed = current * (a - a0);
        let f = f0 + (joule + reversible - SOLID_CAPACITY * (t - t0)) / (MASS * cp);
        let film = 0.5 * (t + f);
        let k = 0.5 + 0.002 * (film - 300.0);
        let eta = viscosity * 0.001 * (1.0 - 0.01 * (film - 300.0));
        let h = 0.664 * (1000.0 * 0.1 / eta).sqrt() * (eta * cp / k).cbrt() * k / 0.1;
        let power = SURFACE * h * (t - f);
        let energy = SOLID_CAPACITY * (t - t0) + MASS * cp * (f - f0) + internal - internal0
            + numerical
            - impressed;
        near(energy, 0.0, 1e-10);
        (
            SOLID_CAPACITY * (t - t0) - joule - reversible + DT * power,
            Reference {
                amplitude: a,
                solid: t,
                fluid: f,
                joule,
                reversible,
                coolant_power: power,
                impressed,
                numerical,
                internal,
                parameter,
            },
        )
    };
    let (mut low, mut high) = (290.0, 310.0);
    assert!(evaluate(low).0 < 0.0 && evaluate(high).0 > 0.0);
    for _ in 0..100 {
        let t = 0.5 * (low + high);
        if evaluate(t).0 > 0.0 {
            high = t;
        } else {
            low = t;
        }
    }
    let (residual, result) = evaluate(0.5 * (low + high));
    near(residual, 0.0, 1e-10);
    result
}

#[test]
fn g1_g3_complete_coolant_changes_retained_3d_field_and_magnetic_thermal_work() {
    let run = |viscosity, cp, resistivity| {
        let card = source::absolute_affine_magnetic_card(
            resistivity,
            0.005,
            1.0,
            0.0004,
            ElectricalPropertyBasis::Resistivity,
        );
        let solid = solid_curve(&card);
        let (fluid_card, fluid_curve) = fluid_material(viscosity, cp);
        let fluid = liquid_law(&fluid_card);
        let thermal = initial_thermal(&solid, &fluid_curve, cp);
        let query = source::absolute_point(300.0);
        let fluid_query = fluid_point(290.0);
        let material = source::bundle_at(&card, &query, ElectricalPropertyBasis::Resistivity);
        let (positions, tetrahedra) = geometry();
        let problem = SpatialVectorMagnetostatics::new(
            positions,
            tetrahedra,
            vec![SpatialMagneticMaterial::Linear(material.clone()); 6],
            assembly_limits(),
            &mut || false,
        )
        .unwrap();
        let state = initial_linear(&problem);
        let diagonal = problem.mesh().edges.binary_search(&[0, 7]).unwrap();
        near(state.potential_wb()[diagonal], AMPLITUDE, 1e-10);
        let regions: Vec<_> = (0..6)
            .map(|tetrahedron| SpatialMagneticRegion::Cell {
                tetrahedron,
                cell: BODY,
            })
            .collect();
        let laws: Vec<_> = (0..6)
            .map(|tetrahedron| SpatialMagneticThermalLaw {
                tetrahedron,
                law: ThermalMagneticMaterial::new(&card, &material, "T", [280.0, 320.0]).unwrap(),
                maximum_flux_density_t: 6.1,
            })
            .collect();
        let cells = [
            (BODY, 0, &card, &solid, &query),
            (FLUID, 1, &fluid_card, &fluid_curve, &fluid_query),
        ]
        .map(|(cell, id, card, curve, context)| MaterialCell {
            cell,
            material_id: MaterialId(id),
            curve: curve.curve(),
            point: ThermalMaterialPoint {
                material: card,
                context,
                temperature_axis: "T",
                liquid_fraction_axis: None,
            },
            conductivity: None,
        });
        let bulk = MaterialNetwork {
            cells: &cells,
            resistors: &[],
            links: &[],
        };
        let electrical = ElectricalContactNetwork {
            bulk: &bulk,
            contacts: &[],
            electrical_bases: &[],
        };
        let heat = MaterialHeatNetwork {
            electrical: &electrical,
            thermal_contacts: &[],
            radiation: None,
        };
        let links = [binding()];
        let admitted = evaluate_source_link_convection(
            &thermal,
            &links,
            1,
            |_, _, t| fluid.at_temperature(t),
            || false,
        )
        .unwrap();
        let original_cells = thermal.cells().clone();
        let original_potential = state.potential_wb().to_vec();
        let gate = fs_exec::CancelGate::new();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        let mut result = [0.0; 3];
        pool.scope(|arena| {
            let cx = fs_exec::Cx::new(
                &gate,
                arena,
                fs_exec::StreamKey {
                    seed: 193,
                    kernel_id: 41,
                    tile: 0,
                    iteration: 0,
                },
                fs_exec::Budget::INFINITE,
                fs_exec::ExecMode::Deterministic,
            );
            assert!(heat.prepare(&cx, &thermal).is_err());
            let prepared = heat
                .prepare_with_link_conductances(&cx, &thermal, admitted.conductances_w_k())
                .unwrap();
            let network = SourceConvectionHeatNetwork {
                heat: &prepared,
                cells: &[],
                maximum_cells: 0,
                resolve_fluid: |_: &SourceConvectiveCell,
                                _: &DepositionSimulation<'_>,
                                t: Temperature| {
                    fluid.at_temperature(t)
                },
            };
            let provider = network.with_links(&links, 1, |_, _, t| fluid.at_temperature(t));
            let (mut accepted, mut field) = (thermal.clone(), state.clone());
            let mut expected_state = [AMPLITUDE, 300.0, 290.0];
            for tick in 1..=2 {
                let expected = reference(
                    expected_state,
                    viscosity,
                    cp,
                    resistivity,
                    tick == 1,
                    false,
                    0.0004,
                );
                let old_cells = accepted.cells().clone();
                let old_potential = field.potential_wb().to_vec();
                let next = heat_assembly_candidate(
                    &provider,
                    &cx,
                    &accepted,
                    f64::from(tick) * DT,
                    &field,
                    &[[0.0; 3]; 6],
                    &[[0.0; 3]; 8],
                    &regions,
                    1e-12,
                    &laws,
                    solve_controls(),
                    field_controls(),
                    assembly_limits(),
                    &[0.0; 2],
                )
                .unwrap();
                let report = next.report();
                let native = &report.coupled.electrical;
                near(
                    next.thermal().temperature_k(BODY).unwrap(),
                    expected.solid,
                    3e-7,
                );
                near(
                    next.thermal().temperature_k(FLUID).unwrap(),
                    expected.fluid,
                    3e-7,
                );
                near(
                    next.field().potential_wb()[diagonal],
                    expected.amplitude,
                    3e-8,
                );
                near(native.heat_j, expected.joule, 3e-7);
                near(native.impressed_work_j, expected.impressed, 3e-7);
                near(native.mechanical_work_j, 0.0, 1e-12);
                near(
                    report.coupled.source.reversible_heat_j,
                    expected.reversible,
                    3e-7,
                );
                near(
                    report.coupled.source.field_internal_energy_j[1],
                    expected.internal,
                    3e-7,
                );
                near(report.coupled.deposited_joule_j, expected.joule, 3e-7);
                near(
                    report.coupled.thermal.stored_energy_change_j,
                    expected.joule + expected.reversible,
                    3e-7,
                );
                near(report.external_heat_j, 0.0, 3e-7);
                near(report.total_energy_residual_j, 0.0, 3e-7);
                let coolant = &report
                    .coupled
                    .source
                    .network
                    .finite_convection
                    .as_ref()
                    .unwrap()
                    .samples()[0];
                near(coolant.heat_into_fluid_w, expected.coolant_power, 3e-6);
                near(
                    coolant.property_temperature_k,
                    0.5 * (expected.solid + expected.fluid),
                    3e-7,
                );
                assert_eq!(
                    coolant.material.material().card_identity(),
                    fluid_card.content_hash()
                );
                assert_eq!(coolant.material.material().query().axes()["p"], 100000.0);
                for conductor in next.field().conductors().iter().flatten() {
                    assert_eq!(conductor.material().card_identity(), card.content_hash());
                    near(
                        conductor.material().query().axes()["T"],
                        expected.solid,
                        3e-7,
                    );
                    near(
                        conductor.resistivity_ohm_m(),
                        resistivity * (1.0 + 0.005 * (expected.solid - 300.0)),
                        3e-10,
                    );
                }
                for local in &report.coupled.source.fields {
                    near(local.temperature_k[1], expected.solid, 3e-7);
                    near(local.flux_density_t[1], 2.0 * expected.amplitude, 3e-8);
                }
                assert!(std::ptr::eq(
                    next.field().problem().mesh(),
                    state.problem().mesh()
                ));
                assert_eq!(accepted.cells(), &old_cells);
                assert_eq!(field.potential_wb(), old_potential);
                assert_eq!(next.thermal().links(), thermal.links());
                result = [
                    expected.solid,
                    expected.fluid,
                    next.field().potential_wb()[diagonal],
                ];
                expected_state = expected.state();
                (accepted, field, _) = next.into_parts();
            }
            let saved_cells = accepted.cells().clone();
            let saved_field = field.potential_wb().to_vec();
            let queries = std::cell::Cell::new(0);
            let refused = SourceConvectionHeatNetwork {
                heat: &prepared,
                cells: &[],
                maximum_cells: 0,
                resolve_fluid: |_: &SourceConvectiveCell,
                                _: &DepositionSimulation<'_>,
                                t: Temperature| {
                    Ok::<_, &'static str>(fluid.at_temperature(t).unwrap())
                },
            };
            let failing = refused.with_links(&links, 1, |_, _, t| {
                queries.set(queries.get() + 1);
                if queries.get() > 1 {
                    Err("spatial coolant source refused")
                } else {
                    Ok(fluid.at_temperature(t).unwrap())
                }
            });
            assert!(matches!(
                heat_assembly_candidate(
                    &failing,
                    &cx,
                    &accepted,
                    3.0 * DT,
                    &field,
                    &[[0.0; 3]; 6],
                    &[[0.0; 3]; 8],
                    &regions,
                    1e-12,
                    &laws,
                    solve_controls(),
                    field_controls(),
                    assembly_limits(),
                    &[0.0; 2]
                ),
                Err(CircuitEnthalpyError::Source(
                    SpatialInductionHeatError::Heat(ConvectionHeatError::Convection(
                        SourceCellConvectionError::Source {
                            cell: FLUID,
                            source: "spatial coolant source refused",
                            ..
                        }
                    ))
                ))
            ));
            assert!(queries.get() > 1);
            let mut insufficient = solve_controls();
            insufficient.maximum_circuit_trials = 1;
            assert!(matches!(
                heat_assembly_candidate(
                    &provider,
                    &cx,
                    &accepted,
                    3.0 * DT,
                    &field,
                    &[[0.0; 3]; 6],
                    &[[0.0; 3]; 8],
                    &regions,
                    1e-12,
                    &laws,
                    insufficient,
                    field_controls(),
                    assembly_limits(),
                    &[0.0; 2]
                ),
                Err(CircuitEnthalpyError::CircuitBudget)
            ));
            gate.request();
            assert!(matches!(
                heat_assembly_candidate(
                    &provider,
                    &cx,
                    &accepted,
                    3.0 * DT,
                    &field,
                    &[[0.0; 3]; 6],
                    &[[0.0; 3]; 8],
                    &regions,
                    1e-12,
                    &laws,
                    solve_controls(),
                    field_controls(),
                    assembly_limits(),
                    &[0.0; 2]
                ),
                Err(CircuitEnthalpyError::Cancelled)
            ));
            assert_eq!(accepted.cells(), &saved_cells);
            assert_eq!(field.potential_wb(), saved_field);
            assert_eq!(thermal.cells(), &original_cells);
            assert_eq!(state.potential_wb(), original_potential);
            assert_eq!(state.time_s(), 0.0);
        });
        result
    };
    let baseline = run(1.0, 4000.0, 0.2);
    let viscous = run(4.0, 4000.0, 0.2);
    let capacity = run(1.0, 8000.0, 0.2);
    let conductor = run(1.0, 4000.0, 0.4);
    assert!(viscous[0] > baseline[0] && viscous[2] < baseline[2]);
    assert!(capacity[0] < baseline[0] && capacity[2] > baseline[2]);
    assert!(conductor[2] < baseline[2] - 0.1);
}

#[test]
fn g1_complete_coolant_keeps_nonlinear_field_numerical_loss_out_of_heat() {
    // Both fixed and affine-mu(T) limits keep native BE numerical loss out of
    // heat; the latter also closes magnetic parameter work and field entropy.
    let run = |viscosity, slope| {
        let (rho, cp) = (0.2, 4000.0);
        let card = source::absolute_affine_magnetic_card(
            rho,
            0.005,
            1.0,
            slope,
            ElectricalPropertyBasis::Resistivity,
        );
        let solid = solid_curve(&card);
        let (fluid_card, fluid_curve) = fluid_material(viscosity, cp);
        let fluid = liquid_law(&fluid_card);
        let thermal = initial_thermal(&solid, &fluid_curve, cp);
        let query = source::absolute_point(300.0);
        let fluid_query = fluid_point(290.0);
        let material = source::bundle_at(&card, &query, ElectricalPropertyBasis::Resistivity);
        let (positions, tetrahedra) = geometry();
        let problem = SpatialNonlinearVectorMagnetostatics::new(
            positions,
            tetrahedra,
            vec![SpatialNonlinearMaterial::Linear(material.clone()); 6],
            assembly_limits(),
            &mut || false,
        )
        .unwrap();
        let diagonal = problem.mesh().edges.binary_search(&[0, 7]).unwrap();
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
        let mut guess = vec![0.0; problem.mesh().edges.len()];
        guess[diagonal] = AMPLITUDE;
        let field = problem
            .solve(
                &[[8.0 * AMPLITUDE; 3]; 6],
                &boundary,
                &guess,
                NonlinearVectorControls {
                    field: vector_controls(),
                    relative_residual_tolerance: 1e-12,
                    maximum_newton_iterations: 40,
                    maximum_backtracks: 30,
                },
                &mut || false,
            )
            .unwrap();
        let thermal_law =
            ThermalMagneticMaterial::new(&card, &material, "T", [280.0, 320.0]).unwrap();
        let laws: Vec<_> = (0..6)
            .map(|tetrahedron| SpatialMagneticThermalLaw {
                tetrahedron,
                law: MagneticConstitutiveThermalLaw::Linear(thermal_law.clone()),
                maximum_flux_density_t: 6.1,
            })
            .collect();
        let conductor = ResolvedConductorState::try_from_material_with_basis(
            material,
            ElectricalPropertyBasis::Resistivity,
        )
        .unwrap();
        let state = SpatialNonlinearInduction::new(
            &field,
            vec![Some(conductor); 6],
            induction_limits(),
            1e-10,
            &mut || false,
        )
        .unwrap();
        let cells = [
            (BODY, 0, &card, &solid, &query),
            (FLUID, 1, &fluid_card, &fluid_curve, &fluid_query),
        ]
        .map(|(cell, id, card, curve, context)| MaterialCell {
            cell,
            material_id: MaterialId(id),
            curve: curve.curve(),
            point: ThermalMaterialPoint {
                material: card,
                context,
                temperature_axis: "T",
                liquid_fraction_axis: None,
            },
            conductivity: None,
        });
        let bulk = MaterialNetwork {
            cells: &cells,
            resistors: &[],
            links: &[],
        };
        let electrical = ElectricalContactNetwork {
            bulk: &bulk,
            contacts: &[],
            electrical_bases: &[],
        };
        let heat = MaterialHeatNetwork {
            electrical: &electrical,
            thermal_contacts: &[],
            radiation: None,
        };
        let regions: Vec<_> = (0..6)
            .map(|tetrahedron| SpatialMagneticRegion::Cell {
                tetrahedron,
                cell: BODY,
            })
            .collect();
        let links = [binding()];
        let admitted = evaluate_source_link_convection(
            &thermal,
            &links,
            1,
            |_, _, t| fluid.at_temperature(t),
            || false,
        )
        .unwrap();
        let gate = fs_exec::CancelGate::new();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        let mut result = [0.0; 2];
        pool.scope(|arena| {
            let cx = fs_exec::Cx::new(
                &gate,
                arena,
                fs_exec::StreamKey {
                    seed: 197,
                    kernel_id: 43,
                    tile: 0,
                    iteration: 0,
                },
                fs_exec::Budget::INFINITE,
                fs_exec::ExecMode::Deterministic,
            );
            let prepared = heat
                .prepare_with_link_conductances(&cx, &thermal, admitted.conductances_w_k())
                .unwrap();
            let network = SourceConvectionHeatNetwork {
                heat: &prepared,
                cells: &[],
                maximum_cells: 0,
                resolve_fluid: |_: &SourceConvectiveCell,
                                _: &DepositionSimulation<'_>,
                                t: Temperature| {
                    fluid.at_temperature(t)
                },
            };
            let provider = network.with_links(&links, 1, |_, _, t| fluid.at_temperature(t));
            let (mut accepted, mut magnetic) = (thermal.clone(), state.clone());
            let mut expected_state = [AMPLITUDE, 300.0, 290.0];
            for tick in 1..=2 {
                let expected = reference(expected_state, viscosity, cp, rho, false, true, slope);
                let before = accepted.cells().clone();
                let old_field = magnetic.potential_wb().to_vec();
                let field_solve = SpatialNonlinearInductionControls {
                    induction: field_controls(),
                    maximum_newton_iterations: 40,
                    maximum_backtracks: 30,
                };
                let next = if slope == 0.0 {
                    nonlinear_heat_assembly_candidate(
                        &provider,
                        &cx,
                        &accepted,
                        f64::from(tick) * DT,
                        &magnetic,
                        &[[0.0; 3]; 6],
                        &[[0.0; 3]; 8],
                        &regions,
                        1e-12,
                        solve_controls(),
                        field_solve,
                        &[0.0; 2],
                    )
                } else {
                    nonlinear_thermal_heat_assembly_candidate(
                        &provider,
                        &cx,
                        &accepted,
                        f64::from(tick) * DT,
                        &magnetic,
                        &[[0.0; 3]; 6],
                        &[[0.0; 3]; 8],
                        &regions,
                        1e-12,
                        &laws,
                        solve_controls(),
                        field_solve,
                        &[0.0; 2],
                    )
                }
                .unwrap();
                let report = next.report();
                let native = &report.coupled.electrical;
                near(
                    next.field().potential_wb()[diagonal],
                    expected.amplitude,
                    3e-8,
                );
                near(
                    next.thermal().temperature_k(BODY).unwrap(),
                    expected.solid,
                    3e-7,
                );
                near(
                    next.thermal().temperature_k(FLUID).unwrap(),
                    expected.fluid,
                    3e-7,
                );
                near(native.heat_j, expected.joule, 3e-7);
                near(native.numerical_dissipation_j, expected.numerical, 3e-7);
                assert!(native.numerical_dissipation_j > 0.1);
                near(native.parameter_work_j, expected.parameter, 3e-7);
                near(
                    report.coupled.source.reversible_heat_j,
                    expected.reversible,
                    3e-7,
                );
                near(
                    report.coupled.source.field_internal_energy_j[1],
                    expected.internal,
                    3e-7,
                );
                near(report.coupled.deposited_joule_j, expected.joule, 3e-7);
                near(
                    report.coupled.thermal.stored_energy_change_j,
                    expected.joule + expected.reversible,
                    3e-7,
                );
                if slope != 0.0 {
                    assert!(expected.reversible.abs() > 1e-3);
                    for field in &report.coupled.source.fields {
                        near(field.temperature_k[1], expected.solid, 3e-7);
                        near(field.flux_density_t[1], 2.0 * expected.amplitude, 3e-8);
                    }
                }
                near(report.external_heat_j, 0.0, 3e-7);
                near(report.total_energy_residual_j, 0.0, 3e-7);
                let sample = &report
                    .coupled
                    .source
                    .network
                    .finite_convection
                    .as_ref()
                    .unwrap()
                    .samples()[0];
                near(sample.heat_into_fluid_w, expected.coolant_power, 3e-6);
                assert_eq!(
                    sample.material.material().card_identity(),
                    fluid_card.content_hash()
                );
                for conductor in next.field().conductors().iter().flatten() {
                    near(
                        conductor.material().query().axes()["T"],
                        expected.solid,
                        3e-7,
                    );
                }
                assert_eq!(accepted.cells(), &before);
                assert_eq!(magnetic.potential_wb(), old_field);
                assert!(std::ptr::eq(next.field().problem(), state.problem()));
                result = [
                    next.thermal().temperature_k(BODY).unwrap(),
                    next.field().potential_wb()[diagonal],
                ];
                expected_state = expected.state();
                (accepted, magnetic, _) = next.into_parts();
            }
        });
        result
    };
    let baseline = run(1.0, 0.0);
    let viscous = run(4.0, 0.0);
    assert!(viscous[0] > baseline[0] && viscous[1] < baseline[1]);
    let thermal = run(1.0, 0.0004);
    let thermal_viscous = run(4.0, 0.0004);
    assert!((thermal[0] - baseline[0]).abs() > 0.01);
    assert!((thermal[1] - baseline[1]).abs() > 1e-5);
    assert!(thermal_viscous[0] > thermal[0] && thermal_viscous[1] < thermal[1]);
}
