//! Actual perturbed enthalpy solves with upstream heating and bypass mixing.
//! The plateau's history seed is independent of the temperature observation.

use fs_airflow::conjugate::{AirSegment, ConjugateConfig, SolidRegionState};
use fs_airflow::graph::thermal::coupled_transport::sensitivity::enthalpy::CoupledEnthalpyLinearization;
use fs_airflow::graph::thermal::coupled_transport::sensitivity::{
    CoupledObjective, CoupledSensitivityError, InterfaceSolveConfig,
};
use fs_airflow::graph::thermal::coupled_transport::{
    CoupledTransportSolution, solve_coupled_transport_iqn,
};
use fs_airflow::graph::thermal::transport::{
    BranchThermalModel, TransportAir, TransportConfig, TransportInlet, TransportNetwork,
};
use fs_airflow::graph::{FixedPressure, GraphBranch, GraphSolveConfig, LossGraph};
use fs_airflow::{LossElement, LossResistance, SourceProvenance, ToleranceBasis};
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::ContentHash;
use fs_conduction::transient::enthalpy::adjoint::EnthalpyStepLinearization;
use fs_conduction::transient::enthalpy::{
    EnthalpyBackwardEuler, EnthalpyBudget, EnthalpyStepConfig, EnthalpyStepSolution,
};
use fs_conduction::{
    AmbientRadiationConfig, AmbientRadiationPatch, ConductionMesh, ConductionProblem,
    ConductivityModel, ConductivityTable, EMISSIVITY_DIMS, LinearConfig,
    SURFACE_EMISSIVITY_PROPERTY, ScalarField, SurfaceEmissivity, ThermalBc, ThermalBoundary,
    ThermalBoundaryBuilder,
};
use fs_couple::iqn_ils::IqnIlsConfig;
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, SelectionPolicy, UncertaintyModel,
};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
use fs_qty::{Density, Pressure, Temperature, VolumetricFlowRate};
use fs_rep_mesh::TetComplex;

const DT: f64 = 0.2;
const NAMES: [&str; 2] = ["first", "last"];
const CARRY: [f64; 4] = [0.0, 0.07, 0.0, -0.015];

fn with_gate<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 29,
                kernel_id: 733,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}
fn config() -> EnthalpyStepConfig {
    let mut c = EnthalpyStepConfig {
        newton: Default::default(),
        max_newton_iterations: 32,
        energy_tolerance_j: 1e-8,
    };
    c.newton.absolute_tolerance = 1e-10;
    c.newton.relative_tolerance = 1e-12;
    c.newton.linear_restart = 4;
    c.newton.max_linear_cycles = 8;
    c.newton.forcing_minimum = 1e-12;
    c.newton.forcing_maximum = 1e-4;
    c
}
fn linear() -> LinearConfig {
    LinearConfig {
        tolerance: 1e-12,
        max_iterations: 64,
        restart: 4,
    }
}
fn gate() -> ConjugateConfig {
    ConjugateConfig {
        temperature_tolerance_k: 1e-10,
        max_iterations: 128,
        balance_tolerance_w: 1e-8,
        balance_relative_tolerance: 1e-10,
        ..ConjugateConfig::default()
    }
}
fn interface() -> InterfaceSolveConfig {
    InterfaceSolveConfig {
        max_iterations: 128,
        absolute_tolerance: 1e-12,
        relative_tolerance: 1e-11,
        relaxation: 1.0,
    }
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 3e-7 * (1.0 + expected.abs()),
        "implicit {actual:.12e}, perturbed solve {expected:.12e}"
    );
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}

#[derive(Clone)]
struct Controls {
    old: [f64; 4],
    source: [f64; 4],
    inlet: [f64; 2],
    htc: [f64; 2],
    flow_scale: f64,
}
impl Default for Controls {
    fn default() -> Self {
        Self {
            old: [600.0, 1800.0, 3400.0, 5000.0],
            source: [12.0, -3.0, 8.0, 5.0],
            inlet: [330.0, 290.0],
            htc: [20.0, 10.0],
            flow_scale: 1.0,
        }
    }
}

fn with_network<R>(cx: &Cx<'_>, c: &Controls, f: impl FnOnce(&TransportNetwork<'_>) -> R) -> R {
    let edge = |name, from, to, q: f64| GraphBranch {
        from,
        to,
        loss: LossElement::new(
            name,
            LossResistance::new(1.0 / (q * q)),
            0.0,
            SourceProvenance::new("analytic fixture", "enthalpy-adjoint"),
            ToleranceBasis::Analytic,
        )
        .unwrap(),
    };
    let flow = LossGraph::new(
        4,
        vec![
            edge("first", 0, 2, 0.005),
            edge("bypass", 1, 2, 0.002),
            edge("last", 2, 3, 0.007),
        ],
    )
    .unwrap()
    .solve(
        &[
            FixedPressure {
                node: 0,
                pressure: Pressure::new(2.0 * c.flow_scale * c.flow_scale),
            },
            FixedPressure {
                node: 1,
                pressure: Pressure::new(2.0 * c.flow_scale * c.flow_scale),
            },
            FixedPressure {
                node: 3,
                pressure: Pressure::new(0.0),
            },
        ],
        GraphSolveConfig {
            max_sweeps: 4096,
            max_node_iterations: 80,
            absolute_flow_tolerance: VolumetricFlowRate::new(1e-13),
            relative_flow_tolerance: 1e-12,
        },
        cx,
    )
    .unwrap();
    let network = TransportNetwork::new(
        cx,
        &flow,
        TransportAir {
            density: Density::new(1.0),
            specific_heat_j_kg_k: 1000.0,
        },
        vec![
            BranchThermalModel::Exchange(vec![AirSegment::new(NAMES[0], 0.5, c.htc[0]).unwrap()]),
            BranchThermalModel::Adiabatic,
            BranchThermalModel::Exchange(vec![
                AirSegment::new(NAMES[1], 1.0 + 3.0_f64.sqrt() / 2.0, c.htc[1]).unwrap(),
            ]),
        ],
        &[
            TransportInlet {
                node: 0,
                temperature: Temperature::new(c.inlet[0]),
            },
            TransportInlet {
                node: 1,
                temperature: Temperature::new(c.inlet[1]),
            },
        ],
        TransportConfig {
            absolute_flow_tolerance: VolumetricFlowRate::new(1e-12),
            relative_flow_tolerance: 1e-10,
            absolute_heat_tolerance_w: 1e-7,
            relative_heat_tolerance: 1e-9,
        },
    )
    .unwrap();
    f(&network)
}

struct Fixture {
    mesh: ConductionMesh,
    curve: EquilibriumEnthalpyPhaseCurve,
    material: ConductivityModel,
    patches: Vec<AmbientRadiationPatch>,
}
impl Fixture {
    fn new(radiation: bool) -> Self {
        let mesh = ConductionMesh::new(
            TetComplex::from_tets(4, vec![[0, 1, 2, 3]]),
            vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ],
        )
        .unwrap();
        let curve = EquilibriumEnthalpyPhaseCurve::try_new(
            ContentHash([0x72; 32]),
            [
                (0.0, 250.0, 0.0),
                (1000.0, 350.0, 0.0),
                (3000.0, 350.0, 1.0),
                (4000.0, 450.0, 1.0),
                (6000.0, 550.0, 1.0),
            ]
            .into_iter()
            .map(|(h, t, f)| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: h,
                temperature_k: t,
                liquid_mass_fraction: f,
                bulk_density_kg_m3: 10.0,
            })
            .collect(),
        )
        .unwrap();
        Self {
            mesh,
            curve,
            material: ConductivityModel::isotropic(
                ConductivityTable::declared_curve(vec![(240.0, 2.0), (560.0, 8.0)]).unwrap(),
            ),
            patches: if radiation {
                vec![radiation_patch()]
            } else {
                Vec::new()
            },
        }
    }
    fn stepper(&self, cx: &Cx<'_>) -> EnthalpyBackwardEuler<'_, '_> {
        EnthalpyBackwardEuler::uniform(
            cx,
            &self.mesh,
            &self.curve,
            10.0,
            EnthalpyBudget {
                max_vertices: 4,
                max_elements: 1,
            },
        )
        .unwrap()
    }
    fn boundary(&self, c: &Controls, refs: &[f64]) -> ThermalBoundary {
        ThermalBoundaryBuilder::new(&self.mesh)
            .region(
                NAMES[0],
                |f| f.centroid[0] == 0.0,
                ThermalBc::robin(c.htc[0], refs[0]).unwrap(),
            )
            .unwrap()
            .region(
                NAMES[1],
                |f| f.centroid[0] > 0.0,
                ThermalBc::robin(c.htc[1], refs[1]).unwrap(),
            )
            .unwrap()
            .finish()
            .unwrap()
    }
    fn problem<'a>(
        &'a self,
        boundary: &'a ThermalBoundary,
        source: &'a ScalarField,
    ) -> ConductionProblem<'a> {
        ConductionProblem {
            mesh: &self.mesh,
            boundary,
            material: &self.material,
            source,
            element_materials: None,
        }
    }
    fn solve(
        &self,
        cx: &Cx<'_>,
        c: &Controls,
        network: &TransportNetwork<'_>,
    ) -> (EnthalpyStepSolution, CoupledTransportSolution) {
        let source = ScalarField::nodal("deposition", 4, c.source.to_vec()).unwrap();
        let stepper = self.stepper(cx);
        let mut retained = None;
        let coupled = solve_coupled_transport_iqn(
            cx,
            network,
            &gate(),
            IqnIlsConfig::default(),
            |cx, refs| {
                let boundary = self.boundary(c, refs);
                // Every nonlinear/air iteration solves the SAME physical step.
                let (solid, fluxes) = if self.patches.is_empty() {
                    let s = stepper
                        .advance(
                            cx,
                            self.problem(&boundary, &source),
                            None,
                            &c.old,
                            DT,
                            config(),
                        )
                        .unwrap();
                    let fluxes = s.robin_fluxes.clone();
                    (s, fluxes)
                } else {
                    let s = stepper
                        .advance_with_ambient_radiation(
                            cx,
                            self.problem(&boundary, &source),
                            None,
                            &c.old,
                            DT,
                            config(),
                            &self.patches,
                            AmbientRadiationConfig {
                                max_iterations: 128,
                                ..AmbientRadiationConfig::default()
                            },
                        )
                        .unwrap();
                    (s.conduction, s.convective_robin_fluxes)
                };
                let states = NAMES
                    .iter()
                    .map(|name| {
                        SolidRegionState::from_robin_flux(
                            fluxes.iter().find(|r| r.region == *name).unwrap(),
                        )
                    })
                    .collect();
                retained = Some(solid);
                Ok(states)
            },
        )
        .unwrap();
        (retained.unwrap(), coupled)
    }
    fn bind(
        &self,
        cx: &Cx<'_>,
        c: &Controls,
        solid: EnthalpyStepSolution,
        refs: &[f64],
    ) -> EnthalpyStepLinearization<'_> {
        let source = ScalarField::nodal("deposition", 4, c.source.to_vec()).unwrap();
        let boundary = self.boundary(c, refs);
        let stepper = self.stepper(cx);
        if self.patches.is_empty() {
            stepper
                .linearize_accepted(
                    cx,
                    self.problem(&boundary, &source),
                    None,
                    &c.old,
                    DT,
                    config(),
                    solid,
                )
                .unwrap()
        } else {
            stepper
                .linearize_accepted_with_ambient_radiation(
                    cx,
                    self.problem(&boundary, &source),
                    None,
                    &c.old,
                    DT,
                    config(),
                    &self.patches,
                    solid,
                    8,
                )
                .unwrap()
                .into_transport()
        }
    }
    fn objective(&self, cx: &Cx<'_>, c: &Controls, o: &CoupledObjective, carry: &[f64]) -> f64 {
        with_network(cx, c, |network| {
            let (solid, coupled) = self.solve(cx, c, network);
            let mut value = dot(&solid.specific_enthalpy_j_kg, carry);
            value += solid
                .temperature
                .iter()
                .zip(&o.nodal_temperatures)
                .map(|(t, w)| w * (t - 350.0))
                .sum::<f64>();
            for ((s, w), q) in coupled
                .solid
                .iter()
                .zip(&o.wall_temperatures)
                .zip(&o.solid_heat_rates)
            {
                value += w * (s.mean_wall_temperature_k - 350.0) + q * s.heat_rate_w;
            }
            let a = &coupled.transport;
            value += dot(&a.reference_temperatures_k, &o.air.references)
                + o.air.wall_heat_rate * a.wall_heat_rate_w
                + o.air.external_heat_gain * a.external_heat_gain_w
                + o.air.heat_imbalance * a.heat_imbalance_w
                + o.air.hydraulic_energy_defect * a.hydraulic_energy_defect_w;
            value += a
                .branches
                .iter()
                .zip(&o.air.branch_outlets)
                .map(|(b, w)| w * b.outlet_temperature_k.unwrap_or(0.0))
                .sum::<f64>();
            value += a
                .node_temperatures_k
                .iter()
                .zip(&o.air.node_temperatures)
                .map(|(t, w)| w * t.unwrap_or(0.0))
                .sum::<f64>();
            value
        })
    }
}

fn radiation_patch() -> AmbientRadiationPatch {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS),
            value: PropertyValue::Scalar {
                value: 0.8,
                dims: EMISSIVITY_DIMS,
            },
            validity: ValidityDomain::unconstrained().with("T", 240.0, 560.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: Vec::new(),
            provenance: Provenance {
                source: "synthetic coupled gradient fixture".into(),
                license: "internal-test-use".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic gray surface".into(),
            phase: "solid".into(),
            process: "numerical reference".into(),
            revision: 0,
        },
        claims,
        Vec::new(),
    )
    .unwrap();
    let emissivity =
        SurfaceEmissivity::from_card(NAMES[0], &card, 350.0, SelectionPolicy::SingleClaimOnly)
            .unwrap();
    AmbientRadiationPatch::new(NAMES[0], emissivity, 300.0).unwrap()
}
fn objective(c: &CoupledEnthalpyLinearization<'_, '_, '_, '_>) -> CoupledObjective {
    let mut o = c.zero_objective();
    o.nodal_temperatures = vec![0.7, -0.4, 0.2, 1.1];
    o.wall_temperatures = vec![0.3, -0.2];
    o.solid_heat_rates = vec![0.03, -0.02];
    o.air.references = vec![0.2, -0.1];
    o.air.branch_outlets[2] = 0.4;
    o.air.node_temperatures[2] = -0.15;
    o.air.wall_heat_rate = 0.01;
    o.air.external_heat_gain = -0.02;
    o
}

#[test]
fn mixed_phase_coupled_pullback_matches_history_source_inlet_and_htc_perturbations() {
    let fixture = Fixture::new(false);
    let c = Controls::default();
    with_gate(&CancelGate::new_clock_free(), |cx| {
        with_network(cx, &c, |network| {
            let (solid, primal) = fixture.solve(cx, &c, network);
            assert_eq!(solid.temperature[1], 350.0);
            assert!(solid.temperature[0] < 350.0 && solid.temperature[2] > 350.0);
            let linearization = fixture.bind(cx, &c, solid, &primal.reference_temperatures_k);
            let response = linearization.robin_response(cx, &NAMES, linear()).unwrap();
            let coupled =
                CoupledEnthalpyLinearization::new(cx, network, &response, &gate()).unwrap();
            let o = objective(&coupled);
            let g = coupled
                .pullback_iqn(cx, &o, &CARRY, interface(), IqnIlsConfig::default())
                .unwrap();
            assert!(g.interface_residual < 1e-10);
            for i in 0..4 {
                let mut plus = c.clone();
                let mut minus = c.clone();
                plus.old[i] += 0.05;
                minus.old[i] -= 0.05;
                close(
                    g.solid.previous_specific_enthalpy[i],
                    (fixture.objective(cx, &plus, &o, &CARRY)
                        - fixture.objective(cx, &minus, &o, &CARRY))
                        / 0.1,
                );
                let mut plus = c.clone();
                let mut minus = c.clone();
                plus.source[i] += 0.2;
                minus.source[i] -= 0.2;
                close(
                    g.solid.source_density[i],
                    (fixture.objective(cx, &plus, &o, &CARRY)
                        - fixture.objective(cx, &minus, &o, &CARRY))
                        / 0.4,
                );
            }
            for i in 0..2 {
                let mut plus = c.clone();
                let mut minus = c.clone();
                plus.inlet[i] += 0.01;
                minus.inlet[i] -= 0.01;
                close(
                    g.inlets[i],
                    (fixture.objective(cx, &plus, &o, &CARRY)
                        - fixture.objective(cx, &minus, &o, &CARRY))
                        / 0.02,
                );
                let mut plus = c.clone();
                let mut minus = c.clone();
                plus.htc[i] *= 1e-4_f64.exp();
                minus.htc[i] *= (-1e-4_f64).exp();
                close(
                    g.log_htc[i],
                    (fixture.objective(cx, &plus, &o, &CARRY)
                        - fixture.objective(cx, &minus, &o, &CARRY))
                        / 2e-4,
                );
            }
            let mut temperature_only = coupled.zero_objective();
            temperature_only.nodal_temperatures[3] = 1.0;
            let total = coupled
                .pullback(cx, &temperature_only, &[0.0; 4], interface())
                .unwrap();
            let frozen = response
                .pullback(
                    cx,
                    &[0.0; 4],
                    &temperature_only.nodal_temperatures,
                    &[0.0; 2],
                    &[0.0; 2],
                )
                .unwrap();
            assert!(
                total
                    .solid
                    .source_density
                    .iter()
                    .zip(&frozen.transport.source_density)
                    .any(|(a, b)| (a - b).abs() > 1e-7)
            );
        });
    });
}

#[test]
fn latent_history_carry_survives_when_all_temperature_seeds_vanish() {
    let fixture = Fixture::new(false);
    let c = Controls {
        old: [2000.0; 4],
        ..Controls::default()
    };
    with_gate(&CancelGate::new_clock_free(), |cx| {
        with_network(cx, &c, |network| {
            let (solid, primal) = fixture.solve(cx, &c, network);
            assert_eq!(solid.temperature, [350.0; 4]);
            let linearization = fixture.bind(cx, &c, solid, &primal.reference_temperatures_k);
            let response = linearization.robin_response(cx, &NAMES, linear()).unwrap();
            let coupled =
                CoupledEnthalpyLinearization::new(cx, network, &response, &gate()).unwrap();
            let mut o = coupled.zero_objective();
            o.nodal_temperatures.fill(1.0);
            let temperature = coupled.pullback(cx, &o, &[0.0; 4], interface()).unwrap();
            assert!(
                temperature
                    .solid
                    .previous_specific_enthalpy
                    .iter()
                    .chain(&temperature.solid.source_density)
                    .chain(&temperature.inlets)
                    .all(|v| *v == 0.0)
            );
            let o = coupled.zero_objective();
            let carry = [0.0, 1.0, 0.0, 0.0];
            let h = coupled
                .pullback_iqn(cx, &o, &carry, interface(), IqnIlsConfig::default())
                .unwrap();
            for (a, b) in h.solid.previous_specific_enthalpy.iter().zip(carry) {
                close(*a, b);
            }
            assert!(h.solid.source_density.iter().all(|v| *v > 0.0));
            let mut plus = c.clone();
            let mut minus = c.clone();
            plus.source[2] += 0.2;
            minus.source[2] -= 0.2;
            close(
                h.solid.source_density[2],
                (fixture.objective(cx, &plus, &o, &carry)
                    - fixture.objective(cx, &minus, &o, &carry))
                    / 0.4,
            );
        })
    });
}

#[test]
fn radiation_keeps_full_state_feedback_but_exposes_only_convection_to_air() {
    let fixture = Fixture::new(true);
    let c = Controls::default();
    with_gate(&CancelGate::new_clock_free(), |cx| {
        with_network(cx, &c, |network| {
            let (solid, primal) = fixture.solve(cx, &c, network);
            let combined_heat = solid.robin_out_w;
            let linearization = fixture.bind(cx, &c, solid, &primal.reference_temperatures_k);
            let response = linearization.robin_response(cx, &NAMES, linear()).unwrap();
            assert_eq!(response.ports().len(), 2);
            for (i, p) in response.ports().iter().enumerate() {
                assert_eq!(p.htc_w_m2_k, c.htc[i]);
                assert_eq!(p.reference_k, primal.reference_temperatures_k[i]);
            }
            let convective_heat: f64 = response.robin_fluxes().iter().map(|r| r.heat_rate_w).sum();
            close(convective_heat, primal.transport.wall_heat_rate_w);
            assert!((combined_heat - convective_heat).abs() > 10.0);
            let coupled =
                CoupledEnthalpyLinearization::new(cx, network, &response, &gate()).unwrap();
            let o = objective(&coupled);
            let g = coupled
                .pullback_iqn(cx, &o, &CARRY, interface(), IqnIlsConfig::default())
                .unwrap();
            let dh = [0.4, -0.2, 0.3, -0.5];
            let dq = [3.0, -2.0, 5.0, 1.0];
            let di = [0.3, -0.4];
            let dl = [0.2, -0.1];
            let predicted = dot(&dh, &g.solid.previous_specific_enthalpy)
                + dot(&dq, &g.solid.source_density)
                + dot(&di, &g.inlets[..2])
                + dot(&dl, &g.log_htc);
            let evaluate = |epsilon: f64| {
                let mut shifted = c.clone();
                for i in 0..4 {
                    shifted.old[i] += epsilon * dh[i];
                    shifted.source[i] += epsilon * dq[i];
                }
                for i in 0..2 {
                    shifted.inlet[i] += epsilon * di[i];
                    shifted.htc[i] *= (epsilon * dl[i]).exp();
                }
                fixture.objective(cx, &shifted, &o, &CARRY)
            };
            close(
                predicted,
                (-evaluate(0.02) + 8.0 * evaluate(0.01) - 8.0 * evaluate(-0.01) + evaluate(-0.02))
                    / 0.12,
            );
        })
    });
}

#[test]
fn uniform_flow_scale_matches_resolved_mixed_air_and_radiation_with_direct_heat_terms() {
    let fixture = Fixture::new(true);
    let c = Controls {
        flow_scale: 1.3,
        ..Controls::default()
    };
    with_gate(&CancelGate::new_clock_free(), |cx| {
        with_network(cx, &c, |network| {
            let (solid, primal) = fixture.solve(cx, &c, network);
            assert_eq!(solid.temperature[1], 350.0);
            let linearization = fixture.bind(cx, &c, solid, &primal.reference_temperatures_k);
            let response = linearization.robin_response(cx, &NAMES, linear()).unwrap();
            let coupled =
                CoupledEnthalpyLinearization::new(cx, network, &response, &gate()).unwrap();
            // Includes explicit air watt objectives as well as mixed-air,
            // temperature, convection and latent-history contributions.
            let o = objective(&coupled);
            let g = coupled
                .pullback_flow_scale_iqn(cx, &o, &CARRY, interface(), IqnIlsConfig::default())
                .unwrap();
            let ordinary = coupled
                .pullback_iqn(cx, &o, &CARRY, interface(), IqnIlsConfig::default())
                .unwrap();
            assert_eq!(g.thermal.solid, ordinary.solid);
            assert_eq!(g.thermal.inlets, ordinary.inlets);
            assert_eq!(g.thermal.log_htc, ordinary.log_htc);
            assert_eq!(g.thermal.interface_adjoint, ordinary.interface_adjoint);
            assert_eq!(
                g.thermal.solid_krylov_iterations,
                ordinary.solid_krylov_iterations
            );
            let stationary = coupled
                .pullback_flow_scale(cx, &o, &CARRY, interface())
                .unwrap();
            close(g.log_flow_scale, stationary.log_flow_scale);
            let evaluate = |epsilon: f64| {
                let mut shifted = c.clone();
                // Quadratic hydraulic losses at pressure scaled by s^2
                // independently produce every branch/external flow times s.
                // Each evaluation resolves h, T, radiation and mixed air from
                // the same original history; no primal wall is frozen.
                shifted.flow_scale *= epsilon.exp();
                fixture.objective(cx, &shifted, &o, &CARRY)
            };
            let epsilon = 0.005;
            close(
                g.log_flow_scale,
                (-evaluate(2.0 * epsilon) + 8.0 * evaluate(epsilon) - 8.0 * evaluate(-epsilon)
                    + evaluate(-2.0 * epsilon))
                    / (12.0 * epsilon),
            );
            assert!(g.log_flow_scale.abs() > 1e-3);
        })
    });
}

#[test]
fn uniform_flow_scale_preserves_latent_history_when_temperature_response_is_zero() {
    let fixture = Fixture::new(true);
    let c = Controls {
        old: [2000.0; 4],
        flow_scale: 0.7,
        ..Controls::default()
    };
    with_gate(&CancelGate::new_clock_free(), |cx| {
        with_network(cx, &c, |network| {
            let (solid, primal) = fixture.solve(cx, &c, network);
            assert_eq!(solid.temperature, [350.0; 4]);
            let linearization = fixture.bind(cx, &c, solid, &primal.reference_temperatures_k);
            let response = linearization.robin_response(cx, &NAMES, linear()).unwrap();
            let coupled =
                CoupledEnthalpyLinearization::new(cx, network, &response, &gate()).unwrap();
            let mut temperature = coupled.zero_objective();
            temperature.nodal_temperatures.fill(1.0);
            let t = coupled
                .pullback_flow_scale(cx, &temperature, &[0.0; 4], interface())
                .unwrap();
            assert_eq!(t.log_flow_scale, 0.0);
            let o = coupled.zero_objective();
            let carry = [0.0, 1.0, 0.0, 0.0];
            let h = coupled
                .pullback_flow_scale_iqn(cx, &o, &carry, interface(), IqnIlsConfig::default())
                .unwrap();
            assert!(h.log_flow_scale.abs() > 1.0);
            for (actual, expected) in h.thermal.solid.previous_specific_enthalpy.iter().zip(carry) {
                close(*actual, expected);
            }
            let epsilon = 1e-4_f64;
            let mut plus = c.clone();
            let mut minus = c.clone();
            plus.flow_scale *= epsilon.exp();
            minus.flow_scale *= (-epsilon).exp();
            close(
                h.log_flow_scale,
                (fixture.objective(cx, &plus, &o, &carry)
                    - fixture.objective(cx, &minus, &o, &carry))
                    / (2.0 * epsilon),
            );
        })
    });
}

#[test]
fn cancellation_and_exhausted_interface_work_publish_no_gradient_and_binding_is_reusable() {
    let fixture = Fixture::new(false);
    let c = Controls::default();
    with_gate(&CancelGate::new_clock_free(), |cx| {
        with_network(cx, &c, |network| {
            let (solid, primal) = fixture.solve(cx, &c, network);
            let linearization = fixture.bind(cx, &c, solid, &primal.reference_temperatures_k);
            let response = linearization.robin_response(cx, &NAMES, linear()).unwrap();
            let coupled =
                CoupledEnthalpyLinearization::new(cx, network, &response, &gate()).unwrap();
            let o = objective(&coupled);
            let baseline = coupled
                .pullback_iqn(cx, &o, &CARRY, interface(), IqnIlsConfig::default())
                .unwrap();
            assert!(matches!(
                coupled.pullback(
                    cx,
                    &o,
                    &CARRY,
                    InterfaceSolveConfig {
                        max_iterations: 1,
                        ..interface()
                    }
                ),
                Err(CoupledSensitivityError::DidNotConverge { .. })
            ));
            assert!(matches!(
                coupled.pullback_flow_scale_iqn(
                    cx,
                    &o,
                    &CARRY,
                    InterfaceSolveConfig {
                        max_iterations: 1,
                        ..interface()
                    },
                    IqnIlsConfig::default(),
                ),
                Err(CoupledSensitivityError::DidNotConverge { .. })
            ));
            let cancelled = CancelGate::new_clock_free();
            cancelled.request();
            with_gate(&cancelled, |cx| {
                assert!(matches!(
                    coupled.pullback_iqn(cx, &o, &CARRY, interface(), IqnIlsConfig::default()),
                    Err(CoupledSensitivityError::Interrupted)
                ));
                assert!(matches!(
                    coupled.pullback_flow_scale(cx, &o, &CARRY, interface()),
                    Err(CoupledSensitivityError::Interrupted)
                ));
            });
            let repeated = coupled
                .pullback_iqn(cx, &o, &CARRY, interface(), IqnIlsConfig::default())
                .unwrap();
            assert_eq!(baseline.solid, repeated.solid);
            assert_eq!(baseline.inlets, repeated.inlets);
            assert_eq!(baseline.log_htc, repeated.log_htc);
            assert_eq!(baseline.interface_adjoint, repeated.interface_adjoint);
        })
    });
}
