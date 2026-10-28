//! Real spatial enthalpy solves against independent closed-enclosure balances.
mod support;
#[path = "enthalpy_enclosure/vacuum.rs"]
mod vacuum;
#[path = "enthalpy_enclosure/endpoint.rs"]
mod endpoint;
#[path = "enthalpy_enclosure/capacity.rs"]
mod capacity;
#[path = "enthalpy_enclosure/pointwise.rs"]
mod pointwise;

use fs_blake3::ContentHash;
use fs_conduction::radiation::{GrayDiffuseEnclosure, RadiationSurface, ViewFactorMatrix};
use fs_conduction::transient::enthalpy::{
    EnthalpyBackwardEuler, EnthalpyBudget, EnthalpyStepConfig,
};
use fs_conduction::{
    AmbientRadiationConfig, ConductionMesh, ConductionProblem, ConductivityModel,
    ConductivityTable, EMISSIVITY_DIMS, LinearConfig, STEFAN_BOLTZMANN_W_M2_K4,
    SURFACE_EMISSIVITY_PROPERTY, ScalarField, SurfaceEmissivity, ThermalBc, ThermalBoundary,
    ThermalBoundaryBuilder, material::LiquidMassFractionConductivity,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, SelectionPolicy, UncertaintyModel,
};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};
use fs_rep_mesh::TetComplex;
use fs_solver::NewtonKrylovConfig;
use support::{with_cancelled_cx, with_cx};

const RHO: f64 = 10.0;
const CP: f64 = 10.0;
const HTC: f64 = 3.0;
const DT: f64 = 0.05;
const EPS: [f64; 2] = [0.7, 0.9];

fn emissivity(name: &str, value: f64) -> SurfaceEmissivity {
    emissivity_in(name, value, 240.0, 560.0)
}

fn emissivity_in(name: &str, value: f64, lower: f64, upper: f64) -> SurfaceEmissivity {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS),
            value: PropertyValue::Scalar {
                value,
                dims: EMISSIVITY_DIMS,
            },
            validity: ValidityDomain::unconstrained().with("T", lower, upper),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: Vec::new(),
            provenance: Provenance {
                source: "synthetic two-patch numerical reference".into(),
                license: "internal-test-use".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic gray".into(),
            phase: "solid".into(),
            process: "numerical reference".into(),
            revision: 0,
        },
        claims,
        Vec::new(),
    )
    .unwrap();
    SurfaceEmissivity::from_card(name, &card, 350.0, SelectionPolicy::SingleClaimOnly).unwrap()
}

struct Fixture {
    mesh: ConductionMesh,
    boundary: ThermalBoundary,
    enclosure: GrayDiffuseEnclosure,
    curve: EquilibriumEnthalpyPhaseCurve,
    material: ConductivityModel,
}

impl Fixture {
    fn new(latent: bool) -> Self {
        let tetra = [
            [0.0; 3],
            [1.0, 0.0, 0.0],
            [0.5, 3.0_f64.sqrt() / 2.0, 0.0],
            [0.5, 3.0_f64.sqrt() / 6.0, (2.0_f64 / 3.0).sqrt()],
        ];
        let mut positions = tetra.to_vec();
        positions.extend(tetra.map(|[x, y, z]| [x + 3.0, y, z]));
        let mesh = ConductionMesh::new(
            TetComplex::from_tets(8, vec![[0, 1, 2, 3], [4, 5, 6, 7]]),
            positions,
        )
        .unwrap();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .region(
                "left",
                |face| face.vertices.iter().all(|v| *v < 4),
                ThermalBc::robin(HTC, 300.0).unwrap(),
            )
            .unwrap()
            .region(
                "right",
                |face| face.vertices.iter().all(|v| *v >= 4),
                ThermalBc::robin(HTC, 300.0).unwrap(),
            )
            .unwrap()
            .finish()
            .unwrap();
        let left = RadiationSurface::new(
            &mesh,
            "left",
            |face| face.vertices.iter().all(|v| *v < 4),
            emissivity("left", EPS[0]),
        )
        .unwrap();
        let right = RadiationSurface::new(
            &mesh,
            "right",
            |face| face.vertices.iter().all(|v| *v >= 4),
            emissivity("right", EPS[1]),
        )
        .unwrap();
        // Deliberately prescribed two-patch enclosure, independent of the
        // bodies' displayed separation. No mesh-derived view-factor claim.
        let factors = ViewFactorMatrix::infinite_parallel_plates(left.area_m2()).unwrap();
        let enclosure = GrayDiffuseEnclosure::new(vec![left, right], factors).unwrap();
        let data = if latent {
            vec![
                (0.0, 250.0, 0.0),
                (1000.0, 350.0, 0.0),
                (3000.0, 350.0, 1.0),
                (4500.0, 500.0, 1.0),
            ]
        } else {
            vec![(0.0, 200.0, 0.0), (4000.0, 600.0, 0.0)]
        };
        let knots = data
            .into_iter()
            .map(|(h, t, f)| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: h,
                temperature_k: t,
                liquid_mass_fraction: f,
                bulk_density_kg_m3: RHO,
            })
            .collect();
        let curve = if latent {
            EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([0x78; 32]), knots)
        } else {
            EquilibriumEnthalpyPhaseCurve::try_single_phase(
                ContentHash([0x79; 32]),
                SolidLiquidPhase::Solid,
                knots,
            )
        }
        .unwrap();
        let material = ConductivityModel::isotropic(
            ConductivityTable::declared_curve(vec![(240.0, 2.0), (560.0, 8.0)]).unwrap(),
        );
        Self {
            mesh,
            boundary,
            enclosure,
            curve,
            material,
        }
    }

    fn problem<'a>(&'a self, source: &'a ScalarField) -> ConductionProblem<'a> {
        ConductionProblem {
            mesh: &self.mesh,
            boundary: &self.boundary,
            material: &self.material,
            element_materials: None,
            source,
        }
    }

    fn engine(&self, cx: &fs_exec::Cx<'_>, phase_k: bool) -> EnthalpyBackwardEuler<'_, '_> {
        let engine = EnthalpyBackwardEuler::uniform(
            cx,
            &self.mesh,
            &self.curve,
            RHO,
            EnthalpyBudget {
                max_vertices: 8,
                max_elements: 2,
            },
        )
        .unwrap();
        if phase_k {
            engine
                .with_phase_conductivity(
                    cx,
                    vec![
                        LiquidMassFractionConductivity::declared(0.7, 2.2).unwrap(),
                        LiquidMassFractionConductivity::declared(1.1, 0.6).unwrap(),
                    ],
                )
                .unwrap()
        } else {
            engine
        }
    }
}

fn step_config() -> EnthalpyStepConfig {
    EnthalpyStepConfig {
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-10,
            relative_tolerance: 1e-12,
            linear_restart: 8,
            max_linear_cycles: 8,
            forcing_minimum: 1e-12,
            forcing_maximum: 1e-3,
            ..NewtonKrylovConfig::default()
        },
        max_newton_iterations: 32,
        energy_tolerance_j: 1e-8,
    }
}

fn radiation_config() -> AmbientRadiationConfig {
    AmbientRadiationConfig {
        max_iterations: 160,
        temperature_tolerance_k: 1e-11,
        balance_tolerance_w: 1e-8,
        balance_relative_tolerance: 0.0,
        relaxation: 0.5,
    }
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "actual {actual:.14e}, independent expected {expected:.14e}"
    );
}

#[test]
fn two_spatial_bodies_match_independent_reflected_exchange_and_energy() {
    let fixture = Fixture::new(false);
    let source = ScalarField::Uniform(0.0);
    with_cx(|cx| {
        let engine = fixture.engine(cx, false);
        for temperatures in [[400.0, 300.0], [350.0, 350.0]] {
            let old: Vec<_> = temperatures
                .into_iter()
                .flat_map(|t| [CP * (t - 200.0); 4])
                .collect();
            let history = old.clone();
            let solve = || {
                engine
                    .advance_with_enclosure_radiation(
                        cx,
                        fixture.problem(&source),
                        None,
                        &old,
                        DT,
                        step_config(),
                        &fixture.enclosure,
                        radiation_config(),
                    )
                    .unwrap()
            };
            let got = solve();
            assert_eq!(got.radiation.radiosity, fixture.enclosure.solve(
                cx, &got.radiation.radiosity.surface_temperatures_k,
            ).unwrap());
            let capacity = RHO * CP * 2.0_f64.sqrt() / 12.0;
            let area = 3.0_f64.sqrt();
            let sum = (capacity * (temperatures[0] + temperatures[1]) + DT * HTC * area * 600.0)
                / (capacity + DT * HTC * area);
            let factor = STEFAN_BOLTZMANN_W_M2_K4 / (1.0 / EPS[0] + 1.0 / EPS[1] - 1.0);
            let (mut lo, mut hi) = (0.0, temperatures[0] - temperatures[1]);
            for _ in 0..80 {
                let delta = (lo + hi) / 2.0;
                let left = (sum + delta) / 2.0;
                let right = (sum - delta) / 2.0;
                let residual = capacity * (delta - (temperatures[0] - temperatures[1]))
                    + DT * HTC * area * delta
                    + 2.0 * DT * area * factor * (left.powi(4) - right.powi(4));
                if residual > 0.0 {
                    hi = delta;
                } else {
                    lo = delta;
                }
            }
            let expected = [(sum + (lo + hi) / 2.0) / 2.0, (sum - (lo + hi) / 2.0) / 2.0];
            for (body, &temperature) in expected.iter().enumerate() {
                for node in 0..4 {
                    close(
                        got.conduction.temperature[4 * body + node],
                        temperature,
                        2e-10,
                    );
                }
            }
            let heat = area * factor * (expected[0].powi(4) - expected[1].powi(4));
            close(got.radiation.radiosity.net_outward_heat_w[0], heat, 2e-9);
            close(got.radiation.radiosity.net_outward_heat_w[1], -heat, 2e-9);
            close(
                got.conduction.stored_energy_change_j,
                -DT * got.convective_out_w,
                1e-9,
            );
            assert!(got.physical_residual_norm_j <= got.physical_residual_tolerance_j);
            assert!(got.physical_energy_residual_j.abs() <= step_config().energy_tolerance_j);
            assert_eq!(
                got.conduction.specific_enthalpy_j_kg,
                solve().conduction.specific_enthalpy_j_kg
            );
            assert_eq!(old, history);
        }
    });
}

#[test]
fn latent_patch_transfers_heat_to_sensible_neighbor_without_fictitious_storage() {
    let fixture = Fixture::new(true);
    let source = ScalarField::Uniform(0.0);
    let old = [2000.0, 2000.0, 2000.0, 2000.0, 500.0, 500.0, 500.0, 500.0];
    with_cx(|cx| {
        let got = fixture
            .engine(cx, true)
            .advance_with_enclosure_radiation(
                cx,
                fixture.problem(&source),
                None,
                &old,
                DT,
                step_config(),
                &fixture.enclosure,
                radiation_config(),
            )
            .unwrap();
        for node in 0..4 {
            assert_eq!(got.conduction.temperature[node], 350.0);
            assert!(got.conduction.specific_enthalpy_j_kg[node] < old[node]);
            assert!(got.conduction.liquid_mass_fraction[node] < 0.5);
        }
        assert!(got.conduction.temperature[4] > 300.0);
        let left_heat = got.radiation.radiosity.net_outward_heat_w[0];
        let mass = RHO * 2.0_f64.sqrt() / 12.0;
        let expected = 2000.0 - DT * (HTC * 3.0_f64.sqrt() * 50.0 + left_heat) / mass;
        close(got.conduction.specific_enthalpy_j_kg[0], expected, 1e-10);
        close(
            got.conduction.stored_energy_change_j,
            -DT * got.convective_out_w,
            1e-9,
        );
    });
}

#[test]
fn complete_h_tangent_and_history_source_pullbacks_match_perturbed_physics() {
    let fixture = Fixture::new(true);
    let old = [600.0, 1800.0, 3400.0, 3700.0, 700.0, 1900.0, 3300.0, 3600.0];
    let q = [12.0, -3.0, 8.0, 5.0, 6.0, 9.0, -2.0, 4.0];
    let goal = [0.7, -0.4, 0.2, 1.1, -0.3, 0.6, 0.9, -0.1];
    let direction = [0.2, -0.5, 0.7, 0.3, -0.9, 0.6, -0.1, 0.8];
    with_cx(|cx| {
        let engine = fixture.engine(cx, true);
        let source = ScalarField::nodal("source", 8, q.to_vec()).unwrap();
        let accepted = engine
            .advance_with_enclosure_radiation(
                cx,
                fixture.problem(&source),
                None,
                &old,
                DT,
                step_config(),
                &fixture.enclosure,
                radiation_config(),
            )
            .unwrap();
        let linear = engine
            .linearize_accepted_with_enclosure_radiation(
                cx,
                fixture.problem(&source),
                None,
                &old,
                DT,
                step_config(),
                &fixture.enclosure,
                radiation_config(),
                accepted.conduction,
                32,
            )
            .unwrap();
        let jv = linear.apply_jacobian(cx, &direction).unwrap();
        let jtu = linear.apply_jacobian_transpose(cx, &goal).unwrap();
        close(dot(&goal, &jv), dot(&jtu, &direction), 1e-12);
        // Standalone audits use fresh one-shot radiosity, independently of the
        // factors retained by the complete reflected-radiation tangent.
        let delta = 0.03;
        let residuals = [-delta, delta].map(|offset| {
            let candidate = linear.primal().specific_enthalpy_j_kg.iter().zip(direction)
                .map(|(&h, v)| h + offset * v).collect::<Vec<_>>();
            engine.evaluate_endpoint_with_enclosure(cx, fixture.problem(&source), None,
                &old, &candidate, DT, step_config(), &fixture.enclosure, radiation_config())
                .unwrap().residual_norm_j
        });
        close((residuals[0] + residuals[1]) / (2.0 * delta), dot(&jv, &jv).sqrt(), 1e-7);
        let seed = linear.temperature_pullback(cx, &goal).unwrap();
        let gradient = linear
            .pullback(
                cx,
                &seed,
                LinearConfig {
                    tolerance: 1e-11,
                    max_iterations: 128,
                    restart: 8,
                },
            )
            .unwrap();
        let objective = |history: &[f64], power: &[f64]| {
            let source = ScalarField::nodal("source", 8, power.to_vec()).unwrap();
            let value = engine
                .advance_with_enclosure_radiation(
                    cx,
                    fixture.problem(&source),
                    None,
                    history,
                    DT,
                    step_config(),
                    &fixture.enclosure,
                    radiation_config(),
                )
                .unwrap();
            dot(&value.conduction.temperature, &goal)
        };
        let perturb = |base: &[f64], delta: f64| {
            base.iter()
                .zip(direction)
                .map(|(a, b)| a + delta * b)
                .collect::<Vec<_>>()
        };
        let old_fd =
            (objective(&perturb(&old, 0.03), &q) - objective(&perturb(&old, -0.03), &q)) / 0.06;
        close(
            dot(&gradient.previous_specific_enthalpy, &direction),
            old_fd,
            3e-7,
        );
        let source_fd =
            (objective(&old, &perturb(&q, 0.2)) - objective(&old, &perturb(&q, -0.2))) / 0.4;
        close(dot(&gradient.source_density, &direction), source_fd, 3e-8);
        // Original convection must survive binding; shifted references are
        // algebraic radiative loads and cannot masquerade as physical inlets.
        let response = linear
            .robin_response(
                cx,
                &["left", "right"],
                LinearConfig {
                    tolerance: 1e-11,
                    max_iterations: 128,
                    restart: 8,
                },
            )
            .unwrap();
        assert!(
            response
                .ports()
                .iter()
                .all(|port| port.reference_k == 300.0 && port.htc_w_m2_k == HTC)
        );
    });
}

#[test]
fn enclosure_refuses_underbudget_forged_endpoint_and_cancelled_trials() {
    let fixture = Fixture::new(true);
    let source = ScalarField::Uniform(0.0);
    let old = [2000.0, 2000.0, 2000.0, 2000.0, 500.0, 500.0, 500.0, 500.0];
    with_cx(|cx| {
        let engine = fixture.engine(cx, false);
        let mut one = radiation_config();
        one.max_iterations = 1;
        assert!(
            engine
                .advance_with_enclosure_radiation(
                    cx,
                    fixture.problem(&source),
                    None,
                    &old,
                    DT,
                    step_config(),
                    &fixture.enclosure,
                    one
                )
                .is_err()
        );
        let accepted = engine
            .advance_with_enclosure_radiation(
                cx,
                fixture.problem(&source),
                None,
                &old,
                DT,
                step_config(),
                &fixture.enclosure,
                radiation_config(),
            )
            .unwrap()
            .conduction;
        assert!(
            engine
                .linearize_accepted_with_enclosure_radiation(
                    cx,
                    fixture.problem(&source),
                    None,
                    &old,
                    DT,
                    step_config(),
                    &fixture.enclosure,
                    radiation_config(),
                    accepted.clone(),
                    31
                )
                .is_err()
        );
        let mut forged = accepted;
        forged.specific_enthalpy_j_kg[0] += 10.0;
        forged.newton.residual_norm = 0.0;
        forged.energy_residual_j = 0.0;
        assert!(
            engine
                .linearize_accepted_with_enclosure_radiation(
                    cx,
                    fixture.problem(&source),
                    None,
                    &old,
                    DT,
                    step_config(),
                    &fixture.enclosure,
                    radiation_config(),
                    forged,
                    32
                )
                .is_err()
        );
        with_cancelled_cx(|cancelled| {
            assert!(
                engine
                    .advance_with_enclosure_radiation(
                        cancelled,
                        fixture.problem(&source),
                        None,
                        &old,
                        DT,
                        step_config(),
                        &fixture.enclosure,
                        radiation_config(),
                    )
                    .is_err()
            );
        });
    });
}
