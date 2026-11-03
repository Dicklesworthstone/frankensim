//! Original source melting reaches spatial fields, heat audits and adjoints.
mod support;
#[path = "material_melting/heterogeneous.rs"]
mod heterogeneous;

use fs_conduction::lumped::{
    BiotGate, LumpedEnthalpyBody, LumpedEnthalpyMarchConfig, solve_lumped_enthalpy,
};
use fs_conduction::material::LiquidMassFractionConductivity;
use fs_conduction::transient::enthalpy::adjoint::EnthalpyAdjointError;
use fs_conduction::transient::enthalpy::material::melting::{
    CoexistenceConductivity, MaterialMeltingSource, MaterialMeltingTransport,
};
use fs_conduction::transient::enthalpy::{EnthalpyBudget, EnthalpyError, EnthalpyStepConfig};
use fs_conduction::{
    ConductionError, ConductionMesh, ConductivityModel, LinearConfig, ScalarField, ThermalBc,
    ThermalBoundaryBuilder,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::phase::melting::MeltingSourceQuery;
use fs_material::phase::{HeatCapacitySourceQuery, SolidLiquidPhase};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QuantitySpec, Temperature};
use fs_rep_mesh::TetComplex;
use fs_solver::NewtonKrylovConfig;
use support::{with_cancelled_cx, with_cx};

const CP: Dims = Dims([2, 0, -2, -1, 0, 0]);
const RHO: Dims = Dims([-3, 1, 0, 0, 0, 0]);
const K: Dims = Dims([1, 1, -3, -1, 0, 0]);
const ENERGY: Dims = Dims([2, 0, -2, 0, 0, 0]);
const PRESSURE: Dims = Dims([-1, 1, -2, 0, 0, 0]);
const MASS: f64 = 2.0 / 24.0;
const OLD: [f64; 4] = [60.0, 180.0, 240.0, 420.0];
const Q: [f64; 4] = [8.0, -3.0, 12.0, 5.0];
const DT: f64 = 0.03;

fn absolute_temperature() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}

fn point(t: f64, phase: Option<f64>, pressure: f64) -> QueryPoint {
    let point = QueryPoint::new()
        .with_quantity("T", absolute_temperature(), t)
        .unwrap()
        .with_quantity("p", QuantitySpec::dimensional(PRESSURE), pressure)
        .unwrap();
    phase.map_or_else(
        || point.clone(),
        |f| {
            point
                .clone()
                .with_quantity("phase", QuantitySpec::dimensional(Dims::NONE), f)
                .unwrap()
        },
    )
}

fn card(liquid_conductivity: bool, conductivity_kink: bool) -> MaterialCard {
    let mut claims = ClaimSet::new();
    let pressure = || {
        ValidityDomain::unconstrained().with_quantity(
            "p",
            QuantitySpec::dimensional(PRESSURE),
            100_000.0,
            200_000.0,
        )
    };
    for (phase, lo, hi, cp, density, ks) in [
        (0.0, 300.0, 350.0, 2.0, 2.0, [2.0, 4.0]),
        (1.0, 350.0, 450.0, 4.0, 1.0, [8.0, 12.0]),
    ] {
        for (name, dims, values) in [
            ("cp", CP, [cp, cp]),
            ("rho", RHO, [density, density]),
            ("k", K, ks),
        ] {
            if name == "k" && phase == 1.0 && !liquid_conductivity {
                continue;
            }
            let knots = if name == "k" && phase == 0.0 && conductivity_kink {
                vec![(lo, values[0]), (325.0, 3.5), (hi, values[1])]
            } else {
                vec![(lo, values[0]), (hi, values[1])]
            };
            claims
                .insert_claim(PropertyClaim {
                    key: PropertyKey::new(name, dims),
                    value: PropertyValue::Curve {
                        abscissa: "T".into(),
                        abscissa_dims: Temperature::DIMS,
                        knots,
                        dims,
                    },
                    validity: pressure()
                        .with_quantity("T", absolute_temperature(), lo, hi)
                        .with_quantity(
                            "phase",
                            QuantitySpec::dimensional(Dims::NONE),
                            phase,
                            phase,
                        ),
                    uncertainty: UncertaintyModel::Unstated,
                    interpolation: InterpolationPolicy::LinearInside,
                    observations: vec![],
                    provenance: Provenance {
                        source: "synthetic analytic melting fixture; not measured material data"
                            .into(),
                        license: "test-only".into(),
                        artifact: None,
                    },
                })
                .unwrap();
        }
    }
    for (name, quantity, value) in [
        ("Tm", absolute_temperature(), 350.0),
        ("fusion", QuantitySpec::dimensional(ENERGY), 200.0),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(name, quantity),
                value: PropertyValue::Scalar {
                    value,
                    dims: quantity.dims(),
                },
                validity: pressure(),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic sharp equilibrium transition".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic two-phase thermal reference".into(),
            phase: "solid/liquid".into(),
            process: "analytic fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn resolve(
    card: &MaterialCard,
    coexistence: CoexistenceConductivity,
    reference: &QueryPoint,
) -> Result<
    MaterialMeltingTransport,
    fs_conduction::transient::enthalpy::material::MaterialEnthalpyError,
> {
    let cp = PropertyKey::new("cp", CP);
    let rho = PropertyKey::new("rho", RHO);
    let k = PropertyKey::new("k", K);
    let points = [
        point(300.0, Some(0.0), 100_000.0),
        point(350.0, Some(0.0), 100_000.0),
        point(350.0, Some(1.0), 100_000.0),
        point(450.0, Some(1.0), 100_000.0),
        point(350.0, None, 100_000.0),
    ];
    let branch = |index| HeatCapacitySourceQuery {
        heat_capacity: &cp,
        density: &rho,
        temperature_axis: "T",
        lower: &points[index],
        upper: &points[index + 1],
        selection: SelectionPolicy::SingleClaimOnly,
    };
    MaterialMeltingTransport::from_isotropic_card(
        card,
        MaterialMeltingSource {
            caloric: MeltingSourceQuery {
                solid: branch(0),
                liquid: branch(2),
                melting_point: &PropertyKey::with_quantity("Tm", absolute_temperature()),
                latent_heat_fusion: &PropertyKey::new("fusion", ENERGY),
                transition: &points[4],
                transition_selection: SelectionPolicy::SingleClaimOnly,
                phase_selector: Some("phase"),
            },
            conductivity: [&k, &k],
            reference_phase: SolidLiquidPhase::Solid,
            reference,
            lower_specific_enthalpy_j_kg: 0.0,
            maximum_interpolation_error_j_kg: 1e-6,
            coexistence,
        },
    )
}

fn model(coexistence: CoexistenceConductivity) -> MaterialMeltingTransport {
    resolve(
        &card(true, false),
        coexistence,
        &point(300.0, Some(0.0), 100_000.0),
    )
    .unwrap()
}

fn mesh() -> ConductionMesh {
    ConductionMesh::new(
        TetComplex::from_tets(4, vec![[0, 1, 2, 3]]),
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ],
    )
    .unwrap()
}
fn budget() -> EnthalpyBudget {
    EnthalpyBudget {
        max_vertices: 4,
        max_elements: 1,
    }
}
fn config() -> EnthalpyStepConfig {
    EnthalpyStepConfig {
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-12,
            relative_tolerance: 1e-12,
            linear_restart: 4,
            max_linear_cycles: 8,
            forcing_minimum: 1e-12,
            forcing_maximum: 1e-4,
            ..NewtonKrylovConfig::default()
        },
        max_newton_iterations: 32,
        energy_tolerance_j: 1e-10,
    }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}

// Independent constitutive arithmetic, with h_s=100, h_l=300, rho_s/rho_l=2.
fn reference(h: f64, law: CoexistenceConductivity) -> (f64, f64, f64) {
    if h < 100.0 {
        return (300.0 + h / 2.0, 0.0, 2.0 + h / 50.0);
    }
    if h >= 300.0 {
        return (350.0 + (h - 300.0) / 4.0, 1.0, 8.0 + (h - 300.0) / 100.0);
    }
    let f = (h - 100.0) / 200.0;
    let k = match law {
        CoexistenceConductivity::Parallel => 4.0 + 8.0 * f / (1.0 + f),
        CoexistenceConductivity::Series => 4.0 * (1.0 + f),
    };
    (350.0, f, k)
}

fn residual(h: &[f64], old: &[f64], q: &[f64], law: CoexistenceConductivity) -> [f64; 4] {
    let t = h.iter().map(|&h| reference(h, law).0).collect::<Vec<_>>();
    let k = h.iter().map(|&h| reference(h, law).2 / 4.0).sum::<f64>();
    let laplace = [
        3.0 * t[0] - t[1] - t[2] - t[3],
        t[1] - t[0],
        t[2] - t[0],
        t[3] - t[0],
    ];
    // Right tetrahedron: V=1/6, grad N0=(-1,-1,-1), grad Ni=e_i.
    // Consistent source load integral(N_i*N_j)=(1+delta_ij)/120.
    std::array::from_fn(|i| {
        MASS * (h[i] - old[i])
            + DT * (k * laplace[i] / 6.0 - (q.iter().sum::<f64>() + q[i]) / 120.0)
    })
}

#[test]
fn original_phase_sources_volume_closure_and_latent_tangents_are_retained() {
    let parallel = model(CoexistenceConductivity::Parallel);
    let series = model(CoexistenceConductivity::Series);
    assert_ne!(parallel.identity(), series.identity());
    near(parallel.reference_density_kg_m3(), 2.0, 0.0);
    assert_eq!(
        parallel.caloric().plateau_enthalpy_bounds_j_kg(),
        [100.0, 300.0]
    );
    for transport in [&parallel, &series] {
        for (phase, table) in transport.conductivity_sources().iter().enumerate() {
            assert!(table.receipts().iter().all(|r| {
                r.query_point
                    .iter()
                    .any(|(axis, value)| axis == "phase" && *value == phase as f64)
            }));
        }
        for h in [40.0, 100.0, 160.0, 260.0, 300.0, 480.0] {
            let got = transport.conductivity_at_specific_enthalpy(h).unwrap();
            near(
                got.conductivity_w_m_k,
                reference(h, transport.coexistence()).2,
                2e-14,
            );
        }
        for h in [40.0, 160.0, 260.0, 480.0] {
            let got = transport.conductivity_at_specific_enthalpy(h).unwrap();
            let epsilon = 0.001;
            let fd = (reference(h + epsilon, transport.coexistence()).2
                - reference(h - epsilon, transport.coexistence()).2)
                / (2.0 * epsilon);
            near(got.derivative_per_j_kg, fd, 2e-11);
        }
        let half = transport.conductivity_at_specific_enthalpy(200.0).unwrap();
        near(half.liquid_volume_fraction, 2.0 / 3.0, 1e-15);
        assert!(transport.conductivity_at_specific_enthalpy(-1.0).is_err());
        assert!(transport.conductivity_at_specific_enthalpy(701.0).is_err());
    }
    assert!(parallel.conductivity_sources()[0].eval(351.0).is_err());
    assert!(parallel.conductivity_sources()[1].eval(349.0).is_err());
}

#[test]
fn mixed_phase_spatial_balance_tangent_and_history_adjoint_match_independent_equations() {
    let mesh = mesh();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let q = ScalarField::nodal("deposition", 4, Q.to_vec()).unwrap();
    for law in [
        CoexistenceConductivity::Parallel,
        CoexistenceConductivity::Series,
    ] {
        let material = model(law);
        with_cx(|cx| {
            let solver = material.prepare_uniform(cx, &mesh, budget()).unwrap();
            let problem = material.problem(&mesh, &boundary, &q);
            let response = solver
                .linearize_step(cx, problem, None, &OLD, DT, config())
                .unwrap();
            let h = &response.primal().specific_enthalpy_j_kg;
            for value in residual(h, &OLD, &Q, law) {
                near(value, 0.0, 2e-10);
            }
            near(
                response.primal().stored_energy_change_j,
                DT * Q.iter().sum::<f64>() / 24.0,
                1e-10,
            );
            for column in 0..4 {
                let mut direction = [0.0; 4];
                direction[column] = 1.0;
                let action = response.apply_jacobian(cx, &direction).unwrap();
                let mut plus = h.clone();
                let mut minus = h.clone();
                plus[column] += 0.01;
                minus[column] -= 0.01;
                let rp = residual(&plus, &OLD, &Q, law);
                let rm = residual(&minus, &OLD, &Q, law);
                for row in 0..4 {
                    near(action[row], (rp[row] - rm[row]) / 0.02, 2e-9);
                }
                if column == 1 || column == 2 {
                    assert!(action[0].abs() > 1e-4);
                }
            }
            let weights = [0.7, -0.4, 0.2, 1.1];
            let seed = response.temperature_pullback(cx, &weights).unwrap();
            assert_eq!([seed[1], seed[2]], [0.0, 0.0]);
            let gradient = response
                .pullback(
                    cx,
                    &seed,
                    LinearConfig {
                        tolerance: 1e-12,
                        max_iterations: 64,
                        restart: 4,
                    },
                )
                .unwrap();
            for column in 0..4 {
                let objective = |delta| {
                    let mut old = OLD;
                    old[column] += delta;
                    let next = solver
                        .advance(cx, problem, None, &old, DT, config())
                        .unwrap();
                    next.temperature
                        .iter()
                        .zip(weights)
                        .map(|(t, w)| (t - 350.0) * w)
                        .sum::<f64>()
                };
                let e = 0.02;
                let fd = (-objective(2.0 * e) + 8.0 * objective(e) - 8.0 * objective(-e)
                    + objective(-2.0 * e))
                    / (12.0 * e);
                near(gradient.previous_specific_enthalpy[column], fd, 1e-8);
            }
        });
    }
}

#[test]
fn uniform_source_melting_matches_original_lumped_owner_without_changing_mass() {
    let mesh = mesh();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let material = model(CoexistenceConductivity::Parallel);
    let source = ScalarField::uniform("uniform deposited heat", 100.0).unwrap();
    let body = LumpedEnthalpyBody::try_new(
        "uniform source-melting body",
        1.0 / 3.0,
        1.5 + 0.5 * 3.0_f64.sqrt(),
        0.0,
        0.0,
        0.1,
        2.0,
        material.caloric().curve(),
    )
    .unwrap();
    with_cx(|cx| {
        let lumped = solve_lumped_enthalpy(
            cx,
            &body,
            BiotGate::corpus_default(),
            LumpedEnthalpyMarchConfig {
                initial_specific_enthalpy_j_kg: 50.0,
                ambient_temperature_k: 300.0,
                radiation_temperature_k: 300.0,
                internal_power_w: 100.0 / 6.0,
                duration_s: 6.0,
                maximum_step_s: 1.0,
                maximum_steps: 6,
                enthalpy_tolerance_j_kg: 1e-9,
            },
        )
        .unwrap();
        let solver = material.prepare_uniform(cx, &mesh, budget()).unwrap();
        let masses = solver.reference_nodal_masses_kg().to_vec();
        let mut h = vec![50.0; 4];
        for step in 1..=6 {
            let next = solver
                .advance(
                    cx,
                    material.problem(&mesh, &boundary, &source),
                    None,
                    &h,
                    1.0,
                    config(),
                )
                .unwrap();
            for vertex in 0..4 {
                near(
                    next.specific_enthalpy_j_kg[vertex],
                    50.0 + 50.0 * step as f64,
                    1e-8,
                );
                near(
                    next.specific_enthalpy_j_kg[vertex],
                    lumped.samples()[step].phase_state.specific_enthalpy_j_kg(),
                    1e-8,
                );
                near(
                    next.temperature[vertex],
                    reference(50.0 + 50.0 * step as f64, material.coexistence()).0,
                    1e-8,
                );
            }
            near(next.stored_energy_change_j, 100.0 / 6.0, 1e-9);
            h = next.specific_enthalpy_j_kg;
        }
        assert_eq!(solver.reference_nodal_masses_kg(), masses);
        // The liquid source density is 1 kg/m3, but storage retains rho0=2.
        near(
            material
                .caloric()
                .curve()
                .state_at_specific_enthalpy(h[0])
                .unwrap()
                .bulk_density_kg_m3(),
            1.0,
            0.0,
        );
        let cool = ScalarField::uniform("removed heat", -100.0).unwrap();
        for _ in 0..6 {
            h = solver
                .advance(
                    cx,
                    material.problem(&mesh, &boundary, &cool),
                    None,
                    &h,
                    1.0,
                    config(),
                )
                .unwrap()
                .specific_enthalpy_j_kg;
        }
        for value in h {
            near(value, 50.0, 1e-8);
        }
    });
}

#[test]
fn boundary_driven_spatial_latent_front_refines_in_time_and_closes_heat() {
    let mesh = mesh();
    let material = model(CoexistenceConductivity::Series);
    let source = ScalarField::uniform("no volumetric heat", 0.0).unwrap();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "400 K bath",
            |face| face.centroid[0] == 0.0,
            ThermalBc::dirichlet(400.0).unwrap(),
        )
        .unwrap()
        .adiabatic_remainder()
        .finish()
        .unwrap();
    with_cx(|cx| {
        let solver = material.prepare_uniform(cx, &mesh, budget()).unwrap();
        let mut errors = Vec::new();
        for steps in [5, 10] {
            let dt = 0.05 / steps as f64;
            let mut h = vec![500.0, 150.0, 500.0, 500.0];
            let mut heat = 0.0;
            let mut front = 2.0 / 3.0;
            for step in 1..=steps {
                let next = solver
                    .advance(
                        cx,
                        material.problem(&mesh, &boundary, &source),
                        None,
                        &h,
                        dt,
                        config(),
                    )
                    .unwrap();
                // Independent P1 balance: Tfree=350, kfree=4*(1+f),
                // kmean=(3*10+kfree)/4. Therefore dh/dt=800+0.5*h.
                let discrete_h = 1750.0 / (1.0 - 0.5 * dt).powi(step) - 1600.0;
                near(next.specific_enthalpy_j_kg[1], discrete_h, 1e-8);
                near(next.temperature[1], 350.0, 0.0);
                let computed_front = 0.5 / (1.0 - next.liquid_mass_fraction[1]);
                assert!(computed_front > front && computed_front < 1.0);
                front = computed_front;
                heat += dt * next.dirichlet_in_w;
                h = next.specific_enthalpy_j_kg;
            }
            near(heat, MASS * (h[1] - 150.0), 1e-9);
            // The continuum-in-time semi-discrete solution is independent of
            // the production stepper; no future phase/front is imposed.
            let exact_h = 1750.0 * 0.025_f64.exp() - 1600.0;
            let exact_front = 0.5 / (1.0 - (exact_h - 100.0) / 200.0);
            errors.push((front - exact_front).abs());
        }
        assert!(errors[1] < 0.55 * errors[0]);
    });
}

#[test]
fn source_validity_bound_conductivity_cancellation_and_budget_refuse_atomically() {
    let reference = point(300.0, Some(0.0), 100_000.0);
    assert!(
        resolve(
            &card(false, false),
            CoexistenceConductivity::Parallel,
            &reference
        )
        .is_err()
    );
    assert!(
        resolve(
            &card(true, false),
            CoexistenceConductivity::Parallel,
            &point(360.0, Some(0.0), 100_000.0)
        )
        .is_err()
    );
    assert!(
        resolve(
            &card(true, false),
            CoexistenceConductivity::Parallel,
            &point(300.0, Some(0.0), 150_000.0)
        )
        .is_err()
    );
    let mesh = mesh();
    let material = model(CoexistenceConductivity::Parallel);
    let source = ScalarField::uniform("no source", 0.0).unwrap();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap();
    with_cx(|cx| {
        assert!(matches!(
            material.prepare_uniform(
                cx,
                &mesh,
                EnthalpyBudget {
                    max_vertices: 3,
                    max_elements: 1
                }
            ),
            Err(EnthalpyError::Budget { .. })
        ));
        let solver = material.prepare_uniform(cx, &mesh, budget()).unwrap();
        let mut problem = material.problem(&mesh, &boundary, &source);
        let other = ConductivityModel::isotropic_declared(1.0).unwrap();
        problem.material = &other;
        let saved = OLD.map(f64::to_bits);
        assert!(matches!(
            solver.advance(cx, problem, None, &OLD, DT, config()),
            Err(EnthalpyError::InvalidInput(_))
        ));
        with_cancelled_cx(|cancelled| {
            assert!(matches!(
                solver.advance(
                    cancelled,
                    material.problem(&mesh, &boundary, &source),
                    None,
                    &OLD,
                    DT,
                    config()
                ),
                Err(EnthalpyError::Conduction(ConductionError::Cancelled { .. }))
            ));
        });
        assert_eq!(OLD.map(f64::to_bits), saved);
        assert!(
            solver
                .with_phase_conductivity(
                    cx,
                    vec![LiquidMassFractionConductivity::declared(1.0, 2.0).unwrap()]
                )
                .is_err()
        );
    });
}

#[test]
fn extra_source_conductivity_corner_refuses_classical_adjoint() {
    let material = resolve(
        &card(true, true),
        CoexistenceConductivity::Parallel,
        &point(300.0, Some(0.0), 100_000.0),
    )
    .unwrap();
    let mesh = mesh();
    let source = ScalarField::uniform("equilibrium", 0.0).unwrap();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap();
    with_cx(|cx| {
        let solver = material.prepare_uniform(cx, &mesh, budget()).unwrap();
        // h=50 gives T=325: smooth caloric branch, but a real source k(T) corner.
        let result = solver.linearize_step(
            cx,
            material.problem(&mesh, &boundary, &source),
            None,
            &[50.0; 4],
            DT,
            config(),
        );
        assert!(matches!(result, Err(EnthalpyAdjointError::MaterialKink)));
    });
}
