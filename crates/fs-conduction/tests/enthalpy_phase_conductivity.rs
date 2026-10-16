//! G1/G3 phase-conductivity residuals and endpoint adjoints on a right tet.

mod support;

use fs_blake3::ContentHash;
use fs_conduction::material::LiquidMassFractionConductivity;
use fs_conduction::transient::enthalpy::adjoint::EnthalpyAdjointError;
use fs_conduction::transient::enthalpy::heterogeneous::{
    HeterogeneousEnthalpyBackwardEuler, ReferenceEnthalpyMaterial,
};
use fs_conduction::transient::enthalpy::{
    EnthalpyBackwardEuler, EnthalpyBudget, EnthalpyError, EnthalpyStepConfig,
};
use fs_conduction::{
    ConductionError, ConductionMesh, ConductionProblem, ConductivityModel, ConductivityTable,
    LinearConfig, ScalarField, ThermalBoundary, ThermalBoundaryBuilder,
};
use fs_exec::Cx;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
use fs_rep_mesh::TetComplex;
use fs_solver::NewtonKrylovConfig;
use support::{with_cancelled_cx, with_cx};

const OLD: [f64; 4] = [60.0, 180.0, 240.0, 420.0];
const SOURCE: [f64; 4] = [8.0, -3.0, 12.0, 5.0];
const WEIGHTS: [f64; 4] = [0.7, -0.4, 0.2, 1.1];
const DT: f64 = 0.03;
const MASS: f64 = 2.0 / 24.0;
const SOLID_MULTIPLIER: f64 = 0.6;
const LIQUID_MULTIPLIER: f64 = 2.4;

struct Fixture {
    mesh: ConductionMesh,
    boundary: ThermalBoundary,
    curve: EquilibriumEnthalpyPhaseCurve,
    material: ConductivityModel,
}

impl Fixture {
    fn new() -> Self {
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
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let curve = EquilibriumEnthalpyPhaseCurve::try_new(
            ContentHash([0x37; 32]),
            [
                (0.0, 300.0, 0.0, 4.0),
                (100.0, 350.0, 0.0, 3.0),
                (300.0, 350.0, 1.0, 2.0),
                (700.0, 450.0, 1.0, 1.0),
            ]
            .into_iter()
            .map(|(h, t, liquid, density)| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: h,
                temperature_k: t,
                liquid_mass_fraction: liquid,
                bulk_density_kg_m3: density,
            })
            .collect(),
        )
        .unwrap();
        Self {
            mesh,
            boundary,
            curve,
            material: scaled_material(1.0),
        }
    }

    fn stepper(&self, cx: &Cx<'_>) -> EnthalpyBackwardEuler<'_, '_> {
        EnthalpyBackwardEuler::uniform(cx, &self.mesh, &self.curve, 2.0, budget()).unwrap()
    }

    fn phase_stepper(&self, cx: &Cx<'_>) -> EnthalpyBackwardEuler<'_, '_> {
        self.stepper(cx)
            .with_phase_conductivity(cx, vec![phase_law()])
            .unwrap()
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

    fn temperature_objective(&self, cx: &Cx<'_>, old: &[f64], q: &[f64]) -> f64 {
        let source = ScalarField::nodal("nodal deposition", 4, q.to_vec()).unwrap();
        let primal = self
            .phase_stepper(cx)
            .advance(cx, self.problem(&source), None, old, DT, config())
            .unwrap();
        primal
            .temperature
            .iter()
            .zip(WEIGHTS)
            .map(|(t, w)| w * (t - 350.0))
            .sum()
    }
}

fn budget() -> EnthalpyBudget {
    EnthalpyBudget {
        max_vertices: 4,
        max_elements: 1,
    }
}

fn phase_law() -> LiquidMassFractionConductivity {
    LiquidMassFractionConductivity::declared(SOLID_MULTIPLIER, LIQUID_MULTIPLIER).unwrap()
}

fn scaled_material(factor: f64) -> ConductivityModel {
    ConductivityModel::isotropic(
        ConductivityTable::declared_curve(vec![(300.0, 2.0 * factor), (450.0, 5.0 * factor)])
            .unwrap(),
    )
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

fn linear() -> LinearConfig {
    LinearConfig {
        tolerance: 1e-12,
        max_iterations: 64,
        restart: 4,
    }
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}

fn five_point(mut evaluate: impl FnMut(f64) -> f64) -> f64 {
    let epsilon = 0.02;
    (-evaluate(2.0 * epsilon) + 8.0 * evaluate(epsilon) - 8.0 * evaluate(-epsilon)
        + evaluate(-2.0 * epsilon))
        / (12.0 * epsilon)
}

// Independent P1 algebra: V=1/6, grad N0=(-1,-1,-1), grad N1=(1,0,0),
// grad N2=(0,1,0), grad N3=(0,0,1). Source integration is consistent:
// integral(N_i N_j)=V*(1+delta_ij)/20. Storage uses rho_ref*V/4, even
// though the phase chart's equilibrium density varies with enthalpy.
fn reference_residual(h: &[f64], old: &[f64], source: &[f64]) -> [f64; 4] {
    let t: Vec<_> = h
        .iter()
        .map(|&value| {
            if value < 100.0 {
                300.0 + 0.5 * value
            } else if value < 300.0 {
                350.0
            } else {
                350.0 + 0.25 * (value - 300.0)
            }
        })
        .collect();
    let liquid_mean = h
        .iter()
        .map(|h| ((h - 100.0) / 200.0).clamp(0.0, 1.0))
        .sum::<f64>()
        / 4.0;
    let scale = SOLID_MULTIPLIER + (LIQUID_MULTIPLIER - SOLID_MULTIPLIER) * liquid_mean;
    let k = scale * (2.0 + 0.02 * (t.iter().sum::<f64>() / 4.0 - 300.0));
    let transport = [
        3.0 * t[0] - t[1] - t[2] - t[3],
        t[1] - t[0],
        t[2] - t[0],
        t[3] - t[0],
    ];
    std::array::from_fn(|i| {
        let load = (source.iter().sum::<f64>() + source[i]) / 120.0;
        MASS * (h[i] - old[i]) + DT * (k * transport[i] / 6.0 - load)
    })
}

#[test]
fn mixed_phase_endpoint_and_latent_jacobian_columns_match_independent_p1_algebra() {
    let fixture = Fixture::new();
    let source = ScalarField::nodal("nodal deposition", 4, SOURCE.to_vec()).unwrap();
    with_cx(|cx| {
        let stepper = fixture.phase_stepper(cx);
        for &mass in stepper.reference_nodal_masses_kg() {
            assert!((mass - MASS).abs() < 1e-15);
        }
        let response = stepper
            .linearize_step(cx, fixture.problem(&source), None, &OLD, DT, config())
            .unwrap();
        let primal = response.primal();
        let h = &primal.specific_enthalpy_j_kg;
        let slopes: Vec<_> = h
            .iter()
            .map(|&value| {
                fixture
                    .curve
                    .temperature_derivative_at_specific_enthalpy(value)
                    .unwrap()
            })
            .collect();
        assert_eq!(slopes, [0.5, 0.0, 0.0, 0.25]);
        for residual in reference_residual(h, &OLD, &SOURCE) {
            assert!(residual.abs() < 2e-10, "P1 residual: {residual}");
        }
        let stored = MASS * h.iter().zip(OLD).map(|(h, old)| h - old).sum::<f64>();
        let supplied = DT * SOURCE.iter().sum::<f64>() / 24.0;
        assert!((stored - supplied).abs() < 1e-10);
        assert!((primal.stored_energy_change_j - stored).abs() < 1e-12);

        for column in 0..4 {
            let mut basis = [0.0; 4];
            basis[column] = 1.0;
            let action = response.apply_jacobian(cx, &basis).unwrap();
            let mut plus = h.clone();
            let mut minus = h.clone();
            plus[column] += 0.01;
            minus[column] -= 0.01;
            let rp = reference_residual(&plus, &OLD, &SOURCE);
            let rm = reference_residual(&minus, &OLD, &SOURCE);
            for row in 0..4 {
                let fd = (rp[row] - rm[row]) / 0.02;
                assert!(
                    (action[row] - fd).abs() < 2e-9,
                    "J[{row},{column}]: action={}, FD={fd}",
                    action[row]
                );
            }
            if [1, 2].contains(&column) {
                // dT/dh is exactly zero here. Any off-diagonal action must
                // come from k_base(Tbar) * ds/df * df/dh, not a thermal slope.
                assert!(action[0].abs() > 1e-4);
                let transport_norm = action
                    .iter()
                    .enumerate()
                    .map(|(row, value)| (value - MASS * basis[row]).abs())
                    .sum::<f64>();
                assert!(transport_norm > 1e-3);
            }
        }
        let u = [0.7, -1.2, 0.6, 0.3];
        let v = [-0.5, 0.8, 1.4, -0.2];
        let jv = response.apply_jacobian(cx, &v).unwrap();
        let jtu = response.apply_jacobian_transpose(cx, &u).unwrap();
        assert!((dot(&u, &jv) - dot(&jtu, &v)).abs() < 1e-13);
        let ju = response.apply_jacobian(cx, &u).unwrap();
        assert!((dot(&u, &jv) - dot(&ju, &v)).abs() > 1e-4);
    });
}

#[test]
fn latent_history_and_nodal_source_pullbacks_match_perturbed_forward_solves() {
    let fixture = Fixture::new();
    let source = ScalarField::nodal("nodal deposition", 4, SOURCE.to_vec()).unwrap();
    with_cx(|cx| {
        let stepper = fixture.phase_stepper(cx);
        let accepted = stepper
            .advance(cx, fixture.problem(&source), None, &OLD, DT, config())
            .unwrap();
        let response = stepper
            .linearize_accepted(
                cx,
                fixture.problem(&source),
                None,
                &OLD,
                DT,
                config(),
                accepted,
            )
            .unwrap();
        let seed = response.temperature_pullback(cx, &WEIGHTS).unwrap();
        assert_eq!(seed[1], 0.0);
        assert_eq!(seed[2], 0.0);
        let gradient = response.pullback(cx, &seed, linear()).unwrap();
        for vertex in 0..4 {
            let history_fd = five_point(|delta| {
                let mut old = OLD;
                old[vertex] += delta;
                fixture.temperature_objective(cx, &old, &SOURCE)
            });
            let source_fd = five_point(|delta| {
                let mut q = SOURCE;
                q[vertex] += delta;
                fixture.temperature_objective(cx, &OLD, &q)
            });
            assert!(
                (gradient.previous_specific_enthalpy[vertex] - history_fd).abs() < 1e-8,
                "history vertex {vertex}: adjoint={}, FD={history_fd}",
                gradient.previous_specific_enthalpy[vertex]
            );
            assert!(
                (gradient.source_density[vertex] - source_fd).abs() < 1e-8,
                "source vertex {vertex}: adjoint={}, FD={source_fd}",
                gradient.source_density[vertex]
            );
        }
        // A latent node now changes neighboring temperatures by changing
        // conductivity, even though its own temperature seed remains zero.
        for vertex in [1, 2] {
            assert!(gradient.previous_specific_enthalpy[vertex].abs() > 1e-5);
            assert!(gradient.source_density[vertex].abs() > 1e-4);
        }
    });
}

#[test]
fn constant_phase_multipliers_match_the_corresponding_base_conductivity() {
    let fixture = Fixture::new();
    let source = ScalarField::nodal("nodal deposition", 4, SOURCE.to_vec()).unwrap();
    with_cx(|cx| {
        for factor in [1.0, 2.0] {
            let phase = fixture
                .stepper(cx)
                .with_phase_conductivity(
                    cx,
                    vec![LiquidMassFractionConductivity::declared(factor, factor).unwrap()],
                )
                .unwrap();
            let material = scaled_material(factor);
            let mut scaled_problem = fixture.problem(&source);
            scaled_problem.material = &material;
            let actual = phase
                .linearize_step(cx, fixture.problem(&source), None, &OLD, DT, config())
                .unwrap();
            let reference = fixture
                .stepper(cx)
                .linearize_step(cx, scaled_problem, None, &OLD, DT, config())
                .unwrap();
            for (actual, reference) in actual
                .primal()
                .specific_enthalpy_j_kg
                .iter()
                .zip(&reference.primal().specific_enthalpy_j_kg)
            {
                assert!((actual - reference).abs() < 1e-9);
            }
            let direction = [0.7, -0.4, 0.2, 1.1];
            let actual_jv = actual.apply_jacobian(cx, &direction).unwrap();
            let reference_jv = reference.apply_jacobian(cx, &direction).unwrap();
            for (actual, reference) in actual_jv.iter().zip(&reference_jv) {
                assert!((actual - reference).abs() < 1e-11);
            }
        }
    });
}

#[test]
fn heterogeneous_owner_preserves_the_phase_law_and_reference_storage() {
    let fixture = Fixture::new();
    let source = ScalarField::nodal("nodal deposition", 4, SOURCE.to_vec()).unwrap();
    with_cx(|cx| {
        let materials = [ReferenceEnthalpyMaterial {
            curve: &fixture.curve,
            reference_density_kg_m3: 2.0,
        }];
        let heterogeneous =
            HeterogeneousEnthalpyBackwardEuler::new(cx, &fixture.mesh, &materials, &[0], budget())
                .unwrap()
                .with_phase_conductivity(cx, vec![phase_law()])
                .unwrap();
        let uniform = fixture.phase_stepper(cx);
        assert_eq!(
            heterogeneous.reference_nodal_masses_kg(),
            uniform.reference_nodal_masses_kg()
        );
        let actual = heterogeneous
            .advance(cx, fixture.problem(&source), None, &OLD, DT, config())
            .unwrap();
        let reference = uniform
            .advance(cx, fixture.problem(&source), None, &OLD, DT, config())
            .unwrap();
        for (actual, reference) in actual
            .specific_enthalpy_j_kg
            .iter()
            .zip(&reference.specific_enthalpy_j_kg)
        {
            assert!((actual - reference).abs() < 1e-9);
        }
    });
}

#[test]
fn fraction_only_kink_refuses_only_when_conductivity_depends_on_phase() {
    let mut fixture = Fixture::new();
    fixture.curve = EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([0x38; 32]),
        [
            (0.0, 300.0, 0.0),
            (100.0, 350.0, 0.0),
            (200.0, 350.0, 0.2),
            (300.0, 350.0, 1.0),
            (700.0, 450.0, 1.0),
        ]
        .into_iter()
        .map(|(h, t, liquid)| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: h,
            temperature_k: t,
            liquid_mass_fraction: liquid,
            bulk_density_kg_m3: 2.0,
        })
        .collect(),
    )
    .unwrap();
    let source = ScalarField::Uniform(0.0);
    with_cx(|cx| {
        // Uniform plateau temperature and zero load preserve h=200 exactly.
        // T_h is zero on both sides, but f_h changes from 0.002 to 0.008.
        assert!(matches!(
            fixture.phase_stepper(cx).linearize_step(
                cx,
                fixture.problem(&source),
                None,
                &[200.0; 4],
                DT,
                config()
            ),
            Err(EnthalpyAdjointError::ChartKink { .. })
        ));
        for stepper in [
            fixture.stepper(cx),
            fixture
                .stepper(cx)
                .with_phase_conductivity(
                    cx,
                    vec![LiquidMassFractionConductivity::declared(2.0, 2.0).unwrap()],
                )
                .unwrap(),
        ] {
            let response = stepper
                .linearize_step(
                    cx,
                    fixture.problem(&source),
                    None,
                    &[200.0; 4],
                    DT,
                    config(),
                )
                .unwrap();
            assert_eq!(response.primal().specific_enthalpy_j_kg, [200.0; 4]);
        }
    });
}

#[test]
fn invalid_phase_laws_and_assignment_shapes_refuse_and_builder_polls_cancellation() {
    for invalid in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(LiquidMassFractionConductivity::declared(invalid, 1.0).is_err());
        assert!(LiquidMassFractionConductivity::declared(1.0, invalid).is_err());
    }
    let fixture = Fixture::new();
    with_cx(|cx| {
        for laws in [vec![], vec![phase_law(); 2]] {
            assert!(
                fixture
                    .stepper(cx)
                    .with_phase_conductivity(cx, laws)
                    .is_err()
            );
        }
        let stepper = fixture.stepper(cx);
        with_cancelled_cx(|cancelled| {
            assert!(matches!(
                stepper.with_phase_conductivity(cancelled, vec![phase_law()]),
                Err(EnthalpyError::Conduction(ConductionError::Cancelled { .. }))
            ));
        });
    });
}
