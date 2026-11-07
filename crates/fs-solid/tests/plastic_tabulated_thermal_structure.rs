//! A heated member yields under a retained physical dead load, then springs back.
//! Independent free-sided Hencky equations check native tetrahedra and sheets;
//! authored curves are synthetic and do not claim experimental formability.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::hyper::plane_stress::PlaneStressSettings;
use fs_material::plastic::finite::thermal::{
    AffineThermalFiniteJ2, FiniteJ2TemperatureSlopes, FiniteJ2ThermalParameters,
    FiniteJ2ThermalTabulatedHardening,
};
use fs_material::plastic::finite::{FiniteJ2Limits, FiniteJ2Parameters};
use fs_material::plastic::tabulated::{HardeningPoint, TabulatedHardening};
use fs_solid::hyper2d::plastic::SheetAssemblyBudget;
use fs_solid::hyper2d::plastic::thermal::FiniteThermoplasticSheet;
use fs_solid::hyper3::plastic::thermal::{
    FiniteThermoplasticControls, FiniteThermoplasticTetSolid,
};
use fs_solid::hyper3::plastic::{FinitePlasticLawField, FinitePlasticStepControls};
use fs_solid::{HyperTetSettings, IncrementSettings, Mesh2, TetAssemblyBudget};

const NODES: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
const VOLUME: f64 = 1.0 / 6.0;
const THICKNESS: f64 = 0.01;
const E: f64 = 9.0 * 16_000.0 * 8_000.0 / (3.0 * 16_000.0 + 8_000.0);
const NU: f64 = (3.0 * 16_000.0 - 2.0 * 8_000.0) / (2.0 * (3.0 * 16_000.0 + 8_000.0));
const EXPANSION: f64 = 3e-4;
const PATH: [(f64, f64); 5] = [
    (300.0, 100.0),
    (550.0, 100.0),
    (600.0, 100.0),
    (300.0, 100.0),
    (300.0, 0.0),
];

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 87,
                kernel_id: 297,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}
fn family() -> AffineThermalFiniteJ2 {
    let curve = TabulatedHardening::new(
        [(0.0, 200.0), (0.05, 280.0), (0.2, 400.0)]
            .map(|(plastic_strain, yield_stress_pa)| HardeningPoint {
                plastic_strain,
                yield_stress_pa,
            })
            .to_vec(),
    )
    .unwrap();
    AffineThermalFiniteJ2::new(
        FiniteJ2ThermalParameters {
            reference_temperature_k: 300.0,
            temperature_range_k: [300.0, 600.0],
            reference: FiniteJ2Parameters {
                shear_modulus_pa: 8_000.0,
                bulk_modulus_pa: 16_000.0,
                yield_stress_pa: 200.0,
                hardening_modulus_pa: 0.0,
            },
            slopes: FiniteJ2TemperatureSlopes {
                shear_modulus_pa_k: 0.0,
                bulk_modulus_pa_k: 0.0,
                yield_stress_pa_k: -0.4,
                hardening_modulus_pa_k: 0.0,
            },
            logarithmic_expansion_per_k: EXPANSION,
            logarithmic_expansion_slope_per_k2: 0.0,
        },
        FiniteJ2Limits {
            minimum_jacobian: 0.25,
            maximum_elastic_log_strain: 0.7,
            maximum_plastic_log_strain: 2.0,
            maximum_accumulated_plastic_strain: 0.2,
        },
    )
    .unwrap()
    .with_tabulated_hardening(
        FiniteJ2ThermalTabulatedHardening::new(curve, vec![-0.4, -0.56, -0.8]).unwrap(),
    )
    .unwrap()
}
fn controls() -> FinitePlasticStepControls {
    FinitePlasticStepControls {
        equilibrium_work_tolerance_j: 1e-6,
        maximum_endpoint_work_defect_j: 800.0,
        yield_tolerance_pa: 1e-5,
    }
}
fn tet_settings() -> HyperTetSettings {
    HyperTetSettings {
        load_steps: 1,
        relative_tolerance: 1e-11,
        ..HyperTetSettings::default()
    }
}
fn sheet_settings() -> IncrementSettings {
    IncrementSettings {
        relative_tolerance: 1e-11,
        ..IncrementSettings::default()
    }
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual {actual:.16e}, expected {expected:.16e}, tolerance {tolerance:e}"
    );
}
fn root(mut lower: f64, mut upper: f64, f: impl Fn(f64) -> f64) -> f64 {
    assert!(f(lower) <= 0.0 && f(upper) >= 0.0);
    for _ in 0..120 {
        let center = 0.5 * lower + 0.5 * upper;
        if f(center) < 0.0 {
            lower = center;
        } else {
            upper = center;
        }
    }
    0.5 * lower + 0.5 * upper
}
fn reference_yield(alpha: f64) -> f64 {
    if alpha <= 0.05 {
        200.0 + 1_600.0 * alpha
    } else {
        280.0 + 800.0 * (alpha - 0.05)
    }
}
fn reference_storage(alpha: f64) -> f64 {
    if alpha <= 0.05 {
        800.0 * alpha * alpha
    } else {
        let remaining = alpha - 0.05;
        2.0 + 80.0 * remaining + 400.0 * remaining * remaining
    }
}
struct Point {
    alpha: f64,
    axial: f64,
    transverse: f64,
    free: f64,
    internal: f64,
    heat: f64,
}
fn independent(temperature: f64, nominal: f64, old_alpha: f64) -> Point {
    let log_theta = EXPANSION * (temperature - 300.0);
    let theta = log_theta.exp();
    let factor = 1.0 - 0.002 * (temperature - 300.0);
    let natural = theta * old_alpha.exp();
    let elastic_axial = if nominal == 0.0 {
        natural
    } else {
        root(natural, natural * 0.04_f64.exp(), |axial| {
            E * (axial.ln() - log_theta - old_alpha) - nominal * axial
        })
    };
    let (alpha, q, axial) = if nominal * elastic_axial <= factor * reference_yield(old_alpha) {
        (old_alpha, nominal * elastic_axial, elastic_axial)
    } else {
        // Absolute plastic strain is the retained monotonic axial history;
        // P=q/lambda, so a fixed nominal load is not a fixed Kirchhoff stress.
        let alpha = root(old_alpha, 0.2, |alpha| {
            let q = factor * reference_yield(alpha);
            q - nominal * theta * (alpha + q / E).exp()
        });
        let q = factor * reference_yield(alpha);
        (alpha, q, theta * (alpha + q / E).exp())
    };
    let transverse = theta * (-0.5 * alpha - NU * q / E).exp();
    let stored = reference_storage(alpha);
    let free = q * q / (2.0 * E) + factor * stored;
    let entropy = EXPANSION * q + 0.002 * stored;
    Point {
        alpha,
        axial,
        transverse,
        free,
        internal: free + temperature * entropy,
        heat: 200.0 * factor * (alpha - old_alpha),
    }
}

#[test]
fn heated_tetrahedron_yields_under_dead_load_and_retains_permanent_shape_after_release() {
    with_cx(|cx| {
        let controls = FiniteThermoplasticControls {
            material: controls(),
            maximum_numerical_kinetic_loss_j: 50.0,
        };
        let supports: Vec<_> = (0..12)
            .filter(|dof| ![3, 7, 11].contains(dof))
            .map(|dof| (dof, 0.0))
            .collect();
        let mut solid = FiniteThermoplasticTetSolid::try_new(
            &NODES,
            &TETS,
            vec![family()],
            &[300.0],
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        let mut alpha = 0.0;
        for (index, (temperature, nominal)) in PATH.into_iter().enumerate() {
            let expected = independent(temperature, nominal, alpha);
            let forces = if nominal == 0.0 {
                Vec::new()
            } else {
                vec![(3, nominal * VOLUME)]
            };
            let step = solid
                .advance(
                    &[temperature],
                    &supports,
                    &forces,
                    [0.0; 3],
                    0.1,
                    tet_settings(),
                    controls,
                    cx,
                )
                .unwrap();
            assert!(step.equilibrium.steps[0].linear_iterations > 0);
            for (node, x) in NODES.iter().enumerate() {
                for axis in 0..3 {
                    let stretch = if axis == 0 {
                        expected.axial
                    } else {
                        expected.transverse
                    };
                    close(
                        step.solid.displacement_m()[node][axis],
                        (stretch - 1.0) * x[axis],
                        3e-8,
                    );
                }
            }
            close(step.solid.element_states()[0].alpha, expected.alpha, 3e-8);
            close(step.solid.element_first_piola_pa()[0][0], nominal, 3e-6);
            for entry in [4, 8] {
                close(step.solid.element_first_piola_pa()[0][entry], 0.0, 3e-6);
            }
            close(step.solid.free_energy_j(), VOLUME * expected.free, 2e-6);
            close(
                step.solid.internal_energy_j(),
                VOLUME * expected.internal,
                2e-6,
            );
            close(
                step.energy.plastic_dissipation_j,
                VOLUME * expected.heat,
                2e-6,
            );
            close(step.energy.balance_error_j, 0.0, 1e-6);
            if index == 1 {
                assert!(expected.alpha > 0.0 && expected.alpha < 0.05);
            }
            if index == 2 {
                assert!(expected.alpha > 0.05);
            }
            if index >= 3 {
                assert_eq!(step.solid.element_states(), solid.element_states());
            }
            alpha = expected.alpha;
            solid = step.solid;
        }
        assert!(solid.displacement_m()[1][0] > 0.05);
        assert!(solid.displacement_m()[2][1] < 0.0);
        let states = solid.element_states().to_vec();
        let shape = solid.displacement_m().to_vec();
        let stored = (solid.free_energy_j(), solid.internal_energy_j());
        assert!(
            solid
                .advance(
                    &[600.01],
                    &supports,
                    &[(3, 100.0 * VOLUME)],
                    [0.0; 3],
                    0.1,
                    tet_settings(),
                    controls,
                    cx
                )
                .is_err()
        );
        assert!(
            solid
                .advance(
                    &[600.0],
                    &supports,
                    &[(3, 400.0 * VOLUME)],
                    [0.0; 3],
                    0.1,
                    tet_settings(),
                    controls,
                    cx
                )
                .is_err()
        );
        assert_eq!(solid.element_states(), states);
        assert_eq!(solid.displacement_m(), shape);
        assert_eq!(solid.temperatures_k(), &[300.0]);
        assert_eq!((solid.free_energy_j(), solid.internal_energy_j()), stored);
    });
}

#[test]
fn heated_triangular_and_quadrilateral_sheets_solve_lateral_motion_thickness_and_springback() {
    with_cx(|cx| {
        let family = family();
        for mesh in [
            Mesh2::triangles(1.0, 1.0, 1, 1),
            Mesh2::quads(1.0, 1.0, 1, 1),
        ] {
            let thickness = vec![THICKNESS; mesh.elems.len()];
            let initial_temperature = vec![300.0; mesh.elems.len()];
            let mut supports = Vec::new();
            for (node, &[x, y]) in mesh.nodes.iter().enumerate() {
                if x == 0.0 {
                    supports.push((2 * node, 0.0));
                }
                if y == 0.0 {
                    supports.push((2 * node + 1, 0.0));
                }
            }
            let loads = |nominal: f64| -> Vec<_> {
                if nominal == 0.0 {
                    return Vec::new();
                }
                mesh.nodes
                    .iter()
                    .enumerate()
                    .filter(|(_, x)| x[0] == 1.0)
                    .map(|(node, _)| (2 * node, nominal * THICKNESS / 2.0))
                    .collect()
            };
            let mut sheet = FiniteThermoplasticSheet::try_new(
                &mesh,
                FinitePlasticLawField::Uniform(&family),
                &thickness,
                &initial_temperature,
                PlaneStressSettings {
                    relative_tolerance: 1e-13,
                    ..PlaneStressSettings::default()
                },
                SheetAssemblyBudget {
                    maximum_nodes: 6,
                    maximum_elements: 2,
                    maximum_free_dofs: 12,
                    minimum_scaled_jacobian: 1e-4,
                },
                cx,
            )
            .unwrap();
            let mut alpha = 0.0;
            for (index, (temperature, nominal)) in PATH.into_iter().enumerate() {
                let expected = independent(temperature, nominal, alpha);
                let step = sheet
                    .advance(
                        &vec![temperature; mesh.elems.len()],
                        &supports,
                        &loads(nominal),
                        0.1,
                        sheet_settings(),
                        controls(),
                        cx,
                    )
                    .unwrap();
                for (node, &[x, y]) in mesh.nodes.iter().enumerate() {
                    close(
                        step.sheet.displacement_m()[node][0],
                        (expected.axial - 1.0) * x,
                        3e-8,
                    );
                    close(
                        step.sheet.displacement_m()[node][1],
                        (expected.transverse - 1.0) * y,
                        3e-8,
                    );
                }
                for state in step.sheet.quadrature_states() {
                    close(state.alpha, expected.alpha, 3e-8);
                }
                for value in step.sheet.quadrature_thickness_m() {
                    close(value, THICKNESS * expected.transverse, 3e-10);
                }
                for piola in step.sheet.quadrature_first_piola_pa() {
                    close(piola[0], nominal, 3e-6);
                    close(piola[4], 0.0, 3e-6);
                    close(piola[8], 0.0, 3e-6);
                }
                close(step.sheet.free_energy_j(), THICKNESS * expected.free, 2e-7);
                close(
                    step.sheet.internal_energy_j(),
                    THICKNESS * expected.internal,
                    2e-7,
                );
                close(
                    step.energy.plastic_dissipation_j,
                    THICKNESS * expected.heat,
                    2e-7,
                );
                close(step.energy.balance_error_j, 0.0, 1e-6);
                if index >= 3 {
                    assert_eq!(step.sheet.quadrature_states(), sheet.quadrature_states());
                }
                alpha = expected.alpha;
                sheet = step.sheet;
            }
            let states = sheet.quadrature_states().to_vec();
            let shape = sheet.displacement_m().to_vec();
            let stored = (
                sheet.free_energy_j(),
                sheet.internal_energy_j(),
                sheet.physical_time_s(),
            );
            assert!(
                sheet
                    .advance(
                        &vec![600.01; mesh.elems.len()],
                        &supports,
                        &loads(100.0),
                        0.1,
                        sheet_settings(),
                        controls(),
                        cx
                    )
                    .is_err()
            );
            assert!(
                sheet
                    .advance(
                        &vec![600.0; mesh.elems.len()],
                        &supports,
                        &loads(400.0),
                        0.1,
                        sheet_settings(),
                        controls(),
                        cx
                    )
                    .is_err()
            );
            assert_eq!(sheet.quadrature_states(), states);
            assert_eq!(sheet.displacement_m(), shape);
            assert_eq!(sheet.temperatures_k(), initial_temperature);
            assert_eq!(
                (
                    sheet.free_energy_j(),
                    sheet.internal_energy_j(),
                    sheet.physical_time_s()
                ),
                stored
            );
        }
    });
}
