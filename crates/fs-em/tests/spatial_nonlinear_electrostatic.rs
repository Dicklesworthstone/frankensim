//! Independent original-law charge, storage and tangent oracles in 3D.
//! Authored constitutive fixtures do not represent qualified material data.
use fs_em::{
    electrostatic::{
        ElectrostaticControls, ElectrostaticError,
        nonlinear::{
            ElectrodeControl, NonlinearElectrostaticControls, NonlinearElectrostaticError,
            spatial::{
                SpatialDielectricMaterial, SpatialNonlinearElectrostatics,
                tangent::DifferentialCapacitanceControls,
            },
        },
    },
    spatial::electrostatic::{
        ElectrostaticVoltageFace, SpatialElectrostaticBoundary, SpatialElectrostatics,
    },
};
use fs_evidence::ValidityDomain;
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    field_storage::{
        ELECTRIC_PERMITTIVITY_DIMS, ELECTRIC_PERMITTIVITY_PROPERTY,
        dielectric::reversible::{
            ELECTRIC_DISPLACEMENT_DIMS, ELECTRIC_FIELD_DIMS,
            REVERSIBLE_ELECTRIC_DISPLACEMENT_PROPERTY, ReversibleDielectricMaterial,
        },
        spatial::{PRINCIPAL_ELECTRIC_PERMITTIVITY_PROPERTIES, SpatialDielectricLaw},
    },
    state_point::{
        MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
        ScalarPropertyRequirement, resolve_material_state_point,
    },
};
use fs_qty::QuantitySpec;
use fs_solver::CheckedCgConfig;

const GAP: f64 = 1e-3;
const WIDTH: f64 = 2e-2;
const DEPTH: f64 = 5e-3;
const AREA: f64 = WIDTH * DEPTH;
const E0: f64 = 1e5;
const D0: f64 = 1e-6;
const EPSILON: f64 = D0 / E0;
const V0: f64 = GAP * E0;
const Q0: f64 = AREA * D0;
const U0: f64 = AREA * GAP * E0 * D0;
const C0: f64 = EPSILON * AREA / GAP;

fn near(actual: f64, expected: f64, scale: f64) {
    let allowance = 3e-8 * scale.abs().max(expected.abs());
    assert!(
        (actual - expected).abs() <= allowance,
        "{actual:e} != {expected:e}, allowance={allowance:e}"
    );
}
fn near_vector(actual: [f64; 3], expected: [f64; 3], scale: f64) {
    for i in 0..3 {
        near(actual[i], expected[i], scale);
    }
}
fn limits() -> WeightedAssemblyLimits {
    WeightedAssemblyLimits {
        max_cells: 1024,
        max_dofs: 8192,
        max_triplets: 36 * 1024,
    }
}
fn linear_controls() -> CheckedCgConfig {
    CheckedCgConfig {
        tolerance: 1e-12,
        max_iterations: 512,
        max_corrections: 2,
    }
}
fn controls() -> NonlinearElectrostaticControls {
    NonlinearElectrostaticControls {
        absolute_residual_tolerance_c: 1e-22,
        relative_residual_tolerance: 1e-12,
        charge_tolerance_c: 1e-20,
        energy_tolerance_j: 1e-16,
        max_newton_iterations: 32,
        max_backtracks: 24,
        max_operator_applications: 8192,
        linear: linear_controls(),
    }
}
fn differential_controls() -> DifferentialCapacitanceControls {
    DifferentialCapacitanceControls {
        maximum_terminals: 8,
        maximum_matrix_entries: 64,
        max_operator_applications: 8192,
        coefficient_tolerance_f: 1e-21,
        linear: linear_controls(),
    }
}
fn scalar(properties: &[(&str, f64)]) -> ResolvedMaterialStatePoint {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for &(name, value) in properties {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, ELECTRIC_PERMITTIVITY_DIMS),
                value: PropertyValue::Scalar {
                    value,
                    dims: ELECTRIC_PERMITTIVITY_DIMS,
                },
                validity: ValidityDomain::unconstrained(),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "authored 3D dielectric fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_new(
                name,
                ELECTRIC_PERMITTIVITY_DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )
            .unwrap(),
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic dielectric".into(),
            phase: "solid".into(),
            process: "3D test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    resolve_material_state_point(
        &card,
        &QueryPoint::new(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn linear(epsilon: f64) -> ResolvedMaterialStatePoint {
    scalar(&[(ELECTRIC_PERMITTIVITY_PROPERTY, epsilon)])
}
fn nonlinear() -> ReversibleDielectricMaterial {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(
                REVERSIBLE_ELECTRIC_DISPLACEMENT_PROPERTY,
                ELECTRIC_DISPLACEMENT_DIMS,
            ),
            value: PropertyValue::Curve {
                abscissa: "E".into(),
                abscissa_dims: ELECTRIC_FIELD_DIMS,
                knots: [(0.0, 0.0), (1.0, 2.0), (3.0, 3.0), (4.0, 7.0)]
                    .map(|(e, d)| (e * E0, d * D0))
                    .to_vec(),
                dims: ELECTRIC_DISPLACEMENT_DIMS,
            },
            validity: ValidityDomain::unconstrained().with_quantity(
                "E",
                QuantitySpec::dimensional(ELECTRIC_FIELD_DIMS),
                0.0,
                4.0 * E0,
            ),
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "authored piecewise-affine absolute D(E)".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic nonlinear dielectric".into(),
            phase: "solid".into(),
            process: "3D test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let context = QueryPoint::new()
        .with_quantity("E", QuantitySpec::dimensional(ELECTRIC_FIELD_DIMS), 0.0)
        .unwrap();
    ReversibleDielectricMaterial::resolve(
        &card,
        &context,
        "E",
        4.0 * E0,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn box_mesh(
    xs: &[f64],
    width: f64,
    depth: f64,
    ny: usize,
    nz: usize,
) -> (Vec<[f64; 3]>, Vec<[u32; 4]>) {
    let mut positions = Vec::new();
    for &x in xs {
        for y in 0..=ny {
            for z in 0..=nz {
                positions.push([
                    x,
                    width * y as f64 / ny as f64,
                    depth * z as f64 / nz as f64,
                ]);
            }
        }
    }
    let index = |x: usize, y: usize, z: usize| ((x * (ny + 1) + y) * (nz + 1) + z) as u32;
    let mut cells = Vec::new();
    for x in 0..xs.len() - 1 {
        for y in 0..ny {
            for z in 0..nz {
                let [a, b, c, d, e, f, g, h] = [
                    index(x, y, z),
                    index(x + 1, y, z),
                    index(x, y + 1, z),
                    index(x + 1, y + 1, z),
                    index(x, y, z + 1),
                    index(x + 1, y, z + 1),
                    index(x, y + 1, z + 1),
                    index(x + 1, y + 1, z + 1),
                ];
                cells.extend([
                    [a, b, d, h],
                    [a, d, c, h],
                    [a, c, g, h],
                    [a, g, e, h],
                    [a, e, f, h],
                    [a, f, b, h],
                ]);
            }
        }
    }
    (positions, cells)
}
fn append(first: &mut (Vec<[f64; 3]>, Vec<[u32; 4]>), second: (Vec<[f64; 3]>, Vec<[u32; 4]>)) {
    let offset = first.0.len() as u32;
    first.0.extend(second.0);
    first
        .1
        .extend(second.1.into_iter().map(|cell| cell.map(|v| v + offset)));
}
fn faces_at(
    problem: &SpatialNonlinearElectrostatics,
    original_positions: &[[f64; 3]],
    x: f64,
) -> Vec<usize> {
    problem
        .exterior_faces()
        .into_iter()
        .filter(|&face| {
            problem.mesh().faces[face]
                .iter()
                .all(|&v| original_positions[v as usize][0] == x)
        })
        .collect()
}

#[test]
fn rotated_nonlinear_plate_keeps_original_storage_charge_and_three_dimensional_tangent() {
    let (original, cells) = box_mesh(&[0.0, 0.5 * GAP, GAP], WIDTH, DEPTH, 2, 2);
    let axes = [
        [1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0],
        [2.0 / 3.0, 1.0 / 3.0, -2.0 / 3.0],
        [-2.0 / 3.0, 2.0 / 3.0, -1.0 / 3.0],
    ];
    let positions = original
        .iter()
        .map(|point| {
            std::array::from_fn(|i| {
                point[0].mul_add(
                    axes[0][i],
                    point[1].mul_add(axes[1][i], point[2] * axes[2][i]),
                )
            })
        })
        .collect();
    let laws = vec![SpatialDielectricMaterial::Reversible(nonlinear()); cells.len()];
    let problem =
        SpatialNonlinearElectrostatics::new(positions, cells, laws, limits(), &mut || false)
            .unwrap();
    let terminals = [
        faces_at(&problem, &original, 0.0),
        faces_at(&problem, &original, GAP),
    ];
    let electrodes = problem
        .prepare_electrodes(&[&terminals[0], &terminals[1]], &mut || false)
        .unwrap();
    let rho = vec![0.0; problem.mesh().tets.len()];
    let initial: Vec<_> = original.iter().map(|point| 1.8 * E0 * point[0]).collect();
    let field = electrodes
        .solve(
            &rho,
            &[
                ElectrodeControl::Voltage(0.0),
                ElectrodeControl::Voltage(2.0 * V0),
            ],
            &initial,
            controls(),
            &mut || false,
        )
        .unwrap();
    assert!(field.report().newton_iterations > 0);
    near(field.terminal_charge_c()[0], -2.5 * Q0, Q0);
    near(field.terminal_charge_c()[1], 2.5 * Q0, Q0);
    near(field.energy_j(), 1.75 * U0, U0);
    near(field.coenergy_j(), 3.25 * U0, U0);
    assert!(
        (field.energy_j() - 0.5 * field.terminal_voltage_v()[1] * field.terminal_charge_c()[1])
            .abs()
            > 0.5 * U0
    );
    for cell in 0..problem.mesh().tets.len() {
        near_vector(
            field.electric_field_v_m()[cell],
            axes[0].map(|v| -2.0 * E0 * v),
            E0,
        );
        near_vector(
            field.displacement_c_m2()[cell],
            axes[0].map(|v| -2.5 * D0 * v),
            D0,
        );
        for i in 0..3 {
            for j in 0..3 {
                let expected =
                    EPSILON * ((if i == j { 1.25 } else { 0.0 }) - 0.75 * axes[0][i] * axes[0][j]);
                near(
                    field.differential_permittivity_f_m()[cell][i][j],
                    expected,
                    EPSILON,
                );
            }
        }
    }
    let differential = field
        .differential_capacitance(differential_controls(), &mut || false)
        .unwrap();
    for (actual, expected) in differential.matrix_f().iter().zip([0.5, -0.5, -0.5, 0.5]) {
        near(*actual, expected * C0, C0);
    }
    let charged = electrodes
        .solve(
            &rho,
            &[
                ElectrodeControl::Voltage(0.0),
                ElectrodeControl::Charge(2.5 * Q0),
            ],
            &vec![0.0; original.len()],
            controls(),
            &mut || false,
        )
        .unwrap();
    near(charged.terminal_voltage_v()[1], 2.0 * V0, V0);
    near(charged.energy_j(), field.energy_j(), U0);
    let delta = 1e-3 * V0;
    let mut charge = [0.0; 2];
    for (slot, sign) in [-1.0, 1.0].into_iter().enumerate() {
        let driven = electrodes
            .solve(
                &rho,
                &[
                    ElectrodeControl::Voltage(0.0),
                    ElectrodeControl::Voltage(2.0 * V0 + sign * delta),
                ],
                field.relative_voltage_v(),
                controls(),
                &mut || false,
            )
            .unwrap();
        charge[slot] = driven.terminal_charge_c()[1];
    }
    near((charge[1] - charge[0]) / (2.0 * delta), 0.5 * C0, C0);
}

#[test]
fn oriented_linear_volume_charge_matches_the_native_three_dimensional_hodge() {
    let (positions, cells) = box_mesh(&[0.0, 0.5 * GAP, GAP], GAP, GAP, 2, 2);
    let names = PRINCIPAL_ELECTRIC_PERMITTIVITY_PROPERTIES;
    let material = scalar(&[
        (names[0], 2.0 * EPSILON),
        (names[1], 5.0 * EPSILON),
        (names[2], 7.0 * EPSILON),
    ]);
    let s = 0.5_f64.sqrt();
    let law =
        SpatialDielectricLaw::orthotropic([[s, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, s]]).unwrap();
    let materials = vec![material.clone(); cells.len()];
    let dielectric_laws = vec![law; cells.len()];
    let problem = SpatialNonlinearElectrostatics::new(
        positions.clone(),
        cells.clone(),
        vec![SpatialDielectricMaterial::Linear { material, law }; cells.len()],
        limits(),
        &mut || false,
    )
    .unwrap();
    let native = SpatialElectrostatics::new_with_laws(
        positions.clone(),
        cells,
        materials,
        dielectric_laws,
        limits(),
        &mut || false,
    )
    .unwrap();
    let terminals = [
        faces_at(&problem, &positions, 0.0),
        faces_at(&problem, &positions, GAP),
    ];
    let electrodes = problem
        .prepare_electrodes(&[&terminals[0], &terminals[1]], &mut || false)
        .unwrap();
    let voltages = [0.0, 0.75 * V0];
    let rho = vec![0.7 * D0 / GAP; problem.mesh().tets.len()];
    let field = electrodes
        .solve(
            &rho,
            &voltages.map(ElectrodeControl::Voltage),
            &vec![0.0; positions.len()],
            controls(),
            &mut || false,
        )
        .unwrap();
    let boundary = SpatialElectrostaticBoundary {
        voltage: terminals
            .iter()
            .zip(voltages)
            .flat_map(|(faces, voltage)| {
                faces.iter().map(move |&face| ElectrostaticVoltageFace {
                    face,
                    voltage_v: [voltage; 3],
                })
            })
            .collect(),
        displacement: vec![],
    };
    let exact = native
        .solve(
            &rho,
            &boundary,
            ElectrostaticControls {
                linear: linear_controls(),
                charge_tolerance_c: 1e-20,
                energy_tolerance_j: 1e-16,
            },
            &mut || false,
        )
        .unwrap();
    for (a, b) in field
        .electric_field_v_m()
        .iter()
        .zip(exact.electric_field_v_m())
    {
        near_vector(*a, *b, E0);
    }
    for (a, b) in field
        .displacement_c_m2()
        .iter()
        .zip(exact.displacement_c_m2())
    {
        near_vector(*a, *b, D0);
    }
    near(field.energy_j(), exact.energy_j(), exact.energy_j());
    near(field.coenergy_j(), exact.energy_j(), exact.energy_j());
    near(
        field.volume_charge_c(),
        0.7 * D0 * GAP * GAP,
        D0 * GAP * GAP,
    );
    for (terminal, trace) in electrodes.terminal_vertices().iter().enumerate() {
        let charge: f64 = trace
            .iter()
            .map(|&v| exact.nodal_boundary_charge_c()[v].unwrap())
            .sum();
        near(field.terminal_charge_c()[terminal], charge, D0 * GAP * GAP);
    }
    assert!(field.maximum_component_charge_residual_c() <= controls().charge_tolerance_c);
}

#[test]
fn floating_middle_metal_keeps_trapped_charge_and_the_condensed_nonlinear_tangent() {
    let mut mesh = box_mesh(&[0.0, 0.5 * GAP, GAP], WIDTH, DEPTH, 2, 2);
    let first_cells = mesh.1.len();
    append(
        &mut mesh,
        box_mesh(&[1.5 * GAP, 2.0 * GAP, 2.5 * GAP], WIDTH, DEPTH, 2, 2),
    );
    let (positions, cells) = mesh;
    let mut materials = vec![SpatialDielectricMaterial::Reversible(nonlinear()); first_cells];
    materials.extend(vec![
        SpatialDielectricMaterial::isotropic(linear(EPSILON));
        cells.len() - first_cells
    ]);
    let problem = SpatialNonlinearElectrostatics::new(
        positions.clone(),
        cells,
        materials,
        limits(),
        &mut || false,
    )
    .unwrap();
    let lower = faces_at(&problem, &positions, 0.0);
    let mut middle = faces_at(&problem, &positions, GAP);
    middle.extend(faces_at(&problem, &positions, 1.5 * GAP));
    let upper = faces_at(&problem, &positions, 2.5 * GAP);
    let electrodes = problem
        .prepare_electrodes(&[&lower, &middle, &upper], &mut || false)
        .unwrap();
    let rho = vec![0.0; problem.mesh().tets.len()];
    let boundary = [
        ElectrodeControl::Voltage(0.0),
        ElectrodeControl::Charge(1.5 * Q0),
        ElectrodeControl::Voltage(3.0 * V0),
    ];
    let field = electrodes
        .solve(
            &rho,
            &boundary,
            &vec![0.0; positions.len()],
            controls(),
            &mut || false,
        )
        .unwrap();
    near(field.terminal_voltage_v()[1], 2.0 * V0, V0);
    for (charge, expected) in field.terminal_charge_c().iter().zip([-2.5, 1.5, 1.0]) {
        near(*charge, expected * Q0, Q0);
    }
    near(field.energy_j(), 2.25 * U0, U0);
    near(field.coenergy_j(), 3.75 * U0, U0);
    assert!(
        field
            .reference_terminal_by_vertex()
            .iter()
            .all(|&reference| reference == 0)
    );
    assert!(field.maximum_floating_charge_residual_c() <= controls().charge_tolerance_c);
    let differential = field
        .differential_capacitance(differential_controls(), &mut || false)
        .unwrap();
    assert_eq!(differential.terminals(), &[0, 2]);
    for (actual, sign) in differential.matrix_f().iter().zip([1.0, -1.0, -1.0, 1.0]) {
        near(*actual, sign * C0 / 3.0, C0);
    }
    let delta = 1e-3 * V0;
    let mut charge = [0.0; 2];
    for (slot, sign) in [-1.0, 1.0].into_iter().enumerate() {
        let driven = electrodes
            .solve(
                &rho,
                &[
                    boundary[0],
                    boundary[1],
                    ElectrodeControl::Voltage(3.0 * V0 + sign * delta),
                ],
                field.relative_voltage_v(),
                controls(),
                &mut || false,
            )
            .unwrap();
        near(driven.terminal_charge_c()[1], 1.5 * Q0, Q0);
        charge[slot] = driven.terminal_charge_c()[2];
    }
    near((charge[1] - charge[0]) / (2.0 * delta), C0 / 3.0, C0);
    let shifted = electrodes
        .solve(
            &rho,
            &[
                ElectrodeControl::Voltage(1e6),
                boundary[1],
                ElectrodeControl::Voltage(1e6 + 3.0 * V0),
            ],
            field.relative_voltage_v(),
            controls(),
            &mut || false,
        )
        .unwrap();
    for (&a, &b) in shifted
        .terminal_charge_c()
        .iter()
        .zip(field.terminal_charge_c())
    {
        near(a, b, Q0);
    }
    near(shifted.energy_j(), field.energy_j(), U0);
}

#[test]
fn disconnected_fixed_charge_fields_keep_tiny_voltages_under_independent_large_biases() {
    let mut mesh = box_mesh(&[0.0, 0.5 * GAP, GAP], WIDTH, DEPTH, 1, 1);
    append(
        &mut mesh,
        box_mesh(&[2.0 * GAP, 2.5 * GAP, 3.0 * GAP], WIDTH, DEPTH, 1, 1),
    );
    let (positions, cells) = mesh;
    let materials = vec![SpatialDielectricMaterial::Reversible(nonlinear()); cells.len()];
    let problem = SpatialNonlinearElectrostatics::new(
        positions.clone(),
        cells,
        materials,
        limits(),
        &mut || false,
    )
    .unwrap();
    let terminals = [0.0, GAP, 2.0 * GAP, 3.0 * GAP].map(|x| faces_at(&problem, &positions, x));
    let traces: Vec<&[usize]> = terminals.iter().map(Vec::as_slice).collect();
    let electrodes = problem.prepare_electrodes(&traces, &mut || false).unwrap();
    let rho = vec![0.0; problem.mesh().tets.len()];
    let huge = 2.0_f64.powi(50);
    let voltage = [1e-6, -3e-6];
    let charge = voltage.map(|v| 2.0 * C0 * v);
    let mut settings = controls();
    settings.absolute_residual_tolerance_c = 1e-30;
    settings.charge_tolerance_c = 1e-28;
    settings.energy_tolerance_j = 1e-34;
    let field = electrodes
        .solve(
            &rho,
            &[
                ElectrodeControl::Voltage(huge),
                ElectrodeControl::Charge(charge[0]),
                ElectrodeControl::Voltage(-huge),
                ElectrodeControl::Charge(charge[1]),
            ],
            &vec![0.0; positions.len()],
            settings,
            &mut || false,
        )
        .unwrap();
    near(field.terminal_voltage_v()[1], voltage[0], voltage[0]);
    near(field.terminal_voltage_v()[3], voltage[1], voltage[1]);
    near(field.terminal_charge_c()[1], charge[0], charge[0]);
    near(field.terminal_charge_c()[3], charge[1], charge[1]);
    assert_eq!(field.reference_voltage_for_terminal_v(1), Some(huge));
    assert_eq!(field.reference_voltage_for_terminal_v(3), Some(-huge));
    assert_eq!(field.reference_terminal_for_terminal(1), Some(0));
    assert_eq!(field.reference_terminal_for_terminal(3), Some(2));
    near(
        field.energy_j(),
        C0 * (voltage[0] * voltage[0] + voltage[1] * voltage[1]),
        C0 * 1e-12,
    );
    assert!(field.maximum_component_charge_residual_c() <= settings.charge_tolerance_c);
    // Changing the sole reference of either fixed-charge capacitor is a gauge
    // change. Its local terminal response must remain zero independently.
    let differential = field
        .differential_capacitance(differential_controls(), &mut || false)
        .unwrap();
    assert_eq!(differential.terminals(), &[0, 2]);
    assert!(
        differential
            .matrix_f()
            .iter()
            .all(|value| value.abs() <= 1e-21)
    );
    assert!(
        electrodes
            .solve(
                &rho,
                &[
                    ElectrodeControl::Voltage(0.0),
                    ElectrodeControl::Charge(charge[0]),
                    ElectrodeControl::Charge(-charge[1]),
                    ElectrodeControl::Charge(charge[1]),
                ],
                &vec![0.0; positions.len()],
                settings,
                &mut || false
            )
            .is_err()
    );
}

#[test]
fn one_grounded_electrode_and_volume_charge_use_exact_tetrahedral_source_moments() {
    let positions = vec![
        [0.0, 0.0, 0.0],
        [GAP, 0.0, 0.0],
        [0.0, WIDTH, 0.0],
        [0.0, 0.0, DEPTH],
    ];
    let problem = SpatialNonlinearElectrostatics::new(
        positions.clone(),
        vec![[0, 1, 2, 3]],
        vec![SpatialDielectricMaterial::Reversible(nonlinear())],
        limits(),
        &mut || false,
    )
    .unwrap();
    let ground = faces_at(&problem, &positions, 0.0);
    let electrodes = problem
        .prepare_electrodes(&[&ground], &mut || false)
        .unwrap();
    // On this original P1 tetrahedron the sole free shape function is x/G.
    // Its exact weak equation is Vcell*D/G = rho*Vcell/4, hence D=rho*G/4.
    // This is a discrete oracle, not a claim that the raw cell field enforces
    // pointwise zero flux on every slanted face of the continuum tetrahedron.
    let rho = 10.0 * D0 / GAP;
    let volume = AREA * GAP / 6.0;
    let field = electrodes
        .solve(
            &[rho],
            &[ElectrodeControl::Voltage(0.0)],
            &[0.0; 4],
            controls(),
            &mut || false,
        )
        .unwrap();
    near(field.relative_voltage_v()[1], 2.0 * V0, V0);
    near_vector(field.electric_field_v_m()[0], [-2.0 * E0, 0.0, 0.0], E0);
    near_vector(field.displacement_c_m2()[0], [-2.5 * D0, 0.0, 0.0], D0);
    near(field.volume_charge_c(), rho * volume, Q0);
    near(field.terminal_charge_c()[0], -rho * volume, Q0);
    near(field.energy_j(), 1.75 * U0 / 6.0, U0);
    near(field.coenergy_j(), 3.25 * U0 / 6.0, U0);
    near(
        field.energy_j() + field.coenergy_j(),
        field.relative_voltage_v()[1] * rho * volume / 4.0,
        U0,
    );
    assert!(field.energy_identity_residual_j().abs() <= controls().energy_tolerance_j);
    assert!(field.maximum_free_charge_residual_c() <= controls().charge_tolerance_c);
    let differential = field
        .differential_capacitance(differential_controls(), &mut || false)
        .unwrap();
    assert_eq!(differential.terminals(), &[0]);
    assert!(differential.matrix_f()[0].abs() <= differential_controls().coefficient_tolerance_f);
}

#[test]
fn source_kinks_endpoints_and_bounded_work_refuse_without_corrupting_reuse() {
    // Integer geometry makes the exact source knots representable without
    // confusing a near-knot numerical field with an exact derivative corner.
    let (positions, cells) = box_mesh(&[0.0, 1.0], 1.0, 1.0, 1, 1);
    let materials = vec![SpatialDielectricMaterial::Reversible(nonlinear()); cells.len()];
    assert!(matches!(
        SpatialNonlinearElectrostatics::new(
            positions.clone(),
            cells.clone(),
            materials.clone(),
            limits(),
            &mut || true
        ),
        Err(NonlinearElectrostaticError::Electric(
            ElectrostaticError::Cancelled
        ))
    ));
    let mut too_small = limits();
    too_small.max_cells = 1;
    assert!(
        SpatialNonlinearElectrostatics::new(
            positions.clone(),
            cells.clone(),
            materials.clone(),
            too_small,
            &mut || false
        )
        .is_err()
    );
    let problem = SpatialNonlinearElectrostatics::new(
        positions.clone(),
        cells,
        materials,
        limits(),
        &mut || false,
    )
    .unwrap();
    let terminals = [
        faces_at(&problem, &positions, 0.0),
        faces_at(&problem, &positions, 1.0),
    ];
    assert!(problem.prepare_electrodes(&[], &mut || false).is_err());
    assert!(
        problem
            .prepare_electrodes(&[&terminals[0], &terminals[0]], &mut || false)
            .is_err()
    );
    let exterior = problem.exterior_faces();
    let internal = (0..problem.mesh().faces.len())
        .find(|face| !exterior.contains(face))
        .unwrap();
    assert!(
        problem
            .prepare_electrodes(&[&terminals[0], &[internal]], &mut || false)
            .is_err()
    );
    let electrodes = problem
        .prepare_electrodes(&[&terminals[0], &terminals[1]], &mut || false)
        .unwrap();
    let rho = vec![0.0; problem.mesh().tets.len()];
    let initial = vec![0.0; positions.len()];
    let mut settings = controls();
    settings.charge_tolerance_c = 1e-16;
    settings.energy_tolerance_j = 1e-12;
    let kink = electrodes
        .solve(
            &rho,
            &[
                ElectrodeControl::Voltage(0.0),
                ElectrodeControl::Voltage(E0),
            ],
            &initial,
            settings,
            &mut || false,
        )
        .unwrap();
    assert!(
        kink.differential_capacitance(differential_controls(), &mut || false)
            .is_err()
    );
    let endpoint = electrodes
        .solve(
            &rho,
            &[
                ElectrodeControl::Voltage(0.0),
                ElectrodeControl::Voltage(4.0 * E0),
            ],
            &initial,
            settings,
            &mut || false,
        )
        .unwrap();
    assert!(
        endpoint
            .differential_capacitance(differential_controls(), &mut || false)
            .is_err()
    );
    assert!(matches!(
        electrodes.solve(
            &rho,
            &[
                ElectrodeControl::Voltage(0.0),
                ElectrodeControl::Voltage(4.01 * E0)
            ],
            &initial,
            settings,
            &mut || false
        ),
        Err(NonlinearElectrostaticError::FieldRange { .. })
    ));
    let charged = [
        ElectrodeControl::Voltage(0.0),
        ElectrodeControl::Charge(2.5 * D0),
    ];
    let mut no_work = settings;
    no_work.max_operator_applications = 0;
    assert!(matches!(
        electrodes.solve(&rho, &charged, &initial, no_work, &mut || false),
        Err(NonlinearElectrostaticError::OperatorBudget)
    ));
    no_work = settings;
    no_work.max_newton_iterations = 0;
    assert!(matches!(
        electrodes.solve(&rho, &charged, &initial, no_work, &mut || false),
        Err(NonlinearElectrostaticError::Nonconvergence(_))
    ));
    for stop_after in [1, 20] {
        let mut calls = 0;
        let cancelled = electrodes.solve(&rho, &charged, &initial, settings, &mut || {
            calls += 1;
            calls >= stop_after
        });
        assert!(matches!(
            cancelled,
            Err(NonlinearElectrostaticError::Electric(
                ElectrostaticError::Cancelled
            ))
        ));
    }
    let accepted = electrodes
        .solve(&rho, &charged, &initial, settings, &mut || false)
        .unwrap();
    near(accepted.terminal_voltage_v()[1], 2.0 * E0, E0);
    near(accepted.terminal_charge_c()[1], 2.5 * D0, D0);
    assert!(
        electrodes
            .solve(
                &rho,
                &[
                    ElectrodeControl::Charge(-2.5 * D0),
                    ElectrodeControl::Charge(2.5 * D0)
                ],
                &initial,
                settings,
                &mut || false
            )
            .is_err()
    );
    let mut bad_initial = initial.clone();
    bad_initial[electrodes.terminal_vertices()[1][0]] = 1.0;
    assert!(
        electrodes
            .solve(&rho, &charged, &bad_initial, settings, &mut || false)
            .is_err()
    );
    let mut differential_budget = differential_controls();
    differential_budget.max_operator_applications = 0;
    assert!(
        accepted
            .differential_capacitance(differential_budget, &mut || false)
            .is_err()
    );
}
