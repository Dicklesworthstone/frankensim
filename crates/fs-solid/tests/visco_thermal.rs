//! Native imposed-expansion histories and their separate mechanical work port.
//! Synthetic source cards and analytic fixtures are software checks, not data.

use fs_blake3::ContentHash;
use fs_evidence::ValidityDomain;
use fs_exec::Cx;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::state_point::{
    INVERSE_TEMPERATURE_DIMS, LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
    MaterialPropertySelection,
};
use fs_material::visco::tensor::{Mandel6, PronyBranch, TensorProny};
use fs_material::visco_tensor::source::{
    ResolvedTensorMaterial, TensorBranchProperties, TensorElasticity, TensorPronyProperties,
};
use fs_qty::{Density, Dims, Pressure, Temperature};
use fs_solid::linear::Formulation;
use fs_solid::linear3::hex::{
    HEX20_REFERENCE_NODES, HEX27_REFERENCE_NODES, HexConnectivity, HexQuadrature,
};
use fs_solid::linear3::source::ThermalExpansionModel;
use fs_solid::linear3::tet10::{TetQuadrature, promote_p1_to_p2};
use fs_solid::linear3::visco::*;
use fs_solid::linear3::{
    TetAssemblyBudget, TetElasticError, TetStaticSolveConfig, TetThermalStrainField,
    TetThermalStrainState,
};

const TET_NODES: [[f64; 3]; 4] =
    [[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
const HEX_NODES: [[f64; 3]; 8] = [
    [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [1.0, 1.0, 1.0], [0.0, 1.0, 1.0],
];
const HEXES: [[usize; 8]; 1] = [[0, 1, 2, 3, 4, 5, 6, 7]];

fn with_cx(f: impl FnOnce(&Cx<'_>, &fs_exec::CancelGate)) {
    let gate = fs_exec::CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate, arena,
            fs_exec::StreamKey { seed: 925, kernel_id: 925, tile: 0, iteration: 0 },
            fs_exec::Budget::INFINITE, fs_exec::ExecMode::Deterministic,
        ), &gate);
    });
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}; tolerance={tolerance:e}");
}
fn close(actual: f64, expected: f64) {
    near(actual, expected, 3e-10 * actual.abs().max(expected.abs()).max(1e-14));
}
fn diagonal(values: [f64; 6]) -> Mandel6 {
    core::array::from_fn(|i| core::array::from_fn(|j| if i == j { values[i] } else { 0.0 }))
}
fn outer(v: [f64; 6], scale: f64) -> Mandel6 {
    core::array::from_fn(|i| core::array::from_fn(|j| scale * v[i] * v[j]))
}
fn dot6(a: [f64; 6], b: [f64; 6]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn action(c: Mandel6, x: [f64; 6]) -> [f64; 6] {
    core::array::from_fn(|i| c[i].iter().zip(x).map(|(a, x)| a * x).sum())
}
fn config() -> TetStaticSolveConfig {
    TetStaticSolveConfig { maximum_iterations: 300, relative_residual_tolerance: 2e-11 }
}
fn dynamic_config() -> ViscoDynamicConfig {
    ViscoDynamicConfig {
        maximum_steps: 100, maximum_free_dofs: 100, solve: config(),
        absolute_energy_tolerance_j: 1e-13, relative_energy_tolerance: 2e-10,
    }
}
fn prepare(
    mesh: SolidMeshRef<'_>, materials: TetViscoField<'_>, fixed: &[usize], cx: &Cx<'_>,
) -> ViscoElements {
    ViscoElements::new(mesh, materials, fixed, TetAssemblyBudget::standard(),
        ViscoStateBudget { max_internal_variables: 100_000 }, cx).unwrap()
}
fn tet(material: &TetViscoMaterial, free: &[usize], cx: &Cx<'_>) -> ViscoElements {
    let fixed: Vec<_> = (0..12).filter(|d| !free.contains(d)).collect();
    prepare(SolidMeshRef::Tet4 { nodes_m: &TET_NODES, elements: &TETS },
        TetViscoField::Uniform(material), &fixed, cx)
}
fn scalar_material(e: f64, b: f64, tau: f64, eta: f64, limit: f64) -> TetViscoMaterial {
    let branches = if b == 0.0 { vec![] } else {
        vec![PronyBranch::new(diagonal([b, 0.0, 0.0, 0.0, 0.0, 0.0]), tau).unwrap()]
    };
    TetViscoMaterial {
        law: TensorProny::new_mandel(diagonal([e; 6]), branches, limit).unwrap()
            .with_kelvin_viscosity_mandel(diagonal([eta, 0.0, 0.0, 0.0, 0.0, 0.0])).unwrap(),
        density_kg_m3: 1000.0, material_state_identity: [92; 32],
    }
}
fn thermal(material: &TetViscoMaterial, free: [f64; 6]) -> TetThermalStrainState {
    TetThermalStrainState::try_new(
        300.0 + 10_000.0 * free[0], 300.0, free, material.material_state_identity,
        ContentHash([93; 32]), 0.025,
    ).unwrap()
}
fn query(t: f64) -> QueryPoint {
    QueryPoint::new().with("T", t).unwrap().with("omega", 0.0).unwrap()
}
fn scalar_claim(name: &str, dims: Dims, value: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Scalar { value, dims },
        validity: ValidityDomain::unconstrained().with("T", 300.0, 400.0)
            .with("omega", 0.0, 1000.0),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![],
        provenance: Provenance {
            source: "synthetic frozen-coefficient thermal-history fixture".into(),
            license: "CC0-1.0".into(), artifact: None,
        },
    }
}
fn source_card() -> (MaterialCard, TensorPronyProperties) {
    let branch = TensorBranchProperties {
        stiffness_upper: core::array::from_fn(|i| format!("memory_{i}")),
        relaxation_time: "tau".into(),
    };
    let mut claims = ClaimSet::new();
    for c in [
        scalar_claim("young_modulus", Pressure::DIMS, 120_000.0),
        scalar_claim("poisson_ratio", Dims::NONE, 0.2),
        scalar_claim("density", Density::DIMS, 27.0),
        scalar_claim("strain_limit", Dims::NONE, 0.025),
        scalar_claim("tau", Dims([0, 0, 1, 0, 0, 0]), 0.2),
    ] { claims.insert_claim(c).unwrap(); }
    let mut index = 0;
    for i in 0..6 {
        for j in i..6 {
            claims.insert_claim(scalar_claim(
                &branch.stiffness_upper[index], Pressure::DIMS,
                if i == j { 15_000.0 } else { 0.0 },
            )).unwrap();
            index += 1;
        }
    }
    let mut alpha = scalar_claim(
        LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, INVERSE_TEMPERATURE_DIMS, 1e-5,
    );
    alpha.value = PropertyValue::Curve {
        abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
        knots: vec![(300.0, 1e-5), (340.0, 3e-5), (400.0, 2e-5)],
        dims: INVERSE_TEMPERATURE_DIMS,
    };
    alpha.interpolation = InterpolationPolicy::LinearInside;
    claims.insert_claim(alpha).unwrap();
    (MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured solid".into(), phase: "solid".into(),
            process: "analytic thermal-history fixture".into(), revision: 0,
        }, claims, vec![],
    ).unwrap(), TensorPronyProperties {
        elasticity: TensorElasticity::Isotropic, strain_limit: "strain_limit".into(),
        branches: vec![branch],
    })
}

#[test]
fn original_integrated_source_drives_free_expansion_on_every_native_element_family() {
    with_cx(|cx, _| {
        let (card, properties) = source_card();
        let source = ResolvedTensorMaterial::resolve(
            &card, &query(300.0), properties, MaterialPropertySelection::SingleClaimOnly,
        ).unwrap();
        let rotation = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let material = TetViscoMaterial::from_source(&source, rotation).unwrap();
        let reference = material.equilibrium_material().unwrap();
        let expansion = ThermalExpansionModel::Isotropic.resolve(
            &card, &query(300.0), &query(400.0), MaterialPropertySelection::SingleClaimOnly,
            &reference, ContentHash([94; 32]), 0.025,
        ).unwrap();
        // Independent piecewise-linear coefficient integral:
        // 40*(1+3)e-5/2 + 60*(3+2)e-5/2 = 0.0023.
        let free = 0.0023;
        for i in 0..3 { close(expansion.thermal_strain().free_strain_mandel()[i], free); }
        assert_eq!(reference.material_state_identity, source.state().identity().0);
        let (tet_nodes, tet_cells) = promote_p1_to_p2(&TET_NODES, &TETS, None, cx).unwrap();
        let physical = |p: [i8; 3]| p.map(|x| 0.5 * (1.0 + f64::from(x)));
        let hex20_nodes = HEX20_REFERENCE_NODES.map(physical);
        let hex27_nodes = HEX27_REFERENCE_NODES.map(physical);
        let hex20 = [core::array::from_fn(|i| i)];
        let hex27 = [core::array::from_fn(|i| i)];
        for (mesh, nodes) in [
            (SolidMeshRef::Tet4 { nodes_m: &TET_NODES, elements: &TETS }, &TET_NODES[..]),
            (SolidMeshRef::Tet10 {
                nodes_m: &tet_nodes, elements: &tet_cells,
                quadrature: TetQuadrature { stiffness_points_per_dir: 3, mass_points_per_dir: 4 },
            }, &tet_nodes[..]),
            (SolidMeshRef::Hex8 {
                nodes_m: &HEX_NODES, elements: &HEXES, formulation: Formulation::Standard,
                stiffness_points_per_dir: 2,
            }, &HEX_NODES[..]),
            (SolidMeshRef::Hex8 {
                nodes_m: &HEX_NODES, elements: &HEXES, formulation: Formulation::BBar,
                stiffness_points_per_dir: 2,
            }, &HEX_NODES[..]),
            (SolidMeshRef::Hex {
                nodes_m: &hex20_nodes, elements: HexConnectivity::Hex20(&hex20),
                quadrature: HexQuadrature { points_per_dir: 3 },
            }, &hex20_nodes[..]),
            (SolidMeshRef::Hex {
                nodes_m: &hex27_nodes, elements: HexConnectivity::Hex27(&hex27),
                quadrature: HexQuadrature { points_per_dir: 3 },
            }, &hex27_nodes[..]),
        ] {
            // Symmetry supports remove rigid modes and permit isotropic expansion.
            let fixed: Vec<_> = nodes.iter().enumerate().flat_map(|(node, p)| {
                (0..3).filter(move |&d| p[d] == 0.0).map(move |d| 3 * node + d)
            }).collect();
            let solid = prepare(mesh, TetViscoField::Uniform(&material), &fixed, cx);
            let mass = solid.equilibrium().total_mass_kg;
            let operator = solid.quasistatic_operator(0.07, cx).unwrap();
            let mut state = solid.state(cx).unwrap();
            let zero = vec![0.0; nodes.len() * 3];
            let response = state.step_quasistatic_with_thermal_strain(
                &operator, &zero, &zero,
                TetThermalStrainField::Uniform(expansion.thermal_strain()), config(), cx,
            ).unwrap();
            for (&actual, &x) in response.displacement_m.iter().zip(nodes.iter().flatten()) {
                near(actual, free * x, 2e-12);
            }
            assert!(response.constitutive.full_force_n.iter().all(|v| v.abs() < 2e-7));
            assert!(response.constitutive.energy.stored_j < 1e-18);
            assert!(response.constitutive.step_dissipated_j < 1e-18);
            assert_eq!(solid.equilibrium().total_mass_kg, mass);
            let accepted = state.thermal_strains().unwrap().to_vec();
            state.step_internal_force(&response.displacement_m, 0.04, cx).unwrap();
            assert_eq!(state.thermal_strains().unwrap(), accepted.as_slice());

            // Numerically equal coefficients at another source state cannot
            // replace the original coefficient identity underneath history.
            let refreshed = source.at_query_point(&card, &query(400.0)).unwrap();
            let other = TetViscoMaterial::from_source(&refreshed, rotation).unwrap();
            let stale = ThermalExpansionModel::Isotropic.resolve(
                &card, &query(300.0), &query(400.0), MaterialPropertySelection::SingleClaimOnly,
                &other.equilibrium_material().unwrap(), ContentHash([94; 32]), 0.025,
            ).unwrap();
            let before = state.energy();
            assert!(matches!(state.step_internal_force_with_thermal_strain(
                &response.displacement_m, 0.04,
                TetThermalStrainField::Uniform(stale.thermal_strain()), cx,
            ), Err(ViscoElementError::Elastic(TetElasticError::ThermalElasticStateMismatch { .. }))));
            assert_eq!(state.energy(), before);
            assert_eq!(state.thermal_strains().unwrap(), accepted.as_slice());
        }
    });
}

#[test]
fn restrained_full_tensor_expansion_relaxes_with_original_branch_storage_and_heat() {
    with_cx(|cx, _| {
        let mut equilibrium = diagonal([1200.0, 900.0, 1000.0, 700.0, 600.0, 500.0]);
        equilibrium[0][3] = 130.0;
        equilibrium[3][0] = 130.0;
        let branch = outer([2.0, -1.0, 0.5, 3.0, -2.0, 1.0], 80.0);
        let eta = outer([1.0, 0.0, -1.0, 2.0, 0.5, -0.5], 5.0);
        let tau = 0.2;
        let material = TetViscoMaterial {
            law: TensorProny::new_mandel(
                equilibrium, vec![PronyBranch::new(branch, tau).unwrap()], 0.025,
            ).unwrap().with_kelvin_viscosity_mandel(eta).unwrap(),
            density_kg_m3: 1000.0, material_state_identity: [95; 32],
        };
        let solid = tet(&material, &[3, 4, 5], cx);
        let mut state = solid.state(cx).unwrap();
        let free = [0.001, -0.0004, 0.0002, 0.0007, -0.0002, 0.0003];
        let target = thermal(&material, free);
        let h: f64 = 0.08;
        let phi = tau / h * (1.0 - (-h / tau).exp());
        let ce = dot6(free, action(equilibrium, free));
        let cb = dot6(free, action(branch, free));
        let cv = dot6(free, action(eta, free));
        let work = (0.5 * ce + cb * tau / h * (1.0 - phi) + cv / h) / 6.0;
        let stored = 0.5 * (ce + cb * phi * phi) / 6.0;
        let stress: [f64; 6] = core::array::from_fn(|i| {
            -action(equilibrium, free)[i] - phi * action(branch, free)[i]
                - action(eta, free)[i] / h
        });
        let response = state.step_internal_force_with_thermal_strain(
            &[0.0; 12], h, TetThermalStrainField::Uniform(&target), cx,
        ).unwrap();
        close(response.full_force_n[3], stress[0] / 6.0);
        close(response.full_force_n[4], stress[3] / (6.0 * 2.0_f64.sqrt()));
        close(response.full_force_n[5], stress[5] / (6.0 * 2.0_f64.sqrt()));
        close(response.energy.stored_j, stored);
        close(response.step_work_j, work);
        close(response.step_eigenstrain_work_j, work);
        assert_eq!(response.step_mechanical_work_j, 0.0);
        close(response.step_dissipated_j, work - stored);
        let hold: f64 = 0.12;
        let decay = (-hold / tau).exp();
        let held = state.step_internal_force(&[0.0; 12], hold, cx).unwrap();
        let next_stored = 0.5 * (ce + cb * (phi * decay).powi(2)) / 6.0;
        close(held.energy.stored_j, next_stored);
        close(held.step_dissipated_j, stored - next_stored);
        assert_eq!(held.step_eigenstrain_work_j, 0.0);
        assert_eq!(held.step_work_j, 0.0);
        assert_eq!(held.step_mechanical_work_j, 0.0);
        close(state.eigenstrain_work_j(), work);
        assert_eq!(state.thermal_strains().unwrap(), &[target]);
    });
}

#[test]
fn heterogeneous_thermal_creep_distinguishes_force_control_from_end_restraint() {
    with_cx(|cx, _| {
        let mut nodes = HEX_NODES.to_vec();
        nodes.extend([[2.0, 0.0, 0.0], [2.0, 1.0, 0.0], [2.0, 0.0, 1.0], [2.0, 1.0, 1.0]]);
        let cells = [HEXES[0], [1, 8, 9, 2, 5, 10, 11, 6]];
        let mut materials = [
            scalar_material(1200.0, 600.0, 0.05, 4.0, 0.025),
            scalar_material(600.0, 300.0, 0.2, 12.0, 0.025),
        ];
        materials[1].material_state_identity = [96; 32];
        let strains = [
            thermal(&materials[0], [0.001, 0.0, 0.0, 0.0, 0.0, 0.0]),
            thermal(&materials[1], [0.0004, 0.0, 0.0, 0.0, 0.0, 0.0]),
        ];
        let h: f64 = 0.03;
        let tangent = |e: f64, b: f64, tau: f64, eta: f64| {
            e + b * tau / h * (1.0 - (-h / tau).exp()) + eta / h
        };
        let a = tangent(1200.0, 600.0, 0.05, 4.0);
        let b = tangent(600.0, 300.0, 0.2, 12.0);
        for formulation in [Formulation::Standard, Formulation::BBar] {
            for restrained in [false, true] {
                let fixed: Vec<_> = (0..36).filter(|&d| {
                    d % 3 != 0 || nodes[d / 3][0] == 0.0
                        || (restrained && nodes[d / 3][0] == 2.0)
                }).collect();
                let solid = prepare(SolidMeshRef::Hex8 {
                    nodes_m: &nodes, elements: &cells, formulation, stiffness_points_per_dir: 2,
                }, TetViscoField::PerElement(&materials), &fixed, cx);
                let operator = solid.quasistatic_operator(h, cx).unwrap();
                let mut state = solid.state(cx).unwrap();
                let mut force = [0.0; 36];
                let applied = if restrained { 0.0 } else { 0.05 };
                for node in [8, 9, 10, 11] { force[3 * node] = applied / 4.0; }
                let response = state.step_quasistatic_with_thermal_strain(
                    &operator, &force, &[0.0; 36], TetThermalStrainField::PerElement(&strains),
                    config(), cx,
                ).unwrap();
                let (middle, end, stress) = if restrained {
                    ((a * 0.001 - b * 0.0004) / (a + b), 0.0,
                        -a * b * (0.001 + 0.0004) / (a + b))
                } else {
                    (0.001 + applied / a, 0.0014 + applied * (1.0 / a + 1.0 / b), applied)
                };
                for (node, p) in nodes.iter().enumerate() {
                    near(response.displacement_m[3 * node],
                        if p[0] == 1.0 { middle } else if p[0] == 2.0 { end } else { 0.0 }, 2e-13);
                }
                let left_reaction: f64 = [0, 3, 4, 7].iter()
                    .map(|node| response.reaction_force_n[3 * node]).sum();
                close(left_reaction, -stress);
                let right_reaction: f64 = [8, 9, 10, 11].iter()
                    .map(|node| response.reaction_force_n[3 * node]).sum();
                near(right_reaction, stress - applied, 2e-11);
                let work = &response.constitutive;
                close(work.step_mechanical_work_j + work.step_eigenstrain_work_j,
                    work.energy.stored_j + work.step_dissipated_j);
                assert!(work.step_dissipated_j > 0.0);
                assert_eq!(state.thermal_strains().unwrap(), strains.as_slice());
            }
        }
    });
}

#[test]
fn thermal_midpoint_motion_matches_scalar_memory_inertia_and_separate_work_ports() {
    with_cx(|cx, _| {
        let (e, b, tau, eta) = (1200.0, 600.0, 0.07, 8.0);
        let h: f64 = 0.02;
        let material = scalar_material(e, b, tau, eta, 0.025);
        let solid = tet(&material, &[3], cx);
        let prepared = solid.dynamic_operator(h, dynamic_config(), cx).unwrap();
        let mut state = prepared.start(solid.state(cx).unwrap(), &[0.0; 12], cx).unwrap();
        // One free vertex of a unit affine tet: V=1/6, M_11=rho/60.
        let (volume, mass) = (1.0 / 6.0, material.density_kg_m3 / 60.0);
        close(solid.equilibrium().mass.get(0, 0), mass);
        let decay = (-h / tau).exp();
        let phi = tau / h * (1.0 - decay);
        let ramp = b * tau / h * (1.0 - phi);
        let mean_tangent = 0.5 * e + ramp + eta / h;
        let (mut u, mut v, mut free, mut branch_stress) = (0.0, 0.0, 0.0, 0.0);
        let (mut heat, mut external_work, mut thermal_work) = (0.0, 0.0, 0.0);
        for step in 0..24 {
            let next_free = if step < 6 { 0.0004 * (step + 1) as f64 }
                else if step < 12 { 0.0024 }
                else if step < 18 { 0.0024 - 0.0004 * (step - 11) as f64 }
                else { 0.0 };
            let force = if step < 9 { 0.03 } else { 0.0 };
            let de_free = next_free - free;
            let old_mechanical = u - free;
            let du = (force + 2.0 * mass * v / h
                - volume * (e * old_mechanical + phi * branch_stress - mean_tangent * de_free))
                / (2.0 * mass / (h * h) + volume * mean_tangent);
            let de = du - de_free;
            let mean_stress = e * (old_mechanical + 0.5 * de)
                + phi * branch_stress + ramp * de + eta * de / h;
            // Independent integral of z(t)^2/(B*tau), z=a+c exp(-t/tau).
            let a = b * tau * de / h;
            let c = branch_stress - a;
            let step_heat = volume * (a * a * h / (b * tau)
                + 2.0 * a * c * (1.0 - decay) / b
                + c * c * (1.0 - decay * decay) / (2.0 * b) + eta * de * de / h);
            branch_stress = decay * branch_stress + b * phi * de;
            u += du;
            v = 2.0 * du / h - v;
            free = next_free;
            heat += step_heat;
            external_work += force * du;
            thermal_work -= volume * mean_stress * de_free;
            let mut applied = [0.0; 12];
            applied[3] = force;
            let target = thermal(&material, [free, 0.0, 0.0, 0.0, 0.0, 0.0]);
            let response = if de_free == 0.0 {
                state.advance(&applied, cx).unwrap()
            } else {
                state.advance_with_thermal_strain(
                    &applied, TetThermalStrainField::Uniform(&target), cx,
                ).unwrap()
            };
            close(state.material_history().displacement_m()[3], u);
            close(state.velocity_mps()[3], v);
            close(response.constitutive.step_mechanical_work_j, volume * mean_stress * du);
            close(response.constitutive.step_eigenstrain_work_j, -volume * mean_stress * de_free);
            close(response.constitutive.step_dissipated_j, step_heat);
            close(response.constitutive.energy.dissipated_j, heat);
            close(state.external_work_j(), external_work);
            close(state.eigenstrain_work_j(), thermal_work);
            let stored = 0.5 * volume * (e * (u - free).powi(2) + branch_stress.powi(2) / b);
            close(response.constitutive.energy.stored_j, stored);
            close(response.kinetic_energy_j, 0.5 * mass * v * v);
            near(external_work + thermal_work - stored - 0.5 * mass * v * v - heat, 0.0, 2e-13);

            if step == 5 {
                let before_u = state.material_history().displacement_m().to_vec();
                let before_v = state.velocity_mps().to_vec();
                let before_energy = state.material_history().energy();
                let before_thermal = state.material_history().thermal_strains().unwrap().to_vec();
                let before_work = state.eigenstrain_work_j();
                let before_clock = (state.time_s(), state.steps());
                let stale = TetThermalStrainState::try_new(
                    350.0, 300.0, [0.001; 6], [255; 32], ContentHash([93; 32]), 0.025,
                ).unwrap();
                assert!(state.advance_with_thermal_strain(
                    &applied, TetThermalStrainField::Uniform(&stale), cx,
                ).is_err());
                assert_eq!(state.material_history().displacement_m(), before_u);
                assert_eq!(state.velocity_mps(), before_v);
                assert_eq!(state.material_history().energy(), before_energy);
                assert_eq!(state.material_history().thermal_strains().unwrap(), before_thermal);
                assert_eq!(state.eigenstrain_work_j(), before_work);
                assert_eq!((state.time_s(), state.steps()), before_clock);
            }
        }
    });
}

#[test]
fn dynamic_thermal_offset_does_not_admit_a_fictitious_clamped_endpoint() {
    with_cx(|cx, _| {
        let material = scalar_material(1200.0, 0.0, 0.1, 0.0, 0.01);
        let solid = tet(&material, &[3], cx);
        let h = (1.0_f64 / 3.0).sqrt();
        let prepared = solid.dynamic_operator(h, dynamic_config(), cx).unwrap();
        let mut state = prepared.start(solid.state(cx).unwrap(), &[0.0; 12], cx).unwrap();
        let target = thermal(&material, [0.015, 0.0, 0.0, 0.0, 0.0, 0.0]);
        // The clamped offset strain -0.015 is outside the 0.01 material domain.
        // The actual coupled motion has u=0.0075 and mechanical strain=-0.0075.
        assert!(material.law.validate_mandel_strain(&[-0.015, 0.0, 0.0, 0.0, 0.0, 0.0]).is_err());
        let step = state.advance_with_thermal_strain(
            &[0.0; 12], TetThermalStrainField::Uniform(&target), cx,
        ).unwrap();
        close(state.material_history().displacement_m()[3], 0.0075);
        close(step.constitutive.energy.stored_j, 0.005625);
        close(step.kinetic_energy_j, 0.005625);
        close(step.constitutive.step_eigenstrain_work_j, 0.01125);
        assert_eq!(step.constitutive.step_dissipated_j, 0.0);
    });
}

#[test]
fn bbar_thermal_history_checks_raw_elastic_strain_and_rolls_back_every_state() {
    with_cx(|cx, gate| {
        let material = scalar_material(50.0, 20.0, 0.1, 0.2, 0.0123);
        let solid = prepare(SolidMeshRef::Hex8 {
            nodes_m: &HEX_NODES, elements: &HEXES, formulation: Formulation::BBar,
            stiffness_points_per_dir: 2,
        }, TetViscoField::Uniform(&material), &[], cx);
        let amplitude = 0.007;
        let displacement: Vec<_> = HEX_NODES.iter()
            .flat_map(|[x, y, _]| [amplitude * x * (2.0 * y - 1.0), 0.0, 0.0]).collect();
        let target = thermal(&material, [-0.004, -0.004, -0.004, 0.0, 0.0, 0.0]);
        let q = [0.5 - 0.5 / 3.0_f64.sqrt(), 0.5 + 0.5 / 3.0_f64.sqrt()];
        for x in q {
            for y in q {
                let d = amplitude * (2.0 * y - 1.0);
                let shear = 2.0_f64.sqrt() * amplitude * x;
                material.law.validate_mandel_strain(&[d, 0.0, 0.0, shear, 0.0, 0.0]).unwrap();
                material.law.validate_mandel_strain(&[
                    2.0 * d / 3.0 + 0.004, -d / 3.0 + 0.004, -d / 3.0 + 0.004,
                    shear, 0.0, 0.0,
                ]).unwrap();
            }
        }
        // Raw total and projected elastic strain pass. Raw elastic strain at
        // the upper Gauss point exceeds 0.0123; projection must not hide it.
        let d = amplitude / 3.0_f64.sqrt();
        assert!(material.law.validate_mandel_strain(&[
            d + 0.004, 0.004, 0.004, 2.0_f64.sqrt() * amplitude * q[1], 0.0, 0.0,
        ]).is_err());
        let mut state = solid.state(cx).unwrap();
        let good: Vec<_> = displacement.iter().map(|v| 0.05 * v).collect();
        let warm = thermal(&material, [0.0003, 0.0, 0.0, 0.0, 0.0, 0.0]);
        state.step_internal_force_with_thermal_strain(
            &good, 0.03, TetThermalStrainField::Uniform(&warm), cx,
        ).unwrap();
        let mut untouched = state.clone();
        let energy = state.energy();
        let work = state.eigenstrain_work_j();
        assert!(matches!(state.step_internal_force_with_thermal_strain(
            &displacement, 0.03, TetThermalStrainField::Uniform(&target), cx,
        ), Err(ViscoElementError::Material(_))));
        let rebased = TetThermalStrainState::try_new(
            303.0, 299.0, warm.free_strain_mandel(), material.material_state_identity,
            ContentHash([93; 32]), 0.025,
        ).unwrap();
        assert!(state.step_internal_force_with_thermal_strain(
            &good, 0.03, TetThermalStrainField::Uniform(&rebased), cx,
        ).is_err());
        assert!(state.step_internal_force_with_thermal_strain(
            &good, 0.03, TetThermalStrainField::PerElement(&[]), cx,
        ).is_err());
        assert_eq!(state.energy(), energy);
        assert_eq!(state.eigenstrain_work_j(), work);
        assert_eq!(state.displacement_m(), good);
        assert_eq!(state.thermal_strains().unwrap(), &[warm]);
        let actual = state.step_internal_force(&good, 0.03, cx).unwrap();
        let expected = untouched.step_internal_force(&good, 0.03, cx).unwrap();
        assert_eq!(actual.full_force_n, expected.full_force_n);
        assert_eq!(actual.full_mean_force_n, expected.full_mean_force_n);
        assert_eq!(actual.energy, expected.energy);
        assert_eq!(actual.step_dissipated_j, expected.step_dissipated_j);
        let before_cancel = state.energy();
        gate.request();
        assert!(matches!(state.step_internal_force_with_thermal_strain(
            &good, 0.03, TetThermalStrainField::Uniform(&warm), cx,
        ), Err(ViscoElementError::Elastic(TetElasticError::Cancelled))));
        assert_eq!(state.energy(), before_cancel);
        assert_eq!(state.thermal_strains().unwrap(), &[warm]);
    });
}
