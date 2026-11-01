//! G1/G3/G4: native terminals must retain the actual 3D magnetic history.
//! These manufactured materials check the equations, not measured-device fidelity.
#[path = "spatial_induction_winding/resistors.rs"]
mod resistors;
use fs_em::{
    spatial::{
        nonlinear::SpatialNonlinearMaterial,
        vector_potential::{
            VectorMagneticBoundary, VectorMagneticControls, VectorPotentialFace,
            induction::{
                SpatialInductionControls, SpatialInductionLimits,
                nonlinear::{
                    SpatialNonlinearInduction, SpatialNonlinearInductionControls,
                    winding::{
                        SpatialNonlinearWindingCircuit, SpatialNonlinearWindingConfig,
                        SpatialNonlinearWindingError, SpatialNonlinearWindings,
                        SpatialWindingConductor,
                    },
                },
            },
            nonlinear::{
                NonlinearVectorControls, SpatialNonlinearVectorField,
                SpatialNonlinearVectorMagnetostatics,
            },
        },
    },
    winding::nonlinear::{NonlinearWindingControls, NonlinearWindingInitial, NonlinearWindingPort},
};
use fs_evidence::ValidityDomain;
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    conductor::{
        ELECTRICAL_CONDUCTIVITY_DIMS, ELECTRICAL_CONDUCTIVITY_PROPERTY, ResolvedConductorState,
    },
    field_storage::{
        MAGNETIC_PERMEABILITY_DIMS, MAGNETIC_PERMEABILITY_PROPERTY,
        nonlinear::{REVERSIBLE_FIELD_STRENGTH_PROPERTY, ReversibleMagneticMaterial},
        permanent::{FIELD_STRENGTH_DIMS, FLUX_DENSITY_DIMS},
    },
    state_point::{
        MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
        resolve_material_state_point,
    },
};
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_qty::QuantitySpec;
use fs_solver::{CheckedCgConfig, NewtonKrylovConfig};

fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.16e} != {b:.16e}");
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.into_iter().zip(b).map(|(x, y)| x * y).sum()
}
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|k| a[k] - b[k])
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn material(saturating: bool, sigma: f64) -> (SpatialNonlinearMaterial, ResolvedConductorState) {
    let (name, dims, value, interpolation, validity) = if saturating {
        (
            REVERSIBLE_FIELD_STRENGTH_PROPERTY,
            FIELD_STRENGTH_DIMS,
            PropertyValue::Curve {
                abscissa: "B".into(),
                abscissa_dims: FLUX_DENSITY_DIMS,
                knots: vec![(0.0, 0.0), (0.01, 0.005), (1.0, 1.985), (100.0, 397.985)],
                dims: FIELD_STRENGTH_DIMS,
            },
            InterpolationPolicy::LinearInside,
            ValidityDomain::unconstrained().with_quantity(
                "B",
                QuantitySpec::dimensional(FLUX_DENSITY_DIMS),
                0.0,
                100.0,
            ),
        )
    } else {
        (
            MAGNETIC_PERMEABILITY_PROPERTY,
            MAGNETIC_PERMEABILITY_DIMS,
            PropertyValue::Scalar {
                value: 2.0,
                dims: MAGNETIC_PERMEABILITY_DIMS,
            },
            InterpolationPolicy::ConstantWithinValidity,
            ValidityDomain::unconstrained(),
        )
    };
    let mut claims = ClaimSet::new();
    for (name, dims, value, interpolation, validity) in [
        (name, dims, value, interpolation, validity),
        (
            ELECTRICAL_CONDUCTIVITY_PROPERTY,
            ELECTRICAL_CONDUCTIVITY_DIMS,
            PropertyValue::Scalar {
                value: sigma,
                dims: ELECTRICAL_CONDUCTIVITY_DIMS,
            },
            InterpolationPolicy::ConstantWithinValidity,
            ValidityDomain::unconstrained(),
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value,
                validity,
                interpolation,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured closed 3D winding".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic coil/core".into(),
            phase: "solid".into(),
            process: "regression".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let query = if saturating {
        QueryPoint::new()
            .with_quantity("B", QuantitySpec::dimensional(FLUX_DENSITY_DIMS), 0.0)
            .unwrap()
    } else {
        QueryPoint::new()
    };
    let bundle = resolve_material_state_point(
        &card,
        &query,
        &[
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::NonNegative)
                .unwrap(),
            ScalarPropertyRequirement::try_new(
                ELECTRICAL_CONDUCTIVITY_PROPERTY,
                ELECTRICAL_CONDUCTIVITY_DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )
            .unwrap(),
        ],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let magnetic = if saturating {
        SpatialNonlinearMaterial::Reversible(
            ReversibleMagneticMaterial::resolve(
                &card,
                &query,
                "B",
                100.0,
                MaterialPropertySelection::SingleClaimOnly,
            )
            .unwrap(),
        )
    } else {
        SpatialNonlinearMaterial::Linear(bundle.clone())
    };
    (
        magnetic,
        ResolvedConductorState::try_from_material(bundle).unwrap(),
    )
}

fn cube(material: SpatialNonlinearMaterial) -> SpatialNonlinearVectorMagnetostatics {
    let mut positions = vec![];
    for k in 0..=2 {
        for j in 0..=2 {
            for i in 0..=2 {
                positions.push([f64::from(i) / 2.0, f64::from(j) / 2.0, f64::from(k) / 2.0]);
            }
        }
    }
    let id = |i, j, k| ((k * 3 + j) * 3 + i) as u32;
    let mut tets = vec![];
    for k in 0..2 {
        for j in 0..2 {
            for i in 0..2 {
                let (a, b, c, d) = (
                    id(i, j, k),
                    id(i + 1, j, k),
                    id(i, j + 1, k),
                    id(i + 1, j + 1, k),
                );
                let (e, f, g, h) = (
                    id(i, j, k + 1),
                    id(i + 1, j, k + 1),
                    id(i, j + 1, k + 1),
                    id(i + 1, j + 1, k + 1),
                );
                tets.extend([
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
    let count = tets.len();
    SpatialNonlinearVectorMagnetostatics::new(
        positions,
        tets,
        vec![material; count],
        WeightedAssemblyLimits {
            max_cells: 1000,
            max_dofs: 4000,
            max_triplets: 100_000,
        },
        &mut || false,
    )
    .unwrap()
}

// J_unit=2 curl(N_center e_z) has a complete closed return and carries -1 A
// through x=1/2,y>1/2: -2 integral N(1/2,1/2,z) dz = -1.
// The second profile adds a separately normalized x-axis loop in series.
fn profiles(p: &SpatialNonlinearVectorMagnetostatics) -> Vec<Vec<[f64; 3]>> {
    let center = p.positions_m().iter().position(|&x| x == [0.5; 3]).unwrap() as u32;
    let mut result = vec![vec![], vec![]];
    for tet in &p.mesh().tets {
        let gradient = if tet.contains(&center) {
            let face: Vec<_> = tet.iter().copied().filter(|&v| v != center).collect();
            let point = p.positions_m()[face[0] as usize];
            let normal = cross(
                sub(p.positions_m()[face[1] as usize], point),
                sub(p.positions_m()[face[2] as usize], point),
            );
            let height = dot(normal, sub(p.positions_m()[center as usize], point));
            normal.map(|v| v / height)
        } else {
            [0.0; 3]
        };
        result[0].push(cross(gradient, [0.0, 0.0, 2.0]));
        result[1].push(cross(gradient, [2.0, 0.0, 2.0]));
    }
    result
}
fn current(profiles: &[Vec<[f64; 3]>], currents: &[f64]) -> Vec<[f64; 3]> {
    (0..profiles[0].len())
        .map(|cell| {
            std::array::from_fn(|k| {
                profiles
                    .iter()
                    .zip(currents)
                    .map(|(p, i)| p[cell][k] * i)
                    .sum()
            })
        })
        .collect()
}
fn static_controls() -> NonlinearVectorControls {
    let cg = CheckedCgConfig {
        tolerance: 1e-12,
        max_iterations: 3000,
        max_corrections: 2,
    };
    NonlinearVectorControls {
        field: VectorMagneticControls {
            field: cg,
            gauge: cg,
            maximum_gauge_product_terms: 2_000_000,
            current_tolerance_a: 2e-8,
            flux_tolerance_wb: 2e-9,
            energy_tolerance_j: 2e-8,
        },
        relative_residual_tolerance: 2e-11,
        maximum_newton_iterations: 50,
        maximum_backtracks: 35,
    }
}
fn controls() -> NonlinearWindingControls<SpatialNonlinearInductionControls> {
    NonlinearWindingControls {
        field: SpatialNonlinearInductionControls {
            induction: SpatialInductionControls {
                relative_tolerance: 2e-11,
                absolute_residual_a: 1e-11,
                restart: 64,
                maximum_cycles: 100,
                maximum_krylov_scalars: 2_000_000,
                current_tolerance_a: 2e-8,
                flux_tolerance_wb: 2e-9,
                energy_tolerance_j: 2e-8,
                maximum_cell_magnetic_peclet: 100.0,
            },
            maximum_newton_iterations: 50,
            maximum_backtracks: 35,
        },
        maximum_field_trials: 200,
        maximum_newton_iterations: 25,
        difference_step_a: 1e-5,
        terminal_voltage_tolerance_v: 2e-7,
        power_tolerance_w: 2e-7,
        energy_tolerance_j: 3e-7,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-9,
            relative_tolerance: 1e-12,
            forcing_minimum: 1e-4,
            linear_restart: 8,
            ..NewtonKrylovConfig::default()
        },
    }
}
fn config() -> SpatialNonlinearWindingConfig {
    SpatialNonlinearWindingConfig {
        max_windings: 4,
        max_profile_values: 100_000,
        current_balance_tolerance_turns: 1e-10,
    }
}
fn solve<'m>(
    p: &'m SpatialNonlinearVectorMagnetostatics,
    j: &[[f64; 3]],
) -> SpatialNonlinearVectorField<'m> {
    let boundary = VectorMagneticBoundary {
        potential: p
            .exterior_faces()
            .into_iter()
            .map(|face| VectorPotentialFace {
                face,
                edge_integrals_wb: [0.0; 3],
            })
            .collect(),
    };
    p.solve(
        j,
        &boundary,
        &vec![0.0; p.mesh().edges.len()],
        static_controls(),
        &mut || false,
    )
    .unwrap()
}
fn history<'m>(
    p: &'m SpatialNonlinearVectorMagnetostatics,
    j: &[[f64; 3]],
    eddy: Vec<Option<ResolvedConductorState>>,
) -> SpatialNonlinearInduction<'m> {
    SpatialNonlinearInduction::new(
        &solve(p, j),
        eddy,
        SpatialInductionLimits {
            maximum_unknowns: 20000,
            maximum_conducting_cells: 1000,
            maximum_assembly_terms: 1_000_000,
        },
        2e-9,
        &mut || false,
    )
    .unwrap()
}
fn sources(
    profiles: &[Vec<[f64; 3]>],
    conductor: &ResolvedConductorState,
    fill: f64,
) -> Vec<SpatialWindingConductor> {
    profiles
        .iter()
        .enumerate()
        .flat_map(|(winding, profile)| {
            profile
                .iter()
                .enumerate()
                .filter(|(_, p)| **p != [0.0; 3])
                .map(move |(cell, _)| SpatialWindingConductor {
                    winding,
                    cell,
                    fill_fraction: fill,
                    conductor: conductor.clone(),
                })
        })
        .collect()
}
fn resistance(profile: &[[f64; 3]], sigma: f64, fill: f64) -> f64 {
    // Every tetrahedron has volume 1/48; integrate |Junit|²/(sigma*fill).
    profile
        .iter()
        .map(|&j| dot(j, j) / (48.0 * sigma * fill))
        .sum()
}
fn initial(currents: Vec<f64>, external_inputs: Vec<f64>) -> NonlinearWindingInitial {
    NonlinearWindingInitial {
        currents_a: currents,
        remaining_fluxes_wb: vec![],
        charges_c: vec![],
        external_inputs,
    }
}
fn pair(p: &SpatialNonlinearVectorMagnetostatics, a: &[f64], profile: &[[f64; 3]]) -> f64 {
    // Independent exact integral of Whitney Nedelec functions against cellwise
    // constant J: int W_ij = volume/4*(grad N_j-grad N_i).
    let mut integral = 0.0;
    for (cell, tet) in p.mesh().tets.iter().enumerate() {
        let gradients: Vec<_> = tet
            .iter()
            .map(|&v| {
                let face: Vec<_> = tet.iter().copied().filter(|&u| u != v).collect();
                let point = p.positions_m()[face[0] as usize];
                let normal = cross(
                    sub(p.positions_m()[face[1] as usize], point),
                    sub(p.positions_m()[face[2] as usize], point),
                );
                let height = dot(normal, sub(p.positions_m()[v as usize], point));
                normal.map(|x| x / height)
            })
            .collect();
        for i in 0..4 {
            for j in i + 1..4 {
                let mut vertices = [tet[i], tet[j]];
                let sign = if vertices[0] < vertices[1] { 1.0 } else { -1.0 };
                vertices.sort();
                let edge = p.mesh().edges.binary_search(&vertices).unwrap();
                integral +=
                    sign * a[edge] * dot(sub(gradients[j], gradients[i]), profile[cell]) / 192.0;
            }
        }
    }
    integral
}

#[test]
fn closed_3d_windings_match_independent_two_loop_backward_euler_equations() {
    let (magnetic, conductor) = material(false, 18.0);
    let p = cube(magnetic);
    let profiles = profiles(&p);
    // Static source solves are independent of the field/circuit interface and
    // its linkage implementation. Energy polarization gives the complete L.
    let e0 = solve(&p, &profiles[0]).stored_energy_j();
    let e1 = solve(&p, &profiles[1]).stored_energy_j();
    let es = solve(&p, &current(&profiles, &[1.0, 1.0])).stored_energy_j();
    let ends = [0.02, 0.03];
    let l = [
        2.0 * e0 + ends[0],
        es - e0 - e1,
        es - e0 - e1,
        2.0 * e1 + ends[1],
    ];
    assert!(l[1].abs() > 1e-5);
    let rw = [
        resistance(&profiles[0], 18.0, 0.5),
        resistance(&profiles[1], 18.0, 0.5),
    ];
    let graph = CircuitGraph {
        node_count: 6,
        branches: vec![
            (1, 0, Branch::VoltageSource { port: 0 }),
            (1, 2, Branch::Resistor { ohms: 0.7 }),
            (2, 0, Branch::Inductor { henries: f64::NAN }),
            (3, 0, Branch::VoltageSource { port: 2 }),
            (3, 4, Branch::Resistor { ohms: 1.1 }),
            (4, 0, Branch::Inductor { henries: f64::NAN }),
            (5, 0, Branch::Inductor { henries: 0.4 }),
            (1, 5, Branch::Resistor { ohms: 0.6 }),
            (5, 0, Branch::Capacitor { farads: 0.3 }),
        ],
        transformers: vec![],
    };
    let mut start = initial(vec![0.0, 0.0], vec![0.0, 17.0, 0.0]);
    start.remaining_fluxes_wb = vec![0.08];
    start.charges_c = vec![0.06];
    let original = SpatialNonlinearWindingCircuit::new(
        history(&p, &current(&profiles, &[0.0, 0.0]), vec![None; 48]),
        graph,
        SpatialNonlinearWindings {
            branch_indices: vec![5, 2],
            profiles_per_m2: profiles.clone(),
            end_inductance_h: ends.to_vec(),
        },
        sources(&profiles, &conductor, 0.5),
        start,
        config(),
        &mut || false,
    )
    .unwrap();
    for i in 0..2 {
        near(original.resistance_ohm()[i], rw[i], 1e-12);
    }
    let mut accepted = original.clone();
    let mut old = [0.0, 0.0];
    let (mut old_phi, mut old_q) = (0.08, 0.06);
    let dt = 0.04;
    for index in 1..=3 {
        let r = [rw[0] + 1.1, rw[1] + 0.7];
        let v = [-0.6, 1.3];
        let (a, d) = (l[0] + dt * r[0], l[3] + dt * r[1]);
        let b = [
            l[0] * old[0] + l[1] * old[1] + dt * v[0],
            l[2] * old[0] + l[3] * old[1] + dt * v[1],
        ];
        let determinant = a * d - l[1] * l[2];
        let expected = [
            (d * b[0] - l[1] * b[1]) / determinant,
            (a * b[1] - l[2] * b[0]) / determinant,
        ];
        // The unrelated native LC storage keeps its trapezoidal law. Its
        // resistor is fed by the ideal source, so it cannot alter the two
        // winding equations above. Solve its two scalar ODEs independently.
        let (lc_a, lc_b, lc_d) = (-dt / 0.6, dt / 0.8, 1.0 + dt / 0.36);
        let lc_rhs = [
            old_phi + dt * old_q / 0.6,
            old_q + dt * 1.3 / 0.6 - dt * old_q / 0.36 - dt * old_phi / 0.8,
        ];
        let lc_determinant = lc_d - lc_a * lc_b;
        let phi = (lc_d * lc_rhs[0] - lc_a * lc_rhs[1]) / lc_determinant;
        let q = (lc_rhs[1] - lc_b * lc_rhs[0]) / lc_determinant;
        let lc_voltage = (q + old_q) / 0.6;
        let lc_heat = dt * (1.3 - lc_voltage).powi(2) / 0.6;
        let lc_supply = dt * 1.3 * (1.3 - lc_voltage) / 0.6;
        let lc_before = 0.5 * (old_phi * old_phi / 0.4 + old_q * old_q / 0.3);
        let lc_after = 0.5 * (phi * phi / 0.4 + q * q / 0.3);
        let candidate = accepted
            .step_to_candidate(
                f64::from(index) * dt,
                &[[0.0; 3]; 27],
                &[1.3, 17.0, -0.6],
                controls(),
                &mut || false,
            )
            .unwrap();
        let next = candidate.circuit();
        let step = candidate.step();
        near(next.state()[next.dae().flux_index[0]], phi, 3e-7);
        near(next.state()[next.dae().charge_index[0]], q, 3e-7);
        for i in 0..2 {
            near(next.currents_a()[i], expected[i], 3e-7);
            near(
                next.flux_linkages_wb()[i],
                l[2 * i] * expected[0] + l[2 * i + 1] * expected[1],
                3e-7,
            );
            near(
                step.terminal_voltage_v[i],
                v[i] - [1.1, 0.7][i] * expected[i],
                3e-7,
            );
            near(
                step.terminal_residual_v[i],
                0.0,
                controls().terminal_voltage_tolerance_v,
            );
        }
        let energy = |i: [f64; 2]| {
            0.5 * (l[0] * i[0] * i[0] + 2.0 * l[1] * i[0] * i[1] + l[3] * i[1] * i[1])
        };
        near(step.stored_before_j, energy(old) + lc_before, 3e-7);
        near(step.stored_after_j, energy(expected) + lc_after, 3e-7);
        near(
            step.heat_j,
            dt * (r[0] * expected[0].powi(2) + r[1] * expected[1].powi(2)) + lc_heat,
            3e-7,
        );
        near(
            step.numerical_dissipation_j,
            energy([expected[0] - old[0], expected[1] - old[1]]),
            3e-7,
        );
        near(
            step.supplied_work_j,
            dt * (v[0] * expected[0] + v[1] * expected[1]) + lc_supply,
            3e-7,
        );
        near(
            step.native.supplied,
            step.native.dissipated + lc_after - lc_before,
            3e-7,
        );
        near(step.field.heat_j, 0.0, 1e-14);
        near(step.energy_residual_j, 0.0, controls().energy_tolerance_j);
        old = expected;
        old_phi = phi;
        old_q = q;
        accepted = candidate.into_parts().0;
    }
    assert_eq!(original.currents_a(), [0.0, 0.0]);
    assert_eq!(original.time_s(), 0.0);
}

fn driven<'m>(
    p: &'m SpatialNonlinearVectorMagnetostatics,
    profile: Vec<[f64; 3]>,
    conductor: &ResolvedConductorState,
    eddy: bool,
) -> SpatialNonlinearWindingCircuit<'m> {
    let conductors = profile
        .iter()
        .map(|&j| (eddy && j == [0.0; 3]).then(|| conductor.clone()))
        .collect();
    let profiles = vec![profile];
    SpatialNonlinearWindingCircuit::new_with_ports(
        history(p, &current(&profiles, &[0.0]), conductors),
        CircuitGraph {
            node_count: 2,
            branches: vec![
                (1, 0, Branch::CurrentSource { port: 0 }),
                (1, 0, Branch::Inductor { henries: f64::NAN }),
            ],
            transformers: vec![],
        },
        SpatialNonlinearWindings {
            branch_indices: vec![1],
            profiles_per_m2: profiles.clone(),
            end_inductance_h: vec![0.03],
        },
        vec![NonlinearWindingPort::VoltageInput {
            initial_voltage_v: 0.0,
            residual_scale_ohm: 1.0,
            current_tolerance_a: 2e-7,
        }],
        sources(&profiles, conductor, 0.6),
        initial(vec![0.0], vec![0.0]),
        config(),
        &mut || false,
    )
    .unwrap()
}

fn psi(b: [f64; 3]) -> f64 {
    let b = dot(b, b).sqrt();
    if b <= 0.01 {
        0.25 * b * b
    } else if b <= 1.0 {
        0.000025 + 0.005 * (b - 0.01) + (b - 0.01).powi(2)
    } else {
        0.985075 + 1.985 * (b - 1.0) + 2.0 * (b - 1.0).powi(2)
    }
}
fn field_strength(b: [f64; 3]) -> [f64; 3] {
    let norm = dot(b, b).sqrt();
    if norm == 0.0 {
        return [0.0; 3];
    }
    let h = if norm <= 0.01 {
        0.5 * norm
    } else if norm <= 1.0 {
        0.005 + 2.0 * (norm - 0.01)
    } else {
        1.985 + 4.0 * (norm - 1.0)
    };
    b.map(|v| v * h / norm)
}

#[test]
fn native_current_drive_uses_saturated_linkage_and_separates_physical_and_numerical_heat() {
    let (magnetic, conductor) = material(true, 18.0);
    let p = cube(magnetic);
    let profile = profiles(&p).remove(0);
    let original = driven(&p, profile.clone(), &conductor, false);
    let rw = resistance(&profile, 18.0, 0.6);
    let mut accepted = original.clone();
    let mut old_current = 0.0;
    let mut old_b = vec![[0.0; 3]; 48];
    let mut old_linkage = 0.0;
    for (target, i) in [(0.1, 1.0), (0.18, 0.4)] {
        let dt = target - accepted.time_s();
        let direct = solve(
            &p,
            &profile.iter().map(|j| j.map(|v| v * i)).collect::<Vec<_>>(),
        );
        let linkage = pair(&p, direct.potential_wb(), &profile);
        let stored_before =
            old_b.iter().copied().map(psi).sum::<f64>() / 48.0 + 0.015 * old_current * old_current;
        let stored_after = direct
            .flux_density_t()
            .iter()
            .copied()
            .map(psi)
            .sum::<f64>()
            / 48.0
            + 0.015 * i * i;
        let numerical = old_b
            .iter()
            .zip(direct.flux_density_t())
            .map(|(&old, &next)| dot(field_strength(next), sub(next, old)) - (psi(next) - psi(old)))
            .sum::<f64>()
            / 48.0
            + 0.015 * (i - old_current).powi(2);
        let voltage = rw * i + (linkage - old_linkage + 0.03 * (i - old_current)) / dt;
        let candidate = accepted
            .step_to_candidate(target, &[[0.0; 3]; 27], &[i], controls(), &mut || false)
            .unwrap();
        let step = candidate.step();
        near(candidate.circuit().currents_a()[0], i, 2e-7);
        near(
            candidate.circuit().flux_linkages_wb()[0],
            linkage + 0.03 * i,
            3e-7,
        );
        near(step.terminal_voltage_v[0], voltage, 3e-7);
        near(step.native.y[0], voltage, 3e-7);
        near(step.native.y[1], -i, 2e-7);
        near(step.native.supplied, 0.0, 3e-7);
        near(step.stored_before_j, stored_before, 3e-7);
        near(step.stored_after_j, stored_after, 3e-7);
        near(step.heat_j, dt * rw * i * i, 3e-7);
        near(step.numerical_dissipation_j, numerical, 3e-7);
        near(step.supplied_work_j, dt * i * voltage, 3e-7);
        near(step.energy_residual_j, 0.0, controls().energy_tolerance_j);
        assert!(numerical > 1e-5);
        assert!(
            direct
                .flux_density_t()
                .iter()
                .any(|&b| dot(b, b).sqrt() > 0.01)
        );
        for cell in &step.cells {
            let j = profile[cell.cell].map(|v| v * i / 0.6);
            near(cell.conductor_volume_m3, 0.6 / 48.0, 1e-14);
            for k in 0..3 {
                near(cell.current_density_a_m2[k], j[k], 2e-6);
                near(cell.electric_field_v_m[k], j[k] / 18.0, 2e-7);
            }
            near(cell.heat_j, dt * dot(j, j) * 0.6 / (48.0 * 18.0), 3e-7);
        }
        old_current = i;
        old_b = direct.flux_density_t().to_vec();
        old_linkage = linkage;
        accepted = candidate.into_parts().0;
    }
    assert_eq!(original.time_s(), 0.0);
    assert_eq!(original.currents_a(), [0.0]);
}

#[test]
fn native_terminals_preserve_saturated_diffusion_history_and_exact_spatial_joule_work() {
    let (magnetic, conductor) = material(true, 18.0);
    let p = cube(magnetic);
    let profile = profiles(&p).remove(0);
    let original = driven(&p, profile.clone(), &conductor, true);
    let mut accepted = original.clone();
    let rw = resistance(&profile, 18.0, 0.6);
    for (target, i) in [(0.05, 1.0), (0.12, 0.4)] {
        let dt = target - accepted.time_s();
        let previous = accepted.currents_a()[0];
        let old_a = accepted.field().potential_wb().to_vec();
        let old_b = accepted
            .field()
            .field(&mut || false)
            .unwrap()
            .flux_density_t()
            .to_vec();
        // Independently drive the original A-phi owner with the declared ideal
        // current. This checks that terminal trials use its accepted history,
        // not a restarted magnetostatic field or a fitted inductance.
        let direct = accepted
            .field()
            .step_candidate(
                &profile.iter().map(|j| j.map(|v| v * i)).collect::<Vec<_>>(),
                &[[0.0; 3]; 27],
                dt,
                controls().field,
                &mut || false,
            )
            .unwrap();
        let candidate = accepted
            .step_to_candidate(target, &[[0.0; 3]; 27], &[i], controls(), &mut || false)
            .unwrap();
        let step = candidate.step();
        for (&a, &b) in candidate
            .circuit()
            .field()
            .potential_wb()
            .iter()
            .zip(direct.state().potential_wb())
        {
            near(a, b, 3e-7);
        }
        let delta: Vec<_> = direct
            .state()
            .potential_wb()
            .iter()
            .zip(&old_a)
            .map(|(a, b)| a - b)
            .collect();
        let voltage = rw * i + (pair(&p, &delta, &profile) + 0.03 * (i - previous)) / dt;
        near(step.terminal_voltage_v[0], voltage, 3e-7);
        near(
            step.field.impressed_work_j,
            i * pair(&p, &delta, &profile),
            3e-7,
        );
        let mut eddy_heat = 0.0;
        for cell in &step.field.cells {
            assert_eq!(profile[cell.cell], [0.0; 3]);
            let e = cell.conductor_electric_field_v_m;
            let j = cell.current_density_a_m2;
            let se = std::array::from_fn(|k| e.iter().map(|v| v[k]).sum());
            let sj = std::array::from_fn(|k| j.iter().map(|v| v[k]).sum());
            let diagonal = (0..4).map(|k| dot(e[k], j[k])).sum::<f64>();
            let heat = dt * cell.volume_m3 / 20.0 * (dot(se, sj) + diagonal);
            near(cell.heat_j, heat, 2e-9);
            eddy_heat += heat;
        }
        assert!(eddy_heat > 1e-8);
        near(step.field.heat_j, eddy_heat, 2e-8);
        near(step.heat_j, eddy_heat + dt * rw * i * i, 3e-7);
        let next_b = &step.field.endpoint_flux_density_t;
        let magnetic_before = old_b.iter().copied().map(psi).sum::<f64>() / 48.0;
        let magnetic_after = next_b.iter().copied().map(psi).sum::<f64>() / 48.0;
        let numerical = old_b
            .iter()
            .zip(next_b)
            .map(|(&old, &next)| dot(field_strength(next), sub(next, old)) - (psi(next) - psi(old)))
            .sum::<f64>()
            / 48.0
            + 0.015 * (i - previous).powi(2);
        near(
            step.stored_before_j,
            magnetic_before + 0.015 * previous * previous,
            3e-7,
        );
        near(step.stored_after_j, magnetic_after + 0.015 * i * i, 3e-7);
        near(step.numerical_dissipation_j, numerical, 3e-7);
        near(step.supplied_work_j, dt * i * voltage, 3e-7);
        near(step.energy_residual_j, 0.0, controls().energy_tolerance_j);
        assert!(
            step.field.maximum_charge_residual_a <= controls().field.induction.current_tolerance_a
        );
        assert!(
            step.field.maximum_ampere_residual_a <= controls().field.induction.current_tolerance_a
        );
        accepted = candidate.into_parts().0;
    }
    let static_field = solve(
        &p,
        &profile
            .iter()
            .map(|j| j.map(|v| 0.4 * v))
            .collect::<Vec<_>>(),
    );
    let history_difference = accepted
        .field()
        .potential_wb()
        .iter()
        .zip(static_field.potential_wb())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0_f64, f64::max);
    assert!(history_difference > 1e-6);
    assert_eq!(original.time_s(), 0.0);
    assert!(original.field().potential_wb().iter().all(|&a| a == 0.0));
}

#[test]
fn terminal_budget_and_late_cancellation_leave_both_histories_unchanged() {
    let (magnetic, conductor) = material(true, 18.0);
    let p = cube(magnetic);
    let profile = profiles(&p).remove(0);
    let original = driven(&p, profile, &conductor, true);
    let old_a = original.field().potential_wb().to_vec();
    let old_native = original.state().to_vec();
    let mut limited = controls();
    limited.maximum_field_trials = 1;
    assert!(matches!(
        original.step_to_candidate(0.1, &[[0.0; 3]; 27], &[0.7], limited, &mut || false),
        Err(SpatialNonlinearWindingError::TrialBudget)
    ));
    let mut calls = 0;
    let success = original
        .step_to_candidate(0.1, &[[0.0; 3]; 27], &[0.7], controls(), &mut || {
            calls += 1;
            false
        })
        .unwrap();
    assert!(calls > 20);
    let mut counter = 0;
    assert!(matches!(
        original.step_to_candidate(0.1, &[[0.0; 3]; 27], &[0.7], controls(), &mut || {
            counter += 1;
            counter == calls - 1
        }),
        Err(SpatialNonlinearWindingError::Cancelled)
    ));
    assert_eq!(original.field().potential_wb(), old_a);
    assert_eq!(original.state(), old_native);
    assert_eq!(original.currents_a(), [0.0]);
    assert_eq!(original.time_s(), 0.0);
    let retry = original
        .step_to_candidate(0.1, &[[0.0; 3]; 27], &[0.7], controls(), &mut || false)
        .unwrap();
    assert_eq!(success.circuit().state(), retry.circuit().state());
    assert_eq!(
        success.circuit().field().potential_wb(),
        retry.circuit().field().potential_wb()
    );
    assert_eq!(
        success.step().heat_j.to_bits(),
        retry.step().heat_j.to_bits()
    );
    assert_eq!(
        success.step().supplied_work_j.to_bits(),
        retry.step().supplied_work_j.to_bits()
    );
}
