//! Original source thermal potentials through native plane continuum solves.
//! Synthetic source data; no calibration or measured-material claim.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::elastic::thermal::{AffineThermoelastic, ThermoelasticSymmetry};
use fs_material::plastic::reduction::ReductionConfig;
use fs_material::state_point::*;
use fs_material::{SmallStrainLaw, tensor};
use fs_qty::{
    Dims, Pressure, QuantitySpec, Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_solid::{
    Formulation, IncrementSettings, Mesh2, PlaneAssemblyBudget, PlaneKind, SmallStrainLawField,
    SmallStrainPlaneProblem,
};

fn temperature_quantity() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new().with_quantity("T", temperature_quantity(), t).unwrap()
}

fn potential(isotropic: bool, axial_only: bool) -> AffineThermoelastic {
    potential_with_expansion_slope(isotropic, axial_only, 1e-8)
}

fn potential_with_expansion_slope(
    isotropic: bool,
    axial_only: bool,
    expansion_slope: f64,
) -> AffineThermoelastic {
    let mut claims = ClaimSet::new();
    let mut insert = |name: &str, dims: Dims, base: f64, slope: f64| {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(280.0, base - 20.0 * slope), (360.0, base + 60.0 * slope)],
                dims,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity("T", temperature_quantity(), 280.0, 360.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic source thermal mechanics".into(),
                license: "test-only".into(),
                artifact: None,
            },
        }).unwrap();
    };
    let inverse_t = Dims([0, 0, 0, -1, 0, 0]);
    let mut elastic_names = Vec::new();
    let mut expansion_names = Vec::new();
    if isotropic {
        insert(YOUNG_MODULUS_PROPERTY, Pressure::DIMS, 12e9, -3e7);
        insert(POISSON_RATIO_PROPERTY, Dims::NONE, 0.23, 0.0);
        insert(LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, inverse_t, 1e-5, expansion_slope);
        elastic_names.extend([
            (YOUNG_MODULUS_PROPERTY, Pressure::DIMS),
            (POISSON_RATIO_PROPERTY, Dims::NONE),
        ]);
        expansion_names.push((LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, inverse_t));
    } else {
        for i in 0..3 {
            for (name, dims, value, slope) in [
                (ORTHOTROPIC_YOUNG_MODULUS_PROPERTIES[i], Pressure::DIMS,
                    [12e9, 5e9, 3e9][i], [-3e7, 1e7, 0.0][i]),
                (ORTHOTROPIC_POISSON_RATIO_PROPERTIES[i], Dims::NONE,
                    if axial_only { 0.0 } else { [0.1, 0.08, 0.15][i] }, 0.0),
                (ORTHOTROPIC_SHEAR_MODULUS_PROPERTIES[i], Pressure::DIMS,
                    [3e9, 2e9, 1e9][i], [0.0, -1e6, 0.0][i]),
            ] {
                insert(name, dims, value, slope);
                elastic_names.push((name, dims));
            }
            let alpha = if axial_only { [1e-5, 0.0, 0.0] } else { [1e-5, -0.4e-5, 2e-5] };
            insert(ORTHOTROPIC_THERMAL_EXPANSION_COEFFICIENT_PROPERTIES[i], inverse_t,
                alpha[i], [expansion_slope, 0.0, 0.0][i]);
            expansion_names.push((ORTHOTROPIC_THERMAL_EXPANSION_COEFFICIENT_PROPERTIES[i], inverse_t));
        }
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic".into(),
            phase: "solid".into(),
            process: "native thermal plane test".into(),
            revision: 0,
        },
        claims,
        vec![],
    ).unwrap();
    let resolve = |names: Vec<(&str, Dims)>| {
        let requirements: Vec<_> = names.into_iter().map(|(name, dims)| {
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::Finite).unwrap()
        }).collect();
        resolve_material_state_point(&card, &point(300.0), &requirements,
            MaterialPropertySelection::SingleClaimOnly).unwrap()
    };
    let elastic = resolve(elastic_names);
    let expansion = resolve(expansion_names);
    let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let (stiffness, expansion_symmetry) = if isotropic {
        (ThermoelasticSymmetry::Isotropic, ThermoelasticSymmetry::Isotropic)
    } else {
        (
            ThermoelasticSymmetry::OrientedOrthotropic {
                principal_to_target: if axial_only { identity } else {
                    tensor::rotation([1.0, 1.0, 2.0].map(|x| x / 6.0_f64.sqrt()), 0.41)
                },
            },
            ThermoelasticSymmetry::OrientedOrthotropic {
                principal_to_target: if axial_only { identity } else {
                    tensor::rotation([0.0, 1.0, 1.0].map(|x| x / 2.0_f64.sqrt()), 0.54)
                },
            },
        )
    };
    AffineThermoelastic::new(
        &card, &elastic, &expansion, stiffness, expansion_symmetry, [280.0, 360.0], 0.02,
    ).unwrap()
}

fn close(a: f64, b: f64, scale: f64) {
    assert!((a - b).abs() <= 2e-9 * scale.max(a.abs()).max(b.abs()), "{a:e} != {b:e}");
}
fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey { seed: 123, kernel_id: 123, tile: 0, iteration: 0 },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ));
    });
}

#[test]
fn g1_frozen_mechanics_retains_original_thermal_potential_and_tangent() {
    for isotropic in [true, false] {
        let potential = potential(isotropic, false);
        for t in [280.0, 300.0, 340.0, 360.0] {
            let source = potential.at_temperature(t).unwrap();
            let law = source.try_mechanical_law().unwrap();
            assert!(std::ptr::eq(law.source(), &source));
            let e = [0.002, -0.0007, 0.0004, 0.0003, -0.0002, 0.0001];
            assert!(law.admits_strain(&e));
            let state = source.state(tensor::to_mandel(&e), 1.0).unwrap();
            let s = law.stress(&e, &());
            let c = law.tangent(&e, &());
            for (a, b) in tensor::to_mandel(&s).into_iter().zip(state.stress_mandel_pa) {
                close(a, b, 1e7);
            }
            for j in 0..6 {
                let h = 1e-6;
                let mut plus = e;
                let mut minus = e;
                plus[j] += h;
                minus[j] -= h;
                let sp = law.stress(&plus, &());
                let sm = law.stress(&minus, &());
                for i in 0..6 {
                    close((sp[i] - sm[i]) / (2.0 * h), c[i][j], 1e10);
                }
                let energy = |x| source.state(tensor::to_mandel(&x), 1.0).unwrap().free_energy_j;
                close((energy(plus) - energy(minus)) / (2.0 * h),
                    s[j] * if j < 3 { 1.0 } else { 2.0 }, 1e7);
            }
        }
    }
    let potential = potential(true, false);
    let hot = potential.at_temperature(360.0).unwrap();
    let law = hot.try_mechanical_law().unwrap();
    assert!(law.admits_strain(&[0.018, 0.018, 0.018, 0.0, 0.0, 0.0]));
    for bad in [
        [0.021, 0.0, 0.0, 0.0, 0.0, 0.0],
        [-0.0195, 0.0, 0.0, 0.0, 0.0, 0.0],
    ] {
        assert!(hot.state(tensor::to_mandel(&bad), 1.0).is_err());
        assert!(!law.admits_strain(&bad));
        assert!(law.stress(&bad, &()).iter().all(|s| s.is_nan()));
        assert!(law.tangent(&bad, &()).iter().flatten().all(|s| s.is_nan()));
    }
}

#[test]
fn g1_native_p1_q1_plane_stress_recovers_free_anisotropic_thermal_expansion() {
    let potential = potential(false, false);
    let hot = potential.at_temperature(350.0).unwrap();
    let law = hot.try_mechanical_law().unwrap();
    let free = tensor::from_mandel(&hot.free_strain_mandel());
    assert!(free[4].abs() > 1e-6 && free[5].abs() > 1e-6,
        "exercise transverse thermal shears");
    with_cx(|cx| {
        for mesh in [Mesh2::triangles(2.0, 1.0, 2, 2), Mesh2::quads(2.0, 1.0, 2, 2)] {
            let right = mesh.nodes.iter().position(|&p| p == [2.0, 0.0]).unwrap();
            let supports = [(0, 0.0), (1, 0.0), (2 * right + 1, 2.0 * free[3])];
            let problem = SmallStrainPlaneProblem {
                mesh: &mesh,
                laws: SmallStrainLawField::Uniform(&law),
                plane: PlaneKind::Stress,
                formulation: Formulation::Standard,
                thickness_m: 0.2,
                reduction: ReductionConfig {
                    max_iterations: 8, stress_tolerance_pa: 1e-6, strain_limit: 0.02,
                },
                rotation_limit_rad: 0.02,
                prescribed_m: &supports,
                nodal_forces_n: &[],
                load_path: &[1.0],
                budget: PlaneAssemblyBudget::default(),
                settings: IncrementSettings::default(),
            };
            let solved = problem.solve(cx).unwrap();
            for (&p, u) in mesh.nodes.iter().zip(&solved.displacement_m) {
                close(u[0], free[0] * p[0] + free[3] * p[1], 0.01);
                close(u[1], free[3] * p[0] + free[1] * p[1], 0.01);
            }
            for (state, stress) in solved.quadrature_states.iter()
                .zip(&solved.quadrature_stress_pa)
            {
                for i in 0..6 {
                    close(state.strain[i], free[i], 0.01);
                    close(stress[i], 0.0, 1e7);
                }
                close(hot.state(tensor::to_mandel(&state.strain), 1.0)
                    .unwrap().free_energy_j, 0.0, 1e-6);
            }
        }
    });
}

#[test]
fn g1_nonuniform_heating_and_dead_load_use_each_original_source_stiffness() {
    let potential = potential(false, true);
    let temperatures = [310.0, 350.0];
    let sources = temperatures.map(|t| potential.at_temperature(t).unwrap());
    let laws: Vec<_> = sources.iter().map(|s| s.try_mechanical_law().unwrap()).collect();
    let mesh = Mesh2::quads(2.0, 1.0, 2, 1);
    let mut supports = Vec::new();
    let mut forces = Vec::new();
    for (node, p) in mesh.nodes.iter().enumerate() {
        supports.push((2 * node + 1, 0.0));
        if p[0] == 0.0 { supports.push((2 * node, 0.0)); }
        if p[0] == 2.0 { forces.push((2 * node, 1e5)); }
    }
    let strain = temperatures.map(|t| {
        let dt = t - 300.0;
        1e-5 * dt + 0.5e-8 * dt * dt + 1e6 / (12e9 - 3e7 * dt)
    });
    with_cx(|cx| {
        let problem = SmallStrainPlaneProblem {
            mesh: &mesh,
            laws: SmallStrainLawField::PerElement(&laws),
            plane: PlaneKind::Stress,
            formulation: Formulation::Standard,
            thickness_m: 0.2,
            reduction: ReductionConfig {
                max_iterations: 8, stress_tolerance_pa: 1e-6, strain_limit: 0.02,
            },
            rotation_limit_rad: 0.02,
            prescribed_m: &supports,
            nodal_forces_n: &forces,
            load_path: &[1.0],
            budget: PlaneAssemblyBudget::default(),
            settings: IncrementSettings::default(),
        };
        let solved = problem.solve(cx).unwrap();
        for (&p, u) in mesh.nodes.iter().zip(&solved.displacement_m) {
            close(u[0], strain[0] * p[0].min(1.0) + strain[1] * (p[0] - 1.0).max(0.0), 0.01);
            close(u[1], 0.0, 0.01);
        }
        for stress in &solved.quadrature_stress_pa {
            close(stress[0], 1e6, 1e6);
            for &s in &stress[1..] { close(s, 0.0, 1e6); }
        }
    });
}

#[test]
fn g1_within_element_thermal_gradient_drives_constrained_plane_and_brick_deformation() {
    use fs_solid::linear3::hex::HEX27_REFERENCE_NODES;
    use fs_solid::plastic3::hex::{HexAssemblyBudget, SmallStrainHexProblem};

    // A transversely restrained strip has the exact Q1 displacement
    // u_x = x(a + b y), u_y = u_z = 0. A linear temperature field gives
    // E(y) = E0 + E1 y and free axial strain a_T + b_T y. Independently
    // minimizing its integrated axial/shear energy determines a and b:
    // [M0 M1; M1 M2 + G L^2/3] [a;b] = [M0 a_T+M1 b_T; M1 a_T+M2 b_T].
    // All integrands have degree <= 3 per coordinate, so native Gauss is exact.
    let length: f64 = 2.0;
    let thickness = 0.2;
    let (mean_rise, gradient) = (20.0, 40.0);
    let (a_t, b_t) = (1e-5 * mean_rise, 1e-5 * gradient);
    let (e0, e1, shear) = (12e9 - 3e7 * mean_rise, -3e7 * gradient, 3e9);
    let (m0, m1, m2) = (e0, e1 / 12.0, e0 / 12.0);
    let bending = m2 - m1 * m1 / m0;
    let b = bending * b_t / (bending + shear * length * length / 3.0);
    let a = a_t + m1 / m0 * (b_t - b);
    let expected_energy = 0.5
        * length
        * thickness
        * (m0 * (a - a_t).powi(2)
            + 2.0 * m1 * (a - a_t) * (b - b_t)
            + m2 * (b - b_t).powi(2)
            + shear * length * length * b * b / 3.0);
    assert!(b.abs() > 5e-5, "the gradient must change the deformation");
    let potential = potential_with_expansion_slope(false, true, 0.0);
    let g = 1.0 / 3.0_f64.sqrt();
    let plane_points = [(-g, -g), (g, -g), (g, g), (-g, g)]
        .map(|(xi, eta)| [0.5 * length * (1.0 + xi), 0.5 * eta]);
    let plane_sources: Vec<_> = plane_points
        .iter()
        .map(|p| {
            potential
                .at_temperature(300.0 + mean_rise + gradient * p[1])
                .unwrap()
        })
        .collect();
    let plane_laws: Vec<_> = plane_sources
        .iter()
        .map(|source| source.try_mechanical_law().unwrap())
        .collect();
    let plane_field = [plane_laws.as_slice()];
    let mut mesh = Mesh2::quads(length, 1.0, 1, 1);
    for p in &mut mesh.nodes {
        p[1] -= 0.5;
    }
    let mut supports = Vec::new();
    for (node, p) in mesh.nodes.iter().enumerate() {
        supports.push((2 * node + 1, 0.0));
        if p[0] == 0.0 {
            supports.push((2 * node, 0.0));
        }
    }
    let expected_stress = |[x, y]: [f64; 2]| {
        [
            (e0 + e1 * y) * (a + b * y - a_t - b_t * y),
            0.0,
            0.0,
            shear * b * x,
            0.0,
            0.0,
        ]
    };
    with_cx(|cx| {
        for plane in [PlaneKind::Stress, PlaneKind::Strain] {
            let problem = SmallStrainPlaneProblem {
                mesh: &mesh,
                laws: SmallStrainLawField::PerQuadrature(&plane_field),
                plane,
                formulation: Formulation::Standard,
                thickness_m: thickness,
                reduction: ReductionConfig {
                    max_iterations: 8,
                    stress_tolerance_pa: 1e-6,
                    strain_limit: 0.02,
                },
                rotation_limit_rad: 0.02,
                prescribed_m: &supports,
                nodal_forces_n: &[],
                load_path: &[1.0],
                budget: PlaneAssemblyBudget::default(),
                settings: IncrementSettings::default(),
            };
            let solved = problem.solve(cx).unwrap();
            for (p, u) in mesh.nodes.iter().zip(&solved.displacement_m) {
                close(u[0], p[0] * (a + b * p[1]), 0.001);
                close(u[1], 0.0, 0.001);
            }
            let mut energy = 0.0;
            for (q, state) in solved.quadrature_states.iter().enumerate() {
                for (actual, expected) in solved.quadrature_stress_pa[q]
                    .iter()
                    .zip(expected_stress(plane_points[q]))
                {
                    close(*actual, expected, 1e6);
                }
                energy += plane_sources[q]
                    .state(tensor::to_mandel(&state.strain), length * thickness / 4.0)
                    .unwrap()
                    .free_energy_j;
            }
            close(energy, expected_energy, 1.0);
            let u = solved.displacement_m.as_flattened();
            let force = problem
                .internal_force(u, &solved.quadrature_states, cx)
                .unwrap();
            for (node, p) in mesh.nodes.iter().enumerate() {
                if p[0] == length {
                    close(force[2 * node], 0.0, 1e6);
                }
            }
            // The pointwise source stiffness must also enter every tangent row.
            let tangent = problem
                .tangent_matrix(u, &solved.quadrature_states, cx)
                .unwrap()
                .to_dense();
            let direction = [0.2, -0.1, -0.3, 0.4, 0.1, 0.2, -0.4, -0.2];
            let h = 1e-6;
            let plus: Vec<_> = u.iter().zip(direction).map(|(v, d)| v + h * d).collect();
            let minus: Vec<_> = u.iter().zip(direction).map(|(v, d)| v - h * d).collect();
            let fp = problem
                .internal_force(&plus, &solved.quadrature_states, cx)
                .unwrap();
            let fm = problem
                .internal_force(&minus, &solved.quadrature_states, cx)
                .unwrap();
            for i in 0..u.len() {
                let kd: f64 = direction
                    .iter()
                    .enumerate()
                    .map(|(j, d)| tangent[i * u.len() + j] * d)
                    .sum();
                close((fp[i] - fm[i]) / (2.0 * h), kd, 1e10);
            }
        }

        let nodes: Vec<_> = HEX27_REFERENCE_NODES[..8]
            .iter()
            .map(|r| {
                [
                    0.5 * length * (1.0 + f64::from(r[0])),
                    0.5 * f64::from(r[1]),
                    0.5 * thickness * (1.0 + f64::from(r[2])),
                ]
            })
            .collect();
        let mut hex_points = Vec::new();
        for xi in [-g, g] {
            for eta in [-g, g] {
                for _ in [-g, g] {
                    hex_points.push([0.5 * length * (1.0 + xi), 0.5 * eta]);
                }
            }
        }
        let sources: Vec<_> = hex_points
            .iter()
            .map(|p| {
                potential
                    .at_temperature(300.0 + mean_rise + gradient * p[1])
                    .unwrap()
            })
            .collect();
        let laws: Vec<_> = sources
            .iter()
            .map(|s| s.try_mechanical_law().unwrap())
            .collect();
        let field = [laws.as_slice()];
        let mut supports = Vec::new();
        for (node, p) in nodes.iter().enumerate() {
            supports.extend([(3 * node + 1, 0.0), (3 * node + 2, 0.0)]);
            if p[0] == 0.0 {
                supports.push((3 * node, 0.0));
            }
        }
        let problem = SmallStrainHexProblem {
            nodes_m: &nodes,
            hexahedra: &[[0, 1, 2, 3, 4, 5, 6, 7]],
            laws: SmallStrainLawField::PerQuadrature(&field),
            formulation: Formulation::Standard,
            strain_limit: 0.02,
            rotation_limit_rad: 0.02,
            prescribed_m: &supports,
            nodal_forces_n: &[],
            load_path: &[1.0],
            budget: HexAssemblyBudget::default(),
            settings: IncrementSettings::default(),
        };
        let solved = problem.solve(cx).unwrap();
        let mut energy = 0.0;
        for (p, u) in nodes.iter().zip(&solved.displacement_m) {
            close(u[0], p[0] * (a + b * p[1]), 0.001);
            close(u[1], 0.0, 0.001);
            close(u[2], 0.0, 0.001);
        }
        for q in 0..8 {
            for (actual, expected) in solved.quadrature_stress_pa[q]
                .iter()
                .zip(expected_stress(hex_points[q]))
            {
                close(*actual, expected, 1e6);
            }
            energy += sources[q]
                .state(
                    tensor::to_mandel(&solved.quadrature_strain[q]),
                    solved.quadrature_volume_m3[q],
                )
                .unwrap()
                .free_energy_j;
        }
        close(energy, expected_energy, 1.0);
    });
}

#[test]
fn g0_original_triangle_and_quad_points_keep_distinct_thermal_stresses_and_exact_coverage() {
    use fs_solid::mesh2::{quad_points, shapes_at};

    let potential = potential_with_expansion_slope(false, true, 0.0);
    with_cx(|cx| {
        for mesh in [
            Mesh2::triangles(1.0, 1.0, 1, 1),
            Mesh2::quads(1.0, 1.0, 1, 1),
        ] {
            let temperatures: Vec<Vec<_>> = mesh
                .elems
                .iter()
                .map(|conn| {
                    quad_points(conn.len())
                        .iter()
                        .map(|&(xi, eta, _)| {
                            let (shape, _, _) = shapes_at(&mesh.nodes, conn, xi, eta);
                            let y: f64 = conn
                                .iter()
                                .zip(shape)
                                .map(|(&node, s)| s * mesh.nodes[node][1])
                                .sum();
                            300.0 + 40.0 * y
                        })
                        .collect()
                })
                .collect();
            let sources: Vec<Vec<_>> = temperatures
                .iter()
                .map(|points| {
                    points
                        .iter()
                        .map(|&t| potential.at_temperature(t).unwrap())
                        .collect()
                })
                .collect();
            let laws: Vec<Vec<_>> = sources
                .iter()
                .map(|points| {
                    points
                        .iter()
                        .map(|s| s.try_mechanical_law().unwrap())
                        .collect()
                })
                .collect();
            let field: Vec<_> = laws.iter().map(Vec::as_slice).collect();
            let supports: Vec<_> = (0..2 * mesh.nodes.len()).map(|d| (d, 0.0)).collect();
            let mut problem = SmallStrainPlaneProblem {
                mesh: &mesh,
                laws: SmallStrainLawField::PerQuadrature(&field),
                plane: PlaneKind::Strain,
                formulation: Formulation::Standard,
                thickness_m: 0.2,
                reduction: ReductionConfig {
                    max_iterations: 8,
                    stress_tolerance_pa: 1e-6,
                    strain_limit: 0.02,
                },
                rotation_limit_rad: 0.02,
                prescribed_m: &supports,
                nodal_forces_n: &[],
                load_path: &[1.0],
                budget: PlaneAssemblyBudget::default(),
                settings: IncrementSettings::default(),
            };
            for formulation in [Formulation::Standard, Formulation::BBar] {
                if formulation == Formulation::BBar && mesh.elems[0].len() == 3 {
                    continue;
                }
                problem.formulation = formulation;
                let solved = problem.solve(cx).unwrap();
                for (stress, &t) in solved
                    .quadrature_stress_pa
                    .iter()
                    .zip(temperatures.iter().flatten())
                {
                    let dt = t - 300.0;
                    close(stress[0], -(12e9 - 3e7 * dt) * 1e-5 * dt, 1e7);
                    for value in &stress[1..] {
                        close(*value, 0.0, 1e7);
                    }
                }
            }
            problem.laws = SmallStrainLawField::PerQuadrature(&field[..field.len() - 1]);
            assert!(
                problem.solve(cx).is_err(),
                "missing elements cannot be broadcast"
            );
            let mut incomplete = field.clone();
            incomplete[0] = &laws[0][..laws[0].len() - 1];
            problem.laws = SmallStrainLawField::PerQuadrature(&incomplete);
            assert!(
                problem.solve(cx).is_err(),
                "missing points cannot be averaged"
            );
        }
    });
}
