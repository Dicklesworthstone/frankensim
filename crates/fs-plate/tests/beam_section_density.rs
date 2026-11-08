//! Independent volume integrals for heterogeneous original beam sections.
use fs_plate::{
    PlateSection,
    shell::{
        ReferenceDensityMoments, ShellMesh, ShellSupport,
        junction::{
            BeamOnShell, ShellComponent, TieBudget, TieSide,
            tie_beams_to_shells_with_reference_inertia,
            tie_prestressed_beams_to_shells_with_reference_inertia, tie_shells,
        },
        stiffened::{
            BeamMass, BeamSection, BeamSectionInertia, ShellBeam,
            assemble_stiffened_shell_with_reference_inertia,
            assemble_stiffened_shell_with_reference_inertia_and_body_acceleration,
            timoshenko::{
                BeamKind, ShearAreas, assemble_beam_chain_with_reference_inertia,
                assemble_prestressed_beam_chain_with_reference_inertia,
                segment_matrices_with_reference_density, segment_matrices_with_reference_inertia,
            },
        },
    },
};
use fs_sparse::Csr;

#[derive(Clone, Copy)]
struct Patch {
    y: [f64; 2],
    z: [f64; 2],
    rho: [f64; 2],
}
fn patches() -> [Patch; 4] {
    [
        Patch {
            y: [-0.04, 0.0],
            z: [-0.02, 0.0],
            rho: [100.0, 700.0],
        },
        Patch {
            y: [0.0, 0.04],
            z: [-0.02, 0.0],
            rho: [900.0, 200.0],
        },
        Patch {
            y: [-0.04, 0.0],
            z: [0.0, 0.02],
            rho: [500.0, 1200.0],
        },
        Patch {
            y: [0.0, 0.04],
            z: [0.0, 0.02],
            rho: [1800.0, 300.0],
        },
    ]
}
fn section() -> BeamSection {
    let area = 0.08 * 0.04;
    BeamSection {
        young_pa: 12e9,
        shear_pa: 0.8e9,
        density_kg_m3: 450.0,
        area_m2: area,
        iy_m4: area * 0.04 * 0.04 / 12.0,
        iz_m4: area * 0.08 * 0.08 / 12.0,
        torsion_m4: 1e-6,
    }
}
fn kinds() -> [BeamKind; 2] {
    [
        BeamKind::EulerBernoulli,
        BeamKind::Timoshenko(ShearAreas {
            y_m2: 0.0027,
            z_m2: 0.0023,
        }),
    ]
}
fn moment(patches: &[Patch], t: f64) -> BeamSectionInertia {
    let integral = |[a, b]: [f64; 2], power: i32| {
        (b.powi(power + 1) - a.powi(power + 1)) / f64::from(power + 1)
    };
    let q: [f64; 6] = core::array::from_fn(|k| {
        let [py, pz] = [[0, 0], [1, 0], [0, 1], [2, 0], [1, 1], [0, 2]][k];
        patches
            .iter()
            .map(|p| ((1.0 - t) * p.rho[0] + t * p.rho[1]) * integral(p.y, py) * integral(p.z, pz))
            .sum()
    });
    BeamSectionInertia::new(q[0], q[1], q[2], q[3], q[4], q[5]).unwrap()
}
fn reversed_y(q: BeamSectionInertia) -> BeamSectionInertia {
    BeamSectionInertia::new(
        q.line_mass_kg_per_m(),
        -q.first_y_kg(),
        q.first_z_kg(),
        q.second_yy_kg_m(),
        -q.second_yz_kg_m(),
        q.second_zz_kg_m(),
    )
    .unwrap()
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    core::array::from_fn(|i| a[i] + b[i])
}
fn rigid_velocity(p: [f64; 3], coordinate: usize) -> [f64; 3] {
    let mut e = [0.0; 3];
    e[coordinate % 3] = 1.0;
    if coordinate < 3 { e } else { cross(e, p) }
}
fn rigid(nodes: &[[f64; 3]], coordinate: usize) -> Vec<f64> {
    nodes
        .iter()
        .flat_map(|&p| {
            let v = rigid_velocity(p, coordinate);
            let mut w = [0.0; 3];
            if coordinate >= 3 {
                w[coordinate - 3] = 1.0;
            }
            [v[0], v[1], v[2], w[0], w[1], w[2]]
        })
        .collect()
}
fn bilinear(m: &[f64], a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .enumerate()
        .map(|(i, x)| {
            x * m[i * b.len()..(i + 1) * b.len()]
                .iter()
                .zip(b)
                .map(|(m, y)| m * y)
                .sum::<f64>()
        })
        .sum()
}
fn sparse(m: &Csr, a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .enumerate()
        .map(|(i, x)| {
            let (columns, values) = m.row(i);
            x * columns
                .iter()
                .zip(values)
                .map(|(&j, y)| y * b[j])
                .sum::<f64>()
        })
        .sum()
}
fn near(a: f64, b: f64, scale: f64) {
    assert!(
        (a - b).abs() <= 3e-11 * scale,
        "{a:e} != {b:e}; scale {scale:e}"
    );
}

// Independently integrate the actual heterogeneous volume: five axial sites
// and two per transverse direction. No production mass coefficients, moments,
// four-point rule, or beam shape functions enter this kinetic-work oracle.
fn volume(
    l: f64,
    patches: &[Patch],
    a: impl Fn(f64, f64, f64) -> [f64; 3],
    b: impl Fn(f64, f64, f64) -> [f64; 3],
) -> f64 {
    let axial = [
        (0.046910077030668, 0.11846344252809454),
        (0.23076534494715845, 0.23931433524968324),
        (0.5, 0.28444444444444444),
        (0.7692346550528415, 0.23931433524968324),
        (0.953089922969332, 0.11846344252809454),
    ];
    let g = 1.0 / 3.0_f64.sqrt();
    let mut out = 0.0;
    for p in patches {
        for (t, wx) in axial {
            let rho = (1.0 - t) * p.rho[0] + t * p.rho[1];
            for eta in [-g, g] {
                for zeta in [-g, g] {
                    let y = 0.5 * ((1.0 - eta) * p.y[0] + (1.0 + eta) * p.y[1]);
                    let z = 0.5 * ((1.0 - zeta) * p.z[0] + (1.0 + zeta) * p.z[1]);
                    out += l
                        * wx
                        * (p.y[1] - p.y[0])
                        * (p.z[1] - p.z[0])
                        * 0.25
                        * rho
                        * dot(a(l * t, y, z), b(l * t, y, z));
                }
            }
        }
    }
    out
}
fn spatial(p0: [f64; 3], axes: [[f64; 3]; 3], x: f64, y: f64, z: f64) -> [f64; 3] {
    core::array::from_fn(|i| p0[i] + x * axes[0][i] + y * axes[1][i] + z * axes[2][i])
}
fn rigid_metric(l: f64, patches: &[Patch], p0: [f64; 3], axes: [[f64; 3]; 3]) -> [[f64; 6]; 6] {
    core::array::from_fn(|a| {
        core::array::from_fn(|b| {
            volume(
                l,
                patches,
                |x, y, z| rigid_velocity(spatial(p0, axes, x, y, z), a),
                |x, y, z| rigid_velocity(spatial(p0, axes, x, y, z), b),
            )
        })
    })
}

#[test]
fn moment_admission_resolves_scaled_and_nearly_singular_physical_inventory() {
    let p = patches();
    let q = moment(&p, 0.0);
    assert!(q.first_y_kg() != 0.0 && q.first_z_kg() != 0.0 && q.second_yz_kg_m() != 0.0);
    // Positive raw moment matrices remain valid even when all floating triple
    // products underflow, overflow, or physical coordinate scales differ.
    for scale in [1e-320, 1e-200, 1.0, 1e200, 1e300] {
        assert!(BeamSectionInertia::new(scale, 0.0, 0.0, scale, 0.0, scale).is_ok());
    }
    let scales = [2.0_f64.powi(300), 2.0_f64.powi(-250), 2.0_f64.powi(150)];
    let base = [2.0, 1.0, -0.5, 3.0, 0.25, 4.0];
    let [a, b, c] = scales;
    assert!(
        BeamSectionInertia::new(
            base[0] * a * a,
            base[1] * a * b,
            base[2] * a * c,
            base[3] * b * b,
            base[4] * b * c,
            base[5] * c * c
        )
        .is_ok()
    );
    let up = f64::from_bits(0.5_f64.to_bits() + 1);
    let down = f64::from_bits(0.5_f64.to_bits() - 1);
    // This rank-two matrix has a strictly positive leading 2x2 minor.
    assert!(BeamSectionInertia::new(1.0, 0.0, 0.5, 1.0, 0.5, 0.5).is_err());
    assert!(BeamSectionInertia::new(1.0, 0.0, 0.5, 1.0, 0.5, up).is_ok());
    assert!(BeamSectionInertia::new(1.0, 0.0, 0.5, 1.0, 0.5, down).is_err());
    assert!(BeamSectionInertia::new(1.0, 2.0, 2.0, 1.0, 2.0, 1.0).is_err());
    for bad in [f64::NAN, f64::INFINITY, -1.0, 0.0] {
        assert!(BeamSectionInertia::new(bad, 0.0, 0.0, 1.0, 0.0, 1.0).is_err());
    }
    assert!(BeamSectionInertia::new(1.0, 0.0, 0.0, 1.0, f64::NAN, 1.0).is_err());
}

#[test]
fn every_rigid_moment_matches_volume_under_orientation_and_endpoint_reversal() {
    let patches = patches();
    let l = 1.3;
    let p0 = [0.4, -0.3, 0.2];
    let moments = [moment(&patches, 0.0), moment(&patches, 1.0)];
    for axes in [
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        [
            [1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0],
            [2.0 / 3.0, -2.0 / 3.0, 1.0 / 3.0],
            [2.0 / 3.0, 1.0 / 3.0, -2.0 / 3.0],
        ],
    ] {
        let p1 = spatial(p0, axes, l, 0.0, 0.0);
        let expected = rigid_metric(l, &patches, p0, axes);
        for kind in kinds() {
            for reverse in [false, true] {
                let (nodes, field) = if reverse {
                    ([p1, p0], [reversed_y(moments[1]), reversed_y(moments[0])])
                } else {
                    ([p0, p1], moments)
                };
                let (_, mass) = segment_matrices_with_reference_inertia(
                    nodes[0],
                    nodes[1],
                    axes[2],
                    section(),
                    kind,
                    field,
                )
                .unwrap();
                for a in 0..6 {
                    for b in 0..6 {
                        near(
                            bilinear(&mass, &rigid(&nodes, a), &rigid(&nodes, b)),
                            expected[a][b],
                            (expected[a][a] * expected[b][b]).sqrt(),
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn independent_volume_keeps_all_nonrigid_translation_rotation_cross_work() {
    let patches = patches();
    let field = [moment(&patches, 0.0), moment(&patches, 1.0)];
    let l = 0.31;
    for kind in kinds() {
        for shear_pa in [8e6, 8e8, 8e10] {
            let section = BeamSection {
                shear_pa,
                ..section()
            };
            let theta_y = [0.7, -0.8, 0.6];
            let theta_z = [-0.4, 0.5, -0.3];
            let gamma = |inertia: f64, area: f64, quadratic: f64| match kind {
                BeamKind::EulerBernoulli => 0.0,
                BeamKind::Timoshenko(_) => {
                    -2.0 * section.young_pa * inertia * quadratic / (section.shear_pa * area)
                }
            };
            let make_v = |v0: f64, t: [f64; 3], g: f64| [v0, t[0] + g, 0.5 * t[1], t[2] / 3.0];
            let v_y = make_v(-0.2, theta_y, gamma(section.iz_m4, 0.0027, theta_y[2]));
            let v_z = make_v(0.1, theta_z, gamma(section.iy_m4, 0.0023, theta_z[2]));
            let evaluate = |c: &[f64], x: f64| c.iter().rev().fold(0.0, |v, c| v * x + c);
            let motion = |x: f64| {
                (
                    [0.4 - 0.7 * x, evaluate(&v_y, x), evaluate(&v_z, x)],
                    [
                        -0.6 + 0.8 * x,
                        -evaluate(&theta_z, x),
                        evaluate(&theta_y, x),
                    ],
                )
            };
            let q: Vec<f64> = [0.0, l]
                .into_iter()
                .flat_map(|x| {
                    let (v, w) = motion(x);
                    [v[0], v[1], v[2], w[0], w[1], w[2]]
                })
                .collect();
            let physical = |x, y, z| {
                let (v, w) = motion(x);
                add(v, cross(w, [0.0, y, z]))
            };
            let exact = volume(l, &patches, physical, physical);
            let (k, m) = segment_matrices_with_reference_inertia(
                [0.0; 3],
                [l, 0.0, 0.0],
                [0.0, 0.0, 1.0],
                section,
                kind,
                field,
            )
            .unwrap();
            near(bilinear(&m, &q, &q), exact, exact);
            assert!(exact > 0.0);
            let centered = field.map(|q| {
                BeamSectionInertia::new(
                    q.line_mass_kg_per_m(),
                    0.0,
                    0.0,
                    q.second_yy_kg_m(),
                    0.0,
                    q.second_zz_kg_m(),
                )
                .unwrap()
            });
            let (k_center, m_center) = segment_matrices_with_reference_inertia(
                [0.0; 3],
                [l, 0.0, 0.0],
                [0.0, 0.0, 1.0],
                section,
                kind,
                centered,
            )
            .unwrap();
            assert_eq!(k, k_center);
            assert!((bilinear(&m, &q, &q) - bilinear(&m_center, &q, &q)).abs() > 1e-4 * exact);
        }
    }
}

#[test]
fn original_scalar_and_p1_density_paths_keep_their_exact_arithmetic() {
    let s = section();
    let z = [0.0, 0.0, 1.0];
    for rho in [[s.density_kg_m3; 2], [700.0; 2], [200.0, 800.0]] {
        let field = rho.map(|r| {
            BeamSectionInertia::from_uniform(&BeamSection {
                density_kg_m3: r,
                ..s
            })
            .unwrap()
        });
        for kind in kinds() {
            let a = segment_matrices_with_reference_density(
                [0.0; 3],
                [0.7, 0.0, 0.0],
                z,
                s,
                kind,
                rho,
                BeamMass::Consistent,
            )
            .unwrap();
            let b = segment_matrices_with_reference_inertia(
                [0.0; 3],
                [0.7, 0.0, 0.0],
                z,
                s,
                kind,
                field,
            )
            .unwrap();
            assert_eq!(a.0.map(f64::to_bits), b.0.map(f64::to_bits));
            assert_eq!(a.1.map(f64::to_bits), b.1.map(f64::to_bits));
        }
    }
}

#[test]
fn native_chain_preserves_section_jumps_and_prestress_retains_inventory() {
    let first = patches();
    let second = first.map(|p| Patch {
        rho: [2.0 * p.rho[1], 0.7 * p.rho[0]],
        ..p
    });
    let field = [
        [moment(&first, 0.0), moment(&first, 1.0)],
        [moment(&second, 0.0), moment(&second, 1.0)],
    ];
    let nodes = [[0.0; 3], [0.4, 0.0, 0.0], [0.4, 0.7, 0.0]];
    let z = [0.0, 0.0, 1.0];
    let (k, m) =
        assemble_beam_chain_with_reference_inertia(&nodes, z, section(), kinds()[1], &field)
            .unwrap();
    let a = rigid_metric(0.4, &first, nodes[0], [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], z]);
    let b = rigid_metric(
        0.7,
        &second,
        nodes[1],
        [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], z],
    );
    for i in 0..6 {
        for j in 0..6 {
            near(
                sparse(&m, &rigid(&nodes, i), &rigid(&nodes, j)),
                a[i][j] + b[i][j],
                ((a[i][i] + b[i][i]) * (a[j][j] + b[j][j])).sqrt(),
            );
        }
    }
    let loaded = assemble_prestressed_beam_chain_with_reference_inertia(
        &nodes,
        z,
        section(),
        kinds()[1],
        &field,
        &[20.0, -5.0],
    )
    .unwrap();
    assert_eq!(m, loaded.1);
    assert_ne!(k, loaded.0);
    assert!(
        assemble_beam_chain_with_reference_inertia(&nodes, z, section(), kinds()[1], &field[..1])
            .is_err()
    );
    assert!(
        assemble_prestressed_beam_chain_with_reference_inertia(
            &nodes,
            z,
            section(),
            kinds()[1],
            &field,
            &[20.0]
        )
        .is_err()
    );
}

#[test]
fn bonded_and_off_node_shell_consumers_retain_full_section_volume_work() {
    let mesh = ShellMesh::new(
        vec![[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0]],
        vec![[0, 1, 3], [0, 3, 2]],
    )
    .unwrap();
    let plate = PlateSection {
        d: [0.5, 0.1, 0.0, 0.1, 0.5, 0.0, 0.0, 0.0, 0.2],
        thickness: 0.02,
        density: 450.0,
    };
    let panels = [[ReferenceDensityMoments::new(500.0, 40.0, 80.0).unwrap(); 3]; 2];
    let patches = patches();
    let field = [moment(&patches, 0.0), moment(&patches, 1.0)];
    let beam = ShellBeam {
        nodes: [0, 1],
        offsets_m: [[0.0, 0.0, 0.05]; 2],
        section_z: [0.0, 0.0, 1.0],
        section: section(),
    };
    let bare = assemble_stiffened_shell_with_reference_inertia(
        &mesh,
        &[plate; 2],
        &panels,
        &[],
        ShellSupport::Free,
        &[],
        &[],
    )
    .unwrap();
    let bonded = assemble_stiffened_shell_with_reference_inertia(
        &mesh,
        &[plate; 2],
        &panels,
        &[],
        ShellSupport::Free,
        &[beam],
        &[field],
    )
    .unwrap();
    let thickness = [0.02; 4];
    let components = [ShellComponent {
        id: 0,
        mesh: &mesh,
        model: &bare,
        nodal_thickness_m: &thickness,
    }];
    let budget = TieBudget {
        max_segments: 100,
        max_slave_dofs: 100,
        gap_tolerance_m: 1e-10,
        angle_tolerance_rad: 1e-8,
    };
    let tied = tie_shells(&components, &[], budget, || false).unwrap();
    let nodes = [[0.0, 0.0, 0.05], [0.4, 0.0, 0.05], [1.0, 0.0, 0.05]];
    let chain = [BeamOnShell {
        beam_nodes_m: &nodes,
        section_z: [0.0, 0.0, 1.0],
        section: section(),
        kind: kinds()[1],
        master: TieSide {
            component: 0,
            nodes: vec![0, 1],
        },
        side: 1,
        centroid_height_m: 0.04,
    }];
    let profiles = [
        [field[0], moment(&patches, 0.4)],
        [moment(&patches, 0.4), field[1]],
    ];
    let coupled = tie_beams_to_shells_with_reference_inertia(
        &tied,
        &components,
        &chain,
        &[&profiles],
        budget,
        || false,
    )
    .unwrap();
    let expected = rigid_metric(
        1.0,
        &patches,
        nodes[0],
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
    );
    for i in 0..6 {
        for j in 0..6 {
            let (a, b) = (rigid(&mesh.nodes, i), rigid(&mesh.nodes, j));
            let shell_work = sparse(&bare.m, &a, &b);
            for model in [&bonded, &coupled.model] {
                near(
                    sparse(&model.m, &a, &b) - shell_work,
                    expected[i][j],
                    (expected[i][i] * expected[j][j]).sqrt(),
                );
            }
        }
    }
    // Supported columns still accelerate the complete inventory. Its full
    // force and torque must retain the transverse density centroid even when
    // the elastic beam endpoint is clamped.
    let (_, body) = assemble_stiffened_shell_with_reference_inertia_and_body_acceleration(
        &mesh,
        &[plate; 2],
        &panels,
        &[0],
        ShellSupport::Clamped,
        &[beam],
        &[field],
    )
    .unwrap();
    let (_, shell_body) = assemble_stiffened_shell_with_reference_inertia_and_body_acceleration(
        &mesh,
        &[plate; 2],
        &panels,
        &[0],
        ShellSupport::Clamped,
        &[],
        &[],
    )
    .unwrap();
    let acceleration = [1.7, -2.4, 8.2];
    let force = body.full_load(acceleration).unwrap();
    let shell_force = shell_body.full_load(acceleration).unwrap();
    assert_eq!(force.len(), 6 * mesh.nodes.len());
    for i in 0..6 {
        let actual: f64 = rigid(&mesh.nodes, i)
            .iter()
            .zip(&force)
            .zip(&shell_force)
            .map(|((v, f), s)| v * (f - s))
            .sum();
        let want: f64 = (0..3).map(|j| expected[i][j] * acceleration[j]).sum();
        near(
            actual,
            want,
            (expected[i][i] * expected[0][0] * dot(acceleration, acceleration)).sqrt(),
        );
    }
    let loaded = tie_prestressed_beams_to_shells_with_reference_inertia(
        &tied,
        &components,
        &chain,
        &[&profiles],
        &[&[20.0, -5.0]],
        budget,
        || false,
    )
    .unwrap();
    assert_eq!(loaded.model.m, coupled.model.m);
    assert_ne!(loaded.model.k, coupled.model.k);
    let lift = coupled.expand_beam(0, &rigid(&mesh.nodes, 4)).unwrap();
    for (a, b) in lift.into_iter().zip(rigid(&nodes, 4)) {
        near(a, b, 1.0);
    }
    assert!(
        tie_beams_to_shells_with_reference_inertia(&tied, &components, &chain, &[], budget, || {
            false
        })
        .is_err()
    );
    assert!(
        tie_beams_to_shells_with_reference_inertia(
            &tied,
            &components,
            &chain,
            &[&profiles[..1]],
            budget,
            || false
        )
        .is_err()
    );
    assert!(
        tie_beams_to_shells_with_reference_inertia(
            &tied,
            &components,
            &chain,
            &[&profiles],
            budget,
            || true
        )
        .is_err()
    );
    assert!(
        assemble_stiffened_shell_with_reference_inertia(
            &mesh,
            &[plate; 2],
            &panels,
            &[],
            ShellSupport::Free,
            &[beam],
            &[]
        )
        .is_err()
    );
}
