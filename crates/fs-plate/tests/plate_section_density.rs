//! Independent through-thickness inventories in the native flat-plate mass.
use fs_plate::{
    AssemblyOptions, EdgeSupport, PlateChart, PlateError, PlateMesh, PlateModel, PlateSection,
    PlateSectionMassMoments, RotaryMass, Stiffener, StiffenerMass, TransverseMass,
};

fn close(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 8e-12 * scale.max(actual.abs()).max(expected.abs()),
        "{actual} != {expected}"
    );
}

fn options() -> AssemblyOptions {
    AssemblyOptions {
        support: EdgeSupport::Clamped,
        pretension: 0.0,
    }
}

fn triangle(thickness: f64) -> PlateChart {
    let mesh =
        PlateMesh::from_unstructured(vec![(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)], vec![[0, 1, 2]])
            .unwrap();
    PlateChart::with_boundary_and_regions(
        mesh,
        PlateSection::isotropic(1000.0, 0.25, thickness, 8.0).unwrap(),
        vec![],
        vec![],
    )
    .unwrap()
}

fn image(chart: &PlateChart, kind: usize) -> PlateModel {
    image_with_beams(chart, kind, &[], StiffenerMass::Lumped)
}

fn image_with_beams(
    chart: &PlateChart,
    kind: usize,
    beams: &[Stiffener],
    beam_mass: StiffenerMass,
) -> PlateModel {
    match kind {
        0 => chart.assemble_with_mass(beams, &options(), TransverseMass::Lumped, beam_mass),
        1 => chart.assemble_with_mass(beams, &options(), TransverseMass::Linear, beam_mass),
        2 => chart.assemble_with_mass(beams, &options(), TransverseMass::EdgeCubic, beam_mass),
        3..=5 => chart.assemble_with_kinetic_mass(
            beams,
            &options(),
            [
                TransverseMass::Lumped,
                TransverseMass::Linear,
                TransverseMass::EdgeCubic,
            ][kind - 3],
            RotaryMass::Dkt,
            beam_mass,
        ),
        6 => chart.assemble_edge_cubic_kinetic_mass(beams, &options(), beam_mass),
        _ => unreachable!(),
    }
    .unwrap()
}

fn action(model: &PlateModel, x: &[f64]) -> Vec<f64> {
    let mut y = vec![0.0; model.free];
    model.m.spmv(x, &mut y);
    y
}

fn product(model: &PlateModel, x: &[f64], y: &[f64]) -> f64 {
    x.iter().zip(action(model, y)).map(|(a, b)| a * b).sum()
}

fn monomial(chart: &PlateChart, powers: [usize; 2]) -> Vec<f64> {
    chart
        .mesh
        .nodes
        .iter()
        .flat_map(|&(x, y)| {
            let at = |a: usize, b: usize| x.powi(a as i32) * y.powi(b as i32);
            [
                at(powers[0], powers[1]),
                if powers[0] > 0 {
                    powers[0] as f64 * at(powers[0] - 1, powers[1])
                } else {
                    0.0
                },
                if powers[1] > 0 {
                    powers[1] as f64 * at(powers[0], powers[1] - 1)
                } else {
                    0.0
                },
            ]
        })
        .collect()
}

fn integral(x: usize, y: usize) -> f64 {
    let factorial = |n: usize| (1..=n).fold(1.0, |v, i| v * i as f64);
    factorial(x) * factorial(y) / factorial(x + y + 2)
}

fn moments(i0: f64, i2: f64) -> PlateSectionMassMoments {
    PlateSectionMassMoments::try_new(i0, i2).unwrap()
}

// Independent Cartesian integration of two DIFFERENT physical fields:
// I0(x,y) = 3 + 6*x + 9*y and I2(x,y) = 0.4 + 0.2*x - 0.1*y.
fn weighted_integral(field: [f64; 3], x: usize, y: usize) -> f64 {
    field[0] * integral(x, y) + field[1] * integral(x + 1, y) + field[2] * integral(x, y + 1)
}

const I0: [f64; 3] = [3.0, 6.0, 9.0];
const I2: [f64; 3] = [0.4, 0.2, -0.1];

fn graded() -> PlateChart {
    triangle(0.4)
        .with_reference_section_mass(vec![[
            moments(3.0, 0.4),
            moments(9.0, 0.6),
            moments(12.0, 0.3),
        ]])
        .unwrap()
}

#[test]
fn independent_moments_reach_every_native_mass_image() {
    let chart = graded();
    let rigid = monomial(&chart, [0, 0]);
    for kind in 0..7 {
        let model = image(&chart, kind);
        assert_eq!(model.k, image(&triangle(0.4), kind).k);
        close(
            product(&model, &rigid, &rigid),
            weighted_integral(I0, 0, 0),
            4.0,
        );
        for p in [[1, 0], [0, 1]] {
            close(
                product(&model, &rigid, &monomial(&chart, p)),
                weighted_integral(I0, p[0], p[1]),
                4.0,
            );
        }
        for i in 0..model.free {
            assert!(model.m.get(i, i) > 0.0);
            for j in 0..model.free {
                assert_eq!(model.m.get(i, j), model.m.get(j, i));
            }
        }
    }

    let powers = [[0, 0], [1, 0], [0, 1], [2, 0], [1, 1], [0, 2]];
    for kind in [1, 4, 5, 6] {
        let model = image(&chart, kind);
        let fields = if kind == 1 || kind == 4 {
            &powers[..3]
        } else {
            &powers[..]
        };
        for a in fields {
            for b in fields {
                let mut expected = weighted_integral(I0, a[0] + b[0], a[1] + b[1]);
                if a[0] * b[0] != 0 {
                    expected +=
                        (a[0] * b[0]) as f64 * weighted_integral(I2, a[0] + b[0] - 2, a[1] + b[1]);
                }
                if a[1] * b[1] != 0 {
                    expected +=
                        (a[1] * b[1]) as f64 * weighted_integral(I2, a[0] + b[0], a[1] + b[1] - 2);
                }
                close(
                    product(&model, &monomial(&chart, *a), &monomial(&chart, *b)),
                    expected,
                    4.0,
                );
            }
        }
    }

    // Full degree-seven transverse integral and degree-five rotary integral.
    // The original edge-cubic reconstruction of nodal x^3 is
    // v = x^3 - x*y + x^2*y + x*y^2. Its gradients are integrated directly.
    let square = |polynomial: &[(f64, usize, usize)], field: [f64; 3]| -> f64 {
        polynomial
            .iter()
            .flat_map(|&(a, x, y)| {
                polynomial
                    .iter()
                    .map(move |&(b, u, v)| a * b * weighted_integral(field, x + u, y + v))
            })
            .sum()
    };
    let transverse = [(1.0, 3, 0), (-1.0, 1, 1), (1.0, 2, 1), (1.0, 1, 2)];
    let cubic_dx = [(3.0, 2, 0), (-1.0, 0, 1), (2.0, 1, 1), (1.0, 0, 2)];
    let cubic_dy = [(-1.0, 1, 0), (1.0, 2, 0), (2.0, 1, 1)];
    let dkt_x = [(3.0, 2, 0), (1.5, 1, 1)];
    let dkt_y = [(1.5, 1, 1)];
    let cubic = [0.0, 0.0, 0.0, 1.0, 3.0, 0.0, 0.0, 0.0, 0.0];
    for (kind, rotation) in [
        (5, square(&dkt_x, I2) + square(&dkt_y, I2)),
        (6, square(&cubic_dx, I2) + square(&cubic_dy, I2)),
    ] {
        close(
            product(&image(&chart, kind), &cubic, &cubic),
            square(&transverse, I0) + rotation,
            4.0,
        );
    }
}

#[test]
fn dense_faces_and_dense_core_have_equal_mass_but_different_native_inertia() {
    // Both symmetric sections occupy -1 <= z <= 1. The core is |z| < 0.5.
    // Exchange densities 2 and 8 between equal-volume core and face layers.
    let layers = |core: f64, faces: f64| {
        let i0 = 2.0 * (core * 0.5 + faces * (1.0 - 0.5));
        let i2 = 2.0 / 3.0 * (core * 0.5_f64.powi(3) + faces * (1.0 - 0.5_f64.powi(3)));
        moments(i0, i2)
    };
    let face = layers(2.0, 8.0);
    let core = layers(8.0, 2.0);
    assert_eq!(face.areal_mass_kg_m2(), core.areal_mass_kg_m2());
    close(face.second_moment_kg() - core.second_moment_kg(), 3.0, 3.0);
    let face_chart = triangle(2.0)
        .with_reference_section_mass(vec![[face; 3]])
        .unwrap();
    let core_chart = triangle(2.0)
        .with_reference_section_mass(vec![[core; 3]])
        .unwrap();
    let rigid = monomial(&face_chart, [0, 0]);
    let tilt = monomial(&face_chart, [1, 0]);
    for kind in 0..7 {
        let a = image(&face_chart, kind);
        let b = image(&core_chart, kind);
        assert_eq!(a.k, b.k);
        close(product(&a, &rigid, &rigid), 5.0, 5.0);
        close(product(&b, &rigid, &rigid), 5.0, 5.0);
        close(
            product(&a, &tilt, &tilt) - product(&b, &tilt, &tilt),
            1.5,
            5.0,
        );
    }
}

#[test]
fn section_inventory_preserves_density_compatibility_and_independent_beams() {
    let rho = [3.0, 9.0, 12.0];
    let h = 0.4;
    let density = triangle(h).with_reference_density(vec![rho]).unwrap();
    let same = triangle(h)
        .with_reference_section_mass(vec![rho.map(|r| moments(r * h, r * h * h * h / 12.0))])
        .unwrap();
    let independent = graded();
    let beam = Stiffener {
        nodes: vec![0, 1],
        e: 1200.0,
        g: 480.0,
        area: 0.04,
        inertia: 0.0002,
        torsion: 0.0001,
        eccentricity: 0.03,
        density: 7.0,
    };
    for kind in 0..7 {
        let original = image(&density, kind);
        let equivalent = image(&same, kind);
        assert_eq!(original.k, equivalent.k);
        for i in 0..original.free {
            for j in 0..original.free {
                close(original.m.get(i, j), equivalent.m.get(i, j), 2.0);
            }
        }
        for beam_mass in [
            StiffenerMass::Lumped,
            StiffenerMass::ConsistentHermite,
            StiffenerMass::ConsistentEccentric,
        ] {
            let old = image_with_beams(&density, kind, std::slice::from_ref(&beam), beam_mass);
            let new = image_with_beams(&independent, kind, std::slice::from_ref(&beam), beam_mass);
            let bare = image(&independent, kind);
            for i in 0..original.free {
                for j in 0..original.free {
                    close(
                        old.m.get(i, j) - original.m.get(i, j),
                        new.m.get(i, j) - bare.m.get(i, j),
                        2.0,
                    );
                }
            }
        }
    }
    // All assembly front doors, including legacy defaults, consume the field.
    assert_eq!(
        independent.assemble(&[], &options()).unwrap().m,
        image(&independent, 0).m
    );
    assert_eq!(
        independent
            .assemble_consistent_transverse_mass(&[], &options())
            .unwrap()
            .m,
        image(&independent, 1).m
    );
    assert_eq!(
        independent
            .assemble_edge_cubic_transverse_mass(&[], &options())
            .unwrap()
            .m,
        image(&independent, 2).m
    );
    let mut changed = independent.clone();
    changed.section = PlateSection::isotropic(2000.0, 0.2, 0.8, 100.0).unwrap();
    for kind in 0..7 {
        assert_eq!(image(&independent, kind).m, image(&changed, kind).m);
    }
    // Reduction uses the same full inventory; no supports are baked into it.
    let full = image(&independent, 6);
    let mut supported = independent.clone();
    supported.boundary_nodes = vec![0, 1];
    let reduced = image(&supported, 6);
    for (i, &r) in reduced.dof_map.iter().enumerate() {
        for (j, &c) in reduced.dof_map.iter().enumerate() {
            if let (Some(r), Some(c)) = (r, c) {
                close(reduced.m.get(r, c), full.m.get(i, j), 4.0);
            }
        }
    }
}

#[test]
fn moment_admission_inventory_selection_and_unrepresentable_mass_refuse() {
    for value in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(PlateSectionMassMoments::try_new(value, 1.0).is_err());
        assert!(PlateSectionMassMoments::try_new(1.0, value).is_err());
    }
    assert!(matches!(
        triangle(0.4).with_reference_section_mass(vec![]),
        Err(PlateError::SectionCount { .. })
    ));
    let selected = graded().with_reference_density(vec![[3.0; 3]]).unwrap();
    assert!(selected.reference_section_mass().is_none());
    assert_eq!(selected.reference_density(), Some([[3.0; 3]].as_slice()));
    let restored = selected
        .with_reference_section_mass(vec![[moments(2.0, 0.1); 3]])
        .unwrap();
    assert!(restored.reference_density().is_none());
    assert!(restored.reference_section_mass().is_some());
    let mut changed = restored.clone();
    changed.mesh.tris.push([0, 1, 2]);
    assert!(matches!(
        changed.assemble(&[], &options()),
        Err(PlateError::SectionCount { .. })
    ));
    let underflow = triangle(0.4)
        .with_reference_section_mass(vec![[moments(1.0, f64::from_bits(1)); 3]])
        .unwrap();
    for kind in 0..7 {
        let result = if kind == 6 {
            underflow.assemble_edge_cubic_kinetic_mass(&[], &options(), StiffenerMass::Lumped)
        } else {
            underflow.assemble_with_kinetic_mass(
                &[],
                &options(),
                [
                    TransverseMass::Lumped,
                    TransverseMass::Linear,
                    TransverseMass::EdgeCubic,
                ][kind % 3],
                if kind < 3 {
                    RotaryMass::Lumped
                } else {
                    RotaryMass::Dkt
                },
                StiffenerMass::Lumped,
            )
        };
        assert!(matches!(result, Err(PlateError::BadSection { .. })));
    }

    // Many finite element inventories meet at one shared node. Their sparse
    // sum must refuse overflow even when each element matrix is representable.
    let count = 64;
    let angle = std::f64::consts::TAU / count as f64;
    let radius = angle.sin().recip().sqrt();
    let mut nodes = vec![(0.0, 0.0)];
    nodes.extend((0..count).map(|i| {
        let theta = i as f64 * angle;
        (radius * theta.cos(), radius * theta.sin())
    }));
    let tris = (0..count)
        .map(|i| [0, i + 1, (i + 1) % count + 1])
        .collect();
    let mesh = PlateMesh::from_unstructured(nodes, tris).unwrap();
    let overflow =
        PlateChart::with_boundary_and_regions(mesh, triangle(0.4).section, vec![], vec![])
            .unwrap()
            .with_reference_section_mass(vec![[moments(f64::MAX, 1.0); 3]; count])
            .unwrap();
    for kind in 0..7 {
        let result = if kind == 6 {
            overflow.assemble_edge_cubic_kinetic_mass(&[], &options(), StiffenerMass::Lumped)
        } else {
            overflow.assemble_with_kinetic_mass(
                &[],
                &options(),
                [
                    TransverseMass::Lumped,
                    TransverseMass::Linear,
                    TransverseMass::EdgeCubic,
                ][kind % 3],
                if kind < 3 {
                    RotaryMass::Lumped
                } else {
                    RotaryMass::Dkt
                },
                StiffenerMass::Lumped,
            )
        };
        assert!(
            matches!(result, Err(PlateError::BadSection { .. })),
            "assembled section-moment overflow must refuse for mass image {kind}"
        );
    }
}
