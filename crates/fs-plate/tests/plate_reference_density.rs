//! G0/G1/G3: reference density placement in the native DKT kinetic images.
use fs_plate::{
    AssemblyOptions, EdgeSupport, PlateChart, PlateError, PlateMesh, PlateModel, PlateSection,
    RotaryMass, SliceOptions, Stiffener, StiffenerMass, TransverseMass, modes,
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
    match kind {
        0 => chart.assemble(&[], &options()),
        1 => chart.assemble_consistent_transverse_mass(&[], &options()),
        2 => chart.assemble_edge_cubic_transverse_mass(&[], &options()),
        3..=5 => chart.assemble_with_kinetic_mass(
            &[],
            &options(),
            [
                TransverseMass::Lumped,
                TransverseMass::Linear,
                TransverseMass::EdgeCubic,
            ][kind - 3],
            RotaryMass::Dkt,
            StiffenerMass::Lumped,
        ),
        6 => chart.assemble_edge_cubic_kinetic_mass(&[], &options(), StiffenerMass::Lumped),
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

// Independent Cartesian integral of rho(x,y) = 3 + 6*x + 9*y.
fn weighted_integral(x: usize, y: usize) -> f64 {
    3.0 * integral(x, y) + 6.0 * integral(x + 1, y) + 9.0 * integral(x, y + 1)
}

#[test]
fn all_native_kinetic_images_retain_mass_and_centroid() {
    let base = triangle(0.4);
    let graded = base
        .clone()
        .with_reference_density(vec![[3.0, 9.0, 12.0]])
        .unwrap();
    assert_eq!(
        graded.reference_density(),
        Some([[3.0, 9.0, 12.0]].as_slice())
    );
    let rigid = monomial(&graded, [0, 0]);
    let tilt_x = monomial(&graded, [1, 0]);
    let tilt_y = monomial(&graded, [0, 1]);
    for kind in 0..7 {
        let model = image(&graded, kind);
        assert_eq!(model.k, image(&base, kind).k);
        close(product(&model, &rigid, &rigid), 1.6, 1.6);
        close(product(&model, &rigid, &tilt_x), 0.55, 1.6);
        close(product(&model, &rigid, &tilt_y), 0.60, 1.6);
        for i in 0..model.free {
            assert!(model.m.get(i, i) > 0.0);
            for j in 0..model.free {
                assert_eq!(model.m.get(i, j), model.m.get(j, i));
            }
        }
        // A constant transverse acceleration uses the same mass action.
        // Its total force and origin moment include the non-diagonal rows.
        let loads = action(&model, &rigid);
        close(loads.iter().step_by(3).sum(), 1.6, 1.6);
        let moment_x: f64 = graded
            .mesh
            .nodes
            .iter()
            .enumerate()
            .map(|(i, &(_, y))| y * loads[3 * i] + loads[3 * i + 2])
            .sum();
        let moment_y: f64 = graded
            .mesh
            .nodes
            .iter()
            .enumerate()
            .map(|(i, &(x, _))| -x * loads[3 * i] - loads[3 * i + 1])
            .sum();
        close(moment_x, 0.60, 1.6);
        close(moment_y, -0.55, 1.6);
    }
}

#[test]
fn consistent_images_reproduce_all_density_weighted_quadratic_kinetic_products() {
    let h = 0.4;
    let chart = triangle(h)
        .with_reference_density(vec![[3.0, 9.0, 12.0]])
        .unwrap();
    let powers = [[0, 0], [1, 0], [0, 1], [2, 0], [1, 1], [0, 2]];
    for kind in [5, 6] {
        let model = image(&chart, kind);
        for a in powers {
            for b in powers {
                let mut rotation = 0.0;
                if a[0] * b[0] != 0 {
                    rotation +=
                        (a[0] * b[0]) as f64 * weighted_integral(a[0] + b[0] - 2, a[1] + b[1]);
                }
                if a[1] * b[1] != 0 {
                    rotation +=
                        (a[1] * b[1]) as f64 * weighted_integral(a[0] + b[0], a[1] + b[1] - 2);
                }
                let expected =
                    h * weighted_integral(a[0] + b[0], a[1] + b[1]) + h * h * h / 12.0 * rotation;
                close(
                    product(&model, &monomial(&chart, a), &monomial(&chart, b)),
                    expected,
                    h,
                );
            }
        }
        // An equal-mean-density plate has the same mass, but different tilt
        // inertia. Averaging rho before assembly would incorrectly pass a
        // total-mass-only test and fail this physical placement check.
        let mean = image(&triangle(h), kind);
        let tilt = monomial(&chart, [1, 0]);
        assert!((product(&model, &tilt, &tilt) - product(&mean, &tilt, &tilt)).abs() > 0.01);
    }
    for kind in [1, 4] {
        let model = image(&chart, kind);
        for a in &powers[..3] {
            for b in &powers[..3] {
                let rotation = (a[0] * b[0] + a[1] * b[1]) as f64 * weighted_integral(0, 0);
                close(
                    product(&model, &monomial(&chart, *a), &monomial(&chart, *b)),
                    h * weighted_integral(a[0] + b[0], a[1] + b[1]) + h * h * h / 12.0 * rotation,
                    h,
                );
            }
        }
    }
    // Nodal data of x^3 produce the existing edge-cubic reconstruction
    // v = x^3 - (1-x-y)*x*y. Its rho*v^2 is genuinely degree seven,
    // beyond the quadratic completeness exercised above. Independent exact
    // Cartesian integration gives int(rho*v^2)=41/280 and
    // int(rho*|grad(v)|^2)=417/140. The original DKT rotation is different:
    // beta = -(3*x^2 + 1.5*x*y, 1.5*x*y), with integral 195/56.
    let cubic = [0.0, 0.0, 0.0, 1.0, 3.0, 0.0, 0.0, 0.0, 0.0];
    for (kind, rotation) in [(2, 99.0 / 8.0), (5, 195.0 / 56.0), (6, 417.0 / 140.0)] {
        close(
            product(&image(&chart, kind), &cubic, &cubic),
            h * (41.0 / 280.0) + h * h * h / 12.0 * rotation,
            h,
        );
    }
}

#[test]
fn unequal_placement_at_equal_mass_changes_the_actual_supported_modes() {
    let mut left = triangle(0.4);
    left.boundary_nodes = vec![0, 1];
    let first = left
        .clone()
        .with_reference_density(vec![[3.0, 9.0, 12.0]])
        .unwrap();
    let second = left.with_reference_density(vec![[12.0, 9.0, 3.0]]).unwrap();
    let a = image(&first, 1);
    let b = image(&second, 1);
    assert_eq!(a.k, b.k);
    assert_eq!(a.dof_map, b.dof_map);
    assert_eq!(a.free, 3);
    // Analytic consistent P1 transverse and row-summed rotary masses at the
    // free vertex. All supplied profiles have total mass 1.6 kg.
    for (model, diagonal) in [(&a, [0.32, 0.008, 0.008]), (&b, [0.20, 0.006, 0.006])] {
        for (i, expected) in diagonal.into_iter().enumerate() {
            close(model.m.get(i, i), expected, expected);
        }
    }
    let opts = SliceOptions {
        mass_diagonal_equilibration: true,
        ..SliceOptions::default()
    };
    let low = modes(&a, (0.0, 1e6), &opts).unwrap();
    let high = modes(&b, (0.0, 1e6), &opts).unwrap();
    assert_eq!(low.expected, 3);
    assert_eq!(high.expected, 3);
    assert_eq!(low.modes.len(), 3);
    assert_eq!(high.modes.len(), 3);
    for (l, h) in low.modes.iter().zip(&high.modes) {
        // Positive difference of the complete reduced mass matrices gives
        // this ordering by the generalized Rayleigh quotient/min-max law.
        assert!(
            l.lambda < 0.8 * h.lambda,
            "{} versus {}",
            l.lambda,
            h.lambda
        );
        close(product(&a, &l.phi, &l.phi), 1.0, 1.0);
        close(product(&b, &h.phi, &h.phi), 1.0, 1.0);
    }
}

#[test]
fn uniform_profiles_preserve_native_arithmetic_and_independent_beam_inertia() {
    let mesh = PlateMesh::rectangle(1.0, 0.7, 1, 1);
    let a = PlateSection::isotropic(1200.0, 0.25, 0.08, 3.0).unwrap();
    let b = PlateSection::isotropic(900.0, 0.2, 0.06, 5.0).unwrap();
    let chart = PlateChart::with_boundary_and_regions(mesh, a, vec![], vec![])
        .unwrap()
        .with_element_sections(vec![a, b])
        .unwrap();
    let field = chart
        .clone()
        .with_reference_density(vec![[3.0; 3], [5.0; 3]])
        .unwrap();
    for kind in 0..7 {
        assert_eq!(image(&chart, kind).m, image(&field, kind).m);
        assert_eq!(image(&chart, kind).k, image(&field, kind).k);
    }
    let beam = Stiffener {
        nodes: vec![0, 3],
        e: 1200.0,
        g: 480.0,
        area: 0.04,
        inertia: 0.0002,
        torsion: 0.0001,
        eccentricity: 0.03,
        density: 7.0,
    };
    let graded = chart
        .clone()
        .with_reference_density(vec![[3.0, 9.0, 12.0], [12.0, 5.0, 3.0]])
        .unwrap();
    for mass in [
        StiffenerMass::Lumped,
        StiffenerMass::ConsistentHermite,
        StiffenerMass::ConsistentEccentric,
    ] {
        let plain = chart
            .assemble_with_mass(
                std::slice::from_ref(&beam),
                &options(),
                TransverseMass::EdgeCubic,
                mass,
            )
            .unwrap();
        let uniform = field
            .assemble_with_mass(
                std::slice::from_ref(&beam),
                &options(),
                TransverseMass::EdgeCubic,
                mass,
            )
            .unwrap();
        assert_eq!(plain.m, uniform.m);
        let old_bare = chart
            .assemble_with_mass(&[], &options(), TransverseMass::EdgeCubic, mass)
            .unwrap();
        let with_beam = graded
            .assemble_with_mass(
                std::slice::from_ref(&beam),
                &options(),
                TransverseMass::EdgeCubic,
                mass,
            )
            .unwrap();
        let bare = graded
            .assemble_with_mass(&[], &options(), TransverseMass::EdgeCubic, mass)
            .unwrap();
        for i in 0..bare.free {
            for j in 0..bare.free {
                close(
                    with_beam.m.get(i, j) - bare.m.get(i, j),
                    plain.m.get(i, j) - old_bare.m.get(i, j),
                    1.0,
                );
            }
        }
    }
    // Updating constitutive scalar density does not silently resample the
    // explicit reference inventory; the caller chose its physical owner.
    let mut changed = graded.clone();
    changed.section.density = 100.0;
    let mut aa = a;
    aa.density = 100.0;
    let mut bb = b;
    bb.density = 200.0;
    changed = changed.with_element_sections(vec![aa, bb]).unwrap();
    for kind in 0..7 {
        assert_eq!(image(&graded, kind).m, image(&changed, kind).m);
    }
}

#[test]
fn reparameterization_and_partition_preserve_physical_quadratic_energy() {
    let mut reference = triangle(0.4)
        .with_reference_density(vec![[3.0, 9.0, 12.0]])
        .unwrap();
    let mut permuted = triangle(0.4);
    permuted.mesh.tris[0] = [1, 2, 0];
    permuted = permuted
        .with_reference_density(vec![[9.0, 12.0, 3.0]])
        .unwrap();
    for kind in 0..7 {
        let original = image(&reference, kind);
        let other = image(&permuted, kind);
        for i in 0..9 {
            for j in 0..9 {
                close(other.m.get(i, j), original.m.get(i, j), 1.6);
            }
        }
    }
    reference.mesh.nodes = vec![(0.2, -0.3), (0.2, 0.7), (-0.8, -0.3)];
    let original = triangle(0.4)
        .with_reference_density(vec![[3.0, 9.0, 12.0]])
        .unwrap();
    let q = [0.2, 0.7, -0.4, -0.8, 0.1, 0.9, 0.5, -0.3, 0.6];
    let rotated: Vec<_> = q
        .chunks_exact(3)
        .flat_map(|v| [v[0], -v[2], v[1]])
        .collect();
    for kind in 0..7 {
        close(
            product(&image(&reference, kind), &rotated, &rotated),
            product(&image(&original, kind), &q, &q),
            1.6,
        );
    }
    // Same Cartesian rho and quadratic velocity over successively partitioned
    // unit squares. The exact integral is independent of the mesh partition.
    let h = 0.2;
    let rectangle_integral = |a: usize, b: usize| {
        let i = |x: usize, y: usize| 1.0 / ((x + 1) * (y + 1)) as f64;
        3.0 * i(a, b) + 6.0 * i(a + 1, b) + 9.0 * i(a, b + 1)
    };
    let expected = h * rectangle_integral(2, 2)
        + h * h * h / 12.0 * (rectangle_integral(2, 0) + rectangle_integral(0, 2));
    for n in [1, 2, 4] {
        let mesh = PlateMesh::rectangle(1.0, 1.0, n, n);
        let density = mesh
            .tris
            .iter()
            .map(|tri| {
                tri.map(|i| {
                    let (x, y) = mesh.nodes[i];
                    3.0 + 6.0 * x + 9.0 * y
                })
            })
            .collect();
        let chart = PlateChart::with_boundary_and_regions(
            mesh,
            PlateSection::isotropic(1000.0, 0.25, h, 1.0).unwrap(),
            vec![],
            vec![],
        )
        .unwrap()
        .with_reference_density(density)
        .unwrap();
        let velocity = monomial(&chart, [1, 1]);
        for kind in [5, 6] {
            close(
                product(&image(&chart, kind), &velocity, &velocity),
                expected,
                expected,
            );
        }
    }
}

#[test]
fn malformed_density_or_stale_coverage_refuses_before_a_model_is_returned() {
    for density in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            triangle(0.4)
                .with_reference_density(vec![[3.0, density, 12.0]])
                .is_err()
        );
    }
    assert!(matches!(
        triangle(0.4).with_reference_density(vec![]),
        Err(PlateError::SectionCount {
            expected: 1,
            actual: 0
        })
    ));
    let original = triangle(0.4)
        .with_reference_density(vec![[3.0, 9.0, 12.0]])
        .unwrap();
    let mut stale = original.clone();
    stale.mesh.tris.push([0, 1, 2]);
    assert!(matches!(
        stale.assemble(&[], &options()),
        Err(PlateError::SectionCount {
            expected: 2,
            actual: 1
        })
    ));
    let field = original.reference_density().unwrap().to_vec();
    assert_eq!(image(&original, 6).m, image(&original, 6).m);
    assert_eq!(original.reference_density().unwrap(), &field);
}
