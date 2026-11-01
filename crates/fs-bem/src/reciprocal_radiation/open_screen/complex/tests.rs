use super::*;

fn surface(offset: Point) -> OpenPanels {
    let mut vertices = Vec::new();
    for j in 0..4 {
        for i in 0..4 {
            let x = (i as f64 - 1.5) * 0.25;
            let y = (j as f64 - 1.5) * 0.25;
            vertices.push([x + offset[0], y + offset[1], 0.125 * x * y + offset[2]]);
        }
    }
    let mut triangles = Vec::new();
    for j in 0..3 {
        for i in 0..3 {
            let a = 4 * j + i;
            triangles.push([a, a + 1, a + 5]);
            triangles.push([a, a + 5, a + 4]);
        }
    }
    OpenPanels::new(vertices, triangles).unwrap()
}
fn options(workers: usize) -> OpenOptions {
    OpenOptions {
        edge_level: 0,
        maximum_edge_level: 0,
        triangle_order: 4,
        angular_refinement: 1,
        maximum_unknowns: 32,
        maximum_cache_bytes: 1_000_000,
        maximum_work: u64::MAX,
        workers,
    }
}
fn limits() -> ScreenComplexLimits {
    ScreenComplexLimits {
        maximum_matrix_entries: 4096,
        maximum_workspace_bytes: 4_000_000,
        maximum_work: u64::MAX,
        maximum_growth_exponent: 8.0,
    }
}
fn near(actual: &[C64], expected: &[C64], tolerance: f64) {
    assert_eq!(actual.len(), expected.len());
    let scale = expected
        .iter()
        .map(|z| z.abs())
        .fold(0.0, f64::max)
        .max(1e-300);
    let error = actual
        .iter()
        .zip(expected)
        .map(|(&a, &b)| (a - b).abs())
        .fold(0.0, f64::max);
    assert!(
        error <= tolerance * scale,
        "relative error {} > {tolerance}",
        error / scale
    );
}

#[test]
fn real_axis_is_the_original_operator_and_uses_the_identical_work_pairing() {
    let surface = surface([0.0; 3]);
    for workers in [1, 3] {
        let owner = PreparedScreen::new(&surface, 1.5, options(workers), || false).unwrap();
        for k in [0.0, 0.2, 1.3] {
            let original = owner.assemble(k, true, &mut || false).unwrap();
            let actual = owner
                .complex_operator(C64::from_re(k), limits(), || false)
                .unwrap();
            assert_eq!(actual.dimension, owner.unknowns());
            for (a, b) in actual.matrix.iter().zip(&original) {
                assert_eq!(
                    (a.re.to_bits(), a.im.to_bits()),
                    (b.re.to_bits(), b.im.to_bits())
                );
            }
            if k == 0.0 {
                assert!(actual.derivative_wavenumber.iter().all(|z| *z == C64::ZERO));
            }
        }
        let rows: Vec<Vec<[f64; 3]>> = vec![
            vec![[1.0; 3]; surface.triangles().len()],
            surface
                .triangles()
                .iter()
                .map(|ids| ids.map(|v| 0.2 + surface.vertices()[v][0]))
                .collect(),
        ];
        let refs: Vec<_> = rows.iter().map(Vec::as_slice).collect();
        let g = owner
            .velocity_port_pairing(&refs, 2 * owner.unknowns(), u64::MAX, || false)
            .unwrap();
        let medium = Medium::air();
        let k = 0.7;
        let source = owner
            .solve_velocity_ports(
                k,
                medium,
                &refs,
                ScreenPortLimits {
                    maximum_ports: 2,
                    maximum_matrix_entries: 4,
                    maximum_jump_entries: 2 * surface.vertices().len(),
                    maximum_work: u64::MAX,
                },
                || false,
            )
            .unwrap();
        let actual = owner
            .complex_operator(C64::from_re(k), limits(), || false)
            .unwrap();
        let factor = fs_la::eigen_complex::lu_complex(&actual.matrix, actual.dimension).unwrap();
        for (j, row) in g.iter().enumerate() {
            let scale = C64::new(0.0, k).scale(medium.density * medium.sound_speed);
            let mut pressure: Vec<_> = row.iter().map(|&x| scale.scale(x)).collect();
            factor.solve(&mut pressure);
            let expected: Vec<_> = owner
                .pressure_dof_map()
                .iter()
                .enumerate()
                .filter_map(|(v, d)| d.map(|_| source.solutions()[j].pressure_jump_pa()[v]))
                .collect();
            near(&pressure, &expected, 5e-13);
            for (i, row) in g.iter().enumerate() {
                let force = row
                    .iter()
                    .zip(&pressure)
                    .fold(C64::ZERO, |sum, (&g, &p)| sum + p.scale(g));
                assert!(
                    (force - source.impedance()[2 * i + j]).abs()
                        <= 5e-13 * source.impedance()[0].abs()
                );
            }
        }
    }
}

#[test]
fn matrix_derivative_is_holomorphic_and_retains_acoustic_reality_symmetry() {
    let surface = surface([0.0; 3]);
    let owner = PreparedScreen::new(&surface, 1.5, options(1), || false).unwrap();
    for k in [C64::from_re(0.6), C64::new(0.7, -0.3), C64::new(-0.8, 0.2)] {
        let operator = owner.complex_operator(k, limits(), || false).unwrap();
        for step in [C64::from_re(1e-5), C64::new(0.0, 1e-5)] {
            let left = owner
                .complex_operator(k - step, limits(), || false)
                .unwrap();
            let right = owner
                .complex_operator(k + step, limits(), || false)
                .unwrap();
            let derivative: Vec<_> = left
                .matrix
                .iter()
                .zip(&right.matrix)
                .map(|(&a, &b)| (b - a) / step.scale(2.0))
                .collect();
            near(&derivative, &operator.derivative_wavenumber, 2e-7);
        }
        let mirror = owner
            .complex_operator(-k.conj(), limits(), || false)
            .unwrap();
        near(
            &mirror.matrix,
            &operator.matrix.iter().map(|z| z.conj()).collect::<Vec<_>>(),
            2e-13,
        );
        near(
            &mirror.derivative_wavenumber,
            &operator
                .derivative_wavenumber
                .iter()
                .map(|z| -z.conj())
                .collect::<Vec<_>>(),
            2e-13,
        );
        let scale = operator.matrix.iter().map(|z| z.abs()).fold(0.0, f64::max);
        for i in 0..operator.dimension {
            for j in 0..operator.dimension {
                assert!(
                    (operator.matrix[i * operator.dimension + j]
                        - operator.matrix[j * operator.dimension + i])
                        .abs()
                        <= 2e-13 * scale
                );
            }
        }
    }
}

#[test]
fn complete_domain_budgets_centered_growth_and_cancellation_are_explicit() {
    let surface = surface([1e8, -1e8, 1e7]);
    let owner = PreparedScreen::new(&surface, 1.5, options(3), || false).unwrap();
    owner
        .admit_complex_domain(0.0, 0.0, limits(), || false)
        .unwrap();
    owner
        .admit_complex_domain(1.4, 0.5, limits(), || false)
        .unwrap();
    let k = C64::new(0.7, -0.3);
    let reference = owner.complex_operator(k, limits(), || false).unwrap();
    assert!(
        reference
            .matrix
            .iter()
            .chain(&reference.derivative_wavenumber)
            .all(|&z| finite(z))
    );
    assert!(
        owner
            .admit_complex_domain(1.6, 0.0, limits(), || false)
            .is_err()
    );
    let mut small = limits();
    small.maximum_growth_exponent = 0.1;
    assert!(
        owner
            .admit_complex_domain(1.0, 0.3, small, || false)
            .is_err()
    );
    for kind in 0..3 {
        let mut short = limits();
        match kind {
            0 => short.maximum_matrix_entries = 2 * owner.unknowns() * owner.unknowns() - 1,
            1 => short.maximum_workspace_bytes = 0,
            _ => short.maximum_work = 0,
        }
        assert!(owner.complex_operator(k, short, || false).is_err());
    }
    assert!(matches!(
        owner.complex_operator(k, limits(), || true),
        Err(HelmholtzError::Cancelled)
    ));
    let mut total = 0;
    owner
        .complex_operator(k, limits(), || {
            total += 1;
            false
        })
        .unwrap();
    let mut count = 0;
    assert!(matches!(
        owner.complex_operator(k, limits(), || {
            count += 1;
            count >= total - 1
        }),
        Err(HelmholtzError::Cancelled)
    ));
    let replay = owner.complex_operator(k, limits(), || false).unwrap();
    assert_eq!(reference.matrix, replay.matrix);
    assert_eq!(
        reference.derivative_wavenumber,
        replay.derivative_wavenumber
    );
}
