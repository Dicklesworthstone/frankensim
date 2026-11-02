//! G0/G1/G3 numerical regressions on original source-owned screen operators.
//! Synthetic caloric data isolate molecular compression; they do not constitute
//! measured atmospheric or open-shell validation. Native execution is separate.
use super::*;
use crate::reciprocal_radiation::open_screen::{ScreenComplexLimits, ScreenComplexValueLimits};
use fs_material::{
    gas::{GasSpec, GasState},
    relaxation::{RelaxingMedium, input::read_relaxation},
};

pub(super) fn exterior(relaxing: bool) -> MolecularExterior {
    let input = if relaxing {
        "declared 1 synthetic 1 vibration 1 900 1 power 101.325 293.15 0 200 500"
    } else {
        "declared 1 synthetic 0"
    };
    let law = read_relaxation(&mut input.split_ascii_whitespace()).unwrap();
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
    let medium = RelaxingMedium::new(gas, vec![(law.species()[0].species.clone(), 1.0)]).unwrap();
    MolecularExterior::new(law.evaluate(&medium).unwrap(), [100.0, 2000.0]).unwrap()
}

pub(super) fn surface() -> OpenPanels {
    let mut vertices = Vec::new();
    for j in 0..4 {
        for i in 0..4 {
            let x = (i as f64 - 1.5) * 0.25;
            let y = (j as f64 - 1.5) * 0.25;
            vertices.push([x, y, 0.125 * x * y]);
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

pub(super) fn options(workers: usize) -> OpenOptions {
    OpenOptions {
        edge_level: 0,
        maximum_edge_level: 0,
        triangle_order: 4,
        angular_refinement: 1,
        maximum_unknowns: 16,
        maximum_cache_bytes: 1_000_000,
        maximum_work: u64::MAX,
        workers,
    }
}

pub(super) fn operator_limits() -> ScreenComplexValueLimits {
    ScreenComplexValueLimits {
        maximum_matrix_entries: 1024,
        maximum_workspace_bytes: 4_000_000,
        maximum_work: u64::MAX,
        maximum_growth_exponent: 8.0,
    }
}

fn derivative_limits() -> ScreenComplexLimits {
    ScreenComplexLimits {
        maximum_matrix_entries: 2048,
        maximum_workspace_bytes: 4_000_000,
        maximum_work: u64::MAX,
        maximum_growth_exponent: 8.0,
    }
}

pub(super) fn limits(surface: &OpenPanels) -> MolecularScreenPortLimits {
    MolecularScreenPortLimits {
        ports: ScreenPortLimits {
            maximum_ports: 2,
            maximum_matrix_entries: 4,
            maximum_jump_entries: 2 * surface.vertices().len(),
            maximum_work: u64::MAX,
        },
        operator: operator_limits(),
        maximum_workspace_bytes: 4_000_000,
    }
}

pub(super) fn area_normal(surface: &OpenPanels, facet: usize) -> (f64, Point) {
    let t = surface.triangles()[facet].map(|v| surface.vertices()[v]);
    let a: Point = core::array::from_fn(|i| t[1][i] - t[0][i]);
    let b: Point = core::array::from_fn(|i| t[2][i] - t[0][i]);
    let cross = [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ];
    let length = cross.iter().map(|x| x * x).sum::<f64>().sqrt();
    (0.5 * length, cross.map(|v| v / length))
}

pub(super) fn velocity_rows(surface: &OpenPanels) -> Vec<Vec<[f64; 3]>> {
    [[0.4, -0.2, 1.0], [0.2, 0.7, 0.4]]
        .iter()
        .enumerate()
        .map(|(mode, direction)| {
            surface
                .triangles()
                .iter()
                .enumerate()
                .map(|(f, ids)| {
                    let (_, normal) = area_normal(surface, f);
                    let projection = (0..3).map(|i| direction[i] * normal[i]).sum::<f64>();
                    ids.map(|v| {
                        let [x, y, _] = surface.vertices()[v];
                        projection
                            * if mode == 0 {
                                1.0 + 0.15 * x
                            } else {
                                0.4 + x - 0.3 * y
                            }
                    })
                })
                .collect()
        })
        .collect()
}

// Independent degree-two quadrature reconstructs BOTH the pressure-test drive
// and its mechanical transpose, without the production P1 mass coefficients.
const WORK_POINTS: [[f64; 3]; 3] = [
    [2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0],
    [1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0],
    [1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0],
];

fn pairing(surface: &OpenPanels, dofs: &[Option<usize>], trace: &[[f64; 3]]) -> Vec<f64> {
    let mut row = vec![0.0; surface.unknowns()];
    for (f, ids) in surface.triangles().iter().enumerate() {
        let (area, _) = area_normal(surface, f);
        for q in WORK_POINTS {
            let v = (0..3).map(|a| q[a] * trace[f][a]).sum::<f64>();
            for a in 0..3 {
                if let Some(i) = dofs[ids[a]] {
                    row[i] += area / 3.0 * v * q[a];
                }
            }
        }
    }
    row
}

pub(super) fn effort(surface: &OpenPanels, trace: &[[f64; 3]], jump: &[C64]) -> C64 {
    let mut value = C64::ZERO;
    for (f, ids) in surface.triangles().iter().enumerate() {
        let (area, _) = area_normal(surface, f);
        for q in WORK_POINTS {
            let v = (0..3).map(|a| q[a] * trace[f][a]).sum::<f64>();
            let p = (0..3).fold(C64::ZERO, |p, a| p + jump[ids[a]].scale(q[a]));
            value = value + p.scale(area / 3.0 * v);
        }
    }
    value
}

pub(super) fn close(actual: &[C64], expected: &[C64], tolerance: f64) {
    assert_eq!(actual.len(), expected.len());
    let scale = expected.iter().map(|z| z.abs()).fold(0.0, f64::max);
    let error = actual
        .iter()
        .zip(expected)
        .map(|(&a, &b)| (a - b).abs())
        .fold(0.0, f64::max);
    assert!(
        error <= tolerance * scale.max(1e-300),
        "error {error}, scale {scale}, tolerance {tolerance}"
    );
}

#[test]
fn empty_reservoir_recovers_original_real_screen_reaction_and_jumps() {
    let surface = surface();
    let material = exterior(false);
    let omega = 650.0;
    let k = material.source().relaxing_wavenumber(omega).unwrap();
    assert_eq!(k.im, 0.0);
    let owner = PreparedScreen::new(&surface, 6.0, options(1), || false).unwrap();
    let rows = velocity_rows(&surface);
    let refs: Vec<_> = rows.iter().map(Vec::as_slice).collect();
    let controls = limits(&surface);
    let original = owner
        .solve_velocity_ports(
            k.re,
            Medium {
                density: material.source().gas().density,
                sound_speed: omega / k.re,
            },
            &refs,
            controls.ports,
            || false,
        )
        .unwrap();
    let actual = owner
        .solve_molecular_velocity_ports(&material, omega, &refs, controls, || false)
        .unwrap();
    assert_eq!(actual.ports(), 2);
    assert_eq!(actual.angular_frequency_rad_s(), omega);
    assert_eq!(actual.wavenumber_per_m(), k);
    assert!(std::ptr::eq(actual.source().source(), material.source()));
    close(actual.impedance(), original.impedance(), 3e-13);
    for (actual, original) in actual.pressure_jumps_pa().iter().zip(original.solutions()) {
        close(actual, original.pressure_jump_pa(), 3e-13);
        for (v, &p) in actual.iter().enumerate() {
            if surface.is_free_edge_vertex(v) {
                assert_eq!(p, C64::ZERO);
            }
        }
    }
    let rates = [C64::new(0.7, -0.2), C64::new(-0.3, 0.6)];
    close(
        &actual.pressure_jump(&rates, || false).unwrap(),
        &original.pressure_jump(&rates, || false).unwrap(),
        5e-13,
    );
    close(
        &actual.resisting_efforts(&rates).unwrap(),
        &original.resisting_efforts(&rates).unwrap(),
        5e-13,
    );
    let work = actual.boundary_work_w(&rates, || false).unwrap();
    let expected = original.input_power_w(&rates).unwrap();
    assert!((work - expected).abs() <= 1e-11 * expected.abs());
    let values = owner
        .complex_matrix(k, operator_limits(), || false)
        .unwrap();
    let native = owner.assemble(k.re, true, &mut || false).unwrap();
    assert_eq!(values, native);
}

#[test]
fn molecular_ports_use_real_momentum_with_full_coherent_work_and_mutual_loading() {
    let surface = surface();
    let material = exterior(true);
    let omega = 650.0;
    let owner = PreparedScreen::new(&surface, 6.0, options(1), || false).unwrap();
    let rows = velocity_rows(&surface);
    let refs: Vec<_> = rows.iter().map(Vec::as_slice).collect();
    let solved = owner
        .solve_molecular_velocity_ports(&material, omega, &refs, limits(&surface), || false)
        .unwrap();
    let k = solved.wavenumber_per_m();
    assert!(k.im > 0.0);
    assert_eq!(k, material.source().relaxing_wavenumber(omega).unwrap());
    assert!(solved.maximum_backward_residual() < 1e-10);
    // The independent derivative-producing owner remains a reference for the
    // new values-only assembly; it also exposes the original pressure rows.
    let h = owner
        .complex_operator(k, derivative_limits(), || false)
        .unwrap();
    let values = owner
        .complex_matrix(k, operator_limits(), || false)
        .unwrap();
    close(&values, &h.matrix, 2e-13);
    let dofs = owner.pressure_dof_map();
    let weak: Vec<_> = rows
        .iter()
        .map(|row| pairing(&surface, dofs, row))
        .collect();
    let rho_omega = material.source().gas().density * omega;
    let physical = C64::new(0.0, rho_omega);
    let false_momentum = C64::new(0.0, 1.0)
        * k.scale(material.source().gas().density * material.source().gas().sound_speed);
    assert!((physical - false_momentum).abs() > 1e-4 * physical.abs());
    let scale = solved
        .impedance()
        .iter()
        .map(|z| z.abs())
        .fold(0.0, f64::max);
    for (j, jump) in solved.pressure_jumps_pa().iter().enumerate() {
        let mut active = vec![C64::ZERO; h.dimension];
        for (v, &dof) in dofs.iter().enumerate() {
            if let Some(i) = dof {
                active[i] = jump[v];
            }
        }
        for i in 0..h.dimension {
            let lhs = (0..h.dimension).fold(C64::ZERO, |sum, a| {
                sum + h.matrix[i * h.dimension + a] * active[a]
            });
            let rhs = physical.scale(weak[j][i]);
            assert!((lhs - rhs).abs() <= 2e-12 * rhs.abs());
        }
        for (i, row) in rows.iter().enumerate() {
            let expected = effort(&surface, row, jump);
            assert!((solved.impedance()[2 * i + j] - expected).abs() <= 3e-13 * scale);
            assert!(
                (solved.impedance()[2 * i + j] - solved.impedance()[2 * j + i]).abs()
                    <= 2e-12 * scale
            );
        }
    }
    assert!(solved.impedance()[1].abs() > 1e-3 * scale);
    let rates = [C64::new(0.7, -0.2), C64::new(-0.3, 0.6)];
    let mut direct: Vec<_> = (0..h.dimension)
        .map(|i| physical * (rates[0].scale(weak[0][i]) + rates[1].scale(weak[1][i])))
        .collect();
    fs_la::eigen_complex::lu_complex(&h.matrix, h.dimension)
        .unwrap()
        .solve(&mut direct);
    let direct: Vec<_> = dofs
        .iter()
        .map(|d| d.map_or(C64::ZERO, |i| direct[i]))
        .collect();
    let combined = solved.pressure_jump(&rates, || false).unwrap();
    close(&combined, &direct, 3e-12);
    let efforts: Vec<_> = rows
        .iter()
        .map(|row| effort(&surface, row, &direct))
        .collect();
    close(&solved.resisting_efforts(&rates).unwrap(), &efforts, 3e-12);
    let work = solved.boundary_work_w(&rates, || false).unwrap();
    let direct_work = 0.5
        * (0..2)
            .map(|i| (rates[i].conj() * efforts[i]).re)
            .sum::<f64>();
    let work_scale = 0.5
        * (0..2)
            .map(|i| rates[i].abs() * efforts[i].abs())
            .sum::<f64>();
    assert!(work > 0.0);
    assert!((work - direct_work).abs() <= 3e-12 * work_scale);
    let diagonal_work = 0.5
        * (0..2)
            .map(|i| rates[i].norm_sq() * solved.impedance()[2 * i + i].re)
            .sum::<f64>();
    assert!((work - diagonal_work).abs() > 1e-4 * work);
    let lossless = exterior(false);
    let lossless = owner
        .solve_molecular_velocity_ports(&lossless, omega, &refs, limits(&surface), || false)
        .unwrap();
    assert!((solved.impedance()[0] - lossless.impedance()[0]).abs() > 1e-5 * scale);
}

#[test]
fn molecular_screen_refuses_unadmitted_sources_resources_and_cancelled_batches() {
    let surface = surface();
    let material = exterior(true);
    let owner = PreparedScreen::new(&surface, 6.0, options(1), || false).unwrap();
    let rows = velocity_rows(&surface);
    let refs: Vec<_> = rows.iter().map(Vec::as_slice).collect();
    let omega = 650.0;
    let controls = limits(&surface);
    for omega in [0.0, f64::NAN, 99.0, 2001.0] {
        assert!(
            owner
                .solve_molecular_velocity_ports(&material, omega, &refs, controls, || false)
                .is_err()
        );
    }
    assert!(
        owner
            .solve_molecular_velocity_ports(&material, omega, &[], controls, || false)
            .is_err()
    );
    let short_row = &rows[0][..rows[0].len() - 1];
    assert!(
        owner
            .solve_molecular_velocity_ports(
                &material,
                omega,
                &[short_row, &rows[1]],
                controls,
                || false
            )
            .is_err()
    );
    for kind in 0..7 {
        let mut short = controls;
        match kind {
            0 => short.ports.maximum_matrix_entries = 3,
            1 => short.ports.maximum_jump_entries -= 1,
            2 => short.ports.maximum_work = 0,
            3 => short.maximum_workspace_bytes = 0,
            4 => short.operator.maximum_matrix_entries = owner.unknowns().pow(2) - 1,
            5 => short.operator.maximum_workspace_bytes = 0,
            _ => short.operator.maximum_growth_exponent = 1e-10,
        }
        assert!(
            owner
                .solve_molecular_velocity_ports(&material, omega, &refs, short, || false)
                .is_err()
        );
    }
    assert!(matches!(
        owner.solve_molecular_velocity_ports(&material, omega, &refs, controls, || true),
        Err(MolecularError::Boundary(HelmholtzError::Cancelled))
    ));
    let mut total = 0;
    let reference = owner
        .solve_molecular_velocity_ports(&material, omega, &refs, controls, || {
            total += 1;
            false
        })
        .unwrap();
    assert!(total > 2);
    let mut polls = 0;
    assert!(matches!(
        owner.solve_molecular_velocity_ports(&material, omega, &refs, controls, || {
            polls += 1;
            polls >= total - 1
        }),
        Err(MolecularError::Boundary(HelmholtzError::Cancelled))
    ));
    let replay = owner
        .solve_molecular_velocity_ports(&material, omega, &refs, controls, || false)
        .unwrap();
    assert_eq!(reference.impedance(), replay.impedance());
    assert_eq!(reference.pressure_jumps_pa(), replay.pressure_jumps_pa());
    assert_eq!(
        reference.maximum_backward_residual(),
        replay.maximum_backward_residual()
    );
    assert!(reference.pressure_jump(&[C64::ONE], || false).is_err());
    assert!(
        reference
            .resisting_efforts(&[C64::ONE, C64::new(f64::NAN, 0.0)])
            .is_err()
    );
    assert!(matches!(
        reference.pressure_jump(&[C64::ONE; 2], || true),
        Err(HelmholtzError::Cancelled)
    ));
    assert_eq!(
        reference
            .boundary_work_w(&[C64::ZERO; 2], || false)
            .unwrap(),
        0.0
    );
}
