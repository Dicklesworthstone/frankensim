//! G1 rotated-tensor heat modes and G0 admission/balance/atomicity.
use super::*;
use crate::thermal::sources::{
    Slab, ThermalSourceError, ThermalSourceLimits, source_step_tensor_mapped,
};
use crate::thermal::tensor::{DiffusivityField, DiffusivityTensor};
use crate::{ThermalLbm, d3q19::ThermalLbm3};
use fs_math::det::{cos, ln, sin, sqrt};

fn mode<const Q: usize>(populations: &[[f64; Q]], basis: &[[f64; 2]]) -> [f64; 2] {
    let mut projection = [0.0; 2];
    for (pop, wave) in populations.iter().zip(basis) {
        let temperature = pop.iter().sum::<f64>() - 1.0;
        for axis in 0..2 {
            projection[axis] += temperature * wave[axis];
        }
    }
    projection
}

fn tensor_wave<const D: usize, const Q: usize>(
    lattice: &AxisLattice<D, Q>,
    shape: [usize; D],
    modes: [i32; D],
    alpha: [[f64; D]; D],
    velocity: [f64; D],
) -> (f64, f64) {
    let tensor = DiffusivityTensor::from_lattice(alpha).unwrap();
    let k: [f64; D] = core::array::from_fn(|axis| {
        core::f64::consts::TAU * f64::from(modes[axis]) / shape[axis] as f64
    });
    let cells: usize = shape.iter().product();
    let basis: Vec<_> = (0..cells)
        .map(|mut cell| {
            let mut phase = 0.0;
            for axis in 0..D {
                phase += k[axis] * (cell % shape[axis]) as f64;
                cell /= shape[axis];
            }
            [cos(phase), sin(phase)]
        })
        .collect();
    let mut populations: Vec<_> = basis
        .iter()
        .map(|wave| lattice.equilibrium(1.0 + 0.1 * wave[0], velocity))
        .collect();
    let mut post = vec![0.0; cells * Q];
    let mut history = FluxHistory::default();
    let mut initial = [0.0; 2];
    let warmup = 32;
    let steps = 256;
    for time in 0..warmup + steps {
        if time == warmup {
            initial = mode(&populations, &basis);
        }
        lattice.collide_tensor(
            &populations,
            &mut post,
            0..cells,
            0.9,
            DiffusivityField::Uniform(&tensor),
            &mut history,
            |_| velocity,
        );
        lattice.stream(&mut populations, &post, shape, None);
    }
    let final_mode = mode(&populations, &basis);
    let norm = initial[0] * initial[0] + initial[1] * initial[1];
    let ratio = [
        (final_mode[0] * initial[0] + final_mode[1] * initial[1]) / norm,
        (final_mode[1] * initial[0] - final_mode[0] * initial[1]) / norm,
    ];
    let amplitude = sqrt(ratio[0] * ratio[0] + ratio[1] * ratio[1]);
    // Independent continuum Fourier symbol, including both cross terms.
    let mut decay = 0.0;
    for row in 0..D {
        for column in 0..D {
            decay += k[row] * alpha[row][column] * k[column];
        }
    }
    let phase = steps as f64 * k.iter().zip(velocity).map(|(k, u)| k * u).sum::<f64>();
    let phase_error = (ratio[1] * cos(phase) - ratio[0] * sin(phase)) / amplitude;
    let heat = populations.iter().flatten().sum::<f64>() / cells as f64;
    assert!((heat - 1.0).abs() < 1e-11, "periodic heat drift: {heat}");
    (-ln(amplitude) / (steps as f64 * decay) - 1.0, phase_error)
}

#[test]
fn g1_rotated_tensor_diffusion_refines_with_and_without_advection() {
    for moving in [false, true] {
        let u = if moving { [0.10, 0.06] } else { [0.0; 2] };
        let v = [u[0], 0.0, u[1]];
        let mut previous = [1.0_f64; 2];
        for n in [24, 48] {
            let results = [
                tensor_wave(&D2Q5, [n, n], [1, 1], [[0.12, 0.04], [0.04, 0.06]], u),
                tensor_wave(
                    &D3Q7,
                    [n, 1, n],
                    [1, 0, 1],
                    [[0.12, 0.0, 0.04], [0.0, 0.09, 0.0], [0.04, 0.0, 0.06]],
                    v,
                ),
            ];
            for (dimension, (error, phase)) in results.into_iter().enumerate() {
                assert!(error.abs() < 0.008, "D={} n={n}: {error}", dimension + 2);
                assert!(
                    phase.abs() < 0.025,
                    "D={} n={n}: phase {phase}",
                    dimension + 2
                );
                if n == 48 {
                    assert!(
                        error.abs() < 0.3 * previous[dimension],
                        "diffusion did not refine"
                    );
                    assert!(phase.abs() < 0.004);
                }
                previous[dimension] = error.abs();
            }
        }
    }
    // Change the sign of the off-diagonal coefficient and wave direction
    // together. The physical decay is unchanged, not a diagonal approximation.
    let reflected = tensor_wave(
        &D2Q5,
        [48, 48],
        [1, -1],
        [[0.12, -0.04], [-0.04, 0.06]],
        [0.10, -0.06],
    );
    let original = tensor_wave(
        &D2Q5,
        [48, 48],
        [1, 1],
        [[0.12, 0.04], [0.04, 0.06]],
        [0.10, 0.06],
    );
    assert!((reflected.0 - original.0).abs() < 2e-10);
    assert!((reflected.1 - original.1).abs() < 2e-10);
}

#[test]
fn g0_tensor_admission_preserves_axes_and_refuses_unresolved_or_nonpassive_laws() {
    let a = [[0.12, -0.04], [-0.04, 0.06]];
    let tensor = DiffusivityTensor::from_lattice(a).unwrap();
    assert_eq!(tensor.lattice_diffusivity(), a);
    for bad in [
        [[0.0, 0.0], [0.0, 0.1]],
        [[0.1, 0.2], [0.2, 0.1]],
        [[0.1, 0.02], [0.01, 0.1]],
        [[f64::NAN, 0.0], [0.0, 0.1]],
        [[f64::INFINITY, 0.0], [0.0, 0.1]],
        [[f64::MIN_POSITIVE, 0.0], [0.0, 0.1]],
        [[f64::MAX, 0.0], [0.0, 0.1]],
    ] {
        assert!(DiffusivityTensor::from_lattice(bad).is_err(), "{bad:?}");
    }
    assert!(DiffusivityTensor::<1>::from_lattice([[0.1]]).is_err());
    assert!(
        DiffusivityTensor::<3>::from_lattice([
            [0.12, 0.02, 0.01],
            [0.02, 0.08, -0.015],
            [0.01, -0.015, 0.06],
        ])
        .is_ok()
    );
}

fn limits(cells: usize) -> ThermalSourceLimits {
    ThermalSourceLimits {
        max_cells: cells,
        temperature_range: [-1.0, 2.0],
        balance_tolerance: 1e-10,
    }
}

fn transaction<S: Slab<D, Q>, const D: usize, const Q: usize>(
    mut actual: S,
    mut reference: S,
    tensor: DiffusivityTensor<D>,
) {
    let (cells, plane) = actual.layout().unwrap();
    let per_cell = vec![tensor; cells];
    let heat: Vec<_> = (0..cells)
        .map(|cell| {
            if cell < plane || cell >= cells - plane {
                0.0
            } else {
                0.0001 * ((cell % 5) as f64 - 2.0)
            }
        })
        .collect();
    let forces: Vec<[f64; D]> =
        vec![core::array::from_fn(|axis| if axis == 0 { 1e-6 } else { 0.0 }); cells];
    for _ in 0..4 {
        let a = source_step_tensor_mapped::<S, D, Q, _>(
            &mut actual,
            &heat,
            &forces,
            limits(cells),
            DiffusivityField::Uniform(&tensor),
            |_, report| Ok(report),
        )
        .unwrap();
        let b = source_step_tensor_mapped::<S, D, Q, _>(
            &mut reference,
            &heat,
            &forces,
            limits(cells),
            DiffusivityField::PerCell(&per_cell),
            |_, report| Ok(report),
        )
        .unwrap();
        assert_eq!(a, b);
        assert!(a.defect.abs() < 1e-10);
        assert_eq!(actual.populations(), reference.populations());
    }
    let before = actual.populations().to_vec();
    let rates = actual.momentum_relaxation().to_vec();
    let failure = source_step_tensor_mapped::<S, D, Q, ()>(
        &mut actual,
        &heat,
        &forces,
        limits(cells),
        DiffusivityField::Uniform(&tensor),
        |_, _| {
            Err(ThermalSourceError::Invalid(
                "reject after both kicks and wall balance",
            ))
        },
    );
    assert!(failure.is_err());
    assert_eq!(actual.populations(), before);
    assert_eq!(actual.momentum_relaxation(), rates);
    let too_many = vec![tensor; cells + 1];
    for field in [&per_cell[..cells - 1], too_many.as_slice()] {
        assert!(
            source_step_tensor_mapped::<S, D, Q, _>(
                &mut actual,
                &heat,
                &forces,
                limits(cells),
                DiffusivityField::PerCell(field),
                |_, report| Ok(report),
            )
            .is_err()
        );
    }
    // A rejected tensor candidate must preserve retained derivative history.
    for _ in 0..3 {
        source_step_tensor_mapped::<S, D, Q, _>(
            &mut actual,
            &heat,
            &forces,
            limits(cells),
            DiffusivityField::Uniform(&tensor),
            |_, report| Ok(report),
        )
        .unwrap();
        source_step_tensor_mapped::<S, D, Q, _>(
            &mut reference,
            &heat,
            &forces,
            limits(cells),
            DiffusivityField::PerCell(&per_cell),
            |_, report| Ok(report),
        )
        .unwrap();
        assert_eq!(actual.populations(), reference.populations());
    }
    // A later scalar interval uses its retained background coefficients.
    actual.advance(&forces).unwrap();
    reference.advance(&forces).unwrap();
    assert_eq!(actual.populations(), reference.populations());
    assert_eq!(actual.momentum_relaxation(), rates);
}

#[test]
fn g0_tensor_steps_share_source_wall_balance_and_atomic_history_in_both_dimensions() {
    let tensor = DiffusivityTensor::from_lattice([[0.12, 0.04], [0.04, 0.06]]).unwrap();
    transaction::<_, 2, 5>(
        ThermalLbm::slab(5, 6, 0.8, 0.9, 1e-5),
        ThermalLbm::slab(5, 6, 0.8, 0.9, 1e-5),
        tensor,
    );
    let tensor = DiffusivityTensor::from_lattice([
        [0.12, 0.02, 0.01],
        [0.02, 0.08, -0.015],
        [0.01, -0.015, 0.06],
    ])
    .unwrap();
    transaction::<_, 3, 7>(
        ThermalLbm3::slab(5, 3, 6, 0.8, 0.9, 1e-5),
        ThermalLbm3::slab(5, 3, 6, 0.8, 0.9, 1e-5),
        tensor,
    );
}

#[test]
fn g0_public_tensor_steps_admit_spatial_laws_and_keep_background_rates() {
    let mut two = ThermalLbm::slab(5, 6, 0.8, 0.9, 0.0);
    let n = two.gpop.len();
    let field: Vec<_> = (0..n)
        .map(|cell| {
            let cross = if cell % 2 == 0 { 0.04 } else { -0.04 };
            DiffusivityTensor::from_lattice([[0.12, cross], [cross, 0.06]]).unwrap()
        })
        .collect();
    for _ in 0..5 {
        let report = two
            .step_with_sources_and_tensor_diffusion(
                &vec![0.0; n],
                &vec![[0.0; 2]; n],
                DiffusivityField::PerCell(&field),
                limits(n),
            )
            .unwrap();
        assert!(report.defect.abs() < 1e-10);
    }
    assert_eq!(two.tau_g, 0.9);
    assert_eq!(two.grid.tau, vec![0.8; n]);

    let mut three = ThermalLbm3::slab(5, 3, 6, 0.8, 0.9, 0.0);
    let n = three.gpop.len();
    let tensor = DiffusivityTensor::from_lattice([
        [0.12, 0.02, 0.01],
        [0.02, 0.08, -0.015],
        [0.01, -0.015, 0.06],
    ])
    .unwrap();
    let mut heat = vec![0.0; n];
    for value in &mut heat[15..n - 15] {
        *value = 0.0001;
    }
    let report = three
        .step_with_sources_and_tensor_diffusion(
            &heat,
            &vec![[0.0; 3]; n],
            DiffusivityField::Uniform(&tensor),
            limits(n),
        )
        .unwrap();
    assert!((report.source - 0.0001 * (n - 30) as f64).abs() < 1e-15);
    assert!(report.defect.abs() < 1e-10);
    assert_eq!(three.tau_g, 0.9);
    assert_eq!(three.grid.tau, vec![0.8; n]);
}
