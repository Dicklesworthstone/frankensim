//! A resistance witness lets its producer read the checked factor rows, and they rebuild the resistance.
use fs_phs::{PortHamiltonian, QuadraticStorage, ResistanceWitness};

fn identity(n: usize) -> Vec<f64> {
    let mut matrix = vec![0.0; n * n];
    for i in 0..n {
        matrix[i * n + i] = 1.0;
    }
    matrix
}

fn source_matrix(n: usize, rows: &[f64]) -> Vec<f64> {
    let mut matrix = vec![0.0; n * n];
    for row in rows.chunks_exact(n) {
        for i in 0..n {
            for j in 0..n {
                matrix[i * n + j] += row[i] * row[j];
            }
        }
    }
    matrix
}

fn system(n: usize, resistance: Vec<f64>) -> PortHamiltonian {
    PortHamiltonian::new(
        n,
        0,
        vec![0.0; n * n],
        resistance,
        Vec::new(),
        Box::new(QuadraticStorage::new(identity(n), n).unwrap()),
    )
    .unwrap()
}

/// The reconstruction gate restated independently of fs-phs: at every entry the Gram of the rows must reproduce
/// `R` within `128 n eps` of `|R_ij|` (of the local product sum where `R_ij` is a structural zero).
fn worst_gate_ratio(matrix: &[f64], rows: &[f64], n: usize) -> f64 {
    let mut worst = 0.0_f64;
    for i in 0..n {
        for j in i..n {
            let (mut gram, mut products) = (0.0, 0.0);
            for k in 0..n {
                let product = rows[k * n + i] * rows[k * n + j];
                gram += product;
                products += product.abs();
            }
            let magnitude = if matrix[i * n + j] != 0.0 { matrix[i * n + j].abs() } else { products };
            if magnitude > 0.0 {
                worst = worst.max((gram - matrix[i * n + j]).abs() / (128.0 * n as f64 * f64::EPSILON * magnitude));
            }
        }
    }
    worst
}

#[test]
fn witness_rows_rebuild_the_resistance_within_the_reconstruction_gate_and_a_bent_resistance_is_refused() {
    let n = 4;
    // Six source squares over four coordinates (taller than n, so the witness is QR-compressed), with
    // opposite-sign couplings and a column that only the last row reaches.
    let rows = [
        1.0, -0.6, 0.0, 0.0, //
        0.3, 1.2, -0.9, 0.0, //
        0.0, 0.0, 0.5, -1.1, //
        0.2, 0.0, 0.0, 0.7, //
        -0.4, 0.8, 0.0, 0.1, //
        0.0, 0.0, 0.3, 0.6,
    ];
    let matrix = source_matrix(n, &rows);
    let witness = ResistanceWitness::try_from_source_rows(n, matrix.clone(), rows.to_vec()).unwrap();
    assert_eq!(witness.rows().len(), n * n);
    let ratio = worst_gate_ratio(&matrix, witness.rows(), n);
    assert!(ratio <= 1.0, "rows rebuild R outside the gate: {ratio:e}");
    // The rows are a square root of R in the physical sense: power sum_k (row_k . effort)^2 = effort^T R effort.
    let effort = [0.7, -1.3, 0.4, 2.1];
    let power: f64 = witness.rows().chunks_exact(n).map(|r| r.iter().zip(&effort).map(|(a, b)| a * b).sum::<f64>().powi(2)).sum();
    let reference: f64 = (0..n).map(|i| (0..n).map(|j| effort[i] * matrix[i * n + j] * effort[j]).sum::<f64>()).sum();
    assert!((power - reference).abs() <= 1e-13 * reference);

    // A system exposes the witness it was given, bit for bit, and none before.
    let bare = system(n, matrix.clone());
    assert!(bare.resistance_witness().is_none());
    let witnessed = bare.with_resistance_witness(witness.clone()).unwrap();
    let exposed = witnessed.resistance_witness().unwrap().rows();
    assert_eq!(
        exposed.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
        witness.rows().iter().map(|x| x.to_bits()).collect::<Vec<_>>()
    );

    // Planted negative: R bent by 1e-9 relative in one diagonal entry is still a valid resistance, but it is not the
    // one the witness was built for, and the rows cannot be rebuilt into it.
    let mut bent = matrix.clone();
    bent[0] *= 1.0 + 1e-9;
    assert!(worst_gate_ratio(&bent, witness.rows(), n) > 1.0);
    assert!(system(n, bent.clone()).with_resistance_witness(witness.clone()).is_err());
    assert!(ResistanceWitness::try_from_source_rows(n, bent, rows.to_vec()).is_err());
}
