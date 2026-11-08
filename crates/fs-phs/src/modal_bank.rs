//! Physical damping on the original mass-normalized modal coordinates.

use crate::linear_exact::positive_heat::factor_resistance;
use crate::{PhsError, PortHamiltonian, QuadraticStorage, ResistanceWitness};
use fs_math::det;

/// One-port mass-normalized modal bank with physical viscous rates [1/s].
///
/// Each coordinate obeys `q'' + c*q' + omega^2*q = drive*u`.
/// Zero frequency admits a free coordinate with its original drag, including
/// exactly zero drag. No damping ratio or artificial stiffness is introduced.
/// See [`modal_bank_ports_with_damping_rates`] for coordinates and admission.
pub fn modal_bank_with_damping_rates(
    omegas: &[f64],
    damping_per_s: &[f64],
    drive: &[f64],
) -> Result<PortHamiltonian, PhsError> {
    modal_bank_ports_with_damping_rates(omegas, damping_per_s, &[drive])
}

/// Mass-normalized modal bank with one physical viscous rate per coordinate.
///
/// The state is interleaved `[q_0, p_0, q_1, p_1, ...]`, with `p=q'` and
/// `H = sum(p_i^2 + omega_i^2*q_i^2)/2`. Each drive slice is one column of B
/// in `p' = -Omega^2*q - C*p + B*u`; output is the conjugate `y=B^T*p`.
/// Supply `c=b/m` for an isolated dimensional mass m with viscous drag b.
/// Frequencies, states and drives must use the same mass-normalized basis.
///
/// This replaces the entire damping law. It never adds a loss floor, rescales
/// motion or derives a damping ratio at zero frequency. The checked loss
/// factor is retained for the existing exact-flow heat and composition owners.
///
/// # Errors
/// Empty/mismatched shapes, nonfinite drives, negative/nonfinite frequencies
/// or rates, unrepresentable nonzero squared frequencies, or dense capacity.
pub fn modal_bank_ports_with_damping_rates(
    omegas: &[f64],
    damping_per_s: &[f64],
    drives: &[&[f64]],
) -> Result<PortHamiltonian, PhsError> {
    assemble(omegas, Damping::Rates(damping_per_s), drives)
}

/// One-port mass-normalized modal bank with complete coupled viscosity [1/s].
/// See [`modal_bank_ports_with_damping_matrix`] for the physical normalization,
/// state order, passivity admission and zero-frequency semantics.
pub fn modal_bank_with_damping_matrix(
    omegas: &[f64],
    damping_per_s: &[f64],
    drive: &[f64],
) -> Result<PortHamiltonian, PhsError> {
    modal_bank_ports_with_damping_matrix(omegas, damping_per_s, &[drive])
}

/// Mass-normalized modal bank with a complete symmetric PSD viscous matrix C.
///
/// `damping_per_s` is row-major in the supplied mode order. The exact original
/// entries, including either off-diagonal sign and every loss nullspace, enter
/// `p' = -Omega^2*q - C*p + B*u`. Its physical power loss is `p^T C p`.
/// State order, ports and zero-frequency coordinates follow
/// [`modal_bank_ports_with_damping_rates`]. No Rayleigh or diagonal fit is made.
///
/// The basis must satisfy `Phi^T M Phi = I`; supply `C = Phi^T C_physical Phi`
/// and the correspondingly projected physical force columns B. With diagonal
/// unnormalized modal masses, `C_ij = C_physical_ij / sqrt(m_i*m_j)` and each
/// drive entry is divided by `sqrt(m_i)`. This constructor does not infer mass
/// from frequency or silently treat a dimensional damping matrix as a rate.
///
/// Exact represented symmetry is required. The existing positive-resistance
/// factor admits the complete matrix and checks its componentwise reconstruction;
/// unresolved negative eigenvalues refuse, without clipping or deleting weak
/// losses. That factor accompanies the unchanged matrix through composition
/// and independently integrated heat. The ordinary pHS structural admission
/// still applies. These are floating-point checks, not interval certificates.
///
/// # Errors
/// The rate constructor's shape/frequency/drive errors, incomplete or nonfinite
/// C, asymmetric C, nonpassive or numerically unresolved C, or dense capacity.
pub fn modal_bank_ports_with_damping_matrix(
    omegas: &[f64],
    damping_per_s: &[f64],
    drives: &[&[f64]],
) -> Result<PortHamiltonian, PhsError> {
    assemble(omegas, Damping::Matrix(damping_per_s), drives)
}

#[derive(Clone, Copy)]
enum Damping<'a> {
    Rates(&'a [f64]),
    Matrix(&'a [f64]),
}

fn assemble(
    omegas: &[f64],
    damping: Damping<'_>,
    drives: &[&[f64]],
) -> Result<PortHamiltonian, PhsError> {
    let nm = omegas.len();
    let m = drives.len();
    let shape = || PhsError::Dimension { what: "physical modal bank extents" };
    let n = nm.checked_mul(2).ok_or_else(shape)?;
    let square = n.checked_mul(n).ok_or_else(shape)?;
    let ports = n.checked_mul(m).ok_or_else(shape)?;
    if nm == 0 || m == 0 || drives.iter().any(|drive| drive.len() != nm) {
        return Err(shape());
    }
    if omegas.iter().any(|&omega| {
        let stiffness = omega * omega;
        omega < 0.0 || !omega.is_finite() || !stiffness.is_finite()
            || (omega != 0.0 && stiffness == 0.0)
    }) || drives.iter().any(|drive| drive.iter().any(|value| !value.is_finite())) {
        return Err(PhsError::Unadmitted {
            what: "finite modal frequencies, represented stiffness and drive projections required",
        });
    }
    match damping {
        Damping::Rates(rates) => {
            if rates.len() != nm { return Err(shape()); }
            if rates.iter().any(|rate| !rate.is_finite() || *rate < 0.0) {
                return Err(PhsError::Unadmitted {
                    what: "finite nonnegative modal damping rates required",
                });
            }
        }
        Damping::Matrix(matrix) => {
            if matrix.len() != nm * nm { return Err(shape()); }
            if matrix.iter().any(|value| !value.is_finite()) {
                return Err(PhsError::Unadmitted { what: "finite modal damping matrix required" });
            }
            for i in 0..nm {
                for j in 0..i {
                    if matrix[i * nm + j] != matrix[j * nm + i] {
                        return Err(PhsError::NotSymmetric { what: "modal damping C" });
                    }
                }
            }
        }
    }
    // Q/R admission uses the existing dense owner. Refuse its work ceiling
    // before allocation rather than reaching its infallible API's assertion.
    fs_la::eigen::admit_jacobi_eigh(n).map_err(|_| PhsError::Unadmitted {
        what: "physical modal bank exceeds dense admission capacity",
    })?;
    let damping_root = match damping {
        Damping::Rates(rates) => {
            let mut root = vec![0.0; nm * nm];
            for i in 0..nm { root[i * nm + i] = det::sqrt(rates[i]); }
            root
        }
        Damping::Matrix(matrix) => factor_resistance(matrix, nm)?,
    };
    let mut q = vec![0.0; square];
    let mut j = vec![0.0; square];
    let mut r = vec![0.0; square];
    let mut root = vec![0.0; square];
    let mut g = vec![0.0; ports];
    for i in 0..nm {
        let (qi, pi) = (2 * i, 2 * i + 1);
        q[qi * n + qi] = omegas[i] * omegas[i];
        q[pi * n + pi] = 1.0;
        j[qi * n + pi] = 1.0;
        j[pi * n + qi] = -1.0;
        for k in 0..nm {
            r[pi * n + 2 * k + 1] = match damping {
                Damping::Rates(rates) => if i == k { rates[i] } else { 0.0 },
                Damping::Matrix(matrix) => matrix[i * nm + k],
            };
            root[pi * n + 2 * k + 1] = damping_root[i * nm + k];
        }
        for (port, drive) in drives.iter().enumerate() {
            g[pi * m + port] = drive[i];
        }
    }
    let witness = ResistanceWitness { matrix: r.clone(), root };
    let storage = Box::new(QuadraticStorage::new(q, n)?);
    PortHamiltonian::new(n, m, j, r, g, storage)?.with_resistance_witness(witness)
}
