//! The joint loaded modal basis: the modes of several structures and air columns, each in its own
//! mass-normalized coordinates, coupled through the added mass of the air they share, solved as one
//! generalized eigenproblem `(K, M + M_a)` and run as one [`ModalStructure`].
//!
//! Inputs: the unloaded angular frequencies `w_k` of `n` modes in mass-normalized coordinates (unit
//! modal mass each, so `K = diag(w_k^2)` and `M = I`), and the symmetric added-mass matrix `M_a` in
//! the same coordinates (kg per kg of modal mass, so without unit). `M + M_a` must be positive
//! definite: fs-modal's dense pencil solve ([`fs_modal::eigh_gen_dense`]) factors it by Cholesky
//! and the basis is refused when that fails. The loaded modes `Phi` satisfy
//! `Phi^T (M + M_a) Phi = I` and `Phi^T K Phi = diag(w'_j^2)`, ascending; the unloaded coordinates
//! are `x = Phi y`.
//!
//! The loss of the coupled set (each mode's own loss on the diagonal, the radiation resistance
//! between modes off it) is one symmetric matrix `C` in the unloaded coordinates [1/s]. In the
//! loaded coordinates it is `Phi^T C Phi` ([`LoadedModalBasis::congruent`]), which is not diagonal
//! in general, so the structure is run undamped and the loss goes to a resistor N-port joined to
//! its modal port by a common-flow junction (fs-system's `elements::Resistance`, which refuses a
//! matrix that is not positive definite, so the definiteness of the whole loss is checked on the
//! sum).
//!
//! Held loss: the structure steps each mode exactly under a force held over the step, and a
//! resistor fed the step's mean velocity returns the force for that mean. For one underdamped mode
//! at `w` with loss `c` the step's determinant (the squared modulus of its eigenvalues) is
//! `(1 - e) / (1 + e)`, with `e = (c dt / 2) sinc^2(w dt / 2)`: the mode's energy envelope decays
//! by that factor per step, while within one step the energy ratio depends on the phase, between
//! `((1 - e) / (1 + e))^2` and 1. The envelope decays at `sinc^2(w dt / 2)` times the continuous
//! rate `c` to first order in `e` (exactly, at that times `atanh(e) / e`);
//! [`LoadedModalBasis::held_loss_factors`] gives the factor for each loaded mode, to be stated with
//! a render, not corrected.
//!
//! No claim: `M_a` and `C` are constants, so a frequency-dependent loading enters at the frequency
//! the caller sampled it at; the basis is dense, `O(n^3)` to build.

use fs_modal::eigh_gen_dense;
use fs_ports::{PortKind, StableId};

use crate::component::{
    ModalAcousticTimeBudget, ModalAttachment, ModalStructure, ModalStructureRefusal,
};

/// The loaded modes of `n` coupled mass-normalized modes.
#[derive(Debug, Clone, PartialEq)]
pub struct LoadedModalBasis {
    n: usize,
    /// Loaded angular frequencies, ascending [rad/s].
    angular_frequencies_rad_s: Vec<f64>,
    /// fs-modal's certified residual of each loaded eigenvalue `w'^2` [1/s^2].
    residuals: Vec<f64>,
    /// `Phi`, row-major: entry `k * n + j` is unloaded coordinate `k` of loaded mode `j`.
    basis: Vec<f64>,
}

impl LoadedModalBasis {
    /// Solve `(diag(w_k^2), I + M_a)` for the unloaded `angular_frequencies_rad_s` and the
    /// row-major `n x n` added mass `added_mass`.
    ///
    /// # Errors
    /// [`ModalStructureRefusal::NoModes`] for no mode; [`ModalStructureRefusal::Frequency`] for an
    /// unloaded frequency that is not positive and finite, or a loaded one (in the loaded order)
    /// whose square came back not positive and finite; [`ModalStructureRefusal::LoadedSize`] for
    /// more modes than fs-la's dense symmetric eigensolver admits, before the matrix is read;
    /// [`ModalStructureRefusal::MatrixShape`], [`ModalStructureRefusal::NonFiniteEntry`] and
    /// [`ModalStructureRefusal::Asymmetric`] for an added mass that is not a finite symmetric
    /// `n x n` matrix; [`ModalStructureRefusal::LoadedMass`] when `I + M_a` is not positive
    /// definite.
    pub fn new(
        angular_frequencies_rad_s: &[f64],
        added_mass: &[f64],
    ) -> Result<Self, ModalStructureRefusal> {
        let n = angular_frequencies_rad_s.len();
        if n == 0 {
            return Err(ModalStructureRefusal::NoModes);
        }
        for (mode, &value) in angular_frequencies_rad_s.iter().enumerate() {
            if !(value > 0.0 && value.is_finite()) {
                return Err(ModalStructureRefusal::Frequency { mode, value });
            }
        }
        // fs-modal's dense solve runs fs-la's Jacobi eigensolver, which panics past its work cap.
        fs_la::eigen::admit_jacobi_eigh(n).map_err(ModalStructureRefusal::LoadedSize)?;
        symmetric(added_mass, n)?;
        let mut stiffness = vec![0.0; n * n];
        let mut mass = added_mass.to_vec();
        for (k, &w) in angular_frequencies_rad_s.iter().enumerate() {
            stiffness[k * n + k] = w * w;
            mass[k * n + k] += 1.0;
        }
        let pairs =
            eigh_gen_dense(&stiffness, &mass, n).map_err(ModalStructureRefusal::LoadedMass)?;
        let mut basis = vec![0.0; n * n];
        for (j, pair) in pairs.iter().enumerate() {
            if !(pair.lambda > 0.0 && pair.lambda.is_finite()) {
                return Err(ModalStructureRefusal::Frequency {
                    mode: j,
                    value: pair.lambda,
                });
            }
            for (k, &value) in pair.phi.iter().enumerate() {
                basis[k * n + j] = value;
            }
        }
        Ok(Self {
            n,
            angular_frequencies_rad_s: pairs.iter().map(|pair| pair.lambda.sqrt()).collect(),
            residuals: pairs.iter().map(|pair| pair.residual).collect(),
            basis,
        })
    }

    /// Number of modes.
    #[must_use]
    pub fn modes(&self) -> usize {
        self.n
    }

    /// Loaded angular frequencies `w'_j`, ascending [rad/s].
    #[must_use]
    pub fn angular_frequencies_rad_s(&self) -> &[f64] {
        &self.angular_frequencies_rad_s
    }

    /// fs-modal's residual `|K phi - w'^2 (M + M_a) phi|` in the inverse-mass norm for each loaded
    /// mode: a bound on the distance from `w'_j^2` to a true eigenvalue of the pencil [1/s^2].
    #[must_use]
    pub fn residuals(&self) -> &[f64] {
        &self.residuals
    }

    /// `Phi`, row-major: entry `k * n + j` is unloaded coordinate `k` of loaded mode `j`.
    #[must_use]
    pub fn basis(&self) -> &[f64] {
        &self.basis
    }

    /// `Phi^T A Phi` for a symmetric row-major `n x n` matrix `A` in the unloaded coordinates, such
    /// as the loss `C`. Each entry above the diagonal is computed once and mirrored, so the result
    /// is symmetric bit for bit.
    ///
    /// # Errors
    /// [`ModalStructureRefusal::MatrixShape`], [`ModalStructureRefusal::NonFiniteEntry`] and
    /// [`ModalStructureRefusal::Asymmetric`] for a matrix that is not finite, symmetric and
    /// `n x n`.
    pub fn congruent(&self, matrix: &[f64]) -> Result<Vec<f64>, ModalStructureRefusal> {
        let n = self.n;
        symmetric(matrix, n)?;
        // `A Phi`, row-major.
        let mut product = vec![0.0; n * n];
        for k in 0..n {
            for j in 0..n {
                product[k * n + j] = (0..n)
                    .map(|l| matrix[k * n + l] * self.basis[l * n + j])
                    .sum();
            }
        }
        let mut out = vec![0.0; n * n];
        for i in 0..n {
            for j in i..n {
                let value: f64 = (0..n)
                    .map(|k| self.basis[k * n + i] * product[k * n + j])
                    .sum();
                out[i * n + j] = value;
                out[j * n + i] = value;
            }
        }
        Ok(out)
    }

    /// An attachment written in the unloaded coordinates, carried to the loaded ones: row `r`'s
    /// weights become `s'_rj = sum_k s_rk Phi_kj`, so the generalized force on loaded mode `j` is
    /// `Phi^T` of the unloaded forces and the row's flow is unchanged.
    ///
    /// # Errors
    /// [`ModalStructureRefusal::Weights`] (attachment index 0) when the weights do not hold
    /// `width` rows of `n` entries, [`ModalStructureRefusal::NonFiniteWeight`] when one is not
    /// finite, given or carried.
    pub fn loaded_attachment(
        &self,
        attachment: ModalAttachment,
    ) -> Result<ModalAttachment, ModalStructureRefusal> {
        self.carried(0, attachment)
    }

    fn carried(
        &self,
        index: usize,
        attachment: ModalAttachment,
    ) -> Result<ModalAttachment, ModalStructureRefusal> {
        let n = self.n;
        let expected = attachment.width.get() * n;
        if attachment.weights.len() != expected {
            return Err(ModalStructureRefusal::Weights {
                attachment: index,
                expected,
                found: attachment.weights.len(),
            });
        }
        if attachment.weights.iter().any(|weight| !weight.is_finite()) {
            return Err(ModalStructureRefusal::NonFiniteWeight { attachment: index });
        }
        let mut weights = vec![0.0; expected];
        for (row, out) in attachment
            .weights
            .chunks_exact(n)
            .zip(weights.chunks_exact_mut(n))
        {
            for (j, slot) in out.iter_mut().enumerate() {
                *slot = (0..n).map(|k| row[k] * self.basis[k * n + j]).sum();
            }
        }
        // Finite weights can overflow when carried.
        if weights.iter().any(|weight| !weight.is_finite()) {
            return Err(ModalStructureRefusal::NonFiniteWeight { attachment: index });
        }
        Ok(ModalAttachment {
            weights,
            ..attachment
        })
    }

    /// The unloaded coordinates `x = Phi y` of the loaded coordinates `loaded` (displacements or
    /// velocities).
    ///
    /// # Panics
    /// If `loaded` does not hold one entry per mode (programmer error).
    #[must_use]
    pub fn unloaded(&self, loaded: &[f64]) -> Vec<f64> {
        let n = self.n;
        assert_eq!(loaded.len(), n, "one loaded coordinate per mode");
        (0..n)
            .map(|k| (0..n).map(|j| self.basis[k * n + j] * loaded[j]).sum())
            .collect()
    }

    /// For each loaded mode, `sinc^2(w'_j dt_s / 2)`: the held scheme's decay rate of the energy
    /// envelope over the continuous rate, to first order in the loss, at the step `dt_s`. Every
    /// factor is NaN for a step that is not positive and finite.
    #[must_use]
    pub fn held_loss_factors(&self, dt_s: f64) -> Vec<f64> {
        if !(dt_s > 0.0 && dt_s.is_finite()) {
            return vec![f64::NAN; self.n];
        }
        self.angular_frequencies_rad_s
            .iter()
            .map(|&w| {
                let x = 0.5 * w * dt_s;
                let sinc = if x > 0.0 { x.sin() / x } else { 1.0 };
                sinc * sinc
            })
            .collect()
    }

    /// The loaded modes as one undamped [`ModalStructure`]: the modal port `modal` of kind `kind`
    /// first (width `n`, identity weights: the loaded coordinates themselves, where a resistor
    /// N-port with `Phi^T C Phi` joins), then `attachments`, each written in the unloaded
    /// coordinates and carried to the loaded ones, then `heat`.
    ///
    /// # Errors
    /// The refusals of [`Self::loaded_attachment`] (with each attachment's index in
    /// `attachments`, a carried weight that is not finite included), and every refusal of
    /// [`ModalStructure::from_modes`].
    pub fn structure(
        &self,
        sample_rate_hz: u32,
        budget: ModalAcousticTimeBudget,
        modal: StableId,
        kind: PortKind,
        attachments: Vec<ModalAttachment>,
    ) -> Result<ModalStructure, ModalStructureRefusal> {
        let n = self.n;
        let mut identity = vec![0.0; n * n];
        for k in 0..n {
            identity[k * n + k] = 1.0;
        }
        let mut all = Vec::with_capacity(attachments.len() + 1);
        all.push(ModalAttachment {
            id: modal,
            kind,
            width: core::num::NonZeroUsize::new(n).expect("a basis has at least one mode"),
            weights: identity,
        });
        for (index, attachment) in attachments.into_iter().enumerate() {
            all.push(self.carried(index, attachment)?);
        }
        let modes: Vec<(f64, f64)> = self
            .angular_frequencies_rad_s
            .iter()
            .map(|&w| (w, 0.0))
            .collect();
        ModalStructure::from_modes(sample_rate_hz, &modes, budget, all)
    }
}

/// A finite row-major `n x n` matrix whose entries equal their mirrors (a signed zero equals its
/// mirror; every later use reads one triangle or the matrix's product).
#[expect(
    clippy::float_cmp,
    reason = "symmetry is tested exactly: an entry must equal its mirror, not lie near it"
)]
fn symmetric(matrix: &[f64], n: usize) -> Result<(), ModalStructureRefusal> {
    if matrix.len() != n * n {
        return Err(ModalStructureRefusal::MatrixShape {
            expected: n * n,
            found: matrix.len(),
        });
    }
    for row in 0..n {
        for column in 0..n {
            if !matrix[row * n + column].is_finite() {
                return Err(ModalStructureRefusal::NonFiniteEntry { row, column });
            }
            // The mirror entry, in an earlier row, was found finite.
            if column < row && matrix[row * n + column] != matrix[column * n + row] {
                return Err(ModalStructureRefusal::Asymmetric { row, column });
            }
        }
    }
    Ok(())
}
