//! Execute the existing modal damping certificate in the original coordinates.
//!
//! Singleton modes keep their certified physical decay rates. Every Full block
//! keeps its complete principal damping matrix, realized by the original passive
//! channel owner. The existing coupled exact solver then advances those islands.
//! This module performs no eigensolve, projection or time integration of its own.
//!
//! The opt-in approximation omits only coupling between certified blocks. Its
//! certificate bounds free-decay poles of the original retained matrix. It does
//! not bound arbitrary forced trajectories, pressure, or heat error. The runtime
//! accounts for work and nonnegative physical loss of the REALIZED operator;
//! omitted coupling is never booked as heat or numerical dissipation.

use super::{ModalPlateViscosity, WindowedPlateViscosity, invalid, owner, poll};
use crate::acoustic_realize::AcousticRealizeError;
use crate::modal_acoustic_time::ModalAcousticTimeModel;
use crate::render::schedule::force::coupled::{
    CoupledModalSystem, MAX_MODAL_CONNECTIONS, ModalColumnConnection, ModalCouplingConfig,
    ModalPreparation,
};
use fs_exec::CancelGate;
use fs_modal::{
    CertifiedModalDamping, ModalDampingBlock, ModalDampingBudget, certify_modal_damping,
    certify_modal_damping_owned_diagonal,
};

#[derive(Clone, Debug)]
struct Cluster {
    modes: Vec<usize>,
    damping: ModalPlateViscosity,
}

/// A passive time-domain realization of one original damping certificate.
///
/// The original mode order, frequencies, force coordinates and observation
/// coordinates are retained. Diagonal entries are physical coefficients [1/s],
/// not extra loss added to the original matrix. Full clusters carry no parallel
/// scalar damping; their complete principal matrices supply all their loss.
///
/// Storage of this image is O(n + sum of squared Full-block sizes), excluding
/// the certificate's own data. In particular, the windowed constructor never
/// forms a global dense damping matrix. Native runtime construction retains its
/// original global connection layout and explicit setup limits.
#[derive(Clone, Debug)]
pub struct CertifiedPlateViscosity {
    certificate: CertifiedModalDamping,
    omega: Vec<f64>,
    diagonal_rates: Vec<f64>,
    clusters: Vec<Cluster>,
}

impl ModalPlateViscosity {
    /// Compile this SAME projected matrix into certified scalar/full blocks.
    ///
    /// When present, the original energy-quotient diagonal evaluator remains the
    /// owner of singleton decay rates. Full blocks use the original matrix with
    /// no refit or additional entry threshold. Their passive factorization and
    /// reconstruction checks belong to ModalPlateViscosity.
    ///
    /// This explicitly selects a free-pole approximation. It does not replace
    /// the unchanged full-matrix connections supplied by Self::connections.
    ///
    /// # Errors
    /// Original certificate or passive-factor refusal, cancellation, or a
    /// frequency vector that does not cover this complete original basis.
    pub fn certify_for_time(
        &self,
        omega: &[f64],
        budget: &ModalDampingBudget,
        gate: &CancelGate,
    ) -> Result<CertifiedPlateViscosity, AcousticRealizeError> {
        poll(gate)?;
        if omega.len() != self.mode_count() {
            return Err(invalid("certified viscosity needs every original modal frequency"));
        }
        let certificate = match self.diagonal_zeta() {
            Some(zeta) => certify_modal_damping_owned_diagonal(
                omega,
                self.matrix_per_s(),
                zeta,
                budget,
            ),
            None => certify_modal_damping(omega, self.matrix_per_s(), budget),
        }
        .map_err(owner)?;
        poll(gate)?;
        let n = self.mode_count();
        CertifiedPlateViscosity::compile(
            omega,
            certificate,
            self.reference_resistance_n_s_m(),
            |i, j| Some(self.matrix_per_s()[n * i + j]),
            gate,
        )
    }
}

impl WindowedPlateViscosity {
    /// Compile certified blocks from the existing near entries and far bounds.
    ///
    /// Only complete principal matrices of Full blocks are formed. Any required
    /// unread entry refuses; there is no new projection, zero fill, wider window
    /// or silently relaxed tolerance. The original far-row bound stays attached
    /// to the resulting certificate.
    ///
    /// # Errors
    /// Original windowed certificate or passive-factor refusal, cancellation,
    /// or a Full block requiring a matrix entry absent from the original window.
    pub fn certify_for_time(
        &self,
        budget: &ModalDampingBudget,
        gate: &CancelGate,
    ) -> Result<CertifiedPlateViscosity, AcousticRealizeError> {
        poll(gate)?;
        let certificate = self.certify(budget).map_err(owner)?;
        poll(gate)?;
        // One N s/m is only a channel factorization scale. It neither changes
        // the projected [1/s] coefficients nor introduces a physical support.
        CertifiedPlateViscosity::compile(
            self.omega_rad_s(),
            certificate,
            1.0,
            |i, j| self.entry(i, j),
            gate,
        )
    }
}

impl CertifiedPlateViscosity {
    fn compile(
        omega: &[f64],
        certificate: CertifiedModalDamping,
        resistance: f64,
        entry: impl Fn(usize, usize) -> Option<f64>,
        gate: &CancelGate,
    ) -> Result<Self, AcousticRealizeError> {
        let mut diagonal_rates = vec![0.0; omega.len()];
        let mut clusters = Vec::new();
        for block in &certificate.blocks {
            poll(gate)?;
            match block {
                ModalDampingBlock::Diagonal { mode, eigenvalue, .. } => {
                    // Preserve the certified centre, including the original
                    // energy-quotient owner's final product and exact zero.
                    let rate = -2.0 * eigenvalue.re;
                    if !rate.is_finite() || rate < 0.0 {
                        return Err(invalid("certified modal decay rate is unrepresentable"));
                    }
                    diagonal_rates[*mode] = rate;
                }
                ModalDampingBlock::Full { modes, .. } => {
                    let size = modes.len();
                    let mut matrix = Vec::with_capacity(size * size);
                    for &i in modes {
                        poll(gate)?;
                        for &j in modes {
                            matrix.push(entry(i, j).ok_or_else(|| {
                                invalid("certified full damping block needs an unread source entry")
                            })?);
                        }
                    }
                    let damping = ModalPlateViscosity::factor(matrix, size, resistance, gate)?;
                    clusters.push(Cluster { modes: modes.clone(), damping });
                }
            }
        }
        poll(gate)?;
        Ok(Self { certificate, omega: omega.to_vec(), diagonal_rates, clusters })
    }

    /// Original retained-matrix free-pole certificate, with its original scope.
    ///
    /// It does not certify the poles after adding another support, contact,
    /// radiation or material-memory coupling to the returned runtime.
    #[must_use]
    pub const fn certificate(&self) -> &CertifiedModalDamping {
        &self.certificate
    }

    /// Original frequencies [rad/s], in the unchanged mechanical coordinate order.
    #[must_use]
    pub fn omega_rad_s(&self) -> &[f64] {
        &self.omega
    }

    /// Physical scalar rates [1/s]. Full-cluster entries are exactly zero:
    /// their passive coupled channels supply their entire material damping.
    #[must_use]
    pub fn diagonal_rates_per_s(&self) -> &[f64] {
        &self.diagonal_rates
    }

    /// Passive channels retained across all Full clusters, including zero
    /// channels for an exactly lossless repeated-frequency cluster.
    #[must_use]
    pub fn coupled_channel_count(&self) -> usize {
        self.clusters.iter().map(|cluster| cluster.damping.channel_count()).sum()
    }

    /// Bind to the existing exact held-force owner without changing coordinates.
    ///
    /// The input must use the original basis/order and matching frequencies,
    /// with every scalar damping rate exactly zero. Frequency equality alone
    /// cannot establish the identity of a degenerate eigenspace; the caller must
    /// supply its original coordinates. States, clock, pressure transfers and
    /// modal budgets are preserved by the existing rate-binding owner.
    ///
    /// Ordinary preparation uses the original exact island decomposition and
    /// scalar fast path, so weak omitted far coupling cannot force one global
    /// exponential. Prepared preparation uses the original explicit dense
    /// budget when larger Full clusters need it. No admission is widened here.
    /// HarmonicOnly refuses because this constructor promises a time runtime.
    ///
    /// Singleton heat is reported as component dissipation and Full-block heat
    /// as connection dissipation by that owner; their sum is this realization's
    /// intrinsic loss. It adds no numerical damping. External work remains on
    /// the caller's complete original force vector.
    ///
    /// # Errors
    /// Incompatible basis/frequencies, already selected damping, cancellation,
    /// or the existing modal, connection, setup, propagation and energy gates.
    pub fn prepare_coupled(
        &self,
        model: ModalAcousticTimeModel,
        config: ModalCouplingConfig,
        preparation: ModalPreparation,
        gate: &CancelGate,
    ) -> Result<CoupledModalSystem, AcousticRealizeError> {
        poll(gate)?;
        let n = self.omega.len();
        if model.modes().len() != n
            || model.modes().iter().zip(&self.omega).any(|(mode, omega)| {
                mode.angular_frequency_rad_s.to_bits() != omega.to_bits()
            })
            || model.damping_rates_per_s().iter().any(|rate| *rate != 0.0)
        {
            return Err(invalid(
                "certified damping requires its complete original basis and zero prior scalar loss",
            ));
        }
        match preparation {
            ModalPreparation::HarmonicOnly(_) => {
                return Err(invalid("certified damping time realization needs time preparation"));
            }
            ModalPreparation::Prepared(budget) => {
                budget.validate().map_err(owner)?;
                if n > budget.max_modes {
                    return Err(invalid("certified damping exceeds the original prepared mode budget"));
                }
            }
            ModalPreparation::Ordinary => {}
        }
        let channels = self.coupled_channel_count();
        // The native owner repeats these checks. Check its allocation/work
        // limits before expanding compact cluster columns to the global layout.
        let extent = channels.checked_add(1)
            .ok_or_else(|| invalid("certified damping connection extent overflows"))?;
        let terms = n.checked_mul(extent).and_then(|x| x.checked_mul(extent))
            .ok_or_else(|| invalid("certified damping setup work overflows"))?;
        if n > config.max_modes
            || config.max_modes > fs_modal::time::MAX_TIME_DOMAIN_ACOUSTIC_MODES
            || channels > config.max_connections
            || config.max_connections > MAX_MODAL_CONNECTIONS
            || terms > config.max_setup_terms
        {
            return Err(invalid("certified damping exceeds the original modal connection/setup limits"));
        }
        let model = model.with_damping_rates(&self.diagonal_rates).map_err(owner)?;
        let mut connections = Vec::with_capacity(channels);
        for cluster in &self.clusters {
            for local in cluster.damping.channel_shapes() {
                poll(gate)?;
                let mut column = vec![0.0; n];
                for (&mode, &value) in cluster.modes.iter().zip(local) {
                    column[mode] = value;
                }
                connections.push(ModalColumnConnection {
                    column,
                    stiffness_n_m: 0.0,
                    damping_n_s_m: cluster.damping.reference_resistance_n_s_m(),
                    rest_extension_m: 0.0,
                });
            }
        }
        poll(gate)?;
        CoupledModalSystem::new_with_columns_preparation(
            vec![model], connections, &[], config, preparation, gate,
        ).map_err(owner)
    }
}

#[cfg(test)]
mod tests;
