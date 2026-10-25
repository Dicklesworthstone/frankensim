//! Physical state transfer between modes on one fixed reference mass space.
//!
//! Shapes, not sorted mode numbers, carry displacement and conjugate momentum.
//! The target basis may permute, change sign, or rotate a degenerate subspace.
//! A changed support/mesh/mass requires a different transfer owner. All storage
//! here is incremental modal storage about the supplied equilibria; it is not
//! the complete prestressed body's thermodynamic energy or a heat debit.

use crate::{ModePair, time::ModalAcousticState};
use fs_exec::Cx;
use fs_sparse::{Csr, DirectOrdering, LdltError, LdltOptions, SymbolicLdlt};

/// Explicit accuracy, state-loss, work and allocation limits for one transfer.
#[derive(Clone, Copy, Debug)]
pub struct ModalTransferControls {
    /// Maximum number of coordinates in the common reference space.
    pub maximum_dofs: usize,
    /// Maximum number of modes in either source basis.
    pub maximum_modes: usize,
    /// Applied separately to each original sparse factorization.
    pub maximum_factor_bytes: usize,
    /// Absolute error permitted in every entry of Phi^T M Phi - I.
    pub maximum_mass_orthogonality_error: f64,
    /// Relative Euclidean residual of K phi = lambda M phi in the target.
    pub maximum_relative_eigen_residual: f64,
    /// Maximum discarded displacement in the mass norm [m sqrt(kg)].
    pub maximum_displacement_residual_m_sqrt_kg: f64,
    /// Maximum kinetic energy in the discarded velocity field [J].
    pub maximum_kinetic_residual_j: f64,
    /// Maximum target elastic energy in the discarded displacement field [J].
    pub maximum_elastic_residual_j: f64,
    /// Signed change of represented incremental storage BEFORE projection.
    pub maximum_abs_parameter_work_j: f64,
    /// Error allowed in full target energy = retained energy + residual energy.
    pub energy_tolerance_j: f64,
}

/// An immutable proposal. Projection residuals are discarded represented state,
/// never physical damping, heat, or a permission to normalize the waveform.
#[derive(Clone, Debug)]
pub struct ModalTransfer {
    /// Mass-normalized coordinates and velocities in the target basis.
    pub states: Vec<ModalAcousticState>,
    /// Storage represented by the original modes and state [J].
    pub previous_modal_energy_j: f64,
    /// Storage of the complete old physical state in the target pencil [J].
    pub target_full_energy_j: f64,
    /// Storage represented by the projected target modal state [J].
    pub retained_modal_energy_j: f64,
    /// Target full energy minus previous modal energy. This is the external
    /// parameter owner's change of incremental storage, not thermodynamic heat.
    pub parameter_work_j: f64,
    /// Mass norm of the discarded displacement field [m sqrt(kg)].
    pub displacement_residual_m_sqrt_kg: f64,
    /// Actual discarded velocity energy under the original mass operator [J].
    pub kinetic_residual_j: f64,
    /// Actual discarded displacement energy under the target stiffness [J].
    pub elastic_residual_j: f64,
    /// Full target storage minus retained and discarded storage [J].
    pub projection_energy_residual_j: f64,
}

/// No accepted motion is changed by a transfer refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum ModalTransferError {
    /// Unusable source, state, coordinate or control data.
    Invalid(&'static str),
    /// Cancellation observed before publication.
    Cancelled,
    /// Original sparse admission/factorization failure.
    Factor(LdltError),
    /// An explicit count or allocation budget was exceeded.
    Budget {
        quantity: &'static str,
        requested: usize,
        cap: usize,
    },
    /// A physical or numerical transfer limit was exceeded.
    Limit {
        quantity: &'static str,
        value: f64,
        limit: f64,
    },
}

impl core::fmt::Display for ModalTransferError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-MODAL-TRANSFER: {self:?}")
    }
}
impl std::error::Error for ModalTransferError {}
impl From<LdltError> for ModalTransferError {
    fn from(error: LdltError) -> Self { Self::Factor(error) }
}

fn poll(cx: &Cx<'_>) -> Result<(), ModalTransferError> {
    cx.checkpoint().map_err(|_| ModalTransferError::Cancelled)
}
fn finite(value: f64) -> Result<f64, ModalTransferError> {
    if value.is_finite() { Ok(value) }
    else { Err(ModalTransferError::Invalid("unrepresentable transfer arithmetic")) }
}
fn limit(quantity: &'static str, value: f64, cap: f64) -> Result<(), ModalTransferError> {
    finite(value)?;
    if value > cap { Err(ModalTransferError::Limit { quantity, value, limit: cap }) }
    else { Ok(()) }
}
fn dot(a: &[f64], b: &[f64], cx: &Cx<'_>) -> Result<f64, ModalTransferError> {
    let mut sum = 0.0;
    for (a, b) in a.chunks(256).zip(b.chunks(256)) {
        poll(cx)?;
        for (a, b) in a.iter().zip(b) { sum += a * b; }
    }
    finite(sum)
}
fn apply(a: &Csr, x: &[f64], cx: &Cx<'_>) -> Result<Vec<f64>, ModalTransferError> {
    poll(cx)?;
    let mut image = vec![0.0; a.nrows()];
    a.spmv(x, &mut image);
    poll(cx)?;
    if image.iter().all(|x| x.is_finite()) { Ok(image) }
    else { Err(ModalTransferError::Invalid("unrepresentable sparse operator application")) }
}
fn energy(a: &Csr, x: &[f64], cx: &Cx<'_>) -> Result<f64, ModalTransferError> {
    let value = finite(0.5 * dot(x, &apply(a, x, cx)?, cx)?)?;
    if value >= 0.0 { Ok(value) }
    else { Err(ModalTransferError::Invalid("negative quadratic storage")) }
}
fn admit_spd(a: &Csr, controls: ModalTransferControls, cx: &Cx<'_>)
    -> Result<(), ModalTransferError>
{
    for row in 0..a.nrows() {
        poll(cx)?;
        let (columns, values) = a.row(row);
        for (&column, &value) in columns.iter().zip(values) {
            if !value.is_finite() || value != a.get(column, row) {
                return Err(ModalTransferError::Invalid(
                    "finite exactly symmetric mass and target stiffness required"));
            }
        }
    }
    let symbolic = SymbolicLdlt::analyze(a, DirectOrdering::Amd)?;
    poll(cx)?;
    let factor = symbolic.factor(a, &LdltOptions {
        max_factor_bytes: Some(controls.maximum_factor_bytes),
        ..LdltOptions::default()
    })?;
    poll(cx)?;
    if factor.inertia().positive != a.nrows() {
        return Err(ModalTransferError::Invalid(
            "positive definite reference mass and stable target stiffness required"));
    }
    Ok(())
}

fn admit_basis(
    mass: &Csr, stiffness: Option<&Csr>, modes: &[ModePair],
    controls: ModalTransferControls, cx: &Cx<'_>,
) -> Result<(), ModalTransferError> {
    for (i, mode) in modes.iter().enumerate() {
        poll(cx)?;
        if mode.phi.len() != mass.nrows() || !mode.lambda.is_finite() || mode.lambda <= 0.0
            || !mode.phi.iter().all(|x| x.is_finite()) {
            return Err(ModalTransferError::Invalid(
                "finite positive-frequency shapes on the complete reference DOF space required"));
        }
        let m_phi = apply(mass, &mode.phi, cx)?;
        for (j, other) in modes[..=i].iter().enumerate() {
            limit("mass orthogonality",
                (dot(&other.phi, &m_phi, cx)? - if i == j { 1.0 } else { 0.0 }).abs(),
                controls.maximum_mass_orthogonality_error)?;
        }
        if let Some(stiffness) = stiffness {
            let k_phi = apply(stiffness, &mode.phi, cx)?;
            let residual: Vec<_> = k_phi.iter().zip(&m_phi)
                .map(|(k, m)| k - mode.lambda * m).collect();
            let scale = finite(dot(&k_phi, &k_phi, cx)?.sqrt()
                + mode.lambda * dot(&m_phi, &m_phi, cx)?.sqrt())?;
            if scale <= 0.0 {
                return Err(ModalTransferError::Invalid("target eigenvector scale is not representable"));
            }
            limit("target eigen residual", dot(&residual, &residual, cx)?.sqrt() / scale,
                controls.maximum_relative_eigen_residual)?;
        }
    }
    Ok(())
}

fn modal_energy(modes: &[ModePair], states: &[ModalAcousticState], cx: &Cx<'_>)
    -> Result<f64, ModalTransferError>
{
    let mut result = 0.0;
    for (mode, state) in modes.iter().zip(states) {
        poll(cx)?;
        let q = state.displacement_m_sqrt_kg;
        let v = state.velocity_m_sqrt_kg_per_s;
        result += 0.5 * (v * v + mode.lambda * q * q);
    }
    finite(result)
}

/// Transfer one accepted physical state to a new modal basis on the same mass.
///
/// `equilibrium_offset_m` is old equilibrium minus new equilibrium in the
/// common reduced-DOF ordering. The source owner must bind both bases, this
/// offset, the previous modal frequencies and the unchanged reference mass to
/// their actual source states. Target eigen residuals are checked independently
/// against the supplied target stiffness; no mode index correspondence is used.
///
/// With mass-orthonormal Phi, the proposal is
/// `q_new = Phi_new^T M (Phi_old q_old + offset)` and
/// `v_new = Phi_new^T M Phi_old v_old`. Since p = M v on this fixed mass
/// space, the second equation projects physical canonical momentum. The full
/// target-pencil energy is computed before projection and then reconciled with
/// the target modal energy plus independently computed residual-field energies.
///
/// Positive-frequency stable pencils and exactly symmetric sparse matrices are
/// required. Changing mass, supports, mesh, rigid modes or constitutive memory
/// requires another admitted transfer. This operation neither advances time
/// nor mutates any source, state or accepted runtime on failure or success.
#[allow(clippy::too_many_arguments)]
pub fn transfer_modal_state(
    mass: &Csr,
    target_stiffness: &Csr,
    previous_modes: &[ModePair],
    target_modes: &[ModePair],
    previous_states: &[ModalAcousticState],
    equilibrium_offset_m: &[f64],
    controls: ModalTransferControls,
    cx: &Cx<'_>,
) -> Result<ModalTransfer, ModalTransferError> {
    poll(cx)?;
    let tolerances = [controls.maximum_mass_orthogonality_error,
        controls.maximum_relative_eigen_residual,
        controls.maximum_displacement_residual_m_sqrt_kg,
        controls.maximum_kinetic_residual_j, controls.maximum_elastic_residual_j,
        controls.maximum_abs_parameter_work_j, controls.energy_tolerance_j];
    if controls.maximum_dofs == 0 || controls.maximum_modes == 0
        || controls.maximum_factor_bytes == 0
        || tolerances.iter().any(|x| !x.is_finite() || *x < 0.0)
        || controls.maximum_mass_orthogonality_error >= 1.0
        || controls.maximum_relative_eigen_residual >= 1.0 {
        return Err(ModalTransferError::Invalid("finite nonnegative transfer limits and nonzero count/allocation budgets required"));
    }
    let n = mass.nrows();
    if n == 0 || mass.ncols() != n || target_stiffness.nrows() != n
        || target_stiffness.ncols() != n || equilibrium_offset_m.len() != n
        || previous_modes.is_empty() || target_modes.is_empty()
        || previous_modes.len() != previous_states.len()
        || previous_modes.len() > n || target_modes.len() > n {
        return Err(ModalTransferError::Invalid("matching nonempty reference space, bases, states and equilibrium offset required"));
    }
    for (quantity, requested, cap) in [
        ("reference DOFs", n, controls.maximum_dofs),
        ("previous modes", previous_modes.len(), controls.maximum_modes),
        ("target modes", target_modes.len(), controls.maximum_modes),
    ] {
        if requested > cap { return Err(ModalTransferError::Budget { quantity, requested, cap }); }
    }
    if equilibrium_offset_m.iter().any(|x| !x.is_finite())
        || previous_states.iter().any(|state| !state.displacement_m_sqrt_kg.is_finite()
            || !state.velocity_m_sqrt_kg_per_s.is_finite()) {
        return Err(ModalTransferError::Invalid("finite previous motion and equilibrium offset required"));
    }
    admit_spd(mass, controls, cx)?;
    admit_spd(target_stiffness, controls, cx)?;
    admit_basis(mass, None, previous_modes, controls, cx)?;
    admit_basis(mass, Some(target_stiffness), target_modes, controls, cx)?;

    let previous_modal_energy_j = modal_energy(previous_modes, previous_states, cx)?;
    let mut u = equilibrium_offset_m.to_vec();
    let mut v = vec![0.0; n];
    for (mode, state) in previous_modes.iter().zip(previous_states) {
        for ((shape, u), v) in mode.phi.chunks(256).zip(u.chunks_mut(256)).zip(v.chunks_mut(256)) {
            poll(cx)?;
            for ((&shape, u), v) in shape.iter().zip(u).zip(v) {
                *u += shape * state.displacement_m_sqrt_kg;
                *v += shape * state.velocity_m_sqrt_kg_per_s;
            }
        }
    }
    let m_u = apply(mass, &u, cx)?;
    let m_v = apply(mass, &v, cx)?;
    let target_full_energy_j = finite(energy(mass, &v, cx)? + energy(target_stiffness, &u, cx)?)?;
    let parameter_work_j = finite(target_full_energy_j - previous_modal_energy_j)?;
    limit("absolute parameter work", parameter_work_j.abs(), controls.maximum_abs_parameter_work_j)?;

    let mut states = Vec::with_capacity(target_modes.len());
    for mode in target_modes {
        states.push(ModalAcousticState {
            displacement_m_sqrt_kg: dot(&mode.phi, &m_u, cx)?,
            velocity_m_sqrt_kg_per_s: dot(&mode.phi, &m_v, cx)?,
        });
    }
    let mut du = u;
    let mut dv = v;
    for (mode, state) in target_modes.iter().zip(&states) {
        for ((shape, du), dv) in mode.phi.chunks(256).zip(du.chunks_mut(256)).zip(dv.chunks_mut(256)) {
            poll(cx)?;
            for ((&shape, du), dv) in shape.iter().zip(du).zip(dv) {
                *du -= shape * state.displacement_m_sqrt_kg;
                *dv -= shape * state.velocity_m_sqrt_kg_per_s;
            }
        }
    }
    let displacement_residual_m_sqrt_kg = finite((2.0 * energy(mass, &du, cx)?).sqrt())?;
    let kinetic_residual_j = energy(mass, &dv, cx)?;
    let elastic_residual_j = energy(target_stiffness, &du, cx)?;
    let retained_modal_energy_j = modal_energy(target_modes, &states, cx)?;
    let projection_energy_residual_j = finite(target_full_energy_j
        - retained_modal_energy_j - kinetic_residual_j - elastic_residual_j)?;
    limit("displacement projection residual", displacement_residual_m_sqrt_kg,
        controls.maximum_displacement_residual_m_sqrt_kg)?;
    limit("kinetic projection residual", kinetic_residual_j, controls.maximum_kinetic_residual_j)?;
    limit("elastic projection residual", elastic_residual_j, controls.maximum_elastic_residual_j)?;
    limit("projection energy residual", projection_energy_residual_j.abs(), controls.energy_tolerance_j)?;
    poll(cx)?;
    Ok(ModalTransfer { states, previous_modal_energy_j, target_full_energy_j,
        retained_modal_energy_j, parameter_work_j, displacement_residual_m_sqrt_kg,
        kinetic_residual_j, elastic_residual_j, projection_energy_residual_j })
}
