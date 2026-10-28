//! Physical state transfer between modes on one reference coordinate space.
//!
//! Shapes, not sorted mode numbers, carry displacement and conjugate momentum.
//! The target basis may permute, change sign, or rotate a degenerate subspace.
//! Changed inertia preserves canonical momentum, not velocity. A changed
//! support/mesh/coordinate map requires a different transfer owner. All storage
//! here is incremental modal storage about the supplied equilibria; it is not
//! the complete prestressed body's thermodynamic energy or a heat debit.

use crate::{ModePair, time::ModalAcousticState};
use fs_exec::Cx;
use fs_sparse::{Coo, Csr, DirectOrdering, LdltError, LdltFactor, LdltOptions, SymbolicLdlt};

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

/// Limits for a change of inertia on the same physical coordinate space.
#[derive(Clone, Copy, Debug)]
pub struct ModalMassTransferControls {
    /// State projection, parameter work, and sparse factorization limits.
    pub transfer: ModalTransferControls,
    /// Maximum relative residual of `M_new v_new = M_old v_old` before
    /// projection: `||M_new v_new - p||_2 / (||M_new v_new||_2 + ||p||_2)`.
    /// Exact zero momentum has zero residual.
    pub maximum_relative_momentum_residual: f64,
}

/// Explicit admission of a target's complete rigid-motion nullspace.
#[derive(Clone, Copy, Debug)]
pub struct ModalRigidTransferControls {
    /// Original projection, work, inertia-solve and factorization budgets.
    pub transfer: ModalMassTransferControls,
    /// Maximum staged entries of the bordered stiffness qualification matrix.
    /// Its conservative count is nnz(K) + 2 * dofs * rigid_modes.
    pub maximum_border_entries: usize,
    /// Relative null residual ||K z||_2 / (||K||_F ||z||_2). A zero K has
    /// exactly zero residual. This does not turn a nonzero eigenvalue into zero.
    pub maximum_relative_rigid_residual: f64,
}

/// A transfer under changed inertia, including the independent mass-solve check.
#[derive(Clone, Debug)]
pub struct ModalMassTransfer {
    /// Projected motion and the complete incremental storage/work accounting.
    pub transfer: ModalTransfer,
    /// Relative momentum residual before any declared modal truncation.
    pub relative_momentum_residual: f64,
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
    /// Actual discarded velocity energy under the target mass operator [J].
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
fn admit_spd(
    a: &Csr,
    controls: ModalTransferControls,
    cx: &Cx<'_>,
) -> Result<LdltFactor, ModalTransferError> {
    for row in 0..a.nrows() {
        poll(cx)?;
        let (columns, values) = a.row(row);
        for (&column, &value) in columns.iter().zip(values) {
            if !value.is_finite() || value != a.get(column, row) {
                return Err(ModalTransferError::Invalid(
                    "finite exactly symmetric mass and target stiffness required",
                ));
            }
        }
    }
    let symbolic = SymbolicLdlt::analyze(a, DirectOrdering::Amd)?;
    poll(cx)?;
    let factor = symbolic.factor(
        a,
        &LdltOptions {
            max_factor_bytes: Some(controls.maximum_factor_bytes),
            ..LdltOptions::default()
        },
    )?;
    poll(cx)?;
    if factor.inertia().positive != a.nrows() {
        return Err(ModalTransferError::Invalid(
            "positive definite reference mass and stable target stiffness required",
        ));
    }
    Ok(factor)
}

fn admit_basis(
    mass: &Csr, stiffness: Option<&Csr>, modes: &[ModePair],
    controls: ModalTransferControls, allow_rigid: bool, cx: &Cx<'_>,
) -> Result<(), ModalTransferError> {
    for (i, mode) in modes.iter().enumerate() {
        poll(cx)?;
        if mode.phi.len() != mass.nrows() || !mode.lambda.is_finite() || mode.lambda < 0.0
            || !allow_rigid && mode.lambda == 0.0
            || !mode.phi.iter().all(|x| x.is_finite()) {
            return Err(ModalTransferError::Invalid(if allow_rigid {
                "finite nonnegative-frequency shapes on the complete reference DOF space required"
            } else {
                "finite positive-frequency shapes on the complete reference DOF space required"
            }));
        }
        let m_phi = apply(mass, &mode.phi, cx)?;
        for (j, other) in modes[..=i].iter().enumerate() {
            limit("mass orthogonality",
                (dot(&other.phi, &m_phi, cx)? - if i == j { 1.0 } else { 0.0 }).abs(),
                controls.maximum_mass_orthogonality_error)?;
        }
        if let Some(stiffness) = stiffness.filter(|_| mode.lambda != 0.0) {
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

// K is positive definite on the M-orthogonal complement of Z iff this
// bordered matrix has inertia (n positive, r negative). The explicit KZ test
// admits the alleged rigid space; no diagonal shift or artificial spring is
// applied to physical dynamics, residual energies or returned eigenvalues.
fn admit_rigid_stiffness(
    mass: &Csr, stiffness: &Csr, modes: &[ModePair],
    controls: ModalRigidTransferControls, cx: &Cx<'_>,
) -> Result<(), ModalTransferError> {
    let rigid: Vec<_> = modes.iter().filter(|m| m.lambda == 0.0).collect();
    if rigid.is_empty() {
        admit_spd(stiffness, controls.transfer.transfer, cx)?;
        return Ok(());
    }
    let n = mass.nrows();
    let count = rigid.len();
    let entries = n.checked_mul(count).and_then(|v| v.checked_mul(2))
        .and_then(|v| v.checked_add(stiffness.nnz()))
        .ok_or(ModalTransferError::Invalid("rigid border extent overflow"))?;
    if entries > controls.maximum_border_entries {
        return Err(ModalTransferError::Budget { quantity: "rigid border entries",
            requested: entries, cap: controls.maximum_border_entries });
    }
    let dim = n.checked_add(count)
        .ok_or(ModalTransferError::Invalid("rigid border dimension overflow"))?;
    let mut norm_k = 0.0_f64;
    let mut scale_k = 0.0_f64;
    for row in 0..n {
        poll(cx)?;
        let (columns, values) = stiffness.row(row);
        for (&col, &value) in columns.iter().zip(values) {
            if !value.is_finite() || value != stiffness.get(col,row) {
                return Err(ModalTransferError::Invalid("rigid stiffness must be finite and exactly symmetric"));
            }
            norm_k = norm_k.hypot(value);
            scale_k = scale_k.max(value.abs());
        }
    }
    finite(norm_k)?;
    // Congruence scaling is confined to the qualification factorization.
    // It prevents the units/magnitude of K and MZ from setting a false pivot
    // floor; the original physical K is used by every transfer energy check.
    if scale_k == 0.0 { scale_k = 1.0; }
    let mut border = Coo::new(dim,dim);
    for row in 0..n {
        poll(cx)?;
        let (columns, values) = stiffness.row(row);
        for (&col, &value) in columns.iter().zip(values) {
            border.push(row,col,value/scale_k);
        }
    }
    for (column,mode) in rigid.iter().enumerate() {
        poll(cx)?;
        let kz = apply(stiffness,&mode.phi,cx)?;
        let norm_z = mode.phi.iter().fold(0.0_f64,|norm,v| norm.hypot(*v));
        let norm_kz = kz.iter().fold(0.0_f64,|norm,v| norm.hypot(*v));
        let relative = if norm_k == 0.0 { 0.0 } else { (norm_kz/norm_k)/norm_z };
        limit("target rigid residual",relative,controls.maximum_relative_rigid_residual)?;
        let mz = apply(mass,&mode.phi,cx)?;
        let scale = mz.iter().fold(0.0_f64,|scale,v| scale.max(v.abs()));
        if scale <= 0.0 || !scale.is_finite() {
            return Err(ModalTransferError::Invalid("rigid constraint has no representable mass image"));
        }
        for (row,value) in mz.iter().enumerate() {
            if row%256 == 0 { poll(cx)?; }
            border.push(row,n+column,value/scale);
            border.push(n+column,row,value/scale);
        }
    }
    poll(cx)?;
    let border = border.assemble();
    let symbolic = SymbolicLdlt::analyze(&border,DirectOrdering::Amd)?;
    poll(cx)?;
    let factor = symbolic.factor(&border,&LdltOptions {
        max_factor_bytes: Some(controls.transfer.transfer.maximum_factor_bytes),
        ..LdltOptions::default()
    })?;
    poll(cx)?;
    if factor.inertia().positive != n || factor.inertia().negative != count {
        return Err(ModalTransferError::Invalid("target stiffness is not stable outside its complete rigid subspace"));
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
    transfer_modal_state_impl(
        mass,
        None,
        target_stiffness,
        previous_modes,
        target_modes,
        previous_states,
        equilibrium_offset_m,
        controls,
        None,
        cx,
    )
    .map(|result| result.transfer)
}

/// Transfer vibration when an explicitly parameterized inertia changes.
///
/// Both SPD masses use the same physical DOFs and inertial frame. In the
/// absence of an external impulse, displacement and canonical momentum are
/// continuous: `u_new = u_old`, `p_new = M_old v_old`, and
/// `v_new = M_new^-1 p_new`. Target modal velocity is therefore
/// `Phi_new^T p_new`, while displacement uses the target mass inner product.
/// This matters even for one mode: retaining velocity under changed inertia
/// would create or destroy momentum.
///
/// The target mass is factored once for both SPD admission and its solve. The
/// proposal separately checks the solve residual and the state lost by modal
/// truncation. Parameter work includes the signed kinetic-storage change from
/// the inertia update and the elastic-storage change; the parameter owner
/// must supply or absorb it. No amplitude normalization or heat is inferred.
///
/// This describes a Hamiltonian parameter update, not accreted/ejected mass,
/// fluid added-mass dynamics, moving coordinates, remeshing or new constraints.
/// Those require their own momentum flux, impulse, or coordinate transfer.
/// Like [`transfer_modal_state`], this function does not advance time or
/// mutate an accepted state, and refuses unstable/zero-frequency pencils.
#[allow(clippy::too_many_arguments)]
pub fn transfer_modal_state_with_mass_change(
    previous_mass: &Csr,
    target_mass: &Csr,
    target_stiffness: &Csr,
    previous_modes: &[ModePair],
    target_modes: &[ModePair],
    previous_states: &[ModalAcousticState],
    equilibrium_offset_m: &[f64],
    controls: ModalMassTransferControls,
    cx: &Cx<'_>,
) -> Result<ModalMassTransfer, ModalTransferError> {
    transfer_modal_state_impl(
        previous_mass,
        Some((target_mass, controls.maximum_relative_momentum_residual)),
        target_stiffness,
        previous_modes,
        target_modes,
        previous_states,
        equilibrium_offset_m,
        controls.transfer,
        None,
        cx,
    )
}

/// Transfer vibration and rigid motion on the same physical coordinate space.
///
/// This explicitly extends [`transfer_modal_state_with_mass_change`] to
/// zero-frequency modes. Supply EVERY target rigid direction as an exactly
/// zero `ModePair::lambda`, together with the target oscillatory modes. Shapes
/// must be M-orthonormal; their K-null residual has its own explicit tolerance.
/// The original sparse LDLT owner qualifies K on their complement using
/// `[K, M Z; Z^T M, 0]`. Missing rigid directions, false null modes and unstable
/// elastic directions refuse. No spring floor, eigenvalue clipping, waveform
/// rescaling or coordinate-based mode matching is introduced.
///
/// Displacement and canonical momentum transfer, changed-inertia work and all
/// truncation limits are the same as the original owner. Pass the same matrix
/// for both masses when inertia is unchanged. A rigid coordinate contributes
/// kinetic storage and zero elastic storage; all motion participates
/// in the original displacement and kinetic residual allowances.
/// This is an instantaneous fixed-frame parameter change, not remeshing,
/// support impulses, finite rotations or a constitutive-memory transfer.
#[allow(clippy::too_many_arguments)]
pub fn transfer_modal_state_with_rigid_modes(
    previous_mass: &Csr,
    target_mass: &Csr,
    target_stiffness: &Csr,
    previous_modes: &[ModePair],
    target_modes: &[ModePair],
    previous_states: &[ModalAcousticState],
    equilibrium_offset_m: &[f64],
    controls: ModalRigidTransferControls,
    cx: &Cx<'_>,
) -> Result<ModalMassTransfer, ModalTransferError> {
    if controls.maximum_border_entries == 0
        || !controls.maximum_relative_rigid_residual.is_finite()
        || !(0.0..1.0).contains(&controls.maximum_relative_rigid_residual) {
        return Err(ModalTransferError::Invalid("positive border capacity and finite rigid residual tolerance in [0,1) required"));
    }
    transfer_modal_state_impl(previous_mass,
        Some((target_mass,controls.transfer.maximum_relative_momentum_residual)),
        target_stiffness,previous_modes,target_modes,previous_states,
        equilibrium_offset_m,controls.transfer.transfer,Some(controls),cx)
}

#[allow(clippy::too_many_arguments)]
fn transfer_modal_state_impl(
    previous_mass: &Csr,
    changed_mass: Option<(&Csr, f64)>,
    target_stiffness: &Csr,
    previous_modes: &[ModePair],
    target_modes: &[ModePair],
    previous_states: &[ModalAcousticState],
    equilibrium_offset_m: &[f64],
    controls: ModalTransferControls,
    rigid: Option<ModalRigidTransferControls>,
    cx: &Cx<'_>,
) -> Result<ModalMassTransfer, ModalTransferError> {
    poll(cx)?;
    let mass = changed_mass.map_or(previous_mass, |(target, _)| target);
    if let Some((_, tolerance)) = changed_mass {
        if !tolerance.is_finite() || !(0.0..1.0).contains(&tolerance) {
            return Err(ModalTransferError::Invalid(
                "finite momentum residual tolerance in [0, 1) required",
            ));
        }
    }
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
    if n == 0
        || mass.ncols() != n
        || previous_mass.nrows() != n
        || previous_mass.ncols() != n
        || target_stiffness.nrows() != n
        || target_stiffness.ncols() != n
        || equilibrium_offset_m.len() != n
        || previous_modes.is_empty()
        || target_modes.is_empty()
        || previous_modes.len() != previous_states.len()
        || previous_modes.len() > n
        || target_modes.len() > n
    {
        return Err(ModalTransferError::Invalid(
            "matching nonempty reference space, bases, states and equilibrium offset required",
        ));
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
    let mass_factor = if changed_mass.is_some() {
        admit_spd(previous_mass, controls, cx)?;
        Some(admit_spd(mass, controls, cx)?)
    } else {
        admit_spd(mass, controls, cx)?;
        None
    };
    if rigid.is_none() { admit_spd(target_stiffness, controls, cx)?; }
    admit_basis(previous_mass, None, previous_modes, controls, rigid.is_some(), cx)?;
    admit_basis(mass, Some(target_stiffness), target_modes, controls, rigid.is_some(), cx)?;
    if let Some(rigid) = rigid {
        admit_rigid_stiffness(mass,target_stiffness,target_modes,rigid,cx)?;
    }

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
    let m_v = apply(previous_mass, &v, cx)?;
    let mut relative_momentum_residual = 0.0;
    if let (Some(factor), Some((_, tolerance))) = (mass_factor, changed_mass) {
        poll(cx)?;
        v = factor.solve(&m_v);
        poll(cx)?;
        let solved_momentum = apply(mass, &v, cx)?;
        let residual: Vec<_> = solved_momentum
            .iter()
            .zip(&m_v)
            .map(|(actual, expected)| actual - expected)
            .collect();
        let scale = finite(
            dot(&solved_momentum, &solved_momentum, cx)?.sqrt() + dot(&m_v, &m_v, cx)?.sqrt(),
        )?;
        let norm = dot(&residual, &residual, cx)?.sqrt();
        relative_momentum_residual = if scale == 0.0 { 0.0 } else { norm / scale };
        limit(
            "momentum solve residual",
            relative_momentum_residual,
            tolerance,
        )?;
    }
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
    Ok(ModalMassTransfer {
        transfer: ModalTransfer {
            states,
            previous_modal_energy_j,
            target_full_energy_j,
            retained_modal_energy_j,
            parameter_work_j,
            displacement_residual_m_sqrt_kg,
            kinetic_residual_j,
            elastic_residual_j,
            projection_energy_residual_j,
        },
        relative_momentum_residual,
    })
}
