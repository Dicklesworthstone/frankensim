//! Full-state Galerkin image of the continuous m=0 modal owner.
//!
//! Lossless assembly can retain original molecular relaxation through
//! [ModalGalerkin::with_molecular_relaxation], without reducing any spatial
//! mode. That attachment supplies real bulk absorption and original thermal
//! memory; it does not introduce wall viscothermal physics.
//!
//! Let C=(R'/R)A be the ORIGINAL Riccati coupling and P=q_dot. Its equations
//! give U=-(S/rho)(q'-Cq). The weak equation is M q_ddot + K q = B U_in,
//! with M=int S/(rho c²) N^T N dx and
//! K=int S/rho [(N'-CN)^T(N'-CN)+N^T diag(gamma²/R²) N] dx.
//! Thus H=(q_dot^T M q_dot+q^T K q)/2 and H_dot=P_in^T U_in.
//! Positive two-point Gauss assembly keeps this structure on the ORIGINAL
//! smooth radius/slope; it creates no staircase or fitted damping.
//!
//! Every axial/radial coordinate is retained. The pressure-release mouth
//! constrains its endpoint potential to zero; a rigid mouth is the natural
//! zero-flow boundary. Canonical state is [q; p=M q_dot], output is inlet
//! modal pressure, input is volume flow INTO the region. A port mouth retains
//! every endpoint coordinate and exposes every modal pressure/flow pair, with
//! the same inward-flow convention on both faces. Both propagating
//! and evanescent behavior follow from the same transverse stiffness.
//!
//! This is a spatial approximation of the existing finite-mode operator.
//! Axial elements and retained modes must be refined separately; quadrature
//! refines with the elements. No 3-D, viscothermal, azimuthal-mode, open-bell
//! radiation or nonlinear-interface accuracy is inferred from energy closure.

use super::riccati::{Basis, ModeCount};
use crate::{DuctError, LossModel, RadiusProfile};
use fs_la::factor::cholesky;
use fs_material::gas::GasState;
use fs_phs::linear_exact::MAX_LINEAR_EXPONENTIAL_BLOCK_DIM;
use fs_phs::{PhsError, PortHamiltonian, QuadraticStorage};

mod molecular;
pub use molecular::ModalMolecularBudget;

/// Actual boundary on every retained far-end mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ModalEnd {
    /// Zero modal flow; the natural boundary of the weak form.
    Rigid,
    /// Zero modal pressure; endpoint potentials are fixed at zero.
    PressureRelease,
    /// Expose every mouth mode as a physical port. Inputs are volume flow
    /// INTO the duct, opposite the axial flow at the mouth; outputs are modal
    /// pressure. Port order is all inlet modes, then all mouth modes.
    Ports,
}

/// Refusal from the physical duct source or the existing structural owner.
#[derive(Debug)]
pub enum ModalTimeError {
    /// Original geometry, gas, mode or cancellation refusal.
    Duct(DuctError),
    /// Original passive-structure / exact-propagation refusal.
    Structure(PhsError),
    /// Invalid dimensions or an unrepresentable assembled operator.
    Invalid(&'static str),
}
impl core::fmt::Display for ModalTimeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Duct(error) => write!(f, "modal time source: {error}"),
            Self::Structure(error) => write!(f, "modal time structure: {error}"),
            Self::Invalid(what) => write!(f, "modal time source: {what}"),
        }
    }
}
impl core::error::Error for ModalTimeError {}
impl From<DuctError> for ModalTimeError {
    fn from(error: DuctError) -> Self {
        Self::Duct(error)
    }
}
impl From<PhsError> for ModalTimeError {
    fn from(error: PhsError) -> Self {
        Self::Structure(error)
    }
}

/// Retained original source and actual consistent spatial operators.
/// Coordinates are ordered by axial node, then radial mode; no eigenmodes
/// of these matrices are selected or discarded.
#[derive(Debug, Clone)]
pub struct ModalGalerkinSource {
    profile: RadiusProfile,
    gas: GasState,
    regional_gases: Option<Vec<GasState>>,
    modes: ModeCount,
    elements: usize,
    end: ModalEnd,
    nodes: Vec<f64>,
    radial_eigenvalues: Vec<f64>,
    mass: Vec<f64>,
    stiffness: Vec<f64>,
    inverse_mass: Vec<f64>,
    inlet_reference: f64,
    mouth_reference: f64,
    molecular: Option<molecular::MolecularSource>,
}
impl ModalGalerkinSource {
    /// Original continuous radius and slope owner.
    pub const fn profile(&self) -> &RadiusProfile {
        &self.profile
    }
    /// Original admitted inlet gas. For a regional source, the other original
    /// states remain available through [`Self::gas_at_element`].
    pub const fn gas(&self) -> &GasState {
        &self.gas
    }
    /// Original gas of one axial element, without averaging its density,
    /// compressibility or source state. Out-of-range indices return None.
    pub fn gas_at_element(&self, element: usize) -> Option<&GasState> {
        if element >= self.elements {
            return None;
        }
        Some(
            self.regional_gases
                .as_ref()
                .map_or(&self.gas, |gases| &gases[element]),
        )
    }
    /// Authored per-element states, when this source was assembled regionally.
    /// The uniform constructor retains its single original state instead.
    pub fn regional_gases(&self) -> Option<&[GasState]> {
        self.regional_gases.as_deref()
    }
    /// Complete retained radial basis, including the plane mode.
    pub const fn modes(&self) -> ModeCount {
        self.modes
    }
    /// Actual number of continuous axial elements.
    pub const fn elements(&self) -> usize {
        self.elements
    }
    /// Authored boundary on every far-end mode.
    pub const fn end(&self) -> ModalEnd {
        self.end
    }
    /// Retained axial nodes [m]. The constrained release endpoint is omitted.
    pub fn nodes_m(&self) -> &[f64] {
        &self.nodes
    }
    /// Original eigenvalues gamma_n from the certified Riccati basis.
    pub fn radial_eigenvalues(&self) -> &[f64] {
        &self.radial_eigenvalues
    }
    /// Number of potential coordinates; the acoustic canonical prefix has twice
    /// this length. Molecular temperature/reservoir blocks follow that prefix.
    pub fn coordinates(&self) -> usize {
        self.nodes.len() * self.modes.0
    }
    /// Row-major consistent M in the full retained spatial basis. With molecular
    /// memory this is the frozen pressure mass, with original geometry retained.
    pub fn mass(&self) -> &[f64] {
        &self.mass
    }
    /// Row-major K, assembled from the original covariant derivative.
    pub fn stiffness(&self) -> &[f64] {
        &self.stiffness
    }
    /// Row-major M inverse; nodal pressure is M inverse times canonical p.
    pub fn inverse_mass(&self) -> &[f64] {
        &self.inverse_mass
    }
    /// Common real inlet modal wave reference rho*c/S [Pa s/m³].
    pub const fn inlet_reference_impedance(&self) -> f64 {
        self.inlet_reference
    }
    /// Common real mouth modal wave reference rho*c/S [Pa s/m³]. This uses
    /// the actual outlet area; it is not an evanescent matching impedance.
    pub const fn mouth_reference_impedance(&self) -> f64 {
        self.mouth_reference
    }
    /// Physical flow-input / pressure-output count. A port mouth contributes
    /// every retained mode, in the same order as the inlet.
    pub fn port_count(&self) -> usize {
        self.modes.0 * if self.end == ModalEnd::Ports { 2 } else { 1 }
    }
}

/// Complete original spatial system and its source description.
pub struct ModalGalerkin {
    source: ModalGalerkinSource,
    system: PortHamiltonian,
}
impl ModalGalerkin {
    /// Retained source and actual spatial matrices.
    pub const fn source(&self) -> &ModalGalerkinSource {
        &self.source
    }
    /// Full physical system; R is exactly zero before molecular attachment.
    pub const fn system(&self) -> &PortHamiltonian {
        &self.system
    }
    /// Transfer both original owners to a physical runtime without refitting.
    pub fn into_parts(self) -> (ModalGalerkinSource, PortHamiltonian) {
        (self.source, self.system)
    }
}

fn cancelled_at(
    cancelled: &mut impl FnMut() -> bool,
    completed: usize,
) -> Result<(), ModalTimeError> {
    if cancelled() {
        return Err(DuctError::Cancelled {
            completed_steps: completed,
        }
        .into());
    }
    Ok(())
}

/// Assemble the same finite-mode continuous operator as mm_profile_impedance.
/// Piecewise-linear axial functions use two positive Gauss nodes per element.
/// The complete state must fit both the caller's ceiling and the existing
/// exact-ZOH dense block capacity; resolution is never silently reduced.
///
/// The zero-frequency plane-potential gauge at a rigid or port mouth remains
/// represented. A release endpoint is a physical Dirichlet constraint. With
/// `ModalEnd::Ports`, inlet and mouth inputs are both flows INTO this region;
/// their outputs are pressures in the original radial basis.
///
/// # Errors
/// Original gas/profile/mode refusals, cancellation, zero resolution, complete
/// state capacity, nonfinite matrix terms, non-SPD mass or structural refusal.
#[allow(clippy::too_many_arguments)]
pub fn assemble_lossless(
    profile: &RadiusProfile,
    gas: &GasState,
    modes: ModeCount,
    elements: usize,
    end: ModalEnd,
    max_states: usize,
    cancelled: impl FnMut() -> bool,
) -> Result<ModalGalerkin, ModalTimeError> {
    assemble(
        profile, gas, None, modes, elements, end, max_states, cancelled,
    )
}

/// Assemble the same full multimodal operator with one original gas state per
/// uniform axial element. Interfaces lie exactly on the authored element grid;
/// refining that grid requires repeating each state's physical interval.
///
/// Each element uses its own `S/rho` and `S/(rho*c^2)` in the SAME stiffness
/// and consistent mass assembly. Continuous pressure potential and the weak
/// flux balance enforce pressure and volume-flow continuity across interfaces,
/// including density contrasts and local transverse cutoffs. No average sound
/// speed, impedance interpolation or fitted resonance correction is used.
///
/// These are prescribed frozen regions at one common static pressure. Mean
/// flow, mixing, thermal/species evolution and moving interfaces are absent.
/// Source states, geometry, every coordinate and the original work limits are
/// retained. Uniform repeated states reproduce the uniform operator exactly.
///
/// # Errors
/// Empty regions, incoherent gas states, unequal static pressures, or any
/// original geometry, resolution, capacity, factorization or cancellation
/// refusal from [`assemble_lossless`].
pub fn assemble_lossless_regional(
    profile: &RadiusProfile,
    gases: &[GasState],
    modes: ModeCount,
    end: ModalEnd,
    max_states: usize,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ModalGalerkin, ModalTimeError> {
    cancelled_at(&mut cancelled, 0)?;
    let inlet = gases.first().ok_or(ModalTimeError::Invalid(
        "regional modal source needs one gas state per axial element",
    ))?;
    assemble(
        profile,
        inlet,
        Some(gases),
        modes,
        gases.len(),
        end,
        max_states,
        cancelled,
    )
}

#[allow(clippy::too_many_arguments)]
fn assemble(
    profile: &RadiusProfile,
    gas: &GasState,
    regional_gases: Option<&[GasState]>,
    modes: ModeCount,
    elements: usize,
    end: ModalEnd,
    max_states: usize,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ModalGalerkin, ModalTimeError> {
    cancelled_at(&mut cancelled, 0)?;
    modes.validate()?;
    let n = modes.0;
    let nodes = elements
        .checked_add(usize::from(end != ModalEnd::PressureRelease))
        .ok_or(ModalTimeError::Invalid("axial node count overflow"))?;
    let coordinates = nodes
        .checked_mul(n)
        .ok_or(ModalTimeError::Invalid("modal coordinate count overflow"))?;
    let states = coordinates
        .checked_mul(2)
        .ok_or(ModalTimeError::Invalid("modal state count overflow"))?;
    let ports = n * if end == ModalEnd::Ports { 2 } else { 1 };
    let block = states
        .checked_add(ports)
        .and_then(|v| v.checked_add(1))
        .and_then(|v| v.checked_mul(2))
        .ok_or(ModalTimeError::Invalid("modal exact-flow block overflow"))?;
    if elements == 0
        || states == 0
        || states > max_states
        || block > MAX_LINEAR_EXPONENTIAL_BLOCK_DIM
    {
        return Err(ModalTimeError::Invalid(
            "complete modal state exceeds the authored/exact-flow capacity",
        ));
    }
    // Lossless gas admission is frequency independent. This query delegates
    // coherence, area and inlet impedance admission to the original owner.
    let inlet_reference =
        crate::segment_wave(gas, profile.inlet_radius(), 1.0, LossModel::Lossless)?
            .characteristic_impedance
            .re;
    if let Some(gases) = regional_gases {
        for (element, state) in gases.iter().enumerate() {
            cancelled_at(&mut cancelled, element)?;
            crate::validate_gas(state)?;
            if state.pressure != gas.pressure {
                return Err(ModalTimeError::Invalid(
                    "regional modal gases need a common static pressure; mean flow is not modeled",
                ));
            }
        }
    }
    let outlet_gas = regional_gases.map_or(gas, |gases| &gases[elements - 1]);
    let mouth_reference = crate::segment_wave(
        outlet_gas,
        profile.outlet_radius(),
        1.0,
        LossModel::Lossless,
    )?
    .characteristic_impedance
    .re;
    let h = profile.length() / elements as f64;
    if !h.is_finite() || h <= 0.0 {
        return Err(ModalTimeError::Invalid(
            "axial element width is unrepresentable",
        ));
    }
    let basis = Basis::new(n);
    let entries = coordinates * coordinates;
    let mut mass = vec![0.0; entries];
    let mut stiffness = vec![0.0; entries];
    let local = 2 * n;
    let mut gradient = vec![0.0; n * local];
    let offset = 0.5 / fs_math::det::sqrt(3.0);
    for element in 0..elements {
        cancelled_at(&mut cancelled, element)?;
        let local_gas = regional_gases.map_or(gas, |gases| &gases[element]);
        for fraction in [0.5 - offset, 0.5 + offset] {
            let x = (element as f64 + fraction) * h;
            let radius = profile.radius(x)?;
            let flare = profile.slope(x)? / radius;
            let area = core::f64::consts::PI * radius * radius;
            let kinetic_weight = 0.5 * h * area / local_gas.density;
            let compressive_weight =
                kinetic_weight / (local_gas.sound_speed * local_gas.sound_speed);
            if !kinetic_weight.is_finite()
                || kinetic_weight <= 0.0
                || !compressive_weight.is_finite()
                || compressive_weight <= 0.0
                || !flare.is_finite()
            {
                return Err(ModalTimeError::Invalid(
                    "modal quadrature weights are unrepresentable",
                ));
            }
            let shape = [1.0 - fraction, fraction];
            let derivative = [-1.0 / h, 1.0 / h];
            for m in 0..n {
                for side in 0..2 {
                    for k in 0..n {
                        gradient[m * local + side * n + k] =
                            (if m == k { derivative[side] } else { 0.0 })
                                - flare * basis.coupling[m * n + k] * shape[side];
                    }
                }
            }
            for u in 0..local {
                let row = element * n + u;
                if row >= coordinates {
                    continue;
                }
                let (side_u, mode_u) = (u / n, u % n);
                for v in 0..=u {
                    let column = element * n + v;
                    if column >= coordinates {
                        continue;
                    }
                    let (side_v, mode_v) = (v / n, v % n);
                    let mut k = 0.0;
                    for m in 0..n {
                        k += gradient[m * local + u] * gradient[m * local + v];
                    }
                    if mode_u == mode_v {
                        let transverse = basis.gammas[mode_u] / radius;
                        k += transverse * transverse * shape[side_u] * shape[side_v];
                        mass[row * coordinates + column] +=
                            compressive_weight * shape[side_u] * shape[side_v];
                    }
                    stiffness[row * coordinates + column] += kinetic_weight * k;
                }
            }
        }
    }
    for i in 0..coordinates {
        for j in 0..i {
            mass[j * coordinates + i] = mass[i * coordinates + j];
            stiffness[j * coordinates + i] = stiffness[i * coordinates + j];
        }
    }
    if mass.iter().chain(&stiffness).any(|v| !v.is_finite()) {
        return Err(ModalTimeError::Invalid(
            "assembled modal storage is nonfinite",
        ));
    }
    cancelled_at(&mut cancelled, elements)?;
    let factor = cholesky(&mass, coordinates)
        .map_err(|_| ModalTimeError::Invalid("consistent modal mass is not positive definite"))?;
    // Triangular inverse followed by a Gram product preserves the inverse
    // mass's symmetry/positive structure without eigenvalue repair.
    let mut inverse_l = vec![0.0; entries];
    for column in 0..coordinates {
        for row in 0..coordinates {
            let mut value = if row == column { 1.0 } else { 0.0 };
            for k in 0..row {
                value -= factor.l(row, k) * inverse_l[k * coordinates + column];
            }
            inverse_l[row * coordinates + column] = value / factor.l(row, row);
        }
    }
    let mut inverse_mass = vec![0.0; entries];
    for i in 0..coordinates {
        for j in 0..=i {
            let mut value = 0.0;
            for k in i..coordinates {
                value += inverse_l[k * coordinates + i] * inverse_l[k * coordinates + j];
            }
            inverse_mass[i * coordinates + j] = value;
            inverse_mass[j * coordinates + i] = value;
        }
    }
    let mut q = vec![0.0; states * states];
    let mut j = vec![0.0; states * states];
    let mut g = vec![0.0; states * ports];
    for row in 0..coordinates {
        j[row * states + coordinates + row] = 1.0;
        j[(coordinates + row) * states + row] = -1.0;
        for column in 0..coordinates {
            q[row * states + column] = stiffness[row * coordinates + column];
            q[(coordinates + row) * states + coordinates + column] =
                inverse_mass[row * coordinates + column];
        }
    }
    for mode in 0..n {
        g[(coordinates + mode) * ports + mode] = 1.0;
        if end == ModalEnd::Ports {
            // Both physical inputs point into the region. Integration by
            // parts gives +P_in*U_in + P_mouth*U_mouth in the power balance.
            g[(2 * coordinates - n + mode) * ports + n + mode] = 1.0;
        }
    }
    if q.iter().any(|v| !v.is_finite()) {
        return Err(ModalTimeError::Invalid(
            "modal canonical storage is unrepresentable",
        ));
    }
    let system = PortHamiltonian::new(
        states,
        ports,
        j,
        vec![0.0; states * states],
        g,
        Box::new(QuadraticStorage::new(q, states)?),
    )?;
    cancelled_at(&mut cancelled, elements)?;
    Ok(ModalGalerkin {
        source: ModalGalerkinSource {
            profile: profile.clone(),
            gas: *gas,
            regional_gases: regional_gases.map(<[GasState]>::to_vec),
            modes,
            elements,
            end,
            nodes: (0..nodes).map(|i| i as f64 * h).collect(),
            radial_eigenvalues: basis.gammas,
            mass,
            stiffness,
            inverse_mass,
            inlet_reference,
            mouth_reference,
            molecular: None,
        },
        system,
    })
}
