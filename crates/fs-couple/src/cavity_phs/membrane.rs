//! Sealed source gas and the original finite-deformation P1 membrane.
//!
//! The membrane owns swept volume, current-area traction and its geometric
//! derivative. The parcel owns the isentrope, pressure, calorics and acoustic
//! compressibility. Their common potential counts gas work once. This is a
//! frozen-composition equilibrium connection, not a heat or time integrator.

mod solve;
pub use solve::{GasMembraneSolution, GasMembraneSolve};
mod vibration;
pub use vibration::GasMembraneVibration;

use fs_material::thermochemical_gas::energy::entropy::{GasEntropyError, GasEntropySolve};
use fs_material::thermochemical_gas::energy::parcel::mechanics::{
    GasPressureResponse, MAX_VOLUME_COORDINATES,
};
use fs_material::thermochemical_gas::energy::parcel::{ClosedGasParcel, GasVolumeStep};
use fs_material::thermochemical_gas::nonideal::energy::parcel::ClosedRealGasParcel;
use fs_solid::membrane::vibration::VibrationError;
use fs_solid::membrane::{MembraneError, MembranePatch};
use fs_solver::{NewtonError, NewtonReport, SolverRunProgress};

/// Original closed material inventory; neither variant substitutes a gas law.
#[derive(Clone, Debug, PartialEq)]
pub enum MembraneGas {
    /// Frozen ideal mixture with its original temperature-dependent calorics.
    Ideal(ClosedGasParcel),
    /// Original homogeneous Peng--Robinson branch, including caloric departures.
    PengRobinson(ClosedRealGasParcel),
}

impl From<ClosedGasParcel> for MembraneGas {
    fn from(gas: ClosedGasParcel) -> Self {
        Self::Ideal(gas)
    }
}
impl From<ClosedRealGasParcel> for MembraneGas {
    fn from(gas: ClosedRealGasParcel) -> Self {
        Self::PengRobinson(gas)
    }
}

impl MembraneGas {
    /// Current physical volume [m3], retained by the original parcel.
    #[must_use]
    pub fn volume_m3(&self) -> f64 {
        match self {
            Self::Ideal(gas) => gas.volume_m3(),
            Self::PengRobinson(gas) => gas.volume_m3(),
        }
    }

    /// Conserved source inventory [kg].
    #[must_use]
    pub fn mass_kg(&self) -> f64 {
        match self {
            Self::Ideal(gas) => gas.mass_kg(),
            Self::PengRobinson(gas) => gas.mass_kg(),
        }
    }

    /// Same-state absolute pressure and adiabatic compressibility.
    pub fn pressure_response(&self) -> Result<GasPressureResponse, GasEntropyError> {
        match self {
            Self::Ideal(gas) => gas.pressure_response(),
            Self::PengRobinson(gas) => gas.pressure_response(),
        }
    }

    fn at_volume(
        &self,
        volume_m3: f64,
        controls: GasEntropySolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Self, GasVolumeStep), GasMembraneError> {
        let mut candidate = self.clone();
        let step = match &mut candidate {
            Self::Ideal(gas) => gas.set_volume_isentropic(volume_m3, controls, cancelled)?,
            Self::PengRobinson(gas) => {
                gas.set_volume_isentropic(volume_m3, controls, cancelled)?
            }
        };
        Ok((candidate, step))
    }

    fn forces(&self, gradient: &[f64], exterior: f64) -> Result<Vec<f64>, GasMembraneError> {
        Ok(match self {
            Self::Ideal(gas) => gas
                .mechanical_load(gas.volume_m3(), gradient, exterior)?
                .generalized_forces()?,
            Self::PengRobinson(gas) => gas
                .mechanical_load(gas.volume_m3(), gradient, exterior)?
                .generalized_forces()?,
        })
    }

    fn tangent(
        &self,
        gradient: &[f64],
        exterior: f64,
        direction: &[f64],
        hessian_action: &[f64],
    ) -> Result<Vec<f64>, GasMembraneError> {
        Ok(match self {
            Self::Ideal(gas) => gas
                .mechanical_load(gas.volume_m3(), gradient, exterior)?
                .tangent_action(direction, hessian_action)?,
            Self::PengRobinson(gas) => gas
                .mechanical_load(gas.volume_m3(), gradient, exterior)?
                .tangent_action(direction, hessian_action)?,
        })
    }
}

/// A refusal publishes neither a new gas state nor a partial equilibrium.
#[derive(Debug)]
pub enum GasMembraneError {
    /// Invalid geometry, controls or unrepresentable physical coordinates.
    Invalid(&'static str),
    /// The original membrane material or geometry refused a trial.
    Membrane(MembraneError),
    /// The original material isentrope or pressure response refused a trial.
    Gas(GasEntropyError),
    /// The original nonlinear solver refused its input.
    Newton(NewtonError),
    /// The original membrane modal owner refused the coupled spectrum.
    Vibration(VibrationError),
    /// The native solve exhausted convergence or globalization.
    Unconverged(NewtonReport),
    /// The native solver stopped without an accepted equilibrium.
    Interrupted(SolverRunProgress),
    /// Combined residual/tangent/final-check evaluation budget exhausted.
    EvaluationBudget,
    /// Cancellation before publication.
    Cancelled,
}
impl core::fmt::Display for GasMembraneError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "sealed gas/membrane: {self:?}")
    }
}
impl std::error::Error for GasMembraneError {}
impl From<MembraneError> for GasMembraneError {
    fn from(error: MembraneError) -> Self {
        match error {
            MembraneError::Cancelled => Self::Cancelled,
            other => Self::Membrane(other),
        }
    }
}
impl From<GasEntropyError> for GasMembraneError {
    fn from(error: GasEntropyError) -> Self {
        match error {
            GasEntropyError::Cancelled => Self::Cancelled,
            other => Self::Gas(other),
        }
    }
}
impl From<VibrationError> for GasMembraneError {
    fn from(error: VibrationError) -> Self {
        match error {
            VibrationError::Membrane(error) => error.into(),
            other => Self::Vibration(other),
        }
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), GasMembraneError> {
    if cancelled() {
        Err(GasMembraneError::Cancelled)
    } else {
        Ok(())
    }
}
fn finite_vectors(values: &[[f64; 3]], n: usize) -> bool {
    values.len() == n && values.iter().flatten().all(|v| v.is_finite())
}

/// An immutable reference gas and a fixed-perimeter membrane cavity chart.
///
/// The supplied gas volume is the actual cavity volume at `reference_m`.
/// Positive membrane orientation points OUT of this gas; the remaining cavity
/// wall is rigid. Absolute exterior pressure is explicit and held constant.
/// The native perimeter check is mandatory even when the initial pressure jump
/// is zero. There must be 1..=256 free Cartesian coordinates; this is the
/// original gas mechanical-load bound, not a membrane mesh-resolution claim.
///
/// A positive volume interval is caller-owned geometric/source admission. It
/// does not prove global injectivity, prevent contact with the rigid enclosure,
/// or resolve nonuniform gas motion. The original membrane and homogeneous
/// frozen-composition gas validity domains remain binding at every trial.
pub struct GasMembrane<'a> {
    patch: &'a MembranePatch,
    gas: MembraneGas,
    reference_m: Vec<[f64; 3]>,
    reference_strain_energy_j: f64,
    fixed: Vec<Option<f64>>,
    free_dofs: Vec<usize>,
    exterior_pressure_pa: f64,
    volume_interval_m3: [f64; 2],
}

/// Same-geometry force, source-gas state and conservative energy evaluation.
/// This is not by itself an equilibrium, stability or time-step certificate.
#[derive(Clone, Debug)]
pub struct GasMembraneEvaluation {
    /// Original reference-node ordering, in metres.
    pub displacement_m: Vec<[f64; 3]>,
    /// Complete source-bearing gas at the actual enclosed volume.
    pub gas: MembraneGas,
    /// Original entropy and boundary-work account from the reference parcel.
    pub gas_volume: GasVolumeStep,
    /// Absolute pressure and same-state acoustic bulk modulus.
    pub pressure: GasPressureResponse,
    /// Internal force minus dead and net-pressure loads [N]. Prescribed rows
    /// are support reactions; all physical pressure loads are retained there.
    pub residual_n: Vec<[f64; 3]>,
    /// Exact fixed-perimeter swept volume from the reference geometry [m3].
    pub swept_volume_m3: f64,
    /// Original membrane stored energy [J].
    pub strain_energy_j: f64,
    /// Held dead loads dotted with the actual displacement change [J].
    pub dead_load_work_j: f64,
    /// Work against the explicit constant exterior pressure [J].
    pub exterior_pressure_work_j: f64,
    /// Delta membrane storage + native gas Delta U + exterior work - dead work.
    /// Gas pressure work is counted through Delta U, never again as p_final*dV.
    pub potential_energy_change_j: f64,
}

struct Trial {
    value: GasMembraneEvaluation,
    gradient: Vec<f64>,
}

impl<'a> GasMembrane<'a> {
    /// Bind the original gas to an admitted reference shape and fixed map.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        patch: &'a MembranePatch,
        gas: impl Into<MembraneGas>,
        reference_m: &[[f64; 3]],
        fixed: &[Option<f64>],
        exterior_pressure_pa: f64,
        volume_interval_m3: [f64; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, GasMembraneError> {
        poll(&mut cancelled)?;
        patch.admit_pressure_perimeter(fixed, &mut cancelled)?;
        let gas = gas.into();
        let [lo, hi] = volume_interval_m3;
        if !exterior_pressure_pa.is_finite()
            || exterior_pressure_pa < 0.0
            || !lo.is_finite()
            || !hi.is_finite()
            || lo <= 0.0
            || lo >= hi
            || gas.volume_m3() < lo
            || gas.volume_m3() > hi
        {
            return Err(GasMembraneError::Invalid(
                "explicit nonnegative exterior pressure and positive volume interval containing the reference gas required",
            ));
        }
        let free_dofs: Vec<_> = (0..fixed.len()).filter(|&i| fixed[i].is_none()).collect();
        if free_dofs.is_empty() || free_dofs.len() > MAX_VOLUME_COORDINATES {
            return Err(GasMembraneError::Invalid("sealed membrane needs 1..=256 free coordinates"));
        }
        let mut model = Self {
            patch,
            gas,
            reference_m: reference_m.to_vec(),
            reference_strain_energy_j: 0.0,
            fixed: fixed.to_vec(),
            free_dofs,
            exterior_pressure_pa,
            volume_interval_m3,
        };
        model.admit_shape(reference_m)?;
        model.reference_strain_energy_j = patch.evaluate(reference_m, &mut cancelled)?.strain_energy_j;
        model.gas.pressure_response()?;
        poll(&mut cancelled)?;
        Ok(model)
    }

    /// Original physical membrane, including its conserved reference mass.
    #[must_use]
    pub const fn membrane(&self) -> &MembranePatch {
        self.patch
    }
    /// Unchanged source inventory defining the reference entropy and volume.
    #[must_use]
    pub const fn reference_gas(&self) -> &MembraneGas {
        &self.gas
    }
    /// Canonical flat node/axis indices for residual and tangent coordinates.
    #[must_use]
    pub fn free_dofs(&self) -> &[usize] {
        &self.free_dofs
    }

    fn admit_shape(&self, displacement: &[[f64; 3]]) -> Result<(), GasMembraneError> {
        if !finite_vectors(displacement, self.patch.mesh().nodes.len())
            || self.fixed.iter().enumerate().any(|(dof, prescribed)| {
                prescribed.is_some_and(|value| displacement[dof / 3][dof % 3] != value)
            })
        {
            return Err(GasMembraneError::Invalid("finite shape with unchanged prescribed coordinates required"));
        }
        Ok(())
    }

    fn expand(&self, coordinates: &[f64], direction: bool) -> Vec<[f64; 3]> {
        let mut u = if direction {
            vec![[0.0; 3]; self.reference_m.len()]
        } else {
            self.reference_m.clone()
        };
        for (&dof, &value) in self.free_dofs.iter().zip(coordinates) {
            u[dof / 3][dof % 3] = value;
        }
        u
    }

    /// Evaluate the original isentrope at this geometry without mutating any
    /// accepted source, membrane or reference state. Loads are held dead forces.
    pub fn evaluate(
        &self,
        displacement_m: &[[f64; 3]],
        dead_load_n: &[[f64; 3]],
        gas: GasEntropySolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<GasMembraneEvaluation, GasMembraneError> {
        Ok(self.trial(displacement_m, dead_load_n, gas, &mut cancelled)?.value)
    }

    fn trial(
        &self,
        displacement: &[[f64; 3]],
        loads: &[[f64; 3]],
        controls: GasEntropySolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Trial, GasMembraneError> {
        poll(cancelled)?;
        self.admit_shape(displacement)?;
        if !finite_vectors(loads, displacement.len()) {
            return Err(GasMembraneError::Invalid("one finite dead force per membrane node required"));
        }
        let membrane = self.patch.evaluate(displacement, &mut *cancelled)?;
        let swept = self.patch.pressure_work_j(&self.reference_m, displacement, 1.0, &mut *cancelled)?;
        let volume = self.gas.volume_m3() + swept;
        if !volume.is_finite()
            || volume < self.volume_interval_m3[0]
            || volume > self.volume_interval_m3[1]
            || (swept != 0.0 && volume == self.gas.volume_m3())
        {
            return Err(GasMembraneError::Invalid("swept cavity volume is outside its interval or unrepresentable"));
        }
        let (gas, gas_volume) = self.gas.at_volume(volume, controls, &mut *cancelled)?;
        let pressure = gas.pressure_response()?;
        let gradient_all = self.patch.pressure_load_n(displacement, 1.0, &mut *cancelled)?;
        let gradient: Vec<_> = self.free_dofs.iter().map(|&i| gradient_all[i / 3][i % 3]).collect();
        let forces = gas.forces(&gradient, self.exterior_pressure_pa)?;
        let net_pressure = pressure.pressure_pa - self.exterior_pressure_pa;
        let mut residual_n = membrane.internal_gradient_n;
        let (mut dead_work, mut correction) = (0.0, 0.0);
        for node in 0..displacement.len() {
            if node % 64 == 0 { poll(cancelled)?; }
            for axis in 0..3 {
                residual_n[node][axis] -= loads[node][axis];
                if self.fixed[3 * node + axis].is_some() {
                    residual_n[node][axis] -= net_pressure * gradient_all[node][axis];
                }
                let work = loads[node][axis] * (displacement[node][axis] - self.reference_m[node][axis]);
                let adjusted = work - correction;
                let next = dead_work + adjusted;
                correction = (next - dead_work) - adjusted;
                dead_work = next;
            }
        }
        for (&dof, force) in self.free_dofs.iter().zip(forces) {
            residual_n[dof / 3][dof % 3] -= force;
        }
        let exterior_work = self.exterior_pressure_pa * swept;
        let potential = (membrane.strain_energy_j - self.reference_strain_energy_j)
            + gas_volume.stored_energy_change_j + exterior_work - dead_work;
        if !finite_vectors(&residual_n, displacement.len())
            || ![dead_work, exterior_work, potential].iter().all(|x| x.is_finite())
        {
            return Err(GasMembraneError::Invalid("coupled force or conservative work is unrepresentable"));
        }
        poll(cancelled)?;
        Ok(Trial {
            value: GasMembraneEvaluation {
                displacement_m: displacement.to_vec(), gas, gas_volume, pressure, residual_n,
                swept_volume_m3: swept, strain_energy_j: membrane.strain_energy_j,
                dead_load_work_j: dead_work, exterior_pressure_work_j: exterior_work,
                potential_energy_change_j: potential,
            },
            gradient,
        })
    }

    /// Exact coupled force-residual tangent in `free_dofs()` order. Direction
    /// entries are Cartesian displacement increments [m]; fixed entries stay
    /// zero. Both follower geometry and the original gas bulk modulus contribute.
    /// This is an equilibrium tangent, not a smooth-mode or stability claim.
    pub fn tangent_apply(
        &self,
        displacement_m: &[[f64; 3]],
        direction: &[f64],
        gas: GasEntropySolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<f64>, GasMembraneError> {
        let zero = vec![[0.0; 3]; self.reference_m.len()];
        let trial = self.trial(displacement_m, &zero, gas, &mut cancelled)?;
        self.tangent_at(&trial, direction, &mut cancelled)
    }

    fn tangent_at(
        &self,
        trial: &Trial,
        direction: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<f64>, GasMembraneError> {
        poll(cancelled)?;
        if direction.len() != self.free_dofs.len() || direction.iter().any(|v| !v.is_finite()) {
            return Err(GasMembraneError::Invalid("one finite tangent direction per free coordinate required"));
        }
        let full_direction = self.expand(direction, true);
        let displacement = &trial.value.displacement_m;
        let internal = self.patch.tangent_apply(displacement, &full_direction, &mut *cancelled)?;
        let hessian = self.patch.pressure_tangent_apply(displacement, &full_direction, 1.0, &mut *cancelled)?;
        let hessian: Vec<_> = self.free_dofs.iter().map(|&i| hessian[i / 3][i % 3]).collect();
        let pressure = trial.value.gas.tangent(
            &trial.gradient, self.exterior_pressure_pa, direction, &hessian,
        )?;
        let action: Vec<_> = self.free_dofs.iter().zip(pressure)
            .map(|(&i, load)| internal[i / 3][i % 3] - load).collect();
        if action.iter().any(|v| !v.is_finite()) {
            return Err(GasMembraneError::Invalid("coupled tangent overflow"));
        }
        poll(cancelled)?;
        Ok(action)
    }
}

#[cfg(test)]
mod tests;
