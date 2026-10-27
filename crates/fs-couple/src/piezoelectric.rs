//! Longitudinal material-driven sensor/actuator on the shared physical stepper.
//!
//! One rod has a fixed root, free transverse stress, uniform end electrodes and
//! the axial displacement shape u(z)=x*z/L. Integrating the original material
//! energy gives H=A*L*w(x/L,Q/A)+p^2/(2m), m=rho*A*L/3. Its derivatives give
//! force, voltage and velocity. The SAME fs-phs descriptor step advances x, p
//! and Q; there is no separate actuator gain, sensor gain or mechanical solver.
//! The low-frequency single-element reduction omits resolved axial waves,
//! transverse deformation, fringing, dielectric loss and depoling.

use fs_material::piezoelectric::{
    LongitudinalPiezoelectric, PiezoelectricError, PiezoelectricResponse,
};
use fs_phs::{
    DescriptorPortHamiltonian, PhsError, PortHamiltonian, StepRecord, Storage, step_descriptor,
};

/// Explicit caller-admitted small-signal operating envelope. These bounds are
/// model restrictions, not invented measured strength or breakdown data.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SmallSignalLimits {
    /// Maximum absolute axial engineering strain, strictly less than one.
    pub max_abs_strain: f64,
    /// Maximum absolute axial electric field [V/m].
    pub max_abs_electric_field_v_m: f64,
}

/// Material, numerical or operating-domain refusal before state publication.
#[derive(Debug)]
pub enum PiezoelectricCouplingError {
    /// Original material source or constitutive evaluation failed.
    Material(PiezoelectricError),
    /// Existing port-Hamiltonian admission or implicit solve failed.
    Native(PhsError),
    /// Invalid geometry, drive, amplitude or unrepresentable arithmetic.
    Invalid(&'static str),
    /// The independent physical work account exceeds the explicit tolerance.
    Balance {
        /// Actual storage change minus mechanical and electrical supply [J].
        residual_j: f64,
        /// Requested absolute work tolerance [J].
        tolerance_j: f64,
    },
    /// Cancellation observed before publishing the joint state.
    Cancelled,
}
impl From<PiezoelectricError> for PiezoelectricCouplingError {
    fn from(error: PiezoelectricError) -> Self {
        Self::Material(error)
    }
}
impl From<PhsError> for PiezoelectricCouplingError {
    fn from(error: PhsError) -> Self {
        Self::Native(error)
    }
}
impl core::fmt::Display for PiezoelectricCouplingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "piezoelectric coupling: {self:?}")
    }
}
impl std::error::Error for PiezoelectricCouplingError {}
type Result<T> = core::result::Result<T, PiezoelectricCouplingError>;
fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn invalid(what: &'static str) -> PiezoelectricCouplingError {
    PiezoelectricCouplingError::Invalid(what)
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(PiezoelectricCouplingError::Cancelled)
    } else {
        Ok(())
    }
}

/// Source-bound geometry with material axis 3 along the rod. Positive extension
/// is tensile; positive terminal voltage is the drop along +3 and positive Q is
/// charge on the electrode at its beginning. Reversing poling changes signed
/// d33 in the material input; a world rotation does not change these local scalars.
#[derive(Clone, Debug, PartialEq)]
pub struct PiezoelectricRod {
    material: LongitudinalPiezoelectric,
    length_m: f64,
    area_m2: f64,
    volume_m3: f64,
    modal_mass_kg: f64,
    limits: SmallSignalLimits,
}
impl PiezoelectricRod {
    /// Derive stiffness, electric storage, coupling and consistent inertia from
    /// the SAME material and reference geometry. All thermodynamic context stays
    /// fixed; live coefficient changes need an explicit thermal/parameter-work owner.
    pub fn new(
        material: LongitudinalPiezoelectric,
        length_m: f64,
        area_m2: f64,
        limits: SmallSignalLimits,
    ) -> Result<Self> {
        let volume_m3 = length_m * area_m2;
        let modal_mass_kg = material.density_kg_m3() * volume_m3 / 3.0;
        if [
            length_m,
            area_m2,
            volume_m3,
            modal_mass_kg,
            1.0 / modal_mass_kg,
            limits.max_abs_strain,
            limits.max_abs_electric_field_v_m,
        ]
        .iter()
        .any(|&x| !positive(x))
            || limits.max_abs_strain >= 1.0
        {
            return Err(invalid(
                "positive finite geometry, inertia and small-signal limits required",
            ));
        }
        let rod = Self {
            material,
            length_m,
            area_m2,
            volume_m3,
            modal_mass_kg,
            limits,
        };
        if [
            rod.short_circuit_stiffness_n_m(),
            rod.open_circuit_stiffness_n_m(),
            rod.axially_clamped_capacitance_f(),
            1.0 / rod.axially_clamped_capacitance_f(),
        ]
        .iter()
        .any(|&x| !positive(x))
            || !rod.coupling_c_m().is_finite()
            || rod.material.strain_charge_coefficient_m_v() != 0.0 && rod.coupling_c_m() == 0.0
        {
            return Err(invalid(
                "derived rod stiffness, capacitance or coupling is unrepresentable",
            ));
        }
        Ok(rod)
    }
    /// Original immutable material snapshot and source receipts.
    pub const fn material(&self) -> &LongitudinalPiezoelectric {
        &self.material
    }
    /// Reference rod length [m].
    pub const fn length_m(&self) -> f64 {
        self.length_m
    }
    /// Reference electrode/cross-sectional area [m^2].
    pub const fn area_m2(&self) -> f64 {
        self.area_m2
    }
    /// Consistent mass for u(z)=x*z/L [kg], not the entire rod as a tip lump.
    pub const fn modal_mass_kg(&self) -> f64 {
        self.modal_mass_kg
    }
    /// Fixed-field axial stiffness [N/m].
    pub fn short_circuit_stiffness_n_m(&self) -> f64 {
        self.material.short_circuit_modulus_pa() * self.area_m2 / self.length_m
    }
    /// Fixed-charge axial stiffness [N/m].
    pub fn open_circuit_stiffness_n_m(&self) -> f64 {
        self.material.open_circuit_modulus_pa() * self.area_m2 / self.length_m
    }
    /// Terminal capacitance with axial extension fixed [F].
    pub fn axially_clamped_capacitance_f(&self) -> f64 {
        self.material.axially_clamped_permittivity_f_m() * self.area_m2 / self.length_m
    }
    /// Reciprocal theta in Q=C*V+theta*x and force=k*x-theta*V [C/m = N/V].
    pub fn coupling_c_m(&self) -> f64 {
        self.material.stress_coefficient_c_m2() * self.area_m2 / self.length_m
    }
    /// Retained caller-admitted amplitude bounds.
    pub const fn limits(&self) -> SmallSignalLimits {
        self.limits
    }

    pub(crate) fn response(&self, extension: f64, charge: f64) -> Result<PiezoelectricResponse> {
        let strain = extension / self.length_m;
        let response = self.material.response(strain, charge / self.area_m2)?;
        if strain.abs() > self.limits.max_abs_strain
            || response.electric_field_v_m.abs() > self.limits.max_abs_electric_field_v_m
        {
            return Err(invalid(
                "rod left its declared small-signal strain or electric-field envelope",
            ));
        }
        Ok(response)
    }
    pub(crate) fn energy(&self, x: &[f64]) -> Result<f64> {
        let response = self.response(x[0], x[2])?;
        let kinetic = 0.5 * x[1] * (x[1] / self.modal_mass_kg);
        let material_energy = self.volume_m3 * response.stored_energy_j_m3;
        let energy = material_energy + kinetic;
        if !energy.is_finite()
            || x[1] != 0.0 && kinetic == 0.0
            || response.stored_energy_j_m3 != 0.0 && material_energy == 0.0
        {
            return Err(invalid("rod energy is not representable"));
        }
        Ok(energy)
    }

    /// Reusable two-port component for the existing fs-phs interconnections.
    /// State is [extension m, momentum kg m/s, charge C]; inputs are [force N,
    /// current A] and outputs [velocity m/s, voltage V]. Original material and
    /// amplitude checks remain in its energy law. Coupled callers own acceptance
    /// of the complete composed step, including each component's source domain.
    pub fn port_hamiltonian(&self) -> Result<PortHamiltonian> {
        Ok(PortHamiltonian::new(
            3,
            2,
            vec![0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            vec![0.0; 9],
            vec![0.0, 0.0, 1.0, 0.0, 0.0, 1.0],
            Box::new(RodStorage(self.clone())),
        )?)
    }

    fn descriptor(&self, voltage_drive: bool) -> Result<DescriptorPortHamiltonian> {
        let n = if voltage_drive { 4 } else { 3 };
        let mut j = vec![0.0; n * n];
        let mut g = vec![0.0; n * 2];
        j[1] = 1.0;
        j[n] = -1.0;
        g[2] = 1.0; // external force into momentum
        if voltage_drive {
            // qdot=i, 0=-dH/dq+V. i is an algebraic multiplier, not storage.
            j[2 * n + 3] = 1.0;
            j[3 * n + 2] = -1.0;
            g[7] = 1.0;
        } else {
            g[5] = 1.0;
        }
        Ok(DescriptorPortHamiltonian::new(
            n,
            3,
            2,
            j,
            vec![0.0; n * n],
            g,
            Box::new(RodStorage(self.clone())),
        )?)
    }
}

#[derive(Clone)]
struct RodStorage(PiezoelectricRod);
impl Storage for RodStorage {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        self.0.energy(x).unwrap_or(f64::NAN)
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        match self.0.response(x[0], x[2]) {
            Ok(response) => {
                out[0] = self.0.area_m2 * response.stress_pa;
                out[1] = x[1] / self.0.modal_mass_kg;
                out[2] = self.0.length_m * response.electric_field_v_m;
            }
            Err(_) => out.fill(f64::NAN),
        }
    }
    fn midpoint_energy_remainder(&self, a: &[f64], b: &[f64], mid: &[f64]) -> Option<f64> {
        Some(if [a, b, mid].iter().all(|x| self.0.energy(x).is_ok()) {
            0.0
        } else {
            f64::NAN
        })
    }
}

/// Held electrical boundary during one mechanical/electrical transaction.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ElectricalDrive {
    /// Current INTO the initial electrode [A]; zero is an open circuit.
    Current(f64),
    /// Ideal voltage along axis 3 [V]; zero is a short circuit. Initial charge
    /// must already satisfy this voltage. Switching impulses are not invented.
    Voltage(f64),
}

/// Explicit physical work and ideal-clamp consistency tolerances.
#[derive(Clone, Copy, Debug)]
pub struct PiezoelectricSolve {
    /// Absolute allowed storage/work residual [J], not a time-accuracy guarantee.
    pub energy_tolerance_j: f64,
    /// Allowed initial/final voltage-constraint residual [V].
    pub voltage_tolerance_v: f64,
}

/// Accepted mechanical and electric state with its immutable material/geometry.
#[derive(Clone, Debug, PartialEq)]
pub struct PiezoelectricBody {
    rod: PiezoelectricRod,
    state: [f64; 3],
}

/// A joint proposed state and the native power-conjugate step that produced it.
pub struct PiezoelectricCandidate {
    body: PiezoelectricBody,
    step: PiezoelectricStep,
}
impl PiezoelectricCandidate {
    /// Complete proposed mechanical/electrical endpoint.
    pub const fn body(&self) -> &PiezoelectricBody {
        &self.body
    }
    /// Native solver and independently reconstructed physical work.
    pub const fn step(&self) -> &PiezoelectricStep {
        &self.step
    }
    /// Publish the joint endpoint and its account together.
    pub fn into_parts(self) -> (PiezoelectricBody, PiezoelectricStep) {
        (self.body, self.step)
    }
}

/// No heat is generated by this reversible small-signal material rung.
pub struct PiezoelectricStep {
    /// Native descriptor result; the voltage-driven case appends algebraic current.
    pub native: StepRecord,
    /// Applied force times actual extension change [J].
    pub mechanical_supply_j: f64,
    /// Current/voltage source work [J]; zero for an open or short circuit.
    pub electrical_supply_j: f64,
    /// Actual stored-energy change minus both source works [J].
    pub energy_residual_j: f64,
}

impl PiezoelectricBody {
    /// Admit a joint state. Positive initial charge, stress and velocity use the
    /// rod's oriented axis/electrode convention. No charge is projected to match
    /// a later ideal voltage source, and no pre-existing energy is discarded.
    pub fn new(
        rod: PiezoelectricRod,
        extension_m: f64,
        velocity_m_s: f64,
        charge_c: f64,
    ) -> Result<Self> {
        let momentum = rod.modal_mass_kg * velocity_m_s;
        let state = [extension_m, momentum, charge_c];
        if !velocity_m_s.is_finite() || velocity_m_s != 0.0 && momentum == 0.0 {
            return Err(invalid("initial velocity or momentum is unrepresentable"));
        }
        rod.energy(&state)?;
        Ok(Self { rod, state })
    }
    /// Original geometry and material binding.
    pub const fn rod(&self) -> &PiezoelectricRod {
        &self.rod
    }
    /// Native [extension, momentum, charge] coordinates in SI units.
    pub const fn state(&self) -> [f64; 3] {
        self.state
    }
    /// Axial tip extension [m].
    pub const fn extension_m(&self) -> f64 {
        self.state[0]
    }
    /// Axial tip velocity [m/s].
    pub fn velocity_m_s(&self) -> f64 {
        self.state[1] / self.rod.modal_mass_kg
    }
    /// Free charge on the electrode at the beginning of axis 3 [C].
    pub const fn charge_c(&self) -> f64 {
        self.state[2]
    }
    /// Accepted source-law stress, electric field and energy density.
    pub fn response(&self) -> Result<PiezoelectricResponse> {
        self.rod.response(self.state[0], self.state[2])
    }
    /// Accepted terminal voltage [V].
    pub fn voltage_v(&self) -> Result<f64> {
        Ok(self.rod.length_m * self.response()?.electric_field_v_m)
    }
    /// Combined material and kinetic energy [J].
    pub fn stored_energy_j(&self) -> Result<f64> {
        self.rod.energy(&self.state)
    }

    /// Evolve both domains through the native descriptor integrator. Cancellation,
    /// unsupported amplitude, incompatible ideal clamp, solver refusal and failed
    /// work balance leave this body unchanged. A compatible boundary may change
    /// between steps without resetting its charge or mechanical history.
    pub fn step_candidate(
        &self,
        force_n: f64,
        electrical: ElectricalDrive,
        duration_s: f64,
        controls: PiezoelectricSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PiezoelectricCandidate> {
        poll(&mut cancelled)?;
        let (drive, voltage_drive) = match electrical {
            ElectricalDrive::Current(value) => (value, false),
            ElectricalDrive::Voltage(value) => (value, true),
        };
        if !force_n.is_finite()
            || !drive.is_finite()
            || !positive(duration_s)
            || !positive(controls.energy_tolerance_j)
            || !positive(controls.voltage_tolerance_v)
        {
            return Err(invalid(
                "finite held drives and positive duration/work/voltage tolerances required",
            ));
        }
        if voltage_drive && (self.voltage_v()? - drive).abs() > controls.voltage_tolerance_v {
            return Err(invalid(
                "ideal voltage clamp is incompatible with initial charge; switching needs an explicit circuit impulse model",
            ));
        }
        let system = self.rod.descriptor(voltage_drive)?;
        let mut x0 = self.state.to_vec();
        if voltage_drive {
            x0.push(self.rod.coupling_c_m() * self.velocity_m_s());
        }
        let inputs = [force_n, drive];
        if voltage_drive {
            // Reuse native algebraic admission, while forbidding even a numerical
            // projection of conserved coordinates during initial consistency.
            let admitted = step_descriptor(&system, &x0, &inputs, 0.0)?;
            if admitted.x[..3] != self.state {
                return Err(invalid(
                    "ideal clamp admission changed a conserved coordinate",
                ));
            }
        }
        poll(&mut cancelled)?;
        let native = step_descriptor(&system, &x0, &inputs, duration_s)?;
        let state = [native.x[0], native.x[1], native.x[2]];
        let body = Self {
            rod: self.rod.clone(),
            state,
        };
        let before = self.stored_energy_j()?;
        let after = body.stored_energy_j()?;
        let mechanical_supply_j = force_n * (state[0] - self.state[0]);
        let electrical_supply_j = if voltage_drive {
            drive * (state[2] - self.state[2])
        } else {
            duration_s * drive * 0.5 * (self.voltage_v()? + body.voltage_v()?)
        };
        let energy_residual_j = (after - before) - mechanical_supply_j - electrical_supply_j;
        if voltage_drive && (body.voltage_v()? - drive).abs() > controls.voltage_tolerance_v {
            return Err(invalid(
                "coupled endpoint violated the held voltage constraint",
            ));
        }
        if [
            before,
            after,
            mechanical_supply_j,
            electrical_supply_j,
            energy_residual_j,
            native.supplied,
            native.delta_h,
            native.dissipated,
        ]
        .iter()
        .any(|x| !x.is_finite())
        {
            return Err(invalid("coupled piezoelectric work is not representable"));
        }
        let tolerance = controls.energy_tolerance_j;
        if energy_residual_j.abs() > tolerance
            || native.balance_residual().abs() > tolerance
            || ((after - before) - native.delta_h).abs() > tolerance
            || (mechanical_supply_j + electrical_supply_j - native.supplied).abs() > tolerance
            || native.dissipated != 0.0
        {
            return Err(PiezoelectricCouplingError::Balance {
                residual_j: energy_residual_j,
                tolerance_j: tolerance,
            });
        }
        poll(&mut cancelled)?;
        Ok(PiezoelectricCandidate {
            body,
            step: PiezoelectricStep {
                native,
                mechanical_supply_j,
                electrical_supply_j,
                energy_residual_j,
            },
        })
    }
}
