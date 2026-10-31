//! Shared isothermal aqueous thermodynamics and selective membranes.
//!
//! Concentrations count independent dissolved particles in mol/m³: a millimolar
//! concentration is 1 mol/m³, and an already enumerated ion is not dissociated a
//! second time. The solvent is incompressible with an explicitly supplied
//! partial molar volume. [`aqueous::AqueousSolution`] keeps finite inventories
//! with either the ideal law or the consistent free energy in [`nonideal`].
//! [`donnan`] separately supplies finite-bath, ideal ionic equilibrium; it must
//! not be combined with a nonideal pressure correction alone. Solute partial
//! volumes, ion reactions and temperature differences are not modeled. Optional
//! steady concentration polarization is provided by [`polarization`].
//! Coefficients need their own measured
//! applicability; raw numerical constructors do not confer material authority.
//!
//! Positive flux is left to right. The ideal selective membrane law is
//! `Jv = Lp (p_left - p_right - RT sum(sigma_i (c_left_i - c_right_i)))`.
//! Solute transport integrates the constant-coefficient Spiegler–Kedem equation
//! through the membrane, retaining both diffusion and solvent drag at arbitrary
//! membrane Peclet number. The exponential fit has nonnegative donor stream
//! coefficients and finite pure-solvent endpoints; no concentration floor is
//! introduced. It is a constitutive evaluation, not a time integrator.
//!
//! References: IUPAC, <https://goldbook.iupac.org/terms/view/O04344> (osmotic
//! pressure); Spiegler and Kedem, Desalination 1 (1966), 311–326,
//! <https://doi.org/10.1016/S0011-9164(00)80018-1> (membrane transport).

use crate::fluid_transport::FluidTransport;
use fs_math::det;
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

/// Exact steady diffusive films coupled to the original selective membrane.
pub mod polarization;

/// Thermodynamically consistent finite-domain solution interactions.
pub mod nonideal;
/// Shared finite aqueous inventories with their selected thermodynamics.
pub mod aqueous;
pub use aqueous::AqueousSolution;
/// Finite-bath ideal ionic partition with fixed charge and explicit trapped particles.
pub mod donnan;

/// Bounded independent-solute layout; an empty layout represents pure solvent.
pub const MAX_OSMOTIC_SOLUTES: usize = 64;

/// Invalid dimensions/layout, nonphysical inputs, or unrepresentable arithmetic.
#[derive(Clone, Debug, PartialEq)]
pub enum OsmosisError {
    /// A finite physical-domain or model requirement was not met.
    Invalid(&'static str),
    /// The concentration, inventory, or output layout differs from the owner.
    Shape {
        /// Number of independent dissolved species required.
        expected: usize,
        /// Number supplied by the caller.
        actual: usize,
    },
    /// A result cannot be represented by finite coherent-SI values.
    Arithmetic(&'static str),
    /// The original fluid owner refused its transport snapshot.
    Fluid(crate::MaterialError),
    /// Nonideal thermodynamics refused its declared concentration support.
    Nonideal(nonideal::NonidealSolutionError),
}

impl core::fmt::Display for OsmosisError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "osmotic material: {self:?}")
    }
}
impl std::error::Error for OsmosisError {}

impl From<nonideal::NonidealSolutionError> for OsmosisError {
    fn from(value: nonideal::NonidealSolutionError) -> Self {
        Self::Nonideal(value)
    }
}

type Result<T> = core::result::Result<T, OsmosisError>;

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn nonnegative(value: f64) -> bool {
    value.is_finite() && value >= 0.0
}
fn finite(value: f64, what: &'static str) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(OsmosisError::Arithmetic(what))
    }
}
fn shape(expected: usize, actual: usize) -> Result<()> {
    if expected == actual {
        Ok(())
    } else {
        Err(OsmosisError::Shape { expected, actual })
    }
}
fn layout(count: usize) -> Result<()> {
    if count <= MAX_OSMOTIC_SOLUTES {
        Ok(())
    } else {
        Err(OsmosisError::Invalid("independent-solute count exceeds 64"))
    }
}
fn thermal_scale(temperature_k: f64) -> Result<f64> {
    if !positive(temperature_k) {
        return Err(OsmosisError::Invalid(
            "temperature must be finite, positive kelvin",
        ));
    }
    finite(R * temperature_k, "RT is not representable")
}

// Fixed-order compensated accumulation avoids losing dilute species against
// dominant solutes, or small driving differences between counter-gradients.
#[derive(Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) {
        let next = self.value + value;
        self.correction += if self.value.abs() >= value.abs() {
            (self.value - next) + value
        } else {
            (value - next) + self.value
        };
        self.value = next;
    }
    fn finish(self) -> Result<f64> {
        finite(
            self.value + self.correction,
            "solute sum is not representable",
        )
    }
}

/// van't Hoff pressure `RT sum(c_i)` [Pa] for independent solute particles.
///
/// `c_i` is mol/m³ of solution in the declared dilute, negligible-solute-volume
/// approximation. No chemical-name, dissociation, or osmotic-coefficient
/// correction is guessed. Empty composition has exactly zero osmotic pressure.
pub fn ideal_dilute_osmotic_pressure_pa(
    temperature_k: f64,
    particle_concentrations_mol_m3: &[f64],
) -> Result<f64> {
    layout(particle_concentrations_mol_m3.len())?;
    let rt = thermal_scale(temperature_k)?;
    let mut total = Sum::default();
    for &concentration in particle_concentrations_mol_m3 {
        if !nonnegative(concentration) {
            return Err(OsmosisError::Invalid(
                "concentrations must be finite, nonnegative mol/m³",
            ));
        }
        total.add(concentration);
    }
    finite(
        rt * total.finish()?,
        "osmotic pressure is not representable",
    )
}

/// IUPAC solvent-activity pressure `-RT ln(a_w)/Vbar_w` [Pa].
///
/// The supplied activity may come from a nonideal solution model or measured
/// data. The solvent partial molar volume is treated as pressure-independent;
/// this helper supplies neither that data nor a species diffusion law. A total
/// activity must not be split among species-specific reflection coefficients.
pub fn osmotic_pressure_from_water_activity_pa(
    temperature_k: f64,
    water_activity: f64,
    solvent_molar_volume_m3_mol: f64,
) -> Result<f64> {
    let rt = thermal_scale(temperature_k)?;
    if !positive(water_activity) || water_activity > 1.0 || !positive(solvent_molar_volume_m3_mol) {
        return Err(OsmosisError::Invalid(
            "activity must be in (0, 1] and solvent molar volume positive",
        ));
    }
    finite(
        rt * (-det::ln(water_activity) / solvent_molar_volume_m3_mol),
        "activity-derived osmotic pressure is not representable",
    )
}

/// Solvent activity from a nonnegative osmotic pressure at fixed molar volume.
/// A result below the positive floating-point range refuses rather than becoming
/// a zero-activity state with an infinite chemical potential.
pub fn water_activity_from_osmotic_pressure(
    temperature_k: f64,
    osmotic_pressure_pa: f64,
    solvent_molar_volume_m3_mol: f64,
) -> Result<f64> {
    let rt = thermal_scale(temperature_k)?;
    if !nonnegative(osmotic_pressure_pa) || !positive(solvent_molar_volume_m3_mol) {
        return Err(OsmosisError::Invalid(
            "osmotic pressure must be nonnegative and solvent molar volume positive",
        ));
    }
    let exponent = finite(
        -(osmotic_pressure_pa / rt) * solvent_molar_volume_m3_mol,
        "activity exponent is not representable",
    )?;
    let activity = det::exp(exponent);
    if positive(activity) && activity <= 1.0 {
        Ok(activity)
    } else {
        Err(OsmosisError::Arithmetic(
            "positive solvent activity is not representable",
        ))
    }
}

/// Finite dilute solution; solute amounts stay fixed when solvent volume changes.
///
/// The moles count independent dissolved particles and retain their caller-owned
/// order. All consumers must use that same order. Temperature is externally held
/// fixed: this state carries no caloric, mechanical-storage, or phase model.
#[derive(Clone, Debug, PartialEq)]
pub struct IdealDiluteSolution {
    temperature_k: f64,
    solvent_molar_volume_m3_mol: f64,
    solvent_volume_m3: f64,
    solute_moles: Box<[f64]>,
}

impl IdealDiluteSolution {
    /// Prepare a finite solution with explicit kelvin, m³/mol, m³ and mol units.
    pub fn new(
        temperature_k: f64,
        solvent_molar_volume_m3_mol: f64,
        solvent_volume_m3: f64,
        solute_moles: Vec<f64>,
    ) -> Result<Self> {
        thermal_scale(temperature_k)?;
        if !positive(solvent_molar_volume_m3_mol) {
            return Err(OsmosisError::Invalid(
                "solvent molar volume must be finite and positive",
            ));
        }
        layout(solute_moles.len())?;
        validate_inventory(temperature_k, solvent_volume_m3, &solute_moles)?;
        Ok(Self {
            temperature_k,
            solvent_molar_volume_m3_mol,
            solvent_volume_m3,
            solute_moles: solute_moles.into_boxed_slice(),
        })
    }

    /// Externally held common solution temperature [K].
    pub const fn temperature_k(&self) -> f64 {
        self.temperature_k
    }
    /// Explicit incompressible solvent partial molar volume [m³/mol].
    pub const fn solvent_molar_volume_m3_mol(&self) -> f64 {
        self.solvent_molar_volume_m3_mol
    }
    /// Finite solvent volume [m³]; dilute solute volume is neglected.
    pub const fn solvent_volume_m3(&self) -> f64 {
        self.solvent_volume_m3
    }
    /// Conserved independent-solute particle amounts [mol], in the original order.
    pub fn solute_moles(&self) -> &[f64] {
        &self.solute_moles
    }

    /// Recompute molar concentrations from current finite inventories, without allocation.
    pub fn concentrations_into(&self, output: &mut [f64]) -> Result<()> {
        shape(self.solute_moles.len(), output.len())?;
        for (concentration, &amount) in output.iter_mut().zip(&self.solute_moles) {
            *concentration = amount / self.solvent_volume_m3;
        }
        Ok(())
    }

    /// Current dilute osmotic pressure [Pa], recomputed after every water exchange.
    pub fn osmotic_pressure_pa(&self) -> Result<f64> {
        let concentration = validate_inventory(
            self.temperature_k,
            self.solvent_volume_m3,
            &self.solute_moles,
        )?;
        finite(
            thermal_scale(self.temperature_k)? * concentration,
            "osmotic pressure is not representable",
        )
    }

    /// Activity consistent with THIS dilute model: `ln(a_w)=-Vbar_w sum(c_i)`.
    /// This is not the distinct finite-mole-fraction ideal-mixture model.
    pub fn water_activity(&self) -> Result<f64> {
        water_activity_from_osmotic_pressure(
            self.temperature_k,
            self.osmotic_pressure_pa()?,
            self.solvent_molar_volume_m3_mol,
        )
    }

    /// New solution after a pure-solvent exchange; all original solute amounts remain.
    /// The returned owned state allocates once to preserve those inventories.
    pub fn with_solvent_volume_m3(&self, solvent_volume_m3: f64) -> Result<Self> {
        validate_inventory(self.temperature_k, solvent_volume_m3, &self.solute_moles)?;
        Ok(Self {
            solvent_volume_m3,
            ..self.clone()
        })
    }

    /// Atomically replace finite inventories in the SAME species layout.
    /// The complete candidate is checked before any field changes; no allocation.
    pub fn replace_inventories(
        &mut self,
        solvent_volume_m3: f64,
        solute_moles: &[f64],
    ) -> Result<()> {
        shape(self.solute_moles.len(), solute_moles.len())?;
        validate_inventory(self.temperature_k, solvent_volume_m3, solute_moles)?;
        self.solvent_volume_m3 = solvent_volume_m3;
        self.solute_moles.copy_from_slice(solute_moles);
        Ok(())
    }
}

fn validate_inventory(temperature_k: f64, volume_m3: f64, amounts: &[f64]) -> Result<f64> {
    if !positive(volume_m3) {
        return Err(OsmosisError::Invalid(
            "finite solvent volume must stay strictly positive",
        ));
    }
    let mut total = Sum::default();
    for &amount in amounts {
        if !nonnegative(amount) {
            return Err(OsmosisError::Invalid(
                "solute amounts must be finite and nonnegative",
            ));
        }
        total.add(finite(
            amount / volume_m3,
            "molar concentration is not representable",
        )?);
    }
    let concentration = total.finish()?;
    finite(
        thermal_scale(temperature_k)? * concentration,
        "osmotic pressure is not representable",
    )?;
    Ok(concentration)
}

/// Thickness-normalized hydraulic membrane permeance [m/(Pa s)].
/// Distinct from intrinsic porous permeability [m²] and head conductivity [m/s].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HydraulicPermeance(f64);

impl HydraulicPermeance {
    /// Admit measured/model-explicit `Lp`; zero is a solvent-impermeable boundary.
    pub fn new(value_m_per_pa_s: f64) -> Result<Self> {
        if nonnegative(value_m_per_pa_s) {
            Ok(Self(value_m_per_pa_s))
        } else {
            Err(OsmosisError::Invalid(
                "hydraulic permeance must be finite and nonnegative m/(Pa s)",
            ))
        }
    }
    /// Numerical SI value [m/(Pa s)].
    pub const fn value(self) -> f64 {
        self.0
    }

    /// Explicit Darcy-slab crosswalk `Lp=k/(mu thickness)` at the original
    /// fluid's current material evaluation. It supplies no solute selectivity,
    /// geometry, temperature, or material provenance beyond those caller inputs.
    pub fn from_darcy_slab(
        intrinsic_permeability_m2: f64,
        thickness_m: f64,
        fluid: &(impl FluidTransport + ?Sized),
    ) -> Result<Self> {
        if !nonnegative(intrinsic_permeability_m2) || !positive(thickness_m) {
            return Err(OsmosisError::Invalid(
                "Darcy crosswalk needs nonnegative intrinsic m² and positive thickness m",
            ));
        }
        let viscosity = fluid
            .transport_state()
            .map_err(OsmosisError::Fluid)?
            .dynamic_viscosity_pa_s();
        Self::new((intrinsic_permeability_m2 / viscosity) / thickness_m)
    }
}

/// Concentration-based solute membrane permeance [m/s], not diffusivity [m²/s].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SolutePermeance(f64);
impl SolutePermeance {
    /// Zero excludes diffusive permeation; solvent drag is independently selective.
    pub fn new(value_m_s: f64) -> Result<Self> {
        if nonnegative(value_m_s) {
            Ok(Self(value_m_s))
        } else {
            Err(OsmosisError::Invalid(
                "solute permeance must be finite and nonnegative m/s",
            ))
        }
    }
    /// Numerical SI value [m/s].
    pub const fn value(self) -> f64 {
        self.0
    }
}

/// Species-specific osmotic reflection coefficient: 0 freely convected, 1 excluded.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReflectionCoefficient(f64);
impl ReflectionCoefficient {
    /// Admit a finite coefficient in the closed physical interval [0, 1].
    pub fn new(value: f64) -> Result<Self> {
        if nonnegative(value) && value <= 1.0 {
            Ok(Self(value))
        } else {
            Err(OsmosisError::Invalid(
                "reflection coefficient must be in [0, 1]",
            ))
        }
    }
    /// Dimensionless coefficient.
    pub const fn value(self) -> f64 {
        self.0
    }
}

/// Prepared constant-coefficient membrane on one ordered independent-solute layout.
#[derive(Clone, Debug, PartialEq)]
pub struct MembraneTransport {
    hydraulic: HydraulicPermeance,
    solutes: Box<[(ReflectionCoefficient, SolutePermeance)]>,
}

/// Common left-to-right solvent result; caller-owned output holds molar solute fluxes.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MembraneFlux {
    /// Solvent volume per membrane area and second [m/s], positive left to right.
    pub volume_flux_m_s: f64,
    /// `RT sum(sigma_i (c_left_i-c_right_i))` [Pa].
    pub effective_osmotic_pressure_drop_pa: f64,
    /// Hydrostatic minus selectively reflected osmotic pressure drop [Pa].
    pub driving_pressure_pa: f64,
}

impl MembraneTransport {
    /// Allocate the fixed species layout once; numerical constitutive inputs
    /// do not imply source support over a later changing solution state.
    pub fn new(
        hydraulic: HydraulicPermeance,
        solutes: Vec<(ReflectionCoefficient, SolutePermeance)>,
    ) -> Result<Self> {
        layout(solutes.len())?;
        Ok(Self {
            hydraulic,
            solutes: solutes.into_boxed_slice(),
        })
    }

    /// Original hydraulic permeance.
    pub const fn hydraulic_permeance(&self) -> HydraulicPermeance {
        self.hydraulic
    }
    /// Original per-species `(reflection, concentration permeance)` layout.
    pub fn solutes(&self) -> &[(ReflectionCoefficient, SolutePermeance)] {
        &self.solutes
    }

    /// Nonnegative donor coefficients [m/s] for `Js=k_left*c_left-k_right*c_right`.
    /// The volume flux has the same left-to-right orientation as [`Self::flux_into`].
    /// This supports positivity-preserving finite-inventory consumers without
    /// reverse-engineering the prepared constitutive law.
    pub fn solute_stream_coefficients(
        &self,
        index: usize,
        volume_flux_m_s: f64,
    ) -> Result<(f64, f64)> {
        let &(reflection, permeability) = self.solutes.get(index).ok_or(OsmosisError::Invalid(
            "solute index is outside the prepared layout",
        ))?;
        let advective = finite(
            (1.0 - reflection.0) * volume_flux_m_s,
            "solute convective velocity is not representable",
        )?;
        let exchange = diffusive_exchange(permeability.0, advective)?;
        if advective >= 0.0 {
            Ok((
                finite(
                    exchange + advective,
                    "donor coefficient is not representable",
                )?,
                exchange,
            ))
        } else {
            Ok((
                exchange,
                finite(
                    exchange - advective,
                    "donor coefficient is not representable",
                )?,
            ))
        }
    }

    /// Evaluate the solvent driving law without evaluating solute transport.
    ///
    /// Pressures may share any consistent gauge. The temperature and species
    /// layout must be common; fluxes are per membrane area, so geometry is
    /// applied exactly once by the caller. Pure water moves toward higher
    /// reflected solute concentration at equal pressure; sufficient applied
    /// pressure reverses that flow. Finite-volume solvers that already obtained
    /// the solute donor coefficients can evaluate their scalar water residual
    /// here without repeating the membrane exponential fits or staging solute fluxes.
    pub fn water_flux(
        &self,
        temperature_k: f64,
        pressure_left_pa: f64,
        pressure_right_pa: f64,
        concentration_left_mol_m3: &[f64],
        concentration_right_mol_m3: &[f64],
    ) -> Result<MembraneFlux> {
        let count = self.solutes.len();
        shape(count, concentration_left_mol_m3.len())?;
        shape(count, concentration_right_mol_m3.len())?;
        let rt = thermal_scale(temperature_k)?;
        if !pressure_left_pa.is_finite() || !pressure_right_pa.is_finite() {
            return Err(OsmosisError::Invalid(
                "both pressures must be finite pascals in a common gauge",
            ));
        }
        let mut reflected = Sum::default();
        for ((&(reflection, _), &left), &right) in self
            .solutes
            .iter()
            .zip(concentration_left_mol_m3)
            .zip(concentration_right_mol_m3)
        {
            if !nonnegative(left) || !nonnegative(right) {
                return Err(OsmosisError::Invalid(
                    "both concentrations must be finite and nonnegative",
                ));
            }
            reflected.add(reflection.0 * (left - right));
        }
        let effective_osmotic_pressure_drop_pa = finite(
            rt * reflected.finish()?,
            "reflected osmotic pressure is not representable",
        )?;
        let driving_pressure_pa = finite(
            (pressure_left_pa - pressure_right_pa) - effective_osmotic_pressure_drop_pa,
            "transmembrane driving pressure is not representable",
        )?;
        let volume_flux_m_s = finite(
            self.hydraulic.0 * driving_pressure_pa,
            "solvent flux is not representable",
        )?;
        Ok(MembraneFlux {
            volume_flux_m_s,
            effective_osmotic_pressure_drop_pa,
            driving_pressure_pa,
        })
    }

    /// Evaluate isothermal solvent and solute fluxes without allocation.
    ///
    /// Uses the same [`Self::water_flux`] solvent law and the fitted solute
    /// transport. All fluxes are left to right per membrane area. Failure leaves
    /// the complete caller-owned `solute_flux_mol_m2_s` output unchanged.
    pub fn flux_into(
        &self,
        temperature_k: f64,
        pressure_left_pa: f64,
        pressure_right_pa: f64,
        concentration_left_mol_m3: &[f64],
        concentration_right_mol_m3: &[f64],
        solute_flux_mol_m2_s: &mut [f64],
    ) -> Result<MembraneFlux> {
        let count = self.solutes.len();
        shape(count, solute_flux_mol_m2_s.len())?;
        let water = self.water_flux(
            temperature_k,
            pressure_left_pa,
            pressure_right_pa,
            concentration_left_mol_m3,
            concentration_right_mol_m3,
        )?;
        let mut output = [0.0; MAX_OSMOTIC_SOLUTES];
        for (index, &(reflection, permeability)) in self.solutes.iter().enumerate() {
            let left = concentration_left_mol_m3[index];
            let right = concentration_right_mol_m3[index];
            let advective = finite(
                (1.0 - reflection.0) * water.volume_flux_m_s,
                "solute convective velocity is not representable",
            )?;
            let exchange = diffusive_exchange(permeability.0, advective)?;
            // This form avoids subtracting almost equal large donor streams
            // for a small concentration gradient at low membrane Peclet number.
            output[index] = finite(
                advective * if advective >= 0.0 { left } else { right } + exchange * (left - right),
                "solute flux is not representable",
            )?;
        }
        solute_flux_mol_m2_s.copy_from_slice(&output[..count]);
        Ok(water)
    }
}

// Symmetric exchange P*B(|u|/P), with B(x)=x/(exp(x)-1). The remaining
// |u| is a one-sided donor stream. This is the exact steady solution of
// Js=u*c-P*dc/ds on normalized membrane coordinate s in [0,1].
fn diffusive_exchange(permeability: f64, advective: f64) -> Result<f64> {
    if !advective.is_finite() {
        return Err(OsmosisError::Invalid("convective flux must be finite"));
    }
    if permeability == 0.0 {
        return Ok(0.0);
    }
    let speed = advective.abs();
    if speed == 0.0 {
        return Ok(permeability);
    }
    let peclet = speed / permeability;
    let exchange = if peclet < 1e-4 {
        let square = peclet * peclet;
        permeability * (1.0 - 0.5 * peclet + square * (1.0 / 12.0 - square / 720.0))
    } else {
        // Negative exponents avoid overflow even at an infinite Peclet ratio.
        speed * det::exp(-peclet) / -det::expm1(-peclet)
    };
    finite(exchange, "diffusive donor coefficient is not representable")
}

#[cfg(test)]
mod tests;
