//! Finite heat at fixed total deformation on original caloric/elastic sources.
//!
//! At fixed F the restraints do no mechanical work. Heat changes the SAME
//! internal energy U=m*h0(T)+V0*(Psi-T*Psi_T), with h0 the explicitly declared
//! stress-free zero-pressure caloric reference. Its Newton capacity is
//! m*Cp(T)-V0*T*Psi_TT, including coefficient/expansion cross terms once.
//! Original Cp integrals are evaluated directly between temperatures, so their
//! arbitrary absolute enthalpy datum cannot erase a small heat transaction.
//!
//! This is a bounded numerical path: positive capacity and increasing U are
//! checked at the requested temperature spacing and every inverse trial. It
//! is not a global monotonicity/stability certificate or a CFL bound. No
//! ThermalLaw implementation, free mechanical motion, phase change, changing
//! reference mass or finite-pressure h-to-u conversion is inferred.

use super::{
    AdiabaticHyperelasticState, CaloricReference, FiniteSolidMaterial,
    FiniteThermoelasticSourceError, FiniteThermoelasticState, ReferenceDeformation,
    ResolvedHeatCapacityCurve, ResolvedNeoHookeanThermalCoefficients, ThermalReference,
    admit_caloric_source, expansion_jet_at_temperature, same_derivative,
};
use fs_material::phase::PhaseStateError;
use fs_thermal::enthalpy::ThermalStorageBasis;
use fs_thermochem::solution::caloric::BinaryCaloricErrorV1;
use fs_thermochem::solution::caloric::temperature::{self, TemperatureControl, TemperatureSample};
use std::sync::Arc;

/// Shared inverse controls plus whole-call bounded temperature-path admission.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteThermoelasticHeatSolve {
    pub thermal: TemperatureControl,
    /// Positive maximum spacing of the initial continuation samples [K].
    pub maximum_temperature_step_k: f64,
    /// 1..=4096 numerical material evaluations across continuation and inversion.
    /// Previously sampled temperatures reuse the exact same immutable state.
    pub maximum_state_evaluations: usize,
    /// Additional allowance for the finite endpoint Clausius bounds [J/K].
    pub entropy_tolerance_j_per_k: f64,
}

#[derive(Debug)]
pub enum FiniteThermoelasticHeatError {
    Invalid(&'static str),
    Source(FiniteThermoelasticSourceError),
    Caloric(PhaseStateError),
    Inverse(BinaryCaloricErrorV1),
    StateBudget {
        completed: usize,
    },
    NonMonotone {
        temperatures_k: [f64; 2],
        energy_change_j: f64,
    },
    ClausiusBounds {
        entropy_change_j_per_k: f64,
        bounds_j_per_k: [f64; 2],
    },
    Cancelled,
}
impl core::fmt::Display for FiniteThermoelasticHeatError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "fixed-deformation thermal transaction: {self:?}")
    }
}
impl std::error::Error for FiniteThermoelasticHeatError {}
impl From<FiniteThermoelasticSourceError> for FiniteThermoelasticHeatError {
    fn from(value: FiniteThermoelasticSourceError) -> Self {
        Self::Source(value)
    }
}
impl From<PhaseStateError> for FiniteThermoelasticHeatError {
    fn from(value: PhaseStateError) -> Self {
        Self::Caloric(value)
    }
}
impl From<BinaryCaloricErrorV1> for FiniteThermoelasticHeatError {
    fn from(value: BinaryCaloricErrorV1) -> Self {
        match value {
            BinaryCaloricErrorV1::Cancelled => Self::Cancelled,
            other => Self::Inverse(other),
        }
    }
}
type Result<T> = core::result::Result<T, FiniteThermoelasticHeatError>;
fn invalid(message: &'static str) -> FiniteThermoelasticHeatError {
    FiniteThermoelasticHeatError::Invalid(message)
}
fn finite(value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("thermal energy or capacity is unrepresentable"))
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(FiniteThermoelasticHeatError::Cancelled)
    } else {
        Ok(())
    }
}

#[derive(Clone, Debug)]
struct Source {
    material: FiniteSolidMaterial,
    coefficients: ResolvedNeoHookeanThermalCoefficients,
    caloric: ResolvedHeatCapacityCurve,
    caloric_reference: CaloricReference,
    before_thermal: ReferenceDeformation,
    total_f: [f64; 9],
    mass_kg: f64,
    volume_m3: f64,
}

/// Accepted temperature and complete native mechanical response at fixed F.
/// Every descendant retains this original immutable material/source inventory.
#[derive(Clone, Debug)]
pub struct FiniteSolidCaloricState {
    source: Arc<Source>,
    potential: FiniteThermoelasticState,
    base_specific_capacity_j_kg_k: f64,
    heat_capacity_j_per_k: f64,
}

/// Atomic signed heat transaction. No endpoint is published on any refusal.
#[derive(Clone, Debug)]
pub struct FiniteThermoelasticHeatBalance {
    final_state: FiniteSolidCaloricState,
    heat_j: f64,
    energy_change_j: f64,
    entropy_change_j_per_k: f64,
    temperature_bracket_k: [f64; 2],
    iterations: usize,
    state_evaluations: usize,
}

impl FiniteSolidMaterial {
    /// Prepare a closed, fully restrained material inventory at its actual T.
    /// The reference mass is explicit; reference density remains the original
    /// material's density, never the caloric chart's temperature-varying metadata.
    ///
    /// The same pinned coefficient interval must fit within the fully solid
    /// caloric and engineering-length paths. Expansion must be affine on its
    /// original complete prepared path; a kink would require its own transition
    /// entropy. Caloric Cp may vary across its original source knots.
    #[allow(clippy::too_many_arguments)]
    pub fn fixed_deformation_caloric_state(
        &self,
        total_f: [f64; 9],
        before_thermal: ReferenceDeformation,
        reference_mass_kg: f64,
        coefficients: &ResolvedNeoHookeanThermalCoefficients,
        caloric: &ResolvedHeatCapacityCurve,
        caloric_reference: CaloricReference,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FiniteSolidCaloricState> {
        poll(&mut cancelled)?;
        let Some(ThermalReference::Expansion(expansion)) = self.thermal_reference() else {
            return Err(invalid(
                "finite heat requires the original prepared engineering-length path",
            ));
        };
        let temperature = expansion.geometry().temperature_k();
        if coefficients.model() != self.law()
            || coefficients.material().card_identity() != self.density().card_identity()
            || coefficients.temperature_k() != temperature
            || !reference_mass_kg.is_finite()
            || reference_mass_kg <= 0.0
        {
            return Err(invalid(
                "same actual material/coefficient snapshot and positive original mass required",
            ));
        }
        admit_caloric_source(self, caloric, temperature, caloric_reference)?;
        let range = coefficients.temperature_range_k();
        let expansion_range = expansion.source().temperature_span_k();
        let knots = caloric.curve().knots();
        if range[0] < expansion_range[0]
            || range[1] > expansion_range[1]
            || range[0] < knots[0].temperature_k
            || range[1] > knots[knots.len() - 1].temperature_k
            || knots.iter().any(|knot| knot.liquid_mass_fraction != 0.0)
        {
            return Err(invalid(
                "original expansion and fully solid calorics must cover the coefficient interval",
            ));
        }
        let slopes = expansion
            .source()
            .at_temperature_with_tangent(temperature)
            .map_err(FiniteThermoelasticSourceError::Expansion)?
            .principal_stretch_derivatives_per_k();
        for knot in expansion.source().temperatures_k() {
            poll(&mut cancelled)?;
            let other = expansion
                .source()
                .at_temperature_with_tangent(knot)
                .map_err(FiniteThermoelasticSourceError::Expansion)?
                .principal_stretch_derivatives_per_k();
            if slopes
                .into_iter()
                .zip(other)
                .any(|(a, b)| !same_derivative(a, b))
            {
                return Err(invalid(
                    "finite heat cannot cross a kink in the original expansion potential",
                ));
            }
        }
        let volume_m3 = finite(reference_mass_kg / self.reference_density_kg_m3())?;
        if volume_m3 <= 0.0 {
            return Err(invalid("positive original volume is unrepresentable"));
        }
        let source = Arc::new(Source {
            material: self.clone(),
            coefficients: coefficients.clone(),
            caloric: caloric.clone(),
            caloric_reference,
            before_thermal,
            total_f,
            mass_kg: reference_mass_kg,
            volume_m3,
        });
        let state = evaluate(&source, temperature)?;
        poll(&mut cancelled)?;
        Ok(state)
    }
}

fn evaluate(source: &Arc<Source>, temperature: f64) -> Result<FiniteSolidCaloricState> {
    let Some(ThermalReference::Expansion(expansion)) = source.material.thermal_reference() else {
        unreachable!("prepared thermal source remains immutable")
    };
    let jet = expansion_jet_at_temperature(expansion, source.before_thermal, temperature, false)?;
    let potential = source
        .coefficients
        .thermoelastic_state_at_temperature(temperature, &source.total_f, jet)
        .map_err(FiniteThermoelasticSourceError::Coefficients)?;
    let base_specific_capacity_j_kg_k = source
        .caloric
        .secant_heat_capacity(temperature, temperature)?;
    let heat_capacity_j_per_k = finite(source.volume_m3.mul_add(
        potential.mechanical_heat_capacity_j_m3_k(),
        source.mass_kg * base_specific_capacity_j_kg_k,
    ))?;
    if heat_capacity_j_per_k <= 0.0 {
        return Err(invalid(
            "combined fixed-deformation capacity must be positive",
        ));
    }
    Ok(FiniteSolidCaloricState {
        source: Arc::clone(source),
        potential,
        base_specific_capacity_j_kg_k,
        heat_capacity_j_per_k,
    })
}

impl FiniteSolidCaloricState {
    /// U at fixed F, not open-stream enthalpy or an arbitrary pressure chart.
    pub const fn storage_basis(&self) -> ThermalStorageBasis {
        ThermalStorageBasis::InternalEnergy
    }
    /// Original anchor snapshot and density; the evolved constitutive response
    /// is returned by `thermodynamics`, not this anchor's held-temperature law.
    pub fn material(&self) -> &FiniteSolidMaterial {
        &self.source.material
    }
    pub fn coefficients(&self) -> &ResolvedNeoHookeanThermalCoefficients {
        &self.source.coefficients
    }
    pub fn caloric_source(&self) -> &ResolvedHeatCapacityCurve {
        &self.source.caloric
    }
    pub fn caloric_reference(&self) -> CaloricReference {
        self.source.caloric_reference
    }
    pub fn reference_mass_kg(&self) -> f64 {
        self.source.mass_kg
    }
    pub fn reference_volume_m3(&self) -> f64 {
        self.source.volume_m3
    }
    pub fn temperature_k(&self) -> f64 {
        self.potential.temperature_k()
    }
    pub const fn thermodynamics(&self) -> &FiniteThermoelasticState {
        &self.potential
    }
    pub const fn heat_capacity_j_per_k(&self) -> f64 {
        self.heat_capacity_j_per_k
    }
    pub const fn base_specific_capacity_j_kg_k(&self) -> f64 {
        self.base_specific_capacity_j_kg_k
    }

    /// Explicit source endpoint without claiming that heat or time has advanced.
    pub fn at_temperature(
        &self,
        temperature: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(&mut cancelled)?;
        let state = evaluate(&self.source, temperature)?;
        poll(&mut cancelled)?;
        Ok(state)
    }

    fn same_source(&self, other: &Self) -> Result<()> {
        if Arc::ptr_eq(&self.source, &other.source) {
            Ok(())
        } else {
            Err(invalid(
                "thermal endpoints require the same original fixed-F inventory",
            ))
        }
    }

    /// Direct original Cp integral plus the change of mechanical Psi-T*Psi_T.
    pub fn internal_energy_change_to(&self, other: &Self) -> Result<f64> {
        self.same_source(other)?;
        let base = self
            .source
            .caloric
            .specific_enthalpy_change(self.temperature_k(), other.temperature_k())?;
        finite(
            self.source.mass_kg.mul_add(
                base,
                self.source.volume_m3
                    * self
                        .potential
                        .internal_energy_increment_to(&other.potential)
                        .map_err(FiniteThermoelasticSourceError::Material)?,
            ),
        )
    }

    /// Original Cp/T integral plus the mechanical entropy change, once each.
    pub fn entropy_change_to(&self, other: &Self) -> Result<f64> {
        self.same_source(other)?;
        let base = self
            .source
            .caloric
            .specific_entropy_change(self.temperature_k(), other.temperature_k())?;
        finite(self.source.mass_kg.mul_add(
            base,
            self.source.volume_m3
                * (other.potential.entropy_j_m3_k() - self.potential.entropy_j_m3_k()),
        ))
    }

    /// Existing native full-compressible acoustic response at this actual T.
    /// A two-sided original source neighborhood is required for this derivative.
    pub fn adiabatic_acoustics(&self) -> Result<AdiabaticHyperelasticState> {
        let [low, high] = self.source.coefficients.temperature_range_k();
        if self.temperature_k() <= low || self.temperature_k() >= high {
            return Err(invalid(
                "adiabatic waves require a two-sided coefficient-source neighborhood",
            ));
        }
        self.potential
            .adiabatic_acoustics(
                self.source.material.reference_density_kg_m3(),
                self.source.material.reference_density_kg_m3() * self.base_specific_capacity_j_kg_k,
            )
            .map_err(FiniteThermoelasticSourceError::Material)
            .map_err(Into::into)
    }

    /// Solve signed heat into the same closed restrained inventory. Continuation
    /// covers the offered bracket at bounded spacing; original sample states are
    /// cached and the shared safeguarded inverse uses the complete actual c_F.
    /// Finite sampling is explicitly not a global stability proof. Zero heat
    /// returns this original state with zero work and no claim elsewhere in the
    /// bracket. Any source, capacity, budget or closure failure is atomic.
    pub fn add_heat(
        &self,
        heat_j: f64,
        controls: FiniteThermoelasticHeatSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FiniteThermoelasticHeatBalance> {
        poll(&mut cancelled)?;
        controls.thermal.validate()?;
        let [lo, hi] = controls.thermal.bracket;
        let range = self.source.coefficients.temperature_range_k();
        if !heat_j.is_finite()
            || lo < range[0]
            || hi > range[1]
            || self.temperature_k() < lo
            || self.temperature_k() > hi
            || !controls.maximum_temperature_step_k.is_finite()
            || controls.maximum_temperature_step_k <= 0.0
            || !(1..=4096).contains(&controls.maximum_state_evaluations)
            || !controls.entropy_tolerance_j_per_k.is_finite()
            || controls.entropy_tolerance_j_per_k < 0.0
        {
            return Err(invalid(
                "finite heat needs an original bracket, temperature spacing and finite work/tolerance budgets",
            ));
        }
        if heat_j == 0.0 {
            return Ok(FiniteThermoelasticHeatBalance {
                final_state: self.clone(),
                heat_j,
                energy_change_j: 0.0,
                entropy_change_j_per_k: 0.0,
                temperature_bracket_k: [self.temperature_k(); 2],
                iterations: 0,
                state_evaluations: 0,
            });
        }
        let divisions = ((hi - lo) / controls.maximum_temperature_step_k)
            .ceil()
            .max(1.0);
        if divisions + 1.0 > controls.maximum_state_evaluations as f64 {
            return Err(FiniteThermoelasticHeatError::StateBudget { completed: 0 });
        }
        let divisions = divisions as usize;
        let mut states = vec![self.clone()];
        let mut evaluations = 0;
        for index in 0..=divisions {
            let t = if index == divisions {
                hi
            } else {
                (hi - lo).mul_add(index as f64 / divisions as f64, lo)
            };
            sample(
                self,
                t,
                &mut states,
                &mut evaluations,
                controls.maximum_state_evaluations,
                &mut cancelled,
            )?;
        }
        let root = temperature::solve(
            controls.thermal,
            |t| -> Result<_> {
                let state = sample(
                    self,
                    t,
                    &mut states,
                    &mut evaluations,
                    controls.maximum_state_evaluations,
                    &mut cancelled,
                )?;
                let change_j = self.internal_energy_change_to(&state)?;
                Ok(TemperatureSample {
                    state,
                    change_j,
                    residual_j: finite(change_j - heat_j)?,
                })
            },
            |state| Some(state.heat_capacity_j_per_k),
        )?;
        let final_state = root.sample.state;
        let change = self.internal_energy_change_to(&final_state)?;
        let entropy = self.entropy_change_to(&final_state)?;
        if final_state.temperature_k() == self.temperature_k()
            || (change - heat_j).abs() > controls.thermal.energy_tolerance_j
        {
            return Err(invalid(
                "nonzero heat must resolve a changed temperature and close actual internal energy",
            ));
        }
        let a = change / self.temperature_k();
        let b = change / final_state.temperature_k();
        let bounds = [a.min(b), a.max(b)];
        let allowance = controls.entropy_tolerance_j_per_k
            + controls.thermal.energy_tolerance_j
                / self.temperature_k().min(final_state.temperature_k())
            + 64.0
                * f64::EPSILON
                * (bounds[0].abs().max(bounds[1].abs())
                    + self.source.volume_m3
                        * (self.potential.entropy_j_m3_k().abs()
                            + final_state.potential.entropy_j_m3_k().abs()));
        if !allowance.is_finite()
            || entropy < bounds[0] - allowance
            || entropy > bounds[1] + allowance
        {
            return Err(FiniteThermoelasticHeatError::ClausiusBounds {
                entropy_change_j_per_k: entropy,
                bounds_j_per_k: bounds,
            });
        }
        poll(&mut cancelled)?;
        Ok(FiniteThermoelasticHeatBalance {
            final_state,
            heat_j,
            energy_change_j: change,
            entropy_change_j_per_k: entropy,
            temperature_bracket_k: root.bracket,
            iterations: root.iterations,
            state_evaluations: evaluations,
        })
    }
}

fn sample(
    initial: &FiniteSolidCaloricState,
    temperature: f64,
    states: &mut Vec<FiniteSolidCaloricState>,
    evaluations: &mut usize,
    maximum: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<FiniteSolidCaloricState> {
    poll(cancelled)?;
    let index = match states.binary_search_by(|state| state.temperature_k().total_cmp(&temperature))
    {
        Ok(index) => return Ok(states[index].clone()),
        Err(index) => index,
    };
    if *evaluations >= maximum {
        return Err(FiniteThermoelasticHeatError::StateBudget {
            completed: *evaluations,
        });
    }
    let state = evaluate(&initial.source, temperature)?;
    *evaluations += 1;
    for (left, right) in index
        .checked_sub(1)
        .map(|i| (&states[i], &state))
        .into_iter()
        .chain(states.get(index).map(|next| (&state, next)))
    {
        let change = left.internal_energy_change_to(right)?;
        if change <= 0.0 {
            return Err(FiniteThermoelasticHeatError::NonMonotone {
                temperatures_k: [left.temperature_k(), right.temperature_k()],
                energy_change_j: change,
            });
        }
    }
    poll(cancelled)?;
    states.insert(index, state.clone());
    Ok(state)
}

impl FiniteThermoelasticHeatBalance {
    pub const fn final_state(&self) -> &FiniteSolidCaloricState {
        &self.final_state
    }
    pub const fn heat_j(&self) -> f64 {
        self.heat_j
    }
    pub const fn internal_energy_change_j(&self) -> f64 {
        self.energy_change_j
    }
    pub const fn entropy_change_j_per_k(&self) -> f64 {
        self.entropy_change_j_per_k
    }
    pub fn energy_residual_j(&self) -> f64 {
        self.energy_change_j - self.heat_j
    }
    pub const fn temperature_bracket_k(&self) -> [f64; 2] {
        self.temperature_bracket_k
    }
    pub const fn iterations(&self) -> usize {
        self.iterations
    }
    pub const fn state_evaluations(&self) -> usize {
        self.state_evaluations
    }
}
