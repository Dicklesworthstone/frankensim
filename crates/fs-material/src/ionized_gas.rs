//! Source-bound collisional DC response of a weakly ionized ideal gas.
//!
//! The original source-bound state supplies every carrier abundance, charge,
//! mass and temperature, including frozen finite-rate populations. A separately
//! sourced momentum-relaxation rate
//! supplies drag against a prescribed massive background. Thus
//! `sigma = sum_s n_s (z_s e)^2 / (m_s nu_s)` and `J = sigma E`.
//! The rates are effective momentum-transfer rates, not total collision counts
//! or freely moving electron-ion collision frequencies. Interspecies drag,
//! magnetic/Hall response, finite-frequency inertia, sheaths, diffusion currents,
//! electrical breakdown and radiation require different constitutive models.
//!
//! Reference for the relaxation law and its frequency dependence: R. Fitzpatrick,
//! *Dielectric constant of a collisional plasma*, equations 1202–1204,
//! <https://farside.ph.utexas.edu/teaching/em/lectures/node103.html>.
//! The DC, local, unmagnetized, weak-drift restriction is part of this law.

use fs_math::det::{exp, ln, sqrt};
use fs_qty::{Dims, Volume};
use fs_thermochem::equilibrium::composition::state::{
    FrozenIonizedGasStateErrorV1, FrozenIonizedGasStateV1,
};
use fs_thermochem::equilibrium::composition::{
    IdealGasEquilibriumSolutionV1, IdealGasEquilibriumV1,
};
use fs_thermochem::{SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K};

mod relaxation;
pub use relaxation::{
    CarrierMomentumRelaxation, CarrierMomentumSource, ElectronNeutralMomentumContribution,
};

/// Effective momentum relaxation rate dimensions, reciprocal seconds.
/// This is a collision rate, without an angular/cyclic-frequency conversion.
pub const MOMENTUM_RELAXATION_RATE_DIMS: Dims = Dims([0, 0, -1, 0, 0, 0]);
/// Exact elementary charge under the SI definition [C].
pub const ELEMENTARY_CHARGE_C: f64 = 1.602_176_634e-19;
const AVOGADRO_PER_MOL: f64 = 6.022_140_76e23;
const BOLTZMANN_J_K: f64 = 1.380_649e-23;
const PLANCK_J_S: f64 = 6.626_070_15e-34;
// CODATA 2022 value. The constitutive law is not an uncertainty enclosure.
const VACUUM_PERMITTIVITY_F_M: f64 = 8.854_187_818_8e-12;

/// Invalid source, unsupported state or failed local transport admission.
#[derive(Clone, Debug, PartialEq)]
pub enum IonizedGasError {
    /// Invalid controls, inconsistent source conditions or arithmetic range loss.
    Invalid(&'static str),
    /// A carrier has no exactly matching momentum-relaxation source.
    MissingRate(SpeciesId),
    /// A logarithmically present charged trace cannot be represented directly.
    UnrepresentableCarrier(SpeciesId),
    /// Original frozen-composition/source admission refusal.
    State(FrozenIonizedGasStateErrorV1),
    /// Cancellation before any snapshot is published.
    Cancelled,
}
impl core::fmt::Display for IonizedGasError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "ionized-gas DC material: {self:?}")
    }
}
impl std::error::Error for IonizedGasError {}

/// Caller-declared constitutive validity limits. All four limits must be
/// positive and strictly less than one; they are explicit domain choices,
/// not experimentally validated defaults or estimates of discretization error.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IonizedGasDcLimits {
    /// Upper bound on total charged-particle mole fraction; the remainder is
    /// the prescribed neutral background for the independent drag closures.
    pub maximum_charged_mole_fraction: f64,
    /// Upper bound on `z_max^2 e^2/(4 pi eps0 a kT_min)`, using total charged
    /// density, Wigner-Seitz radius a and the coldest present carrier temperature.
    pub maximum_coulomb_coupling: f64,
    /// Upper bound on `n_s [h/sqrt(2 pi m_s kT)]^3` for every charged species.
    /// No degeneracy/spin factor is inferred from a species name.
    pub maximum_thermal_phase_space_density: f64,
    /// Actual drift speed must be below this fraction of `sqrt(kT/m_s)`
    /// for every carrier. Checked again on the electric field from the PDE.
    pub maximum_drift_thermal_ratio: f64,
}

/// One canonical charge carrier at the accepted thermodynamic state.
#[derive(Clone, Debug, PartialEq)]
pub struct IonizedGasCarrier {
    species: SpeciesId,
    charge_number: i128,
    number_density_per_m3: f64,
    mass_kg: f64,
    temperature_k: f64,
    is_electron: bool,
    conductivity_s_m: f64,
    maximum_electric_field_v_m: f64,
    relaxation: CarrierMomentumRelaxation,
}
impl IonizedGasCarrier {
    /// Original chemical species identity.
    pub const fn species(&self) -> &SpeciesId {
        &self.species
    }
    /// Signed integer charge in elementary-charge units.
    pub const fn charge_number(&self) -> i128 {
        self.charge_number
    }
    /// Actual carrier number density [1/m³], including its equilibrium fraction.
    pub const fn number_density_per_m3(&self) -> f64 {
        self.number_density_per_m3
    }
    /// Original source molar mass divided by the exact Avogadro constant [kg].
    pub const fn mass_kg(&self) -> f64 {
        self.mass_kg
    }
    /// Actual kinetic temperature [K]: Te for the structural electron row,
    /// Th for ions. A finite-rate composition is not replaced by equilibrium.
    pub const fn temperature_k(&self) -> f64 {
        self.temperature_k
    }
    /// True only for the original nucleus-free charge-minus-one source row.
    pub const fn is_electron(&self) -> bool {
        self.is_electron
    }
    /// This carrier's positive contribution to total conductivity [S/m].
    pub const fn conductivity_s_m(&self) -> f64 {
        self.conductivity_s_m
    }
    /// This carrier's weak-drift limit on the actual electric-field magnitude.
    pub const fn maximum_electric_field_v_m(&self) -> f64 {
        self.maximum_electric_field_v_m
    }
    /// Original effective momentum-relaxation source and complete state receipt.
    pub const fn relaxation(&self) -> &CarrierMomentumRelaxation {
        &self.relaxation
    }
}

/// Electrical work per unit volume supplied to the two thermal pools [W/m³].
/// The static weak-drift closure has no separately evolved drift kinetic energy.
/// Electron/heavy thermal exchange is a separate conservative energy process.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IonizedGasJouleHeating {
    electron_power_density_w_m3: f64,
    heavy_power_density_w_m3: f64,
}
impl IonizedGasJouleHeating {
    pub const fn electron_power_density_w_m3(&self) -> f64 {
        self.electron_power_density_w_m3
    }
    pub const fn heavy_power_density_w_m3(&self) -> f64 {
        self.heavy_power_density_w_m3
    }
    pub fn total_power_density_w_m3(&self) -> f64 {
        self.electron_power_density_w_m3 + self.heavy_power_density_w_m3
    }
}

/// Immutable thermochemical and collisional material snapshot for the existing
/// spatial current, Joule-heating and circuit-resistance owners. Construct a
/// new snapshot from actual populations and collision sources when they change.
#[derive(Clone, Debug, PartialEq)]
pub struct IonizedGasConductivity {
    state: FrozenIonizedGasStateV1,
    equilibrium: Option<IdealGasEquilibriumSolutionV1>,
    carriers: Vec<IonizedGasCarrier>,
    limits: IonizedGasDcLimits,
    conductivity_s_m: f64,
    maximum_electric_field_v_m: f64,
    charged_mole_fraction: f64,
    coulomb_coupling: f64,
    maximum_thermal_phase_space_density: f64,
}
impl IonizedGasConductivity {
    /// Resolve all charged particles from the SAME source-bound equilibrium.
    /// Every present ion/electron needs one rate; neutral or unused extra rates
    /// and duplicate identities refuse. This consumes no neutral-gas viscosity
    /// or Wilke/Eucken transport approximation. The caller supplies effective
    /// rates appropriate to the actual mixture and prescribed background.
    pub fn new(
        chemistry: &IdealGasEquilibriumV1,
        equilibrium: IdealGasEquilibriumSolutionV1,
        rates: Vec<CarrierMomentumRelaxation>,
        limits: IonizedGasDcLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, IonizedGasError> {
        let invalid = IonizedGasError::Invalid;
        poll(&mut cancelled)?;
        if !equilibrium.is_bound_to(chemistry) {
            return Err(invalid(
                "equilibrium must retain this exact elemental, charge and source-model owner",
            ));
        }
        let mut amounts = Vec::with_capacity(equilibrium.species().len());
        for species in equilibrium.species() {
            poll(&mut cancelled)?;
            if species.log_mole_fraction.is_some()
                && (!positive(species.mole_fraction) || !positive(species.amount.value()))
            {
                return Err(IonizedGasError::UnrepresentableCarrier(
                    species.species.clone(),
                ));
            }
            amounts.push(species.amount);
        }
        let t = equilibrium.temperature();
        let p = equilibrium.pressure().value();
        let volume = equilibrium.total_amount().value()
            * (UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * t.value() / p);
        let charge_tolerance =
            equilibrium.residuals().relative_charge_residual + 128.0 * f64::EPSILON;
        let state = FrozenIonizedGasStateV1::new(
            chemistry,
            &amounts,
            Volume::new(volume),
            t,
            t,
            charge_tolerance,
            &mut cancelled,
        )
        .map_err(IonizedGasError::State)?;
        if (state.pressure().value() / p - 1.0).abs() > 512.0 * f64::EPSILON {
            return Err(invalid(
                "equilibrium-to-frozen conversion changed physical pressure",
            ));
        }
        // The original accepted pressure remains the exact source-query
        // authority; reconstruction from NRT/V may differ by roundoff only.
        Self::build(state, Some(equilibrium), p, rates, limits, &mut cancelled)
    }

    /// Preserve actual frozen populations from a finite-rate or prescribed
    /// state, with independent electron and heavy temperatures. No equilibrium
    /// projection, Saha population or additional reaction heat is introduced.
    pub fn from_state(
        state: FrozenIonizedGasStateV1,
        rates: Vec<CarrierMomentumRelaxation>,
        limits: IonizedGasDcLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, IonizedGasError> {
        let pressure = state.pressure().value();
        Self::build(state, None, pressure, rates, limits, &mut cancelled)
    }

    fn build(
        state: FrozenIonizedGasStateV1,
        equilibrium: Option<IdealGasEquilibriumSolutionV1>,
        source_pressure_pa: f64,
        mut rates: Vec<CarrierMomentumRelaxation>,
        limits: IonizedGasDcLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, IonizedGasError> {
        let invalid = IonizedGasError::Invalid;
        poll(cancelled)?;
        for x in [
            limits.maximum_charged_mole_fraction,
            limits.maximum_coulomb_coupling,
            limits.maximum_thermal_phase_space_density,
            limits.maximum_drift_thermal_ratio,
        ] {
            if !positive(x) || x >= 1.0 {
                return Err(invalid("transport limits must lie strictly inside (0,1)"));
            }
        }
        let chemistry = state.thermodynamics();
        if rates.len() > state.amounts().len() {
            return Err(invalid("more collision sources than original species"));
        }
        rates.sort_by(|a, b| a.species().cmp(b.species()));
        if rates.windows(2).any(|w| w[0].species() == w[1].species()) {
            return Err(invalid("duplicate carrier momentum-relaxation sources"));
        }
        let mut total_density = 0.0;
        for i in 0..state.amounts().len() {
            poll(cancelled)?;
            total_density += state
                .number_density_per_m3(i)
                .ok_or(invalid("missing species density"))?;
        }
        if !positive(total_density) {
            return Err(invalid("unrepresentable ideal particle density"));
        }
        let rate_count = rates.len();
        let mut carriers = Vec::with_capacity(rate_count);
        let mut rates = rates.into_iter();
        let mut conductivity = 0.0;
        let mut charged_density = 0.0;
        let mut largest_charge = 0.0_f64;
        let mut minimum_temperature = f64::INFINITY;
        let mut maximum_phase_space = 0.0_f64;
        let mut maximum_field = f64::INFINITY;
        for (index, model) in chemistry.models().iter().enumerate() {
            poll(cancelled)?;
            let z = chemistry
                .charge_vector()
                .get(index)
                .ok_or(invalid("missing original charge"))?;
            let n = state
                .number_density_per_m3(index)
                .ok_or(invalid("missing carrier density"))?;
            if z == 0 || n == 0.0 {
                continue;
            }
            let relaxation = rates
                .next()
                .filter(|rate| rate.species() == model.species())
                .ok_or_else(|| IonizedGasError::MissingRate(model.species().clone()))?;
            if !relaxation.matches_state(&state, source_pressure_pa) {
                return Err(invalid(
                    "collision source must retain the actual temperatures, pressure and any bound populations",
                ));
            }
            let t = state
                .species_temperature(index)
                .ok_or(invalid("missing carrier temperature"))?
                .value();
            let kt = BOLTZMANN_J_K * t;
            let m = model.molar_mass().value() / AVOGADRO_PER_MOL;
            let charge = (z as f64).abs() * ELEMENTARY_CHARGE_C;
            let nu = relaxation.rate_per_s();
            let sigma = (n * charge) * (charge / m) / nu;
            let thermal_speed = sqrt(kt / m);
            let max_e = (nu * m / charge) * thermal_speed * limits.maximum_drift_thermal_ratio;
            // Logarithms avoid falsely accepting a vanishing lambda^3 after
            // an intermediate cube underflows for a massive species.
            let log_lambda = ln(PLANCK_J_S) - 0.5 * ln(2.0 * core::f64::consts::PI * m * kt);
            let phase_space = exp(ln(n) + 3.0 * log_lambda);
            if !positive(n)
                || !positive(m)
                || !positive(sigma)
                || !positive(max_e)
                || !phase_space.is_finite()
            {
                return Err(invalid(
                    "carrier density, conductivity or validity scale is unrepresentable",
                ));
            }
            conductivity += sigma;
            charged_density += n;
            largest_charge = largest_charge.max((z as f64).abs());
            minimum_temperature = minimum_temperature.min(t);
            maximum_phase_space = maximum_phase_space.max(phase_space);
            maximum_field = maximum_field.min(max_e);
            carriers.push(IonizedGasCarrier {
                species: model.species().clone(),
                charge_number: z,
                number_density_per_m3: n,
                mass_kg: m,
                temperature_k: t,
                is_electron: Some(index) == state.electron_index(),
                conductivity_s_m: sigma,
                maximum_electric_field_v_m: max_e,
                relaxation,
            });
        }
        if carriers.is_empty() || carriers.len() != rate_count {
            return Err(invalid(
                "one rate is required for every present charged species and no others",
            ));
        }
        let charged_fraction = charged_density / total_density;
        let inverse_radius =
            exp((ln(4.0 * core::f64::consts::PI / 3.0) + ln(charged_density)) / 3.0);
        let coupling = (largest_charge * ELEMENTARY_CHARGE_C).powi(2) * inverse_radius
            / (4.0
                * core::f64::consts::PI
                * VACUUM_PERMITTIVITY_F_M
                * BOLTZMANN_J_K
                * minimum_temperature);
        if !positive(conductivity)
            || !positive(coupling)
            || charged_fraction > limits.maximum_charged_mole_fraction
            || coupling > limits.maximum_coulomb_coupling
            || maximum_phase_space > limits.maximum_thermal_phase_space_density
        {
            return Err(invalid(
                "state exceeds declared neutral-background, weak-coupling or classical-gas domain",
            ));
        }
        poll(cancelled)?;
        Ok(Self {
            state,
            equilibrium,
            carriers,
            limits,
            conductivity_s_m: conductivity,
            maximum_electric_field_v_m: maximum_field,
            charged_mole_fraction: charged_fraction,
            coulomb_coupling: coupling,
            maximum_thermal_phase_space_density: maximum_phase_space,
        })
    }
    /// Original accepted LTE receipt when constructed from equilibrium. Frozen
    /// finite-rate populations have no equilibrium receipt.
    pub const fn equilibrium(&self) -> Option<&IdealGasEquilibriumSolutionV1> {
        self.equilibrium.as_ref()
    }
    /// Exact frozen source bank, populations, both temperatures and calorics.
    pub const fn state(&self) -> &FrozenIonizedGasStateV1 {
        &self.state
    }
    /// Canonical present charge carriers with their independent collision sources.
    pub fn carriers(&self) -> &[IonizedGasCarrier] {
        &self.carriers
    }
    /// Sum of the actual carrier conductivities [S/m].
    pub const fn conductivity_s_m(&self) -> f64 {
        self.conductivity_s_m
    }
    /// The caller's retained constitutive-domain choices.
    pub const fn limits(&self) -> IonizedGasDcLimits {
        self.limits
    }
    /// Sum of all present charged-particle mole fractions, including electrons.
    pub const fn charged_mole_fraction(&self) -> f64 {
        self.charged_mole_fraction
    }
    /// Conservative charge-coupling screen defined by [`IonizedGasDcLimits`].
    pub const fn coulomb_coupling(&self) -> f64 {
        self.coulomb_coupling
    }
    /// Largest carrier `n lambda_thermal³` at the accepted state.
    pub const fn maximum_thermal_phase_space_density(&self) -> f64 {
        self.maximum_thermal_phase_space_density
    }
    /// The tightest carrier weak-drift limit, checked on the actual field.
    pub const fn maximum_electric_field_v_m(&self) -> f64 {
        self.maximum_electric_field_v_m
    }
    /// Partition the actual local `E.J` into electron and heavy-particle work
    /// using the original carriers, without equipartition or a species-name rule.
    /// Multiply by the same field-cell volume before supplying external powers
    /// to a conservative thermal-pool step. Do not also add the unsplit total.
    pub fn joule_heating_w_m3(
        &self,
        electric_field_v_m: [f64; 2],
    ) -> Result<IonizedGasJouleHeating, IonizedGasError> {
        if !self.admits_electric_field(electric_field_v_m) {
            return Err(IonizedGasError::Invalid(
                "Joule partition requires an admitted actual field",
            ));
        }
        let [x, y] = electric_field_v_m;
        let mut powers = [0.0, 0.0];
        for carrier in &self.carriers {
            let s = carrier.conductivity_s_m;
            let work = x * (s * x) + y * (s * y);
            if !work.is_finite() || work < 0.0 {
                return Err(IonizedGasError::Invalid(
                    "carrier electrical work is unrepresentable",
                ));
            }
            powers[usize::from(carrier.is_electron)] += work;
        }
        let total = powers[0] + powers[1];
        if !total.is_finite() || (electric_field_v_m != [0.0; 2] && total == 0.0) {
            return Err(IonizedGasError::Invalid(
                "total electrical work is unrepresentable",
            ));
        }
        Ok(IonizedGasJouleHeating {
            electron_power_density_w_m3: powers[1],
            heavy_power_density_w_m3: powers[0],
        })
    }
    /// A static, unmagnetized, spatially local field is admissible only within
    /// every carrier's declared weak-drift limit. This is not a breakdown model.
    pub fn admits_electric_field(&self, electric_field_v_m: [f64; 2]) -> bool {
        let [x, y] = electric_field_v_m;
        let scale = self.maximum_electric_field_v_m;
        x.is_finite() && y.is_finite() && (x / scale).hypot(y / scale) <= 1.0
    }
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), IonizedGasError> {
    if cancelled() {
        Err(IonizedGasError::Cancelled)
    } else {
        Ok(())
    }
}
