//! Source-bound collisional DC response of a weakly ionized ideal gas.
//!
//! The original Gibbs-equilibrium owner supplies every carrier abundance,
//! charge, mass and temperature. A separately sourced momentum-relaxation rate
//! supplies drag against a prescribed massive background. Thus
//! `sigma = sum_s n_s (z_s e)^2 / (m_s nu_s)` and `J = sigma E`.
//! The rates are effective momentum-transfer rates, not total collision counts
//! or freely moving electron-ion collision frequencies. Interspecies drag,
//! magnetic/Hall response, finite-frequency inertia, sheaths, diffusion currents,
//! non-LTE breakdown and radiation require different constitutive models.
//!
//! Reference for the relaxation law and its frequency dependence: R. Fitzpatrick,
//! *Dielectric constant of a collisional plasma*, equations 1202–1204,
//! <https://farside.ph.utexas.edu/teaching/em/lectures/node103.html>.
//! The DC, local, unmagnetized, weak-drift restriction is part of this law.

use crate::state_point::ResolvedMaterialStatePoint;
use fs_matdb::PropertyKey;
use fs_math::det::{exp, ln, sqrt};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, Pressure, QuantitySpec, Temperature};
use fs_thermochem::SpeciesId;
use fs_thermochem::equilibrium::composition::{
    IdealGasEquilibriumSolutionV1, IdealGasEquilibriumV1,
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
    /// Cancellation before any snapshot is published.
    Cancelled,
}
impl core::fmt::Display for IonizedGasError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "ionized-gas DC material: {self:?}")
    }
}
impl std::error::Error for IonizedGasError {}

/// One effective carrier/background momentum-relaxation source at an actual
/// material state. The complete original card identity and usage receipts stay
/// attached. Source applicability includes its authored mixture/background;
/// this declaration does not infer a collision cross section from a species name.
#[derive(Clone, Debug, PartialEq)]
pub struct CarrierMomentumRelaxation {
    species: SpeciesId,
    material: ResolvedMaterialStatePoint,
    property: PropertyKey,
    rate_per_s: f64,
}
impl CarrierMomentumRelaxation {
    /// Adopt an explicitly selected positive scalar rate from an existing
    /// source-resolved bundle. Coordinates `T` [K] and `p` [Pa] must be supplied
    /// and later match the original equilibrium state exactly. Field/rate axes
    /// cannot be silently frozen as environmental conditions.
    pub fn new(
        species: SpeciesId,
        material: ResolvedMaterialStatePoint,
        property: PropertyKey,
    ) -> Result<Self, IonizedGasError> {
        let invalid = IonizedGasError::Invalid;
        let source = material.property_by_key(&property).ok_or(invalid(
            "explicit carrier momentum-relaxation source is absent",
        ))?;
        if property.quantity() != QuantitySpec::dimensional(MOMENTUM_RELAXATION_RATE_DIMS)
            || !source.is_static_at_fixed_external_conditions()
            || !positive(source.value_si())
            || !material.state_coordinate("T").is_some_and(positive)
            || !material.state_coordinate("p").is_some_and(positive)
        {
            return Err(invalid(
                "require a positive static momentum-relaxation source and actual T,p",
            ));
        }
        let rate_per_s = source.value_si();
        for (axis, dims, kind) in [
            ("T", Temperature::DIMS, QuantityKind::AbsoluteTemperature),
            ("p", Pressure::DIMS, QuantityKind::Pressure),
        ] {
            let Some(quantity) = material.axis_quantities().get(axis) else {
                return Err(invalid(
                    "collision state requires typed T and p coordinates",
                ));
            };
            if *quantity != QuantitySpec::dimensional(dims)
                && *quantity != QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
            {
                return Err(invalid(
                    "collision state T,p coordinate meanings do not match thermodynamics",
                ));
            }
        }
        Ok(Self {
            species,
            material,
            property,
            rate_per_s,
        })
    }
    /// Exact source-bound carrier identity.
    pub const fn species(&self) -> &SpeciesId {
        &self.species
    }
    /// All original material sources, query coordinates and usage receipts.
    pub const fn material(&self) -> &ResolvedMaterialStatePoint {
        &self.material
    }
    /// Original selected property, without a name-based fallback.
    pub const fn property(&self) -> &PropertyKey {
        &self.property
    }
    /// Effective momentum relaxation rate [1/s] at the retained state.
    pub const fn rate_per_s(&self) -> f64 {
        self.rate_per_s
    }
}

/// Caller-declared constitutive validity limits. All four limits must be
/// positive and strictly less than one; they are explicit domain choices,
/// not experimentally validated defaults or estimates of discretization error.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IonizedGasDcLimits {
    /// Upper bound on total charged-particle mole fraction; the remainder is
    /// the prescribed neutral background for the independent drag closures.
    pub maximum_charged_mole_fraction: f64,
    /// Upper bound on `z_max^2 e^2/(4 pi eps0 a kT)`, using total charged
    /// density and Wigner-Seitz radius a; keeps this ideal mixture weakly coupled.
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

/// Immutable thermochemical and collisional material snapshot for the existing
/// spatial current, Joule-heating and circuit-resistance owners. Re-equilibrate
/// and resolve collision sources again when the physical state changes.
#[derive(Clone, Debug, PartialEq)]
pub struct IonizedGasConductivity {
    equilibrium: IdealGasEquilibriumSolutionV1,
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
        mut rates: Vec<CarrierMomentumRelaxation>,
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
        if rates.len() > equilibrium.species().len() {
            return Err(invalid("more collision sources than equilibrium species"));
        }
        rates.sort_by(|a, b| a.species.cmp(&b.species));
        if rates.windows(2).any(|w| w[0].species == w[1].species) {
            return Err(invalid("duplicate carrier momentum-relaxation sources"));
        }
        let t = equilibrium.temperature().value();
        let p = equilibrium.pressure().value();
        let kt = BOLTZMANN_J_K * t;
        let total_density = p / kt;
        if !positive(total_density) {
            return Err(invalid("unrepresentable ideal particle density"));
        }
        let rate_count = rates.len();
        let mut carriers = Vec::with_capacity(rate_count);
        let mut rates = rates.into_iter();
        let mut conductivity = 0.0;
        let mut charged_fraction = 0.0;
        let mut largest_charge = 0.0_f64;
        let mut maximum_phase_space = 0.0_f64;
        let mut maximum_field = f64::INFINITY;
        for (index, species) in equilibrium.species().iter().enumerate() {
            poll(&mut cancelled)?;
            let z = chemistry
                .charge_vector()
                .get(index)
                .ok_or(invalid("missing original charge"))?;
            if z == 0 || species.log_mole_fraction.is_none() {
                continue;
            }
            if !positive(species.mole_fraction) {
                return Err(IonizedGasError::UnrepresentableCarrier(
                    species.species.clone(),
                ));
            }
            let relaxation = rates
                .next()
                .filter(|rate| rate.species == species.species)
                .ok_or_else(|| IonizedGasError::MissingRate(species.species.clone()))?;
            if relaxation.material.state_coordinate("T") != Some(t)
                || relaxation.material.state_coordinate("p") != Some(p)
            {
                return Err(invalid(
                    "collision source and equilibrium must share actual T,p",
                ));
            }
            let n = species.mole_fraction * total_density;
            let m = chemistry.models()[index].molar_mass().value() / AVOGADRO_PER_MOL;
            let charge = (z as f64).abs() * ELEMENTARY_CHARGE_C;
            let nu = relaxation.rate_per_s;
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
            charged_fraction += species.mole_fraction;
            largest_charge = largest_charge.max((z as f64).abs());
            maximum_phase_space = maximum_phase_space.max(phase_space);
            maximum_field = maximum_field.min(max_e);
            carriers.push(IonizedGasCarrier {
                species: species.species.clone(),
                charge_number: z,
                number_density_per_m3: n,
                mass_kg: m,
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
        let inverse_radius = exp(ln(4.0 * core::f64::consts::PI / 3.0) / 3.0
            + (ln(charged_fraction) + ln(total_density)) / 3.0);
        let coupling = (largest_charge * ELEMENTARY_CHARGE_C).powi(2) * inverse_radius
            / (4.0 * core::f64::consts::PI * VACUUM_PERMITTIVITY_F_M * kt);
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
        poll(&mut cancelled)?;
        Ok(Self {
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
    /// Original accepted composition, temperature, calorics and source receipts.
    pub const fn equilibrium(&self) -> &IdealGasEquilibriumSolutionV1 {
        &self.equilibrium
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
