//! First Chapman–Enskog transport from explicitly authored Lennard-Jones data.
//!
//! Neutral nonpolar dilute gases, Lorentz–Berthelot pair combining rules,
//! Neufeld–Janzen–Aziz collision integrals, Wilke viscosity mixing, and an
//! explicitly estimated Eucken conductivity closure. Molecular geometry,
//! rotational relaxation, polar corrections, Soret/Dufour coupling, dense-gas
//! corrections and experimental validation are not inferred from species names.
//!
//! References: Neufeld, Janzen and Aziz, J. Chem. Phys. 57, 1100–1102 (1972),
//! DOI 10.1063/1.1678363, Table I (including the Omega22 sinusoidal term);
//! <https://cantera.org/3.2/cxx/d5/de4/GasTransport_8cpp_source.html>;
//! <https://cantera.org/3.2/cxx/dd/d16/MixTransport_8cpp_source.html>.
//!
//! With molecular mass `m=M/N_A`, reduced pair mass `m_ij`, and collision
//! diameter `sigma`, the SI first-approximation laws are
//! `mu=5 sqrt(m k_B T/pi)/(16 sigma^2 Omega22)` and
//! `D_ij=3 k_B T sqrt(2 k_B T/(pi m_ij))/(16 p sigma_ij^2 Omega11)`.
//! The explicit Eucken estimate is `lambda_i=mu_i (Cp_i+5 R/4)/M_i`,
//! evaluated from the original NASA-9 molar heat capacity at the actual T.
//! Conductivity mixing averages its mole-weighted arithmetic and harmonic
//! means, as in Cantera's mixture transport implementation. In particular,
//! this reproduces the monatomic first-approximation Prandtl number 2/3,
//! without imposing a constant Prandtl number on a polyatomic mixture.

use crate::equilibrium::composition::IdealGasEquilibriumV1;
use crate::kinetics::{KineticsErrorV1, admit_caloric_interval};
use crate::{
    MolarThermalQuantityV1, SpeciesId, ThermochemErrorV1, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
};
use fs_qty::{DynViscosity, Length, Pressure, Qty, Temperature};

/// Fixed operation-tree version, including the complete Omega22 correlation.
pub const KINETIC_GAS_TRANSPORT_VERSION_V1: u32 = 1;
/// Explicit bound on component and quadratic pair work.
pub const MAX_KINETIC_TRANSPORT_SPECIES_V1: usize = 128;
/// Maximum opt-in source-fit Cp deficit, one part per million of `5R/2`.
/// This admission-policy ceiling is not an uncertainty bound or a correction
/// to the original source. A larger physical deficit needs another closure.
pub const MAX_TRANSLATIONAL_CP_DEFICIT_ALLOWANCE_V1: MolarThermalQuantityV1 =
    MolarThermalQuantityV1::new(2.5 * R * 1e-6);
/// Coherent-SI conductivity, W/(m K).
pub type ThermalConductivityQuantityV1 = Qty<1, 1, -3, -1, 0, 0>;
const BOLTZMANN: f64 = 1.380_649e-23;
const AVOGADRO: f64 = 6.022_140_76e23;
const PI: f64 = core::f64::consts::PI;

/// One original source declaration for a neutral 12–6 collision potential.
#[derive(Clone, Debug, PartialEq)]
pub struct LennardJonesSpeciesV1 {
    species: SpeciesId,
    diameter: Length,
    well_depth_over_k: Temperature,
}
impl LennardJonesSpeciesV1 {
    /// Retain finite positive collision diameter [m] and epsilon/kB [K].
    /// Molecular mass comes exclusively from the original thermodynamic bank.
    pub fn new(
        species: SpeciesId,
        diameter: Length,
        well_depth_over_k: Temperature,
    ) -> Result<Self, KineticTransportErrorV1> {
        if !positive(diameter.value()) || !positive(well_depth_over_k.value()) {
            return Err(KineticTransportErrorV1::Invalid(
                "collision diameter or well depth",
            ));
        }
        Ok(Self {
            species,
            diameter,
            well_depth_over_k,
        })
    }
    /// Original canonical chemical identity.
    pub const fn species(&self) -> &SpeciesId {
        &self.species
    }
    /// Original collision diameter in metres.
    pub const fn diameter(&self) -> Length {
        self.diameter
    }
    /// Original well depth divided by Boltzmann's constant, in kelvin.
    pub const fn well_depth_over_k(&self) -> Temperature {
        self.well_depth_over_k
    }
}

/// Caller-declared physical domain, intersected with original caloric support.
///
/// Every admitted species and pair must also satisfy 0.3 <= T* <= 100, the
/// original collision-correlation range. These explicit pressure bounds are a
/// source declaration, not a certificate that the gas is sufficiently dilute.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct KineticTransportDomainV1 {
    /// Minimum physical temperature [K].
    pub minimum_temperature_k: f64,
    /// Maximum physical temperature [K].
    pub maximum_temperature_k: f64,
    /// Minimum physical pressure [Pa].
    pub minimum_pressure_pa: f64,
    /// Maximum physical pressure [Pa].
    pub maximum_pressure_pa: f64,
}

/// Atomic source/property evaluation refusal; no partial result is returned.
#[derive(Clone, Debug, PartialEq)]
pub enum KineticTransportErrorV1 {
    /// Invalid dimensions, identity, authored data, composition or controls.
    Invalid(&'static str),
    /// Physical T or p exceeds the explicitly declared source domain.
    OutsideDomain,
    /// A pure or pair reduced temperature exceeds the collision fit domain.
    CollisionDomain {
        /// Evaluated dimensionless temperature, T/(epsilon/kB).
        reduced_temperature: f64,
    },
    /// The original NASA-9 model refused its common caloric path.
    CaloricPath(KineticsErrorV1),
    /// Original thermodynamic source evaluation refused.
    Thermochemistry {
        /// Canonical species index.
        species: usize,
        /// Original source refusal.
        source: ThermochemErrorV1,
    },
    /// Positive physical property or required arithmetic is unrepresentable.
    Arithmetic(&'static str),
    /// A bounded workspace/model reservation failed.
    Allocation,
    /// Cancellation before complete evaluation.
    Cancelled,
}
impl core::fmt::Display for KineticTransportErrorV1 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "dilute-gas kinetic transport: {self:?}")
    }
}
impl std::error::Error for KineticTransportErrorV1 {}

/// Source-bound collision model with precomputed composition-independent pairs.
#[derive(Clone, Debug, PartialEq)]
pub struct KineticGasTransportV1 {
    thermodynamics: IdealGasEquilibriumV1,
    model_id: String,
    model_version: u32,
    species: Vec<LennardJonesSpeciesV1>,
    domain: KineticTransportDomainV1,
    translational_cp_deficit_allowance: MolarThermalQuantityV1,
    viscosity_scales: Vec<f64>,
    diffusion_scales: Vec<f64>,
    pair_well_depths: Vec<f64>,
    mass_quarter_roots: Vec<f64>,
    wilke_denominators: Vec<f64>,
}
impl KineticGasTransportV1 {
    /// Bind one complete collision-data bank to the exact original species.
    /// Input collision records are canonicalized by identity. Missing, duplicate
    /// or extra species are refused; no air/default coefficients are supplied.
    pub fn new(
        thermodynamics: IdealGasEquilibriumV1,
        model_id: impl Into<String>,
        model_version: u32,
        mut species: Vec<LennardJonesSpeciesV1>,
        domain: KineticTransportDomainV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, KineticTransportErrorV1> {
        poll(&mut cancelled)?;
        let n = thermodynamics.models().len();
        let model_id = model_id.into();
        if !(1..=MAX_KINETIC_TRANSPORT_SPECIES_V1).contains(&n)
            || species.len() != n
            || model_id.is_empty()
            || model_id.len() > 128
            || model_version == 0
            || !model_id
                .bytes()
                .all(|c| c.is_ascii_alphanumeric() || b"-_.:/".contains(&c))
            || !positive(domain.minimum_temperature_k)
            || !positive(domain.maximum_temperature_k)
            || domain.minimum_temperature_k >= domain.maximum_temperature_k
            || !positive(domain.minimum_pressure_pa)
            || !positive(domain.maximum_pressure_pa)
            || domain.minimum_pressure_pa > domain.maximum_pressure_pa
        {
            return Err(KineticTransportErrorV1::Invalid(
                "model identity, complete species axis or physical domain",
            ));
        }
        species.sort_unstable_by(|a, b| a.species.cmp(&b.species));
        if species
            .iter()
            .zip(thermodynamics.models())
            .any(|(a, b)| a.species != *b.species())
        {
            return Err(KineticTransportErrorV1::Invalid(
                "collision and thermodynamic species identities differ",
            ));
        }
        let mut viscosity_scales = filled(n, 0.0)?;
        let mut diffusion_scales = filled(n * n, 0.0)?;
        let mut pair_well_depths = filled(n * n, 0.0)?;
        let mut mass_quarter_roots = filled(n * n, 0.0)?;
        let mut wilke_denominators = filled(n * n, 0.0)?;
        for (i, (data, source)) in species.iter().zip(thermodynamics.models()).enumerate() {
            poll(&mut cancelled)?;
            admit_caloric_interval(
                source,
                Temperature::new(domain.minimum_temperature_k),
                Temperature::new(domain.maximum_temperature_k),
                &mut cancelled,
            )
            .map_err(|e| match e {
                KineticsErrorV1::Cancelled => KineticTransportErrorV1::Cancelled,
                other => KineticTransportErrorV1::CaloricPath(other),
            })?;
            admit_reduced_temperature(
                domain.minimum_temperature_k / data.well_depth_over_k.value(),
            )?;
            admit_reduced_temperature(
                domain.maximum_temperature_k / data.well_depth_over_k.value(),
            )?;
            let molecular_mass = source.molar_mass().value() / AVOGADRO;
            let diameter_squared = data.diameter.value() * data.diameter.value();
            viscosity_scales[i] = (5.0 / 16.0)
                * fs_math::det::sqrt(molecular_mass * BOLTZMANN / PI)
                / diameter_squared;
            require_positive(viscosity_scales[i], "pure viscosity scale")?;
            for (j, other) in species.iter().enumerate() {
                let other_mass = thermodynamics.models()[j].molar_mass().value();
                let mass = source.molar_mass().value();
                mass_quarter_roots[i * n + j] =
                    fs_math::det::sqrt(fs_math::det::sqrt(other_mass / mass));
                wilke_denominators[i * n + j] = fs_math::det::sqrt(8.0 * (1.0 + mass / other_mass));
                require_positive(mass_quarter_roots[i * n + j], "Wilke molecular-mass ratio")?;
                require_positive(
                    wilke_denominators[i * n + j],
                    "Wilke molecular-mass denominator",
                )?;
                if j <= i {
                    continue;
                }
                let other_molecular_mass = other_mass / AVOGADRO;
                let (smaller, larger) = (
                    molecular_mass.min(other_molecular_mass),
                    molecular_mass.max(other_molecular_mass),
                );
                let reduced_mass = smaller / (1.0 + smaller / larger);
                let diameter = 0.5 * data.diameter.value() + 0.5 * other.diameter.value();
                let well_low = data
                    .well_depth_over_k
                    .value()
                    .min(other.well_depth_over_k.value());
                let well_high = data
                    .well_depth_over_k
                    .value()
                    .max(other.well_depth_over_k.value());
                let well_depth = well_high * fs_math::det::sqrt(well_low / well_high);
                let scale = (3.0 / 16.0)
                    * BOLTZMANN
                    * fs_math::det::sqrt((2.0 / PI) * (BOLTZMANN / reduced_mass))
                    / (diameter * diameter);
                require_positive(scale, "binary diffusion scale")?;
                require_positive(well_depth, "pair well depth")?;
                admit_reduced_temperature(domain.minimum_temperature_k / well_depth)?;
                admit_reduced_temperature(domain.maximum_temperature_k / well_depth)?;
                diffusion_scales[i * n + j] = scale;
                diffusion_scales[j * n + i] = scale;
                pair_well_depths[i * n + j] = well_depth;
                pair_well_depths[j * n + i] = well_depth;
            }
        }
        poll(&mut cancelled)?;
        Ok(Self {
            thermodynamics,
            model_id,
            model_version,
            species,
            domain,
            translational_cp_deficit_allowance: MolarThermalQuantityV1::new(0.0),
            viscosity_scales,
            diffusion_scales,
            pair_well_depths,
            mass_quarter_roots,
            wilke_denominators,
        })
    }
    /// Exact original elemental and NASA-9 source owner.
    pub const fn thermodynamics(&self) -> &IdealGasEquilibriumV1 {
        &self.thermodynamics
    }
    /// Authored transport model identity.
    pub fn model_id(&self) -> &str {
        &self.model_id
    }
    /// Authored transport data version, distinct from the evaluator version.
    pub const fn model_version(&self) -> u32 {
        self.model_version
    }
    /// Explicitly declared physical property domain.
    pub const fn domain(&self) -> KineticTransportDomainV1 {
        self.domain
    }
    /// Explicitly admit a bounded source-fit deficit below the classical `5R/2`
    /// molar Cp limit. The default is zero, retaining the original strict
    /// admission with its separate 32-epsilon arithmetic tolerance.
    ///
    /// The caller must justify this finite nonnegative allowance [J/(mol K)]
    /// against the original fit on the declared domain. It may not exceed
    /// [`MAX_TRANSLATIONAL_CP_DEFICIT_ALLOWANCE_V1`]. Every actual source Cp is
    /// still checked, including absent species; a larger deficit is refused.
    /// This control never clips Cp, edits the source or changes Eucken's use of
    /// the original Cp. It is part of this model's explicit declaration and
    /// equality, not empirical validation or a source-uncertainty certificate.
    pub fn with_translational_cp_deficit_allowance(
        mut self,
        allowance: MolarThermalQuantityV1,
    ) -> Result<Self, KineticTransportErrorV1> {
        if !allowance.value().is_finite()
            || allowance.value() < 0.0
            || allowance.value() > MAX_TRANSLATIONAL_CP_DEFICIT_ALLOWANCE_V1.value()
        {
            return Err(KineticTransportErrorV1::Invalid(
                "translational Cp source-fit deficit allowance",
            ));
        }
        self.translational_cp_deficit_allowance = allowance;
        Ok(self)
    }
    /// Authored source-fit Cp deficit budget [J/(mol K)]; zero is strict mode.
    pub const fn translational_cp_deficit_allowance(&self) -> MolarThermalQuantityV1 {
        self.translational_cp_deficit_allowance
    }
    /// Number of canonical species, including currently absent species.
    pub fn species_count(&self) -> usize {
        self.species.len()
    }
    /// Original collision records in the thermodynamic owner's canonical order.
    pub fn species_data(&self) -> &[LennardJonesSpeciesV1] {
        &self.species
    }
    /// Allocate bounded scratch once, to be reused across physical states.
    pub fn workspace(&self) -> Result<KineticTransportWorkspaceV1, KineticTransportErrorV1> {
        let n = self.species_count();
        Ok(KineticTransportWorkspaceV1 {
            pure_viscosities: filled(n, 0.0)?,
            pure_conductivities: filled(n, 0.0)?,
            binary_diffusivities: filled(n * n, 0.0)?,
        })
    }
    /// Evaluate all coefficients at one actual common T,p and mole composition.
    ///
    /// Fractions must be finite and nonnegative, with unit sum to floating-point
    /// summation precision. They are never normalized, clipped or seeded.
    /// Every original species source is evaluated, including possible products.
    /// Successful calls allocate nothing; a borrowed result exists only after
    /// all values and the final cancellation poll have been admitted.
    pub fn evaluate<'a>(
        &self,
        temperature: Temperature,
        pressure: Pressure,
        mole_fractions: &[f64],
        workspace: &'a mut KineticTransportWorkspaceV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<KineticTransportView<'a>, KineticTransportErrorV1> {
        poll(&mut cancelled)?;
        let n = self.species_count();
        if mole_fractions.len() != n
            || workspace.pure_viscosities.len() != n
            || workspace.pure_conductivities.len() != n
            || workspace.binary_diffusivities.len() != n * n
            || mole_fractions.iter().any(|x| !x.is_finite() || *x < 0.0)
        {
            return Err(KineticTransportErrorV1::Invalid(
                "canonical fraction or workspace dimensions",
            ));
        }
        let fraction_sum = sum(mole_fractions.iter().copied());
        if !fraction_sum.is_finite() || (fraction_sum - 1.0).abs() > 32.0 * f64::EPSILON * n as f64
        {
            return Err(KineticTransportErrorV1::Invalid(
                "mole fractions do not close to unity",
            ));
        }
        let (t, p) = (temperature.value(), pressure.value());
        if !t.is_finite()
            || !p.is_finite()
            || !(self.domain.minimum_temperature_k..=self.domain.maximum_temperature_k).contains(&t)
            || !(self.domain.minimum_pressure_pa..=self.domain.maximum_pressure_pa).contains(&p)
        {
            return Err(KineticTransportErrorV1::OutsideDomain);
        }
        let sqrt_t = fs_math::det::sqrt(t);
        for i in 0..n {
            poll(&mut cancelled)?;
            let omega = collision_integral_22(t / self.species[i].well_depth_over_k.value())?;
            let mu = self.viscosity_scales[i] * (sqrt_t / omega);
            require_positive(mu, "species viscosity")?;
            let source = &self.thermodynamics.models()[i];
            let cp = source
                .evaluate_properties(temperature)
                .map_err(|source| KineticTransportErrorV1::Thermochemistry { species: i, source })?
                .cp()
                .value();
            // Keep floating-point tolerance separate from an explicitly
            // declared source-fit deficit. The Eucken closure below always
            // consumes the original Cp, without a translational floor.
            if cp < 2.5 * R * (1.0 - 32.0 * f64::EPSILON)
                && 2.5 * R - cp > self.translational_cp_deficit_allowance.value()
            {
                return Err(KineticTransportErrorV1::Invalid(
                    "source Cp is below the classical translational limit",
                ));
            }
            let conductivity = (mu / source.molar_mass().value()) * (cp + 1.25 * R);
            require_positive(conductivity, "Eucken species conductivity")?;
            workspace.pure_viscosities[i] = mu;
            workspace.pure_conductivities[i] = conductivity;
            workspace.binary_diffusivities[i * n + i] = 0.0;
            for j in i + 1..n {
                let omega = collision_integral_11(t / self.pair_well_depths[i * n + j])?;
                let d = ((self.diffusion_scales[i * n + j] * sqrt_t) * (t / p)) / omega;
                require_positive(d, "binary diffusion coefficient")?;
                workspace.binary_diffusivities[i * n + j] = d;
                workspace.binary_diffusivities[j * n + i] = d;
            }
        }
        let mut mixed_mu = Compensated::default();
        let mut conductivity_arithmetic = Compensated::default();
        let mut conductivity_reciprocal = Compensated::default();
        let mut maximum_d = 0.0_f64;
        for i in 0..n {
            poll(&mut cancelled)?;
            maximum_d = maximum_d.max(
                workspace.binary_diffusivities[i * n..(i + 1) * n]
                    .iter()
                    .copied()
                    .fold(0.0, f64::max),
            );
            if mole_fractions[i] == 0.0 {
                continue;
            }
            let mut denominator = Compensated::default();
            for j in 0..n {
                if mole_fractions[j] == 0.0 {
                    continue;
                }
                let ratio = workspace.pure_viscosities[i] / workspace.pure_viscosities[j];
                require_positive(ratio, "Wilke pure-viscosity ratio")?;
                let base = 1.0 + fs_math::det::sqrt(ratio) * self.mass_quarter_roots[i * n + j];
                let phi = (base * base) / self.wilke_denominators[i * n + j];
                require_positive(phi, "Wilke pair factor")?;
                denominator.add(mole_fractions[j] * phi);
            }
            require_positive(denominator.value(), "Wilke denominator")?;
            mixed_mu.add(mole_fractions[i] * (workspace.pure_viscosities[i] / denominator.value()));
            conductivity_arithmetic.add(mole_fractions[i] * workspace.pure_conductivities[i]);
            conductivity_reciprocal.add(mole_fractions[i] / workspace.pure_conductivities[i]);
        }
        let mu = mixed_mu.value();
        require_positive(mu, "mixture viscosity")?;
        require_positive(
            conductivity_arithmetic.value(),
            "arithmetic conductivity mixture",
        )?;
        require_positive(
            conductivity_reciprocal.value(),
            "reciprocal conductivity mixture",
        )?;
        let conductivity =
            0.5 * conductivity_arithmetic.value() + 0.5 / conductivity_reciprocal.value();
        require_positive(conductivity, "mixture conductivity")?;
        poll(&mut cancelled)?;
        Ok(KineticTransportView {
            temperature,
            pressure,
            dynamic_viscosity: DynViscosity::new(mu),
            thermal_conductivity: ThermalConductivityQuantityV1::new(conductivity),
            maximum_binary_diffusivity: maximum_d,
            workspace,
        })
    }
}

/// O(N²) reusable property buffers; successful evaluations do not reallocate.
#[derive(Debug)]
pub struct KineticTransportWorkspaceV1 {
    pure_viscosities: Vec<f64>,
    pure_conductivities: Vec<f64>,
    binary_diffusivities: Vec<f64>,
}

/// Complete common-state property evaluation borrowing reusable storage.
#[derive(Clone, Copy, Debug)]
pub struct KineticTransportView<'a> {
    temperature: Temperature,
    pressure: Pressure,
    dynamic_viscosity: DynViscosity,
    thermal_conductivity: ThermalConductivityQuantityV1,
    maximum_binary_diffusivity: f64,
    workspace: &'a KineticTransportWorkspaceV1,
}
impl KineticTransportView<'_> {
    /// Actual common source temperature.
    pub const fn temperature(&self) -> Temperature {
        self.temperature
    }
    /// Actual common source pressure.
    pub const fn pressure(&self) -> Pressure {
        self.pressure
    }
    /// Wilke mixture dynamic viscosity [Pa s].
    pub const fn dynamic_viscosity(&self) -> DynViscosity {
        self.dynamic_viscosity
    }
    /// Eucken/arithmetic-harmonic mixture conductivity [W/(m K)].
    pub const fn thermal_conductivity(&self) -> ThermalConductivityQuantityV1 {
        self.thermal_conductivity
    }
    /// Canonical row-major symmetric binary diffusivity matrix [m²/s].
    /// Diagonal entries are exactly zero because no interspecies pair exists.
    pub fn binary_diffusivities(&self) -> &[f64] {
        &self.workspace.binary_diffusivities
    }
    /// Largest binary coefficient [m²/s], zero for a pure one-species bank.
    pub const fn maximum_binary_diffusivity(&self) -> f64 {
        self.maximum_binary_diffusivity
    }
    /// Canonical pure-species viscosities [Pa s] at the same temperature.
    pub fn pure_viscosities(&self) -> &[f64] {
        &self.workspace.pure_viscosities
    }
    /// Canonical Eucken pure-species conductivities [W/(m K)].
    pub fn pure_conductivities(&self) -> &[f64] {
        &self.workspace.pure_conductivities
    }
}

fn admit_reduced_temperature(t: f64) -> Result<(), KineticTransportErrorV1> {
    if !t.is_finite() || !(0.3..=100.0).contains(&t) {
        return Err(KineticTransportErrorV1::CollisionDomain {
            reduced_temperature: t,
        });
    }
    Ok(())
}
// Evaluate only the needed Table I integral: the O(N^2) binary diffusion loop
// must not also evaluate the viscosity exponentials and oscillatory correction.
fn collision_integral_11(t: f64) -> Result<f64, KineticTransportErrorV1> {
    admit_reduced_temperature(t)?;
    let ln = fs_math::det::ln(t);
    let omega11 = 1.06036 * fs_math::det::exp(-0.15610 * ln)
        + 0.19300 * fs_math::det::exp(-0.47635 * t)
        + 1.03587 * fs_math::det::exp(-1.52996 * t)
        + 1.76474 * fs_math::det::exp(-3.89411 * t);
    require_positive(omega11, "diffusion collision integral")?;
    Ok(omega11)
}
// The full viscosity fit includes its oscillatory correction, omitted in some
// lower-accuracy textbook transcriptions.
fn collision_integral_22(t: f64) -> Result<f64, KineticTransportErrorV1> {
    admit_reduced_temperature(t)?;
    let ln = fs_math::det::ln(t);
    let omega22 = 1.16145 * fs_math::det::exp(-0.14874 * ln)
        + 0.52487 * fs_math::det::exp(-0.77320 * t)
        + 2.16178 * fs_math::det::exp(-2.43787 * t)
        - 0.0006435
            * fs_math::det::exp(0.14874 * ln)
            * fs_math::det::sin(18.0323 * fs_math::det::exp(-0.76830 * ln) - 7.27371);
    require_positive(omega22, "viscosity collision integral")?;
    Ok(omega22)
}
#[cfg(test)]
fn collision_integrals(t: f64) -> Result<(f64, f64), KineticTransportErrorV1> {
    Ok((collision_integral_11(t)?, collision_integral_22(t)?))
}
fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn require_positive(x: f64, name: &'static str) -> Result<(), KineticTransportErrorV1> {
    if positive(x) {
        Ok(())
    } else {
        Err(KineticTransportErrorV1::Arithmetic(name))
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), KineticTransportErrorV1> {
    if cancelled() {
        Err(KineticTransportErrorV1::Cancelled)
    } else {
        Ok(())
    }
}
fn filled<T: Clone>(n: usize, value: T) -> Result<Vec<T>, KineticTransportErrorV1> {
    let mut v = Vec::new();
    v.try_reserve_exact(n)
        .map_err(|_| KineticTransportErrorV1::Allocation)?;
    v.resize(n, value);
    Ok(v)
}
#[derive(Default)]
struct Compensated {
    value: f64,
    correction: f64,
}
impl Compensated {
    fn add(&mut self, value: f64) {
        let adjusted = value - self.correction;
        let next = self.value + adjusted;
        self.correction = (next - self.value) - adjusted;
        self.value = next;
    }
    fn value(&self) -> f64 {
        self.value
    }
}
fn sum(values: impl IntoIterator<Item = f64>) -> f64 {
    let mut total = Compensated::default();
    for value in values {
        total.add(value);
    }
    total.value()
}

#[cfg(test)]
mod tests;
