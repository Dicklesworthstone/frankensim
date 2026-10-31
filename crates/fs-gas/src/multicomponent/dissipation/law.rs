use super::{DissipationError, positive};
use fs_qty::{DynViscosity, Pressure, Temperature};
use fs_thermochem::equilibrium::composition::IdealGasEquilibriumV1;
pub use fs_thermochem::transport::kinetic::ThermalConductivityQuantityV1;
use fs_thermochem::transport::kinetic::{KineticGasTransportV1, KineticTransportWorkspaceV1};
use fs_thermochem::{MassSpecificEnergyQuantityV1, SpeciesId};

/// Common-state transport coefficients in coherent SI.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TransportCoefficients {
    /// Mixture shear viscosity [Pa s].
    pub dynamic_viscosity: DynViscosity,
    /// Authored bulk viscosity [Pa s]; zero explicitly selects Stokes' hypothesis.
    pub bulk_viscosity: DynViscosity,
    /// Thermal conductivity at zero diffusive species flux [W/(m K)].
    /// With Soret/Dufour coupling, this excludes the reciprocal heat flux.
    pub thermal_conductivity: ThermalConductivityQuantityV1,
}
impl TransportCoefficients {
    pub(super) fn validate(self) -> Result<(), DissipationError> {
        if ![
            self.dynamic_viscosity.value(),
            self.bulk_viscosity.value(),
            self.thermal_conductivity.value(),
        ]
        .iter()
        .all(|x| x.is_finite() && *x >= 0.0)
        {
            return Err(DissipationError::Invalid(
                "nonnegative finite molecular transport coefficients required",
            ));
        }
        Ok(())
    }
    pub(super) fn longitudinal_viscosity(self) -> Result<f64, DissipationError> {
        let value = (4.0 / 3.0) * self.dynamic_viscosity.value() + self.bulk_viscosity.value();
        if !value.is_finite() {
            return Err(DissipationError::Invalid("longitudinal viscosity overflow"));
        }
        Ok(value)
    }
}

/// Source-domain declaration for an explicitly constant transport law.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TransportDomain {
    /// Smallest admitted temperature [K].
    pub minimum_temperature_k: f64,
    /// Largest admitted temperature [K].
    pub maximum_temperature_k: f64,
    /// Smallest admitted absolute pressure [Pa].
    pub minimum_pressure_pa: f64,
    /// Largest admitted absolute pressure [Pa].
    pub maximum_pressure_pa: f64,
}
impl TransportDomain {
    fn valid(self) -> bool {
        [
            self.minimum_temperature_k,
            self.maximum_temperature_k,
            self.minimum_pressure_pa,
            self.maximum_pressure_pa,
        ]
        .iter()
        .all(|x| positive(*x))
            && self.minimum_temperature_k <= self.maximum_temperature_k
            && self.minimum_pressure_pa <= self.maximum_pressure_pa
    }
    fn admits(self, t: Temperature, p: Pressure) -> bool {
        (self.minimum_temperature_k..=self.maximum_temperature_k).contains(&t.value())
            && (self.minimum_pressure_pa..=self.maximum_pressure_pa).contains(&p.value())
    }
}

/// An original-source law evaluated at the actual common face state.
///
/// Implementations retain immutable source/law identity and may mutate only
/// reusable workspaces. They write every `n*n` binary-diffusion entry in canonical
/// source order, with zero diagonal and positive symmetric off-diagonals. The
/// caller checks these physical values before using any flux. This trait does
/// not grant validation authority to an empirical or caller-supplied law.
pub trait TransportLaw {
    /// Exact original caloric/species bank used by this transport law.
    fn source(&self) -> &IdealGasEquilibriumV1;
    /// Evaluate actual T,p and mole composition; all output diffusivities are m²/s.
    fn evaluate(
        &mut self,
        temperature: Temperature,
        pressure: Pressure,
        mole_fractions: &[f64],
        binary_diffusivities: &mut [f64],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<TransportCoefficients, DissipationError>;

    /// Optional source-owned heats of transport `H_i` [J/kg] at the actual
    /// common face state, in canonical species order. Return `true` only after
    /// writing every finite coefficient. These are kinetic cross coefficients,
    /// distinct from the original thermodynamic partial enthalpies.
    ///
    /// The barycentric closure uses the balanced thermal force
    /// `x_i M_i (H_i-sum_j Y_j H_j)/(RT) grad(ln T)` and reciprocal reduced
    /// heat flux `sum_i H_i j_i`. A common shift of every H has no effect.
    /// The accompanying conductivity must use the zero-diffusive-flux
    /// convention. The default explicitly omits this additional closure;
    /// neither gas names nor collision/thermodynamic data infer these values.
    fn heats_of_transport(
        &mut self,
        _temperature: Temperature,
        _pressure: Pressure,
        _mole_fractions: &[f64],
        output_j_kg: &mut [f64],
        _cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<bool, DissipationError> {
        output_j_kg.fill(0.0);
        Ok(false)
    }
}

/// Explicit constant Soret/Dufour cross coefficients on a declared source domain.
///
/// Wraps an original molecular law without estimating missing transport data.
/// The author declares that its conductivity is measured/modelled at zero
/// diffusive species flux and supplies independent heats of transport in J/kg,
/// not partial enthalpies or heat capacities. The shared face owner applies
/// both cross effects together. A state-dependent provider can instead implement
/// [`TransportLaw::heats_of_transport`] directly.
#[derive(Clone, Debug, PartialEq)]
pub struct SoretDufourTransport<L: TransportLaw> {
    law: L,
    model_id: String,
    model_version: u32,
    domain: TransportDomain,
    heats: Vec<(SpeciesId, MassSpecificEnergyQuantityV1)>,
}
impl<L: TransportLaw> SoretDufourTransport<L> {
    /// Bind authored coefficients to the exact underlying canonical species.
    /// Positive, negative and zero heats are admitted; no universal sign of
    /// thermal separation is assumed. The complete coupled Onsager law remains
    /// positive semidefinite for positive pair D and nonnegative conductivity.
    pub fn new(
        law: L,
        model_id: impl Into<String>,
        model_version: u32,
        domain: TransportDomain,
        mut heats: Vec<(SpeciesId, MassSpecificEnergyQuantityV1)>,
    ) -> Result<Self, DissipationError> {
        let model_id = model_id.into();
        heats.sort_by(|a, b| a.0.cmp(&b.0));
        if model_id.is_empty()
            || model_id.len() > 128
            || model_version == 0
            || !domain.valid()
            || heats.len() != law.source().models().len()
            || heats
                .iter()
                .zip(law.source().models())
                .any(|((id, h), model)| id != model.species() || !h.value().is_finite())
        {
            return Err(DissipationError::Invalid(
                "Soret/Dufour identity, species, heats or domain",
            ));
        }
        Ok(Self {
            law,
            model_id,
            model_version,
            domain,
            heats,
        })
    }
    /// Underlying original-source molecular law.
    pub const fn law(&self) -> &L {
        &self.law
    }
    /// Authored thermal cross-coefficient source identity.
    pub fn model_id(&self) -> &str {
        &self.model_id
    }
    /// Authored thermal cross-coefficient version.
    pub const fn model_version(&self) -> u32 {
        self.model_version
    }
    /// Validity of this constant thermal cross-coefficient approximation.
    pub const fn domain(&self) -> TransportDomain {
        self.domain
    }
    /// Original coefficients, with typed J/kg values and canonical species IDs.
    pub fn heats(&self) -> &[(SpeciesId, MassSpecificEnergyQuantityV1)] {
        &self.heats
    }
}
impl<L: TransportLaw> TransportLaw for SoretDufourTransport<L> {
    fn source(&self) -> &IdealGasEquilibriumV1 {
        self.law.source()
    }
    fn evaluate(
        &mut self,
        t: Temperature,
        p: Pressure,
        x: &[f64],
        binary: &mut [f64],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<TransportCoefficients, DissipationError> {
        if !self.domain.admits(t, p) {
            return Err(DissipationError::Invalid(
                "Soret/Dufour state outside source domain",
            ));
        }
        self.law.evaluate(t, p, x, binary, cancelled)
    }
    fn heats_of_transport(
        &mut self,
        t: Temperature,
        p: Pressure,
        x: &[f64],
        out: &mut [f64],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<bool, DissipationError> {
        if cancelled() {
            return Err(DissipationError::Cancelled);
        }
        if !self.domain.admits(t, p) || x.len() != self.heats.len() || out.len() != self.heats.len()
        {
            return Err(DissipationError::Invalid(
                "Soret/Dufour state or output axis",
            ));
        }
        for (output, (_, heat)) in out.iter_mut().zip(&self.heats) {
            *output = heat.value();
        }
        Ok(true)
    }
}

/// Caller-authored constant coefficients on an explicit physical domain.
///
/// Useful for analytic transport fixtures and genuinely constant-coefficient
/// approximations. No gas labels or thermodynamic values infer diffusivities.
#[derive(Clone, Debug, PartialEq)]
pub struct ConstantTransport {
    source: IdealGasEquilibriumV1,
    model_id: String,
    model_version: u32,
    domain: TransportDomain,
    coefficients: TransportCoefficients,
    binary_diffusivities: Vec<f64>,
}
impl ConstantTransport {
    /// Admit one complete canonical pair matrix and a declared constant law.
    pub fn new(
        source: IdealGasEquilibriumV1,
        model_id: impl Into<String>,
        model_version: u32,
        domain: TransportDomain,
        coefficients: TransportCoefficients,
        binary_diffusivities: Vec<f64>,
    ) -> Result<Self, DissipationError> {
        let model_id = model_id.into();
        if model_id.is_empty() || model_id.len() > 128 || model_version == 0 || !domain.valid() {
            return Err(DissipationError::Invalid(
                "constant molecular law identity or domain",
            ));
        }
        coefficients.validate()?;
        validate_binary(&binary_diffusivities, source.models().len())?;
        Ok(Self {
            source,
            model_id,
            model_version,
            domain,
            coefficients,
            binary_diffusivities,
        })
    }
    /// Authored law identity.
    pub fn model_id(&self) -> &str {
        &self.model_id
    }
    /// Authored law version.
    pub const fn model_version(&self) -> u32 {
        self.model_version
    }
    /// Explicit validity of the constant approximation.
    pub const fn domain(&self) -> TransportDomain {
        self.domain
    }
    /// Original coherent-SI constant coefficients.
    pub const fn coefficients(&self) -> TransportCoefficients {
        self.coefficients
    }
}
impl TransportLaw for ConstantTransport {
    fn source(&self) -> &IdealGasEquilibriumV1 {
        &self.source
    }
    fn evaluate(
        &mut self,
        t: Temperature,
        p: Pressure,
        x: &[f64],
        out: &mut [f64],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<TransportCoefficients, DissipationError> {
        if cancelled() {
            return Err(DissipationError::Cancelled);
        }
        if !self.domain.admits(t, p)
            || x.len() != self.source.models().len()
            || out.len() != self.binary_diffusivities.len()
        {
            return Err(DissipationError::Invalid(
                "constant molecular law state or output axis",
            ));
        }
        out.copy_from_slice(&self.binary_diffusivities);
        Ok(self.coefficients)
    }
}

/// Allocation-free adapter to the original dilute-gas collision-property owner.
pub struct KineticTransport {
    model: KineticGasTransportV1,
    workspace: KineticTransportWorkspaceV1,
    bulk_viscosity: DynViscosity,
}
impl KineticTransport {
    /// Bind a collision model and explicit nonnegative bulk viscosity.
    /// `0 Pa s` selects the Stokes approximation rather than an inferred value.
    pub fn new(
        model: KineticGasTransportV1,
        bulk_viscosity: DynViscosity,
    ) -> Result<Self, DissipationError> {
        if !bulk_viscosity.value().is_finite() || bulk_viscosity.value() < 0.0 {
            return Err(DissipationError::Invalid("bulk viscosity"));
        }
        let workspace = model.workspace()?;
        Ok(Self {
            model,
            workspace,
            bulk_viscosity,
        })
    }
    /// Original source-owned collision model and versioned data.
    pub const fn model(&self) -> &KineticGasTransportV1 {
        &self.model
    }
    /// Explicit bulk-viscosity declaration [Pa s].
    pub const fn bulk_viscosity(&self) -> DynViscosity {
        self.bulk_viscosity
    }
}
impl TransportLaw for KineticTransport {
    fn source(&self) -> &IdealGasEquilibriumV1 {
        self.model.thermodynamics()
    }
    fn evaluate(
        &mut self,
        t: Temperature,
        p: Pressure,
        x: &[f64],
        out: &mut [f64],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<TransportCoefficients, DissipationError> {
        if out.len() != self.model.species_count() * self.model.species_count() {
            return Err(DissipationError::Invalid("kinetic transport output axis"));
        }
        let view = self
            .model
            .evaluate(t, p, x, &mut self.workspace, || cancelled())?;
        out.copy_from_slice(view.binary_diffusivities());
        Ok(TransportCoefficients {
            dynamic_viscosity: view.dynamic_viscosity(),
            bulk_viscosity: self.bulk_viscosity,
            thermal_conductivity: view.thermal_conductivity(),
        })
    }
}

pub(crate) fn validate_binary(d: &[f64], n: usize) -> Result<f64, DissipationError> {
    if d.len() != n * n {
        return Err(DissipationError::Invalid("binary diffusivity dimensions"));
    }
    let mut maximum = 0.0_f64;
    for i in 0..n {
        if d[i * n + i] != 0.0 {
            return Err(DissipationError::Invalid("binary diffusion diagonal"));
        }
        for j in i + 1..n {
            if !positive(d[i * n + j]) || d[i * n + j] != d[j * n + i] {
                return Err(DissipationError::Invalid(
                    "positive symmetric binary diffusivities required",
                ));
            }
            maximum = maximum.max(d[i * n + j]);
        }
    }
    Ok(maximum)
}
