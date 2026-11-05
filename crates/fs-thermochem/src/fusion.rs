//! Bounded Maxwellian nuclear-fusion laws, independent of chemical reactions.
//!
//! Bosch–Hale (1992), equation (14), table VII, supplies separate D–T,
//! D–D neutron, D–D proton and D–He3 reactivities. The input is one common
//! Maxwellian ion temperature; rates are in m³/s per particle pair. The
//! identical-reactant factor belongs to event counting, not to reactivity.
//! Nuclear reactions conserve baryon number and electric charge. They change
//! chemical elements and rest mass, so chemical elemental/mass-conservation
//! certificates must not be applied to them.
//!
//! This module owns no burn integrator, confinement, beam distribution,
//! screening, neutron transport, charged-product slowing or ignition claim.
//! Reaction Q values use nominal CODATA 2022 bare-nucleus masses. Product
//! energy partition is explicitly leading nonrelativistic two-body kinematics.

use fs_qty::Temperature;

pub use crate::electron_collision::{BOLTZMANN_J_PER_K, ELECTRON_MASS_KG, ELECTRONVOLT_J};

/// Explicit binding of a D–T plasma to its original gas thermodynamics.
pub mod plasma;
mod reactivity;

#[cfg(test)]
mod tests;

/// Primary source of the reactivity equation, coefficients and fit domains.
pub const BOSCH_HALE_SOURCE_URL_V1: &str = "https://doi.org/10.1088/0029-5515/32/4/I07";
/// Primary source of the nominal bare-nucleus masses, CODATA 2022.
pub const NUCLEAR_MASS_SOURCE_URL_V1: &str =
    "https://physics.nist.gov/cuu/Constants/Table/allascii.txt";
/// Exact SI speed of light [m/s].
pub const SPEED_OF_LIGHT_M_S: f64 = 299_792_458.0;
/// Exact SI unit conversion, up to binary64 rounding: kelvin per keV.
pub const KELVIN_PER_KEV: f64 = 1000.0 * ELECTRONVOLT_J / BOLTZMANN_J_PER_K;

/// A bare nucleus or free neutron, separate from chemical species identities.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum NuclearSpeciesV1 {
    /// Hydrogen-1 nucleus, charge +e.
    Proton,
    /// Hydrogen-2 nucleus, charge +e.
    Deuterium,
    /// Hydrogen-3 nucleus, charge +e.
    Tritium,
    /// Helium-3 nucleus, charge +2e.
    Helium3,
    /// Helium-4 nucleus (alpha particle), charge +2e.
    Helium4,
    /// Free neutron, electrically neutral.
    Neutron,
}

impl NuclearSpeciesV1 {
    /// Fixed supported species order; no isotope inference from chemical names.
    pub const ALL: [Self; 6] = [
        Self::Proton,
        Self::Deuterium,
        Self::Tritium,
        Self::Helium3,
        Self::Helium4,
        Self::Neutron,
    ];

    /// Integer baryon number (protons plus neutrons).
    #[must_use]
    pub const fn baryon_number(self) -> u8 {
        match self {
            Self::Proton | Self::Neutron => 1,
            Self::Deuterium => 2,
            Self::Tritium | Self::Helium3 => 3,
            Self::Helium4 => 4,
        }
    }

    /// Positive nuclear charge in units of the elementary charge.
    #[must_use]
    pub const fn charge_number(self) -> u8 {
        match self {
            Self::Proton | Self::Deuterium | Self::Tritium => 1,
            Self::Helium3 | Self::Helium4 => 2,
            Self::Neutron => 0,
        }
    }

    /// Nominal CODATA 2022 bare-nucleus/free-neutron rest mass [kg].
    /// These are not atomic masses; bound electrons are absent. Source
    /// uncertainty is not propagated by the fixed-coefficient model.
    #[must_use]
    pub const fn mass_kg(self) -> f64 {
        match self {
            Self::Proton => 1.672_621_925_95e-27,
            Self::Deuterium => 3.343_583_776_8e-27,
            Self::Tritium => 5.007_356_751_2e-27,
            Self::Helium3 => 5.006_412_786_2e-27,
            Self::Helium4 => 6.644_657_345_0e-27,
            Self::Neutron => 1.674_927_500_56e-27,
        }
    }
}

/// Supported, separately evaluated two-body fusion channels.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum FusionReactionV1 {
    /// D + T → He4 + n.
    DeuteriumTritium,
    /// D + D → He3 + n.
    DeuteriumDeuteriumNeutron,
    /// D + D → T + p.
    DeuteriumDeuteriumProton,
    /// D + He3 → He4 + p.
    DeuteriumHelium3,
}

impl FusionReactionV1 {
    /// Fixed channel order. The two D–D branches are independent rate laws.
    pub const ALL: [Self; 4] = [
        Self::DeuteriumTritium,
        Self::DeuteriumDeuteriumNeutron,
        Self::DeuteriumDeuteriumProton,
        Self::DeuteriumHelium3,
    ];

    /// Reactant slots; repeated D slots explicitly represent two consumed ions.
    #[must_use]
    pub const fn reactants(self) -> [NuclearSpeciesV1; 2] {
        use NuclearSpeciesV1::{Deuterium, Helium3, Tritium};
        match self {
            Self::DeuteriumTritium => [Deuterium, Tritium],
            Self::DeuteriumDeuteriumNeutron | Self::DeuteriumDeuteriumProton => {
                [Deuterium, Deuterium]
            }
            Self::DeuteriumHelium3 => [Deuterium, Helium3],
        }
    }

    /// Product slots in the same order as the product energy partition.
    #[must_use]
    pub const fn products(self) -> [NuclearSpeciesV1; 2] {
        use NuclearSpeciesV1::{Helium3, Helium4, Neutron, Proton, Tritium};
        match self {
            Self::DeuteriumTritium => [Helium4, Neutron],
            Self::DeuteriumDeuteriumNeutron => [Helium3, Neutron],
            Self::DeuteriumDeuteriumProton => [Tritium, Proton],
            Self::DeuteriumHelium3 => [Helium4, Proton],
        }
    }

    /// Net nuclei created per event; D depletion in each D–D branch is −2.
    #[must_use]
    pub fn stoichiometric_coefficient(self, species: NuclearSpeciesV1) -> i8 {
        let reactants = self.reactants();
        let products = self.products();
        i8::from(products[0] == species) + i8::from(products[1] == species)
            - i8::from(reactants[0] == species)
            - i8::from(reactants[1] == species)
    }

    /// Pair-counting multiplier for `R = factor * n_a * n_b * <sigma v>`.
    /// Number densities are particles/m³. The D–D factor is 1/2, and its
    /// deuterium sink is twice this event rate for each branch.
    #[must_use]
    pub const fn reactant_pair_factor(self) -> f64 {
        match self {
            Self::DeuteriumDeuteriumNeutron | Self::DeuteriumDeuteriumProton => 0.5,
            Self::DeuteriumTritium | Self::DeuteriumHelium3 => 1.0,
        }
    }

    /// Published inclusive ion-temperature interval [keV], table VII.
    #[must_use]
    pub const fn temperature_range_kev(self) -> [f64; 2] {
        match self {
            Self::DeuteriumHelium3 => [0.5, 190.0],
            _ => [0.2, 100.0],
        }
    }

    /// Published inclusive ion-temperature interval converted to kelvin.
    #[must_use]
    pub fn temperature_range(self) -> [Temperature; 2] {
        let [min, max] = self.temperature_range_kev();
        [
            Temperature::new(min * KELVIN_PER_KEV),
            Temperature::new(max * KELVIN_PER_KEV),
        ]
    }

    /// Rest-mass defect converted to nuclear energy release [J/event].
    /// This energy is not automatically deposited in any plasma population.
    #[must_use]
    pub fn q_j(self) -> f64 {
        let [a, b] = self.reactants();
        let [c, d] = self.products();
        ((a.mass_kg() + b.mass_kg()) - (c.mass_kg() + d.mass_kg()))
            * SPEED_OF_LIGHT_M_S
            * SPEED_OF_LIGHT_M_S
    }

    /// Leading nonrelativistic share of relative kinetic energy for each
    /// product, in `products()` order. Equal and opposite product momenta
    /// give the inverse-mass split. The share of center-of-mass translation
    /// is the complementary fraction, neglecting mass-defect boost corrections.
    #[must_use]
    pub fn product_relative_energy_fractions(self) -> [f64; 2] {
        let [a, b] = self.products();
        let first = b.mass_kg() / (a.mass_kg() + b.mass_kg());
        [first, 1.0 - first]
    }

    /// Cold-reactant Q carried by charged products [J/event], using the
    /// declared nonrelativistic partition. Slowing/deposition is consumer-owned.
    #[must_use]
    pub fn charged_q_j(self) -> f64 {
        self.q_j() - self.neutron_q_j()
    }

    /// Cold-reactant Q carried by a neutron [J/event]. Neutron-free channels
    /// return zero. Transport and deposited/escaped energy are consumer-owned.
    #[must_use]
    pub fn neutron_q_j(self) -> f64 {
        if self.products()[1] == NuclearSpeciesV1::Neutron {
            self.q_j() * self.product_relative_energy_fractions()[1]
        } else {
            0.0
        }
    }

    /// Evaluate the bounded Bosch–Hale Maxwellian fit and analytic derivative.
    /// Temperature is in kelvin; reactivity is m³/s, and its derivative m³/(s K).
    /// Outside the published fit interval evaluation fails, without clipping.
    pub fn reactivity(
        self,
        ion_temperature: Temperature,
    ) -> Result<FusionReactivityV1, FusionErrorV1> {
        reactivity::evaluate(self, ion_temperature)
    }

    /// Evaluate the fitted rate and mean kinetic energies of one reacting pair.
    ///
    /// Both reactants share the same Maxwellian temperature and have no relative
    /// bulk drift. Reaction weighting changes their mean relative energy to
    /// `kT [3/2 + d ln(rate)/d ln(T)]`; their independent center-of-mass
    /// translation contributes another `3kT/2`. Leading nonrelativistic two-body
    /// kinematics partitions Q plus relative energy by inverse product mass and
    /// center-of-mass energy by product mass, neglecting mass-defect boost
    /// corrections. The product sum is Q plus the removed reacting-ion energy.
    ///
    /// These are energies per event: the identical-pair factor remains in event
    /// counting. No deposition split, charged-product slowing or neutron escape
    /// closure is supplied. The rate-fit residual is not an error certificate
    /// for its derivative or these energy moments.
    pub fn thermal_event(
        self,
        ion_temperature: Temperature,
    ) -> Result<FusionThermalEventV1, FusionErrorV1> {
        let reactivity = self.reactivity(ion_temperature)?;
        let relative = reactivity.mean_reacting_relative_energy_j();
        let translation = 1.5 * BOLTZMANN_J_PER_K * ion_temperature.value();
        let removed = relative + translation;
        let fraction = self.product_relative_energy_fractions()[0];
        let q = self.q_j();
        let first = fraction * (q + relative) + (1.0 - fraction) * translation;
        let second = q + removed - first;
        if [relative, translation, removed, first, second]
            .iter()
            .any(|&e| !e.is_finite() || e <= 0.0)
        {
            return Err(FusionErrorV1::Arithmetic);
        }
        Ok(FusionThermalEventV1 {
            reactivity,
            reacting_ion_energy_j: removed,
            product_kinetic_energies_j: [first, second],
        })
    }
}

/// Refused temperature domain or nonrepresentable evaluation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum FusionErrorV1 {
    /// Temperature must be finite and strictly positive.
    InvalidTemperature,
    /// The source fit does not cover this finite positive temperature.
    OutsideTemperatureDomain {
        /// Requested reaction channel.
        reaction: FusionReactionV1,
        /// Requested temperature in kelvin.
        temperature_k: f64,
        /// Inclusive source interval in kelvin.
        range_k: [f64; 2],
    },
    /// A required rate or moment could not be represented as finite binary64.
    Arithmetic,
}

impl core::fmt::Display for FusionErrorV1 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "nuclear fusion: {self:?}")
    }
}
impl std::error::Error for FusionErrorV1 {}

/// One allocation-free Maxwellian rate evaluation, with its analytic slope.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FusionReactivityV1 {
    ion_temperature: Temperature,
    value_m3_s: f64,
    derivative_m3_s_k: f64,
    logarithmic_temperature_derivative: f64,
}

impl FusionReactivityV1 {
    /// Common Maxwellian ion temperature used for this evaluation.
    #[must_use]
    pub const fn temperature(self) -> Temperature {
        self.ion_temperature
    }

    /// Particle-pair reactivity `<sigma v>` in m³/s, before pair counting.
    #[must_use]
    pub const fn value_m3_s(self) -> f64 {
        self.value_m3_s
    }

    /// Analytic `d<sigma v>/dT` in m³/(s K).
    #[must_use]
    pub const fn derivative_m3_s_k(self) -> f64 {
        self.derivative_m3_s_k
    }

    /// Dimensionless `d ln<sigma v>/d ln T` of the original fitted law.
    #[must_use]
    pub const fn logarithmic_temperature_derivative(self) -> f64 {
        self.logarithmic_temperature_derivative
    }

    /// Reaction-weighted mean incoming relative kinetic energy [J/event].
    /// For a nonrelativistic Maxwellian integral, differentiating
    /// `T^(-3/2) integral(sigma(E) E exp(-E/kT) dE)` gives
    /// `<E_rel> = kT [3/2 + d ln<sigma v>/d ln T]`.
    /// This applies the identity to the fitted rate; its derivative's physical
    /// error is not bounded by the published rate-fit residual.
    #[must_use]
    pub fn mean_reacting_relative_energy_j(self) -> f64 {
        BOLTZMANN_J_PER_K
            * self.ion_temperature.value()
            * (1.5 + self.logarithmic_temperature_derivative)
    }
}

/// One allocation-free thermal event from the same original fitted reactivity.
/// Product kinetic energies include their Q and thermal shares, before any
/// transport or thermalization. Their order is [`FusionReactionV1::products`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FusionThermalEventV1 {
    reactivity: FusionReactivityV1,
    reacting_ion_energy_j: f64,
    product_kinetic_energies_j: [f64; 2],
}

impl FusionThermalEventV1 {
    /// Original Maxwellian reactivity and analytic slope at the requested T.
    #[must_use]
    pub const fn reactivity(self) -> FusionReactivityV1 {
        self.reactivity
    }

    /// Mean total kinetic energy removed with the two reacting ions [J/event].
    /// Includes relative motion and center-of-mass translation; it is not the
    /// unweighted mean kinetic energy of two randomly selected thermal ions.
    #[must_use]
    pub const fn reacting_ion_energy_j(self) -> f64 {
        self.reacting_ion_energy_j
    }

    /// Mean outgoing kinetic energies [J/event], in the channel's product order.
    #[must_use]
    pub const fn product_kinetic_energies_j(self) -> [f64; 2] {
        self.product_kinetic_energies_j
    }
}
