//! Ground-state single photoionization, A + photon -> I + e.
//!
//! A source table describes one ionization event per absorbed photon, ending
//! in the same ground-state ion as the retained atomic channel. Total optical
//! extinction, excitation and multielectron yields are not interchangeable
//! with this cross section. The first sample declares the physical onset;
//! it must not merely be the lowest available sample of a truncated source.
//!
//! For discrete photon number fluxes Phi_g, the local rate per neutral is
//! sum_g sigma(E_g)*Phi_g. The absorbed energy is E_g per event, the original
//! source chemical cost is chi/NA, and their difference is photoelectron
//! kinetic energy in the stationary-heavy-target approximation. Ion recoil,
//! secondary ionization, stimulated recombination and transport are not supplied.
//! A thermal consumer must explicitly assume photoelectron thermalization.
//!
//! Source tables retain their interpolation and finite energy interval. The
//! optional H I, He I and He II presets directly evaluate the published
//! Verner et al. fit (1996), <https://arxiv.org/abs/astro-ph/9601009>.
//! He I support stops before excited residual-ion channels can open; the
//! total fit cannot specify their individual energy or population yields.
//! No high-energy continuation or spectral quadrature is inferred.

use super::AtomicIonizationV1;
use crate::electron_collision::{
    AVOGADRO_PER_MOL, ELECTRON_MASS_KG, ELECTRONVOLT_J, MAX_CROSS_SECTION_POINTS_V1,
};
use fs_matdb::Provenance;
use fs_math::det::{exp, ln, sqrt};
use fs_qty::{Area, Energy};

/// Maximum discrete photon groups in one local evaluation.
pub const MAX_PHOTON_GROUPS_V1: usize = 4096;

/// One original photon-ionization sample, with explicit SI quantities.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhotoionizationCrossSectionPointV1 {
    /// Photon energy [J], not frequency or electronvolts.
    pub photon_energy: Energy,
    /// Single-ionization event cross section [m²].
    pub cross_section: Area,
}

/// Monochromatic group whose weight is already integrated over its bandwidth.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhotonNumberFluxV1 {
    /// Representative/discrete photon energy [J].
    pub photon_energy: Energy,
    /// Photon number flux [m⁻² s⁻¹], not energy flux or spectral density.
    ///
    /// For a beam this is its ordinary directional number flux. For a local
    /// angular distribution it is the angle-integrated number intensity,
    /// equivalently c times photon number density, not the net vector flux.
    /// No bandwidth or 4*pi factor is inserted by the evaluator.
    pub number_flux_per_m2_s: f64,
}

/// Invalid source/field, missing source coverage, arithmetic range or cancellation.
#[derive(Clone, Debug, PartialEq)]
pub enum PhotoionizationErrorV1 {
    /// Invalid physical data, provenance or bounded source size.
    Invalid(&'static str),
    /// Photon energy exceeds the source's declared upper limit.
    OutsideEnergyDomain,
    /// A nonnegative rate or energy moment is outside binary64 range.
    Arithmetic(&'static str),
    /// Cancellation observed before a complete result was available.
    Cancelled,
}
impl core::fmt::Display for PhotoionizationErrorV1 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "photoionization: {self:?}")
    }
}
impl std::error::Error for PhotoionizationErrorV1 {}

#[derive(Clone, Debug, PartialEq)]
enum PhotoionizationLaw {
    Table(std::sync::Arc<[PhotoionizationCrossSectionPointV1]>),
    Verner {
        fit: VernerFit,
        maximum_energy: Energy,
    },
}

/// Literal parameters from Verner et al. (1996), Eq.(1), Table 1. The energy
/// and area scales retain the published eV and Mb convention until evaluation.
#[derive(Clone, Copy, Debug, PartialEq)]
struct VernerFit {
    e0: f64,
    sigma0: f64,
    ya: f64,
    p: f64,
    yw: f64,
    y0: f64,
    y1: f64,
}
const VERNER_H_I: VernerFit = VernerFit {
    e0: 0.4298,
    sigma0: 54_750.0,
    ya: 32.88,
    p: 2.963,
    yw: 0.0,
    y0: 0.0,
    y1: 0.0,
};
const VERNER_HE_I: VernerFit = VernerFit {
    e0: 13.61,
    sigma0: 949.2,
    ya: 1.469,
    p: 3.188,
    yw: 2.039,
    y0: 0.4434,
    y1: 2.136,
};
const VERNER_HE_II: VernerFit = VernerFit {
    e0: 1.720,
    sigma0: 13_690.0,
    ya: 32.88,
    p: 2.963,
    yw: 0.0,
    y0: 0.0,
    y1: 0.0,
};

/// Immutable photon source on the original validated atomic species bank.
#[derive(Clone, Debug, PartialEq)]
pub struct AtomicPhotoionizationV1 {
    channel: AtomicIonizationV1,
    provenance: Provenance,
    law: PhotoionizationLaw,
    threshold: Energy,
}
impl AtomicPhotoionizationV1 {
    /// Bind a complete piecewise-linear event table to the original atom/ion.
    ///
    /// The first sample is the declared onset and may have NONZERO cross
    /// section: a photoionization edge need not vanish like electron impact.
    /// Onset must be at least the original source gap chi/NA. A discrepancy
    /// below that gap of at most 64*epsilon relative is admitted as roundoff;
    /// only the effective onset is rounded upward to chi/NA. Original samples
    /// are retained unchanged. Above the final sample every query refuses.
    pub fn new(
        channel: AtomicIonizationV1,
        provenance: Provenance,
        points: Vec<PhotoionizationCrossSectionPointV1>,
    ) -> Result<Self, PhotoionizationErrorV1> {
        if provenance.source.trim().is_empty() || provenance.license.trim().is_empty() {
            return Err(PhotoionizationErrorV1::Invalid(
                "source citation and usage terms",
            ));
        }
        if !(2..=MAX_CROSS_SECTION_POINTS_V1).contains(&points.len())
            || points.iter().any(|p| {
                !p.photon_energy.value().is_finite()
                    || p.photon_energy.value() <= 0.0
                    || !nonnegative(p.cross_section.value())
            })
            || points
                .windows(2)
                .any(|p| p[0].photon_energy.value() >= p[1].photon_energy.value())
        {
            return Err(PhotoionizationErrorV1::Invalid(
                "ordered SI photon cross-section table",
            ));
        }
        let gap = channel.ionization_energy_j_mol() / AVOGADRO_PER_MOL;
        let first = points[0].photon_energy.value();
        let last = points[points.len() - 1].photon_energy.value();
        if !gap.is_finite()
            || gap <= 0.0
            || gap - first > 64.0 * f64::EPSILON * gap
            || last < gap
        {
            return Err(PhotoionizationErrorV1::Invalid(
                "photon threshold below source ionization gap",
            ));
        }
        Ok(Self {
            channel,
            provenance,
            law: PhotoionizationLaw::Table(points.into()),
            threshold: Energy::new(first.max(gap)),
        })
    }

    /// Published Verner et al. (1996) H I ground-state fit, Table 1 / Eq.(1).
    ///
    /// Original support is 13.60..=50000 eV. The rounded 13.60 eV onset is
    /// retained; it does not replace the original source gap in energy accounts.
    /// The gap must lie within the table's 0.005 eV rounding allowance below
    /// that onset. Source charges, the H elemental row, electronic weights
    /// (2,1), and proton/atom masses are checked without reading species names.
    /// The 1e-3 relative mass allowance is an identity-admission policy, not
    /// fit uncertainty or an isotope correction; D/T masses do not qualify.
    ///
    /// Parameters: E0=0.4298 eV, sigma0=54750 Mb, ya=32.88, P=2.963,
    /// yw=y0=y1=0; 1 Mb=1e-22 m². Sources:
    /// <https://www.pa.uky.edu/~verner/dima/photo/photo.dat>,
    /// <https://arxiv.org/abs/astro-ph/9601009>.
    /// This is a nonrelativistic fitted cross section, not an exact hydrogen
    /// wavefunction calculation or a certified physical-error bound.
    pub fn verner_hydrogen(
        channel: AtomicIonizationV1,
    ) -> Result<Self, PhotoionizationErrorV1> {
        let bank = channel.thermodynamics();
        let matrix = bank.elemental_matrix();
        let h = matrix.elements().iter().position(|element| element.as_str() == "H")
            .ok_or(PhotoionizationErrorV1::Invalid("Verner H I requires the H elemental row"))?;
        let indices = [channel.neutral_index(), channel.ion_index(), channel.electron_index()];
        let weights = channel.statistical_weights();
        if matrix.get(h, indices[0]) != Some(1)
            || matrix.get(h, indices[1]) != Some(1)
            || matrix.get(h, indices[2]) != Some(0)
            || indices.map(|i| bank.charge_vector().get(i)) != [Some(0), Some(1), Some(-1)]
            || weights.neutral != 2
            || weights.ion != 1
        {
            return Err(PhotoionizationErrorV1::Invalid("Verner H I ground-state identities"));
        }
        // CODATA 2022; the electron mass is the shared collision constant.
        // https://physics.nist.gov/cuu/pdf/wall_2022.pdf
        const PROTON_MASS_KG: f64 = 1.672_621_925_95e-27;
        for (i, expected) in indices.into_iter().zip([
            PROTON_MASS_KG + ELECTRON_MASS_KG, PROTON_MASS_KG, ELECTRON_MASS_KG,
        ]) {
            let mass = bank.models()[i].molar_mass().value() / AVOGADRO_PER_MOL;
            if (mass / expected - 1.0).abs() > 1e-3 {
                return Err(PhotoionizationErrorV1::Invalid("Verner H I source masses"));
            }
        }
        let gap = channel.ionization_energy_j_mol() / AVOGADRO_PER_MOL;
        let edge = 13.60 * ELECTRONVOLT_J;
        if gap < 13.595 * ELECTRONVOLT_J || gap - edge > 64.0 * f64::EPSILON * edge {
            return Err(PhotoionizationErrorV1::Invalid("Verner H I source ionization gap"));
        }
        Ok(Self {
            channel,
            provenance: Provenance {
                source: "Verner, Ferland, Korista & Yakovlev (1996), ApJ 465, 487, Eq.(1), Table 1 H I; https://www.pa.uky.edu/~verner/dima/photo/photo.dat; https://doi.org/10.1086/177435".into(),
                license: "Published numerical fit parameters; independently implemented equation".into(),
                artifact: None,
            },
            law: PhotoionizationLaw::Verner {
                fit: VERNER_H_I,
                maximum_energy: Energy::new(50_000.0 * ELECTRONVOLT_J),
            },
            threshold: Energy::new(edge.max(gap)),
        })
    }

    /// Verner et al. (1996), Eq.(1) / Table 1 He I, on its ground-state channel.
    /// The original rounded 24.59 eV onset and all seven fit parameters remain
    /// unchanged. Structural He, charges (0,+1,-1), weights (1,2), He-4 masses
    /// and a source gap in 24.585..=24.59 eV are required; names do not select it.
    ///
    /// The published total fit extends to 50 keV, but this single-ground-state
    /// owner admits only up to 65.39 eV. He+ 2p_1/2 opens at about 65.40042 eV
    /// from neutral He; above it the total cross section cannot assign the
    /// residual ion's excitation or the electron's heat. 65.39 eV stays below
    /// that opening even at the admitted source gap's lower rounding limit.
    /// The first He+ excitation is 329179.2940 cm^-1 (about 40.81303 eV):
    /// <https://physics.nist.gov/PhysRefData/Handbook/Tables/heliumtable6.htm>.
    /// Unresolved higher-energy tails are not replaced with zero or a guessed
    /// branching ratio. These are smoothed fit values, not resolved resonances.
    pub fn verner_helium(channel: AtomicIonizationV1) -> Result<Self, PhotoionizationErrorV1> {
        Self::verner_helium_stage(channel, false)
    }

    /// Verner et al. (1996), Eq.(1) / Table 1 He II, with published support
    /// 54.42..=50000 eV. Requires structural He, charges (+1,+2,-1), ground
    /// weights (2,1), He-4 masses and a source gap in 54.415..=54.42 eV.
    /// The stripped upper ion has no unresolved excited electronic channel.
    /// The rounded fit onset never replaces the original gap in energy accounts.
    pub fn verner_helium_ion(
        channel: AtomicIonizationV1,
    ) -> Result<Self, PhotoionizationErrorV1> {
        Self::verner_helium_stage(channel, true)
    }

    fn verner_helium_stage(
        channel: AtomicIonizationV1,
        ionized: bool,
    ) -> Result<Self, PhotoionizationErrorV1> {
        let bank = channel.thermodynamics();
        let matrix = bank.elemental_matrix();
        let he = matrix.elements().iter().position(|element| element.as_str() == "He")
            .ok_or(PhotoionizationErrorV1::Invalid("Verner helium requires the He elemental row"))?;
        let indices = [channel.lower_index(), channel.upper_index(), channel.electron_index()];
        let weights = channel.statistical_weights();
        let (charge, expected_weights, edge_ev, fit, maximum_ev, label) = if ionized {
            (1, [2, 1], 54.42, VERNER_HE_II, 50_000.0, "He II")
        } else {
            (0, [1, 2], 24.59, VERNER_HE_I, 65.39, "He I")
        };
        if matrix.get(he, indices[0]) != Some(1)
            || matrix.get(he, indices[1]) != Some(1)
            || matrix.get(he, indices[2]) != Some(0)
            || indices.map(|i| bank.charge_vector().get(i))
                != [Some(charge), Some(charge + 1), Some(-1)]
            || [weights.neutral, weights.ion] != expected_weights
        {
            return Err(PhotoionizationErrorV1::Invalid("Verner helium ground-state identities"));
        }
        // CODATA 2022 alpha mass; the 1e-3 relative window is the same identity
        // policy as H I, not a fit-error bound or a He-3 isotope correction.
        // https://physics.nist.gov/cuu/Constants/Table/allascii.txt
        const ALPHA_MASS_KG: f64 = 6.644_657_345_0e-27;
        let upper_mass = ALPHA_MASS_KG + if ionized { 0.0 } else { ELECTRON_MASS_KG };
        for (index, expected) in indices.into_iter().zip([
            upper_mass + ELECTRON_MASS_KG, upper_mass, ELECTRON_MASS_KG,
        ]) {
            let mass = bank.models()[index].molar_mass().value() / AVOGADRO_PER_MOL;
            if (mass / expected - 1.0).abs() > 1e-3 {
                return Err(PhotoionizationErrorV1::Invalid("Verner helium source masses"));
            }
        }
        let gap = channel.ionization_energy_j_mol() / AVOGADRO_PER_MOL;
        let edge = edge_ev * ELECTRONVOLT_J;
        if gap < (edge_ev - 0.005) * ELECTRONVOLT_J
            || gap - edge > 64.0 * f64::EPSILON * edge
        {
            return Err(PhotoionizationErrorV1::Invalid("Verner helium source ionization gap"));
        }
        Ok(Self {
            channel,
            provenance: Provenance {
                source: format!("Verner, Ferland, Korista & Yakovlev (1996), ApJ 465, 487, Eq.(1), Table 1 {label}; https://www.pa.uky.edu/~verner/dima/photo/photo.dat; https://doi.org/10.1086/177435; He I ground-state support below the NIST He II 2p threshold: https://physics.nist.gov/PhysRefData/Handbook/Tables/heliumtable6.htm"),
                license: "Published numerical fit parameters; independently implemented equation".into(),
                artifact: None,
            },
            law: PhotoionizationLaw::Verner {
                fit,
                maximum_energy: Energy::new(maximum_ev * ELECTRONVOLT_J),
            },
            threshold: Energy::new(edge.max(gap)),
        })
    }

    /// Original atom/ion/electron identities, calorics and ground-state checks.
    pub const fn channel(&self) -> &AtomicIonizationV1 {
        &self.channel
    }
    /// Original photon source citation, usage terms and artifact identity.
    pub const fn provenance(&self) -> &Provenance {
        &self.provenance
    }
    /// Effective physical onset [J]; zero cross section below it.
    pub const fn threshold(&self) -> Energy {
        self.threshold
    }
    /// Maximum of the immutable law inside its declared energy interval.
    pub(super) fn maximum_cross_section_in_domain(&self) -> Result<Area, PhotoionizationErrorV1> {
        match &self.law {
            PhotoionizationLaw::Table(points) => Ok(Area::new(
                points
                    .iter()
                    .map(|p| p.cross_section.value())
                    .fold(0.0, f64::max),
            )),
            // All three literal fits decrease above their admitted onset.
            // H I/He II have y=x>31 and c=5.5-P/2=4.0185, so
            // d ln(sigma)/d ln(x) < 2x/(x-1)-c < 0.
            // For He I, c=3.906, b=y1, w=yw, x>0. The derivative of
            // ((x-1)^2+w^2)*(x^2+b^2)^(-c/2) has sign -(x*Q(x)+2b^2),
            // Q=1.906*x^2-5.812*x+11.020285026. Q's leading term is
            // positive and discriminant -50.239309038224 < 0. The remaining
            // (1+sqrt(y/ya))^-P factor also decreases. Thus this is an
            // analytic maximum, not a sampled or hydrogen-only assumption.
            PhotoionizationLaw::Verner { .. } => self.cross_section(self.threshold),
        }
    }

    /// Last supported original source energy [J].
    pub fn maximum_energy(&self) -> Energy {
        match &self.law {
            PhotoionizationLaw::Table(points) => points[points.len() - 1].photon_energy,
            PhotoionizationLaw::Verner { maximum_energy, .. } => *maximum_energy,
        }
    }
    /// Unmodified original samples; absent for the analytic Verner fits.
    pub fn points(&self) -> Option<&[PhotoionizationCrossSectionPointV1]> {
        match &self.law {
            PhotoionizationLaw::Table(points) => Some(points.as_ref()),
            PhotoionizationLaw::Verner { .. } => None,
        }
    }

    /// Allocation-free evaluation: O(log(samples)) for tables, O(1) for the fit.
    /// No source law permits upper extrapolation.
    pub fn cross_section(&self, energy: Energy) -> Result<Area, PhotoionizationErrorV1> {
        let e = energy.value();
        if !nonnegative(e) {
            return Err(PhotoionizationErrorV1::Invalid(
                "nonnegative finite photon energy",
            ));
        }
        if e > self.maximum_energy().value() {
            return Err(PhotoionizationErrorV1::OutsideEnergyDomain);
        }
        if e < self.threshold.value() {
            return Ok(Area::new(0.0));
        }
        let sigma = match &self.law {
            PhotoionizationLaw::Verner { fit, .. } => {
                let x = (e / ELECTRONVOLT_J) / fit.e0 - fit.y0;
                let y = if fit.y1 == 0.0 { x } else { sqrt(x * x + fit.y1 * fit.y1) };
                let exponent = (0.5 * fit.p - 5.5) * ln(y)
                    - fit.p * ln(1.0 + sqrt(y / fit.ya));
                (fit.sigma0 * 1e-22) * ((x - 1.0) * (x - 1.0) + fit.yw * fit.yw)
                    * exp(exponent)
            }
            PhotoionizationLaw::Table(points) => {
                let upper = points.partition_point(|p| p.photon_energy.value() < e);
                if points[upper].photon_energy.value() == e {
                    return Ok(points[upper].cross_section);
                }
                let left = points[upper - 1];
                let right = points[upper];
                let fraction = (e - left.photon_energy.value())
                    / (right.photon_energy.value() - left.photon_energy.value());
                let a = left.cross_section.value();
                let b = right.cross_section.value();
                // Avoid cancellation on decreasing edges.
                if a <= b {
                    a + fraction * (b - a)
                } else {
                    b + (1.0 - fraction) * (a - b)
                }
            }
        };
        if !nonnegative(sigma) {
            return Err(PhotoionizationErrorV1::Arithmetic(
                "photon cross section",
            ));
        }
        Ok(Area::new(sigma))
    }

    /// Local optically thin response to explicitly supplied discrete photons.
    ///
    /// This reads the local field; it neither attenuates photons nor consumes
    /// neutral inventory. A finite cell owner must close both ledgers. With
    /// neutral molar concentration c_A, q=c_A*rate is mol/(m³ s), with source
    /// [-q,+q,+q] in [A,I,e] order. At fixed local field, dq/dc_A=rate and
    /// dq/dc_I=dq/dc_e=0; the bank's canonical indices are on channel().
    /// Per-neutral powers multiply c_A*NA to give W/m³. Electron heating is a
    /// sensible-energy source; changing electron population also changes its
    /// heat capacity and must not be mistaken for an extra energy debit.
    ///
    /// Input order defines the compensated sums. Evaluation allocates nothing,
    /// polls cancellation every 32 groups and before publication, and refuses
    /// unsupported energies even for zero-weight groups. Empty input is dark.
    pub fn evaluate(
        &self,
        groups: &[PhotonNumberFluxV1],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PhotoionizationEvaluationV1, PhotoionizationErrorV1> {
        if cancelled() {
            return Err(PhotoionizationErrorV1::Cancelled);
        }
        if groups.len() > MAX_PHOTON_GROUPS_V1 {
            return Err(PhotoionizationErrorV1::Invalid("photon group count"));
        }
        let gap = self.channel.ionization_energy_j_mol() / AVOGADRO_PER_MOL;
        let mut sums = [0.0; 3];
        let mut corrections = [0.0; 3];
        for (index, group) in groups.iter().enumerate() {
            if index != 0 && index % 32 == 0 && cancelled() {
                return Err(PhotoionizationErrorV1::Cancelled);
            }
            let energy = group.photon_energy.value();
            if !energy.is_finite()
                || energy <= 0.0
                || !nonnegative(group.number_flux_per_m2_s)
            {
                return Err(PhotoionizationErrorV1::Invalid(
                    "positive photon energy and nonnegative flux",
                ));
            }
            let sigma = self.cross_section(group.photon_energy)?.value();
            let rate = sigma * group.number_flux_per_m2_s;
            let values = [
                rate,
                rate * energy,
                rate * (energy - gap).max(0.0),
            ];
            if values.iter().any(|value| !value.is_finite()) {
                return Err(PhotoionizationErrorV1::Arithmetic(
                    "photon rate or energy moment",
                ));
            }
            for j in 0..3 {
                let corrected = values[j] - corrections[j];
                let next = sums[j] + corrected;
                corrections[j] = (next - sums[j]) - corrected;
                sums[j] = next;
            }
        }
        let chemical = gap * sums[0];
        if sums.iter().any(|value| !nonnegative(*value)) || !nonnegative(chemical) {
            return Err(PhotoionizationErrorV1::Arithmetic(
                "accumulated photon moments",
            ));
        }
        if cancelled() {
            return Err(PhotoionizationErrorV1::Cancelled);
        }
        Ok(PhotoionizationEvaluationV1 {
            rate_per_neutral_per_s: sums[0],
            absorbed_power_per_neutral_w: sums[1],
            chemical_power_per_neutral_w: chemical,
            electron_heating_per_neutral_w: sums[2],
        })
    }
}

fn nonnegative(value: f64) -> bool {
    value.is_finite() && value >= 0.0
}

/// Local one-photon event and energy moments, before finite inventory evolution.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhotoionizationEvaluationV1 {
    /// Single-ionization probability rate per neutral atom [s⁻¹].
    pub rate_per_neutral_per_s: f64,
    /// Energy removed from the local photon field per neutral atom [W].
    pub absorbed_power_per_neutral_w: f64,
    /// Energy placed in original ion/electron formation levels per neutral [W].
    pub chemical_power_per_neutral_w: f64,
    /// Excess photoelectron kinetic energy per neutral atom [W].
    pub electron_heating_per_neutral_w: f64,
}
