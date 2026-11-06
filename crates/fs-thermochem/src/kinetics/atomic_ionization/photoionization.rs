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
//! The tabulated source, interpolation and finite energy interval remain
//! explicit; no hydrogen fit, high-energy tail or spectral quadrature is
//! invented. Published atomic photoionization sources include Verner et al.,
//! ApJ 465, 487 (1996), <https://arxiv.org/abs/astro-ph/9601009>.

use super::AtomicIonizationV1;
use crate::electron_collision::{AVOGADRO_PER_MOL, MAX_CROSS_SECTION_POINTS_V1};
use fs_matdb::Provenance;
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
    /// Photon energy exceeds the last original source sample.
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

/// Immutable photon source on the original validated atomic species bank.
#[derive(Clone, Debug, PartialEq)]
pub struct AtomicPhotoionizationV1 {
    channel: AtomicIonizationV1,
    provenance: Provenance,
    points: std::sync::Arc<[PhotoionizationCrossSectionPointV1]>,
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
            points: points.into(),
            threshold: Energy::new(first.max(gap)),
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
    /// Last supported original source energy [J].
    pub fn maximum_energy(&self) -> Energy {
        self.points[self.points.len() - 1].photon_energy
    }
    /// Unmodified original interpolation samples.
    pub fn points(&self) -> &[PhotoionizationCrossSectionPointV1] {
        &self.points
    }

    /// Allocation-free O(log(samples)) evaluation with no upper extrapolation.
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
        let upper = self.points.partition_point(|p| p.photon_energy.value() < e);
        if self.points[upper].photon_energy.value() == e {
            return Ok(self.points[upper].cross_section);
        }
        let left = self.points[upper - 1];
        let right = self.points[upper];
        let fraction = (e - left.photon_energy.value())
            / (right.photon_energy.value() - left.photon_energy.value());
        let a = left.cross_section.value();
        let b = right.cross_section.value();
        // Preserve nonnegativity and avoid cancellation on decreasing edges.
        let sigma = if a <= b {
            a + fraction * (b - a)
        } else {
            b + (1.0 - fraction) * (a - b)
        };
        if !nonnegative(sigma) {
            return Err(PhotoionizationErrorV1::Arithmetic(
                "interpolated photon cross section",
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
