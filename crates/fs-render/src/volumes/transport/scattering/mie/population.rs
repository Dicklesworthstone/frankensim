//! Independent particles: sum cross sections and intensities, never amplitudes.
//!
//! For number fractions `f_i`, the mean cross section is `sum(f_i C_i)` and
//! the particle-mass coefficient is `sum(f_i C_i) / sum(f_i m_i)`. The phase
//! and its first moment use weights `f_i C_sca,i / sum(f_i C_sca,i)`.
//! See Frisvad, Christensen & Jensen (2007), §2.2, equations 27--32:
//! <https://cseweb.ucsd.edu/~henrik/papers/lorenz_mie_theory/computing_scattering_properties_using_lorenz_mie_theory.pdf>.
//! Independent-particle cross sections can feed the existing multiple-scattering
//! radiative transport. Near-field coupling, coherent interparticle scattering,
//! aggregates and concentration-dependent constitutive corrections are excluded.

use super::distribution::{Compensated, normalize_number_bins};
use super::{
    MieAngularSolution, MieBudget, MieCrossSections, MieError, MiePhaseTable, MieSolution,
    MieSphere, poll, prepare_phase_table,
};
use fs_exec::Cx;

/// Hard bound on retained particle populations and radius quadrature nodes.
pub const MAX_MIE_POPULATION_BINS: usize = 256;

/// Preparation limits for one wavelength, including every particle bin.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MiePopulationBudget {
    /// Per-sphere recurrence and multipole limits.
    pub sphere: MieBudget,
    /// Maximum number of positive-number bins to solve.
    pub max_bins: usize,
    /// Maximum retained multipole pairs across all particle bins.
    pub max_total_multipoles: usize,
}

impl Default for MiePopulationBudget {
    fn default() -> Self {
        Self {
            sphere: MieBudget::default(),
            max_bins: MAX_MIE_POPULATION_BINS,
            max_total_multipoles: MAX_MIE_POPULATION_BINS * 160,
        }
    }
}

/// A number-weighted population of homogeneous spheres in one real host.
///
/// Radii, composition and material density may differ between bins. Weights
/// are relative particle counts, not mass, volume or scattering fractions.
/// Multiplying every weight by the same factor leaves the population unchanged.
/// Absolute concentration is supplied separately by the transport density field.
#[derive(Clone, Debug, PartialEq)]
pub struct MiePopulation {
    bins: Vec<(MieSphere, f64)>,
    host_index_n: f64,
    mean_particle_mass_kg: f64,
}

impl MiePopulation {
    /// Admit 1--256 bins with finite nonnegative number weights and at least
    /// one positive weight. Zero-number bins are discarded; all retained bins
    /// must have exactly the same nonabsorbing host index. Positive number or
    /// mass contributions that underflow are refused, not silently discarded.
    pub fn try_new(bins: &[(MieSphere, f64)]) -> Result<Self, MieError> {
        let bins = normalize_number_bins(bins)?;
        let host_index_n = bins[0].0.host_index_n;
        let mut mass = Compensated::default();
        for &(sphere, fraction) in &bins {
            if sphere.host_index_n != host_index_n {
                return Err(MieError::Parameters);
            }
            mass.add(positive_product(fraction, sphere.particle_mass_kg()?)?);
        }
        let mean_particle_mass_kg = mass.total();
        if !mean_particle_mass_kg.is_finite() || mean_particle_mass_kg <= 0.0 {
            return Err(MieError::Numerical);
        }
        Ok(Self {
            bins,
            host_index_n,
            mean_particle_mass_kg,
        })
    }

    /// Retained spheres and normalized number fractions, in declaration order.
    #[must_use]
    pub fn bins(&self) -> &[(MieSphere, f64)] {
        &self.bins
    }

    /// Common real nonabsorbing host index.
    #[must_use]
    pub const fn host_index_n(&self) -> f64 {
        self.host_index_n
    }

    /// Number-mean particle mass, kg. Number density times this value equals
    /// the dispersed-particle mass density expected by a scattering grid.
    #[must_use]
    pub const fn mean_particle_mass_kg(&self) -> f64 {
        self.mean_particle_mass_kg
    }

    /// Solve every retained particle at one vacuum wavelength in nanometres.
    /// Cancellation is polled between bins and inside each native sphere solve.
    /// The total multipole budget is passed down before each solve, so a
    /// partially prepared population never escapes a budget refusal.
    pub fn at_wavelength(
        &self,
        cx: &Cx<'_>,
        wavelength_nm: f64,
        budget: MiePopulationBudget,
    ) -> Result<MiePopulationSolution, MieError> {
        poll(cx)?;
        if self.bins.len() > budget.max_bins {
            return Err(MieError::Budget);
        }
        let mut members = Vec::new();
        members
            .try_reserve_exact(self.bins.len())
            .map_err(|_| MieError::Budget)?;
        let (mut scattering, mut absorption) = (Compensated::default(), Compensated::default());
        let (mut total_multipoles, mut largest_multipoles) = (0, 0);
        for &(sphere, fraction) in &self.bins {
            poll(cx)?;
            let remaining = budget
                .max_total_multipoles
                .checked_sub(total_multipoles)
                .ok_or(MieError::Budget)?;
            let solution = sphere.at_wavelength(
                cx,
                wavelength_nm,
                MieBudget {
                    max_multipoles: budget.sphere.max_multipoles.min(remaining),
                    max_recurrence_order: budget.sphere.max_recurrence_order,
                },
            )?;
            let sections = solution.cross_sections();
            let scattered = positive_product(fraction, sections.scattering_m2)?;
            scattering.add(scattered);
            absorption.add(positive_product(fraction, sections.absorption_m2)?);
            total_multipoles += solution.multipoles();
            largest_multipoles = largest_multipoles.max(solution.multipoles());
            members.push((solution, scattered));
        }
        let (scattering_m2, absorption_m2) = (scattering.total(), absorption.total());
        let mut moment = Compensated::default();
        if scattering_m2 > 0.0 && scattering_m2.is_finite() {
            for (solution, weight) in &mut members {
                let raw = *weight;
                *weight /= scattering_m2;
                if raw > 0.0 && *weight == 0.0 {
                    return Err(MieError::Numerical);
                }
                moment.add(*weight * solution.cross_sections().mean_cosine);
            }
        }
        let sections = if members.len() == 1 {
            // Preserve the native monodisperse coefficients, including their
            // stable Q/rho/r arithmetic and exact bits.
            members[0].0.cross_sections()
        } else {
            MieCrossSections {
                scattering_m2,
                absorption_m2,
                mass_scattering_m2_kg: scattering_m2 / self.mean_particle_mass_kg,
                mass_absorption_m2_kg: absorption_m2 / self.mean_particle_mass_kg,
                mean_cosine: moment.total(),
            }
        };
        if [
            sections.scattering_m2,
            sections.absorption_m2,
            sections.mass_scattering_m2_kg,
            sections.mass_absorption_m2_kg,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v < 0.0)
            || (scattering_m2 > 0.0 && sections.mass_scattering_m2_kg == 0.0)
            || (absorption_m2 > 0.0 && sections.mass_absorption_m2_kg == 0.0)
            || !(-1.0..=1.0).contains(&sections.mean_cosine)
            || !sections.extinction_m2().is_finite()
            || !sections.transport_m2().is_finite()
        {
            return Err(MieError::Numerical);
        }
        poll(cx)?;
        Ok(MiePopulationSolution {
            wavelength_nm,
            mean_particle_mass_kg: self.mean_particle_mass_kg,
            members,
            total_multipoles,
            largest_multipoles,
            sections,
        })
    }
}

/// Prepared independent-particle population at one vacuum wavelength.
///
/// Full-series angle evaluation allocates nothing. Rendering should use
/// [`Self::phase_table`] or [`MiePopulation::visible_spectrum`] to avoid
/// re-evaluating all particles at path vertices.
#[derive(Clone, Debug, PartialEq)]
pub struct MiePopulationSolution {
    wavelength_nm: f64,
    mean_particle_mass_kg: f64,
    members: Vec<(MieSolution, f64)>,
    total_multipoles: usize,
    largest_multipoles: usize,
    sections: MieCrossSections,
}

impl MiePopulationSolution {
    /// Vacuum wavelength, nm.
    #[must_use]
    pub const fn wavelength_nm(&self) -> f64 {
        self.wavelength_nm
    }
    /// Number-mean particle mass, kg.
    #[must_use]
    pub const fn mean_particle_mass_kg(&self) -> f64 {
        self.mean_particle_mass_kg
    }
    /// Largest member series order, which determines the angular polynomial
    /// degree needed to resolve the entire incoherent population.
    #[must_use]
    pub const fn multipoles(&self) -> usize {
        self.largest_multipoles
    }
    /// Retained multipole pairs across all bins, for work accounting.
    #[must_use]
    pub const fn total_multipoles(&self) -> usize {
        self.total_multipoles
    }
    /// Cross sections are number means per particle; mass coefficients divide
    /// those means by the number-mean particle mass, not by a mean radius.
    /// Zero-scattering populations use `mean_cosine = 0` as a storage convention;
    /// they have no angular law and angle evaluation returns `NoScattering`.
    #[must_use]
    pub const fn cross_sections(&self) -> MieCrossSections {
        self.sections
    }

    /// Full scattering-weighted, unpolarized phase density in sr⁻¹.
    /// Work is bounded by 256 members, each within the native Mie domain.
    pub fn eval_cos(&self, cosine: f64) -> Result<f64, MieError> {
        if !cosine.is_finite() || cosine.abs() > 1.0 + 2e-10 {
            return Err(MieError::Parameters);
        }
        if self.sections.scattering_m2 == 0.0 {
            return Err(MieError::NoScattering);
        }
        let mut density = Compensated::default();
        for (solution, weight) in &self.members {
            if *weight > 0.0 {
                density.add(*weight * solution.eval_cos(cosine)?);
            }
        }
        let value = density.total();
        if value.is_finite() && value >= 0.0 {
            Ok(value)
        } else {
            Err(MieError::Numerical)
        }
    }

    /// Prepare the population's normalized native angular table, with the
    /// same sampled/refinement diagnostics as [`MieSolution::phase_table`].
    pub fn phase_table(
        &self,
        cx: &Cx<'_>,
        tolerance: f64,
        max_knots: usize,
    ) -> Result<MiePhaseTable, MieError> {
        prepare_phase_table(self, cx, tolerance, max_knots)
    }
}

impl MieAngularSolution for MiePopulationSolution {
    fn cross_sections(&self) -> MieCrossSections {
        self.cross_sections()
    }
    fn multipoles(&self) -> usize {
        self.multipoles()
    }
    fn eval_cos(&self, cosine: f64) -> Result<f64, MieError> {
        self.eval_cos(cosine)
    }
}

fn positive_product(left: f64, right: f64) -> Result<f64, MieError> {
    let product = left * right;
    if product.is_finite() && !(left > 0.0 && right > 0.0 && product == 0.0) {
        Ok(product)
    } else {
        Err(MieError::Numerical)
    }
}

#[cfg(test)]
mod tests;
