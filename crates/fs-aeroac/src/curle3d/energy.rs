//! Quadratic perturbation energy of the ORIGINAL coherent pressure/velocity.
//!
//! For a homogeneous linear fluid, E=p²/(2 K_s)+rho |u|²/2 and laboratory
//! transport is J=p u+U E. Contracting the source-free linearized pressure and
//! momentum equations gives partial_t E+div J=0 for constant mean velocity U.
//! All retained perturbation velocity contributes, including convected vortical
//! motion: this is not exclusively radiated sound power or total fluid energy.
//!
//! U is the original MEAN FLUID velocity, never source or microphone velocity.
//! Native adapters bind it from the accepted frame. These constant-work queries
//! add no flow/propagation model or Mach gate. Original material, linearization,
//! frame and source-domain admission remain on their owners. Quadratic outputs
//! can have twice the field bandwidth; sampled spectra and averaging must
//! resolve that bandwidth. No plane-wave conversion or source independence is
//! inferred; coherent fields must be summed before this reduction.

use super::CompactFieldStatistics;
use super::statistics::{invalid, nonfinite};
use crate::AeroacError;
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::dd::Dd;

/// Density and laboratory transport of the supplied perturbation fields.
/// Values are instantaneous or time averaged according to the input owner.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AcousticEnergyTransport {
    pub energy_density_j_m3: f64,
    /// Pressure work p*u (or its original time average), excluding mean advection.
    pub pressure_work_flux_w_m2: [f64; 3],
    /// Pressure work plus U*E in fixed laboratory Cartesian axes [W/m²].
    /// Flux relative to a moving control surface additionally subtracts V_surface*E.
    pub energy_flux_w_m2: [f64; 3],
}

/// Reduce one co-located, already coherently assembled real pressure/velocity.
/// The snapshot's original adiabatic modulus is retained, including for liquids.
/// Scaled DD squares avoid forming p*p, rho*c*c or rho*|u|² before scaling.
pub fn instantaneous(
    medium: FluidAcousticState,
    mean_flow_velocity_m_s: [f64; 3],
    pressure_pa: f64,
    particle_velocity_m_s: [f64; 3],
) -> Result<AcousticEnergyTransport, AeroacError> {
    if !pressure_pa.is_finite() || particle_velocity_m_s.iter().any(|v| !v.is_finite()) {
        return Err(nonfinite("perturbation energy pressure or velocity"));
    }
    let half = Dd::from_f64(0.5);
    let pressure =
        Dd::from_f64(pressure_pa) / Dd::from_f64(medium.adiabatic_bulk_modulus_pa()).sqrt();
    let density_scale = Dd::from_f64(medium.density_kg_m3()).sqrt();
    let mut energy = (half * pressure) * pressure;
    for velocity in particle_velocity_m_s {
        let scaled = density_scale * Dd::from_f64(velocity);
        energy = energy + (half * scaled) * scaled;
    }
    finish(
        energy,
        particle_velocity_m_s.map(|u| pressure_pa * u),
        mean_flow_velocity_m_s,
    )
}

/// Reduce ORIGINAL physical second moments without another phasor half-factor.
/// The statistics owner already applied its peak/RMS/bin averaging convention.
/// Covariance validity and source coherence remain that owner's responsibility;
/// this query admits the nonnegative real variances and finite pressure work
/// needed for the isotropic-fluid energy reduction, without reconstructing fields.
pub fn time_averaged(
    medium: FluidAcousticState,
    mean_flow_velocity_m_s: [f64; 3],
    field: &CompactFieldStatistics,
) -> Result<AcousticEnergyTransport, AeroacError> {
    let pressure = field.mean_square_pressure_pa2;
    if !pressure.is_finite() || pressure < 0.0 {
        return Err(invalid(
            "perturbation energy requires nonnegative pressure variance",
        ));
    }
    let modulus = medium.adiabatic_bulk_modulus_pa();
    // Preserve the established harmonic owner arithmetic in its normal range;
    // do not halve a subnormal variance before division by a small modulus.
    let mut energy = if pressure >= f64::MIN_POSITIVE {
        Dd::from_f64((0.5 * pressure) / modulus)
    } else {
        (Dd::from_f64(pressure) / Dd::from_f64(modulus)) * Dd::from_f64(0.5)
    };
    let density = medium.density_kg_m3();
    for axis in 0..3 {
        let variance = field.velocity_covariance_m2_s2[axis][axis];
        if !variance.re.is_finite() || variance.re < 0.0 || variance.im != 0.0 {
            return Err(invalid(
                "perturbation energy requires real nonnegative velocity variances",
            ));
        }
        let (a, b) = if density >= variance.re {
            (0.5 * density, variance.re)
        } else {
            (density, 0.5 * variance.re)
        };
        energy = energy + Dd::from_f64(a) * Dd::from_f64(b);
    }
    finish(energy, field.active_intensity_w_m2, mean_flow_velocity_m_s)
}

fn finish(
    energy: Dd,
    pressure_work_flux_w_m2: [f64; 3],
    mean_flow_velocity_m_s: [f64; 3],
) -> Result<AcousticEnergyTransport, AeroacError> {
    let energy_density_j_m3 = energy.to_f64();
    if !energy_density_j_m3.is_finite()
        || energy_density_j_m3 < 0.0
        || pressure_work_flux_w_m2
            .iter()
            .chain(&mean_flow_velocity_m_s)
            .any(|v| !v.is_finite())
    {
        return Err(nonfinite("perturbation energy, pressure work or mean flow"));
    }
    let energy_flux_w_m2 = core::array::from_fn(|axis| {
        if mean_flow_velocity_m_s[axis] == 0.0 {
            pressure_work_flux_w_m2[axis]
        } else {
            (Dd::from_f64(pressure_work_flux_w_m2[axis])
                + Dd::from_f64(mean_flow_velocity_m_s[axis]) * energy)
                .to_f64()
        }
    });
    if energy_flux_w_m2.iter().any(|v| !v.is_finite()) {
        return Err(nonfinite("transported perturbation energy flux"));
    }
    Ok(AcousticEnergyTransport {
        energy_density_j_m3,
        pressure_work_flux_w_m2,
        energy_flux_w_m2,
    })
}
