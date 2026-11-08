//! Reusable linear polarization response of one original absorbing film.
//!
//! Two independent Cartesian incident fields prepare the Maxwell response
//! once. Later Jones/analyzer queries combine their anchored E/H amplitudes
//! before tensor absorption. This retains elliptical interference and both
//! reflected internal modes without another material query or boundary solve.

use super::super::super::slab::kernel::{finite, flux};
use super::{
    AbsorbingSlabResponse, C0, Incidence, OrientedUniaxialSlab, PreparedAbsorbingSlabField,
    ReflectionError, Result, Vector, positive,
};
use fs_math::{c64::C64, dd::Dd, det};

/// Two fixed-size Maxwell responses at one original frequency and direction.
/// The stored input basis is global tangential `(Ex, Ey)`, in peak V/m.
/// Its amplitudes need not be orthonormal in incident power at oblique angles.
pub struct PreparedAbsorbingSlabPolarization<'f, 's, 'm> {
    basis: [PreparedAbsorbingSlabField<'f, 's, 'm>; 2],
    incident_phase_index: f64,
    incident_permeability_h_m: f64,
}

impl<'f, 's, 'm> PreparedAbsorbingSlabPolarization<'f, 's, 'm> {
    pub const fn source(&self) -> &'f OrientedUniaxialSlab<'s, 'm> {
        self.basis[0].source
    }

    pub const fn omega_rad_s(&self) -> f64 {
        self.basis[0].response.omega
    }

    /// Resolved phase index of the original lossless incident half-space.
    pub const fn incident_phase_index(&self) -> f64 {
        self.incident_phase_index
    }

    /// Original incident magnetic permeability [H/m]. Together with phase
    /// index this lets a radiance consumer admit its physical exterior medium.
    pub const fn incident_permeability_h_m(&self) -> f64 {
        self.incident_permeability_h_m
    }

    /// Reuse the original Maxwell solution for arbitrary finite nonzero Jones
    /// amplitudes. All four internal waves combine coherently BEFORE loss.
    /// The returned field retains the same exact source and volume-integral API.
    /// No allocations, material evaluations, or boundary solves occur here.
    pub fn prepare_field(
        &self,
        electric_v_m: [C64; 2],
    ) -> Result<PreparedAbsorbingSlabField<'f, 's, 'm>> {
        Ok(PreparedAbsorbingSlabField {
            source: self.source(),
            response: combine_responses(
                [&self.basis[0].response, &self.basis[1].response],
                electric_v_m,
            )?,
        })
    }
}

impl<'s, 'm> OrientedUniaxialSlab<'s, 'm> {
    /// Prepare the two global tangential polarization basis fields at one
    /// frequency. The original positive incidence cosine is retained even
    /// when a near-grazing tangential index would round to the light cone.
    /// `plane_axis` is the global transverse propagation direction; reversing
    /// it reverses the transverse wavevector without reorienting the crystal.
    pub fn prepare_absorbing_polarization(
        &self,
        omega_rad_s: f64,
        incidence_cosine: f64,
        plane_axis: [f64; 2],
    ) -> Result<PreparedAbsorbingSlabPolarization<'_, 's, 'm>> {
        let omega = positive(omega_rad_s)?;
        let source = self.entrance_interface();
        let incident_permeability_h_m = source.incident_mu_h_m;
        let incident_phase_index = positive(
            (C0 * det::sqrt(incident_permeability_h_m))
                * det::sqrt(source.incident.lossless_permittivity(omega)?),
        )?;
        let basis = [
            self.prepare_absorbing(
                omega,
                Incidence::Cosine(incidence_cosine),
                plane_axis,
                [C64::ONE, C64::ZERO],
            )?,
            self.prepare_absorbing(
                omega,
                Incidence::Cosine(incidence_cosine),
                plane_axis,
                [C64::ZERO, C64::ONE],
            )?,
        ];
        Ok(PreparedAbsorbingSlabPolarization {
            basis,
            incident_phase_index,
            incident_permeability_h_m,
        })
    }
}

fn bad(message: &'static str) -> ReflectionError {
    ReflectionError::Invalid(message)
}

fn combine_responses(
    basis: [&AbsorbingSlabResponse; 2],
    coefficients: [C64; 2],
) -> Result<AbsorbingSlabResponse> {
    if coefficients.iter().any(|&v| !finite(v)) || coefficients == [C64::ZERO; 2] {
        return Err(bad(
            "film polarization requires finite nonzero incident amplitudes",
        ));
    }
    let [a, b] = basis;
    // Both basis fields originate from the same private preparation above.
    // Keep the numerical seam honest if an internal caller ever changes that.
    if a.thickness_m != b.thickness_m
        || a.omega != b.omega
        || a.plane != b.plane
        || a.axis != b.axis
        || a.loss_f_m != b.loss_f_m
        || a.tangential_wavenumber_rad_m != b.tangential_wavenumber_rad_m
        || a.incident_normal_wavenumber_rad_m != b.incident_normal_wavenumber_rad_m
        || a.exit_normal_wavenumber_rad_m != b.exit_normal_wavenumber_rad_m
        || a.right_incident_electric_v_m != [C64::ZERO; 3]
        || b.right_incident_electric_v_m != [C64::ZERO; 3]
        || a.right_incident_magnetic_a_m != [C64::ZERO; 3]
        || b.right_incident_magnetic_a_m != [C64::ZERO; 3]
        || (0..2).any(|j| {
            a.forward[j].normal_wavenumber_rad_m != b.forward[j].normal_wavenumber_rad_m
                || a.backward[j].normal_wavenumber_rad_m != b.backward[j].normal_wavenumber_rad_m
        })
    {
        return Err(bad(
            "film polarization basis must share its original mode geometry",
        ));
    }
    let mut response = *a;
    response.incident_electric_v_m = combine(
        a.incident_electric_v_m,
        b.incident_electric_v_m,
        coefficients,
    )?;
    response.incident_magnetic_a_m = combine(
        a.incident_magnetic_a_m,
        b.incident_magnetic_a_m,
        coefficients,
    )?;
    response.reflected_electric_v_m = combine(
        a.reflected_electric_v_m,
        b.reflected_electric_v_m,
        coefficients,
    )?;
    response.reflected_magnetic_a_m = combine(
        a.reflected_magnetic_a_m,
        b.reflected_magnetic_a_m,
        coefficients,
    )?;
    response.transmitted_electric_v_m = combine(
        a.transmitted_electric_v_m,
        b.transmitted_electric_v_m,
        coefficients,
    )?;
    response.transmitted_magnetic_a_m = combine(
        a.transmitted_magnetic_a_m,
        b.transmitted_magnetic_a_m,
        coefficients,
    )?;
    for j in 0..2 {
        response.forward[j].electric_v_m = combine(
            a.forward[j].electric_v_m,
            b.forward[j].electric_v_m,
            coefficients,
        )?;
        response.forward[j].magnetic_a_m = combine(
            a.forward[j].magnetic_a_m,
            b.forward[j].magnetic_a_m,
            coefficients,
        )?;
        response.backward[j].electric_v_m = combine(
            a.backward[j].electric_v_m,
            b.backward[j].electric_v_m,
            coefficients,
        )?;
        response.backward[j].magnetic_a_m = combine(
            a.backward[j].magnetic_a_m,
            b.backward[j].magnetic_a_m,
            coefficients,
        )?;
    }
    let incoming = positive(flux(
        response.incident_electric_v_m,
        response.incident_magnetic_a_m,
    ))?;
    let reflected = -flux(
        response.reflected_electric_v_m,
        response.reflected_magnetic_a_m,
    ) / incoming;
    let mut transmitted = flux(
        response.transmitted_electric_v_m,
        response.transmitted_magnetic_a_m,
    ) / incoming;
    let tolerance = 8192.0 * f64::EPSILON;
    if response.exit_normal_wavenumber_rad_m.re == 0.0 {
        if !transmitted.is_finite() || transmitted.abs() > tolerance {
            return Err(bad(
                "film polarization nonpropagating exit carries normal power",
            ));
        }
        transmitted = 0.0;
    }
    let absorbed = (Dd::ONE - Dd::from_f64(reflected) - Dd::from_f64(transmitted)).to_f64();
    if !reflected.is_finite()
        || !transmitted.is_finite()
        || !absorbed.is_finite()
        || reflected < 0.0
        || transmitted < 0.0
        || absorbed < -tolerance
        || response.loss_f_m == [0.0; 2] && absorbed.abs() > tolerance
    {
        return Err(bad(
            "film polarization violates passive normal-power balance",
        ));
    }
    response.incident_normal_power_w_m2 = incoming;
    response.reflected_power_fraction = reflected;
    response.transmitted_power_fraction = transmitted;
    response.absorbed_power_fraction = if response.loss_f_m == [0.0; 2] {
        0.0
    } else {
        absorbed.max(0.0)
    };
    Ok(response)
}

fn combine(a: Vector, b: Vector, coefficients: [C64; 2]) -> Result<Vector> {
    let product = |x: C64, y: C64| {
        let [xr, xi, yr, yi] = [x.re, x.im, y.re, y.im].map(Dd::from_f64);
        [xr * yr - xi * yi, xr * yi + xi * yr]
    };
    let result = core::array::from_fn(|j| {
        let x = product(a[j], coefficients[0]);
        let y = product(b[j], coefficients[1]);
        C64::new((x[0] + y[0]).to_f64(), (x[1] + y[1]).to_f64())
    });
    if result.iter().any(|&v| !finite(v)) {
        return Err(bad("film polarization fields are unrepresentable"));
    }
    Ok(result)
}

#[cfg(test)]
mod tests;
