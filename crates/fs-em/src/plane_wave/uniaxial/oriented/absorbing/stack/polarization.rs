//! Cached polarization-dependent loss in the original ordered crystal layers.
//!
//! Two global tangential basis fields share one material-frequency resolution.
//! Each later query combines only the selected layer's electric eigenwaves
//! before evaluating its original tensor loss. The other layers' reflections
//! remain present in those waves; no stack copy or new boundary solve is needed.

use super::super::super::slab::kernel::{finite, flux};
use super::super::slab::polarization::combine;
use super::{
    AbsorbingStackLayerField, AbsorbingStackResponse, Incidence, OrientedUniaxialStack,
    PreparedAbsorbingStackField, Result, bad, kernel, poll, unit_axis,
};
use fs_math::c64::C64;

/// Two retained Maxwell responses at one original frequency and direction.
/// Input coefficients are global tangential `(Ex, Ey)` peak phasors [V/m].
/// Their physical incident-power metric is retained at oblique incidence.
pub struct PreparedAbsorbingStackPolarization<'f, 's, 'm> {
    basis: [PreparedAbsorbingStackField<'f, 's, 'm>; 2],
    incident_phase_index: f64,
    incident_permeability_h_m: f64,
}

impl<'f, 's, 'm> PreparedAbsorbingStackPolarization<'f, 's, 'm> {
    pub const fn source(&self) -> &'f OrientedUniaxialStack<'s, 'm> {
        self.basis[0].source
    }

    pub const fn omega_rad_s(&self) -> f64 {
        self.basis[0].response.omega
    }

    pub const fn incident_phase_index(&self) -> f64 {
        self.incident_phase_index
    }

    pub const fn incident_permeability_h_m(&self) -> f64 {
        self.incident_permeability_h_m
    }

    /// Positive inward normal power for finite nonzero Jones coefficients.
    /// E/H combine before the Poynting contraction, preserving the original
    /// incident medium and propagation direction without a source query.
    pub fn incident_normal_power_w_m2(&self, electric_v_m: [C64; 2]) -> Result<f64> {
        incident_power(
            [&self.basis[0].response, &self.basis[1].response],
            electric_v_m,
        )
    }

    /// Exact dielectric work in one original layer prism. Depths are LOCAL
    /// to that layer's entrance and area is parallel to the stack faces.
    /// Four coherent electric waves are combined before the analytic loss
    /// integral. This fixed-work query allocates nothing and never resolves
    /// materials or traverses another layer.
    pub fn layer_absorbed_power_w(
        &self,
        layer_index: usize,
        electric_v_m: [C64; 2],
        area_m2: f64,
        start_depth_m: f64,
        end_depth_m: f64,
    ) -> Result<f64> {
        let layer = |i: usize| {
            self.basis[i]
                .response
                .layers
                .get(layer_index)
                .ok_or_else(|| bad("unknown absorbing stack polarization layer"))
        };
        layer_power(
            [layer(0)?, layer(1)?],
            electric_v_m,
            area_m2,
            start_depth_m,
            end_depth_m,
        )
    }
}

impl<'s, 'm> OrientedUniaxialStack<'s, 'm> {
    /// Cache the two Cartesian incident polarizations of the original stack.
    /// The original positive incidence cosine is retained near grazing.
    /// `plane_axis` is the global transverse propagation direction; each
    /// material keeps its independently authored global optic axis.
    ///
    /// Both retained responses and one elimination workspace are admitted
    /// against the original preparation budget before any material allocation.
    /// Cancellation returns no partially prepared polarization owner.
    pub fn prepare_absorbing_polarization(
        &self,
        omega_rad_s: f64,
        incidence_cosine: f64,
        plane_axis: [f64; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedAbsorbingStackPolarization<'_, 's, 'm>> {
        poll(&mut cancelled)?;
        let workspace_bytes = kernel::workspace_bytes(self.layers().len(), 2)?
            .checked_add(core::mem::size_of::<
                PreparedAbsorbingStackPolarization<'_, 's, 'm>,
            >())
            .ok_or_else(|| bad("absorbing stack polarization workspace overflows"))?;
        if workspace_bytes > self.preparation_limits().maximum_workspace_bytes {
            return Err(bad(
                "absorbing stack polarization exceeds its original preparation budget",
            ));
        }
        let plane = unit_axis(plane_axis).map_err(bad)?;
        let data = self.resolve_absorbing_frequency(omega_rad_s, &mut cancelled)?;
        let incidence = Incidence::Cosine(incidence_cosine);
        let first = kernel::solve(
            &data,
            incidence,
            [C64::ONE, C64::ZERO],
            plane,
            &mut cancelled,
        )?;
        let second = kernel::solve(
            &data,
            incidence,
            [C64::ZERO, C64::ONE],
            plane,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        Ok(PreparedAbsorbingStackPolarization {
            basis: [
                PreparedAbsorbingStackField {
                    source: self,
                    response: first,
                },
                PreparedAbsorbingStackField {
                    source: self,
                    response: second,
                },
            ],
            incident_phase_index: data.ni,
            incident_permeability_h_m: data.mui,
        })
    }
}

fn admit_coefficients(coefficients: [C64; 2]) -> Result<()> {
    if coefficients.iter().any(|&v| !finite(v)) || coefficients == [C64::ZERO; 2] {
        return Err(bad(
            "stack polarization requires finite nonzero incident amplitudes",
        ));
    }
    Ok(())
}

fn incident_power(basis: [&AbsorbingStackResponse; 2], coefficients: [C64; 2]) -> Result<f64> {
    admit_coefficients(coefficients)?;
    let [a, b] = basis;
    let electric = combine(
        a.incident_electric_v_m,
        b.incident_electric_v_m,
        coefficients,
    )?;
    let magnetic = combine(
        a.incident_magnetic_a_m,
        b.incident_magnetic_a_m,
        coefficients,
    )?;
    let power = flux(electric, magnetic);
    if !power.is_finite() || power <= 0.0 {
        return Err(bad("stack polarization incident power is unrepresentable"));
    }
    Ok(power)
}

fn layer_power(
    basis: [&AbsorbingStackLayerField; 2],
    coefficients: [C64; 2],
    area: f64,
    start: f64,
    end: f64,
) -> Result<f64> {
    admit_coefficients(coefficients)?;
    let [a, b] = basis;
    if a.entrance_z_m != b.entrance_z_m
        || a.exit_z_m != b.exit_z_m
        || a.thickness_m != b.thickness_m
        || a.omega != b.omega
        || a.axis != b.axis
        || a.loss_f_m != b.loss_f_m
        || (0..2).any(|j| {
            a.forward[j].normal_wavenumber_rad_m != b.forward[j].normal_wavenumber_rad_m
                || a.backward[j].normal_wavenumber_rad_m != b.backward[j].normal_wavenumber_rad_m
        })
    {
        return Err(bad(
            "stack polarization basis must share its original layer geometry",
        ));
    }
    let mut layer = *a;
    // Only E enters dielectric work. The fixed-size temporary is consumed by
    // this loss query; magnetic amplitudes need no recombination or evaluation.
    for j in 0..2 {
        layer.forward[j].electric_v_m = combine(
            a.forward[j].electric_v_m,
            b.forward[j].electric_v_m,
            coefficients,
        )?;
        layer.backward[j].electric_v_m = combine(
            a.backward[j].electric_v_m,
            b.backward[j].electric_v_m,
            coefficients,
        )?;
    }
    layer.absorbed_power_w(area, start, end)
}

#[cfg(test)]
mod tests;
