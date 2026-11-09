//! Cached polarization-dependent fields and loss in ordered crystal layers.
//!
//! Two global tangential basis fields share one material-frequency resolution.
//! Later queries combine anchored eigenwaves before propagation and tensor
//! loss. Independent positive factors add second moments, never phasors.
//! The other layers' reflections remain present in those waves; no stack
//! copy or new boundary solve is needed.

use super::super::super::slab::kernel::{finite, flux};
use super::super::slab::polarization::combine;
use super::super::{AbsorbingCrystalField, AbsorbingCrystalStatistics, kernel as crystal};
use super::{
    AbsorbingStackLayerField, AbsorbingStackResponse, FieldSum, Incidence, OrientedUniaxialStack,
    PreparedAbsorbingStackField, Result, StackRegion, WaveGeometry, bad, kernel, poll, unit_axis,
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

    /// Full Cartesian peak E/H, signed cycle-mean Poynting vector and local
    /// tensor heat for arbitrary nonzero Jones amplitudes in global Ex/Ey.
    /// The selected layer's four anchored waves combine BEFORE propagation;
    /// this preserves ellipticity, cross-polarized reflection and standing
    /// waves. Exterior observations retain their original side of each face.
    /// Fixed work and storage: no allocation, stack traversal or source solve.
    pub fn field_at(
        &self,
        electric_v_m: [C64; 2],
        region: StackRegion,
        point_m: [f64; 3],
    ) -> Result<AbsorbingCrystalField> {
        field_at(
            [&self.basis[0].response, &self.basis[1].response],
            electric_v_m,
            region,
            point_m,
        )
    }

    /// Partially polarized illumination at this prepared frequency/direction.
    /// Each Jones factor includes its square-root statistical weight in peak
    /// V/m and is explicitly uncorrelated with the others, even at the same
    /// frequency. Factors add coherency, signed flux and heat; their complex
    /// fields never interfere. No bandwidth, phase averaging or RMS factor
    /// is inferred. Empty/zero factors are dark and still validate geometry.
    /// Work is O(factors), with fixed storage and cancellation per factor.
    pub fn field_statistics(
        &self,
        factors_v_m: &[[C64; 2]],
        maximum_factors: usize,
        region: StackRegion,
        point_m: [f64; 3],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AbsorbingCrystalStatistics> {
        field_statistics(
            [&self.basis[0].response, &self.basis[1].response],
            factors_v_m,
            maximum_factors,
            region,
            point_m,
            &mut cancelled,
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
    combined_layer(basis, coefficients, false)?.absorbed_power_w(area, start, end)
}

fn combined_layer(
    basis: [&AbsorbingStackLayerField; 2],
    coefficients: [C64; 2],
    include_magnetic: bool,
) -> Result<AbsorbingStackLayerField> {
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
    // Heat-only queries retain their original electric-only cost.
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
        if include_magnetic {
            layer.forward[j].magnetic_a_m = combine(
                a.forward[j].magnetic_a_m,
                b.forward[j].magnetic_a_m,
                coefficients,
            )?;
            layer.backward[j].magnetic_a_m = combine(
                a.backward[j].magnetic_a_m,
                b.backward[j].magnetic_a_m,
                coefficients,
            )?;
        }
    }
    Ok(layer)
}

fn admit_point(
    response: &AbsorbingStackResponse,
    region: StackRegion,
    point: [f64; 3],
) -> Result<()> {
    if point.iter().any(|v| !v.is_finite()) {
        return Err(bad("finite absorbing stack field point required"));
    }
    let valid = if let StackRegion::Layer(index) = region {
        let layer = response
            .layers
            .get(index)
            .ok_or_else(|| bad("unknown absorbing stack polarization layer"))?;
        (layer.entrance_z_m..=layer.exit_z_m).contains(&point[2])
    } else if region == StackRegion::Transmitted {
        point[2] >= response.thickness_m
    } else {
        point[2] <= 0.0
    };
    if valid {
        Ok(())
    } else {
        Err(bad("absorbing stack point is outside its selected region"))
    }
}

fn field_at(
    basis: [&AbsorbingStackResponse; 2],
    coefficients: [C64; 2],
    region: StackRegion,
    point: [f64; 3],
) -> Result<AbsorbingCrystalField> {
    admit_coefficients(coefficients)?;
    let [a, b] = basis;
    admit_point(a, region, point)?;
    if let StackRegion::Layer(index) = region {
        let layer = combined_layer([&a.layers[index], &b.layers[index]], coefficients, true)?;
        let mut sum = layer.field_sum(point)?;
        layer.add_waves(&mut sum, a.tangential_wavenumber_rad_m, a.plane)?;
        let field = sum.finish().map_err(bad)?;
        let heat = crystal::heating(field.electric_v_m, layer.axis, layer.loss_f_m, layer.omega)
            .map_err(bad)?;
        return Ok(AbsorbingCrystalField {
            field,
            absorbed_power_density_w_m3: heat,
        });
    }
    let mut sum = FieldSum::at_point(point).map_err(bad)?;
    let mut add = |ea, eb, ha, hb, normal, anchor| {
        sum.add_wave(
            &combine(ea, eb, coefficients)?,
            &combine(ha, hb, coefficients)?,
            WaveGeometry {
                tangential_wavenumber_rad_m: a.tangential_wavenumber_rad_m,
                incidence_plane_axis: a.plane,
                normal_wavenumber_rad_m: normal,
                anchor_z_m: anchor,
            },
        )
        .map_err(bad)
    };
    if matches!(
        region,
        StackRegion::Incident | StackRegion::TotalIncidentSide
    ) {
        add(
            a.incident_electric_v_m,
            b.incident_electric_v_m,
            a.incident_magnetic_a_m,
            b.incident_magnetic_a_m,
            C64::from_re(a.incident_normal_wavenumber_rad_m),
            0.0,
        )?;
    }
    if matches!(
        region,
        StackRegion::Reflected | StackRegion::TotalIncidentSide
    ) {
        add(
            a.reflected_electric_v_m,
            b.reflected_electric_v_m,
            a.reflected_magnetic_a_m,
            b.reflected_magnetic_a_m,
            C64::from_re(-a.incident_normal_wavenumber_rad_m),
            0.0,
        )?;
    }
    if region == StackRegion::Transmitted {
        add(
            a.transmitted_electric_v_m,
            b.transmitted_electric_v_m,
            a.transmitted_magnetic_a_m,
            b.transmitted_magnetic_a_m,
            a.exit_normal_wavenumber_rad_m,
            a.thickness_m,
        )?;
    }
    Ok(AbsorbingCrystalField {
        field: sum.finish().map_err(bad)?,
        absorbed_power_density_w_m3: 0.0,
    })
}

fn field_statistics(
    basis: [&AbsorbingStackResponse; 2],
    factors: &[[C64; 2]],
    maximum_factors: usize,
    region: StackRegion,
    point: [f64; 3],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<AbsorbingCrystalStatistics> {
    poll(cancelled)?;
    if factors.len() > maximum_factors {
        return Err(bad(
            "stack polarization exceeds its independent-factor budget",
        ));
    }
    admit_point(basis[0], region, point)?;
    let mut sum = crystal::StatisticsSum::new();
    for &coefficients in factors {
        poll(cancelled)?;
        if coefficients != [C64::ZERO; 2] {
            sum.add(field_at(basis, coefficients, region, point)?);
        }
    }
    let result = sum.finish().map_err(|error| match error {
        crystal::IlluminationError::Cancelled => super::AngularSpectrumError::Cancelled,
        crystal::IlluminationError::Electromagnetic(message) => bad(message),
    })?;
    poll(cancelled)?;
    Ok(result)
}

#[cfg(test)]
mod tests;
