//! Finite rigid-frame porous layers in the existing pressure/volume-flow chain.
//!
//! A layer can transmit as well as reflect sound. Its JCA(L) density and bulk
//! modulus come from the original material owner; the existing stable two-port
//! supplies propagation, interfaces and composition with other layers or gaps.
//! The effective material coefficients already include open porosity: volume
//! flow is superficial velocity times the FULL exposed face area, with no
//! additional porosity factor. All port phasors use exp(-i omega t).
//!
//! This is a homogeneous, normally incident equivalent-fluid layer with a
//! motionless frame and open-pore pressure/volume-flow continuity. The caller
//! must justify the plane-wave and rigid-frame regimes and retain the source
//! support. Finite lateral diffraction, frame waves, contact skins, mean flow,
//! anisotropy and transient temperature evolution are outside this adapter.

use crate::compose::TwoPort;
use crate::fluid::FluidDuctError;
use crate::{DuctError, Segment, SegmentWave};
use fs_material::fluid_acoustics::thermoviscous::ThermoviscousFluid;
use fs_material::porous::{PorousFluidLaw, PorousFluidResponse};
use fs_math::c64::C64;

fn bad(what: &'static str) -> DuctError {
    DuctError::BadParameter { what }
}
fn finite(value: C64) -> bool {
    value.re.is_finite() && value.im.is_finite()
}

/// One source-owned, frequency-specific porous slab, with both faces exposed.
/// The prepared scattering block is reusable without evaluating the material
/// again. Connecting blocks requires the same positive real port reference.
#[derive(Debug, Clone)]
pub struct PorousLayer<'a> {
    material: &'a PorousFluidLaw,
    response: PorousFluidResponse,
    thickness_m: f64,
    area_m2: f64,
    two_port: TwoPort,
}

/// Coherent, cycle-averaged power of a porous slab under two-sided excitation.
/// Pressure phasors are peak amplitudes, so power contains the factor 1/2.
/// Right-port volume flow points INTO the layer, opposite to the axial flow.
#[derive(Debug, Clone, Copy)]
pub struct PorousLayerPower {
    /// Outgoing pressure-wave amplitudes at the left and right faces [Pa].
    pub outgoing_pressure_pa: [C64; 2],
    /// Total pressure at each face [Pa], including both incident fields.
    pub pressure_pa: [C64; 2],
    /// Volume flow into each face [m3/s], on the full exposed area.
    pub inward_volume_flow_m3_s: [C64; 2],
    /// Total power in the two prescribed incoming waves [W].
    pub incident_power_w: f64,
    /// Total power in the two coherently combined outgoing waves [W].
    pub outgoing_power_w: f64,
    /// Net acoustic power delivered to the layer [W]. This is the signed
    /// incoming-minus-outgoing balance, retaining floating-point roundoff;
    /// it is never clipped or confused with transmission to the other face.
    /// The rigid-frame JCA(L) model dissipates it through pore viscosity and
    /// heat exchange; this result does not partition heat between constituents.
    pub absorbed_power_w: f64,
}

impl<'a> PorousLayer<'a> {
    /// Prepare one bounded JCA(L) material query and a stable propagation block.
    /// `area_m2` is the full face area, not the open-pore area. The positive
    /// real `reference_impedance` is in Pa s/m3 and must match connected blocks.
    /// For plane waves in the same surrounding fluid on both faces, choosing
    /// its rho*c/area makes the incident/outgoing waves physical plane waves.
    ///
    /// Only decaying exponentials are used, even for an acoustically opaque
    /// slab. No transfer matrix containing exponentially growing terms is
    /// formed. Additional rigid-wall diffusion must not be added to the same
    /// pore losses. Preparation has fixed work; composition uses the existing
    /// block-budget and cancellation controls in [`crate::compose`].
    pub fn try_new<M: ThermoviscousFluid + ?Sized>(
        material: &'a PorousFluidLaw,
        fluid: &M,
        omega: f64,
        thickness_m: f64,
        area_m2: f64,
        reference_impedance: f64,
    ) -> Result<Self, FluidDuctError> {
        if !thickness_m.is_finite() || thickness_m <= 0.0 || !area_m2.is_finite() || area_m2 <= 0.0
        {
            return Err(
                bad("porous layer needs positive finite thickness and full face area").into(),
            );
        }
        TwoPort::through(reference_impedance)?;
        let response = material.effective_fluid(fluid, omega)?;
        let (specific, wavenumber) = response.characteristic();
        let specific = specific.conj();
        let characteristic = specific.scale(1.0 / area_m2);
        // Uniform plane-wave propagation depends only on area. This equivalent
        // radius supplies that area to the original cylindrical wave owner;
        // it does not add a circular-pore or circular-wall loss approximation.
        let radius = (area_m2 / core::f64::consts::PI).sqrt();
        if !finite(characteristic)
            || characteristic.abs() == 0.0
            || !radius.is_finite()
            || radius <= 0.0
        {
            return Err(bad(
                "porous layer volume impedance or equivalent radius is unrepresentable",
            )
            .into());
        }
        let two_port = TwoPort::from_wave(
            &Segment::Cylinder {
                radius,
                length: thickness_m,
            },
            &SegmentWave {
                wavenumber: wavenumber.conj(),
                characteristic_impedance: characteristic,
                specific_impedance: specific,
                // Pore losses belong to JCA(L), not an additional tube model.
                shear_number: f64::INFINITY,
            },
            reference_impedance,
        )?;
        Ok(Self {
            material,
            response,
            thickness_m,
            area_m2,
            two_port,
        })
    }

    pub const fn material(&self) -> &'a PorousFluidLaw {
        self.material
    }
    pub const fn material_response(&self) -> &PorousFluidResponse {
        &self.response
    }
    pub const fn thickness_m(&self) -> f64 {
        self.thickness_m
    }
    pub const fn area_m2(&self) -> f64 {
        self.area_m2
    }

    /// Existing stable pressure/volume-flow block, suitable for duct chains,
    /// porous multilayers, inviscid air gaps and explicit terminal loads.
    pub const fn two_port(&self) -> TwoPort {
        self.two_port
    }

    /// Evaluate arbitrary coherent incoming pressure waves at both faces.
    /// Each outgoing wave is formed BEFORE its power is evaluated, preserving
    /// the interference term that controls two-sided absorption. A one-sided
    /// absorption coefficient is `absorbed_power_w / incident_power_w` with the
    /// other incident amplitude zero; `1-|reflection|^2` would count transmitted
    /// sound as absorption. This fixed-size query performs no allocations.
    pub fn power_for_incident_pressure(
        &self,
        incident_pressure_pa: [C64; 2],
    ) -> Result<PorousLayerPower, DuctError> {
        if !incident_pressure_pa.into_iter().all(finite) {
            return Err(bad("porous layer incident pressure must be finite"));
        }
        let s = self.two_port.scattering;
        let incident = incident_pressure_pa;
        let outgoing = [
            s[0] * incident[0] + s[1] * incident[1],
            s[2] * incident[0] + s[3] * incident[1],
        ];
        let reference = self.two_port.reference_impedance;
        let pressure = std::array::from_fn(|i| incident[i] + outgoing[i]);
        let flow = std::array::from_fn(|i| (incident[i] - outgoing[i]).scale(1.0 / reference));
        let amplitude_scale = 1.0 / reference.sqrt();
        let incident_power = 0.5
            * incident
                .into_iter()
                .map(|p| p.scale(amplitude_scale).norm_sq())
                .sum::<f64>();
        let outgoing_power = 0.5
            * outgoing
                .into_iter()
                .map(|p| p.scale(amplitude_scale).norm_sq())
                .sum::<f64>();
        let absorbed = incident_power - outgoing_power;
        if !outgoing.into_iter().all(finite)
            || !pressure.into_iter().all(finite)
            || !flow.into_iter().all(finite)
            || !incident_power.is_finite()
            || !outgoing_power.is_finite()
            || !absorbed.is_finite()
        {
            return Err(bad("porous layer field or power is unrepresentable"));
        }
        Ok(PorousLayerPower {
            outgoing_pressure_pa: outgoing,
            pressure_pa: pressure,
            inward_volume_flow_m3_s: flow,
            incident_power_w: incident_power,
            outgoing_power_w: outgoing_power,
            absorbed_power_w: absorbed,
        })
    }
}

#[cfg(test)]
mod tests;
