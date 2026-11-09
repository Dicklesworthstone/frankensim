//! Source-bound solid jackets in the original rectangular acoustic energy equation.
//!
//! Nested rectangular material layers and the gas share one conforming Q1
//! temperature field. Full in-plane conductivity tensors retain transverse,
//! lateral and corner conduction; no local-normal or isothermal-solid reduction
//! is made. Contacts are perfect. The four outer faces have independent baths.
//! The existing gas molecular Cp, no-slip shear and stable acoustic TwoPort are
//! reused. All fields are peak exp(-i omega t) phasors about the same uniform T0.
//!
//! Sources explicitly declare their in-plane material axes and caloric
//! constraint. This is a frozen linear tangent, not phase change, thermoelastic
//! motion, contact resistance, mean heating evolution or axial heat diffusion.

use super::{
    DiffusionField, RectangularAdmission, RectangularSection, RectangularThermalField,
    RectangularThermalLimits, RectangularThermalWallError as Error, admit_tail, amplitudes, bad,
    finite, poll, positive,
};
use crate::{DuctError, SegmentWave, compose::TwoPort};
use fs_material::{
    fluid_acoustics::thermoviscous::{ThermoviscousFluid, ThermoviscousFluidState},
    gas::GasState,
    relaxation::RelaxationState,
    thermal::transport::orthotropic::{
        OrthotropicThermalTransportSample, PreparedOrthotropicThermalTransport,
    },
};
use fs_math::{c64::C64, det};

mod physics;

/// One complete rectangular shell, ordered from the gas toward the exterior.
/// Every corner belongs to the outermost layer reached along either axis.
/// The declared material basis uses [short, long] gas coordinates.
#[derive(Clone, Copy, Debug)]
pub struct RectangularSolidLayer<'a> {
    pub source: &'a PreparedOrthotropicThermalTransport<2>,
    /// Positive increments on [short zero/end, long zero/end] faces [m].
    pub thickness_m: [f64; 4],
    /// Positive cells through each corresponding increment. Gas cells are
    /// specified separately by RectangularThermalLimits::cells.
    pub cells: [usize; 4],
}

/// Outer-face condition at the original common background temperature.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum RectangularJacketBath {
    Isothermal,
    Adiabatic,
    /// Per actual outer surface area; no gas-film conductance is inferred.
    Conductance {
        conductance_w_m2_k: f64,
    },
}
impl RectangularJacketBath {
    fn conductance(self) -> Result<f64, Error> {
        let g = match self {
            Self::Conductance { conductance_w_m2_k } => conductance_w_m2_k,
            _ => 0.0,
        };
        if g.is_finite() && g >= 0.0 {
            Ok(g)
        } else {
            Err(bad("rectangular jacket bath requires finite nonnegative conductance").into())
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RectangularSolidThermalField {
    pub layer: usize,
    pub coordinates_m: [f64; 2],
    pub temperature_k: C64,
    pub temperature_gradient_k_m: [C64; 2],
    /// Cartesian components of the original -K grad(T) [W/m2].
    pub heat_flux_w_m2: [C64; 2],
    pub conduction_dissipation_w_m3: f64,
    pub thermal_storage_j_m3: f64,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RectangularJacketBoundaryField {
    pub coordinates_m: [f64; 2],
    pub temperature_k: C64,
    /// Original G*T law, or one-sided Q1 tensor flux for an isothermal face.
    pub wall_to_bath_heat_flux_w_m2: C64,
    /// Independent one-sided Q1 tensor flux; refine for a local flux claim.
    pub reconstructed_heat_flux_w_m2: C64,
    pub bath_dissipation_w_m2: f64,
}
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct RectangularConjugateDissipation {
    pub viscous_w_m: f64,
    pub gas_conduction_w_m: f64,
    pub molecular_w_m: f64,
    pub solid_conduction_w_m: f64,
    pub baths_w_m: [f64; 4],
    pub total_w_m: f64,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RectangularSolidThermalObservation {
    pub conduction_dissipation_w_m: f64,
    pub thermal_storage_j_m: f64,
    /// Original C_vol integral(T') over this layer [J/m], not quadratic storage.
    pub oscillatory_thermal_energy_j_m: C64,
}

/// Immutable source-bound gas and all original solid cells. Preparation uses
/// the native bounded sparse factor. Point queries use two axis binary searches
/// and one Q1 cell; integrated observations are O(1), allocate nothing and
/// perform no material query. Shear points retain the original Fourier cost.
#[derive(Debug)]
pub struct PreparedRectangularConduction<'g, 'm, 's, M: ?Sized> {
    section: &'g RectangularSection,
    source: &'m M,
    layers: Vec<RectangularSolidLayer<'s>>,
    samples: Vec<OrthotropicThermalTransportSample<2>>,
    state: ThermoviscousFluidState,
    omega: f64,
    temperature: f64,
    baths: [RectangularJacketBath; 4],
    wave: SegmentWave,
    viscous: DiffusionField,
    velocity_scale: C64,
    thermal: physics::JacketField,
    loss_roots: [f64; 8],
    molecular_root: f64,
    layer_loss_roots: Vec<f64>,
    layer_storage_roots: Vec<f64>,
    relative_port_work_defect: f64,
    peak_payload_bytes: usize,
}

impl RectangularSection {
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_conjugate_jacket_fields<'g, 'm, 's>(
        &'g self,
        gas: &'m GasState,
        omega: f64,
        layers: &[RectangularSolidLayer<'s>],
        baths: [RectangularJacketBath; 4],
        admission: RectangularAdmission,
        limits: RectangularThermalLimits,
    ) -> Result<PreparedRectangularConduction<'g, 'm, 's, GasState>, Error> {
        self.prepare_conjugate_jacket_fields_controlled(
            gas,
            omega,
            layers,
            baths,
            admission,
            limits,
            || false,
        )
    }
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_conjugate_jacket_fields_controlled<'g, 'm, 's>(
        &'g self,
        gas: &'m GasState,
        omega: f64,
        layers: &[RectangularSolidLayer<'s>],
        baths: [RectangularJacketBath; 4],
        admission: RectangularAdmission,
        limits: RectangularThermalLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedRectangularConduction<'g, 'm, 's, GasState>, Error> {
        poll(&mut cancelled, 0)?;
        let state = gas.thermoviscous_state()?;
        prepare(
            self,
            gas,
            gas,
            state,
            C64::from_re(state.transport().specific_heat_j_kg_k()),
            omega,
            layers,
            baths,
            admission,
            limits,
            &mut cancelled,
        )
    }
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_molecular_conjugate_jacket_fields<'g, 'm, 's>(
        &'g self,
        material: &'m RelaxationState,
        omega: f64,
        layers: &[RectangularSolidLayer<'s>],
        baths: [RectangularJacketBath; 4],
        admission: RectangularAdmission,
        limits: RectangularThermalLimits,
    ) -> Result<PreparedRectangularConduction<'g, 'm, 's, RelaxationState>, Error> {
        self.prepare_molecular_conjugate_jacket_fields_controlled(
            material,
            omega,
            layers,
            baths,
            admission,
            limits,
            || false,
        )
    }
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_molecular_conjugate_jacket_fields_controlled<'g, 'm, 's>(
        &'g self,
        material: &'m RelaxationState,
        omega: f64,
        layers: &[RectangularSolidLayer<'s>],
        baths: [RectangularJacketBath; 4],
        admission: RectangularAdmission,
        limits: RectangularThermalLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedRectangularConduction<'g, 'm, 's, RelaxationState>, Error> {
        poll(&mut cancelled, 0)?;
        let gas = material.gas();
        let state = gas.thermoviscous_state()?;
        let cp = if material.processes().is_empty() {
            C64::from_re(state.transport().specific_heat_j_kg_k())
        } else {
            material.caloric_response(omega)?.cp_j_kg_k
        };
        prepare(
            self,
            material,
            gas,
            state,
            cp,
            omega,
            layers,
            baths,
            admission,
            limits,
            &mut cancelled,
        )
    }
}

#[allow(clippy::too_many_arguments)]
fn prepare<'g, 'm, 's, M: ?Sized>(
    section: &'g RectangularSection,
    source: &'m M,
    gas: &GasState,
    state: ThermoviscousFluidState,
    cp: C64,
    omega: f64,
    layers: &[RectangularSolidLayer<'s>],
    baths: [RectangularJacketBath; 4],
    admission: RectangularAdmission,
    limits: RectangularThermalLimits,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<PreparedRectangularConduction<'g, 'm, 's, M>, Error> {
    if layers.is_empty()
        || layers.len() > 8
        || !positive(omega)
        || !positive(admission.maximum_relative_tail)
        || admission.maximum_relative_tail >= 1.0
        || !positive(admission.maximum_transverse_phase)
        || admission.maximum_transverse_phase >= core::f64::consts::PI
        || !finite(cp)
        || cp.re <= 0.0
        || cp.im < 0.0
    {
        return Err(bad("invalid original rectangular conjugate source or admission").into());
    }
    let k0 = omega / state.acoustic().sound_speed_m_s();
    if !positive(k0) || k0 * section.long_m > admission.maximum_transverse_phase {
        return Err(bad("conjugate gas exceeds its original plane-pressure envelope").into());
    }
    let fourier_bytes = section
        .terms()
        .checked_mul(128)
        .and_then(|v| v.checked_add(8192))
        .ok_or_else(|| bad("conjugate shear payload extent overflow"))?;
    let thermal_bytes = limits
        .maximum_payload_bytes
        .checked_sub(fourier_bytes)
        .ok_or_else(|| bad("conjugate payload cannot retain its original shear field"))?;
    let mut samples = Vec::with_capacity(layers.len());
    let mut numerical = Vec::with_capacity(layers.len());
    for (index, layer) in layers.iter().enumerate() {
        poll(cancelled, index)?;
        let sample = layer.source.at_temperature(gas.temperature)?;
        numerical.push(physics::Layer {
            thickness: layer.thickness_m,
            cells: layer.cells,
            conductivity: sample.conductivity_w_m_k(),
            capacity: sample.density_kg_m3() * sample.specific_heat_j_kg_k(),
        });
        samples.push(sample);
    }
    let transport = state.transport();
    let thermal = physics::prepare(
        [section.short_m, section.long_m],
        transport.thermal_conductivity_w_m_k(),
        cp.conj().scale(gas.density),
        omega,
        &numerical,
        baths,
        RectangularThermalLimits {
            maximum_payload_bytes: thermal_bytes,
            ..limits
        },
        cancelled,
    )?;
    let diffusion = section.diffusion(omega, transport.kinematic_viscosity_m2_s(), cancelled)?;
    admit_tail(diffusion, admission.maximum_relative_tail)?;
    let viscous = section.prepare_diffusion_fields(
        C64::new(0.0, omega / transport.kinematic_viscosity_m2_s()),
        cancelled,
    )?;
    let velocity_scale = diffusion.mean.conj().recip().scale(1.0 / section.area_m2);
    let molecular_root = det::sqrt(omega * gas.density * cp.im / (2.0 * gas.temperature));
    let scale = 1.0 / (2.0 * gas.temperature);
    let mut loss_roots = [0.0; 8];
    loss_roots[0] = velocity_scale.abs()
        * det::sqrt(transport.dynamic_viscosity_pa_s() * 0.5)
        * viscous.gradient_norm;
    loss_roots[1] = det::sqrt(thermal.forms[0].conduction * scale);
    loss_roots[2] = molecular_root * det::sqrt(thermal.forms[0].field);
    let layer_loss_roots: Vec<_> = thermal.forms[1..]
        .iter()
        .map(|f| det::sqrt(f.conduction * scale))
        .collect();
    let layer_storage_roots: Vec<_> = thermal.forms[1..]
        .iter()
        .zip(&numerical)
        .map(|(f, layer)| det::sqrt(layer.capacity * f.field * scale * 0.5))
        .collect();
    loss_roots[3] = det::sqrt(thermal.forms[1..].iter().map(|f| f.conduction).sum::<f64>() * scale);
    for side in 0..4 {
        loss_roots[4 + side] =
            det::sqrt(baths[side].conductance()? * thermal.bath_forms[side] * scale);
    }
    let compliance =
        C64::from_re(1.0 / (gas.density * gas.specific_gas_constant * gas.temperature))
            - thermal.mean_gas.scale(1.0 / gas.temperature);
    let s = C64::new(0.0, omega);
    let series = s.scale(gas.density / section.area_m2) / diffusion.mean;
    let shunt = s.scale(section.area_m2) * compliance;
    let (wavenumber, characteristic_impedance) =
        crate::helmholtz_wave_from_telegraph(series, shunt);
    let wave = SegmentWave {
        wavenumber,
        characteristic_impedance,
        specific_impedance: characteristic_impedance.scale(section.area_m2),
        shear_number: 0.5
            * section.short_m
            * det::sqrt(omega / transport.kinematic_viscosity_m2_s()),
    };
    let heat = loss_roots[1..].iter().map(|v| v * v).sum::<f64>();
    let port_heat = 0.5 * shunt.re;
    let relative_port_work_defect = if heat == 0.0 && port_heat == 0.0 {
        0.0
    } else {
        (heat - port_heat).abs() / (heat + port_heat.abs()).max(f64::MIN_POSITIVE)
    };
    if ![
        compliance,
        wavenumber,
        characteristic_impedance,
        wave.specific_impedance,
        velocity_scale,
    ]
    .into_iter()
    .all(finite)
        || compliance.re <= 0.0
        || wavenumber.re <= 0.0
        || wavenumber.im < 0.0
        || characteristic_impedance.re <= 0.0
        || k0.max(wavenumber.abs()) * section.long_m > admission.maximum_transverse_phase
        || !relative_port_work_defect.is_finite()
        || relative_port_work_defect > limits.maximum_relative_residual
        || loss_roots
            .into_iter()
            .chain(layer_storage_roots.iter().copied())
            .any(|v| !v.is_finite() || v < 0.0)
    {
        return Err(bad("original conjugate wave or physical port work exceeds admission").into());
    }
    let peak_payload_bytes = thermal.peak_payload_bytes + fourier_bytes;
    poll(cancelled, section.terms())?;
    Ok(PreparedRectangularConduction {
        section,
        source,
        layers: layers.to_vec(),
        samples,
        state,
        omega,
        temperature: gas.temperature,
        baths,
        wave,
        viscous,
        velocity_scale,
        thermal,
        loss_roots,
        molecular_root,
        layer_loss_roots,
        layer_storage_roots,
        relative_port_work_defect,
        peak_payload_bytes,
    })
}

impl<'g, 'm, 's, M: ?Sized> PreparedRectangularConduction<'g, 'm, 's, M> {
    pub const fn source(&self) -> &'m M {
        self.source
    }
    pub const fn section(&self) -> &'g RectangularSection {
        self.section
    }
    pub fn layers(&self) -> &[RectangularSolidLayer<'s>] {
        &self.layers
    }
    pub fn solid_samples(&self) -> &[OrthotropicThermalTransportSample<2>] {
        &self.samples
    }
    pub const fn state(&self) -> &ThermoviscousFluidState {
        &self.state
    }
    pub const fn angular_frequency_rad_s(&self) -> f64 {
        self.omega
    }
    pub const fn background_temperature_k(&self) -> f64 {
        self.temperature
    }
    pub const fn wave(&self) -> &SegmentWave {
        &self.wave
    }
    pub fn mean_gas_temperature_per_pressure_k_pa(&self) -> C64 {
        self.thermal.mean_gas.conj()
    }
    pub fn thermal_axes_m(&self) -> [&[f64]; 2] {
        [&self.thermal.axes[0], &self.thermal.axes[1]]
    }
    pub fn relative_backward_residual(&self) -> f64 {
        self.thermal.relative_residual
    }
    pub fn relative_heat_balance_defect(&self) -> f64 {
        self.thermal.relative_heat_defect
    }
    pub const fn relative_thermal_work_defect(&self) -> f64 {
        self.relative_port_work_defect
    }
    pub const fn peak_payload_bytes(&self) -> usize {
        self.peak_payload_bytes
    }
    pub fn thermal_factor_nonzeros(&self) -> usize {
        self.thermal.factor_nonzeros
    }
    pub fn two_port(&self, length_m: f64, reference: f64) -> Result<TwoPort, DuctError> {
        TwoPort::from_uniform_wave(&self.wave, self.section.area_m2, length_m, reference)
    }
    pub fn gas_field_at(
        &self,
        coordinates_m: [f64; 2],
        pressure_pa: C64,
        volume_flow_m3_s: C64,
    ) -> Result<RectangularThermalField, Error> {
        self.gas_field_at_controlled(coordinates_m, pressure_pa, volume_flow_m3_s, || false)
    }
    pub fn gas_field_at_controlled(
        &self,
        coordinates_m: [f64; 2],
        pressure_pa: C64,
        volume_flow_m3_s: C64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RectangularThermalField, Error> {
        poll(&mut cancelled, 0)?;
        amplitudes(pressure_pa, volume_flow_m3_s)?;
        let t = self
            .thermal
            .at(0, coordinates_m)?
            .map(|v| v.conj() * pressure_pa);
        let v = self
            .viscous
            .at(
                coordinates_m[0] / self.section.short_m,
                coordinates_m[1],
                self.section.long_m,
                &mut cancelled,
            )?
            .map(|v| v.conj() * self.velocity_scale * volume_flow_m3_s);
        let transport = self.state.transport();
        let viscous = (v[1].norm_sq() + v[2].norm_sq()) * transport.dynamic_viscosity_pa_s() * 0.5;
        let conduction = (t[1].norm_sq() + t[2].norm_sq()) * transport.thermal_conductivity_w_m_k()
            / (2.0 * self.temperature);
        let molecular = t[0].scale(self.molecular_root).norm_sq();
        if !v.into_iter().chain(t).all(finite)
            || [viscous, conduction, molecular]
                .into_iter()
                .any(|v| !v.is_finite())
        {
            return Err(
                bad("original conjugate gas field or local work is unrepresentable").into(),
            );
        }
        Ok(RectangularThermalField {
            coordinates_m,
            pressure_pa,
            axial_velocity_m_s: v[0],
            transverse_velocity_gradient_per_s: [v[1], v[2]],
            gas_temperature_k: t[0],
            transverse_temperature_gradient_k_m: [t[1], t[2]],
            viscous_dissipation_w_m3: viscous,
            conduction_dissipation_w_m3: conduction,
            molecular_dissipation_w_m3: molecular,
        })
    }
    /// Explicit layer selects the physical side of a material interface.
    pub fn solid_field_at(
        &self,
        layer: usize,
        coordinates_m: [f64; 2],
        pressure_pa: C64,
    ) -> Result<RectangularSolidThermalField, Error> {
        amplitudes(pressure_pa, C64::ZERO)?;
        let sample = self
            .samples
            .get(layer)
            .ok_or_else(|| bad("original rectangular solid layer required"))?;
        let t = self
            .thermal
            .at(layer + 1, coordinates_m)?
            .map(|v| v.conj() * pressure_pa);
        let k = sample.conductivity_w_m_k();
        let gradient = [t[1], t[2]];
        let flux =
            core::array::from_fn(|i| -(gradient[0].scale(k[i][0]) + gradient[1].scale(k[i][1])));
        let conduction = physics::quadratic(k, gradient)? / (2.0 * self.temperature);
        let storage = sample.density_kg_m3() * sample.specific_heat_j_kg_k()
            / (4.0 * self.temperature)
            * t[0].norm_sq();
        if !t.into_iter().chain(flux).all(finite) || !conduction.is_finite() || !storage.is_finite()
        {
            return Err(bad("original solid field, heat or storage is unrepresentable").into());
        }
        Ok(RectangularSolidThermalField {
            layer,
            coordinates_m,
            temperature_k: t[0],
            temperature_gradient_k_m: gradient,
            heat_flux_w_m2: flux,
            conduction_dissipation_w_m3: conduction,
            thermal_storage_j_m3: storage,
        })
    }
    pub fn outer_wall_observation(
        &self,
        side: usize,
        fraction: f64,
        pressure_pa: C64,
    ) -> Result<RectangularJacketBoundaryField, Error> {
        if side >= 4 || !fraction.is_finite() || !(0.0..=1.0).contains(&fraction) {
            return Err(bad("original outer face and fraction in [0,1] required").into());
        }
        let extent: [[f64; 2]; 2] = core::array::from_fn(|d| {
            [
                self.thermal.axes[d][0],
                *self.thermal.axes[d].last().unwrap(),
            ]
        });
        let axis = usize::from(side >= 2);
        let mut point = [0.0; 2];
        point[axis] = extent[axis][side % 2];
        point[1 - axis] = if fraction == 1.0 {
            extent[1 - axis][1]
        } else {
            extent[1 - axis][0] + fraction * (extent[1 - axis][1] - extent[1 - axis][0])
        };
        let field = self.solid_field_at(self.layers.len() - 1, point, pressure_pa)?;
        let reconstructed =
            field.heat_flux_w_m2[axis].scale(if side % 2 == 0 { -1.0 } else { 1.0 });
        let g = self.baths[side].conductance()?;
        let flux = if self.baths[side] == RectangularJacketBath::Isothermal {
            reconstructed
        } else {
            field.temperature_k.scale(g)
        };
        let loss = field
            .temperature_k
            .scale(det::sqrt(g / (2.0 * self.temperature)))
            .norm_sq();
        if !finite(flux) || !loss.is_finite() {
            return Err(bad("outer jacket bath field is unrepresentable").into());
        }
        Ok(RectangularJacketBoundaryField {
            coordinates_m: point,
            temperature_k: field.temperature_k,
            wall_to_bath_heat_flux_w_m2: flux,
            reconstructed_heat_flux_w_m2: reconstructed,
            bath_dissipation_w_m2: loss,
        })
    }
    pub fn layer_observation(
        &self,
        layer: usize,
        pressure_pa: C64,
    ) -> Result<RectangularSolidThermalObservation, Error> {
        amplitudes(pressure_pa, C64::ZERO)?;
        let sample = self
            .samples
            .get(layer)
            .ok_or_else(|| bad("original rectangular solid layer required"))?;
        let loss = pressure_pa.scale(self.layer_loss_roots[layer]).norm_sq();
        let storage = pressure_pa.scale(self.layer_storage_roots[layer]).norm_sq();
        let energy = self.thermal.forms[layer + 1].temperature_integral.conj()
            * pressure_pa
            * C64::from_re(sample.density_kg_m3() * sample.specific_heat_j_kg_k());
        if !loss.is_finite() || !storage.is_finite() || !finite(energy) {
            return Err(bad("solid integrated heat or storage is unrepresentable").into());
        }
        Ok(RectangularSolidThermalObservation {
            conduction_dissipation_w_m: loss,
            thermal_storage_j_m: storage,
            oscillatory_thermal_energy_j_m: energy,
        })
    }
    /// Total signed oscillatory heat to the four baths per axial metre. Strong
    /// isothermal faces use their summed ORIGINAL weak reaction; no corner
    /// reaction is arbitrarily distributed among pointwise face observations.
    pub fn total_wall_to_bath_heat_flux_w_m(&self, pressure_pa: C64) -> Result<C64, Error> {
        amplitudes(pressure_pa, C64::ZERO)?;
        let flux = self.thermal.bath_flux.conj() * pressure_pa;
        if !finite(flux) {
            return Err(bad("total original bath heat flux is unrepresentable").into());
        }
        Ok(flux)
    }
    pub fn dissipation_per_axial_length(
        &self,
        pressure_pa: C64,
        volume_flow_m3_s: C64,
    ) -> Result<RectangularConjugateDissipation, Error> {
        amplitudes(pressure_pa, volume_flow_m3_s)?;
        let values = core::array::from_fn::<_, 8, _>(|i| {
            (if i == 0 {
                volume_flow_m3_s
            } else {
                pressure_pa
            })
            .scale(self.loss_roots[i])
            .norm_sq()
        });
        let total = values.iter().sum::<f64>();
        if values.into_iter().any(|v| !v.is_finite()) || !total.is_finite() {
            return Err(bad("conjugate integrated physical work is unrepresentable").into());
        }
        Ok(RectangularConjugateDissipation {
            viscous_w_m: values[0],
            gas_conduction_w_m: values[1],
            molecular_w_m: values[2],
            solid_conduction_w_m: values[3],
            baths_w_m: [values[4], values[5], values[6], values[7]],
            total_w_m: total,
        })
    }
}

#[cfg(test)]
mod tests;
