//! Coherent sound and volume-integrated heat in the original annular section.
//!
//! The native two-port supplies both face fields. Oppositely anchored waves
//! retain interference without growing exponentials; the shared positive wave
//! integrals then integrate the ORIGINAL radial shear/thermal loss forms.
//! No attenuation fit, face-power subtraction, radial remesh or axial sampling
//! replaces those forms. All phasors use exp(-i omega t) and peak amplitudes.
//!
//! Geometry, material state and frequency remain frozen. The complete section's
//! low-reduced-frequency and fully developed assumptions still apply; finite
//! length does not add entrance fields, moving walls or a heat-partition law.

use super::field::{AnnularField, PreparedAnnularFields};
use super::{bad, finite};
use crate::DuctError;
use crate::compose::{TwoPort, TwoPortPower, pressure_wave_integral_weights};
use fs_material::fluid_acoustics::thermoviscous::ThermoviscousFluid;
use fs_material::gas::GasState;
use fs_math::{c64::C64, det};

/// Complete axial pressure/flow at one original cross section.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AnnularAxialField {
    pub pressure_pa: C64,
    /// Positive from the segment's left face toward its right face [m³/s].
    pub volume_flow_m3_s: C64,
    /// Signed cycle-mean acoustic power in that same direction [W].
    pub acoustic_power_w: f64,
}

/// Independently integrated irreversible radial work over a finite volume.
/// This is heat generated in the fluid model, without wall/fluid partition.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct AnnularDissipation {
    pub viscous_power_w: f64,
    pub thermal_power_w: f64,
    pub absorbed_power_w: f64,
}

/// One coherent finite segment, bound to its original radial/material owner.
/// Preparation and every query use fixed work, no allocation or material query.
#[derive(Debug)]
pub struct AnnularExcitation<'f, 'g, 'm, M: ThermoviscousFluid + ?Sized> {
    fields: &'f PreparedAnnularFields<'g, 'm, M>,
    length_m: f64,
    two_port: TwoPort,
    boundary: TwoPortPower,
    forward_left: C64,
    backward_right: C64,
    impedance: C64,
    wavenumber: C64,
}

impl<'g, 'm, M: ThermoviscousFluid + ?Sized> PreparedAnnularFields<'g, 'm, M> {
    /// Excite a positive finite length through the original stable two-port.
    /// The reference is real pressure/volume-flow impedance [Pa s/m³], shared
    /// with the connected network. Both incident pressures combine coherently.
    /// Axial and radial queries reuse this one preparation and its source state.
    pub fn prepare_excitation(
        &self,
        length_m: f64,
        reference_impedance: f64,
        incident_pressure_pa: [C64; 2],
    ) -> Result<AnnularExcitation<'_, 'g, 'm, M>, DuctError> {
        let two_port = self.wave().two_port(length_m, reference_impedance)?;
        let boundary = two_port.power_for_incident_pressure(incident_pressure_pa)?;
        // The average-wave owner's specific impedance and actual area are the
        // same normalization used by TwoPort::from_uniform_wave.
        let impedance =
            self.wave().wave.specific_impedance / C64::from_re(self.section().area_m2());
        let amplitudes = std::array::from_fn::<_, 2, _>(|side| {
            // Face flows point INWARD: this selects F at x=0 and B at x=L.
            (boundary.pressure_pa[side] + impedance * boundary.inward_volume_flow_m3_s[side])
                .scale(0.5)
        });
        if !amplitudes.into_iter().all(finite) {
            return Err(bad("annular internal wave amplitudes are unrepresentable"));
        }
        Ok(AnnularExcitation {
            fields: self,
            length_m,
            two_port,
            boundary,
            forward_left: amplitudes[0],
            backward_right: amplitudes[1],
            impedance,
            wavenumber: self.wave().wave.wavenumber,
        })
    }
}

impl<'f, 'g, 'm, M: ThermoviscousFluid + ?Sized> AnnularExcitation<'f, 'g, 'm, M> {
    pub const fn fields(&self) -> &'f PreparedAnnularFields<'g, 'm, M> {
        self.fields
    }
    pub const fn length_m(&self) -> f64 {
        self.length_m
    }
    pub const fn two_port(&self) -> TwoPort {
        self.two_port
    }
    pub const fn boundary_power(&self) -> &TwoPortPower {
        &self.boundary
    }

    fn admit_depth(&self, x: f64) -> Result<(), DuctError> {
        if !x.is_finite() || x < 0.0 || x > self.length_m {
            Err(bad(
                "annular axial position must lie inside the original segment",
            ))
        } else {
            Ok(())
        }
    }

    /// Complete coherent pressure and volume flow, including standing waves.
    /// An opaque segment still evaluates only decaying exponentials.
    pub fn at_depth(&self, depth_m: f64) -> Result<AnnularAxialField, DuctError> {
        self.admit_depth(depth_m)?;
        let forward = self.forward_left * travel(self.wavenumber, depth_m)?;
        let backward = self.backward_right * travel(self.wavenumber, self.length_m - depth_m)?;
        let pressure_pa = forward + backward;
        let volume_flow_m3_s = (forward - backward) / self.impedance;
        let acoustic_power_w = 0.5 * (pressure_pa * volume_flow_m3_s.conj()).re;
        if !finite(pressure_pa) || !finite(volume_flow_m3_s) || !acoustic_power_w.is_finite() {
            return Err(bad("annular axial field or power is unrepresentable"));
        }
        Ok(AnnularAxialField {
            pressure_pa,
            volume_flow_m3_s,
            acoustic_power_w,
        })
    }

    /// Original radial P1 velocity, thermal response and heat at this depth.
    pub fn field_in_cell(
        &self,
        cell: usize,
        fraction: f64,
        depth_m: f64,
    ) -> Result<AnnularField, DuctError> {
        let axial = self.at_depth(depth_m)?;
        self.fields
            .field_in_cell(cell, fraction, axial.pressure_pa, axial.volume_flow_m3_s)
    }

    /// Exact axial integral over the COMPLETE original annulus, in watts.
    /// Positive forms retain coherent node losses even for subwavelength cells.
    /// This is independent of the signed incident-minus-outgoing diagnostic.
    pub fn dissipation_between(
        &self,
        start_m: f64,
        end_m: f64,
    ) -> Result<AnnularDissipation, DuctError> {
        self.integrate(self.fields.total_loss_roots(), start_m, end_m)
    }

    /// The same finite-length integral in ONE original radial ring.
    /// Its physical volume is fields().cell_area_m2(cell) * (end_m-start_m).
    /// Neither factor of volume nor peak-to-mean factor should be added again
    /// when these watts become a native thermal cell's source.
    pub fn cell_dissipation_between(
        &self,
        cell: usize,
        start_m: f64,
        end_m: f64,
    ) -> Result<AnnularDissipation, DuctError> {
        self.integrate(self.fields.cell_loss_roots(cell)?, start_m, end_m)
    }

    fn integrate(
        &self,
        roots: [C64; 2],
        start_m: f64,
        end_m: f64,
    ) -> Result<AnnularDissipation, DuctError> {
        self.admit_depth(start_m)?;
        self.admit_depth(end_m)?;
        if end_m < start_m {
            return Err(bad("annular heat interval endpoints must be ordered"));
        }
        if start_m == end_m {
            return Ok(AnnularDissipation::default());
        }
        let forward = self.forward_left * travel(self.wavenumber, start_m)?;
        let backward = self.backward_right * travel(self.wavenumber, self.length_m - end_m)?;
        let (plus, minus) = pressure_wave_integral_weights(self.wavenumber, end_m - start_m)?;
        let sum = forward + backward;
        let difference = forward - backward;
        // Positive weights are shared with the native porous owner, including
        // its small-phase series. Apply original radial loss roots BEFORE the
        // squared norm; no new axial quadrature or face subtraction is needed.
        let plus_root = det::sqrt(plus);
        let minus_root = det::sqrt(minus);
        let viscous_power_w = ((difference / self.impedance) * roots[0])
            .scale(plus_root)
            .norm_sq()
            + ((sum / self.impedance) * roots[0])
                .scale(minus_root)
                .norm_sq();
        let thermal_power_w = (sum * roots[1]).scale(plus_root).norm_sq()
            + (difference * roots[1]).scale(minus_root).norm_sq();
        let absorbed_power_w = viscous_power_w + thermal_power_w;
        if ![viscous_power_w, thermal_power_w, absorbed_power_w]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(bad("annular volume-integrated heat is unrepresentable"));
        }
        Ok(AnnularDissipation {
            viscous_power_w,
            thermal_power_w,
            absorbed_power_w,
        })
    }
}

impl AnnularExcitation<'_, '_, '_, GasState> {
    /// Original ideal-gas temperature perturbation in kelvin, at this depth.
    pub fn gas_temperature_in_cell(
        &self,
        cell: usize,
        fraction: f64,
        depth_m: f64,
    ) -> Result<C64, DuctError> {
        let axial = self.at_depth(depth_m)?;
        self.fields
            .gas_temperature_in_cell(cell, fraction, axial.pressure_pa)
    }
}

fn travel(k: C64, distance: f64) -> Result<C64, DuctError> {
    let attenuation = det::exp(-k.im * distance);
    if attenuation == 0.0 {
        return Ok(C64::ZERO);
    }
    let phase = k.re * distance;
    if !phase.is_finite() {
        return Err(bad("annular surviving wave phase is unrepresentable"));
    }
    let (sin, cos) = det::sin_cos(phase);
    Ok(C64::new(cos, sin).scale(attenuation))
}
