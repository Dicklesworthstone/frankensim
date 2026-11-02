//! Interior pressure, flow and heating from the original equivalent-fluid law.
//!
//! In exp(-i omega t), p'=i omega rho_eff v and v'=i omega p/K_eff.
//! Hence -d Re(p v*)/dx / 2 = omega/2 [Im(rho_eff)|v|²
//! + Im(1/K_eff)|p|²]. These are viscous and thermal dissipation per TOTAL
//! specimen volume; porosity is already in the original material coefficients.

use super::{DuctError, PorousLayer, PorousLayerPower, bad, finite};
use fs_math::{c64::C64, det};

/// Cycle-averaged interior fields and local irreversible acoustic losses.
/// Pressure and velocity are peak phasors; positive flow points left to right.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PorousLayerField {
    /// Total acoustic pressure [Pa].
    pub pressure_pa: C64,
    /// Superficial velocity on the full face area, not pore velocity [m/s].
    pub superficial_velocity_m_s: C64,
    /// Axial volume flow [m3/s].
    pub volume_flow_m3_s: C64,
    /// Cycle-mean pore-viscous dissipation per total specimen volume [W/m3].
    pub viscous_dissipation_w_m3: f64,
    /// Cycle-mean thermal-relaxation dissipation per specimen volume [W/m3].
    pub thermal_dissipation_w_m3: f64,
}

/// Exact cell-integrated acoustic dissipation within the homogeneous model.
/// The mechanism split does not assign heat separately to fluid and skeleton.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PorousLayerDissipation {
    /// Pore-viscous dissipation in the requested slab cell [W].
    pub viscous_power_w: f64,
    /// Thermal-relaxation dissipation in the same cell [W].
    pub thermal_power_w: f64,
    /// Sum of viscous and thermal power [W], suitable for a thermal cell source.
    pub absorbed_power_w: f64,
}

/// A coherent excitation prepared once for repeated, allocation-free interior
/// field and thermal-cell queries. The original layer and source remain bound.
#[derive(Debug, Clone)]
pub struct PorousLayerExcitation<'layer, 'material> {
    layer: &'layer PorousLayer<'material>,
    boundary: PorousLayerPower,
    forward_left: C64,
    backward_right: C64,
    impedance: C64,
    wavenumber: C64,
    viscous_weight: f64,
    thermal_weight: f64,
}

impl<'material> PorousLayer<'material> {
    /// Prepare the exact interior solution driven by both incident waves.
    ///
    /// Forward/backward amplitudes are anchored at opposite faces. Only
    /// decaying exponentials are evaluated, including opaque layers. The same
    /// boundary solution as [`Self::power_for_incident_pressure`] is retained;
    /// no material re-query, additional loss law or spatial quadrature occurs.
    pub fn prepare_excitation(
        &self,
        incident_pressure_pa: [C64; 2],
    ) -> Result<PorousLayerExcitation<'_, 'material>, DuctError> {
        let boundary = self.power_for_incident_pressure(incident_pressure_pa)?;
        let (z, k) = self.response.characteristic();
        let impedance = z.conj();
        let amplitudes = std::array::from_fn::<_, 2, _>(|i| {
            // Both face flows point into the slab: at the right face this
            // combination selects the left-travelling internal pressure.
            (boundary.pressure_pa[i]
                + impedance * boundary.inward_volume_flow_m3_s[i].scale(1.0 / self.area_m2))
            .scale(0.5)
        });
        let omega_half = 0.5 * self.response.angular_frequency_rad_s();
        let viscous_weight = -omega_half * self.response.dynamic_density_kg_m3().im;
        let thermal_weight = -omega_half * self.response.dynamic_bulk_modulus_pa().recip().im;
        if !amplitudes.into_iter().all(finite)
            || [viscous_weight, thermal_weight]
                .into_iter()
                .any(|v| !v.is_finite() || v < 0.0)
        {
            return Err(bad(
                "porous interior amplitudes or loss weights are unrepresentable",
            ));
        }
        Ok(PorousLayerExcitation {
            layer: self,
            boundary,
            forward_left: amplitudes[0],
            backward_right: amplitudes[1],
            impedance,
            wavenumber: k.conj(),
            viscous_weight,
            thermal_weight,
        })
    }
}

impl<'layer, 'material> PorousLayerExcitation<'layer, 'material> {
    /// Original source-owned homogeneous layer and coherent fluid snapshot.
    pub const fn layer(&self) -> &'layer PorousLayer<'material> {
        self.layer
    }

    /// The unchanged incoming/outgoing pressure and surface-power balance.
    pub const fn boundary_power(&self) -> &PorousLayerPower {
        &self.boundary
    }

    fn admit_depth(&self, depth: f64) -> Result<(), DuctError> {
        if !depth.is_finite() || depth < 0.0 || depth > self.layer.thickness_m {
            Err(bad("porous field depth must lie inside the physical layer"))
        } else {
            Ok(())
        }
    }

    /// Fields and heating at a depth measured from the left exposed face.
    /// Thermal relaxation and viscous drag retain their original JCA(L) weights.
    /// Their sum is local acoustic power converted to heat, not an independently
    /// qualified fluid/skeleton temperature model.
    pub fn at_depth(&self, depth_m: f64) -> Result<PorousLayerField, DuctError> {
        self.admit_depth(depth_m)?;
        let forward = self.forward_left * travel(self.wavenumber, depth_m)?;
        let backward =
            self.backward_right * travel(self.wavenumber, self.layer.thickness_m - depth_m)?;
        let pressure = forward + backward;
        let velocity = (forward - backward) / self.impedance;
        let flow = velocity.scale(self.layer.area_m2);
        let viscous = self.viscous_weight * velocity.norm_sq();
        let thermal = self.thermal_weight * pressure.norm_sq();
        if ![pressure, velocity, flow].into_iter().all(finite)
            || !viscous.is_finite()
            || !thermal.is_finite()
        {
            return Err(bad("porous interior field or heating is unrepresentable"));
        }
        Ok(PorousLayerField {
            pressure_pa: pressure,
            superficial_velocity_m_s: velocity,
            volume_flow_m3_s: flow,
            viscous_dissipation_w_m3: viscous,
            thermal_dissipation_w_m3: thermal,
        })
    }

    /// Integrate each physical loss mechanism over `[start_m, end_m]`, including
    /// the full exposed area. The coherent standing-wave cross term is retained.
    ///
    /// This is an analytic, constant-work integral, including cells much smaller
    /// than a wavelength or attenuation length. Positive quadratic forms retain
    /// losses near pressure/velocity nodes; subtracting nearly equal face powers
    /// would lose them. Equal endpoints return exact zero. No clipping or
    /// allocation is used. Cell partitions can therefore feed a thermal mesh
    /// without a separate acoustic sampling/convergence axis.
    pub fn dissipation_between(
        &self,
        start_m: f64,
        end_m: f64,
    ) -> Result<PorousLayerDissipation, DuctError> {
        self.admit_depth(start_m)?;
        self.admit_depth(end_m)?;
        if end_m < start_m {
            return Err(bad("porous absorption cell endpoints must be ordered"));
        }
        if start_m == end_m {
            return Ok(PorousLayerDissipation {
                viscous_power_w: 0.0,
                thermal_power_w: 0.0,
                absorbed_power_w: 0.0,
            });
        }
        let forward = self.forward_left * travel(self.wavenumber, start_m)?;
        let backward =
            self.backward_right * travel(self.wavenumber, self.layer.thickness_m - end_m)?;
        let (plus, minus) = integral_weights(self.wavenumber, end_m - start_m)?;
        let sum = forward + backward;
        let difference = forward - backward;
        let pressure_integral = plus * sum.norm_sq() + minus * difference.norm_sq();
        let velocity_integral = minus * (sum / self.impedance).norm_sq()
            + plus * (difference / self.impedance).norm_sq();
        let viscous = self.layer.area_m2 * (self.viscous_weight * velocity_integral);
        let thermal = self.layer.area_m2 * (self.thermal_weight * pressure_integral);
        let absorbed = viscous + thermal;
        if [viscous, thermal, absorbed]
            .into_iter()
            .any(|v| !v.is_finite() || v < 0.0)
        {
            return Err(bad("porous cell dissipation is unrepresentable"));
        }
        Ok(PorousLayerDissipation {
            viscous_power_w: viscous,
            thermal_power_w: thermal,
            absorbed_power_w: absorbed,
        })
    }
}

fn travel(k: C64, distance: f64) -> Result<C64, DuctError> {
    let attenuation = det::exp(-k.im * distance);
    if attenuation == 0.0 {
        return Ok(C64::ZERO);
    }
    let angle = k.re * distance;
    if !angle.is_finite() {
        return Err(bad("porous surviving wave phase is unrepresentable"));
    }
    let (sin, cos) = det::sin_cos(angle);
    Ok(C64::new(cos, sin).scale(attenuation))
}

// For forward amplitude F at the cell's left and backward B at its right,
// integral |p|² = (G+H)/2 |F+B|² + (G-H)/2 |F-B|², where
// G=(1-exp(-2 alpha d))/(2 alpha), H=exp(-alpha d) d sinc(Re(k)d).
// The velocity integral swaps the two nonnegative weights and divides by |Z|².
fn integral_weights(k: C64, length: f64) -> Result<(f64, f64), DuctError> {
    let a = k.im * length;
    let b = k.re * length;
    let decay = det::exp(-a);
    let g = if a == 0.0 {
        length
    } else if a < 0.5 {
        length * (-det::expm1(-2.0 * a) / (2.0 * a))
    } else {
        (-0.5 * det::expm1(-2.0 * a)) / k.im
    };
    let minus = if a <= 0.125 && b <= 0.125 {
        // exp(-a) [sinh(a)/a - sin(b)/b], without subtracting two
        // nearly equal units. Terms through x^10 suffice in this small domain.
        let aa = a * a;
        let bb = b * b;
        let hyperbolic = aa
            * (1.0 / 6.0
                + aa * (1.0 / 120.0
                    + aa * (1.0 / 5040.0 + aa * (1.0 / 362880.0 + aa / 39916800.0))));
        let circular = bb
            * (1.0 / 6.0
                - bb * (1.0 / 120.0
                    - bb * (1.0 / 5040.0 - bb * (1.0 / 362880.0 - bb / 39916800.0))));
        length * (0.5 * decay * (hyperbolic + circular))
    } else {
        let h = if decay == 0.0 {
            0.0
        } else if b == 0.0 {
            length * decay
        } else if !b.is_finite() {
            return Err(bad("porous cell interference phase is unrepresentable"));
        } else {
            (length * (det::sin(b) / b)) * decay
        };
        0.5 * g - 0.5 * h
    };
    let plus = g - minus;
    if [plus, minus].into_iter().any(|v| !v.is_finite() || v < 0.0) {
        return Err(bad("porous cell loss integral is unrepresentable"));
    }
    Ok((plus, minus))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn g1_integrals_retain_lossless_and_subroundoff_node_limits() {
        assert_eq!(integral_weights(C64::ZERO, 0.3).unwrap(), (0.3, 0.0));
        for alpha in [0.0, 0.7] {
            let k = C64::new(3.0, alpha);
            let d = 1e-9;
            let (plus, minus) = integral_weights(k, d).unwrap();
            let leading = (k.re * k.re + k.im * k.im) * d * d * d / 12.0;
            assert!((minus / leading - 1.0).abs() < 1e-9);
            assert!((plus / d - 1.0).abs() < 1e-9);
            assert!(minus > 0.0);
        }
    }
}
