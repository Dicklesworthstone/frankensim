//! Ideal thermal/mechanical contacts between original homogeneous NSF regions.
//! All four fields cross every interior interface; only the two external ends
//! are thermally terminated. Reference junctions change power coordinates,
//! never coefficients or thermodynamic pressure into total traction.

use super::{
    BulkThermalCylinder, BulkThermalError, BulkThermoviscousEquilibrium, ComposeOptions, DuctError,
    Segment, ThermalEndCondition, ThermalNetwork, bad, redheffer_star,
};
use fs_math::{c64::C64, det};

/// One authored homogeneous region, in left-to-right physical order. The
/// material snapshot retains independent shear/bulk viscosity and conductivity.
#[derive(Clone, Copy, Debug)]
pub struct BulkThermalRegion {
    pub material: BulkThermoviscousEquilibrium,
    /// Positive finite physical axial length [m].
    pub length_m: f64,
}

impl BulkThermalCylinder {
    /// Compose full NSF regions at one constant radius and frequency. Interior
    /// contacts preserve total traction, volume flow, temperature perturbation
    /// and heat power. Each region must have exactly the same reference T0;
    /// thermal conditions are applied ONCE, at the two external ends.
    ///
    /// `maximum_blocks` bounds authored regions; `maximum_condition` bounds the
    /// existing local and Redheffer solves. The first region's positive real
    /// references are used for both external ports. Separate endpoint material
    /// maps reconstruct pressure, which need not be continuous across different
    /// longitudinal viscosities. No material or loss coefficients are averaged.
    ///
    /// These are ideal zero-thickness contacts in a frozen equilibrium. The
    /// caller establishes mechanical/phase compatibility. No contact thermal
    /// resistance, capillarity, phase transfer, mixing or evolving interface is
    /// supplied; the original slip-wall, constant-area NSF scope is retained.
    pub fn from_regions(
        regions: &[BulkThermalRegion],
        radius_m: f64,
        omega_rad_s: f64,
        ends: [ThermalEndCondition; 2],
        options: ComposeOptions,
    ) -> Result<Self, BulkThermalError> {
        Self::from_regions_controlled(regions, radius_m, omega_rad_s, ends, options, || false)
    }

    /// Same complete composition with cancellation between bounded region
    /// preparations and inside the native reference/junction solves. A refusal
    /// publishes no partially connected or thermally terminated cylinder.
    pub fn from_regions_controlled(
        regions: &[BulkThermalRegion],
        radius_m: f64,
        omega_rad_s: f64,
        ends: [ThermalEndCondition; 2],
        options: ComposeOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, BulkThermalError> {
        options.validate()?;
        let first = regions.first().ok_or(DuctError::EmptyDuct)?;
        if regions.len() > options.maximum_blocks {
            return Err(bad("regional thermal cylinder exceeds its declared region cap").into());
        }
        let temperature = first.material.temperature_k();
        for (index, region) in regions.iter().enumerate() {
            poll(&mut cancelled, index)?;
            Segment::Cylinder {
                radius: radius_m,
                length: region.length_m,
            }
            .validate()?;
            if region.material.temperature_k() != temperature {
                return Err(bad(
                    "regional thermal contacts require one common equilibrium temperature",
                )
                .into());
            }
        }
        let mut network = ThermalNetwork::new(
            &first.material,
            radius_m,
            first.length_m,
            omega_rad_s,
            options.maximum_condition,
        )?;
        for (index, region) in regions.iter().enumerate().skip(1) {
            poll(&mut cancelled, index)?;
            let right = ThermalNetwork::new(
                &region.material,
                radius_m,
                region.length_m,
                omega_rad_s,
                options.maximum_condition,
            )?;
            network.append(right, options.maximum_condition, &mut cancelled)?;
        }
        poll(&mut cancelled, regions.len())?;
        Ok(Self::from_network(
            network,
            ends,
            options.maximum_condition,
        )?)
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool, completed_steps: usize) -> Result<(), DuctError> {
    if cancelled() {
        Err(DuctError::Cancelled { completed_steps })
    } else {
        Ok(())
    }
}

impl ThermalNetwork {
    // Shared by classical and molecular regions. Join BEFORE thermal end
    // elimination: resetting local heat or temperature would remove the
    // thermal contact mode and change the physical compliance and loss.
    pub(super) fn append(
        &mut self,
        mut right: Self,
        limit: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), DuctError> {
        if self.temperature_k != right.temperature_k {
            return Err(bad(
                "regional thermal contacts require one common equilibrium temperature",
            ));
        }
        right.rebase(self.references, limit, cancelled)?;
        let (scattering, condition) =
            redheffer_star(&self.scattering, &right.scattering, 2, limit, cancelled)?;
        self.scattering.copy_from_slice(&scattering);
        self.endpoints[1] = right.endpoints[1];
        self.condition = self.condition.max(right.condition).max(condition);
        Ok(())
    }

    // Change both EXTERNAL power references through lossless physical-field
    // continuity. Keep the original endpoint pressure maps in THEIR references.
    fn rebase(
        &mut self,
        references: [f64; 2],
        limit: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), DuctError> {
        if references == self.references {
            return Ok(());
        }
        let left = reference_junction(references, self.references)?;
        let right = reference_junction(self.references, references)?;
        let (first, c1) = redheffer_star(&left, &self.scattering, 2, limit, &mut *cancelled)?;
        let (result, c2) = redheffer_star(&first, &right, 2, limit, cancelled)?;
        self.scattering.copy_from_slice(&result);
        self.references = references;
        self.condition = self.condition.max(c1).max(c2);
        Ok(())
    }
}

// For each conjugate pair, e_left=e_right and j_left=-j_right. With real
// references Ra,Rb, r=(Rb-Ra)/(Ra+Rb), t=2sqrt(RaRb)/(Ra+Rb). Scaled roots
// avoid overflowing Ra+Rb or Ra*Rb. This junction is reciprocal and lossless.
fn reference_junction(left: [f64; 2], right: [f64; 2]) -> Result<[C64; 16], DuctError> {
    if left
        .iter()
        .chain(&right)
        .any(|r| !r.is_finite() || *r <= 0.0)
    {
        return Err(bad(
            "thermal junction references must be finite and positive",
        ));
    }
    let mut s = [C64::ZERO; 16];
    for i in 0..2 {
        let a = det::sqrt(left[i]);
        let b = det::sqrt(right[i]);
        let ratio = a.min(b) / a.max(b);
        let denominator = 1.0 + ratio * ratio;
        let reflection = ((1.0 - ratio) * (1.0 + ratio) / denominator)
            * if right[i] >= left[i] { 1.0 } else { -1.0 };
        let transmission = 2.0 * ratio / denominator;
        s[i * 4 + i] = C64::from_re(reflection);
        s[i * 4 + i + 2] = C64::from_re(transmission);
        s[(i + 2) * 4 + i] = C64::from_re(transmission);
        s[(i + 2) * 4 + i + 2] = C64::from_re(-reflection);
    }
    Ok(s)
}

#[cfg(test)]
mod tests {
    use super::*;
    use ThermalEndCondition::{Adiabatic, Isothermal};
    use fs_material::fluid_acoustics::{
        FluidAcousticState, thermoviscous::ThermoviscousFluidState,
    };
    use fs_material::fluid_transport::{FluidTransportState, ThermalFluidTransportState};

    const RADIUS: f64 = 0.02;
    fn material(alpha: f64, scale: f64, temperature: f64) -> BulkThermoviscousEquilibrium {
        let (rho, c, cp) = (2.0 * scale, 10.0 * det::sqrt(scale), 8.0);
        let transport =
            FluidTransportState::try_new(rho, cp, 0.6 * scale, 0.003 * scale * scale).unwrap();
        let beta_s = 1.0 / (rho * c * c);
        let beta_t = beta_s + alpha * alpha * temperature / (rho * cp);
        let fluid = ThermoviscousFluidState::try_new(
            FluidAcousticState::try_new(rho, c).unwrap(),
            transport,
            cp * beta_s / beta_t,
        )
        .unwrap();
        BulkThermoviscousEquilibrium::from_states(
            fluid,
            ThermalFluidTransportState::try_new(transport, temperature, Some(alpha)).unwrap(),
            0.002 / scale,
        )
        .unwrap()
    }
    fn region(alpha: f64, scale: f64, length_m: f64) -> BulkThermalRegion {
        BulkThermalRegion {
            material: material(alpha, scale, 12.5),
            length_m,
        }
    }
    fn near(a: C64, b: C64, relative: f64) {
        assert!(
            (a - b).abs() <= relative * b.abs().max(1e-100),
            "{a:?} != {b:?}"
        );
    }
    fn prepare(
        regions: &[BulkThermalRegion],
        omega: f64,
        ends: [ThermalEndCondition; 2],
    ) -> BulkThermalCylinder {
        BulkThermalCylinder::from_regions(regions, RADIUS, omega, ends, ComposeOptions::default())
            .unwrap()
    }

    #[test]
    fn g1_regional_adiabatic_compliance_retains_internal_heat_exchange() {
        let area = core::f64::consts::PI * RADIUS * RADIUS;
        let omega = 1e-5;
        for alpha in [0.08, -0.08, -0.04] {
            let regions = [region(0.08, 1.0, 0.04), region(alpha, 1.0, 0.04)];
            let (mut isothermal, mut expansion, mut capacity, mut reset) = (0.0, 0.0, 0.0, 0.0);
            for region in &regions {
                let volume = area * region.length_m;
                let state = region.material.absorption().state();
                isothermal += region.material.isothermal_compressibility_pa_inv() * volume;
                expansion += region.material.isobaric_expansion_per_k() * volume;
                capacity += state.acoustic().density_kg_m3()
                    * state.transport().specific_heat_j_kg_k()
                    * volume;
                reset += state.acoustic().adiabatic_compressibility_per_pa() * volume;
            }
            // Uniform leading pressure/temperature and zero TOTAL entropy
            // flux give this compliance. Local entropy is NOT held fixed:
            // the regions exchange heat through their internal contact.
            let adiabatic = isothermal - 12.5 * expansion * expansion / capacity;
            for ends in [[Adiabatic, Adiabatic], [Isothermal, Adiabatic]] {
                let response = prepare(&regions, omega, ends)
                    .response(None, ComposeOptions::default())
                    .unwrap();
                let observed = C64::ONE
                    / (C64::new(0.0, -omega)
                        * response.mechanical.inlet_traction_impedance_pa_s_m3);
                near(
                    observed,
                    C64::from_re(if ends[0] == Adiabatic {
                        adiabatic
                    } else {
                        isothermal
                    }),
                    2e-6,
                );
                assert_eq!(response.mechanical.outlet_over_inlet_flow, C64::ZERO);
            }
            if alpha == -0.08 {
                assert!(
                    (adiabatic / reset - 2.0).abs() < 2e-15,
                    "separate adiabatic segment resets would lose half the compliance"
                );
            }
        }
    }

    // Independent dimensional mass, momentum, enthalpy and constitutive
    // equations, in [traction,theta,volume flow,heat power]. No scattering or
    // reference-junction algebra enters this through-region field oracle.
    fn primitive(
        m: &BulkThermoviscousEquilibrium,
        y: [C64; 4],
        area: f64,
        omega: f64,
    ) -> ([C64; 4], f64, C64) {
        let state = m.absorption().state();
        let rho = state.acoustic().density_kg_m3();
        let tr = state.transport();
        let eta = (4.0 / 3.0) * tr.dynamic_viscosity_pa_s() + m.absorption().bulk_viscosity_pa_s();
        let (alpha, t, beta, kappa, cp) = (
            m.isobaric_expansion_per_k(),
            m.temperature_k(),
            m.isothermal_compressibility_pa_inv(),
            tr.thermal_conductivity_w_m_k(),
            tr.specific_heat_j_kg_k(),
        );
        let pressure = (y[0] + C64::new(0.0, -omega * eta * alpha * t) * y[1])
            / C64::new(1.0, -omega * eta * beta);
        let ux = C64::new(0.0, omega) * (pressure.scale(beta) - y[1].scale(alpha * t));
        let tx = y[3].scale(-1.0 / (kappa * area * t));
        (
            [
                C64::new(0.0, omega * rho / area) * y[2],
                tx,
                ux.scale(area),
                C64::new(0.0, omega * area)
                    * (y[1].scale(rho * cp * t) - pressure.scale(alpha * t)),
            ],
            0.5 * area * (eta * ux.norm_sq() + kappa * t * tx.norm_sq()),
            pressure,
        )
    }

    #[test]
    fn g1_regional_fields_and_heat_close_original_nsf_across_material_contacts() {
        let regions = [region(0.08, 1.0, 0.03), region(-0.06, 1.7, 0.025)];
        let area = core::f64::consts::PI * RADIUS * RADIUS;
        let omega = 40.0;
        for ends in [[Adiabatic, Adiabatic], [Isothermal, Adiabatic]] {
            let line = prepare(&regions, omega, ends);
            let response = line
                .response(
                    Some(C64::new(1.3, 0.2).scale(line.references()[0])),
                    ComposeOptions::default(),
                )
                .unwrap();
            let m = response.mechanical;
            let inlet = [
                m.inlet_traction_impedance_pa_s_m3,
                response
                    .inlet_temperature_over_inlet_flow_k_s_m3
                    .scale(1.0 / 12.5),
                C64::ONE,
                response.inlet_heat_power_over_inlet_flow_pa,
            ];
            let outlet = [
                m.outlet_traction_over_inlet_flow_pa_s_m3,
                response
                    .outlet_temperature_over_inlet_flow_k_s_m3
                    .scale(1.0 / 12.5),
                m.outlet_over_inlet_flow,
                response.outlet_heat_power_over_inlet_flow_pa,
            ];
            near(
                m.inlet_pressure_over_flow_pa_s_m3,
                primitive(&regions[0].material, inlet, area, omega).2,
                2e-12,
            );
            near(
                m.outlet_pressure_over_inlet_flow_pa_s_m3,
                primitive(&regions[1].material, outlet, area, omega).2,
                2e-12,
            );
            let (mut y, mut heat, mut scale) = (inlet, 0.0, inlet.map(|z| z.abs()));
            for (index, region) in regions.iter().enumerate() {
                let steps = 1024;
                let dx = region.length_m / steps as f64;
                let mut local_heat = 0.0;
                for step in 0..=steps {
                    let (k1, source, _) = primitive(&region.material, y, area, omega);
                    local_heat += source
                        * if step == 0 || step == steps {
                            1.0
                        } else if step % 2 == 0 {
                            2.0
                        } else {
                            4.0
                        };
                    for i in 0..4 {
                        scale[i] = scale[i].max(y[i].abs());
                    }
                    if step == steps {
                        break;
                    }
                    let advance =
                        |k: [C64; 4], h: f64| core::array::from_fn(|i| y[i] + k[i].scale(h));
                    let k2 = primitive(&region.material, advance(k1, dx * 0.5), area, omega).0;
                    let k3 = primitive(&region.material, advance(k2, dx * 0.5), area, omega).0;
                    let k4 = primitive(&region.material, advance(k3, dx), area, omega).0;
                    y = core::array::from_fn(|i| {
                        y[i] + (k1[i] + k2[i].scale(2.0) + k3[i].scale(2.0) + k4[i]).scale(dx / 6.0)
                    });
                }
                heat += local_heat * dx / 3.0;
                if index == 0 {
                    // Ideal contact preserves THESE four fields, not pressure.
                    let left = primitive(&regions[0].material, y, area, omega).2;
                    let right = primitive(&regions[1].material, y, area, omega).2;
                    assert!((left - right).abs() > 1e-8 * left.abs());
                }
            }
            for i in 0..4 {
                assert!(
                    (y[i] - outlet[i]).abs() < 2e-9 * scale[i],
                    "field {i}: {:?} != {:?}",
                    y[i],
                    outlet[i]
                );
            }
            let flux = |z: [C64; 4]| 0.5 * (z[0] * z[2].conj() + z[1] * z[3].conj()).re;
            let work = flux(inlet) - flux(outlet);
            assert!(
                heat > 0.0 && (work / heat - 1.0).abs() < 2e-9,
                "work {work}, heat {heat}"
            );
        }
    }

    #[test]
    fn g3_regional_reference_changes_splitting_and_opaque_regrouping_are_physical() {
        let options = ComposeOptions::default();
        let ends = [Adiabatic, Adiabatic];
        let original = region(0.08, 1.0, 0.04);
        let mut rebased = ThermalNetwork::new(
            &original.material,
            RADIUS,
            original.length_m,
            100.0,
            options.maximum_condition,
        )
        .unwrap();
        let reference =
            BulkThermalCylinder::from_network(rebased.clone(), ends, options.maximum_condition)
                .unwrap();
        rebased
            .rebase(
                [rebased.references[0] * 9.0, rebased.references[1] * 0.04],
                options.maximum_condition,
                &mut || false,
            )
            .unwrap();
        let changed =
            BulkThermalCylinder::from_network(rebased, ends, options.maximum_condition).unwrap();
        let load = Some(C64::new(1.2, 0.3).scale(reference.references()[0]));
        let a = reference.response(load, options).unwrap();
        let b = changed.response(load, options).unwrap();
        near(
            a.mechanical.inlet_traction_impedance_pa_s_m3,
            b.mechanical.inlet_traction_impedance_pa_s_m3,
            2e-11,
        );
        near(
            a.mechanical.outlet_pressure_over_inlet_flow_pa_s_m3,
            b.mechanical.outlet_pressure_over_inlet_flow_pa_s_m3,
            2e-11,
        );
        near(
            a.outlet_temperature_over_inlet_flow_k_s_m3,
            b.outlet_temperature_over_inlet_flow_k_s_m3,
            2e-11,
        );
        for length in [0.04, 10000.0] {
            let single = BulkThermalCylinder::new(
                &original.material,
                RADIUS,
                length,
                1000.0,
                ends,
                options.maximum_condition,
            )
            .unwrap();
            let part = BulkThermalRegion {
                length_m: length / 4.0,
                ..original
            };
            let split = prepare(&[part; 4], 1000.0, ends);
            for (a, b) in single.scattering().iter().zip(split.scattering()) {
                assert!((*a - b).abs() < 2e-11);
            }
            let second = region(-0.06, 1.7, length / 4.0);
            let half = prepare(&[part, second], 1000.0, [Isothermal, Adiabatic]);
            let whole = prepare(&[part, second, part, second], 1000.0, ends);
            let (joined, _) = redheffer_star(
                &half.scattering(),
                &half.scattering(),
                2,
                options.maximum_condition,
                || false,
            )
            .unwrap();
            for (a, b) in whole.scattering().iter().zip(joined) {
                assert!((*a - b).abs() < 2e-11);
            }
        }
    }

    #[test]
    fn g0_regional_cylinder_admits_complete_same_temperature_chains_only() {
        let valid = region(0.08, 1.0, 0.04);
        let ends = [Adiabatic, Adiabatic];
        let options = ComposeOptions::default();
        assert!(BulkThermalCylinder::from_regions(&[], RADIUS, 100.0, ends, options).is_err());
        assert!(
            BulkThermalCylinder::from_regions(
                &[valid; 2],
                RADIUS,
                100.0,
                ends,
                ComposeOptions {
                    maximum_blocks: 1,
                    ..options
                }
            )
            .is_err()
        );
        let wrong = BulkThermalRegion {
            material: material(-0.06, 1.7, 13.0),
            ..valid
        };
        assert!(
            BulkThermalCylinder::from_regions(&[valid, wrong], RADIUS, 100.0, ends, options)
                .is_err()
        );
        let mut polls = 0;
        let result = BulkThermalCylinder::from_regions_controlled(
            &[valid; 2],
            RADIUS,
            100.0,
            ends,
            options,
            || {
                polls += 1;
                polls >= 3
            },
        );
        assert!(matches!(
            result,
            Err(BulkThermalError::Duct(DuctError::Cancelled { .. }))
        ));
        for length in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(
                BulkThermalCylinder::from_regions(
                    &[
                        valid,
                        BulkThermalRegion {
                            length_m: length,
                            ..valid
                        }
                    ],
                    RADIUS,
                    100.0,
                    ends,
                    options
                )
                .is_err()
            );
        }
    }
}
