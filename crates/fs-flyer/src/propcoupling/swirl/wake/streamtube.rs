//! Material annular streamtubes from the original native blade solution.
//!
//! For annulus i, m_i=rho*(V+a_i)*A_i and ell_i=Q_i/m_i stay fixed.
//! U_i(s)=V+[1+s/sqrt(s²+R²)]*a_i, A_i(s)=m_i/(rho*U_i).
//! Cumulative areas give contiguous annuli around an unloaded cylindrical hub.
//! With b_i=r_inner² and U_i'=dU_i/ds, continuity gives
//! v_r=[U_i*b_i'-U_i'*(r²-b_i)]/(2r). At either interface v_r/U_i=dr/ds;
//! therefore no mass crosses a material interface even when axial velocity
//! jumps. Specific angular momentum gives v_theta=ell_i/r.
//!
//! These are steady axisymmetric SOURCE streamtubes, not an Euler pressure
//! solution, exterior-flow solution, turbulent mixing law or post-wing wake.
//! The prescribed axial development retains the parent's finite 2R window.
//! Meridional kinetic power is visible, not debited from swirl or called heat.
//! See NASA CR-3887, sections 7.1 (axial induced velocity) and 7.2 (swirl),
//! for the distinction between streamtube contraction and circulation transport:
//! <https://ntrs.nasa.gov/citations/19850014058>.

use super::{AnnularWake, CrossSection, Refusal, Solution, WakeControls, array, det, reject, vec};

/// One physical section of an original native annular streamtube.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialAnnulus {
    pub inner_radius_m: f64,
    pub outer_radius_m: f64,
    pub axial_speed_mps: f64,
    /// Derivatives dr/ds of the two material surfaces.
    pub radius_derivatives: [f64; 2],
    pub mass_flux_kg_s: f64,
    pub angular_momentum_flux_nm: f64,
    pub resolved_swirl_power_w: f64,
    /// Resolved outward transport of radial kinetic energy [W].
    /// No radial-pressure or total-energy closure is inferred from this number.
    pub radial_kinetic_power_w: f64,
    /// Axial kinetic transport above the same native incident V [W].
    pub axial_kinetic_power_increment_w: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct Transport {
    induction_mps: f64,
    mass_flux_kg_s: f64,
    specific_angular_momentum_m2_s: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub(super) struct Streamtubes {
    cells: Vec<Transport>,
    maximum_radius_m: f64,
    maximum_power_defect_w: f64,
}

#[derive(Clone, Copy)]
struct Cell {
    inner_squared: f64,
    width_squared: f64,
    inner_squared_derivative: f64,
    outer_squared_derivative: f64,
    axial: f64,
    axial_derivative: f64,
    transport: Transport,
}

fn finite(values: &[f64]) -> Result<(), Refusal> {
    if values.iter().all(|x| x.is_finite()) {
        Ok(())
    } else {
        Err(reject(
            "swirl-wake-invalid",
            "unrepresentable material streamtube",
        ))
    }
}

fn axial_factor(radius: f64, distance: f64) -> (f64, f64) {
    let x = distance / radius;
    let root = det::sqrt(1.0 + x * x);
    (1.0 + x / root, 1.0 / (radius * root * root * root))
}

// Stable positive moments of 1/(1+e*t), t in [-1,1]. For narrow annuli the
// convergent series avoids subtracting nearly equal logarithms. M1=-e*M2 and
// M2-M1²/M0=2*M2/M0 give a positive completed-square radial-power integral.
fn radial_moments(inner_squared: f64, width_squared: f64) -> (f64, f64, f64) {
    let midpoint = inner_squared + 0.5 * width_squared;
    let e = 0.5 * width_squared / midpoint;
    let e2 = e * e;
    if e <= 0.125 {
        let mut term = 1.0;
        let mut m2 = 0.0;
        for n in 0..16 {
            m2 += 2.0 * term / (2 * n + 3) as f64;
            term *= e2;
        }
        (2.0 + e2 * m2, m2, e)
    } else {
        let outer = det::sqrt(inner_squared + width_squared);
        let inner = det::sqrt(inner_squared);
        let m0 = 2.0 * (det::ln(outer) - det::ln(inner)) / e;
        (m0, (m0 - 2.0) / e2, e)
    }
}

impl Cell {
    fn powers(self, incident: f64) -> Result<[f64; 3], Refusal> {
        let t = self.transport;
        let midpoint = self.inner_squared + 0.5 * self.width_squared;
        let (m0, m2, e) = radial_moments(self.inner_squared, self.width_squared);
        let ell = t.specific_angular_momentum_m2_s;
        let swirl = 0.25 * t.mass_flux_kg_s * (ell * ell / midpoint) * m0;
        let k = 0.25 * self.axial_derivative * self.width_squared;
        let h = 0.5 * self.axial * self.inner_squared_derivative - k;
        let center = h + k * (e * m2 / m0);
        let radial =
            0.25 * t.mass_flux_kg_s / midpoint * (m0 * center * center + k * k * (2.0 * m2 / m0));
        let axial = 0.5 * t.mass_flux_kg_s * (self.axial - incident) * (self.axial + incident);
        finite(&[swirl, radial, axial])?;
        Ok([swirl, radial, axial])
    }

    fn report(self, wake: &AnnularWake) -> Result<MaterialAnnulus, Refusal> {
        let inner = det::sqrt(self.inner_squared);
        let outer = det::sqrt(self.inner_squared + self.width_squared);
        let [swirl, radial, axial] = self.powers(wake.incident_axial_mps)?;
        let result = MaterialAnnulus {
            inner_radius_m: inner,
            outer_radius_m: outer,
            axial_speed_mps: self.axial,
            radius_derivatives: [
                0.5 * self.inner_squared_derivative / inner,
                0.5 * self.outer_squared_derivative / outer,
            ],
            mass_flux_kg_s: self.transport.mass_flux_kg_s,
            angular_momentum_flux_nm: wake.spin_sign
                * self.transport.mass_flux_kg_s
                * self.transport.specific_angular_momentum_m2_s,
            resolved_swirl_power_w: swirl,
            radial_kinetic_power_w: radial,
            axial_kinetic_power_increment_w: axial,
        };
        finite(&result.radius_derivatives)?;
        Ok(result)
    }
}

impl Streamtubes {
    pub(super) fn new(
        wake: &AnnularWake,
        solution: &Solution,
        controls: WakeControls,
    ) -> Result<Self, Refusal> {
        let mut cells = Vec::with_capacity(wake.cells.len());
        for (source, station) in wake.cells.iter().zip(&solution.stations) {
            let induction = station.induction.axial_mps;
            let speed = wake.incident_axial_mps + induction;
            let mass_flux = wake.density * speed * source.area();
            let specific = source.torque_nm / mass_flux;
            finite(&[induction, speed, mass_flux, specific])?;
            if speed <= 0.0 || mass_flux <= 0.0 {
                return Err(reject(
                    "swirl-wake-outside-domain",
                    "material annulus requires positive native disk throughflow",
                ));
            }
            cells.push(Transport {
                induction_mps: induction,
                mass_flux_kg_s: mass_flux,
                specific_angular_momentum_m2_s: specific,
            });
        }
        let mut result = Self {
            cells,
            maximum_radius_m: 0.0,
            maximum_power_defect_w: 0.0,
        };
        // All U_i are affine in the SAME monotone axial factor f(s).
        // Bound each reciprocal area independently over f in [1,1+2/sqrt(5)].
        // Summing these positive bounds is conservative even for mixed-sign
        // induction; assuming the cumulative radius is monotone would be wrong.
        let far_factor = axial_factor(wake.radius_m, 2.0 * wake.radius_m).0;
        let hub_squared = wake.cells[0].inner_radius_m * wake.cells[0].inner_radius_m;
        let mut minimum_inner_squared = hub_squared;
        let mut maximum_inner_squared = hub_squared;
        let mut native_power = 0.0;
        for (t, source) in result.cells.iter().zip(&wake.cells) {
            let near = wake.incident_axial_mps + t.induction_mps;
            let far = wake.incident_axial_mps + far_factor * t.induction_mps;
            let (low, high) = (near.min(far), near.max(far));
            if low <= 0.0 || !high.is_finite() {
                return Err(reject(
                    "swirl-wake-outside-domain",
                    "material annulus reverses throughflow inside the observation window",
                ));
            }
            let volume_flux_over_pi = t.mass_flux_kg_s / (core::f64::consts::PI * wake.density);
            let minimum_width = volume_flux_over_pi / high;
            let maximum_width = volume_flux_over_pi / low;
            // P(U,b)=pi*rho*U*ell²/2*ln(1+C/(U*b)) increases with U and
            // decreases with b. Independent U/b interval endpoints bound every
            // interior cell power; opposite cell defects never cancel.
            let power_at = |speed: f64, inner: f64| -> Result<f64, Refusal> {
                let width = volume_flux_over_pi / speed;
                let cell = Cell {
                    inner_squared: inner,
                    width_squared: width,
                    inner_squared_derivative: 0.0,
                    outer_squared_derivative: 0.0,
                    axial: speed,
                    axial_derivative: 0.0,
                    transport: *t,
                };
                Ok(cell.powers(wake.incident_axial_mps)?[0])
            };
            let lower_power = power_at(low, maximum_inner_squared)?;
            let upper_power = power_at(high, minimum_inner_squared)?;
            result.maximum_power_defect_w += (lower_power - source.native_swirl_power_w)
                .abs()
                .max((upper_power - source.native_swirl_power_w).abs());
            native_power += source.native_swirl_power_w;
            minimum_inner_squared += minimum_width;
            maximum_inner_squared += maximum_width;
            finite(&[
                minimum_inner_squared,
                maximum_inner_squared,
                result.maximum_power_defect_w,
                native_power,
            ])?;
        }
        result.maximum_radius_m = det::sqrt(maximum_inner_squared);
        let budget =
            controls.absolute_power_tolerance_w + controls.relative_power_tolerance * native_power;
        if !budget.is_finite() || result.maximum_power_defect_w > budget {
            return Err(reject(
                "swirl-wake-power-budget",
                &format!(
                    "material-annulus projection bound {} W exceeds admitted {} W",
                    result.maximum_power_defect_w, budget
                ),
            ));
        }
        // Admit physical endpoint geometry, including representable annular width.
        for distance in [0.0, 2.0 * wake.radius_m] {
            for cell in result.sections(wake, distance)? {
                cell?.report(wake)?;
            }
        }
        Ok(result)
    }

    pub(super) fn maximum_radius_m(&self) -> f64 {
        self.maximum_radius_m
    }

    pub(super) fn maximum_power_defect_w(&self) -> f64 {
        self.maximum_power_defect_w
    }

    fn sections<'a>(
        &'a self,
        wake: &'a AnnularWake,
        distance: f64,
    ) -> Result<impl Iterator<Item = Result<Cell, Refusal>> + 'a, Refusal> {
        if !distance.is_finite() || distance < 0.0 || distance > 2.0 * wake.radius_m {
            return Err(reject(
                "swirl-wake-outside-domain",
                "material cross section outside 0..=2R",
            ));
        }
        let (factor, derivative) = axial_factor(wake.radius_m, distance);
        let hub = wake.cells[0].inner_radius_m;
        Ok(self
            .cells
            .iter()
            .scan((hub * hub, 0.0), move |state, &transport| {
                let cell = (|| {
                    let axial = wake.incident_axial_mps + factor * transport.induction_mps;
                    let axial_derivative = derivative * transport.induction_mps;
                    let width =
                        transport.mass_flux_kg_s / (core::f64::consts::PI * wake.density * axial);
                    let width_derivative = -width * (axial_derivative / axial);
                    let outer_squared = state.0 + width;
                    let outer_derivative = state.1 + width_derivative;
                    finite(&[
                        axial,
                        axial_derivative,
                        width,
                        width_derivative,
                        outer_squared,
                        outer_derivative,
                    ])?;
                    if axial <= 0.0 || width <= 0.0 || outer_squared <= state.0 {
                        return Err(reject(
                            "swirl-wake-outside-domain",
                            "unresolved material annular width or throughflow",
                        ));
                    }
                    let cell = Cell {
                        inner_squared: state.0,
                        width_squared: width,
                        inner_squared_derivative: state.1,
                        outer_squared_derivative: outer_derivative,
                        axial,
                        axial_derivative,
                        transport,
                    };
                    *state = (outer_squared, outer_derivative);
                    Ok(cell)
                })();
                Some(cell)
            }))
    }

    pub(super) fn material_annuli(
        &self,
        wake: &AnnularWake,
        distance: f64,
    ) -> Result<Vec<MaterialAnnulus>, Refusal> {
        self.sections(wake, distance)?
            .map(|cell| cell?.report(wake))
            .collect()
    }

    pub(super) fn cross_section(
        &self,
        wake: &AnnularWake,
        distance: f64,
    ) -> Result<CrossSection, Refusal> {
        let mut flux = CrossSection {
            axial_speed_mps: 0.0,
            mass_flux_kg_s: 0.0,
            angular_momentum_flux_nm: 0.0,
            resolved_swirl_power_w: 0.0,
            native_swirl_power_w: 0.0,
            power_defect_w: 0.0,
            absolute_cell_power_defect_w: 0.0,
        };
        let mut area = 0.0;
        for (cell, source) in self.sections(wake, distance)?.zip(&wake.cells) {
            let cell = cell?;
            let power = cell.powers(wake.incident_axial_mps)?[0];
            area += core::f64::consts::PI * cell.width_squared;
            flux.mass_flux_kg_s += cell.transport.mass_flux_kg_s;
            flux.angular_momentum_flux_nm += wake.spin_sign * source.torque_nm;
            flux.resolved_swirl_power_w += power;
            flux.native_swirl_power_w += source.native_swirl_power_w;
            flux.absolute_cell_power_defect_w += (power - source.native_swirl_power_w).abs();
        }
        flux.axial_speed_mps = flux.mass_flux_kg_s / (wake.density * area);
        flux.power_defect_w = flux.resolved_swirl_power_w - flux.native_swirl_power_w;
        finite(&[
            flux.axial_speed_mps,
            flux.mass_flux_kg_s,
            flux.angular_momentum_flux_nm,
            flux.resolved_swirl_power_w,
            flux.native_swirl_power_w,
            flux.power_defect_w,
            flux.absolute_cell_power_defect_w,
        ])?;
        Ok(flux)
    }

    pub(super) fn relative_velocity_at(
        &self,
        wake: &AnnularWake,
        point: [f64; 3],
    ) -> Result<[f64; 3], Refusal> {
        finite(&point)?;
        let offset = vec(point) - vec(wake.center_m);
        let axial_position = offset.dot(vec(wake.axis));
        let radial = offset - vec(wake.axis).scale(axial_position);
        let radius_squared = radial.dot(radial);
        finite(&[axial_position, radius_squared])?;
        if axial_position > 0.0 {
            return Ok([0.0; 3]);
        }
        let distance = -axial_position;
        let hub = wake.cells[0].inner_radius_m;
        if radius_squared < hub * hub
            || radius_squared >= self.maximum_radius_m * self.maximum_radius_m
        {
            return Ok([0.0; 3]);
        }
        if distance <= 0.0 || distance >= 2.0 * wake.radius_m {
            return Err(reject(
                "swirl-wake-outside-domain",
                "material-wake surface outside 0<s<2R",
            ));
        }
        for cell in self.sections(wake, distance)? {
            let cell = cell?;
            if radius_squared >= cell.inner_squared
                && radius_squared < cell.inner_squared + cell.width_squared
            {
                // Physical fluid v_r points outward; material-relative velocity
                // gets its NEGATIVE. Downstream fluid is along -thrust axis,
                // hence its induced axial material-relative component is +axis.
                let radial_coefficient = 0.5
                    * (cell.axial * cell.inner_squared_derivative
                        - cell.axial_derivative * (radius_squared - cell.inner_squared))
                    / radius_squared;
                let swirl_coefficient =
                    wake.spin_sign * cell.transport.specific_angular_momentum_m2_s / radius_squared;
                let result = array(
                    vec(wake.axis).scale(cell.axial - wake.incident_axial_mps)
                        - radial.scale(radial_coefficient)
                        - vec(wake.axis).cross(radial).scale(swirl_coefficient),
                );
                finite(&result)?;
                return Ok(result);
            }
        }
        Ok([0.0; 3])
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::propcoupling::swirl::wake::{Annulus, WakeTransport};
    use fs_airscrew::swirl::{Induction, Power, Station};

    #[test]
    fn mixed_native_induction_bounds_nonmonotone_radius_and_swirl_power() {
        // Prescribed analytic annuli isolate the transport law. Full native
        // blade solves and actual mechanics are exercised by swirl_wake tests.
        let radii = [0.25, det::sqrt(0.8125), 1.0];
        let inductions = [5.0, -4.0];
        let torques = [0.8, -0.3];
        let mut source = Vec::new();
        let mut stations = Vec::new();
        for i in 0..2 {
            let radius = 0.5 * (radii[i] + radii[i + 1]);
            let mass = 1.225
                * (10.0 + inductions[i])
                * core::f64::consts::PI
                * (radii[i + 1] * radii[i + 1] - radii[i] * radii[i]);
            let tangential = torques[i] / (2.0 * mass * radius);
            source.push(Annulus {
                inner_radius_m: radii[i],
                outer_radius_m: radii[i + 1],
                torque_nm: torques[i],
                native_swirl_power_w: tangential * torques[i] / radius,
            });
            stations.push(Station {
                r_over_r: radius,
                induction: Induction {
                    axial_mps: inductions[i],
                    tangential_mps: tangential,
                },
                evaluations: 0,
                alpha_rad: 0.0,
                prandtl_f: 1.0,
                thrust_n: 0.0,
                torque_nm: torques[i],
                axial_momentum_flux_n: 0.0,
                angular_momentum_flux_nm: torques[i],
                axial_residual_n: 0.0,
                angular_residual_nm: 0.0,
            });
        }
        let solution = Solution {
            thrust_n: 0.0,
            torque_nm: 0.5,
            ct: 0.0,
            cq: 0.0,
            j: 0.0,
            power: Power::default(),
            stations,
        };
        let mut wake = AnnularWake {
            center_m: [0.0; 3],
            axis: [1.0, 0.0, 0.0],
            spin_sign: 1.0,
            radius_m: 1.0,
            density: 1.225,
            incident_axial_mps: 10.0,
            disk_induction_mps: 0.0,
            cells: source,
            maximum_power_defect_w: 0.0,
            streamtubes: None,
        };
        let tubes = Streamtubes::new(
            &wake,
            &solution,
            WakeControls {
                relative_power_tolerance: 1.0,
                ..WakeControls::default()
            },
        )
        .unwrap();
        let bound = tubes.maximum_radius_m();
        wake.maximum_power_defect_w = tubes.maximum_power_defect_w();
        wake.streamtubes = Some(tubes);
        assert_eq!(wake.transport(), WakeTransport::MaterialStreamtubes);
        let initial = wake.material_annuli(0.0).unwrap()[1].outer_radius_m;
        assert!(wake.material_annuli(0.1).unwrap()[1].outer_radius_m < initial);
        assert!(wake.material_annuli(2.0).unwrap()[1].outer_radius_m > initial);
        let mass = wake.cross_section(0.0).unwrap().mass_flux_kg_s;
        for step in 0..=100 {
            let s = 2.0 * step as f64 / 100.0;
            let actual = wake.cross_section(s).unwrap();
            assert_eq!(actual.mass_flux_kg_s, mass);
            assert!((actual.angular_momentum_flux_nm - 0.5).abs() < 1e-15);
            assert!(actual.absolute_cell_power_defect_w <= wake.maximum_power_defect_w + 1e-14);
            assert!(wake.material_annuli(s).unwrap()[1].outer_radius_m <= bound);
        }
    }
}
