//! Certified compact radiators from a thin plate.
//!
//! Geometry + orthotropic section go in; `fs-plate` assembles the DKT
//! pencil and `fs-modal` returns inertia-certified eigenpairs. A
//! guitar top, a bulkhead, and a panel are the same object.

use crate::acoustic_realize::AcousticRealizeError;
use fs_material::elastic::OrthotropicElastic;
use fs_plate::{
    AssemblyOptions, EdgeSupport, PlateError, PlateMesh, PlateSection, assemble, modes,
};
use fs_scenario::{RadiatingPlate, ThinPlate};

/// One driven compact radiator harvested from a certified plate mode.
#[derive(Debug, Clone)]
pub struct CompactBody {
    /// Radiating monopole area [m²] (`∫ φ dA`).
    pub area_m2: f64,
    /// Modal mass [kg] (1 for a mass-normalized mode).
    pub mass_kg: f64,
    /// Angular frequency [rad/s].
    pub omega: f64,
    /// Viscous damping ratio.
    pub zeta: f64,
    y: f64,
    v: f64,
}

impl CompactBody {
    /// From a caller-supplied compact monopole.
    ///
    /// # Errors
    /// Non-physical parameters.
    pub fn from_radiator(spec: RadiatingPlate) -> Result<Self, AcousticRealizeError> {
        if !(spec.area_m2 > 0.0
            && spec.mass_kg > 0.0
            && spec.frequency_hz > 0.0
            && spec.damping_ratio >= 0.0)
        {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "compact radiator parameters must be physical",
            });
        }
        Ok(Self {
            area_m2: spec.area_m2,
            mass_kg: spec.mass_kg,
            omega: core::f64::consts::TAU * spec.frequency_hz,
            zeta: spec.damping_ratio,
            y: 0.0,
            v: 0.0,
        })
    }

    /// Advance under a generalized force and return acceleration.
    pub fn drive(&mut self, force_n: f64, dt: f64) -> f64 {
        let acc = force_n / self.mass_kg
            - 2.0 * self.zeta * self.omega * self.v
            - self.omega * self.omega * self.y;
        self.v += dt * acc;
        self.y += dt * self.v;
        acc
    }

    /// Compact monopole pressure at distance `listener_m`.
    #[must_use]
    pub fn radiate(&self, acc: f64, rho: f64, listener_m: f64) -> f64 {
        rho * self.area_m2 * acc / (4.0 * core::f64::consts::PI * listener_m)
    }

    /// Volume velocity of the monopole [m³/s].
    #[must_use]
    pub fn volume_velocity(&self) -> f64 {
        self.area_m2 * self.v
    }

    /// Drive and radiate in one step.
    pub fn drive_and_radiate(&mut self, force_n: f64, dt: f64, rho: f64, listener_m: f64) -> f64 {
        let acc = self.drive(force_n, dt);
        self.radiate(acc, rho, listener_m)
    }
}

/// Certified compact radiators of a thin plate.
///
/// # Errors
/// Section, mesh, or modal-window refusals.
pub fn certified_radiators(plate: ThinPlate) -> Result<Vec<CompactBody>, AcousticRealizeError> {
    if plate.n_modes == 0 {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "thin plate needs at least one mode",
        });
    }
    let law = OrthotropicElastic::new(
        [plate.e1_pa, plate.e2_pa, plate.e2_pa],
        [plate.nu12, plate.nu12, plate.nu12],
        [plate.g12_pa, plate.g12_pa, plate.g12_pa],
        1.0,
    )
    .map_err(|_| AcousticRealizeError::InvalidDescription {
        what: "plate elastic constants refused",
    })?;
    let section = PlateSection::orthotropic(&law, plate.thickness_m, plate.density_kg_m3)
        .map_err(map_plate)?;
    let (nx, ny) = (5, 4);
    let mesh = PlateMesh::rectangle(plate.length_m, plate.width_m, nx, ny);
    let boundary = PlateMesh::rectangle_boundary(nx, ny);
    let model = assemble(
        &mesh,
        &section,
        &boundary,
        &[],
        &AssemblyOptions {
            pretension: 0.0,
            support: EdgeSupport::SimplySupported,
        },
    )
    .map_err(map_plate)?;
    let omega11 = ss_omega11(&section, plate.length_m, plate.width_m);
    let lo = (0.25 * omega11).powi(2);
    let hi = (omega11 * (plate.n_modes as f64 + 2.0) * 4.0).powi(2);
    let report = modes(&model, (lo, hi), &fs_modal::SliceOptions::default()).map_err(map_plate)?;
    if report.modes.is_empty() {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "plate modal window returned no certified modes",
        });
    }
    let n_keep = plate.n_modes.min(report.modes.len());
    let nn = mesh.node_count();
    let nodal_area = plate.length_m * plate.width_m / nn as f64;
    let mut out = Vec::with_capacity(n_keep);
    for pair in report.modes.iter().take(n_keep) {
        let omega = pair.lambda.max(0.0).sqrt();
        if !(omega > 0.0) {
            continue;
        }
        let mut monopole = 0.0;
        let mut drive = 0.0;
        let mut drive_w = 0.0;
        let mid_y = 0.5 * plate.width_m;
        for (inode, &(x, y)) in mesh.nodes.iter().enumerate() {
            let Some(reduced) = model.dof_map[3 * inode] else {
                continue;
            };
            let phi = pair.phi.get(reduced).copied().unwrap_or(0.0);
            monopole += phi * nodal_area;
            if x <= plate.length_m / (nx as f64).max(1.0) {
                let w = 1.0 / (1.0 + (y - mid_y).abs());
                drive += phi * w;
                drive_w += w;
            }
        }
        let phi_drive = if drive_w > 0.0 { drive / drive_w } else { 1.0 };
        let mass = if phi_drive.abs() > 1.0e-12 {
            1.0 / phi_drive.abs()
        } else {
            1.0
        };
        out.push(CompactBody {
            area_m2: monopole.abs().max(1.0e-8),
            mass_kg: mass,
            omega,
            zeta: plate.damping_ratio,
            y: 0.0,
            v: 0.0,
        });
    }
    if out.is_empty() {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "no usable plate radiators after harvesting",
        });
    }
    Ok(out)
}

fn ss_omega11(section: &PlateSection, a: f64, b: f64) -> f64 {
    let pi = core::f64::consts::PI;
    let d = section.d[0].abs().max(1.0e-18);
    let rho_h = (section.density * section.thickness).max(1.0e-18);
    pi * pi * (1.0 / (a * a) + 1.0 / (b * b)) * (d / rho_h).sqrt()
}

fn map_plate(_err: PlateError) -> AcousticRealizeError {
    AcousticRealizeError::InvalidDescription {
        what: "thin-plate modal solve refused",
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn simply_supported_first_mode_is_near_the_analytic_value() {
        let plate = ThinPlate {
            length_m: 0.20,
            width_m: 0.15,
            thickness_m: 0.002,
            density_kg_m3: 7800.0,
            e1_pa: 200e9,
            e2_pa: 200e9,
            nu12: 0.3,
            g12_pa: 200e9 / (2.0 * 1.3),
            damping_ratio: 0.01,
            n_modes: 2,
        };
        let bodies = certified_radiators(plate).expect("modes");
        let section = PlateSection::isotropic(200e9, 0.3, 0.002, 7800.0).expect("section");
        let want = ss_omega11(&section, 0.20, 0.15);
        let got = bodies[0].omega;
        assert!(
            (got - want).abs() / want < 0.20,
            "certified ω={got:.1} vs SS analytic {want:.1}"
        );
    }
}
