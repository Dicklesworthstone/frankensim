//! Linearized heave plus two principal-axis rotations of one rigid backing.
//! Inertias are about the declared centre of mass, not contact-effective masses.
//! The centre's horizontal translations and yaw remain constrained. Springs and
//! resistances are diagonal in these physical coordinates; no full stand tensor
//! or finite rotation is inferred. Every patch uses its actual XY lever arm.
use super::{CarrierPort, Error, MountSite, MountedFelt};
use fs_couple::render::plate::impact::BodyPotential;
use fs_couple::render::plate::impact::felt::mount::modal::{ModalMountSite, ModalMountedFelt};
use fs_couple::modal_acoustic_time::ModalAcousticState;
use std::io::Write;

#[derive(Debug, Clone, Copy)]
pub struct Support {
    pub centre_xy_m: [f64; 2],
    pub inertia_kg_m2: [f64; 2],
    pub stiffness_nm_rad: [f64; 2],
    pub resistance_nm_s_rad: [f64; 2],
}
impl Support {
    pub fn validate(&self) -> Result<(), Error> {
        if self.centre_xy_m.iter().any(|x| !x.is_finite()) {
            return Err("rocking mount needs a finite centre of mass".into());
        }
        for i in 0..2 {
            let inertia = self.inertia_kg_m2[i];
            let stiffness = self.stiffness_nm_rad[i];
            let resistance = self.resistance_nm_s_rad[i];
            let inverse_root = 1.0 / inertia.sqrt();
            let frequency_squared = stiffness / inertia;
            let drag = resistance / inertia;
            if !inertia.is_finite() || inertia <= 0.0
                || !stiffness.is_finite() || stiffness < 0.0
                || !resistance.is_finite() || resistance < 0.0
                || !inverse_root.is_finite() || inverse_root <= 0.0
                || !frequency_squared.is_finite() || !drag.is_finite()
                || (stiffness > 0.0 && frequency_squared == 0.0)
                || (resistance > 0.0 && drag == 0.0) {
                return Err("rocking mount needs positive principal inertias and representable passive coefficients".into());
            }
        }
        Ok(())
    }
    /// d(vertical backing motion)/d(mass-normalized angular coordinate).
    /// Right-handed motion is u_z + (y-y_cm)*theta_x - (x-x_cm)*theta_y.
    fn vertical_row(&self, xy_m: [f64; 2]) -> Result<[f64; 2], Error> {
        let levers = [xy_m[1] - self.centre_xy_m[1], self.centre_xy_m[0] - xy_m[0]];
        let result = [levers[0] / self.inertia_kg_m2[0].sqrt(),
            levers[1] / self.inertia_kg_m2[1].sqrt()];
        if result.iter().any(|v| !v.is_finite())
            || (0..2).any(|i| levers[i] != 0.0 && result[i] == 0.0) {
            return Err("rocking mount point lever cannot be represented".into());
        }
        Ok(result)
    }
    pub fn assemble(&self, structural: usize, base: MountedFelt,
        sites: &[MountSite], positions: &[[f64; 2]]) -> Result<MountedFelt, Error>
    {
        self.validate()?;
        if positions.len() != sites.len() { return Err("missing rocking mount point positions".into()); }
        let MountedFelt { mut body, port, .. } = base;
        let BodyPotential::Linear(omega) = &mut body.potential else {
            return Err("rocking adapter requires the existing translating carrier".into());
        };
        if omega.len() != 1 { return Err("rocking adapter requires one original translation".into()); }
        for i in 0..2 {
            omega.push((self.stiffness_nm_rad[i] / self.inertia_kg_m2[i]).sqrt());
            // Reference-aligned, at rest in rotation. Translation retains the
            // original explicit initial condition. No hidden angular launch.
            body.initial.push(ModalAcousticState::default());
            body.damping_per_s.push(self.resistance_nm_s_rad[i] / self.inertia_kg_m2[i]);
        }
        let mut rows = Vec::with_capacity(sites.len());
        for (site, &xy) in sites.iter().zip(positions) {
            let angular = self.vertical_row(xy)?;
            let sign = site.compression_per_carrier_m;
            rows.push(ModalMountSite { pad: site.pad.clone(), support_weights: vec![
                sign * port.inverse_sqrt_mass, sign * angular[0], sign * angular[1]] });
        }
        let support = ModalMountedFelt::new(structural, body, &rows)?;
        Ok(MountedFelt { body: support.body, pads: support.pads, port })
    }
}

/// Observe accepted support motion without touching material or acoustic states.
#[derive(Debug, Clone, Copy)]
pub struct Port {
    pub translation: CarrierPort,
    pub inverse_sqrt_inertia: Option<[f64; 2]>,
}
impl Port {
    pub fn motion(&self, state: &[f64]) -> Result<[f64; 2], Error> {
        Ok(self.translation.motion(state)?)
    }
    pub fn angular_motion(&self, state: &[f64]) -> Result<Option<[[f64; 2]; 2]>, Error> {
        let Some(weights) = self.inverse_sqrt_inertia else { return Ok(None); };
        let first = self.translation.coordinate.checked_add(1).and_then(|i| i.checked_mul(2))
            .ok_or("rocking observation address overflow")?;
        let end = first.checked_add(4).ok_or("rocking observation extent overflow")?;
        let values = state.get(first..end).ok_or("short rocking mount state")?;
        let mut result = [[0.0; 2]; 2];
        for i in 0..2 {
            if !weights[i].is_finite() || weights[i] <= 0.0 { return Err("invalid rotational scale".into()); }
            result[i] = [values[2*i] * weights[i], values[2*i+1] * weights[i]];
        }
        if result.iter().flatten().any(|v| !v.is_finite()) { return Err("nonfinite rocking mount motion".into()); }
        Ok(Some(result))
    }
    pub fn header(&self, out: &mut impl Write) -> Result<(), Error> {
        super::header(out)?;
        if self.inverse_sqrt_inertia.is_some() {
            write!(out, ",mount_rotation_x_rad,mount_angular_velocity_x_rad_s,mount_rotation_y_rad,mount_angular_velocity_y_rad_s")?;
        }
        Ok(())
    }
}

#[cfg(test)]
#[path = "rocking_tests.rs"]
mod tests;
