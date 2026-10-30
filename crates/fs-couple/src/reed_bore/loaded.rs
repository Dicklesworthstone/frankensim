//! A moving aperture closed against the complete body's same-window pressure.
//!
//! The body supplies its exact held-volume-flow relation `p_bar = free + B q`.
//! For `B > 0` this is the existing characteristic junction with impedance B
//! and incoming pressure free/2. The original scalar contact solve therefore
//! closes pressure, changing aperture and swept-face flow together. No body
//! state is probed, reduced, advanced or overwritten by a reed trial.
//!
//! The original mass/spring/damper block evolves exactly under the interval's
//! held pressure/contact force; contact uses fs-dcontact's potential secant.
//! The aperture is evaluated at the midpoint gap, so nonlinear time accuracy
//! still requires refinement. This is not an exact nonlinear continuous flow.

use super::{
    AcousticRealizeError, BeatingReed, OpeningPropagator, reed_parameters_valid,
    reed_pressure_face, reed_structural, reed_swept_flow,
};
use crate::bernoulli_aperture::{
    moving,
    paired::{BladeState, ExactBlade},
};
use fs_dcontact::{AffineContactStep, ContactStorage, Obstacle};

/// Original reed/contact storage and independently evaluated work and heat [J].
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct LoadedReedAccount {
    pub(crate) energy_before_j: f64,
    pub(crate) energy_after_j: f64,
    pub(crate) delta_h_j: f64,
    pub(crate) viscous_heat_j: f64,
    pub(crate) contact_heat_j: f64,
    pub(crate) jet_heat_j: f64,
    pub(crate) mouth_work_j: f64,
    pub(crate) bore_work_j: f64,
    pub(crate) residual_j: f64,
}

/// One accepted candidate; the caller publishes it with the complete body.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct LoadedReedStep {
    pub(crate) opening_m: f64,
    pub(crate) velocity_m_s: f64,
    pub(crate) mean_velocity_m_s: f64,
    pub(crate) pressure_pa: f64,
    pub(crate) volume_flow_m3_s: f64,
    pub(crate) jet_flow_m3_s: f64,
    pub(crate) account: LoadedReedAccount,
}

/// Quasistatic valve on the same held-flow boundary, with no mechanical store.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct QuasistaticLoadStep {
    pub(crate) pressure_pa: f64,
    pub(crate) volume_flow_m3_s: f64,
    pub(crate) jet_heat_j: f64,
    pub(crate) mouth_work_j: f64,
    pub(crate) bore_work_j: f64,
    pub(crate) residual_j: f64,
}

/// Close the original quasistatic pressure-dependent aperture against the
/// complete body's mean held-flow response. The existing strict characteristic
/// root retains its closure, branch selection, clamp and flow dead zone.
#[allow(clippy::too_many_arguments)]
pub(crate) fn solve_quasistatic_load(
    reed: BeatingReed,
    rho: f64,
    dt: f64,
    mouth_pressure_pa: f64,
    free_pressure_pa: f64,
    mean_impedance: f64,
    previous_flow_m3_s: f64,
) -> Result<QuasistaticLoadStep, AcousticRealizeError> {
    if !reed_parameters_valid(reed)
        || reed.mass_kg != 0.0
        || ![
            rho,
            dt,
            mouth_pressure_pa,
            free_pressure_pa,
            mean_impedance,
            previous_flow_m3_s,
        ]
        .iter()
        .all(|x| x.is_finite())
        || rho <= 0.0
        || dt <= 0.0
        || mean_impedance <= 0.0
    {
        return Err(refuse(
            "quasistatic loaded valve requires finite zero-mass parameters and positive gas, clock and response",
        ));
    }
    let incoming = 0.5 * free_pressure_pa;
    let guess = incoming + mean_impedance * previous_flow_m3_s;
    let outgoing = super::solve_reed_wave(
        reed,
        rho,
        mean_impedance,
        0.0,
        incoming,
        mouth_pressure_pa,
        guess,
        0.0,
    )?;
    let pressure_pa = outgoing + incoming;
    let drop = mouth_pressure_pa - pressure_pa;
    // Read the original constitutive law rather than dividing a rounded
    // pressure difference by a possibly tiny body response.
    let volume_flow_m3_s = super::aperture_of(reed).volume_flow(drop, rho);
    let jet_heat_j = dt * drop * volume_flow_m3_s;
    let mouth_work_j = dt * mouth_pressure_pa * volume_flow_m3_s;
    let bore_work_j = dt * pressure_pa * volume_flow_m3_s;
    let residual_j = jet_heat_j + bore_work_j - mouth_work_j;
    let boundary_residual = pressure_pa - (free_pressure_pa + mean_impedance * volume_flow_m3_s);
    let pressure_scale = pressure_pa.abs()
        + free_pressure_pa.abs()
        + mouth_pressure_pa.abs()
        + (mean_impedance * volume_flow_m3_s).abs();
    let work_scale = jet_heat_j.abs() + mouth_work_j.abs() + bore_work_j.abs();
    if ![
        pressure_pa,
        volume_flow_m3_s,
        jet_heat_j,
        mouth_work_j,
        bore_work_j,
        residual_j,
        boundary_residual,
        pressure_scale,
        work_scale,
    ]
    .iter()
    .all(|x| x.is_finite())
        || jet_heat_j < 0.0
        || residual_j.abs() > 64.0 * f64::EPSILON * work_scale.max(f64::MIN_POSITIVE)
        || boundary_residual.abs() > 64.0 * f64::EPSILON * pressure_scale.max(f64::MIN_POSITIVE)
    {
        return Err(refuse(
            "quasistatic loaded valve physical work or current-window flow balance refused",
        ));
    }
    Ok(QuasistaticLoadStep {
        pressure_pa,
        volume_flow_m3_s,
        jet_heat_j,
        mouth_work_j,
        bore_work_j,
        residual_j,
    })
}

/// Fixed-clock mechanics and original contact law. A successful warm solve
/// allocates no memory; exponentials and owned storage are prepared once.
pub(crate) struct PreparedReedLoad {
    reed: BeatingReed,
    rho: f64,
    dt: f64,
    stiffness: f64,
    face: f64,
    blade: ExactBlade,
    contact: ContactStorage,
}

fn refuse(what: &'static str) -> AcousticRealizeError {
    AcousticRealizeError::Reed { what }
}

impl PreparedReedLoad {
    pub(crate) fn new(
        reed: BeatingReed,
        rho: f64,
        dt: f64,
        lay: Obstacle,
    ) -> Result<Self, AcousticRealizeError> {
        if !reed_parameters_valid(reed)
            || reed.mass_kg <= 0.0
            || !rho.is_finite()
            || rho <= 0.0
            || !dt.is_finite()
            || dt <= 0.0
        {
            return Err(refuse(
                "loaded reed requires finite massive mechanics, gas and positive duration",
            ));
        }
        AffineContactStep::new(
            &lay,
            reed.rest_opening_m,
            crate::unilateral_contact::distributed::MAX_APERTURE_CONTACT_POINTS,
        )
        .map_err(|error| AcousticRealizeError::Nonlinear(error.to_string()))?;
        let (stiffness, damping) = reed_structural(reed);
        let face = reed_pressure_face(reed);
        if !face.is_finite() || face <= 0.0 {
            return Err(refuse(
                "loaded reed pressure area must be finite and positive",
            ));
        }
        let blade = ExactBlade::from_mechanics(reed.mass_kg, damping, stiffness, dt)?;
        // ContactStorage remains the owner of the physical potential. Its
        // probe reads the original admitted law without allocating a trial.
        let zero = fs_phs::QuadraticStorage::new(vec![0.0; 4], 2)
            .map_err(|error| AcousticRealizeError::Nonlinear(error.to_string()))?;
        let contact = ContactStorage::new(Box::new(zero), 1, vec![lay])
            .map_err(|error| AcousticRealizeError::Nonlinear(error.to_string()))?;
        Ok(Self {
            reed,
            rho,
            dt,
            stiffness,
            face,
            blade,
            contact,
        })
    }

    fn energy(&self, opening: f64, velocity: f64) -> f64 {
        let displacement = opening - self.reed.rest_opening_m;
        0.5 * self.reed.mass_kg * velocity * velocity
            + 0.5 * self.stiffness * displacement * displacement
            + self.contact.probe(&[opening, 0.0]).contact_energy
    }

    /// Solve against `p_bar = free_pressure_pa + mean_impedance * q`.
    /// q is the held body volume flow and includes the moving face. A zero or
    /// negative driving-point response is refused, never clipped or replaced
    /// by a delayed boundary. No state is mutated on success or refusal.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn step(
        &self,
        opening_m: f64,
        velocity_m_s: f64,
        mouth_pressure_pa: f64,
        free_pressure_pa: f64,
        mean_impedance: f64,
    ) -> Result<LoadedReedStep, AcousticRealizeError> {
        if ![
            opening_m,
            velocity_m_s,
            mouth_pressure_pa,
            free_pressure_pa,
            mean_impedance,
        ]
        .iter()
        .all(|x| x.is_finite())
            || mean_impedance <= 0.0
        {
            return Err(refuse(
                "loaded reed requires finite state and a positive mean pressure response",
            ));
        }
        let contact = AffineContactStep::new(
            &self.contact.obstacles()[0],
            opening_m,
            crate::unilateral_contact::distributed::MAX_APERTURE_CONTACT_POINTS,
        )
        .map_err(|error| AcousticRealizeError::Nonlinear(error.to_string()))?;
        let incoming = 0.5 * free_pressure_pa;
        let next = super::contact::step_with_contact_response(
            self.reed,
            self.rho,
            mean_impedance,
            incoming,
            mouth_pressure_pa,
            opening_m,
            velocity_m_s,
            self.dt,
            0.0,
            |gap| Ok(gap.max(0.0)),
            None,
            OpeningPropagator::Exact {
                blade: &self.blade,
                external_n: 0.0,
                memory: &[],
            },
            |gap, swept| {
                moving::characteristic_pressure(
                    self.reed.width_m,
                    gap,
                    self.rho,
                    mean_impedance,
                    incoming,
                    mouth_pressure_pa,
                    swept,
                )
            },
            |after, velocity| {
                contact
                    .response(after, velocity)
                    .map_err(|error| AcousticRealizeError::Nonlinear(error.to_string()))
            },
        )?;
        let pressure_pa = next.outgoing + incoming;
        let drop = mouth_pressure_pa - pressure_pa;
        let jet_flow_m3_s = moving::volume_flow(
            self.reed.width_m,
            f64::midpoint(opening_m, next.opening).max(0.0),
            drop,
            self.rho,
        );
        let swept = reed_swept_flow(self.face, next.mean_velocity);
        // Reconstruct from physical contributions, avoiding division of a
        // pressure difference by a possibly tiny mean body impedance.
        let volume_flow_m3_s = jet_flow_m3_s + swept;
        let response = contact
            .response(next.opening, next.mean_velocity)
            .map_err(|error| AcousticRealizeError::Nonlinear(error.to_string()))?;
        let start = BladeState {
            displacement_m: opening_m - self.reed.rest_opening_m,
            velocity_m_s,
        };
        let held_force = response.force - self.face * drop;
        let viscous_heat_j = self.blade.viscous_loss_j(start, held_force, &[]);
        let contact_heat_j = self.dt * response.dissipated_power;
        let jet_heat_j = self.dt * drop * jet_flow_m3_s;
        let mouth_work_j = self.dt * mouth_pressure_pa * volume_flow_m3_s;
        let bore_work_j = self.dt * pressure_pa * volume_flow_m3_s;
        let energy_before_j = self.energy(opening_m, velocity_m_s);
        let energy_after_j = self.energy(next.opening, next.velocity);
        // Stable quadratic endpoint change plus the original contact energy.
        let delta_h_j = self.reed.mass_kg
            * (next.velocity - velocity_m_s)
            * f64::midpoint(next.velocity, velocity_m_s)
            + self.stiffness
                * (next.opening - opening_m)
                * (f64::midpoint(next.opening, opening_m) - self.reed.rest_opening_m)
            + self.contact.probe(&[next.opening, 0.0]).contact_energy
            - self.contact.probe(&[opening_m, 0.0]).contact_energy;
        let residual_j =
            delta_h_j + viscous_heat_j + contact_heat_j + jet_heat_j + bore_work_j - mouth_work_j;
        let scale = energy_before_j
            + energy_after_j
            + viscous_heat_j.abs()
            + contact_heat_j.abs()
            + jet_heat_j.abs()
            + mouth_work_j.abs()
            + bore_work_j.abs();
        let boundary_residual =
            pressure_pa - (free_pressure_pa + mean_impedance * volume_flow_m3_s);
        let pressure_scale = pressure_pa.abs()
            + free_pressure_pa.abs()
            + mouth_pressure_pa.abs()
            + mean_impedance * (jet_flow_m3_s.abs() + swept.abs());
        if ![
            pressure_pa,
            volume_flow_m3_s,
            jet_flow_m3_s,
            energy_before_j,
            energy_after_j,
            delta_h_j,
            viscous_heat_j,
            contact_heat_j,
            jet_heat_j,
            mouth_work_j,
            bore_work_j,
            residual_j,
            scale,
            pressure_scale,
            boundary_residual,
        ]
        .iter()
        .all(|x| x.is_finite())
            || viscous_heat_j < 0.0
            || contact_heat_j < 0.0
            || jet_heat_j < 0.0
            || residual_j.abs() > 2e-11 * scale.max(f64::MIN_POSITIVE)
            || boundary_residual.abs() > 2e-11 * pressure_scale.max(f64::MIN_POSITIVE)
        {
            return Err(refuse(
                "loaded reed physical work or current-window flow balance refused",
            ));
        }
        Ok(LoadedReedStep {
            opening_m: next.opening,
            velocity_m_s: next.velocity,
            mean_velocity_m_s: next.mean_velocity,
            pressure_pa,
            volume_flow_m3_s,
            jet_flow_m3_s,
            account: LoadedReedAccount {
                energy_before_j,
                energy_after_j,
                delta_h_j,
                viscous_heat_j,
                contact_heat_j,
                jet_heat_j,
                mouth_work_j,
                bore_work_j,
                residual_j,
            },
        })
    }
}

#[cfg(test)]
mod tests;
