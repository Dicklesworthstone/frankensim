//! Two-sided shell motion and sound in the original molecular exterior.
//!
//! The same complete K/C, exact facet-local P1 work map and mechanical LU are
//! used by the nondispersive owner. The gas changes H(k(omega)); the acoustic
//! drive remains i*rho*omega*G^T*v. Pressure work includes molecular heat in
//! the exterior and is not mislabeled as distant radiated flux.
use super::*;
use crate::bem_loaded_modes::MolecularModesPower;
use fs_bem::helmholtz::molecular::MolecularExterior;
use fs_bem::reciprocal_radiation::open_screen::{
    MolecularScreenPortLimits, MolecularScreenVelocityPortSolution,
    PreparedMolecularScreenReceivers,
};

/// One original boundary batch and complete mechanical response admission.
/// Receiver quadrature has its own explicit row/work/depth limits.
#[derive(Debug, Clone, Copy)]
pub struct MolecularOpenScreenHarmonicOptions {
    pub boundary: MolecularScreenPortLimits,
    pub maximum_load_cases: usize,
    /// Same live mechanical payload as OpenScreenHarmonicOptions, including
    /// result conversion. Original source and boundary payloads are separate.
    pub maximum_response_bytes: usize,
    /// 32*modes^3 + 256*(cases+1)*modes^2, separate from boundary work.
    pub maximum_work: u64,
    pub physical_residual_tolerance: f64,
    pub power_relative_tolerance: f64,
}

impl MolecularOpenScreenHarmonicOptions {
    fn mechanical(self) -> OpenScreenHarmonicOptions {
        OpenScreenHarmonicOptions {
            boundary: self.boundary.ports,
            maximum_load_cases: self.maximum_load_cases,
            maximum_response_bytes: self.maximum_response_bytes,
            maximum_work: self.maximum_work,
            physical_residual_tolerance: self.physical_residual_tolerance,
            power_relative_tolerance: self.power_relative_tolerance,
        }
    }
}

/// Actual force-scaled shell motion in its original mass-normalized basis.
#[derive(Debug, Clone)]
pub struct MolecularOpenScreenHarmonicResponse {
    pub angular_frequency_rad_s: f64,
    pub modal_displacement: Vec<C64>,
    pub modal_velocity: Vec<C64>,
    /// Fluid force ON the shell, -Z*v, in original mode order.
    pub fluid_reaction: Vec<C64>,
    /// Kinetic and elastic shell energy only, with no inferred fluid storage.
    pub mean_structural_energy_j: f64,
    pub mechanical_residual: f64,
    pub maximum_boundary_residual: f64,
    /// Input = original viscosity + signed total exterior pressure-jump work.
    /// This screen has no separate sheet loss or mechanical Maxwell arms, so
    /// wall_w and relaxation_w are zero; gas memory is in exterior_w once.
    pub power: MolecularModesPower,
}

/// One complete molecular boundary bank and one mechanical factorization for
/// all authored force cases. No source modes or mutual loads are discarded.
pub struct MolecularOpenScreenHarmonicBatch<'p, 's> {
    radiation: MolecularScreenVelocityPortSolution<'p, 's>,
    responses: Vec<MolecularOpenScreenHarmonicResponse>,
}

/// Adaptive pressure/particle-velocity rows bound to this loaded source batch.
/// Every case uses its actual solved modal velocity without a new quadrature.
pub struct PreparedMolecularOpenScreenReceivers<'b, 'p, 's> {
    receiver: PreparedMolecularScreenReceivers<'b, 'p, 's>,
    responses: &'b [MolecularOpenScreenHarmonicResponse],
}

impl OpenScreenLoadedModes {
    /// Bind the original shell mechanics using the supplied evaluated gas's
    /// equilibrium density/speed. The frequency response below retains the
    /// complete molecular law; no effective complex sound speed is fabricated.
    pub fn from_shell_molecular(
        shell: &ShellReduction,
        exterior: &MolecularExterior,
        viscosity: Option<&ModalPlateViscosity>,
        options: OpenScreenModesOptions,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, BemModesError> {
        let gas = exterior.source().gas();
        Self::from_shell(
            shell,
            Medium {
                density: gas.density,
                sound_speed: gas.sound_speed,
            },
            viscosity,
            options,
            cancelled,
        )
    }
}

impl<'s> PreparedOpenScreenModes<'s> {
    /// Solve the original source against H(k(omega)) and real i*rho*omega
    /// forcing. Exact P1 G supplies BOTH the drive and the resisting effort.
    /// All force cases share the complete self/mutual Z and mechanical LU.
    /// The supplied gas must match this source's equilibrium density/speed.
    ///
    /// This is a frozen homogeneous, inviscid-momentum molecular exterior on
    /// both faces. It is not a cavity-backed or finite-thickness shell, a CQ
    /// reaction realization, a nonlinear trajectory, or a loaded pole solve.
    /// Original mesh/angular/growth bounds and signed work remain unrepaired.
    pub fn harmonic_molecular_response_batch<'p>(
        &'p self,
        exterior: &MolecularExterior,
        omega: f64,
        forces: &[&[C64]],
        options: MolecularOpenScreenHarmonicOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularOpenScreenHarmonicBatch<'p, 's>, BemModesError> {
        poll(&mut cancelled)?;
        let mechanical = options.mechanical();
        self.source
            .admit_harmonic_batch(omega, forces, mechanical)?;
        let gas = exterior.source().gas();
        if gas.density.to_bits() != self.source.medium.density.to_bits()
            || gas.sound_speed.to_bits() != self.source.medium.sound_speed.to_bits()
        {
            return Err(bad(
                "molecular open-screen gas differs from the source equilibrium medium",
            ));
        }
        let rows: Vec<_> = self.source.weights.iter().map(Vec::as_slice).collect();
        let radiation = self
            .screen
            .solve_molecular_velocity_ports(
                exterior,
                omega,
                &rows,
                options.boundary,
                &mut cancelled,
            )
            .map_err(BemModesError::MolecularScreen)?;
        let responses = self.source.solve_harmonic_batch(
            omega,
            forces,
            mechanical,
            radiation.impedance(),
            radiation.maximum_backward_residual(),
            |velocity| radiation.resisting_efforts(velocity),
            &mut cancelled,
        )?;
        let responses = responses
            .into_iter()
            .map(|r| MolecularOpenScreenHarmonicResponse {
                angular_frequency_rad_s: r.angular_frequency_rad_s,
                modal_displacement: r.modal_displacement,
                modal_velocity: r.modal_velocity,
                fluid_reaction: r.fluid_reaction,
                mean_structural_energy_j: r.mean_structural_energy_j,
                mechanical_residual: r.mechanical_residual,
                maximum_boundary_residual: r.maximum_boundary_residual,
                power: MolecularModesPower {
                    input_w: r.power.input_w,
                    kelvin_w: r.power.kelvin_w,
                    relaxation_w: 0.0,
                    exterior_w: r.power.radiated_w,
                    wall_w: 0.0,
                    balance_w: r.power.balance_w,
                    apparent_power_scale_w: r.power.apparent_power_scale_w,
                    relative_balance: r.power.relative_balance,
                },
            })
            .collect();
        poll(&mut cancelled)?;
        Ok(MolecularOpenScreenHarmonicBatch {
            radiation,
            responses,
        })
    }
}

impl<'p, 's> MolecularOpenScreenHarmonicBatch<'p, 's> {
    pub fn responses(&self) -> &[MolecularOpenScreenHarmonicResponse] {
        &self.responses
    }
    pub fn radiation(&self) -> &MolecularScreenVelocityPortSolution<'p, 's> {
        &self.radiation
    }
    pub fn pressure_jump(
        &self,
        case: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<C64>, BemModesError> {
        poll(&mut cancelled)?;
        let response = self
            .responses
            .get(case)
            .ok_or_else(|| bad("molecular open-screen load-case index is out of range"))?;
        self.radiation
            .pressure_jump(&response.modal_velocity, cancelled)
            .map_err(BemModesError::ScreenBoundary)
    }
    pub fn prepare_receivers<'b>(
        &'b self,
        points_m: &[[f64; 3]],
        separation_m: f64,
        options: near_field::Options,
        maximum_row_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedMolecularOpenScreenReceivers<'b, 'p, 's>, BemModesError> {
        poll(&mut cancelled)?;
        let receiver = self
            .radiation
            .prepare_receivers(
                points_m,
                separation_m,
                options,
                maximum_row_bytes,
                &mut cancelled,
            )
            .map_err(BemModesError::ScreenBoundary)?;
        Ok(PreparedMolecularOpenScreenReceivers {
            receiver,
            responses: &self.responses,
        })
    }
}

impl PreparedMolecularOpenScreenReceivers<'_, '_, '_> {
    pub fn points_m(&self) -> &[[f64; 3]] {
        self.receiver.points_m()
    }
    pub fn kernel_evaluations(&self) -> usize {
        self.receiver.kernel_evaluations()
    }
    /// Pressure [Pa] and particle velocity [m/s] on either side, including
    /// their absolute quadrature estimates for the fixed solved P1 jump.
    /// Those estimates do not bound boundary discretization or material error.
    pub fn evaluate(
        &self,
        case: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<near_field::VectorObservation, BemModesError> {
        poll(&mut cancelled)?;
        let response = self
            .responses
            .get(case)
            .ok_or_else(|| bad("molecular open-screen load-case index is out of range"))?;
        self.receiver
            .evaluate(&response.modal_velocity, cancelled)
            .map_err(BemModesError::ScreenBoundary)
    }
}

#[cfg(test)]
mod tests;
