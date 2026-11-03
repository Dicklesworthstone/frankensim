//! Complete small-vibration pencil of the SAME sealed-gas equilibrium.
//!
//! The finite gas supplies its original local isentropic stiffness, including
//! temperature-dependent calorics or real-gas departures. Reference inertia,
//! material admission, eigensolve and work-conjugate projection remain with the
//! original membrane/modal owners. No extra gas spring, mass or damping is fit.

use super::*;
use fs_solid::membrane::vibration::{MembraneModes, VibrationSettings};

/// Accepted coupled equilibrium and all positive mass-normalized discrete modes.
///
/// The pencil already contains the gas compressibility and geometric pressure
/// tangent. Modal pressure is an observation of that same gas; applying it again
/// through the pressure port would double-count its restoring force. The port
/// remains available for an independent incremental external pressure load.
///
/// This local homogeneous-gas model omits gas inertia, acoustic wave propagation,
/// thermal exchange and damping. It does not advance a finite-amplitude gas or
/// membrane trajectory or establish continuum stability/mesh convergence.
#[derive(Debug)]
pub struct GasMembraneVibration<'a> {
    solution: GasMembraneSolution,
    basis: MembraneModes<'a>,
}

impl<'a> GasMembrane<'a> {
    /// Solve this original coupled equilibrium, then lower its COMPLETE tangent
    /// through the original membrane reference mass and dense modal owner.
    ///
    /// K = K_membrane - (p-pe) Hess(V) - (dp/dV) grad(V) grad(V)^T.
    /// The single accepted native gas state supplies every tangent column; the
    /// modal lowering performs no additional gas inversion. The full positive
    /// discrete spectrum, source-sheet linearization and native modal residual
    /// gates apply. `settings` bounds the dense leaf to at most 128 free DOFs.
    ///
    /// The accepted volume and temperature must lie strictly inside the declared
    /// volume interval and entropy-solve bracket. As for the original isentrope,
    /// the caller must supply a continuous homogeneous source branch around the
    /// equilibrium. These local checks do not certify an unqueried source domain,
    /// cross a caloric discontinuity, or establish global phase stability.
    pub fn solve_vibration(
        &self,
        dead_load_n: &[[f64; 3]],
        initial_displacement_m: &[[f64; 3]],
        solve: GasMembraneSolve,
        settings: VibrationSettings,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<GasMembraneVibration<'a>, GasMembraneError> {
        poll(&mut cancelled)?;
        let n = self.free_dofs.len();
        settings.validate(n)?;
        let GasMembraneSolution {
            equilibrium,
            newton,
            evaluations,
            gas_trials,
        } = self.solve(dead_load_n, initial_displacement_m, solve, &mut cancelled)?;
        let volume = equilibrium.gas.volume_m3();
        let temperature = match &equilibrium.gas {
            MembraneGas::Ideal(gas) => gas.state().state().temperature,
            MembraneGas::PengRobinson(gas) => gas.state().temperature_k(),
        };
        if volume <= self.volume_interval_m3[0]
            || volume >= self.volume_interval_m3[1]
            || temperature <= solve.gas.minimum_temperature_k
            || temperature >= solve.gas.maximum_temperature_k
        {
            return Err(GasMembraneError::Invalid(
                "two-sided gas modes need equilibrium strictly inside volume and temperature brackets",
            ));
        }
        let patch: &'a MembranePatch = self.patch;
        let gradient_all =
            patch.pressure_load_n(&equilibrium.displacement_m, 1.0, &mut cancelled)?;
        let gradient = self
            .free_dofs
            .iter()
            .map(|&i| gradient_all[i / 3][i % 3])
            .collect();
        // Move the accepted state into one immutable trial. Do not reconstruct
        // the isentrope for each column or linearize a different physical state.
        let trial = Trial {
            value: equilibrium,
            gradient,
        };
        let mut stiffness = vec![0.0; n * n];
        let mut direction = vec![0.0; n];
        for column in 0..n {
            poll(&mut cancelled)?;
            direction[column] = 1.0;
            let action = self.tangent_at(&trial, &direction, &mut cancelled)?;
            direction[column] = 0.0;
            for (row, value) in action.into_iter().enumerate() {
                stiffness[row * n + column] = value;
            }
        }
        let basis = patch.modes_from_tangent(
            &trial.value.displacement_m,
            &self.free_dofs,
            &stiffness,
            settings,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        Ok(GasMembraneVibration {
            solution: GasMembraneSolution {
                equilibrium: trial.value,
                newton,
                evaluations,
                gas_trials,
            },
            basis,
        })
    }
}

impl<'a> GasMembraneVibration<'a> {
    /// Same accepted source-bearing gas, shape, potential and nonlinear report.
    /// Evaluation/gas-trial counts describe the equilibrium solve; its final
    /// gas state is reused without inversion during modal assembly.
    pub fn solution(&self) -> &GasMembraneSolution {
        &self.solution
    }

    /// Original reference mass, full spectrum, nodal projection/reconstruction
    /// and external incremental pressure/volume port at the accepted geometry.
    pub fn basis(&self) -> &MembraneModes<'a> {
        &self.basis
    }

    /// Original local adiabatic dp/dV [Pa/m3] of the accepted native gas.
    pub fn pressure_volume_derivative_pa_m3(&self) -> f64 {
        self.solution
            .equilibrium
            .pressure
            .pressure_volume_derivative_pa_m3
    }

    /// First-order enclosed-volume change [m3] for q [m sqrt(kg)]. This is a
    /// linear observation, not a new finite-volume parcel or geometry proposal.
    pub fn volume_perturbation_m3(&self, coordinates: &[f64]) -> Result<f64, GasMembraneError> {
        let coupling = self.basis.volume_coupling();
        if coordinates.len() != coupling.len() || coordinates.iter().any(|q| !q.is_finite()) {
            return Err(GasMembraneError::Invalid(
                "one finite value per gas/membrane modal coordinate required",
            ));
        }
        let (mut volume, mut correction) = (0.0, 0.0);
        for (&b, &q) in coupling.iter().zip(coordinates) {
            let adjusted = b * q - correction;
            let next = volume + adjusted;
            correction = (next - volume) - adjusted;
            volume = next;
        }
        if !volume.is_finite() {
            return Err(GasMembraneError::Invalid(
                "modal swept volume is unrepresentable",
            ));
        }
        Ok(volume)
    }

    /// First-order absolute-gas-pressure perturbation [Pa]: (dp/dV) B.q.
    /// Exterior pressure is held constant, so this is also the pressure-jump
    /// perturbation. Its restoring force is ALREADY in the eigenvalue pencil.
    pub fn pressure_perturbation_pa(&self, coordinates: &[f64]) -> Result<f64, GasMembraneError> {
        let pressure =
            self.pressure_volume_derivative_pa_m3() * self.volume_perturbation_m3(coordinates)?;
        if !pressure.is_finite() {
            return Err(GasMembraneError::Invalid(
                "modal gas pressure is unrepresentable",
            ));
        }
        Ok(pressure)
    }
}
