//! Curl of the original conservative viscous acceleration.
//!
//! Evaluate exactly the configured molecular/SGS face stress, then apply the
//! owner's finite-volume divergence, divide by each cell's accepted density,
//! and differentiate that acceleration. Variable viscosity, bulk stress,
//! density gradients and physical wall/bath tractions are retained. Replacing
//! this with nu*laplacian(omega) is valid only under narrower assumptions.
//!
//! This observation never applies an additional stress or heat source. The
//! stress adds to the inviscid material-rate terms; advection is already
//! included in D(omega)/Dt. Curl reconstruction has the same smooth-field,
//! resolution and nonperiodic-interior restrictions as the native inviscid
//! observation. It is not a wall-sheet or shock derivative.

use super::super::{
    CartesianGrid, Error, cartesian, evolve,
    vorticity::{self as inviscid, VorticityBudget3d},
};
use super::{
    Euler3d, MolecularError, MolecularTransport3d, TransportLaw, allocated, poll, spatial,
};

/// Resolved inviscid production and the configured stress contribution.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct StressVorticityBudget3d {
    /// Original EOS, velocity-gradient and pressure/density budget.
    pub inviscid: VorticityBudget3d,
    /// curl((div(tau))/rho) [1/s^2] from the actual conservative face stress.
    ///
    /// Includes any explicitly selected SGS and wall traction model. Molecular
    /// and modeled stress are not falsely reported as separately observed terms.
    pub stress: [f64; 3],
}

impl StressVorticityBudget3d {
    /// Resolved inviscid plus configured stress contribution to D(omega)/Dt.
    ///
    /// Imposed body-force curl is excluded. Fourier conduction changes later
    /// thermodynamics but contributes no separate instantaneous momentum curl.
    pub fn resolved_material_rate(self) -> [f64; 3] {
        let rate = self.inviscid.inviscid_material_rate();
        std::array::from_fn(|a| rate[a] + self.stress[a])
    }
}

/// Additional reusable scratch for native stress acceleration and its curl.
///
/// Source primitives, velocity gradients and face fluxes reuse the prepared
/// molecular operator's existing storage. Only a completely evaluated output
/// can be borrowed; failed queries may change private scratch.
pub struct StressVorticityWorkspace3d {
    grid: CartesianGrid,
    acceleration: Vec<[f64; 3]>,
    output: Vec<StressVorticityBudget3d>,
}

impl StressVorticityWorkspace3d {
    /// Allocated numerical capacity, excluding fixed metadata/allocator overhead
    /// and the molecular operator's already prepared source/face workspaces.
    pub fn workspace_bytes(&self) -> usize {
        self.acceleration.capacity() * core::mem::size_of::<[f64; 3]>()
            + self.output.capacity() * core::mem::size_of::<StressVorticityBudget3d>()
    }
}

impl<L: TransportLaw> MolecularTransport3d<L> {
    /// Required additional numerical bytes for acceleration and combined output.
    pub fn required_vorticity_workspace_bytes(&self) -> Result<usize, MolecularError> {
        self.grid
            .len()
            .checked_mul(
                core::mem::size_of::<[f64; 3]>() + core::mem::size_of::<StressVorticityBudget3d>(),
            )
            .ok_or(MolecularError::Invalid(
                "stress vorticity workspace size overflow",
            ))
    }

    /// Prepare a reusable native stress/vorticity observation within a byte cap.
    pub fn vorticity_workspace(
        &self,
        maximum_bytes: usize,
    ) -> Result<StressVorticityWorkspace3d, MolecularError> {
        inviscid::admit_axes::<Error>(&self.grid, self.boundaries)?;
        if self.required_vorticity_workspace_bytes()? > maximum_bytes {
            return Err(MolecularError::Invalid(
                "stress vorticity workspace exceeds its byte budget",
            ));
        }
        let work = StressVorticityWorkspace3d {
            grid: self.grid,
            acceleration: allocated(self.grid.len(), [0.0; 3])?,
            output: allocated(self.grid.len(), StressVorticityBudget3d::default())?,
        };
        if work.workspace_bytes() > maximum_bytes {
            return Err(MolecularError::Invalid(
                "stress vorticity workspace capacity exceeds its byte budget",
            ));
        }
        Ok(work)
    }

    /// Observe the accepted field through the original complete stress operator.
    ///
    /// The configured molecular law, optional SGS, no-slip/Navier-slip walls,
    /// and enabled held reservoir exchange supply exactly the face stresses
    /// used by advancement. Each cell's original EOS is decoded once. Held
    /// accepted reservoir traces are used; no time-dependent provider is sampled.
    ///
    /// Queries allocate no field storage and leave accepted cells, physical
    /// time and published interval loads unchanged. Transport-law and closure
    /// evaluation may update their reusable caches. Cancellation/source refusal
    /// publishes no partial result; the workspace can be reused after failure.
    pub fn vorticity<'a>(
        &mut self,
        flow: &Euler3d,
        work: &'a mut StressVorticityWorkspace3d,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<&'a [StressVorticityBudget3d], MolecularError> {
        poll(&mut cancelled)?;
        inviscid::admit_axes::<Error>(&self.grid, self.boundaries)?;
        if self.gas != flow.gas || self.grid != flow.grid || self.boundaries != flow.boundaries {
            return Err(MolecularError::Invalid(
                "molecular vorticity operator and flow differ",
            ));
        }
        if work.grid != self.grid {
            return Err(MolecularError::Invalid(
                "stress vorticity workspace and physical grid differ",
            ));
        }
        spatial(
            &self.gas,
            &self.grid,
            self.boundaries,
            &self.walls,
            self.reservoir_exchange.then(|| flow.reservoirs.accepted()),
            &flow.states,
            &mut self.law,
            &mut self.work.primitive,
            &mut self.work.gradient,
            &mut self.work.viscosity,
            &mut self.work.flux,
            &mut self.work.binary,
            #[cfg(feature = "experimental-sgs")]
            self.sgs.as_mut(),
            &mut cancelled,
        )?;
        // Molecular mass flux is zero, so du/dt = (d(rho*u)/dt)/rho.
        // Use the original divergence and sign convention, not a second
        // constitutive reconstruction or an increment of the accepted field.
        cartesian::divergence::<3, Error, _>(
            &self.grid,
            3,
            1.0,
            |axis, face, component| {
                evolve::flux_array(self.work.flux[axis][face], axis)[component + 1]
            },
            &mut cancelled,
            |cell, rate, _| {
                let acceleration: [f64; 3] =
                    std::array::from_fn(|d| rate[d] / flow.states[cell].density);
                if !acceleration.into_iter().all(f64::is_finite) {
                    return Err(Error::Invalid(
                        "unrepresentable physical stress acceleration",
                    ));
                }
                work.acceleration[cell] = acceleration;
                Ok(())
            },
        )?;
        for (cell, output) in work.output.iter_mut().enumerate() {
            if cell % 64 == 0 {
                poll(&mut cancelled)?;
            }
            let column = |axis| {
                inviscid::derivative(&self.grid, self.boundaries, cell, axis, |i| {
                    work.acceleration[i]
                })
            };
            let [x, y, z] = [column(0), column(1), column(2)];
            let result = StressVorticityBudget3d {
                inviscid: inviscid::budget::<Error>(
                    &self.grid,
                    self.boundaries,
                    &self.work.primitive,
                    cell,
                )?,
                stress: [y[2] - z[1], z[0] - x[2], x[1] - y[0]],
            };
            if !result
                .stress
                .into_iter()
                .chain(result.resolved_material_rate())
                .all(f64::is_finite)
            {
                return Err(MolecularError::Invalid(
                    "unrepresentable physical stress vorticity rate",
                ));
            }
            *output = result;
        }
        poll(&mut cancelled)?;
        Ok(&work.output)
    }
}
