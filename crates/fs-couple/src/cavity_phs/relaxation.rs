//! Finite-rate water compression in the native lumped Helmholtz network.
//!
//! The retained material supplies the original density, equilibrium modulus,
//! positive compression arm and caller-declared relaxation time. The generic
//! Helmholtz owner supplies neck inertia, end corrections and reciprocal
//! pressure/volume-flow ports. Its ordinary quadratic propagator advances the
//! neck and compression memory together; this adapter adds no integrator.
//!
//! State is [q, pi, z]: excess compression volume [m^3], neck momentum [Pa s],
//! and energy-normalized relaxed compression [sqrt(J)]. With fixed cavity
//! volume V, q_dot = U_in - U_neck_out, delta_mass = rho0*q and
//! delta_rho = rho0*q/V. These are linearized open-cavity continuity variables;
//! q is not a moving-wall volume for a fixed-mass nonlinear EOS.
//!
//! The cavity and neck are homogeneous about one original coexistence state.
//! This small-signal, spatially uniform model requires wavelengths large
//! compared with the cavity and neck. The caller owns the relaxation-rate,
//! phase-fraction and homogenization domain. No phase inventory, thermal state,
//! viscosity, radiation resistance or additional wall loss is inferred.
//! The explicit neck resistance represents a separate admitted passive load;
//! both declared end corrections remain in the original neck inertia.

use crate::cavity_modal::neck_ends;
use fs_material::phase::water::acoustics::relaxation::RelaxingRegionalWater;
use fs_math::det;
use fs_phs::linear_exact::{LinearSplitBudget, PreparedExactZoh};
use fs_phs::{PhsError, PortHamiltonian, QuadraticStorage, RelaxationBranch};
use fs_scenario::HelmholtzCavity;

/// Original material and geometry beside their three-state reciprocal network.
pub struct RelaxingWaterCavity {
    material: RelaxingRegionalWater,
    geometry: HelmholtzCavity,
    system: PortHamiltonian,
    neck_inertance_pa_s2_m3: f64,
    neck_resistance_pa_s_m3: f64,
    memory_scale: f64,
}

/// Read-only physical observations of the accepted native state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RelaxingWaterCavityObservation {
    pub pressure_increment_pa: f64,
    /// Positive flow leaves the cavity through its neck.
    pub neck_outflow_m3_s: f64,
    /// Retained separately even when adding it to the reference mass rounds away.
    pub excess_mass_kg: f64,
    /// First-order density response at fixed cavity volume.
    pub density_increment_kg_m3: f64,
    /// Reference density plus its first-order perturbation.
    pub linearized_density_kg_m3: f64,
    /// Complete compression, neck kinetic and material-memory storage.
    pub stored_energy_j: f64,
    /// Instantaneous loss already present in the native pHS resistance.
    pub compression_dissipation_w: f64,
    /// Instantaneous loss of the explicit neck load, counted once.
    pub neck_dissipation_w: f64,
}

fn invalid(what: &'static str) -> PhsError {
    PhsError::Unadmitted { what }
}

impl RelaxingWaterCavity {
    /// Bind the original thermal density and bulk endpoints to explicit geometry.
    ///
    /// The material's relaxation time is unchanged. Resistance is acoustic
    /// pressure/volume-flow resistance [Pa s/m^3]; zero gives an unloaded neck.
    /// Native geometry, positivity and representability admissions remain active.
    pub fn try_new(
        material: RelaxingRegionalWater,
        geometry: HelmholtzCavity,
        neck_resistance_pa_s_m3: f64,
    ) -> Result<Self, PhsError> {
        let acoustic = material.equilibrium().state();
        let native = fs_phs::helmholtz_resonator_flow(
            geometry.volume_m3,
            geometry.neck_radius_m,
            geometry.neck_length_m,
            acoustic.density_kg_m3(),
            acoustic.sound_speed_m_s(),
            neck_resistance_pa_s_m3,
            neck_ends(&geometry),
        )?;
        // Read the native inertia chart; do not recalculate its geometry.
        let inverse_mass = match native.storage().exact_quadratic_matrix() {
            Some([_, _, _, value]) if value.is_finite() && *value > 0.0 => *value,
            _ => return Err(invalid("water cavity neck inertia is unrepresentable")),
        };
        let neck_inertance_pa_s2_m3 = inverse_mass.recip();
        let equilibrium = acoustic.adiabatic_bulk_modulus_pa() / geometry.volume_m3;
        let strength = material.arm().delta_bulk_modulus_pa / geometry.volume_m3;
        let reference_mass = acoustic.density_kg_m3() * geometry.volume_m3;
        if ![
            neck_inertance_pa_s2_m3,
            equilibrium,
            strength,
            reference_mass,
        ]
        .into_iter()
        .all(|value| value.is_finite() && value > 0.0)
        {
            return Err(invalid(
                "water cavity storage or reference mass is unrepresentable",
            ));
        }
        // Keep the original source modulus directly, including its rounding.
        // The native generic constructor otherwise forms it again from rho*c^2.
        let (j, r, g) = native.structure();
        let system = PortHamiltonian::new(
            2,
            1,
            j.to_vec(),
            r.to_vec(),
            g.to_vec(),
            Box::new(QuadraticStorage::new(
                vec![equilibrium, 0.0, 0.0, inverse_mass],
                2,
            )?),
        )?
        .with_relaxation_branches(vec![RelaxationBranch {
            projection: vec![1.0, 0.0],
            stiffness: strength,
            relaxation_time_s: material.arm().relaxation_time_s,
        }])?;
        Ok(Self {
            material,
            geometry,
            system,
            neck_inertance_pa_s2_m3,
            neck_resistance_pa_s_m3,
            memory_scale: det::sqrt(strength),
        })
    }

    /// Complete original water source and its authored compression timescale.
    pub const fn material(&self) -> &RelaxingRegionalWater {
        &self.material
    }

    pub const fn geometry(&self) -> &HelmholtzCavity {
        &self.geometry
    }

    /// Native density-derived neck inertia, including both original end corrections.
    pub const fn neck_inertance_pa_s2_m3(&self) -> f64 {
        self.neck_inertance_pa_s2_m3
    }

    /// One input is volume flow into the cavity; its output is pressure.
    /// This is the original pHS consumed by native interconnections and solvers.
    pub fn system(&self) -> &PortHamiltonian {
        &self.system
    }

    /// Move the original material and its network together into an interconnection.
    /// Native transformers preserve pressure/flow and mechanical-force/velocity power.
    pub fn into_parts(self) -> (RelaxingRegionalWater, PortHamiltonian) {
        (self.material, self.system)
    }

    /// Prepare the complete fixed-clock native propagator once. Its step_into
    /// reuses scratch and returns integrated input work and positive dissipation.
    pub fn prepare(
        &self,
        duration_s: f64,
        budget: LinearSplitBudget,
    ) -> Result<PreparedExactZoh<'_>, PhsError> {
        PreparedExactZoh::new(&self.system, duration_s, budget)
    }

    /// Form a complete native state. None initializes relaxed memory at q;
    /// Some(0.0) represents an instantaneous compression with frozen history.
    /// An explicit memory value uses the existing energy-normalized z coordinate.
    pub fn initial_state(
        &self,
        compression_volume_m3: f64,
        neck_outflow_m3_s: f64,
        memory_sqrt_j: Option<f64>,
    ) -> Result<[f64; 3], PhsError> {
        let state = [
            compression_volume_m3,
            -self.neck_inertance_pa_s2_m3 * neck_outflow_m3_s,
            memory_sqrt_j.unwrap_or(self.memory_scale * compression_volume_m3),
        ];
        self.observe(&state)?;
        Ok(state)
    }

    /// Observe the same storage gradient used by the native dynamics.
    ///
    /// The reported density is a linear acoustic perturbation. Positivity here
    /// is a necessary endpoint check, not a nonlinear EOS or a whole-trajectory
    /// phase-domain certificate. No observation mutates the source or its history.
    pub fn observe(&self, state: &[f64]) -> Result<RelaxingWaterCavityObservation, PhsError> {
        if state.len() != 3 || state.iter().any(|value| !value.is_finite()) {
            return Err(invalid(
                "water cavity observation needs its complete finite state",
            ));
        }
        let mut effort = [0.0; 3];
        self.system.storage().gradient(state, &mut effort);
        let rho = self.material.source().bulk_density_kg_m3();
        let density_increment_kg_m3 = rho * (state[0] / self.geometry.volume_m3);
        let result = RelaxingWaterCavityObservation {
            pressure_increment_pa: effort[0],
            neck_outflow_m3_s: -effort[1],
            excess_mass_kg: rho * state[0],
            density_increment_kg_m3,
            linearized_density_kg_m3: rho + density_increment_kg_m3,
            stored_energy_j: self.system.hamiltonian(state),
            compression_dissipation_w: effort[2] * effort[2]
                / self.material.arm().relaxation_time_s,
            neck_dissipation_w: self.neck_resistance_pa_s_m3 * effort[1] * effort[1],
        };
        if ![
            result.pressure_increment_pa,
            result.neck_outflow_m3_s,
            result.excess_mass_kg,
            result.density_increment_kg_m3,
            result.linearized_density_kg_m3,
            result.stored_energy_j,
            result.compression_dissipation_w,
            result.neck_dissipation_w,
        ]
        .into_iter()
        .all(f64::is_finite)
            || result.linearized_density_kg_m3 <= 0.0
            || result.stored_energy_j < 0.0
            || result.compression_dissipation_w < 0.0
            || result.neck_dissipation_w < 0.0
        {
            return Err(invalid("water cavity observation is unrepresentable"));
        }
        Ok(result)
    }
}

#[cfg(test)]
mod tests;
