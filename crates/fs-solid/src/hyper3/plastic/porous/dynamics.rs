//! Finite porous motion through the original momentum and compaction owners.
//!
//! The consistent P1 reference mass and backward-Euler tangent enter the same
//! native finite-element Newton solve as the porous incremental potential.
//! Every trial starts from the accepted non-isochoric plastic metric. A complete
//! endpoint retains displacement, velocity, plastic volume, hardening and mass
//! only after all material, momentum and work gates pass.
//!
//! Original mass is independent of plastic volume: current density is rho0/J,
//! never rho0/Jp. Physical plastic dissipation, signed material endpoint-work
//! defect, pressure endpoint-work correction and backward-Euler numerical
//! kinetic loss have separate accounts and budgets. Numerical kinetic loss is
//! not heat; timestep refinement is still necessary for time accuracy even when
//! the work identity closes.
//!
//! Density is an explicit original mass inventory. `FinitePorousMaterial`
//! currently exposes no sourced density; borrowing a source-resolved law keeps
//! its card and context but makes no density-calibration claim. Conditions stay
//! fixed: this path adds no evolving porosity, densification, thermal coupling,
//! contact, mass transfer or pressure stabilization. Held current-area pressure
//! uses the native surface force, analytic follower tangent and exact declared
//! path work in the same momentum solve; it is never frozen to reference area.
//! Sealed cavities use that same solve with original gas storage and acoustic
//! derivatives. Their homogeneous gas closure resolves neither gas momentum
//! nor pressure waves, and transfers no plastic heat into the gas.

use std::sync::Arc;

use fs_exec::Cx;
use fs_material::plastic::finite::porous::{FinitePorous, FinitePorousMaterial};
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;

use super::cavity::FinitePorousCavityEnergy;
use super::{
    FinitePorousPressureEnergy, FinitePorousTetError, FinitePorousTetSolid, FinitePorousTetStep,
};
use crate::hyper3::plastic::{
    AcceptedMotion, FinitePlasticTetEnergy, checkpoint, energy_gate, invalid,
};
use crate::hyper3::pressure::cavity::{AdiabaticCavity, AdiabaticCavityGas, TetCavityGasLoad};
use crate::hyper3::pressure::{ClosedTetPressureSurface, TetPressure};
use crate::hyper3::{HyperTetSettings, HyperTetSolution, inertia_error};
use crate::inertia::{BackwardEulerInertia, P1ReferenceMass};

/// Original separate material and numerical-inertia admission controls.
pub use crate::hyper3::plastic::dynamics::FinitePlasticMotionControls;
/// Original physical kinetic-work account, separate from material heat.
pub use crate::inertia::InertialEnergy;

/// Accepted porous configuration, complete history, velocity and original mass.
#[derive(Debug)]
pub struct FinitePorousTetMotion<'a, L: FinitePorousMaterial = FinitePorous> {
    solid: FinitePorousTetSolid<'a, L>,
    velocity_mps: Vec<[f64; 3]>,
    reference_density_kg_m3: Arc<[f64]>,
    reference_nodal_density_kg_m3: Option<Arc<[f64]>>,
    mass: Arc<P1ReferenceMass>,
}

/// Complete admitted physical interval; the previously accepted body is intact.
#[derive(Debug)]
pub struct FinitePorousMotionStep<'a, L: FinitePorousMaterial = FinitePorous> {
    /// Complete next state, including non-isochoric plastic metric and mass.
    pub motion: FinitePorousTetMotion<'a, L>,
    /// Native solve evidence and support reactions including consistent inertia.
    /// `strain_energy_j` remains the material's reduced incremental potential.
    pub equilibrium: HyperTetSolution,
    /// Physical storage/yield loss and signed material integration defect.
    pub material_energy: FinitePlasticTetEnergy,
    /// Kinetic change, momentum work and separate numerical kinetic loss.
    pub inertia: InertialEnergy,
    /// Present for held current-area pressure, with the original pressure
    /// owner's exact work and signed endpoint correction. The combined external
    /// defect accounts only for material/pressure endpoint integration,
    /// separately from inertial work and numerical kinetic loss.
    pub pressure: Option<FinitePorousPressureEnergy>,
    /// External work minus storage change, plastic dissipation, kinetic change,
    /// numerical kinetic loss and the signed external endpoint-work defect.
    /// Computed from the native momentum and independent kinetic identities to
    /// retain their precision instead of subtracting large absolute energies.
    /// A cavity-coupled result supplies the complete external defect in its
    /// enclosing `FinitePorousCavitiesMotionStep`.
    pub balance_error_j: f64,
}

/// Complete admitted motion and independently sealed gas endpoints.
///
/// All gas receipts retain the original source order. The whole-system work
/// identity is dead + support + held-pressure + summed exterior work = material
/// storage change + summed gas storage changes + physical plastic loss + kinetic
/// change + numerical kinetic loss + `external_work_defect_j` +
/// `step.balance_error_j`. Gas compression work is storage, not material heat.
///
/// Continue through `step.motion.advance_with_cavities` using the same original
/// undeformed-volume gas references. Deformed gas endpoints are output evidence,
/// never replacement reference parcels.
#[derive(Debug)]
pub struct FinitePorousCavitiesMotionStep<
    'a,
    L: FinitePorousMaterial = FinitePorous,
    G = ClosedGasParcel,
> {
    /// Complete material/momentum state; equilibrium includes all cavity loads.
    pub step: FinitePorousMotionStep<'a, L>,
    /// Native accepted gas states with independent inventories and tangents.
    pub cavities: Vec<TetCavityGasLoad<G>>,
    /// Stable original-source gas, exterior, and endpoint work for each cavity.
    pub cavity_energy: Vec<FinitePorousCavityEnergy>,
    /// Material endpoint defect plus held-pressure and individual cavity
    /// corrections. Kinetic numerical loss retains its own separate account.
    pub external_work_defect_j: f64,
}

fn validate_nodal_vectors(
    node_count: usize,
    values: &[[f64; 3]],
    invalid_input: &'static str,
    cx: &Cx<'_>,
) -> Result<(), FinitePorousTetError> {
    checkpoint(cx)?;
    if values.len() != node_count {
        return Err(invalid(invalid_input).into());
    }
    for (node, vector) in values.iter().enumerate() {
        if node % 256 == 0 {
            checkpoint(cx)?;
        }
        if vector.iter().any(|value| !value.is_finite()) {
            return Err(invalid(invalid_input).into());
        }
    }
    Ok(())
}

impl<'a, L: FinitePorousMaterial> FinitePorousTetMotion<'a, L> {
    /// Bind a virgin or previously compacted solid to fixed mass and velocity.
    ///
    /// One finite positive density per original tetrahedron defines reference
    /// mass, regardless of accepted total/plastic volume. The solid retains its
    /// original material owners and all accepted histories. The original
    /// consistent mass matrix is prepared once and shared by later endpoints.
    ///
    /// # Errors
    /// Invalid density/velocity, unrepresentable reference mass or cancellation.
    pub fn try_new(
        solid: FinitePorousTetSolid<'a, L>,
        reference_density_kg_m3: &[f64],
        velocity_mps: &[[f64; 3]],
        cx: &Cx<'_>,
    ) -> Result<Self, FinitePorousTetError> {
        checkpoint(cx)?;
        if velocity_mps.len() != solid.state.displacement_m.len()
            || velocity_mps
                .iter()
                .flatten()
                .any(|value| !value.is_finite())
        {
            return Err(
                invalid("one finite initial velocity vector per original node required").into(),
            );
        }
        let mass = P1ReferenceMass::new(
            solid
                .state
                .elements
                .iter()
                .map(|element| (element.nodes, element.volume)),
            reference_density_kg_m3,
            solid.state.displacement_m.as_flattened().len(),
            cx,
        )
        .map_err(inertia_error)?;
        checkpoint(cx)?;
        Ok(Self {
            solid,
            velocity_mps: velocity_mps.to_vec(),
            reference_density_kg_m3: Arc::from(reference_density_kg_m3),
            reference_nodal_density_kg_m3: None,
            mass: Arc::new(mass),
        })
    }

    /// Bind positive nodal density through exact P1 reference-mass integration.
    ///
    /// The native mass owner integrates `rho_ref * N_a * N_b`, retaining the
    /// original field and its full matrix across every compaction step. Cell
    /// means are available for volume/density queries but do not replace this
    /// mass distribution or its coupling to moving supports.
    ///
    /// # Errors
    /// Invalid nodal density/velocity, unrepresentable mass, or cancellation.
    pub fn try_new_with_nodal_density(
        solid: FinitePorousTetSolid<'a, L>,
        reference_nodal_density_kg_m3: &[f64],
        velocity_mps: &[[f64; 3]],
        cx: &Cx<'_>,
    ) -> Result<Self, FinitePorousTetError> {
        validate_nodal_vectors(
            solid.state.nodes_m.len(),
            velocity_mps,
            "one finite initial velocity vector per original node required",
            cx,
        )?;
        let (mass, means) = P1ReferenceMass::from_nodal_density(
            solid
                .state
                .elements
                .iter()
                .map(|element| (element.nodes, element.volume)),
            reference_nodal_density_kg_m3,
            solid.state.displacement_m.as_flattened().len(),
            cx,
        )
        .map_err(inertia_error)?;
        checkpoint(cx)?;
        Ok(Self {
            solid,
            velocity_mps: velocity_mps.to_vec(),
            reference_density_kg_m3: Arc::from(means),
            reference_nodal_density_kg_m3: Some(Arc::from(reference_nodal_density_kg_m3)),
            mass: Arc::new(mass),
        })
    }

    /// Accepted finite configuration, source binding and full compaction history.
    pub fn solid(&self) -> &FinitePorousTetSolid<'a, L> {
        &self.solid
    }

    /// Accepted velocity in original node order, m/s.
    pub fn velocity_mps(&self) -> &[[f64; 3]] {
        &self.velocity_mps
    }

    /// Original cell density, or exact cell means of a retained P1 nodal field.
    pub fn reference_density_kg_m3(&self) -> &[f64] {
        &self.reference_density_kg_m3
    }

    /// Complete original nodal density field, if that constructor was selected.
    pub fn reference_nodal_density_kg_m3(&self) -> Option<&[f64]> {
        self.reference_nodal_density_kg_m3.as_deref()
    }

    /// Original mass, unchanged by deformation or permanent plastic compaction.
    pub fn reference_mass_kg(&self) -> f64 {
        self.mass.reference_mass_kg()
    }

    /// Density over each current cell: retained reference density divided by J.
    /// The plastic Jacobian never rescales the mass inventory.
    ///
    /// # Errors
    /// Unrepresentable current density or cancellation.
    pub fn current_density_kg_m3(&self, cx: &Cx<'_>) -> Result<Vec<f64>, FinitePorousTetError> {
        self.solid
            .current_density_kg_m3(&self.reference_density_kg_m3, cx)
    }

    /// Project a nodal acceleration field through this body's consistent mass.
    ///
    /// Returned original DOF/force pairs are `M*a` in newtons, including uniform
    /// gravity on graded material. Feed them to `advance` as held nodal loads;
    /// do not also apply the same acceleration as a body-force density. Evolving
    /// external-field potentials and tangents remain with their owning model.
    pub fn reference_acceleration_load_n(
        &self,
        acceleration_m_s2: &[[f64; 3]],
        cx: &Cx<'_>,
    ) -> Result<Vec<(usize, f64)>, FinitePorousTetError> {
        validate_nodal_vectors(
            self.velocity_mps.len(),
            acceleration_m_s2,
            "one finite acceleration vector per original node required",
            cx,
        )?;
        let forces = self
            .mass
            .action(acceleration_m_s2.as_flattened(), cx)
            .map_err(inertia_error)?;
        let mut loads = Vec::with_capacity(forces.len());
        for (dof, force) in forces.into_iter().enumerate() {
            if dof % 256 == 0 {
                checkpoint(cx)?;
            }
            if force != 0.0 {
                loads.push((dof, force));
            }
        }
        checkpoint(cx)?;
        Ok(loads)
    }

    /// Advance one physical interval under held dead loads and absolute supports.
    ///
    /// Body force is force per original reference volume. Free bodies require
    /// no artificial restraints. Prescribed motion participates in the full
    /// consistent mass action and support reactions. Every trial borrows the
    /// same accepted history; any failed solve or gate leaves `self` unchanged.
    ///
    /// # Errors
    /// Invalid timestep/controls, native solve/material refusal, work or kinetic
    /// loss budget exceedance, unrepresentable motion or cancellation.
    #[allow(clippy::too_many_arguments)]
    pub fn advance(
        &self,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        dt_s: f64,
        settings: HyperTetSettings,
        controls: FinitePlasticMotionControls,
        cx: &Cx<'_>,
    ) -> Result<FinitePorousMotionStep<'a, L>, FinitePorousTetError> {
        self.advance_loaded(
            prescribed_m,
            nodal_forces_n,
            body_force_n_m3,
            dt_s,
            settings,
            controls,
            None,
            cx,
        )
    }

    /// Advance momentum and retained compaction under held current-area pressure.
    ///
    /// The exact original closed surface, all supported pressure fields, current
    /// normals and areas, and the complete pressure derivative enter the native
    /// backward-Euler solve. Positive differential pressure compresses the body.
    /// The field is held over this interval; time-varying pressure requires
    /// resolved successive intervals. This prescribes a reservoir load and does
    /// not add a fluid equation of state or evolve a sealed cavity.
    ///
    /// The pressure receipt retains exact declared-path work separately from
    /// its endpoint-force work. The material work budget bounds the SUM of the
    /// absolute material and pressure endpoint defects; the kinetic-loss budget
    /// remains independent. Neither numerical defect is physical plastic heat.
    /// Free bodies need no artificial supports, and prescribed-motion reactions
    /// include pressure and the complete original consistent-mass action.
    ///
    /// # Errors
    /// The same transactional gates as `advance`, plus exact surface/pressure
    /// admission. Any refused trial, final work gate or cancellation leaves
    /// displacement, velocity, history and original mass unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_pressure<'p>(
        &self,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        dt_s: f64,
        surface: &ClosedTetPressureSurface,
        pressure: impl Into<TetPressure<'p>>,
        settings: HyperTetSettings,
        controls: FinitePlasticMotionControls,
        cx: &Cx<'_>,
    ) -> Result<FinitePorousMotionStep<'a, L>, FinitePorousTetError> {
        self.advance_loaded(
            prescribed_m,
            nodal_forces_n,
            body_force_n_m3,
            dt_s,
            settings,
            controls,
            Some((surface, pressure.into())),
            cx,
        )
    }

    /// Advance solid momentum and retained compaction with independently sealed gases.
    ///
    /// Every cavity retains its original reference geometry and undeformed gas
    /// parcel. Current pressure, follower geometry, original acoustic derivatives,
    /// porous tangents and the full consistent M/dt² act in one native solve.
    /// The original ideal or Peng--Robinson gas owns its caloric and phase-domain
    /// admission. This homogeneous isentropic closure neglects gas momentum and
    /// pressure waves; physical plastic loss is not transferred to gas heat.
    ///
    /// Direct gas energy increments and compact geometric work are evaluated
    /// between accepted endpoints, without subtracting reference-relative
    /// potentials or rebuilding reference parcels. Every cavity keeps its own
    /// mass, species, entropy and exterior pressure. Additional held pressure,
    /// dead/body loads and absolute moving supports retain their existing roles.
    ///
    /// One physical interval is accepted only after material, cavity, momentum
    /// and work gates pass. The sum of absolute material, held-pressure and
    /// individual cavity endpoint defects uses the material budget; numerical
    /// kinetic loss uses its independent budget. An empty cavity slice follows
    /// the original pressure/dead-load dynamic path.
    ///
    /// # Errors
    /// Original gas-domain/source/geometry refusal, native solve refusal,
    /// unrepresentable timestep or work, an exceeded material/kinetic budget,
    /// or cancellation. Every refusal preserves the original motion and every
    /// borrowed gas reference, including after a converged coupled solve.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_cavities<G: AdiabaticCavityGas>(
        &self,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        dt_s: f64,
        cavities: &[AdiabaticCavity<'_, G>],
        pressure: Option<(&ClosedTetPressureSurface, TetPressure<'_>)>,
        settings: HyperTetSettings,
        controls: FinitePlasticMotionControls,
        cx: &Cx<'_>,
    ) -> Result<FinitePorousCavitiesMotionStep<'a, L, G>, FinitePorousTetError> {
        let inertia = self.prepare_inertia(dt_s, controls, cx)?;
        let (coupled, motion) = self.solid.advance_cavities_loaded(
            prescribed_m,
            nodal_forces_n,
            body_force_n_m3,
            cavities,
            pressure,
            settings,
            controls.material,
            Some(&inertia),
            cx,
        )?;
        let step = self.accept_motion(
            coupled.step,
            motion.expect("inertial cavity endpoint returns native motion"),
            controls,
            cx,
        )?;
        checkpoint(cx)?;
        Ok(FinitePorousCavitiesMotionStep {
            step,
            cavities: coupled.cavities,
            cavity_energy: coupled.cavity_energy,
            external_work_defect_j: coupled.external_work_defect_j,
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn advance_loaded(
        &self,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        dt_s: f64,
        settings: HyperTetSettings,
        controls: FinitePlasticMotionControls,
        pressure: Option<(&ClosedTetPressureSurface, TetPressure<'_>)>,
        cx: &Cx<'_>,
    ) -> Result<FinitePorousMotionStep<'a, L>, FinitePorousTetError> {
        let inertia = self.prepare_inertia(dt_s, controls, cx)?;
        let (step, motion) = self.solid.advance_loaded(
            prescribed_m,
            nodal_forces_n,
            body_force_n_m3,
            settings,
            controls.material,
            Some(&inertia),
            pressure,
            cx,
        )?;
        self.accept_motion(
            step,
            motion.expect("inertial endpoint returns native motion"),
            controls,
            cx,
        )
    }

    fn prepare_inertia<'s>(
        &'s self,
        dt_s: f64,
        controls: FinitePlasticMotionControls,
        cx: &Cx<'_>,
    ) -> Result<BackwardEulerInertia<'s>, FinitePorousTetError> {
        checkpoint(cx)?;
        if !controls.maximum_numerical_kinetic_loss_j.is_finite()
            || controls.maximum_numerical_kinetic_loss_j < 0.0
        {
            return Err(
                invalid("finite nonnegative numerical kinetic-loss budget required").into(),
            );
        }
        Ok(BackwardEulerInertia::new(
            &self.mass,
            self.solid.state.displacement_m.as_flattened(),
            self.velocity_mps.as_flattened(),
            dt_s,
            cx,
        )
        .map_err(inertia_error)?)
    }

    fn accept_motion(
        &self,
        step: FinitePorousTetStep<'a, L>,
        motion: AcceptedMotion,
        controls: FinitePlasticMotionControls,
        cx: &Cx<'_>,
    ) -> Result<FinitePorousMotionStep<'a, L>, FinitePorousTetError> {
        let (velocity_mps, inertia) = motion;
        energy_gate(
            "porous backward-Euler numerical kinetic loss",
            inertia.numerical_loss_j,
            controls.maximum_numerical_kinetic_loss_j,
        )?;
        let balance_error_j = step.energy.equilibrium_work_residual_j + inertia.balance_error_j;
        energy_gate(
            "total finite porous inertial work",
            balance_error_j,
            controls.material.equilibrium_work_tolerance_j,
        )?;
        checkpoint(cx)?;
        Ok(FinitePorousMotionStep {
            motion: Self {
                solid: step.solid,
                velocity_mps,
                reference_density_kg_m3: Arc::clone(&self.reference_density_kg_m3),
                reference_nodal_density_kg_m3: self
                    .reference_nodal_density_kg_m3
                    .as_ref()
                    .map(Arc::clone),
                mass: Arc::clone(&self.mass),
            },
            equilibrium: step.equilibrium,
            material_energy: step.energy,
            inertia,
            pressure: step.pressure,
            balance_error_j,
        })
    }
}
