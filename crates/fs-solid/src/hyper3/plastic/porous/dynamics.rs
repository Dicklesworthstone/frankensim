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
//! defect and backward-Euler numerical kinetic loss have separate accounts and
//! budgets. Numerical kinetic loss is not heat; timestep refinement is still
//! necessary for time accuracy even when the work identity closes.
//!
//! Density is an explicit original mass inventory. `FinitePorousMaterial`
//! currently exposes no sourced density; borrowing a source-resolved law keeps
//! its card and context but makes no density-calibration claim. Conditions stay
//! fixed: this path adds no evolving porosity, densification, thermal coupling,
//! contact, follower loads, mass transfer or pressure stabilization.

use std::sync::Arc;

use fs_exec::Cx;
use fs_material::plastic::finite::porous::{FinitePorous, FinitePorousMaterial};

use super::{FinitePorousTetError, FinitePorousTetSolid};
use crate::hyper3::plastic::{FinitePlasticTetEnergy, checkpoint, energy_gate, invalid};
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
    /// Combined material/momentum and independent kinetic-identity residual.
    pub balance_error_j: f64,
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

    /// Immutable original cellwise density, kg/m³, in reference element order.
    pub fn reference_density_kg_m3(&self) -> &[f64] {
        &self.reference_density_kg_m3
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
    ) -> Result<Fini…2037 tokens truncated…       let dissipated = if index == 0 { Y * alpha / 6.0 } else { 0.0 };
            near(step.solid.stored_energy_j(), stored, 4e-7);
            near(step.energy.stored_after_j, stored, 4e-7);
            near(step.energy.plastic_dissipation_j, dissipated, 3e-7);
            near(
                step.energy.incremental_potential_j,
                stored + dissipated,
                4e-7,
            );
            near(step.equilibrium.strain_energy_j, stored + dissipated, 4e-7);
            near(step.energy.equilibrium_work_residual_j, 0.0, 2e-7);
            near(
                step.energy.dead_load_endpoint_work_j + step.energy.support_endpoint_work_j,
                stored - old_storage + dissipated + step.energy.endpoint_work_defect_j,
                8e-7,
            );
            if index > 0 {
                assert_eq!(
                    step.solid.element_states(),
                    old_state,
                    "elastic unload/reload must retain compaction"
                );
                near(step.energy.plastic_dissipation_j, 0.0, 0.0);
            }
            heat += step.energy.plastic_dissipation_j;
            solid = step.solid;
        }
        assert!(plastic_j < 0.55);
        assert!(unloaded < 0.82);
        near(heat, Y * alpha / 6.0, 4e-7);
        assert_eq!(solid.reference_nodes_m(), NODES);
        assert_eq!(solid.reference_tetrahedra(), TETS);
    });
}

#[test]
fn heterogeneous_elements_and_rigid_motion_preserve_distinct_plastic_volumes() {
    with_cx(|cx, _| {
        // A shared face joins reference volumes 1/6 and 1/3. The second
        // material is pressure independent and cannot compact plastically.
        let nodes = [NODES[0], NODES[1], NODES[2], NODES[3], [1.0; 3]];
        let tets = [[0, 1, 2, 3], [4, 2, 1, 3]];
        let materials = [law(A), law(0.0)];
        let solid = FinitePorousTetSolid::try_new(
            &nodes,
            &tets,
            FinitePlasticLawField::PerElement(&materials),
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        let stretch = 0.85;
        let f = [stretch, 0.0, 0.0, 0.0, stretch, 0.0, 0.0, 0.0, stretch];
        let mut bounded = controls();
        // A finite rotation of a prestressed body has a signed endpoint
        // quadrature error. Permit it explicitly, while checking heat is zero.
        bounded.maximum_endpoint_work_defect_j = 100.0;
        let compacted = solid
            .advance(
                &affine(&nodes, f, [0.0; 3]),
                &[],
                [0.0; 3],
                HyperTetSettings::default(),
                bounded,
                cx,
            )
            .unwrap();
        let (alpha, jp, tau) = hydrostatic(stretch);
        near(compacted.solid.element_states()[0].alpha, alpha, 3e-8);
        near(
            materials[0]
                .plastic_jacobian(&compacted.solid.element_states()[0])
                .unwrap(),
            jp,
            3e-8,
        );
        assert_eq!(
            compacted.solid.element_states()[1],
            materials[1].initial_state()
        );
        near(
            compacted.energy.plastic_dissipation_j,
            Y * alpha / 6.0,
            4e-7,
        );
        let elastic_tau = 3.0 * K * stretch.ln();
        let storage = (0.5 * tau * tau / K + 0.5 * H * alpha * alpha) / 6.0
            + 0.5 * elastic_tau * elastic_tau / K / 3.0;
        near(compacted.solid.stored_energy_j(), storage, 5e-7);
        let (sin, cos) = 0.7_f64.sin_cos();
        let rotation = [cos, -sin, 0.0, sin, cos, 0.0, 0.0, 0.0, 1.0];
        let rotated_f = rotation.map(|value| stretch * value);
        let rotated = compacted
            .solid
            .advance(
                &affine(&nodes, rotated_f, [0.4, -0.7, 1.2]),
                &[],
                [0.0; 3],
                HyperTetSettings::default(),
                bounded,
                cx,
            )
            .unwrap();
        assert_eq!(
            rotated.solid.element_states(),
            compacted.solid.element_states()
        );
        near(rotated.energy.plastic_dissipation_j, 0.0, 0.0);
        near(rotated.solid.stored_energy_j(), storage, 5e-7);
        assert!(rotated.energy.endpoint_work_defect_j < -1.0);
        for element in 0..2 {
            let old_p = compacted.solid.element_first_piola_pa()[element];
            let new_p = rotated.solid.element_first_piola_pa()[element];
            for i in 0..9 {
                let expected = (0..3)
                    .map(|j| rotation[3 * (i / 3) + j] * old_p[3 * j + i % 3])
                    .sum();
                near(new_p[i], expected, 5e-7);
                near(
                    rotated.solid.element_deformation()[element][i],
                    rotated_f[i],
                    3e-8,
                );
            }
        }
    });
}

#[test]
fn refused_work_and_cancellation_leave_compacted_state_retryable() {
    with_cx(|cx, gate| {
        let material = law(A);
        let initial = FinitePorousTetSolid::try_new(
            &NODES,
            &TETS,
            FinitePlasticLawField::Uniform(&material),
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        let f = [0.9, 0.0, 0.0, 0.0, 0.9, 0.0, 0.0, 0.0, 0.9];
        let solid = initial
            .advance(
                &affine(&NODES, f, [0.0; 3]),
                &[],
                [0.0; 3],
                HyperTetSettings::default(),
                controls(),
                cx,
            )
            .unwrap()
            .solid;
        let before_state = solid.element_states().to_vec();
        let before_f = solid.element_deformation().to_vec();
        let before_u = solid.displacement_m().to_vec();
        let before_p = solid.element_first_piola_pa().to_vec();
        let before_storage = solid.stored_energy_j();
        let target = affine(
            &NODES,
            [0.8, 0.0, 0.0, 0.0, 0.8, 0.0, 0.0, 0.0, 0.8],
            [0.0; 3],
        );
        let mut strict = controls();
        strict.maximum_endpoint_work_defect_j = 0.0;
        let refused = solid
            .advance(
                &target,
                &[],
                [0.0; 3],
                HyperTetSettings::default(),
                strict,
                cx,
            )
            .unwrap_err();
        assert!(matches!(
            refused,
            FinitePorousTetError::Solid(FinitePlasticTetError::EnergyGate {
                quantity: "porous signed endpoint integration defect",
                ..
            })
        ));
        let retry = solid
            .advance(
                &target,
                &[],
                [0.0; 3],
                HyperTetSettings::default(),
                controls(),
                cx,
            )
            .unwrap();
        let repeat = solid
            .advance(
                &target,
                &[],
                [0.0; 3],
                HyperTetSettings::default(),
                controls(),
                cx,
            )
            .unwrap();
        assert_eq!(retry.equilibrium, repeat.equilibrium);
        assert_eq!(retry.energy, repeat.energy);
        assert_eq!(retry.solid.element_states(), repeat.solid.element_states());
        assert!(retry.solid.element_states()[0].alpha > before_state[0].alpha);
        gate.request();
        let cancelled = solid
            .advance(
                &target,
                &[],
                [0.0; 3],
                HyperTetSettings::default(),
                controls(),
                cx,
            )
            .unwrap_err();
        assert!(matches!(
            cancelled,
            FinitePorousTetError::Solid(FinitePlasticTetError::Native(HyperTetError::Cancelled))
        ));
        assert_eq!(solid.element_states(), before_state);
        assert_eq!(solid.element_deformation(), before_f);
        assert_eq!(solid.displacement_m(), before_u);
        assert_eq!(solid.element_first_piola_pa(), before_p);
        near(solid.stored_energy_j(), before_storage, 0.0);
    });
}

#[test]
fn law_count_and_physical_endpoint_controls_are_admitted_before_state_changes() {
    with_cx(|cx, _| {
        let material = law(A);
        let no_laws: [FinitePorous; 0] = [];
        let refused = FinitePorousTetSolid::try_new(
            &NODES,
            &TETS,
            FinitePlasticLawField::PerElement(&no_laws),
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap_err();
        assert!(matches!(
            refused,
            FinitePorousTetError::Solid(FinitePlasticTetError::Native(
                HyperTetError::LawCountMismatch {
                    expected: 1,
                    got: 0
                }
            ))
        ));
        let solid = FinitePorousTetSolid::try_new(
            &NODES,
            &TETS,
            FinitePlasticLawField::Uniform(&material),
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        let multiple = HyperTetSettings {
            load_steps: 2,
            ..HyperTetSettings::default()
        };
        assert!(
            solid
                .advance(&SUPPORTS, &[], [0.0; 3], multiple, controls(), cx)
                .is_err()
        );
        let mut nonfinite = controls();
        nonfinite.yield_tolerance_pa = f64::NAN;
        assert!(
            solid
                .advance(
                    &SUPPORTS,
                    &[],
                    [0.0; 3],
                    HyperTetSettings::default(),
                    nonfinite,
                    cx
                )
                .is_err()
        );
        assert_eq!(solid.element_states(), [material.initial_state()]);
        assert_eq!(solid.displacement_m(), [[0.0; 3]; 4]);
    });
}
