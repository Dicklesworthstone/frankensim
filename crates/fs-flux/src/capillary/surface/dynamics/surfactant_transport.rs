//! Moving insoluble material surfaces with conservative diffusion and mechanics.
//!
//! A half diffusion step on the initial areas, the original mechanical step,
//! and a half diffusion step on the accepted areas form a symmetric composition.
//! Both stages consume the SAME amounts and free-energy law. Native trajectory
//! and amount arrays publish together only after every stage and the combined
//! mass/energy account pass. The existing trajectory history is never cloned.
//!
//! The graph and its physical conductances are held through the interval. A
//! changing spatial metric or diffusivity must be supplied by its actual owner;
//! this module does not infer transport distances from triangle centroids.
//! Temporal and spatial refinement remain necessary. Diffusion loss is an
//! isothermal free-energy loss, not an automatic heat or internal-energy source.
//! Surfactant inertia is neglected: material reference nodal masses stay fixed
//! while the coating carries interfacial free energy and changes its stress.

use super::super::surfactant::{
    SurfactantFace, SurfactantLaw,
    diffusion::{
        SurfaceDiffusionControls, SurfaceDiffusionEdge, SurfaceDiffusionReport,
        SurfacePartitionDiffusionEdge,
        adaptive::{
            AdaptiveSurfaceDiffusionControls, AdaptiveSurfaceDiffusionError,
            AdaptiveSurfaceDiffusionReport, diffuse_partitioned_surface_surfactant_adaptive,
            diffuse_surface_surfactant_adaptive,
        },
        diffuse_partitioned_surface_surfactant, diffuse_surface_surfactant,
    },
};
use super::*;

/// One homogeneous isothermal coating and its held physical diffusion graph.
#[derive(Clone, Copy, Debug)]
pub struct SurfaceSurfactantTransport<'a> {
    /// One held-temperature, common chemical free-energy law for all faces.
    pub law: SurfactantLaw,
    /// Symmetric authored conductances in fixed canonical edge order.
    pub edges: &'a [SurfaceDiffusionEdge],
    /// Substep, work, allocation and physical-balance limits for each half step.
    pub controls: SurfaceDiffusionControls,
}

/// Held chemical laws and a partition-aware graph on the original material faces.
#[derive(Clone, Copy, Debug)]
pub struct SurfacePartitionedSurfactantTransport<'a> {
    /// One original constitutive law per face. Connected patches share the
    /// same family and thermal energy; references/capacities and clean tensions
    /// may differ. Ideal references use one common chemical-energy datum.
    pub laws: &'a [SurfactantLaw],
    /// Symmetric authored molar conductances in fixed canonical edge order.
    pub edges: &'a [SurfacePartitionDiffusionEdge],
    /// Substep, work, allocation and physical-balance limits for each half step.
    pub controls: SurfaceDiffusionControls,
}

/// A homogeneous coating with controlled internal diffusion splitting error.
#[derive(Clone, Copy, Debug)]
pub struct AdaptiveSurfaceSurfactantTransport<'a> {
    pub law: SurfactantLaw,
    pub edges: &'a [SurfaceDiffusionEdge],
    /// Trial and pair budgets cover BOTH diffusion halves, including rejected
    /// and coarse evaluations. Concentration tolerances control each accepted
    /// internal interval, not the mechanical/composition time error.
    pub controls: AdaptiveSurfaceDiffusionControls,
}

/// Original per-face laws with adaptively resolved partition diffusion.
#[derive(Clone, Copy, Debug)]
pub struct AdaptiveSurfacePartitionedSurfactantTransport<'a> {
    pub laws: &'a [SurfactantLaw],
    pub edges: &'a [SurfacePartitionDiffusionEdge],
    /// One combined work budget and physical allowance for both halves.
    /// Mechanical and diffusion/motion composition errors still require
    /// independent refinement of the caller's macro time step.
    pub controls: AdaptiveSurfaceDiffusionControls,
}

/// Which held-geometry diffusion interval refused the atomic composition.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SurfaceDiffusionPhase {
    BeforeMotion,
    AfterMotion,
}

#[derive(Clone, Copy)]
enum Transport<'a> {
    Uniform(SurfaceSurfactantTransport<'a>),
    Partitioned(SurfacePartitionedSurfactantTransport<'a>),
    AdaptiveUniform(AdaptiveSurfaceSurfactantTransport<'a>),
    AdaptivePartitioned(AdaptiveSurfacePartitionedSurfactantTransport<'a>),
}

struct DiffusionStageReport {
    physical: SurfaceDiffusionReport,
    adaptive: Option<AdaptiveSurfaceDiffusionReport>,
}

#[derive(Default)]
struct AdaptiveWork {
    trials: usize,
    pair_exchanges: usize,
}

impl AdaptiveWork {
    fn stage(
        &mut self,
        controls: AdaptiveSurfaceDiffusionControls,
        phase: SurfaceDiffusionPhase,
        advance: impl FnOnce(
            AdaptiveSurfaceDiffusionControls,
        )
            -> Result<AdaptiveSurfaceDiffusionReport, AdaptiveSurfaceDiffusionError>,
    ) -> Result<DiffusionStageReport, SurfaceDynamicsError> {
        let refusal = |error| SurfaceDynamicsError::AdaptiveDiffusion { phase, error };
        // Preserve the source's initial control validation when the original
        // trial cap is zero. Exhaustion after a successful half is a budget
        // refusal, not an invalid-control error from passing zero downstream.
        if self.trials > 0 && self.trials == controls.maximum_trials {
            return Err(refusal(AdaptiveSurfaceDiffusionError::TrialBudget {
                completed_trials: self.trials,
            }));
        }
        let mut remaining = controls;
        remaining.maximum_trials = controls.maximum_trials.checked_sub(self.trials).ok_or(
            SurfaceDynamicsError::Invalid("adaptive trial accounting exceeded its budget"),
        )?;
        remaining.maximum_pair_exchanges = controls
            .maximum_pair_exchanges
            .checked_sub(self.pair_exchanges)
            .ok_or(SurfaceDynamicsError::Invalid(
                "adaptive pair accounting exceeded its budget",
            ))?;
        let report = match advance(remaining) {
            Ok(report) => report,
            Err(AdaptiveSurfaceDiffusionError::TrialBudget { completed_trials }) => {
                return Err(refusal(AdaptiveSurfaceDiffusionError::TrialBudget {
                    completed_trials: self.trials.checked_add(completed_trials).ok_or(
                        SurfaceDynamicsError::Invalid("adaptive trial count overflow"),
                    )?,
                }));
            }
            Err(AdaptiveSurfaceDiffusionError::PairExchangeBudget {
                completed,
                required,
            }) => {
                return Err(refusal(AdaptiveSurfaceDiffusionError::PairExchangeBudget {
                    completed: self.pair_exchanges.checked_add(completed).ok_or(
                        SurfaceDynamicsError::Invalid("adaptive pair count overflow"),
                    )?,
                    required,
                }));
            }
            Err(error) => return Err(refusal(error)),
        };
        self.trials =
            self.trials
                .checked_add(report.trials)
                .ok_or(SurfaceDynamicsError::Invalid(
                    "adaptive trial count overflow",
                ))?;
        self.pair_exchanges = self
            .pair_exchanges
            .checked_add(report.evaluated_pair_exchanges)
            .ok_or(SurfaceDynamicsError::Invalid(
                "adaptive pair count overflow",
            ))?;
        Ok(DiffusionStageReport {
            physical: report.diffusion,
            adaptive: Some(report),
        })
    }
}

impl Transport<'_> {
    fn diffuse(
        self,
        areas: &[f64],
        amounts: &mut [f64],
        duration_s: f64,
        phase: SurfaceDiffusionPhase,
        work: &mut AdaptiveWork,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DiffusionStageReport, SurfaceDynamicsError> {
        let physical = match self {
            Self::Uniform(transport) => diffuse_surface_surfactant(
                transport.law,
                areas,
                amounts,
                transport.edges,
                duration_s,
                transport.controls,
                cancelled,
            )?,
            Self::Partitioned(transport) => diffuse_partitioned_surface_surfactant(
                transport.laws,
                areas,
                amounts,
                transport.edges,
                duration_s,
                transport.controls,
                cancelled,
            )?,
            Self::AdaptiveUniform(transport) => {
                return work.stage(transport.controls, phase, |remaining| {
                    diffuse_surface_surfactant_adaptive(
                        transport.law,
                        areas,
                        amounts,
                        transport.edges,
                        duration_s,
                        remaining,
                        cancelled,
                    )
                });
            }
            Self::AdaptivePartitioned(transport) => {
                return work.stage(transport.controls, phase, |remaining| {
                    diffuse_partitioned_surface_surfactant_adaptive(
                        transport.laws,
                        areas,
                        amounts,
                        transport.edges,
                        duration_s,
                        remaining,
                        cancelled,
                    )
                });
            }
        };
        Ok(DiffusionStageReport {
            physical,
            adaptive: None,
        })
    }
}

/// Accepted composition, including independent substage and complete accounts.
#[derive(Clone, Debug)]
pub struct SurfaceSurfactantReport {
    /// Original midpoint mechanical response at the intermediate inventory.
    pub mechanics: SurfaceDynamicsReport,
    /// Initial-geometry half diffusion step.
    pub diffusion_before: SurfaceDiffusionReport,
    /// Solved-final-geometry half diffusion step.
    pub diffusion_after: SurfaceDiffusionReport,
    /// Surface storage before any diffusion or motion [J].
    pub surface_free_energy_before_j: f64,
    /// Surface storage after motion and both diffusion stages [J].
    pub surface_free_energy_after_j: f64,
    /// Sum of constitutive increments, avoiding subtraction of large totals [J].
    pub surface_free_energy_change_j: f64,
    /// Isothermal free-energy loss from both diffusion stages [J].
    pub diffusion_free_energy_dissipation_j: f64,
    /// Compensated amount change over the entire accepted step [mol].
    pub amount_defect_mol: f64,
    /// Delta K + delta E_surface + viscous loss + diffusion loss - external work [J].
    pub energy_residual_j: f64,
}

/// One physical coupled account and the actual adaptive work in both halves.
#[derive(Clone, Debug)]
pub struct AdaptiveSurfaceSurfactantReport {
    pub coupled: SurfaceSurfactantReport,
    pub diffusion_before: AdaptiveSurfaceDiffusionReport,
    pub diffusion_after: AdaptiveSurfaceDiffusionReport,
    /// All comparisons in both halves, including rejected comparisons.
    pub total_trials: usize,
    /// Coarse, accepted and rejected pair applications across both halves.
    pub total_evaluated_pair_exchanges: usize,
}

fn adaptive_report(
    coupled: SurfaceSurfactantReport,
    stages: [Option<AdaptiveSurfaceDiffusionReport>; 2],
    energy_increase_tolerance_j: f64,
) -> Result<AdaptiveSurfaceSurfactantReport, SurfaceDynamicsError> {
    let [Some(diffusion_before), Some(diffusion_after)] = stages else {
        return Err(SurfaceDynamicsError::Invalid(
            "adaptive transport requires both stage diagnostics",
        ));
    };
    // Numerical energy allowances, like amount allowances, are not doubled
    // just because the same physical diffusion interval straddles mechanics.
    if -coupled.diffusion_free_energy_dissipation_j > energy_increase_tolerance_j {
        return Err(SurfaceDiffusionError::FreeEnergyIncrease {
            change_j: -coupled.diffusion_free_energy_dissipation_j,
        }
        .into());
    }
    let total_trials = diffusion_before
        .trials
        .checked_add(diffusion_after.trials)
        .ok_or(SurfaceDynamicsError::Invalid(
            "adaptive trial count overflow",
        ))?;
    let total_evaluated_pair_exchanges = diffusion_before
        .evaluated_pair_exchanges
        .checked_add(diffusion_after.evaluated_pair_exchanges)
        .ok_or(SurfaceDynamicsError::Invalid(
            "adaptive pair count overflow",
        ))?;
    Ok(AdaptiveSurfaceSurfactantReport {
        coupled,
        diffusion_before,
        diffusion_after,
        total_trials,
        total_evaluated_pair_exchanges,
    })
}

impl CapillarySurfaceDynamics<'_> {
    /// Evolve the native mechanical state and its insoluble coating atomically.
    ///
    /// Material advection and concentration-dependent stress use the original
    /// exact finite-area force and tangent. Diffusion changes the same face
    /// amounts on the actual initial/final geometry. Its graph is homogeneous
    /// in temperature and chemical law. Heterogeneous coatings use
    /// [`Self::advance_partitioned_surfactant_diffusive`].
    ///
    /// The symmetric composition is second order for smooth admitted states and
    /// a held graph, with the diffusion substep limit controlling edge splitting.
    /// Packing, negative tension, any solve/balance failure or cancellation
    /// leaves BOTH the original trajectory and coating unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_surfactant_diffusive(
        &self,
        state: &mut SecondOrderState,
        amounts_mol: &mut [f64],
        transport: SurfaceSurfactantTransport<'_>,
        viscosity: SurfaceViscosityField<'_>,
        external_forces_n: &[Vector],
        dt_s: f64,
        controls: SurfaceDynamicsControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<SurfaceSurfactantReport, SurfaceDynamicsError> {
        self.advance_with_surfactant_transport(
            state,
            amounts_mol,
            Transport::Uniform(transport),
            viscosity,
            external_forces_n,
            dt_s,
            controls,
            cancelled,
            |report, _| Ok(report),
        )
    }

    /// Diffuse and move a coating with distinct material affinities or site
    /// capacities through the same atomic mechanical composition.
    ///
    /// Each face retains its original law during both diffusion stages and the
    /// mechanical trial. The actual initial/final areas set concentrations and
    /// site capacities. The same laws supply chemical potentials, capillary
    /// stress and every surface-energy increment. Temperature, references and
    /// conductances remain held physical inputs throughout this step.
    ///
    /// The partitioned diffusion admission rules apply independently to each
    /// conducting connection. Any refusal, including final-geometry diffusion,
    /// publishes neither trajectory nor amounts.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_partitioned_surfactant_diffusive(
        &self,
        state: &mut SecondOrderState,
        amounts_mol: &mut [f64],
        transport: SurfacePartitionedSurfactantTransport<'_>,
        viscosity: SurfaceViscosityField<'_>,
        external_forces_n: &[Vector],
        dt_s: f64,
        controls: SurfaceDynamicsControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<SurfaceSurfactantReport, SurfaceDynamicsError> {
        self.advance_with_surfactant_transport(
            state,
            amounts_mol,
            Transport::Partitioned(transport),
            viscosity,
            external_forces_n,
            dt_s,
            controls,
            cancelled,
            |report, _| Ok(report),
        )
    }

    /// Move a homogeneous coating while adaptively resolving diffusion on the
    /// actual initial and solved final areas. Both halves share one trial/pair
    /// budget, and only accepted fine trajectories contribute physical loss.
    ///
    /// This controls internal diffusion edge splitting. It does not estimate
    /// the mechanical or diffusion/motion composition error; refine `dt_s`
    /// separately. Any stage refusal publishes neither mechanics nor amounts.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_surfactant_diffusive_adaptive(
        &self,
        state: &mut SecondOrderState,
        amounts_mol: &mut [f64],
        transport: AdaptiveSurfaceSurfactantTransport<'_>,
        viscosity: SurfaceViscosityField<'_>,
        external_forces_n: &[Vector],
        dt_s: f64,
        controls: SurfaceDynamicsControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<AdaptiveSurfaceSurfactantReport, SurfaceDynamicsError> {
        self.advance_with_surfactant_transport(
            state,
            amounts_mol,
            Transport::AdaptiveUniform(transport),
            viscosity,
            external_forces_n,
            dt_s,
            controls,
            cancelled,
            |report, stages| {
                adaptive_report(
                    report,
                    stages,
                    transport.controls.energy_increase_tolerance_j,
                )
            },
        )
    }

    /// Adaptively resolve heterogeneous partition diffusion through the same
    /// original mechanical trial, original per-face laws and atomic account.
    ///
    /// Initial/final geometry supplies each half's capacities; authored laws
    /// and conductances remain held. Shared work limits include every coarse
    /// and rejected evaluation. The temporal estimator controls internal
    /// diffusion splitting only, not macro mechanical/composition accuracy.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_partitioned_surfactant_diffusive_adaptive(
        &self,
        state: &mut SecondOrderState,
        amounts_mol: &mut [f64],
        transport: AdaptiveSurfacePartitionedSurfactantTransport<'_>,
        viscosity: SurfaceViscosityField<'_>,
        external_forces_n: &[Vector],
        dt_s: f64,
        controls: SurfaceDynamicsControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<AdaptiveSurfaceSurfactantReport, SurfaceDynamicsError> {
        self.advance_with_surfactant_transport(
            state,
            amounts_mol,
            Transport::AdaptivePartitioned(transport),
            viscosity,
            external_forces_n,
            dt_s,
            controls,
            cancelled,
            |report, stages| {
                adaptive_report(
                    report,
                    stages,
                    transport.controls.energy_increase_tolerance_j,
                )
            },
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn advance_with_surfactant_transport<R>(
        &self,
        state: &mut SecondOrderState,
        amounts_mol: &mut [f64],
        transport: Transport<'_>,
        viscosity: SurfaceViscosityField<'_>,
        external_forces_n: &[Vector],
        dt_s: f64,
        controls: SurfaceDynamicsControls,
        mut cancelled: impl FnMut() -> bool,
        finish: impl FnOnce(
            SurfaceSurfactantReport,
            [Option<AdaptiveSurfaceDiffusionReport>; 2],
        ) -> Result<R, SurfaceDynamicsError>,
    ) -> Result<R, SurfaceDynamicsError> {
        if cancelled() {
            return Err(SurfaceDynamicsError::Cancelled);
        }
        if state.q.len() != 3 * self.masses_kg.len()
            || amounts_mol.len() != self.faces.len()
            || !dt_s.is_finite()
            || dt_s <= 0.0
        {
            return Err(SurfaceDynamicsError::Invalid(
                "complete mechanical coordinates, one amount per face and a positive finite step required",
            ));
        }
        let (maximum_faces, maximum_edges, amount_tolerance_mol, edges, laws) = match transport {
            Transport::Uniform(data) => (
                data.controls.maximum_faces,
                data.controls.maximum_edges,
                data.controls.amount_tolerance_mol,
                data.edges.len(),
                None,
            ),
            Transport::Partitioned(data) => (
                data.controls.maximum_faces,
                data.controls.maximum_edges,
                data.controls.amount_tolerance_mol,
                data.edges.len(),
                Some(data.laws),
            ),
            Transport::AdaptiveUniform(data) => (
                data.controls.maximum_faces,
                data.controls.maximum_edges,
                data.controls.amount_tolerance_mol,
                data.edges.len(),
                None,
            ),
            Transport::AdaptivePartitioned(data) => (
                data.controls.maximum_faces,
                data.controls.maximum_edges,
                data.controls.amount_tolerance_mol,
                data.edges.len(),
                Some(data.laws),
            ),
        };
        if laws.is_some_and(|laws| laws.len() != self.faces.len()) {
            return Err(SurfaceDiffusionError::InvalidInput(
                "one original constitutive law required per face",
            )
            .into());
        }
        if self.faces.len() > maximum_faces || edges > maximum_edges {
            return Err(SurfaceDiffusionError::WorkBudget.into());
        }
        let before_areas = face_areas(self.faces, &state.q, &mut cancelled)?;
        let mut candidate_amounts = Vec::new();
        candidate_amounts
            .try_reserve_exact(amounts_mol.len())
            .map_err(|_| SurfaceDynamicsError::Allocation)?;
        candidate_amounts.extend_from_slice(amounts_mol);
        let mut work = AdaptiveWork::default();
        let before = transport.diffuse(
            &before_areas,
            &mut candidate_amounts,
            0.5 * dt_s,
            SurfaceDiffusionPhase::BeforeMotion,
            &mut work,
            &mut cancelled,
        )?;
        let (candidate, mechanics) = {
            let mut coating_faces = Vec::new();
            let coating = match transport {
                Transport::Uniform(SurfaceSurfactantTransport { law, .. })
                | Transport::AdaptiveUniform(AdaptiveSurfaceSurfactantTransport { law, .. }) => {
                    SurfactantField::Uniform {
                        law,
                        amounts_mol: &candidate_amounts,
                    }
                }
                Transport::Partitioned(SurfacePartitionedSurfactantTransport { laws, .. })
                | Transport::AdaptivePartitioned(AdaptiveSurfacePartitionedSurfactantTransport {
                    laws,
                    ..
                }) => {
                    coating_faces
                        .try_reserve_exact(self.faces.len())
                        .map_err(|_| SurfaceDynamicsError::Allocation)?;
                    for (face, (&law, &amount_mol)) in
                        laws.iter().zip(&candidate_amounts).enumerate()
                    {
                        if face % 64 == 0 && cancelled() {
                            return Err(SurfaceDynamicsError::Cancelled);
                        }
                        coating_faces.push(SurfactantFace { law, amount_mol });
                    }
                    SurfactantField::PerFace(&coating_faces)
                }
            };
            self.surface_trial(
                state,
                SurfaceForces::Insoluble(coating),
                viscosity,
                external_forces_n,
                dt_s,
                controls,
                &mut cancelled,
            )?
        };
        let after_areas = face_areas(self.faces, &candidate.q, &mut cancelled)?;
        let after = transport.diffuse(
            &after_areas,
            &mut candidate_amounts,
            0.5 * dt_s,
            SurfaceDiffusionPhase::AfterMotion,
            &mut work,
            &mut cancelled,
        )?;
        let diffusion_before = before.physical;
        let diffusion_after = after.physical;

        let mut amount_defect = Sum::default();
        for (face, (&before, &after)) in amounts_mol.iter().zip(&candidate_amounts).enumerate() {
            if face % 64 == 0 && cancelled() {
                return Err(SurfaceDynamicsError::Cancelled);
            }
            amount_defect.add(after - before);
        }
        let amount_defect = amount_defect.value();
        if !amount_defect.is_finite() {
            return Err(SurfaceDiffusionError::Unrepresentable.into());
        }
        if amount_defect.abs() > amount_tolerance_mol {
            return Err(SurfaceDiffusionError::AmountBalance {
                defect_mol: amount_defect,
            }
            .into());
        }
        let mut surface_change = Sum::default();
        for change in [
            diffusion_before.free_energy_change_j,
            mechanics.surface.energy_change_j,
            diffusion_after.free_energy_change_j,
        ] {
            surface_change.add(change);
        }
        let surface_change = surface_change.value();
        let diffusion_loss =
            diffusion_before.free_energy_dissipation_j + diffusion_after.free_energy_dissipation_j;
        let mut energy = Sum::default();
        for value in [
            mechanics.kinetic_change_j,
            surface_change,
            mechanics.surface_viscous_dissipation_j,
            diffusion_loss,
            -mechanics.external_work_j,
        ] {
            energy.add(value);
        }
        let energy = energy.value();
        if !surface_change.is_finite()
            || !diffusion_loss.is_finite()
            || !energy.is_finite()
            || energy.abs() > controls.energy_tolerance_j
        {
            return Err(SurfaceDynamicsError::Balance {
                impulse_n_s: mechanics.impulse_residual_n_s,
                kinematic_m: mechanics.kinematic_residual_m,
                energy_j: energy,
            });
        }
        if cancelled() {
            return Err(SurfaceDynamicsError::Cancelled);
        }
        let report = SurfaceSurfactantReport {
            surface_free_energy_before_j: diffusion_before.free_energy_before_j,
            surface_free_energy_after_j: diffusion_after.free_energy_after_j,
            surface_free_energy_change_j: surface_change,
            diffusion_free_energy_dissipation_j: diffusion_loss,
            amount_defect_mol: amount_defect,
            energy_residual_j: energy,
            mechanics,
            diffusion_before,
            diffusion_after,
        };
        let result = finish(report, [before.adaptive, after.adaptive])?;
        publish_state(state, candidate);
        amounts_mol.copy_from_slice(&candidate_amounts);
        Ok(result)
    }
}

fn face_areas(
    faces: &[[usize; 3]],
    coordinates: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<f64>, SurfaceDynamicsError> {
    let mut areas = Vec::new();
    areas
        .try_reserve_exact(faces.len())
        .map_err(|_| SurfaceDynamicsError::Allocation)?;
    for (face, &indices) in faces.iter().enumerate() {
        if face % 64 == 0 && cancelled() {
            return Err(SurfaceDynamicsError::Cancelled);
        }
        let vertices = indices.map(|i| core::array::from_fn(|d| coordinates[3 * i + d]));
        areas.push(geometry(vertices)?.area);
    }
    Ok(areas)
}
