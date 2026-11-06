//! Atomic first-order composition of prescribed advection and molecular mixing.
//!
//! Both operators act over the SAME physical interval: donor-cell advection is
//! followed by the original convex-split diffusion sweep. This Lie composition
//! is first order. Refining its macro interval is necessary to control splitting
//! error; internal advection CFL substeps do not establish coupled accuracy.
//! The final chemical/interfacial energy can increase under prescribed stirring.
//! Only the diffusion stage carries the original dissipation inequality; the
//! complete report does not identify energy change with heat or mechanical work.

/// Temporal refinement of the complete stirring and demixing operator.
pub mod adaptive;

use super::advection::{
    BinaryMiscibilityAdvection, BinaryMiscibilityAdvectionControls, BinaryMiscibilityAdvectionStep,
};
use super::{
    BinaryMiscibilityError as Error, BinaryMiscibilityNetwork, BinaryMiscibilitySolve,
    BinaryMiscibilityState, BinaryMiscibilityStep, Result, Sum, finite, poll,
};

/// Original stage reports and independently audited complete-interval balances.
#[derive(Clone, Debug, PartialEq)]
pub struct BinaryMiscibilityTransportStep {
    pub duration_s: f64,
    pub advection: BinaryMiscibilityAdvectionStep,
    pub diffusion: BinaryMiscibilityStep,
    /// Net component transfer in canonical face order, including both operators.
    pub transferred_component_moles: Vec<[f64; 2]>,
    pub component_residual_mol: [f64; 2],
    pub maximum_capacity_residual_mol: f64,
    pub free_energy_before_j: f64,
    pub free_energy_after_j: f64,
    /// Direct original-to-final increment, without subtraction of bulk totals.
    pub free_energy_change_j: f64,
    /// Total direct increment minus the sum of the two direct stage increments.
    /// Reported roundoff accounting residual, not a dissipation bound.
    pub stage_energy_residual_j: f64,
}

/// No failed stage or final audit publishes a partially transported field.
#[derive(Clone, Debug, PartialEq)]
pub struct BinaryMiscibilityTransportCandidate {
    state: BinaryMiscibilityState,
    step: BinaryMiscibilityTransportStep,
}
impl BinaryMiscibilityTransportCandidate {
    pub const fn state(&self) -> &BinaryMiscibilityState {
        &self.state
    }
    pub const fn step(&self) -> &BinaryMiscibilityTransportStep {
        &self.step
    }
    pub fn into_parts(self) -> (BinaryMiscibilityState, BinaryMiscibilityTransportStep) {
        (self.state, self.step)
    }
}

impl BinaryMiscibilityNetwork {
    /// Advect then diffuse over one physical interval, retaining the same model,
    /// graph, component inventories and compensated low words throughout.
    /// The smaller stage conservation allowance applies ONCE to the complete
    /// interval. Stage errors, cancellation or a final balance failure leave
    /// `initial` unchanged and return no intermediate state or partial receipt.
    pub fn transport_candidate(
        &self,
        initial: &BinaryMiscibilityState,
        flow: &BinaryMiscibilityAdvection,
        duration_s: f64,
        advection_controls: BinaryMiscibilityAdvectionControls,
        diffusion_controls: BinaryMiscibilitySolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BinaryMiscibilityTransportCandidate> {
        poll(&mut cancelled)?;
        let (mut advected, advection) = self
            .advect_candidate(
                initial,
                flow,
                duration_s,
                advection_controls,
                &mut cancelled,
            )?
            .into_parts();
        // Operator splitting applies two generators over one physical interval.
        // Only this private trial's clock is reset; its inventories are retained.
        advected.elapsed_s = initial.elapsed_s;
        let (state, diffusion) = self
            .step_candidate(&advected, duration_s, diffusion_controls, &mut cancelled)?
            .into_parts();
        let allowance = advection_controls
            .amount_balance_tolerance_mol
            .min(diffusion_controls.amount_balance_tolerance_mol);
        let mut balances = [Sum::default(); 2];
        let mut capacity_residual: f64 = 0.0;
        for (cell, (old, new)) in initial.amounts.iter().zip(&state.amounts).enumerate() {
            poll(&mut cancelled)?;
            for species in 0..2 {
                balances[species].add(finite(
                    new[species] - old[species],
                    "combined component balance",
                )?)?;
            }
            capacity_residual = capacity_residual.max(finite(
                (new[0] + new[1] - self.owner.capacities_mol[cell]).abs(),
                "combined capacity residual",
            )?);
        }
        let balances = [balances[0].value()?, balances[1].value()?];
        if balances.iter().any(|x| x.abs() > allowance) || capacity_residual > allowance {
            return Err(Error::Conservation {
                component_residual_mol: balances,
                capacity_residual_mol: capacity_residual,
            });
        }
        let mut transferred_component_moles = Vec::with_capacity(self.owner.faces.len());
        for (advective, &diffusive) in advection
            .transferred_component_moles
            .iter()
            .zip(&diffusion.transferred_a_mol)
        {
            poll(&mut cancelled)?;
            let mut pair = [Sum::default(); 2];
            for species in 0..2 {
                pair[species].add(advective[species])?;
                pair[species].add(if species == 0 { diffusive } else { -diffusive })?;
            }
            transferred_component_moles.push([pair[0].value()?, pair[1].value()?]);
        }
        // Immutable stage endpoints already supply these same boundary energies.
        let before = finite(advection.free_energy_before_j, "transport initial energy")?;
        let after = finite(diffusion.free_energy_after_j, "transport final energy")?;
        let change = self.free_energy_change(initial, &state, &mut cancelled)?;
        let mut stage_residual = Sum::default();
        stage_residual.add(change)?;
        stage_residual.add(-advection.free_energy_change_j)?;
        stage_residual.add(-diffusion.free_energy_change_j)?;
        let stage_energy_residual_j = stage_residual.value()?;
        poll(&mut cancelled)?;
        Ok(BinaryMiscibilityTransportCandidate {
            state,
            step: BinaryMiscibilityTransportStep {
                duration_s,
                advection,
                diffusion,
                transferred_component_moles,
                component_residual_mol: balances,
                maximum_capacity_residual_mol: capacity_residual,
                free_energy_before_j: before,
                free_energy_after_j: after,
                free_energy_change_j: change,
                stage_energy_residual_j,
            },
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::miscibility::advection::BinaryMiscibilityAdvectionFace;
    use crate::miscibility::{BinaryMiscibilityCell, BinaryMiscibilityFace};
    use fs_qty::{AmountConcentration, Temperature};
    use fs_thermochem::SpeciesId;
    use fs_thermochem::solution::{BinaryRegularSolutionV1, gradient::BinarySquareGradientV1};

    fn fixture(
        diffusivity: f64,
        volumes: [f64; 2],
        rate: f64,
    ) -> (BinaryMiscibilityNetwork, BinaryMiscibilityAdvection) {
        let pair = BinaryRegularSolutionV1::new(
            [SpeciesId::new("A").unwrap(), SpeciesId::new("B").unwrap()],
            0.0,
            0.0,
            Temperature::new(250.0),
            Temperature::new(600.0),
        )
        .unwrap();
        let model = BinarySquareGradientV1::new(
            pair,
            Temperature::new(300.0),
            AmountConcentration::new(1.0),
            0.0,
        )
        .unwrap();
        let net = BinaryMiscibilityNetwork::new(
            model,
            volumes
                .into_iter()
                .map(|volume_m3| BinaryMiscibilityCell { volume_m3 })
                .collect(),
            vec![BinaryMiscibilityFace {
                cells: [0, 1],
                area_m2: 1.0,
                distance_m: 1.0,
                diffusivity_m2_s: diffusivity,
            }],
            1 << 20,
            || false,
        )
        .unwrap();
        let flow = net
            .prepare_advection(
                vec![BinaryMiscibilityAdvectionFace {
                    forward_m3_s: rate,
                    reverse_m3_s: rate,
                }],
                0.0,
                || false,
            )
            .unwrap();
        (net, flow)
    }
    fn state(net: &BinaryMiscibilityNetwork, x: [f64; 2]) -> BinaryMiscibilityState {
        net.state_from_mole_fractions(x.map(|x| [x, 1.0 - x]).to_vec(), || false)
            .unwrap()
    }
    fn near(actual: f64, expected: f64, tolerance: f64) {
        assert!(
            (actual - expected).abs() <= tolerance,
            "{actual:e} != {expected:e}"
        );
    }
    fn transport(
        net: &BinaryMiscibilityNetwork,
        state: &BinaryMiscibilityState,
        flow: &BinaryMiscibilityAdvection,
        dt: f64,
    ) -> BinaryMiscibilityTransportCandidate {
        net.transport_candidate(
            state,
            flow,
            dt,
            BinaryMiscibilityAdvectionControls::default(),
            BinaryMiscibilitySolve::default(),
            || false,
        )
        .unwrap()
    }

    #[test]
    fn zero_flow_is_original_diffusion_with_one_physical_clock() {
        let (net, flow) = fixture(0.7, [1.0; 2], 0.0);
        let initial = state(&net, [0.6, 0.4]);
        let direct = net
            .step_candidate(&initial, 0.25, BinaryMiscibilitySolve::default(), || false)
            .unwrap();
        let combined = transport(&net, &initial, &flow, 0.25);
        assert_eq!(combined.state(), direct.state());
        assert_eq!(&combined.step().diffusion, direct.step());
        assert_eq!(
            combined.step().advection.transferred_component_moles,
            [[0.0; 2]]
        );
        let next = transport(&net, combined.state(), &flow, 0.25);
        assert_eq!(next.state().elapsed_s(), 0.5);
    }

    #[test]
    fn zero_diffusivity_is_original_advection_without_a_second_time_increment() {
        let (net, flow) = fixture(0.0, [1.0; 2], 0.25);
        let initial = state(&net, [0.6, 0.4]);
        let direct = net
            .advect_candidate(
                &initial,
                &flow,
                0.7,
                BinaryMiscibilityAdvectionControls::default(),
                || false,
            )
            .unwrap();
        let combined = transport(&net, &initial, &flow, 0.7);
        assert_eq!(combined.state(), direct.state());
        assert_eq!(&combined.step().advection, direct.step());
        assert_eq!(combined.step().diffusion.transferred_a_mol, [0.0]);
        assert_eq!(combined.state().elapsed_s(), 0.7);
    }

    #[test]
    fn combined_small_mode_decay_converges_under_macro_time_refinement() {
        let (net, flow) = fixture(0.7, [1.0; 2], 0.3);
        let amplitude = 1e-4;
        // Independently linearized ideal two-cell mode: a'=-2(Q+D)a.
        let expected = 0.5 + amplitude * (-2.0_f64).exp();
        let mut errors = Vec::new();
        for steps in [5, 20, 80] {
            let mut current = state(&net, [0.5 + amplitude, 0.5 - amplitude]);
            for _ in 0..steps {
                current = transport(&net, &current, &flow, 1.0 / steps as f64)
                    .into_parts()
                    .0;
            }
            errors.push((current.mole_fractions()[0][0] - expected).abs());
            near(current.elapsed_s(), 1.0, 3e-14);
        }
        assert!(errors[1] < errors[0] / 3.0, "{errors:?}");
        assert!(errors[2] < errors[1] / 3.0, "{errors:?}");
        assert!(errors[2] < 2e-7, "{errors:?}");
    }

    #[test]
    fn final_receipts_close_species_capacities_and_stage_energy_changes() {
        let (net, flow) = fixture(0.7, [1.0, 2.0], 0.4);
        let initial = state(&net, [0.2, 0.8]);
        let combined = transport(&net, &initial, &flow, 0.4);
        let report = combined.step();
        for species in 0..2 {
            let transfer = report.transferred_component_moles[0][species];
            near(
                initial.component_moles()[0][species] - transfer,
                combined.state().component_moles()[0][species],
                2e-14,
            );
            near(
                initial.component_moles()[1][species] + transfer,
                combined.state().component_moles()[1][species],
                2e-14,
            );
            near(report.component_residual_mol[species], 0.0, 2e-14);
        }
        near(report.maximum_capacity_residual_mol, 0.0, 2e-14);
        near(
            report.free_energy_before_j,
            net.free_energy_j(&initial, || false).unwrap(),
            1e-11,
        );
        near(
            report.free_energy_after_j,
            net.free_energy_j(combined.state(), || false).unwrap(),
            1e-11,
        );
        near(report.stage_energy_residual_j, 0.0, 2e-10);
        near(
            report.free_energy_change_j,
            report.advection.free_energy_change_j + report.diffusion.free_energy_change_j,
            2e-10,
        );
        near(combined.state().elapsed_s(), 0.4, 1e-15);
    }

    #[test]
    fn diffusion_failure_and_late_cancellation_keep_the_original_state() {
        let (net, flow) = fixture(1.0, [1.0; 2], 0.1);
        let initial = state(&net, [0.9, 0.1]);
        let original = initial.clone();
        let mut advection_polls = 0;
        net.advect_candidate(
            &initial,
            &flow,
            1.0,
            BinaryMiscibilityAdvectionControls::default(),
            || {
                advection_polls += 1;
                false
            },
        )
        .unwrap();
        assert!(matches!(
            net.transport_candidate(
                &initial,
                &flow,
                1.0,
                BinaryMiscibilityAdvectionControls::default(),
                BinaryMiscibilitySolve {
                    maximum_root_iterations: 1,
                    ..BinaryMiscibilitySolve::default()
                },
                || false,
            ),
            Err(Error::Root { .. })
        ));
        let mut polls = 0;
        assert!(matches!(
            net.transport_candidate(
                &initial,
                &flow,
                1.0,
                BinaryMiscibilityAdvectionControls::default(),
                BinaryMiscibilitySolve::default(),
                || {
                    polls += 1;
                    polls > advection_polls + 8
                },
            ),
            Err(Error::Cancelled)
        ));
        assert!(polls > advection_polls);
        assert_eq!(initial, original);
    }
}
