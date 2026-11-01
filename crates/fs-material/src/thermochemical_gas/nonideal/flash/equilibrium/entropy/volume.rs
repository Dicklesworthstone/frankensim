//! Original-source entropy/volume equilibrium for a conserved mixture feed.
//!
//! At each trial pressure the shared PS owner recovers the original temperature
//! and phase amounts. The existing caloric pressure bracket then closes volume.
//! A reference sample is evaluated once; trials reuse it without rebuilding full
//! source receipts or subtracting large absolute entropy/energy references.
//! No phase-dependent gamma, frozen heat capacity or acoustic closure is used.

use super::super::volume::{
    MAX_ENERGY_VOLUME_PRESSURE_EVALUATIONS, MAX_ENERGY_VOLUME_PT_EVALUATIONS, PressureSample,
    VolumeRootControls, Work, solve_volume_bracket,
};
use super::{PhaseEntropyControls, ResolvedThermochemicalPhaseState};
use crate::thermochemical_gas::nonideal::ThermochemicalPengRobinsonModel;
use crate::thermochemical_gas::nonideal::flash::{
    E,
    enthalpy::{CaloricCoordinate, CaloricRootControls, EnthalpyWork, Target, spacing},
    poll, positive,
};
use fs_qty::{Pressure, Temperature};
use fs_thermochem::eos::peng_robinson::phase_screen::PengRobinsonPhaseScreenV1;

type EntropyTarget = Target<PengRobinsonPhaseScreenV1>;

/// Shared hard ceiling on pressure inversions in one caloric/volume solve.
pub const MAX_ENTROPY_VOLUME_PRESSURE_EVALUATIONS: usize = MAX_ENERGY_VOLUME_PRESSURE_EVALUATIONS;
/// Shared hard ceiling on all original PT selections, including the reference.
pub const MAX_ENTROPY_VOLUME_PT_EVALUATIONS: usize = MAX_ENERGY_VOLUME_PT_EVALUATIONS;

/// Explicit source rectangle, entropy/volume accuracy and nested work limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhaseEntropyVolumeControls {
    /// Increasing positive pressures enclosing a positive-to-negative volume residual.
    pub pressure_bracket: [Pressure; 2],
    /// Absolute mass-specific volume residual tolerance [m3/kg].
    pub specific_volume_tolerance_m3_kg: f64,
    /// Pressure correction or remaining enclosing-bracket tolerance [Pa].
    pub pressure_tolerance_pa: f64,
    /// Whole-call number of original PS inversions, including both endpoints.
    pub maximum_pressure_evaluations: usize,
    /// Whole-call original PT selections across all nested solves and the reference.
    pub maximum_pt_evaluations: usize,
    /// Original entropy tolerance, source temperature interval and phase policy.
    pub phase: PhaseEntropyControls,
}
impl PhaseEntropyVolumeControls {
    /// Supply the pressure bracket and original entropy/phase controls explicitly.
    #[must_use]
    pub fn new(pressure_bracket: [Pressure; 2], phase: PhaseEntropyControls) -> Self {
        Self {
            pressure_bracket,
            specific_volume_tolerance_m3_kg: 1e-12,
            pressure_tolerance_pa: 1e-3,
            maximum_pressure_evaluations: 64,
            maximum_pt_evaluations: 4096,
            phase,
        }
    }

    /// Admit every physical and work limit before source evaluation.
    pub fn validate(self) -> Result<(), E> {
        self.phase.validate()?;
        let [lo, hi] = self.pressure_bracket.map(|p| p.value());
        if !positive(lo)
            || !positive(hi)
            || lo >= hi
            || !positive(self.specific_volume_tolerance_m3_kg)
            || !positive(self.pressure_tolerance_pa)
            || self.maximum_pressure_evaluations == 0
            || self.maximum_pressure_evaluations > MAX_ENTROPY_VOLUME_PRESSURE_EVALUATIONS
            || self.maximum_pt_evaluations == 0
            || self.maximum_pt_evaluations > MAX_ENTROPY_VOLUME_PT_EVALUATIONS
        {
            return Err(E::Invalid("invalid bounded entropy/volume flash controls"));
        }
        Ok(())
    }

    /// Conservative total cubic EOS ceiling, distinct from actual PT calls.
    pub fn maximum_eos_evaluation_bound(self) -> Result<usize, E> {
        self.validate()?;
        self.phase
            .phase_options()
            .maximum_eos_evaluation_bound()
            .map_err(E::Phase)?
            .checked_mul(self.maximum_pt_evaluations)
            .ok_or(E::Invalid("entropy/volume EOS work bound overflow"))
    }
}

/// Accepted original-source S,V state and separately reconstructed energy changes.
#[derive(Clone, Debug, PartialEq)]
pub struct PhaseEntropyVolumeSolution {
    state: ResolvedThermochemicalPhaseState,
    entropy_residual_j_kg_k: f64,
    volume_residual_m3_kg: f64,
    accepted_increment_j_kg_k: Option<f64>,
    enthalpy_change_j_kg: Option<f64>,
    internal_energy_change_j_kg: Option<f64>,
    pressure_evaluations: usize,
    pt_evaluations: usize,
    eos_evaluation_upper_bound: usize,
}
impl PhaseEntropyVolumeSolution {
    /// Complete accepted state with original material and phase evidence.
    #[must_use]
    pub const fn state(&self) -> &ResolvedThermochemicalPhaseState {
        &self.state
    }
    /// Move the state without another source or phase evaluation.
    #[must_use]
    pub fn into_state(self) -> ResolvedThermochemicalPhaseState {
        self.state
    }
    /// Source entropy coordinate minus requested target [J/(kg K)].
    #[must_use]
    pub const fn entropy_residual_j_kg_k(&self) -> f64 {
        self.entropy_residual_j_kg_k
    }
    /// Actual original phase volume per mass minus target [m3/kg].
    #[must_use]
    pub const fn volume_residual_m3_kg(&self) -> f64 {
        self.volume_residual_m3_kg
    }
    /// Integrated original-source Delta s, only for the increment API.
    #[must_use]
    pub const fn accepted_increment_j_kg_k(&self) -> Option<f64> {
        self.accepted_increment_j_kg_k
    }
    /// Integrated original-source Delta h, only for the increment API.
    #[must_use]
    pub const fn enthalpy_change_j_kg(&self) -> Option<f64> {
        self.enthalpy_change_j_kg
    }
    /// Independently integrated original-source Delta u, only for the increment API.
    #[must_use]
    pub const fn internal_energy_change_j_kg(&self) -> Option<f64> {
        self.internal_energy_change_j_kg
    }
    /// Actual pressure/entropy inversions, including bracket endpoints.
    #[must_use]
    pub const fn pressure_evaluations(&self) -> usize {
        self.pressure_evaluations
    }
    /// Actual original PT selections, including the single reference selection.
    #[must_use]
    pub const fn pt_evaluations(&self) -> usize {
        self.pt_evaluations
    }
    /// Original screen work plus optional split-work ceilings; an upper bound.
    #[must_use]
    pub const fn eos_evaluation_upper_bound(&self) -> usize {
        self.eos_evaluation_upper_bound
    }
}

impl ThermochemicalPengRobinsonModel {
    /// Recover pressure, temperature and phases at specified specific entropy
    /// and volume. The declared source rectangle must admit every sampled state.
    /// Finite phase-screen evidence is retained without a global minimum or
    /// unrestricted uniqueness claim. Absolute entropy must resolve its tolerance.
    pub fn phase_equilibrium_entropy_volume(
        &self,
        specific_volume_m3_kg: f64,
        specific_entropy_j_kg_k: f64,
        controls: PhaseEntropyVolumeControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PhaseEntropyVolumeSolution, E> {
        self.admit_entropy_volume(specific_volume_m3_kg, controls, None, &mut cancelled)?;
        if !specific_entropy_j_kg_k.is_finite()
            || spacing(specific_entropy_j_kg_k) > controls.phase.entropy_tolerance_j_kg_k
        {
            return Err(E::Invalid(
                "absolute entropy cannot resolve tolerance; use source increments",
            ));
        }
        self.solve_entropy_volume(
            specific_volume_m3_kg,
            Target::Absolute(specific_entropy_j_kg_k),
            controls,
            Work::default(),
            &mut cancelled,
        )
    }

    /// Change entropy and volume at fixed complete species inventory.
    /// Zero Delta s gives the reversible adiabatic endpoint. Energy changes
    /// integrate the original NASA laws before phase-departure differences,
    /// so independent species entropy/formation datums cancel before rounding.
    /// This closure supplies no piston inertia, heat rate, phase kinetics or work
    /// quadrature. At zero increment and identical volume, an in-bracket reference
    /// returns exactly after one PT selection and no pressure inversion.
    #[allow(clippy::too_many_arguments)]
    pub fn phase_equilibrium_entropy_volume_increment(
        &self,
        initial_temperature: Temperature,
        initial_pressure: Pressure,
        final_specific_volume_m3_kg: f64,
        change_j_kg_k: f64,
        controls: PhaseEntropyVolumeControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PhaseEntropyVolumeSolution, E> {
        self.admit_entropy_volume(
            final_specific_volume_m3_kg,
            controls,
            Some(initial_temperature),
            &mut cancelled,
        )?;
        if !positive(initial_pressure.value()) || !change_j_kg_k.is_finite() {
            return Err(E::Invalid(
                "finite entropy increment and positive reference pressure required",
            ));
        }
        let mut work = Work::default();
        let initial = self.volume_caloric_sample(
            initial_temperature,
            initial_pressure,
            controls.phase.phase_options(),
            controls.maximum_pt_evaluations,
            &mut work,
            &mut cancelled,
        )?;
        self.solve_entropy_volume(
            final_specific_volume_m3_kg,
            Target::Increment {
                initial,
                requested: change_j_kg_k,
            },
            controls,
            work,
            &mut cancelled,
        )
    }

    fn admit_entropy_volume(
        &self,
        volume: f64,
        controls: PhaseEntropyVolumeControls,
        initial: Option<Temperature>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), E> {
        controls.validate()?;
        if !positive(volume) || spacing(volume) > controls.specific_volume_tolerance_m3_kg {
            return Err(E::Invalid(
                "positive representable specific volume required",
            ));
        }
        self.admit_caloric_flash(
            controls.pressure_bracket[0],
            controls.phase.temperature_bracket,
            initial,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn entropy_volume_pressure(
        &self,
        pressure: f64,
        volume: f64,
        target: &EntropyTarget,
        controls: PhaseEntropyVolumeControls,
        scalar: CaloricRootControls,
        work: &mut Work,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PressureSample, E> {
        poll(cancelled)?;
        if work.pressure >= controls.maximum_pressure_evaluations {
            return Err(E::PressureWorkLimit {
                pressure_evaluations: work.pressure,
            });
        }
        work.pressure += 1;
        let root = self
            .solve_flash_coordinate(
                Pressure::new(pressure),
                target,
                CaloricCoordinate::Entropy,
                scalar,
                EnthalpyWork::default(),
                &mut |t, p, c| {
                    self.volume_caloric_sample(
                        t,
                        p,
                        controls.phase.phase_options(),
                        controls.maximum_pt_evaluations,
                        work,
                        c,
                    )
                },
                cancelled,
            )
            .map_err(|error| match error {
                E::WorkLimit { .. } => E::WorkLimit {
                    flash_evaluations: work.pt,
                },
                other => other,
            })?;
        let actual_volume = root.sample.molar_volume / self.mechanics().molar_mass().value();
        if !positive(actual_volume)
            || matches!(target, Target::Absolute(_))
                && spacing(root.coordinate) > scalar.residual_tolerance
        {
            return Err(E::Invalid(
                "entropy/volume original-state resolution is insufficient",
            ));
        }
        Ok(PressureSample {
            sample: root.sample,
            coordinate: root.coordinate,
            volume_residual: actual_volume - volume,
            entropy_change: matches!(target, Target::Increment { .. }).then_some(root.coordinate),
        })
    }

    fn solve_entropy_volume(
        &self,
        volume: f64,
        target: EntropyTarget,
        controls: PhaseEntropyVolumeControls,
        mut work: Work,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PhaseEntropyVolumeSolution, E> {
        let scalar = match &target {
            Target::Absolute(_) => controls.phase.root_controls(),
            Target::Increment { requested, .. } => {
                controls.phase.increment_root_controls(*requested)?
            }
        };
        let mut volume_tolerance = controls.specific_volume_tolerance_m3_kg;
        if let Target::Increment { initial, requested } = &target {
            let change = volume - initial.molar_volume / self.mechanics().molar_mass().value();
            if *requested == 0.0
                && change == 0.0
                && initial.temperature.value() >= controls.phase.temperature_bracket[0].value()
                && initial.temperature.value() <= controls.phase.temperature_bracket[1].value()
                && initial.pressure.value() >= controls.pressure_bracket[0].value()
                && initial.pressure.value() <= controls.pressure_bracket[1].value()
            {
                return self.finish_entropy_volume(
                    PressureSample {
                        sample: initial.clone(),
                        coordinate: 0.0,
                        volume_residual: 0.0,
                        entropy_change: Some(0.0),
                    },
                    &target,
                    work,
                    cancelled,
                );
            }
            if change != 0.0 {
                volume_tolerance = volume_tolerance.min(0.25 * change.abs());
                if !positive(volume_tolerance) || spacing(volume) > volume_tolerance {
                    return Err(E::Invalid(
                        "entropy/volume change cannot resolve volume tolerance",
                    ));
                }
            }
        }
        let accepted = solve_volume_bracket(
            VolumeRootControls {
                pressure_bracket: controls.pressure_bracket,
                volume_tolerance,
                pressure_tolerance: controls.pressure_tolerance_pa,
                coordinate_target: target.requested(),
                coordinate_tolerance: scalar.residual_tolerance,
            },
            &mut |p, c| {
                self.entropy_volume_pressure(p, volume, &target, controls, scalar, &mut work, c)
            },
            cancelled,
        )?;
        self.finish_entropy_volume(accepted, &target, work, cancelled)
    }

    fn finish_entropy_volume(
        &self,
        accepted: PressureSample,
        target: &EntropyTarget,
        work: Work,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PhaseEntropyVolumeSolution, E> {
        poll(cancelled)?;
        let (enthalpy, internal_energy) = match target {
            Target::Absolute(_) => (None, None),
            Target::Increment { initial, .. } => (
                Some(self.caloric_enthalpy_change(initial, &accepted.sample)?),
                Some(self.caloric_internal_energy_change(initial, &accepted.sample)?),
            ),
        };
        let residual = accepted.coordinate - target.requested();
        let state = self.publish_caloric_phase(accepted.sample, cancelled)?;
        poll(cancelled)?;
        Ok(PhaseEntropyVolumeSolution {
            state,
            entropy_residual_j_kg_k: residual,
            volume_residual_m3_kg: accepted.volume_residual,
            accepted_increment_j_kg_k: accepted.entropy_change,
            enthalpy_change_j_kg: enthalpy,
            internal_energy_change_j_kg: internal_energy,
            pressure_evaluations: work.pressure,
            pt_evaluations: work.pt,
            eos_evaluation_upper_bound: work.eos_bound,
        })
    }
}
