//! Spatial vector hysteresis and resolved stationary eddy-current diffusion.
//!
//! Backward Euler minimizes the original frozen-history magnetic potential
//! plus `delta A . C_sigma delta A / (2 dt) - load . delta A`. The existing
//! floating-conductor projection, Newton/FEEC/checked-CG solve and exact P1
//! Joule integral retain ownership. No phenomenological core-loss term is added.
//!
//! Every trial starts from the same accepted vector polarizations. The physical
//! balance is `delta W + Q_pinning + Q_Joule + D_numerical = W_supplied`.
//! The endpoint numerical remainder is never heat. This is first-order,
//! fixed-geometry, isotropic magnetoquasistatics with fixed boundary traces and
//! explicit balanced impressed currents; no motion, terminal circuit, changing
//! magnetic coefficients, anisotropic fit or experimental validation is implied.

use super::{
    EddyCellHeat, EddyConductor, EddyError, EddyHeat, EddySection, admit_algebraic,
    algebraic_blocks, bind, finite, heat_cells, load,
    nonlinear::{BackwardConduction, NonlinearEddyControls, admit_controls},
    poll, positive, regions,
};
use crate::planar::{
    PlanarMagnetostatics,
    nonlinear::hysteretic::{
        PlanarHystereticControls, PlanarHystereticMagnetostatics, PlanarHystereticState,
        PlanarHystereticStep,
    },
};
use fs_material::{
    field_storage::admit_source_temperature_refresh, state_point::ResolvedMaterialStatePoint,
};
use fs_rep_mesh::{Metric2, TriComplex2};

/// Original induction and material/field refusals, including cancellation.
pub use super::nonlinear::NonlinearEddyError as HystereticEddyError;
type Result<T> = core::result::Result<T, HystereticEddyError>;

/// The magnetic work budget also bounds the complete induction account and
/// agreement between the weak conductivity operator and physical Joule heat.
#[derive(Clone, Copy, Debug)]
pub struct HystereticEddyControls {
    pub magnetic: PlanarHystereticControls,
    pub current_balance_tolerance_a: f64,
}

fn controls(value: HystereticEddyControls) -> Result<()> {
    admit_controls(NonlinearEddyControls {
        nonlinear: value.magnetic.nonlinear,
        energy_tolerance_j: value.magnetic.energy_tolerance_j,
        current_balance_tolerance_a: value.current_balance_tolerance_a,
    })?;
    Ok(())
}

impl EddySection for PlanarHystereticMagnetostatics {
    fn mesh(&self) -> &TriComplex2 {
        self.mesh()
    }
    fn thickness_m(&self) -> f64 {
        match self.mesh().metric() {
            Metric2::Planar { thickness } => thickness,
            _ => unreachable!("the magnetic owner admits only a planar metric"),
        }
    }
    fn admit_conductor(
        &self,
        triangle: usize,
        material: &ResolvedMaterialStatePoint,
    ) -> core::result::Result<(), EddyError> {
        let original = &self.materials()[triangle];
        if original.material().card_identity() != material.card_identity()
            || original.material().material() != material.material()
        {
            return Err(EddyError::Invalid(
                "magnetic fit and conductor must retain the same original material card",
            ));
        }
        // Resolve only the original fitted requirements. The electrical bundle
        // may have different properties, but cannot replace magnetic evidence.
        let checked = original.at_query_point(&material.query())?;
        admit_source_temperature_refresh(original.material(), checked.material())?;
        for ((axis, before), (_, after)) in original
            .material()
            .query_point()
            .iter()
            .zip(checked.material().query_point())
        {
            if before.to_bits() != after.to_bits() {
                original
                    .admit_temperature_interval(axis, [before.min(*after), before.max(*after)])?;
            }
        }
        Ok(())
    }
}

/// Accepted field, every vector polarization, conductor sources and one clock.
/// Cloning makes a checkpoint; proposing a step never mutates this state.
#[derive(Clone)]
pub struct PlanarHystereticEddy<'m> {
    magnetic: PlanarHystereticState<'m>,
    conductors: Vec<EddyConductor>,
    insulators: Vec<usize>,
    boundary: Vec<(usize, f64)>,
    components: Vec<usize>,
    algebraic: Vec<Vec<usize>>,
    impressed_jz_a_m2: Vec<f64>,
    time_s: f64,
}

/// The original magnetic and electrical receipts, with one physical heat sum.
#[derive(Clone, Debug)]
pub struct HystereticEddyStep {
    pub magnetic: PlanarHystereticStep,
    /// Exact spatial conductivity heat, retaining each original material source.
    pub cells: Vec<EddyCellHeat>,
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    pub joule_heat_j: f64,
    pub hysteresis_heat_j: f64,
    /// Q_Joule + Q_pinning, deposited once by the finite thermal receiver.
    pub heat_j: f64,
    pub supplied_work_j: f64,
    /// Endpoint integration loss, excluded from heat_j and thermal deposition.
    pub numerical_dissipation_j: f64,
    /// Stable delta W + physical heat + numerical loss - supplied work.
    pub energy_residual_j: f64,
    pub region_currents_a: Vec<(u32, f64)>,
}

pub struct PlanarHystereticEddyCandidate<'m> {
    state: PlanarHystereticEddy<'m>,
    step: HystereticEddyStep,
}
impl<'m> PlanarHystereticEddyCandidate<'m> {
    pub fn state(&self) -> &PlanarHystereticEddy<'m> {
        &self.state
    }
    pub fn step(&self) -> &HystereticEddyStep {
        &self.step
    }
    pub fn into_parts(self) -> (PlanarHystereticEddy<'m>, HystereticEddyStep) {
        (self.state, self.step)
    }
}

impl<'m> PlanarHystereticEddy<'m> {
    /// Start a transient from an actual prepared magnetic history. A released
    /// source can immediately drive conductor currents; insulating/floating
    /// algebraic constraints must already hold without an initialization impulse.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        initial: &PlanarHystereticState<'m>,
        conductors: Vec<EddyConductor>,
        insulating_triangles: &[usize],
        boundary: &[(usize, f64)],
        time_s: f64,
        limits: HystereticEddyControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        Self::new_with_impressed_current(
            initial,
            conductors,
            insulating_triangles,
            boundary,
            &vec![0.0; initial.problem().mesh().faces().len()],
            time_s,
            limits,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub fn new_with_impressed_current(
        initial: &PlanarHystereticState<'m>,
        conductors: Vec<EddyConductor>,
        insulating_triangles: &[usize],
        boundary: &[(usize, f64)],
        impressed_jz_a_m2: &[f64],
        time_s: f64,
        limits: HystereticEddyControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(cancelled)?;
        controls(limits)?;
        if !time_s.is_finite() || time_s < 0.0 {
            return Err(EddyError::Invalid("finite nonnegative initial time required").into());
        }
        let problem = initial.problem();
        let prescribed =
            PlanarMagnetostatics::admit_mesh_boundary(problem.mesh(), boundary, cancelled)?;
        if prescribed
            .iter()
            .zip(initial.potential_wb_m())
            .any(|(value, old)| value.is_some_and(|value| value != *old))
        {
            return Err(
                EddyError::Invalid("fixed induction boundary must match accepted field").into(),
            );
        }
        let (conductors, _) = bind(problem, conductors, None, cancelled)?;
        let components = regions(problem, &conductors, insulating_triangles, cancelled)?;
        let algebraic = algebraic_blocks(problem, &conductors, &prescribed, cancelled)?;
        let loads = load(
            problem,
            impressed_jz_a_m2,
            &components,
            limits.current_balance_tolerance_a,
            cancelled,
        )?;
        let magnetic = problem.magnetic_load(initial, cancelled)?;
        admit_algebraic(
            problem,
            &magnetic,
            &loads,
            &algebraic,
            limits.current_balance_tolerance_a,
            cancelled,
        )?;
        poll(cancelled)?;
        Ok(Self {
            magnetic: initial.clone(),
            conductors,
            insulators: insulating_triangles.to_vec(),
            boundary: boundary.to_vec(),
            components,
            algebraic,
            impressed_jz_a_m2: impressed_jz_a_m2.to_vec(),
            time_s,
        })
    }

    pub fn problem(&self) -> &'m PlanarHystereticMagnetostatics {
        self.magnetic.problem()
    }
    pub fn magnetic_state(&self) -> &PlanarHystereticState<'m> {
        &self.magnetic
    }
    pub fn conductors(&self) -> &[EddyConductor] {
        &self.conductors
    }
    pub fn insulating_triangles(&self) -> &[usize] {
        &self.insulators
    }
    pub fn potential_wb_m(&self) -> &[f64] {
        self.magnetic.potential_wb_m()
    }
    pub fn impressed_current_density_a_m2(&self) -> &[f64] {
        &self.impressed_jz_a_m2
    }
    pub fn energy_j(&self) -> f64 {
        self.magnetic.stored_energy_j()
    }
    pub fn time_s(&self) -> f64 {
        self.time_s
    }

    /// Resolve an actual backward-Euler field increment with frozen accepted
    /// magnetic memory and source-pinned endpoint conductivity. Each electrical
    /// region remains floating; impressed current has an explicit balanced return.
    /// Constitutive support, algebraic constraints, current and energy balances
    /// all pass before the complete candidate can be retained.
    pub fn step_to_candidate(
        &self,
        target_s: f64,
        conductors: Vec<EddyConductor>,
        impressed_jz_a_m2: &[f64],
        limits: HystereticEddyControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PlanarHystereticEddyCandidate<'m>> {
        poll(cancelled)?;
        controls(limits)?;
        let dt = target_s - self.time_s;
        if !target_s.is_finite() || !positive(dt) {
            return Err(
                EddyError::Invalid("positive representable induction interval required").into(),
            );
        }
        let problem = self.problem();
        let (conductors, conduction) =
            bind(problem, conductors, Some(&self.conductors), cancelled)?;
        let loads = load(
            problem,
            impressed_jz_a_m2,
            &self.components,
            limits.current_balance_tolerance_a,
            cancelled,
        )?;
        let quadratic = BackwardConduction {
            accepted: self.potential_wb_m(),
            conduction: &conduction,
            dt,
        };
        let (magnetic, magnetic_step) = problem
            .candidate_loaded(
                &self.magnetic,
                &loads,
                &self.boundary,
                self.potential_wb_m(),
                Some(&quadratic),
                limits.magnetic,
                cancelled,
            )?
            .into_parts();
        let magnetic_load = problem.magnetic_load(&magnetic, cancelled)?;
        admit_algebraic(
            problem,
            &magnetic_load,
            &loads,
            &self.algebraic,
            limits.current_balance_tolerance_a,
            cancelled,
        )?;
        let increments = magnetic
            .potential_wb_m()
            .iter()
            .zip(self.potential_wb_m())
            .map(|(next, old)| finite(next - old))
            .collect::<core::result::Result<Vec<_>, _>>()?;
        let EddyHeat {
            cells,
            heat_j: joule_heat_j,
            region_currents_a,
        } = heat_cells(
            &conduction,
            &conductors,
            &increments,
            dt,
            problem.thickness_m(),
            limits.current_balance_tolerance_a,
            cancelled,
        )?;
        let hysteresis_heat_j = magnetic_step.irreversible_heat_j;
        let heat_j = finite(joule_heat_j + hysteresis_heat_j)?;
        let numerical_dissipation_j = magnetic_step.numerical_dissipation_j;
        let supplied_work_j = magnetic_step.supplied_work_j;
        let energy_residual_j = finite(
            magnetic_step.stored_energy_change_j + heat_j + numerical_dissipation_j
                - supplied_work_j,
        )?;
        if (joule_heat_j - magnetic_step.quadratic_work_j).abs()
            > limits.magnetic.energy_tolerance_j
            || energy_residual_j.abs() > limits.magnetic.energy_tolerance_j
        {
            return Err(EddyError::Balance(
                "hysteretic storage, pinning heat, spatial Joule heat, numerical loss and supplied work do not close",
            ).into());
        }
        let step = HystereticEddyStep {
            stored_before_j: self.energy_j(),
            stored_after_j: magnetic.stored_energy_j(),
            magnetic: magnetic_step,
            cells,
            joule_heat_j,
            hysteresis_heat_j,
            heat_j,
            supplied_work_j,
            numerical_dissipation_j,
            energy_residual_j,
            region_currents_a,
        };
        let state = Self {
            magnetic,
            conductors,
            insulators: self.insulators.clone(),
            boundary: self.boundary.clone(),
            components: self.components.clone(),
            algebraic: self.algebraic.clone(),
            impressed_jz_a_m2: impressed_jz_a_m2.to_vec(),
            time_s: target_s,
        };
        poll(cancelled)?;
        Ok(PlanarHystereticEddyCandidate { state, step })
    }
}
