//! Advected material shear memory with finite, locally partitioned friction heat.
//!
//! A declared uniform-pressure brush attaches the native two-region Iwan strains
//! to surface A material coordinates. Translation clips the original parcels;
//! surviving strains and kinetic branches are copied exactly, never averaged.
//! Fresh material enters unloaded. Elastic energy leaving the footprint is an
//! explicit mechanical boundary flux, not frictional heat or numerical damping.
//!
//! Transport precedes the local return maps (first-order Lie splitting). The
//! resulting piecewise-constant work density is integrated over each complete
//! linear thermal footprint path. This is a convergent time discretization of
//! the declared brush model, not exact finite-step residence-time integration.
//! Refine time independently of the original finite thermal cells. The pressure
//! is uniform; a reference ellipse supplies the native shear law's area and
//! total load, not an elliptical pressure field. Torsion, changing normal load,
//! adhesion, recirculated outflow memory and native body dynamics are absent.

use super::*;
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{
    NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, SolverRunProgress,
};
use fs_tribo::partial_slip::{
    GeneralizedWorkOwnership, PartialSlipCheckpoint, PartialSlipError, PartialSlipKinematics,
    PartialSlipLaw, PartialSlipState, PartialSlipStep, TangentFrame,
};
use std::cell::Cell;
use std::convert::Infallible;

/// One unchanged physical material parcel, measured in surface A coordinates.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialSlipParcel {
    pub span_m: [f64; 2],
    state: PartialSlipState,
}
impl MaterialSlipParcel {
    pub const fn state(&self) -> &PartialSlipState {
        &self.state
    }
}

/// Complete accepted history on one unchanged strip route and clock.
#[derive(Clone)]
pub struct MovingPartialSlip {
    route: MovingFrictionRoute,
    reference: PartialSlipCheckpoint,
    parcels: Vec<MaterialSlipParcel>,
    position: MovingContactPosition,
    time_s: f64,
    frame: TangentFrame,
}

/// Candidate parcel geometry. Each region lies within one ORIGINAL A cell;
/// a source may therefore query that material's actual endpoint temperature.
#[derive(Clone, Debug)]
pub struct MovingSlipRegion {
    pub material_span_m: [f64; 2],
    pub surface_a_cell: CellId,
    pub area_fraction: f64,
    state: PartialSlipState,
}
impl MovingSlipRegion {
    pub const fn accepted_state(&self) -> &PartialSlipState {
        &self.state
    }
}

/// Source laws at the candidate thermal endpoint, in region order. Only the
/// original law's static/kinetic coefficients may change. A source that needs
/// an interface temperature must resolve its spatial applicability itself;
/// this owner never invents an average of unequal interface temperatures.
pub struct MovingPartialSlipSource<P> {
    pub laws: Vec<PartialSlipLaw>,
    pub conductivity_w_m_k: BTreeMap<CellId, f64>,
    pub transport: StepCoefficients<P>,
}

#[derive(Clone, Copy, Debug)]
pub struct MovingPartialSlipControls {
    pub maximum_regions: usize,
    pub maximum_cells: usize,
    pub maximum_trials: usize,
    pub maximum_newton_iterations: usize,
    pub difference_step_j: f64,
    /// Numerical coordinate scale, never an additional heat capacity.
    pub temperature_scale_j_per_k: f64,
    pub heat: FrictionHeatTolerances,
    pub energy_tolerance_j: f64,
    pub newton: NewtonKrylovConfig,
}

#[derive(Debug)]
pub enum MovingPartialSlipError<E = Infallible> {
    Invalid(&'static str),
    Thermal(DepositionError),
    Contact(PartialSlipError),
    Source(E),
    Newton(NewtonError),
    NotConverged(NewtonReport),
    Interrupted(SolverRunProgress),
    TrialBudget,
    RegionBudget,
}
impl<E> From<DepositionError> for MovingPartialSlipError<E> {
    fn from(e: DepositionError) -> Self {
        Self::Thermal(e)
    }
}
impl<E> From<PartialSlipError> for MovingPartialSlipError<E> {
    fn from(e: PartialSlipError) -> Self {
        Self::Contact(e)
    }
}
impl<E> From<NewtonError> for MovingPartialSlipError<E> {
    fn from(e: NewtonError) -> Self {
        Self::Newton(e)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for MovingPartialSlipError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "moving material partial slip: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for MovingPartialSlipError<E> {}

/// Native response is normalized to the reference patch. Multiply every
/// extensive work/force by area_fraction; strains and utilization are unchanged.
pub struct MovingSlipRegionStep {
    pub region: MovingSlipRegion,
    pub response: PartialSlipStep,
}
pub struct MovingPartialSlipReport<P> {
    pub regions: Vec<MovingSlipRegionStep>,
    pub heat: MovingFrictionHeatReport,
    pub source: P,
    pub nonlinear: NewtonReport,
    pub trials: usize,
    pub force_on_a_n: [f64; 3],
    /// Moment of the distributed shear about StripFrame::origin_m, averaged
    /// over the same linear path as heat deposition. No local spin law is added.
    pub moment_on_a_about_origin_nm: [f64; 3],
    pub relative_slip_m: f64,
    pub mechanical_work_j: f64,
    pub storage_before_j: f64,
    pub storage_after_j: f64,
    pub elastic_outflow_j: f64,
    pub numerical_dissipation_j: f64,
    /// Delta(U_contact + H_thermal) + D_numerical + U_out - W - Q_other.
    pub total_energy_residual_j: f64,
}

pub struct MovingPartialSlipCandidate<'m, P> {
    contact: MovingPartialSlip,
    thermal: DepositionSimulation<'m>,
    report: MovingPartialSlipReport<P>,
}
impl<'m, P> MovingPartialSlipCandidate<'m, P> {
    pub const fn contact(&self) -> &MovingPartialSlip {
        &self.contact
    }
    pub const fn thermal(&self) -> &DepositionSimulation<'m> {
        &self.thermal
    }
    pub const fn report(&self) -> &MovingPartialSlipReport<P> {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        MovingPartialSlip,
        DepositionSimulation<'m>,
        MovingPartialSlipReport<P>,
    ) {
        (self.contact, self.thermal, self.report)
    }
}

struct Staged<P> {
    regions: Vec<MovingSlipRegionStep>,
    source: P,
    work: f64,
    storage: f64,
    numerical: f64,
    force: [f64; 3],
    moment: [f64; 3],
}

fn stored(law: &PartialSlipLaw, state: &PartialSlipState) -> Result<f64, DepositionError> {
    let p = law.parameters();
    let weights = [
        1.0 - p.partial_slip_hardening_fraction,
        p.partial_slip_hardening_fraction,
    ];
    let mut energy = 0.0;
    for (i, strain) in state.region_elastic_displacement_m().iter().enumerate() {
        energy = finite(
            energy
                + 0.5
                    * p.tangential_stiffness_n_per_m
                    * weights[i]
                    * finite(strain[0] * strain[0] + strain[1] * strain[1])?,
        )?;
    }
    Ok(energy)
}
fn compatible<E>(
    reference: &PartialSlipLaw,
    law: &PartialSlipLaw,
) -> Result<(), MovingPartialSlipError<E>> {
    let p = reference.parameters();
    let mut q = law.parameters();
    q.static_mu = p.static_mu;
    q.kinetic_mu = p.kinetic_mu;
    if p != q || reference.source_id() != law.source_id() || reference.model_id() != law.model_id()
    {
        return Err(MovingPartialSlipError::Invalid(
            "moving history must retain its original source and storage law",
        ));
    }
    Ok(())
}
fn numerical<E>(e: FiniteDifferenceError<MovingPartialSlipError<E>>) -> MovingPartialSlipError<E> {
    match e {
        FiniteDifferenceError::Source(e) => e,
        FiniteDifferenceError::Invalid(s) => MovingPartialSlipError::Invalid(s),
    }
}

impl MovingPartialSlip {
    /// Begin with the supplied uniform native preload, including both regional
    /// strains and kinetic branches. Its nominal area must equal the strip's
    /// area. Uniform normal pressure and fixed reference stiffness are explicit
    /// constitutive assumptions; no measured calibration is inferred.
    pub fn new(
        thermal: &DepositionSimulation<'_>,
        route: MovingFrictionRoute,
        reference: PartialSlipCheckpoint,
        position: MovingContactPosition,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, MovingPartialSlipError> {
        checkpoint()?;
        route.validate(thermal, &mut checkpoint)?;
        let frame = TangentFrame::new(route.frame.normal, route.frame.tangent)?;
        let patch = reference.normal_patch();
        let area =
            core::f64::consts::PI * patch.semi_axis_longitudinal_m() * patch.semi_axis_lateral_m();
        if !area.is_finite()
            || !close(area, route.footprint_area_m2())
            || reference.state().region_elastic_twist_rad() != [0.0; 2]
            || reference
                .state()
                .reference_frame()
                .is_some_and(|f| f != frame)
            || position
                .body_offsets_m
                .iter()
                .chain(std::iter::once(&position.footprint_offset_m))
                .any(|x| !x.is_finite())
        {
            return Err(MovingPartialSlipError::Invalid(
                "uniform shear brush requires the same reference area/frame and zero torsional history",
            ));
        }
        reference
            .law()
            .restore_checkpoint(patch, reference.interface(), &reference)?;
        let span = route
            .footprint_m
            .map(|x| x + position.footprint_offset_m - position.body_offsets_m[0]);
        for side in 0..2 {
            let surfaces = &route.surfaces[side];
            let offset = position.footprint_offset_m - position.body_offsets_m[side];
            if route.footprint_m[0] + offset < surfaces[0].face.span_m[0]
                || route.footprint_m[1] + offset > surfaces.last().unwrap().face.span_m[1]
            {
                return Err(MovingPartialSlipError::Invalid(
                    "initial footprint leaves the original material surface",
                ));
            }
        }
        let parcels = vec![MaterialSlipParcel {
            span_m: span,
            state: reference.state().clone(),
        }];
        Ok(Self {
            route,
            reference,
            parcels,
            position,
            time_s: thermal.time_s(),
            frame,
        })
    }
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub const fn position(&self) -> MovingContactPosition {
        self.position
    }
    pub fn parcels(&self) -> &[MaterialSlipParcel] {
        &self.parcels
    }
    pub const fn reference(&self) -> &PartialSlipCheckpoint {
        &self.reference
    }
    pub fn stored_energy_j(&self) -> Result<f64, DepositionError> {
        let length = self.route.footprint_m[1] - self.route.footprint_m[0];
        self.parcels.iter().try_fold(0.0, |sum, p| {
            finite(
                sum + (p.span_m[1] - p.span_m[0]) / length
                    * stored(self.reference.law(), &p.state)?,
            )
        })
    }

    /// Accept transport, original regional memory and all finite thermal cells
    /// together. Body offsets determine actual signed longitudinal slip; the
    /// footprint offset independently determines material residence. A region's
    /// source may use its exact original A-cell temperature. No averaged
    /// interface temperature is supplied as a constitutive material coordinate.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn step_candidate<'m, P, E>(
        &self,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        end: MovingContactPosition,
        ownership: &GeneralizedWorkOwnership,
        controls: MovingPartialSlipControls,
        energy_guess_j: &[f64],
        interface_temperature_guess_k: f64,
        resolve: impl Fn(
            &DepositionSimulation<'m>,
            &[MovingSlipRegion],
            &[SweptContactPatch],
            &[f64],
        ) -> Result<MovingPartialSlipSource<P>, E>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MovingPartialSlipCandidate<'m, P>, MovingPartialSlipError<E>> {
        let cancel = RefCell::new(cancelled);
        let check = || {
            if (*cancel.borrow_mut())() {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        };
        check()?;
        controls.newton.validate()?;
        let dt = target_s - self.time_s;
        if thermal.time_s() != self.time_s
            || !dt.is_finite()
            || dt <= 0.0
            || controls.maximum_regions == 0
            || controls.maximum_cells == 0
            || thermal.cells().len() > controls.maximum_cells
            || controls.maximum_trials == 0
            || controls.maximum_newton_iterations == 0
            || energy_guess_j.len() != thermal.cells().len()
            || energy_guess_j.iter().any(|x| !x.is_finite())
            || [
                controls.difference_step_j,
                controls.temperature_scale_j_per_k,
                controls.energy_tolerance_j,
                interface_temperature_guess_k,
            ]
            .iter()
            .any(|x| !x.is_finite() || *x <= 0.0)
        {
            return Err(MovingPartialSlipError::Invalid(
                "synchronized finite bounded moving contact solve required",
            ));
        }
        let path = MovingContactPath {
            start: self.position,
            end,
        };
        let heat_step = MovingFrictionHeatStep::new(thermal, target_s, &self.route, path, &check)?;
        let shift = finite(end.footprint_offset_m - end.body_offsets_m[0])?;
        let footprint = self.route.footprint_m.map(|x| x + shift);
        let length = self.route.footprint_m[1] - self.route.footprint_m[0];
        let old = [
            self.parcels[0].span_m[0],
            self.parcels.last().unwrap().span_m[1],
        ];
        let mut transported = Vec::new();
        let mut outflow = 0.0;
        for p in &self.parcels {
            check()?;
            let lo = p.span_m[0].max(footprint[0]);
            let hi = p.span_m[1].min(footprint[1]);
            let retained = (hi - lo).max(0.0);
            outflow = finite(
                outflow
                    + ((p.span_m[1] - p.span_m[0]) - retained) / length
                        * stored(self.reference.law(), &p.state)?,
            )?;
            if retained > 0.0 {
                transported.push(MaterialSlipParcel {
                    span_m: [lo, hi],
                    state: p.state.clone(),
                });
            }
        }
        if footprint[0] < old[0] {
            transported.insert(
                0,
                MaterialSlipParcel {
                    span_m: [footprint[0], old[0].min(footprint[1])],
                    state: PartialSlipState::zero(),
                },
            );
        }
        if footprint[1] > old[1] {
            transported.push(MaterialSlipParcel {
                span_m: [old[1].max(footprint[0]), footprint[1]],
                state: PartialSlipState::zero(),
            });
        }
        // Split only at actual material boundaries. No interpolation or spring
        // averaging occurs, including when a warm parcel revisits a thermal cell.
        let mut regions = Vec::new();
        let mut first_cell = 0;
        for parcel in transported {
            while first_cell < self.route.surfaces[0].len()
                && self.route.surfaces[0][first_cell].face.span_m[1] <= parcel.span_m[0]
            {
                first_cell += 1;
            }
            for face in self.route.surfaces[0].iter().skip(first_cell) {
                check()?;
                if face.face.span_m[0] >= parcel.span_m[1] {
                    break;
                }
                let span = [
                    parcel.span_m[0].max(face.face.span_m[0]),
                    parcel.span_m[1].min(face.face.span_m[1]),
                ];
                if span[1] <= span[0] {
                    continue;
                }
                if regions.len() >= controls.maximum_regions {
                    return Err(MovingPartialSlipError::RegionBudget);
                }
                regions.push(MovingSlipRegion {
                    material_span_m: span,
                    surface_a_cell: face.face.cell,
                    area_fraction: (span[1] - span[0]) / length,
                    state: parcel.state.clone(),
                });
            }
        }
        let before = self.stored_energy_j()?;
        let after_transport = regions.iter().try_fold(0.0, |sum, r| {
            finite(sum + r.area_fraction * stored(self.reference.law(), &r.state)?)
        })?;
        if !close(regions.iter().map(|r| r.area_fraction).sum(), 1.0)
            || finite(after_transport + outflow - before)?.abs() > controls.energy_tolerance_j
        {
            return Err(MovingPartialSlipError::Invalid(
                "material transport failed original area/storage conservation",
            ));
        }
        // Geometry-only quadrature is prepared once outside all Newton trials.
        let faces: [BTreeMap<_, _>; 2] = std::array::from_fn(|side| {
            self.route.surfaces[side]
                .iter()
                .map(|p| (p.face.cell, p.face.span_m))
                .collect()
        });
        let mut trajectories = Vec::with_capacity(heat_step.patches().len());
        for patch in heat_step.patches() {
            check()?;
            trajectories.push([
                edges(
                    faces[0][&patch.cells[0]],
                    path.start.body_offsets_m[0],
                    path.end.body_offsets_m[0],
                )?,
                edges(
                    faces[1][&patch.cells[1]],
                    path.start.body_offsets_m[1],
                    path.end.body_offsets_m[1],
                )?,
            ]);
        }
        let mut deposition = Vec::with_capacity(regions.len());
        let mut covered = vec![0.0; heat_step.patches().len()];
        for region in &regions {
            check()?;
            let foot = edges(
                region.material_span_m.map(|x| x - shift),
                path.start.footprint_offset_m,
                path.end.footprint_offset_m,
            )?;
            let fb = bounds(foot);
            let area = region.area_fraction * self.route.footprint_area_m2();
            let mut weights = Vec::new();
            let mut sum = 0.0;
            for (index, &[a, b]) in trajectories.iter().enumerate() {
                check()?;
                let ab = bounds(a);
                let bb = bounds(b);
                if ab[1] <= fb[0] || ab[0] >= fb[1] || bb[1] <= fb[0] || bb[0] >= fb[1] {
                    continue;
                }
                let overlap = finite(
                    self.route.frame.width_m
                        * swept_overlap([a[0], a[1], b[0], b[1], foot[0], foot[1]])?,
                )?;
                if overlap > 0.0 {
                    sum = finite(sum + overlap)?;
                    covered[index] = finite(covered[index] + overlap)?;
                    weights.push((index, overlap / area));
                }
            }
            if (sum - area).abs()
                > 1024.0 * f64::EPSILON * area * (trajectories.len() as f64).sqrt()
            {
                return Err(MovingPartialSlipError::Invalid(
                    "regional material work lost swept deposition area",
                ));
            }
            deposition.push(weights);
        }
        for (area, patch) in covered.into_iter().zip(heat_step.patches()) {
            if (area - patch.mean_area_m2).abs()
                > 1024.0
                    * f64::EPSILON
                    * self.route.footprint_area_m2()
                    * (regions.len() as f64).sqrt()
            {
                return Err(MovingPartialSlipError::Invalid(
                    "regional deposition overlaps or omits original contact coverage",
                ));
            }
        }
        let slip = finite(
            (end.body_offsets_m[0] - self.position.body_offsets_m[0])
                - (end.body_offsets_m[1] - self.position.body_offsets_m[1]),
        )?;
        let kinematics = PartialSlipKinematics {
            creepage: [slip / dt, 0.0],
            rolling_speed_mps: 1.0,
            torsional_spin_rad_per_s: 0.0,
            contact_spin_rad_per_s: 0.0,
            dt_s: dt,
        };
        let cells = heat_step.cell_ids().len();
        let patches = heat_step.patches().len();
        let trials = Cell::new(0);
        let evaluate = |x: &[f64]| {
            check()?;
            if trials.get() >= controls.maximum_trials {
                return Err(MovingPartialSlipError::TrialBudget);
            }
            trials.set(trials.get() + 1);
            if x.len() != cells + patches {
                return Err(MovingPartialSlipError::Invalid(
                    "incomplete moving contact iterate",
                ));
            }
            let temperatures: Vec<_> = x[cells..]
                .iter()
                .map(|v| interface_temperature_guess_k + v / controls.temperature_scale_j_per_k)
                .collect();
            heat_step.evaluate(
                &x[..cells],
                &temperatures,
                |endpoint, patches, temperatures| {
                    let MovingPartialSlipSource {
                        laws,
                        conductivity_w_m_k,
                        transport,
                    } = resolve(endpoint, &regions, patches, temperatures)
                        .map_err(MovingPartialSlipError::Source)?;
                    if laws.len() != regions.len() {
                        return Err(MovingPartialSlipError::Invalid(
                            "one original source law per material region required",
                        ));
                    }
                    let StepCoefficients {
                        sources_w,
                        conductances_w_k,
                        payload,
                    } = transport;
                    let mut staged = Staged {
                        regions: Vec::with_capacity(regions.len()),
                        source: payload,
                        work: 0.0,
                        storage: 0.0,
                        numerical: 0.0,
                        force: [0.0; 3],
                        moment: [0.0; 3],
                    };
                    let mut heat = vec![0.0; patches.len()];
                    for (i, (region, law)) in regions.iter().zip(laws).enumerate() {
                        check()?;
                        compatible(self.reference.law(), &law)?;
                        let response = law.advance(
                            self.reference.normal_patch(),
                            self.reference.interface(),
                            self.frame,
                            kinematics,
                            ownership,
                            &region.state,
                        )?;
                        let a = region.area_fraction;
                        staged.work = finite(
                            staged.work + a * response.generalized_work.work_into_interface_j,
                        )?;
                        staged.storage = finite(staged.storage + a * response.stored_energy_j)?;
                        staged.numerical = finite(
                            staged.numerical
                                + a * response.generalized_work.numerical_dissipation_j,
                        )?;
                        for (f, r) in staged.force.iter_mut().zip(response.wrench.force_n) {
                            *f = finite(*f + a * r)?;
                        }
                        let centre =
                            f64::midpoint(region.material_span_m[0], region.material_span_m[1])
                                - shift
                                + f64::midpoint(
                                    self.position.footprint_offset_m,
                                    end.footprint_offset_m,
                                );
                        let arm = self.route.frame.tangent.map(|t| centre * t);
                        let force = response.wrench.force_n.map(|f| a * f);
                        for (axis, moment) in staged.moment.iter_mut().enumerate() {
                            let j = (axis + 1) % 3;
                            let k = (axis + 2) % 3;
                            *moment = finite(*moment + arm[j] * force[k] - arm[k] * force[j])?;
                        }
                        let work = a * response.dissipation.heat_j;
                        for &(patch, weight) in &deposition[i] {
                            heat[patch] = finite(heat[patch] + weight * work)?;
                        }
                        staged.regions.push(MovingSlipRegionStep {
                            region: region.clone(),
                            response,
                        });
                    }
                    Ok(MovingFrictionInput {
                        dissipated_j: heat,
                        conductivity_w_m_k,
                        transport: StepCoefficients {
                            sources_w,
                            conductances_w_k,
                            payload: staged,
                        },
                    })
                },
                &check,
            )
        };
        let problem = FallibleFiniteDifference::new(
            cells + patches,
            controls.difference_step_j,
            |x: &[f64]| evaluate(x).map(|trial| trial.residual_j().to_vec()),
        )
        .map_err(numerical)?;
        let mut initial = energy_guess_j.to_vec();
        initial.resize(cells + patches, 0.0);
        let construction = NewtonKrylovState::new(&problem, initial, controls.newton);
        if let Some(error) = problem.take_failure() {
            return Err(numerical(error));
        }
        let mut solver = construction?;
        let result =
            solver.run_controlled(&problem, controls.maximum_newton_iterations, &mut || {
                (*cancel.borrow_mut())() || problem.has_failed()
            });
        if let Some(error) = problem.take_failure() {
            return Err(numerical(error));
        }
        check()?;
        if result.progress != SolverRunProgress::Complete {
            return Err(MovingPartialSlipError::Interrupted(result.progress));
        }
        if !result.report.converged {
            return Err(MovingPartialSlipError::NotConverged(result.report));
        }
        let trial = evaluate(&solver.x)?;
        let staged = trial.payload();
        let total = finite(
            staged.storage - before
                + trial.report().thermal.stored_energy_change_j
                + staged.numerical
                + outflow
                - staged.work
                - trial.report().external_heat_j,
        )?;
        if total.abs() > controls.energy_tolerance_j {
            return Err(MovingPartialSlipError::Invalid(
                "moving brush and original thermal storage do not balance",
            ));
        }
        let (thermal, staged, heat) = trial.accept(controls.heat, &check)?;
        let mut parcels: Vec<MaterialSlipParcel> = Vec::new();
        for r in &staged.regions {
            check()?;
            // Coalesce only exactly identical original histories. This has zero
            // numerical diffusion and does not conflate distinct kinetic branches.
            if let Some(last) = parcels.last_mut() {
                if last.span_m[1] == r.region.material_span_m[0]
                    && last.state == r.response.next_state
                {
                    last.span_m[1] = r.region.material_span_m[1];
                    continue;
                }
            }
            parcels.push(MaterialSlipParcel {
                span_m: r.region.material_span_m,
                state: r.response.next_state.clone(),
            });
        }
        check()?;
        let mut contact = self.clone();
        contact.parcels = parcels;
        contact.position = end;
        contact.time_s = target_s;
        Ok(MovingPartialSlipCandidate {
            contact,
            thermal,
            report: MovingPartialSlipReport {
                regions: staged.regions,
                heat,
                source: staged.source,
                nonlinear: result.report,
                trials: trials.get(),
                force_on_a_n: staged.force,
                moment_on_a_about_origin_nm: staged.moment,
                relative_slip_m: slip,
                mechanical_work_j: staged.work,
                storage_before_j: before,
                storage_after_j: staged.storage,
                elastic_outflow_j: outflow,
                numerical_dissipation_j: staged.numerical,
                total_energy_residual_j: total,
            },
        })
    }
}
