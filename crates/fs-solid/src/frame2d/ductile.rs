//! Local ductile initiation on the original distributed-plasticity frame.
//!
//! Every accepted native load, displacement or Newmark increment contributes
//! its original accumulated equivalent plastic strain and physical axial
//! stress. Section resultants, total strain and stress relative to a backstress
//! are not substitutes. A failed block publishes neither frame motion nor any
//! initiation history; refine that block from the unchanged accepted state.
//!
//! Gauges name original member/section/fiber points. Unobserved points are not
//! qualified by these gauges. Small axial strains, fixed actual temperatures,
//! rate-independent J2/Voce/tabulated/Chaboche mechanics and endpoint-held
//! initiation quadrature define the scope. No alpha is invented for empirical
//! Menegotto--Pinto steel or Mander concrete. This gate supplies no progressive
//! softening, local buckling, necking, fracture energy or separated members.

use super::{DynamicHistory, Frame2d, Rayleigh, SectionModel, StaticStep};
use crate::SolidError;
use crate::ductile::DuctileMaterialBinding;
use crate::fiber::chaboche::{ChabocheFiberLaw, ChabocheFiberState};
use crate::fiber::material::FiberMaterial;
use crate::fiber::{Fiber, FiberLaw, FiberState};
use fs_exec::Cx;
use fs_material::ductile::{DuctileError, DuctileState, DuctileStep};
use fs_material::plastic::chaboche::axial::{UniaxialChaboche, UniaxialChabocheState};
use fs_material::state_point::ResolvedMaterialStatePoint;
use fs_qty::semantic::{QuantityKind, ValueForm};
use fs_qty::{QuantitySpec, Temperature};
use std::collections::BTreeSet;

/// Immutable calibration at one original fiber integration point.
#[derive(Debug)]
pub struct FrameDuctileGauge<'a, L> {
    /// Original member index, in frame insertion order.
    pub element: usize,
    /// Original Gauss--Lobatto section index within that member.
    pub section: usize,
    /// Original material-fiber index within that section.
    pub fiber: usize,
    /// Immutable local plastic-strain initiation calibration.
    pub law: &'a L,
    /// Actual fixed specimen temperature, in kelvin.
    pub temperature_k: f64,
    /// Original material binding; required when the initiation law is sourced.
    /// Its native coefficients must match the selected fiber exactly.
    pub material_source: Option<&'a FiberMaterial>,
}
impl<L> Copy for FrameDuctileGauge<'_, L> {}
impl<L> Clone for FrameDuctileGauge<'_, L> {
    fn clone(&self) -> Self {
        *self
    }
}

/// Physical observation domain and explicit retained-state/work limits.
#[derive(Clone, Copy, Debug)]
pub struct FrameDuctileConfig {
    /// Maximum number of retained original fiber points.
    pub maximum_gauges: usize,
    /// Maximum number of accepted physical increments across the entire path.
    pub maximum_accepted_steps: u64,
    /// Positive axial-strain limit no greater than 0.05.
    pub maximum_absolute_axial_strain: f64,
}

/// Read-only quantities at an accepted original fiber point.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FrameDuctileSample {
    /// Total axial strain in the native small-strain fiber model.
    pub axial_strain: f64,
    /// Signed physical axial Cauchy stress [Pa], before any backstress shift.
    pub stress_pa: f64,
    /// Irreversible dimensionless plastic path length, including reverse flow.
    pub equivalent_plastic_strain: f64,
}

/// A refusal publishes none of the proposed mechanical or initiation state.
#[derive(Debug)]
pub enum FrameDuctileError {
    /// Invalid point, physical clock, history or model-domain request.
    Invalid(&'static str),
    /// Declared point or physical-increment budget would be exceeded.
    Budget,
    /// Refusal from the unchanged native equilibrium or momentum owner.
    Mechanics(SolidError),
    /// Original material initiation law refused a local physical condition.
    Point {
        /// Index in the path's retained gauge order.
        gauge: usize,
        /// Original material-law error.
        source: DuctileError,
    },
    /// Missing or incompatible original native material binding.
    SourceMismatch {
        /// Index in the path's retained gauge order.
        gauge: usize,
        /// Material identity, coefficient or physical-condition mismatch.
        what: &'static str,
    },
    /// Initiation with the smallest endpoint-held plastic-increment fraction.
    /// Fractions at different fibers do not establish physical event ordering.
    /// The entire proposed block is refused, including its earlier increments.
    Initiated {
        /// Index in the path's retained gauge order; ties use the first gauge.
        gauge: usize,
        /// Zero-based native increment within the refused block.
        step: usize,
        /// Physical endpoints bracketing that increment [s], not an exact event time.
        interval_s: [f64; 2],
        /// Original local initiation proposal, including any utilization overshoot.
        candidate: DuctileStep,
    },
    /// Cancellation observed before publishing a complete block.
    Cancelled,
}
impl core::fmt::Display for FrameDuctileError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native frame ductile-initiation gate refused: {self:?}")
    }
}
impl std::error::Error for FrameDuctileError {}
impl From<SolidError> for FrameDuctileError {
    fn from(value: SolidError) -> Self {
        if value == SolidError::Cancelled {
            Self::Cancelled
        } else {
            Self::Mechanics(value)
        }
    }
}

/// Original frame, accepted local histories and one physical clock.
/// Cloning retains the complete native constitutive and inertial checkpoint.
#[derive(Debug)]
pub struct FrameDuctilePath<'a, L: DuctileMaterialBinding> {
    frame: Frame2d,
    gauges: Vec<FrameDuctileGauge<'a, L>>,
    points: Vec<DuctileState>,
    samples: Vec<FrameDuctileSample>,
    config: FrameDuctileConfig,
    time_s: f64,
    accepted_steps: u64,
}
impl<L: DuctileMaterialBinding> Clone for FrameDuctilePath<'_, L> {
    fn clone(&self) -> Self {
        Self {
            frame: self.frame.clone(),
            gauges: self.gauges.clone(),
            points: self.points.clone(),
            samples: self.samples.clone(),
            config: self.config,
            time_s: self.time_s,
            accepted_steps: self.accepted_steps,
        }
    }
}

impl<'a, L: DuctileMaterialBinding> FrameDuctilePath<'a, L> {
    /// Attach virgin initiation only to fibers with no prior plastic flow.
    /// An elastic preload is permitted; prior inelastic damage cannot be
    /// reconstructed from an endpoint. Duplicate point bindings refuse.
    pub fn new(
        frame: Frame2d,
        gauges: Vec<FrameDuctileGauge<'a, L>>,
        initial_time_s: f64,
        config: FrameDuctileConfig,
        cx: &Cx<'_>,
    ) -> Result<Self, FrameDuctileError> {
        checkpoint(cx)?;
        frame.require_untimed_materials()?;
        if gauges.is_empty() || gauges.len() > config.maximum_gauges {
            return Err(FrameDuctileError::Budget);
        }
        if !initial_time_s.is_finite()
            || !config.maximum_absolute_axial_strain.is_finite()
            || config.maximum_absolute_axial_strain <= 0.0
            || config.maximum_absolute_axial_strain > 0.05
        {
            return Err(FrameDuctileError::Invalid(
                "finite initial time and a positive axial-strain limit at most 0.05 required",
            ));
        }
        let mut locations = BTreeSet::new();
        let mut samples = Vec::with_capacity(gauges.len());
        for (index, gauge) in gauges.iter().enumerate() {
            checkpoint(cx)?;
            if !locations.insert((gauge.element, gauge.section, gauge.fiber)) {
                return Err(FrameDuctileError::Invalid("duplicate original fiber gauge"));
            }
            let (fiber, _, _) = point(&frame, gauge)?;
            bind(index, gauge, &fiber.law)?;
            let sample = sample(&frame, gauge, config.maximum_absolute_axial_strain)?;
            if sample.equivalent_plastic_strain != 0.0 {
                return Err(FrameDuctileError::Invalid(
                    "virgin ductility cannot reconstruct previously accumulated plastic history",
                ));
            }
            // Even an elastic starting point must admit its actual temperature.
            DuctileState::virgin()
                .trial(
                    gauge.law,
                    0.0,
                    1.0,
                    gauge.temperature_k,
                    &[sample.stress_pa, 0.0, 0.0, 0.0, 0.0, 0.0],
                )
                .map_err(|source| FrameDuctileError::Point {
                    gauge: index,
                    source,
                })?;
            samples.push(sample);
        }
        checkpoint(cx)?;
        Ok(Self {
            points: vec![DuctileState::virgin(); gauges.len()],
            frame,
            gauges,
            samples,
            config,
            time_s: initial_time_s,
            accepted_steps: 0,
        })
    }

    /// Complete immutable accepted native mechanics.
    pub const fn frame(&self) -> &Frame2d {
        &self.frame
    }
    /// Original calibration and physical-point bindings, in retained order.
    pub fn gauges(&self) -> &[FrameDuctileGauge<'a, L>] {
        &self.gauges
    }
    /// Accepted initiation states in the same order as [`Self::gauges`].
    pub fn points(&self) -> &[DuctileState] {
        &self.points
    }
    /// Actual accepted native fiber observations, without a second return map.
    pub fn samples(&self) -> &[FrameDuctileSample] {
        &self.samples
    }
    /// Accepted absolute physical clock [s]. Native history records retain
    /// their original block-local time coordinates.
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    /// Physical increments accepted since this path was constructed.
    pub const fn accepted_steps(&self) -> u64 {
        self.accepted_steps
    }

    /// Native load control; as in `Frame2d::static_load`, each call ramps from
    /// zero to the supplied target. Each accepted increment spans the explicit
    /// `duration_per_step_s`; a load factor never supplies physical time.
    pub fn static_load(
        &mut self,
        load: &[f64],
        steps: usize,
        duration_per_step_s: f64,
        cx: &Cx<'_>,
    ) -> Result<Vec<StaticStep>, FrameDuctileError> {
        self.advance(steps, duration_per_step_s, cx, |frame, accepted| {
            frame.static_load_observed(load, steps, &mut || checkpoint(cx), accepted)
        })
    }

    /// Observe every original converged displacement target, including reversals.
    /// Each target has the same explicit physical duration; use successive calls
    /// for unequally spaced samples.
    pub fn displacement_control(
        &mut self,
        pattern: &[f64],
        control: usize,
        targets: &[f64],
        duration_per_step_s: f64,
        cx: &Cx<'_>,
    ) -> Result<Vec<StaticStep>, FrameDuctileError> {
        self.advance(targets.len(), duration_per_step_s, cx, |frame, accepted| {
            frame.displacement_control_observed(
                pattern,
                control,
                targets,
                &mut || checkpoint(cx),
                accepted,
            )
        })
    }

    /// Native Newmark inertia and initiation use exactly the same `dt`.
    /// The first ground sample is the existing endpoint and consumes no time
    /// or ductility. Adjacent blocks repeat that sample, as the native owner
    /// requires; its original Rayleigh reference and physical losses survive.
    #[allow(clippy::too_many_arguments)]
    pub fn newmark_with_load(
        &mut self,
        ground: &[f64],
        dt: f64,
        dir: (f64, f64),
        damping: Rayleigh,
        sustained_load: &[f64],
        cx: &Cx<'_>,
    ) -> Result<DynamicHistory, FrameDuctileError> {
        self.advance(ground.len().saturating_sub(1), dt, cx, |frame, accepted| {
            frame.newmark_with_load_observed(
                ground,
                dt,
                dir,
                damping,
                sustained_load,
                &mut || checkpoint(cx),
                accepted,
            )
        })
    }

    fn advance<T>(
        &mut self,
        steps: usize,
        duration_s: f64,
        cx: &Cx<'_>,
        run: impl FnOnce(
            &mut Frame2d,
            &mut dyn FnMut(&Frame2d) -> Result<(), SolidError>,
        ) -> Result<T, SolidError>,
    ) -> Result<T, FrameDuctileError> {
        checkpoint(cx)?;
        if !duration_s.is_finite() || duration_s <= 0.0 {
            return Err(FrameDuctileError::Invalid(
                "positive finite physical duration required",
            ));
        }
        let additional = u64::try_from(steps).map_err(|_| FrameDuctileError::Budget)?;
        let accepted_steps = self
            .accepted_steps
            .checked_add(additional)
            .filter(|&n| n <= self.config.maximum_accepted_steps)
            .ok_or(FrameDuctileError::Budget)?;
        let mut frame = self.frame.clone();
        let mut points = self.points.clone();
        let mut samples = self.samples.clone();
        let mut time_s = self.time_s;
        let mut observed = 0;
        let mut failure = None;
        let result = run(&mut frame, &mut |accepted| {
            let result = (|| {
                checkpoint(cx)?;
                let next_time_s = time_s + duration_s;
                if !next_time_s.is_finite() || next_time_s <= time_s {
                    return Err(FrameDuctileError::Invalid(
                        "physical clock increment is not representable",
                    ));
                }
                let mut initiation: Option<(usize, DuctileStep)> = None;
                for (index, gauge) in self.gauges.iter().enumerate() {
                    if index % 256 == 0 {
                        checkpoint(cx)?;
                    }
                    let endpoint =
                        sample(accepted, gauge, self.config.maximum_absolute_axial_strain)?;
                    let candidate = points[index]
                        .trial(
                            gauge.law,
                            endpoint.equivalent_plastic_strain,
                            duration_s,
                            gauge.temperature_k,
                            &[endpoint.stress_pa, 0.0, 0.0, 0.0, 0.0, 0.0],
                        )
                        .map_err(|source| FrameDuctileError::Point {
                            gauge: index,
                            source,
                        })?;
                    if candidate.state.initiated()
                        && initiation.as_ref().is_none_or(|(_, prior)| {
                            candidate.initiation_fraction.unwrap_or(0.0)
                                < prior.initiation_fraction.unwrap_or(0.0)
                        })
                    {
                        initiation = Some((index, candidate));
                    }
                    points[index] = candidate.state;
                    samples[index] = endpoint;
                }
                if let Some((gauge, candidate)) = initiation {
                    return Err(FrameDuctileError::Initiated {
                        gauge,
                        step: observed,
                        interval_s: [time_s, next_time_s],
                        candidate,
                    });
                }
                time_s = next_time_s;
                observed += 1;
                Ok(())
            })();
            result.map_err(|error| {
                failure = Some(error);
                SolidError::MaterialRefused {
                    what: "native frame ductility observation refused".into(),
                }
            })
        });
        let result = result.map_err(|error| failure.unwrap_or_else(|| error.into()))?;
        checkpoint(cx)?;
        if observed != steps {
            return Err(FrameDuctileError::Invalid(
                "native accepted-step count differs from ductile observation",
            ));
        }
        self.frame = frame;
        self.points = points;
        self.samples = samples;
        self.time_s = time_s;
        self.accepted_steps = accepted_steps;
        Ok(result)
    }
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), SolidError> {
    cx.checkpoint().map_err(|_| SolidError::Cancelled)
}

fn point<'f, L>(
    frame: &'f Frame2d,
    gauge: &FrameDuctileGauge<'_, L>,
) -> Result<(&'f Fiber, &'f FiberState, f64), FrameDuctileError> {
    let invalid = || {
        FrameDuctileError::Invalid("gauge must name an original member, fiber section and fiber")
    };
    let element = frame.elements.get(gauge.element).ok_or_else(invalid)?;
    let SectionModel::Fiber(section) = element.sections.get(gauge.section).ok_or_else(invalid)?
    else {
        return Err(invalid());
    };
    let fiber = section.fibers.get(gauge.fiber).ok_or_else(invalid)?;
    let state = section.states.get(gauge.fiber).ok_or_else(invalid)?;
    let &(eps0, curvature) = element.commit.d.get(gauge.section).ok_or_else(invalid)?;
    Ok((fiber, state, eps0 - fiber.y * curvature))
}

fn axial_chaboche<const N: usize>(
    law: &UniaxialChaboche<N>,
    state: &UniaxialChabocheState<N>,
    strain: f64,
) -> (f64, f64) {
    let youngs = law.backbone().backbone().backbone().elastic.youngs;
    (
        youngs * (strain - state.plastic_strain()),
        state.equivalent_plastic_strain(),
    )
}

fn sample<L>(
    frame: &Frame2d,
    gauge: &FrameDuctileGauge<'_, L>,
    maximum_strain: f64,
) -> Result<FrameDuctileSample, FrameDuctileError> {
    let (fiber, state, strain) = point(frame, gauge)?;
    if !strain.is_finite() || strain.abs() > maximum_strain {
        return Err(FrameDuctileError::Invalid(
            "accepted fiber strain exceeds the declared small-strain domain",
        ));
    }
    let (stress, alpha) = match (&fiber.law, state) {
        (FiberLaw::J2(law), FiberState::J2(state)) => (
            law.backbone().elastic.youngs * (strain - state.plastic_strain()),
            state.equivalent_plastic_strain(),
        ),
        (FiberLaw::Voce(law), FiberState::Voce(state)) => (
            law.backbone().backbone().elastic.youngs * (strain - state.plastic_strain()),
            state.equivalent_plastic_strain(),
        ),
        (FiberLaw::Tabulated(law), FiberState::Tabulated(state)) => (
            law.backbone().backbone().elastic.youngs * (strain - state.plastic_strain()),
            state.equivalent_plastic_strain(),
        ),
        (FiberLaw::Chaboche(law), FiberState::Chaboche(state)) => {
            axial_chaboche(law, state, strain)
        }
        (FiberLaw::ChabocheBranches(law), FiberState::ChabocheBranches(state)) => {
            match (law, state) {
                (ChabocheFiberLaw::Two(l), ChabocheFiberState::Two(s)) => {
                    axial_chaboche(l, s, strain)
                }
                (ChabocheFiberLaw::Three(l), ChabocheFiberState::Three(s)) => {
                    axial_chaboche(l, s, strain)
                }
                (ChabocheFiberLaw::Four(l), ChabocheFiberState::Four(s)) => {
                    axial_chaboche(l, s, strain)
                }
                (ChabocheFiberLaw::Five(l), ChabocheFiberState::Five(s)) => {
                    axial_chaboche(l, s, strain)
                }
                (ChabocheFiberLaw::Six(l), ChabocheFiberState::Six(s)) => {
                    axial_chaboche(l, s, strain)
                }
                (ChabocheFiberLaw::Seven(l), ChabocheFiberState::Seven(s)) => {
                    axial_chaboche(l, s, strain)
                }
                (ChabocheFiberLaw::Eight(l), ChabocheFiberState::Eight(s)) => {
                    axial_chaboche(l, s, strain)
                }
                _ => {
                    return Err(FrameDuctileError::Invalid(
                        "native Chaboche law and history branch counts differ",
                    ));
                }
            }
        }
        _ => {
            return Err(FrameDuctileError::Invalid(
                "ductility requires native rate-independent J2/Voce/tabulated/Chaboche plastic history",
            ));
        }
    };
    if !stress.is_finite() || !alpha.is_finite() || alpha < 0.0 {
        return Err(FrameDuctileError::Invalid(
            "accepted fiber stress or plastic history is not representable",
        ));
    }
    Ok(FrameDuctileSample {
        axial_strain: strain,
        stress_pa: stress,
        equivalent_plastic_strain: alpha,
    })
}

fn same_law(left: &FiberLaw, right: &FiberLaw) -> bool {
    match (left, right) {
        (FiberLaw::J2(a), FiberLaw::J2(b)) => a == b,
        (FiberLaw::Voce(a), FiberLaw::Voce(b)) => a == b,
        (FiberLaw::Tabulated(a), FiberLaw::Tabulated(b)) => a == b,
        (FiberLaw::Chaboche(a), FiberLaw::Chaboche(b)) => a == b,
        (FiberLaw::ChabocheBranches(a), FiberLaw::ChabocheBranches(b)) => a == b,
        _ => false,
    }
}

fn source_temperature_matches(material: &ResolvedMaterialStatePoint, temperature_k: f64) -> bool {
    let mut found = false;
    for (axis, quantity) in material.axis_quantities() {
        let absolute = *quantity == QuantitySpec::dimensional(Temperature::DIMS)
            || quantity.semantic_type().is_some_and(|semantic| {
                semantic.kind() == QuantityKind::AbsoluteTemperature
                    && semantic.form() == ValueForm::Static
            });
        if absolute {
            found = true;
            if material.state_coordinate(axis) != Some(temperature_k) {
                return false;
            }
        }
    }
    found
}

fn bind<L: DuctileMaterialBinding>(
    index: usize,
    gauge: &FrameDuctileGauge<'_, L>,
    native: &FiberLaw,
) -> Result<(), FrameDuctileError> {
    let mismatch = |what| FrameDuctileError::SourceMismatch { gauge: index, what };
    let source = match gauge.material_source {
        Some(binding) => {
            if !same_law(binding.law(), native) {
                return Err(mismatch(
                    "material binding differs from the original native fiber law",
                ));
            }
            binding.source()
        }
        None => None,
    };
    match source {
        Some(source) => {
            if !source_temperature_matches(source, gauge.temperature_k)
                || !gauge
                    .law
                    .matches_material_temperature(source, gauge.temperature_k)
            {
                return Err(mismatch(
                    "calibration differs from the original fiber card, state, fixed context or actual temperature",
                ));
            }
        }
        None if !gauge.law.permits_unbound_material() => {
            return Err(mismatch(
                "source-bound initiation requires the original source-bound fiber material",
            ));
        }
        None => {}
    }
    Ok(())
}
