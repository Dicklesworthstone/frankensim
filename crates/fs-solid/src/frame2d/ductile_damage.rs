//! Progressive axial ductility in the original force-based frame.
//!
//! Every prepared damage fiber consumes one explicit physical interval per
//! accepted native increment. Newton iterations remain disposable trials from
//! the same committed history. The original section stress and full algorithmic
//! tangent enter equilibrium, including negative tangents during softening.
//!
//! The receiving band is the original Gauss--Lobatto tributary length `w*L`.
//! Its work normalization is local; it does not prove mesh-independent crack
//! trajectories. This one-dimensional law attenuates tension and compression.
//! It does not create separated members, contact or a compression-closure law.
//! Singular sections refuse without a residual stiffness floor. Returned energy
//! covers progressive fibers only, not other members, inertia or Rayleigh loss.

use super::{DynamicHistory, Frame2d, Rayleigh, SectionModel, StaticStep};
use crate::SolidError;
use crate::fiber::{FiberLaw, FiberState};
use fs_exec::Cx;
use fs_exec::reduce::Compensated;
use fs_material::ductile::evolution::DuctileEvolutionError;
use fs_material::ductile::uniaxial::UniaxialDuctileJ2State;

/// Retained-point, physical-step and small-strain limits.
#[derive(Clone, Copy, Debug)]
pub struct FrameDuctileDamageConfig {
    /// Maximum number of progressive fibers across all integration sections.
    pub maximum_points: usize,
    /// Maximum number of accepted physical increments across this path.
    pub maximum_accepted_steps: u64,
    /// Positive bound on absolute native axial strain, at most 0.05.
    pub maximum_absolute_axial_strain: f64,
}

/// Original receiving quadrature point and its material observation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FrameDuctileDamagePoint {
    /// Member index in original insertion order.
    pub element: usize,
    /// Gauss--Lobatto section index in that member.
    pub section: usize,
    /// Original fiber index in that section.
    pub fiber: usize,
    /// Original quadrature tributary length and required material band [m].
    pub characteristic_length_m: f64,
    /// Original tributary fiber volume [m3].
    pub reference_volume_m3: f64,
    /// Accepted native total axial strain.
    pub axial_strain: f64,
    /// Accepted degraded physical axial stress [Pa].
    pub stress_pa: f64,
    /// Accepted native equivalent plastic path length, including reversals.
    pub equivalent_plastic_strain: f64,
    /// Accepted irreversible attenuation, between zero and one.
    pub damage: f64,
}

/// Local work of the progressive fibers over one accepted block, in joules.
/// No value is automatically deposited as heat or added to the native dynamic
/// history's work account. In particular, numerical loss is not physical heat.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct FrameDuctileDamageEnergy {
    /// Stored material energy at the beginning of the block.
    pub initial_stored_j: f64,
    /// Stored material energy at the end of the block.
    pub final_stored_j: f64,
    /// Sum of original endpoint stress times total strain increment and volume.
    pub mechanical_work_j: f64,
    /// Degraded native plastic dissipation.
    pub plastic_dissipation_j: f64,
    /// Released stored energy due to irreversible damage.
    pub damage_release_j: f64,
    /// Backward-Euler material integration loss, accounted separately.
    pub numerical_dissipation_j: f64,
    /// Work minus storage change and the three independently retained losses.
    pub residual_j: f64,
}

/// Original mechanical output plus the progressive fibers' separate account.
#[derive(Clone, Debug)]
pub struct FrameDuctileDamageBlock<T> {
    /// Output of the unchanged native load, displacement or momentum solve.
    pub mechanics: T,
    /// Material work over exactly the accepted physical increments.
    pub energy: FrameDuctileDamageEnergy,
    /// Absolute physical beginning and end of the block [s].
    pub interval_s: [f64; 2],
}

/// Refusal leaves the complete input frame, damage, motion and clock unchanged.
#[derive(Debug)]
pub enum FrameDuctileDamageError {
    /// Invalid quadrature, physical domain or unsupported mixed timed material.
    Invalid(&'static str),
    /// Point or accepted-physical-step limit would be exceeded.
    Budget,
    /// Native equilibrium or momentum failure, including singular sections.
    Mechanics(SolidError),
    /// Original material refused interval preparation.
    Material {
        /// Point index in the retained deterministic order.
        point: usize,
        /// Original constitutive refusal.
        source: DuctileEvolutionError,
    },
    /// Cancellation before publication of a complete physical block.
    Cancelled,
}
impl core::fmt::Display for FrameDuctileDamageError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "progressive frame ductility refused: {self:?}")
    }
}
impl std::error::Error for FrameDuctileDamageError {}
impl From<SolidError> for FrameDuctileDamageError {
    fn from(error: SolidError) -> Self {
        if error == SolidError::Cancelled {
            Self::Cancelled
        } else {
            Self::Mechanics(error)
        }
    }
}

/// Complete native frame with progressive fibers and an explicit physical clock.
/// Cloning preserves the full accepted constitutive and inertial checkpoint.
#[derive(Clone, Debug)]
pub struct FrameDuctileDamagePath {
    frame: Frame2d,
    points: Vec<FrameDuctileDamagePoint>,
    config: FrameDuctileDamageConfig,
    time_s: f64,
    accepted_steps: u64,
}

impl FrameDuctileDamagePath {
    /// Admit prepared progressive fibers at their original receiving lengths.
    /// Other untimed fibers and elastic sections are allowed; other timed
    /// families require their own coupled physical owner and are refused.
    /// Each progressive fiber retains its prepared, fixed physical temperature.
    pub fn new(
        frame: Frame2d,
        initial_time_s: f64,
        config: FrameDuctileDamageConfig,
        cx: &Cx<'_>,
    ) -> Result<Self, FrameDuctileDamageError> {
        checkpoint(cx)?;
        if !initial_time_s.is_finite()
            || !config.maximum_absolute_axial_strain.is_finite()
            || config.maximum_absolute_axial_strain <= 0.0
            || config.maximum_absolute_axial_strain > 0.05
        {
            return Err(FrameDuctileDamageError::Invalid(
                "finite clock and positive axial-strain bound at most 0.05 required",
            ));
        }
        let mut points = Vec::new();
        for (element, member) in frame.elements.iter().enumerate() {
            checkpoint(cx)?;
            for (section, model) in member.sections.iter().enumerate() {
                let SectionModel::Fiber(model) = model else {
                    continue;
                };
                if model.fibers.len() != model.states.len() {
                    return Err(FrameDuctileDamageError::Invalid(
                        "original fiber and state counts differ",
                    ));
                }
                for (fiber, (definition, state)) in
                    model.fibers.iter().zip(&model.states).enumerate()
                {
                    if fiber % 256 == 0 {
                        checkpoint(cx)?;
                    }
                    let FiberLaw::DuctileJ2(law) = &definition.law else {
                        if definition.law.requires_physical_duration() {
                            return Err(FrameDuctileDamageError::Invalid(
                                "mixed timed material requires a coupled physical owner",
                            ));
                        }
                        continue;
                    };
                    if points.len() >= config.maximum_points {
                        return Err(FrameDuctileDamageError::Budget);
                    }
                    let FiberState::DuctileJ2(state) = state else {
                        return Err(FrameDuctileDamageError::Invalid(
                            "progressive fiber needs its original damage history",
                        ));
                    };
                    let length = member.w[section] * member.l0;
                    let volume = length * definition.area;
                    let selected = law.settings().characteristic_length_m;
                    if !length.is_finite()
                        || length <= 0.0
                        || !volume.is_finite()
                        || volume <= 0.0
                        || (selected - length).abs() > 32.0 * f64::EPSILON * selected.max(length)
                        || !law.owns_state(state)
                    {
                        return Err(FrameDuctileDamageError::Invalid(
                            "damage model must retain original history and the w*L receiving band",
                        ));
                    }
                    let mut point = FrameDuctileDamagePoint {
                        element,
                        section,
                        fiber,
                        characteristic_length_m: length,
                        reference_volume_m3: volume,
                        axial_strain: 0.0,
                        stress_pa: 0.0,
                        equivalent_plastic_strain: 0.0,
                        damage: 0.0,
                    };
                    observe(&frame, &mut point, config.maximum_absolute_axial_strain)?;
                    points.push(point);
                }
            }
        }
        if points.is_empty() {
            return Err(FrameDuctileDamageError::Invalid(
                "at least one original progressive fiber required",
            ));
        }
        checkpoint(cx)?;
        Ok(Self {
            frame,
            points,
            config,
            time_s: initial_time_s,
            accepted_steps: 0,
        })
    }

    /// Read-only complete accepted native mechanics.
    pub const fn frame(&self) -> &Frame2d {
        &self.frame
    }
    /// Original points and their accepted physical observations.
    pub fn points(&self) -> &[FrameDuctileDamagePoint] {
        &self.points
    }
    /// Absolute accepted physical time [s]. Native histories retain block-local times.
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    /// Number of accepted physical increments since construction.
    pub const fn accepted_steps(&self) -> u64 {
        self.accepted_steps
    }

    /// Original zero-to-target load ramp with one explicit duration per increment.
    pub fn static_load(
        &mut self,
        load: &[f64],
        steps: usize,
        duration_per_step_s: f64,
        cx: &Cx<'_>,
    ) -> Result<FrameDuctileDamageBlock<Vec<StaticStep>>, FrameDuctileDamageError> {
        self.advance(steps, duration_per_step_s, cx, |frame, accepted| {
            frame.static_load_observed(load, steps, &mut || checkpoint(cx), accepted)
        })
    }

    /// Follow original bordered equilibrium through displacement targets,
    /// including regular descending branches where the native solver converges.
    pub fn displacement_control(
        &mut self,
        pattern: &[f64],
        control: usize,
        targets: &[f64],
        duration_per_step_s: f64,
        cx: &Cx<'_>,
    ) -> Result<FrameDuctileDamageBlock<Vec<StaticStep>>, FrameDuctileDamageError> {
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

    /// Original Newmark motion with the same `dt` in the damage criterion.
    /// The first ground sample represents the already accepted endpoint.
    #[allow(clippy::too_many_arguments)]
    pub fn newmark_with_load(
        &mut self,
        ground: &[f64],
        dt: f64,
        dir: (f64, f64),
        damping: Rayleigh,
        sustained_load: &[f64],
        cx: &Cx<'_>,
    ) -> Result<FrameDuctileDamageBlock<DynamicHistory>, FrameDuctileDamageError> {
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
    ) -> Result<FrameDuctileDamageBlock<T>, FrameDuctileDamageError> {
        checkpoint(cx)?;
        if !duration_s.is_finite() || duration_s <= 0.0 {
            return Err(FrameDuctileDamageError::Invalid(
                "positive finite physical duration required",
            ));
        }
        let additional = u64::try_from(steps).map_err(|_| FrameDuctileDamageError::Budget)?;
        let count = self
            .accepted_steps
            .checked_add(additional)
            .filter(|&n| n <= self.config.maximum_accepted_steps)
            .ok_or(FrameDuctileDamageError::Budget)?;
        let mut frame = self.frame.clone();
        for (index, point) in self.points.iter().enumerate() {
            if index % 256 == 0 {
                checkpoint(cx)?;
            }
            let SectionModel::Fiber(section) =
                &mut frame.elements[point.element].sections[point.section]
            else {
                unreachable!()
            };
            let FiberLaw::DuctileJ2(law) = &mut section.fibers[point.fiber].law else {
                unreachable!()
            };
            *law = law.with_duration(duration_s).map_err(|source| {
                FrameDuctileDamageError::Material {
                    point: index,
                    source,
                }
            })?;
        }
        let mut points = self.points.clone();
        let initial_stored_j = stored(&frame, &points, cx)?;
        let mut sums = [Compensated::zero(); 4];
        let mut time_s = self.time_s;
        let mut observed = 0usize;
        let mut failure = None;
        let result = run(&mut frame, &mut |accepted| {
            let result = (|| {
                checkpoint(cx)?;
                let next_time = time_s + duration_s;
                if !next_time.is_finite() || next_time <= time_s {
                    return Err(FrameDuctileDamageError::Invalid(
                        "physical clock increment is not representable",
                    ));
                }
                for (index, point) in points.iter_mut().enumerate() {
                    if index % 256 == 0 {
                        checkpoint(cx)?;
                    }
                    let state =
                        observe(accepted, point, self.config.maximum_absolute_axial_strain)?;
                    let energy = state.last_step_energy();
                    for (sum, density) in sums.iter_mut().zip([
                        energy.mechanical_work_j_m3,
                        energy.plastic_dissipation_j_m3,
                        energy.damage_release_j_m3,
                        energy.numerical_dissipation_j_m3,
                    ]) {
                        *sum = sum.accumulate(point.reference_volume_m3 * density);
                        if !sum.value().is_finite() {
                            return Err(FrameDuctileDamageError::Invalid(
                                "progressive fiber work is unrepresentable",
                            ));
                        }
                    }
                }
                time_s = next_time;
                observed += 1;
                Ok(())
            })();
            result.map_err(|error| {
                failure = Some(error);
                SolidError::MaterialRefused {
                    what: "progressive frame observation refused".into(),
                }
            })
        });
        let mechanics = result.map_err(|error| failure.unwrap_or_else(|| error.into()))?;
        checkpoint(cx)?;
        if observed != steps {
            return Err(FrameDuctileDamageError::Invalid(
                "native accepted-step count differs from physical observation",
            ));
        }
        let final_stored_j = stored(&frame, &points, cx)?;
        let [
            mechanical_work_j,
            plastic_dissipation_j,
            damage_release_j,
            numerical_dissipation_j,
        ] = sums.map(|s| s.value());
        let residual_j = mechanical_work_j
            - (final_stored_j - initial_stored_j)
            - plastic_dissipation_j
            - damage_release_j
            - numerical_dissipation_j;
        if !residual_j.is_finite() {
            return Err(FrameDuctileDamageError::Invalid(
                "material energy residual is unrepresentable",
            ));
        }
        let energy = FrameDuctileDamageEnergy {
            initial_stored_j,
            final_stored_j,
            mechanical_work_j,
            plastic_dissipation_j,
            damage_release_j,
            numerical_dissipation_j,
            residual_j,
        };
        let interval_s = [self.time_s, time_s];
        self.frame = frame;
        self.points = points;
        self.time_s = time_s;
        self.accepted_steps = count;
        Ok(FrameDuctileDamageBlock {
            mechanics,
            energy,
            interval_s,
        })
    }
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), SolidError> {
    cx.checkpoint().map_err(|_| SolidError::Cancelled)
}

fn state<'a>(frame: &'a Frame2d, point: &FrameDuctileDamagePoint) -> &'a UniaxialDuctileJ2State {
    let SectionModel::Fiber(section) = &frame.elements[point.element].sections[point.section]
    else {
        unreachable!()
    };
    let FiberState::DuctileJ2(state) = &section.states[point.fiber] else {
        unreachable!()
    };
    state
}

fn observe<'a>(
    frame: &'a Frame2d,
    point: &mut FrameDuctileDamagePoint,
    maximum_strain: f64,
) -> Result<&'a UniaxialDuctileJ2State, FrameDuctileDamageError> {
    let member = &frame.elements[point.element];
    let SectionModel::Fiber(section) = &member.sections[point.section] else {
        unreachable!()
    };
    let accepted = state(frame, point);
    let (eps0, curvature) = member.commit.d[point.section];
    let strain = eps0 - section.fibers[point.fiber].y * curvature;
    if !strain.is_finite()
        || strain.abs() > maximum_strain
        || strain != accepted.strain()
        || !accepted.is_valid()
    {
        return Err(FrameDuctileDamageError::Invalid(
            "original accepted fiber strain/history leaves the declared domain",
        ));
    }
    point.axial_strain = strain;
    point.stress_pa = accepted.stress_pa();
    point.equivalent_plastic_strain = accepted.equivalent_plastic_strain();
    point.damage = accepted.damage();
    Ok(accepted)
}

fn stored(
    frame: &Frame2d,
    points: &[FrameDuctileDamagePoint],
    cx: &Cx<'_>,
) -> Result<f64, FrameDuctileDamageError> {
    let mut sum = Compensated::zero();
    for (index, point) in points.iter().enumerate() {
        if index % 256 == 0 {
            checkpoint(cx)?;
        }
        sum = sum.accumulate(point.reference_volume_m3 * state(frame, point).stored_energy_j_m3());
    }
    let value = sum.value();
    if value.is_finite() && value >= 0.0 {
        Ok(value)
    } else {
        Err(FrameDuctileDamageError::Invalid(
            "progressive fiber storage is unrepresentable",
        ))
    }
}

/// Invert the actual indefinite section tangent without a stiffness floor.
/// Scaling each row avoids squaring a large elastic modulus in the determinant.
pub(super) fn section_flexibility(t: [[f64; 2]; 2]) -> Result<[[f64; 2]; 2], SolidError> {
    let refused = || {
        SolidError::MaterialRefused { what: "singular progressive section requires a different physical topology or continuation model".into() }
    };
    if t.iter().flatten().any(|x| !x.is_finite()) {
        return Err(refused());
    }
    let a = t[0][0].abs().max(t[0][1].abs());
    let b = t[1][0].abs().max(t[1][1].abs());
    if a == 0.0 || b == 0.0 {
        return Err(refused());
    }
    let n = [[t[0][0] / a, t[0][1] / a], [t[1][0] / b, t[1][1] / b]];
    let det = n[0][0] * n[1][1] - n[0][1] * n[1][0];
    if det == 0.0 || !det.is_finite() {
        return Err(refused());
    }
    let inverse = [
        [(n[1][1] / det) / a, (-n[0][1] / det) / b],
        [(-n[1][0] / det) / a, (n[0][0] / det) / b],
    ];
    if inverse.iter().flatten().any(|x| !x.is_finite()) {
        return Err(refused());
    }
    Ok(inverse)
}
