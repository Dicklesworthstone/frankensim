//! Supplied heat closes the original caloric and native expansion-work account.
//!
//! Planning owns the resolved endpoint sources; advancing borrows them into the
//! existing native volume step. No self-referential liquid or replacement EOS.

use super::*;
use crate::conjugate::{FlowFace, InternalPressureFace};
use fs_material::liquid_temperature::{LiquidCaloricEndpoint, LiquidCaloricError};

mod balance;

/// One authored heat interval and its original source/pressure accuracy limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreeSurfaceSlabHeatControls {
    pub volume: FreeSurfaceSlabVolumeControls,
    /// Explicit finite bracket for NET specific enthalpy [J/kg]. It selects
    /// the admissible source branch; no uniqueness outside it is claimed.
    pub net_specific_enthalpy_bracket_j_kg: [f64; 2],
    /// Maximum |represented delta_H + correction-pressure work - heat| [J].
    pub maximum_heat_defect_j: f64,
    /// Maximum actual native work minus the scalar response's work [J].
    pub maximum_pressure_work_defect_j: f64,
    /// Total original-caloric scalar evaluations, including both endpoints.
    pub maximum_source_evaluations: usize,
}

/// Original source, physical-domain or bounded heat-balance refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum FreeSurfaceSlabHeatError {
    Native(Error),
    Caloric(LiquidCaloricError),
    Invalid(&'static str),
    /// A missing sign bracket, exhausted work, or stalled floating-point
    /// interval never publishes a best-effort thermal or liquid state.
    Balance {
        evaluations: usize,
        residual_j: f64,
    },
}

impl core::fmt::Display for FreeSurfaceSlabHeatError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Native(e) => write!(f, "thermal slab heat: {e}"),
            Self::Caloric(e) => write!(f, "thermal slab caloric endpoint: {e}"),
            Self::Invalid(e) => write!(f, "thermal slab heat: {e}"),
            Self::Balance {
                evaluations,
                residual_j,
            } => write!(
                f,
                "thermal slab heat did not close in {evaluations} source evaluations (residual {residual_j} J)"
            ),
        }
    }
}
impl std::error::Error for FreeSurfaceSlabHeatError {}
impl From<Error> for FreeSurfaceSlabHeatError {
    fn from(error: Error) -> Self {
        Self::Native(error)
    }
}
impl From<LiquidCaloricError> for FreeSurfaceSlabHeatError {
    fn from(error: LiquidCaloricError) -> Self {
        Self::Caloric(error)
    }
}

/// Heat and pressure work are distinct physical channels, evaluated on the
/// original mass and chart. Numerical pressure/remap losses stay in the
/// accompanying volume report and are never deposited as thermal energy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreeSurfaceSlabHeatReport {
    pub supplied_heat_j: f64,
    pub net_specific_enthalpy_increment_j_kg: f64,
    pub represented_enthalpy_increment_j: f64,
    pub pressure_correction_work_j: f64,
    pub heat_balance_defect_j: f64,
    pub scalar_heat_balance_defect_j: f64,
    pub pressure_response_defect_j: f64,
    /// Native face-dual inertia for unit outward interface speed [kg].
    pub surface_speed_inertia_kg: f64,
    pub source_evaluations: usize,
}

/// A source-bound heat endpoint. Its borrow prevents the original liquid from
/// changing while this plan or any returned endpoint liquid is alive.
pub struct FreeSurfaceSlabHeatPlan<'a, 's, 'm> {
    original: &'a MaterialFreeSurfaceEnthalpy3<'s, 'm>,
    endpoint: LiquidCaloricEndpoint<'m>,
    anchored_face: Face3,
    controls: FreeSurfaceSlabHeatControls,
    supplied_heat_j: f64,
    solution: balance::Solution,
    surface_speed_inertia_kg: f64,
}

/// Accepted native moving liquid plus its independently checked heat closure.
#[derive(Clone)]
pub struct FreeSurfaceSlabHeatStep<'s, 'm> {
    pub volume: FreeSurfaceSlabVolumeStep<'s, 'm>,
    pub report: FreeSurfaceSlabHeatReport,
}

impl<'a, 's, 'm> FreeSurfaceSlabHeatPlan<'a, 's, 'm> {
    /// The owned endpoint retains the original caloric datum, source claims,
    /// external conditions and any density-derived expansion metadata.
    #[must_use]
    pub const fn endpoint(&self) -> &LiquidCaloricEndpoint<'m> {
        &self.endpoint
    }

    /// Planned original-mass enthalpy change [J], before native readback.
    #[must_use]
    pub const fn net_enthalpy_increment_j(&self) -> f64 {
        self.solution.enthalpy_increment_j
    }

    /// Planned mechanical correction work [J]. The accepted native step
    /// independently measures this amount before it can be returned.
    #[must_use]
    pub const fn required_pressure_work_j(&self) -> f64 {
        self.solution.work_j
    }

    /// Execute the ONE actual moving-interface pressure/remap step, then check
    /// its measured work and represented enthalpy against the supplied heat.
    /// All original native source, speed, mass, momentum, geometry and energy
    /// gates remain active. No rejected proposal changes the original owner.
    pub fn advance<'next>(
        &'next self,
        gate: &CancelGate,
        cx: &Cx<'_>,
    ) -> Result<FreeSurfaceSlabHeatStep<'next, 'next>, FreeSurfaceSlabHeatError>
    where
        'm: 'next,
    {
        checkpoint(gate, cx)?;
        let volume = self.original.advance_uniform_slab_volume(
            self.anchored_face,
            self.solution.increment_j_kg,
            self.endpoint.source(),
            self.endpoint.law(),
            self.controls.volume,
            gate,
            cx,
        )?;
        let actual = &volume.report;
        let enthalpy = Dd::from_f64(actual.requested_enthalpy_increment_j)
            + Dd::from_f64(actual.enthalpy_defect_j);
        let work = actual.required_mechanical_work_j;
        let work_defect = value(Dd::from_f64(work) - Dd::from_f64(self.solution.work_j));
        let heat_defect = value(enthalpy + Dd::from_f64(work) - Dd::from_f64(self.supplied_heat_j));
        if !work_defect.is_finite()
            || !heat_defect.is_finite()
            || work_defect.abs() > self.controls.maximum_pressure_work_defect_j
            || heat_defect.abs() > self.controls.maximum_heat_defect_j
        {
            return Err(FreeSurfaceSlabHeatError::Invalid(
                "actual native pressure work and caloric storage exceed the supplied-heat budget",
            ));
        }
        checkpoint(gate, cx)?;
        let report = FreeSurfaceSlabHeatReport {
            supplied_heat_j: self.supplied_heat_j,
            net_specific_enthalpy_increment_j_kg: self.solution.increment_j_kg,
            represented_enthalpy_increment_j: value(enthalpy),
            pressure_correction_work_j: work,
            heat_balance_defect_j: heat_defect,
            scalar_heat_balance_defect_j: self.solution.residual_j,
            pressure_response_defect_j: work_defect,
            surface_speed_inertia_kg: self.surface_speed_inertia_kg,
            source_evaluations: self.solution.evaluations,
        };
        Ok(FreeSurfaceSlabHeatStep { volume, report })
    }
}

impl<'s, 'm> MaterialFreeSurfaceEnthalpy3<'s, 'm> {
    /// Convert supplied heat into a conservative uniform thermal-volume step.
    ///
    /// The source-split balance is Q = M*delta_h + W, where W is the ORIGINAL
    /// pressure owner's correction work. The fixed-pressure caloric h already
    /// includes ambient p*dV; it is not charged a second time. The pressure
    /// correction remains mechanical and does not change chart pressure.
    /// Contraction can require positive work; retained expansion can return
    /// signed work. The original projection's numerical loss and explicit
    /// native/face/affine mapping discrepancies never become heat.
    ///
    /// A uniform slab has an affine volume-source velocity. On the original
    /// frozen cut metric, the native pressure work is B*u*(u-u0), where B is
    /// the native face-dual inertia and u is the outward free-surface speed.
    /// We read B from the EXISTING pressure owner's energy audit on a unit-speed
    /// affine field whose target outflow is already satisfied (zero pressure
    /// correction). No new Poisson solver or independently fitted mass is used.
    /// B is measured once. Bounded scalar trials then use the original chart
    /// and analytic specific-volume derivative; they do not rerun a 3-D solve
    /// or requery every material property. The accepted native step independently
    /// checks this response, the complete heat balance, and all original gates.
    ///
    /// The explicit bracket selects the source branch and must stay pure liquid.
    /// An unsupported phase, source, geometry, work or cancellation refuses.
    /// The same slab/source-split scope as advance_uniform_slab_volume applies:
    /// this does not add viscous evolution, arbitrary flow, or solid melting.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn prepare_uniform_slab_heat<'a>(
        &'a self,
        anchored_face: Face3,
        supplied_heat_j: f64,
        controls: FreeSurfaceSlabHeatControls,
        gate: &CancelGate,
        cx: &Cx<'_>,
    ) -> Result<FreeSurfaceSlabHeatPlan<'a, 's, 'm>, FreeSurfaceSlabHeatError> {
        checkpoint(gate, cx)?;
        self.admit_uniform_slab_sources(supplied_heat_j)?;
        for budget in [
            controls.maximum_heat_defect_j,
            controls.maximum_pressure_work_defect_j,
        ] {
            if !budget.is_finite() || budget < 0.0 {
                return Err(FreeSurfaceSlabHeatError::Invalid(
                    "invalid slab heat accuracy budget",
                ));
            }
        }
        let dt = self.flow.scaling.timestep_s();
        let dx = self.flow.scaling.cell_size_m();
        if controls.volume.volume.projection.duration_s != dt {
            return Err(FreeSurfaceSlabHeatError::Invalid(
                "slab heat must retain the native timestep",
            ));
        }
        let axis = anchored_face as usize / 2;
        let upper = anchored_face as usize % 2 == 1;
        let size = self.flow.lattice.domain_cells();
        let count = size[0]
            .checked_mul(size[1])
            .and_then(|n| n.checked_mul(size[2]))
            .filter(|n| *n > 0 && *n <= controls.volume.volume.max_cells)
            .ok_or(FreeSurfaceSlabHeatError::Invalid(
                "slab heat exceeds its original grid budget",
            ))?;
        if self
            .volume_slab
            .as_ref()
            .is_some_and(|v| v.anchored_face != anchored_face)
        {
            return Err(FreeSurfaceSlabHeatError::Invalid(
                "slab heat changed the anchored face",
            ));
        }
        self.flow
            .lattice
            .admit_volume_slab(
                axis,
                upper,
                self.volume_slab.as_ref().map(|v| &v.profile),
                gate,
            )
            .map_err(Error::from)?;
        let columns = size
            .iter()
            .enumerate()
            .filter(|(a, _)| *a != axis)
            .try_fold(1_usize, |n, (_, &width)| n.checked_mul(width))
            .ok_or(FreeSurfaceSlabHeatError::Invalid(
                "slab heat cross section overflows",
            ))?;
        let masses = self.flow.lattice.mass_snapshot();
        let mass_unit = self.flow.scaling.mass_unit_kg();
        let mut mass = Dd::ZERO;
        let mut fills = BTreeMap::new();
        let mut depth = 0;
        let mut occupied = 0;
        let mut h0 = None;
        for (&location, &native_mass) in &masses {
            checkpoint(gate, cx)?;
            if !native_mass.is_finite() || native_mass < 0.0 {
                return Err(FreeSurfaceSlabHeatError::Invalid(
                    "slab heat has invalid native inventory",
                ));
            }
            if native_mass == 0.0 {
                continue;
            }
            let coordinates = FreeSurface3::transport_coordinates(location).ok_or(
                FreeSurfaceSlabHeatError::Invalid("slab heat has unlocated mass"),
            )?;
            let h = *self
                .h
                .get(&location)
                .ok_or(FreeSurfaceSlabHeatError::Invalid(
                    "slab heat lost original enthalpy",
                ))?;
            if h0.is_some_and(|old: f64| old.to_bits() != h.to_bits()) {
                return Err(FreeSurfaceSlabHeatError::Invalid(
                    "slab heat needs uniform original enthalpy",
                ));
            }
            h0 = Some(h);
            let layer = if upper {
                size[axis] - 1 - coordinates[axis]
            } else {
                coordinates[axis]
            };
            let sample = self
                .flow
                .sample_at_si(coordinates[0], coordinates[1], coordinates[2])
                .map_err(Error::from)?
                .ok_or(FreeSurfaceSlabHeatError::Invalid(
                    "slab heat lost a native sample",
                ))?;
            let fill = sample.fill_fraction;
            if !fill.is_finite()
                || fill <= 0.0
                || fill > 1.0
                || fills.get(&layer).is_some_and(|&v| v != fill)
            {
                return Err(FreeSurfaceSlabHeatError::Invalid(
                    "slab heat needs planar positive layer fills",
                ));
            }
            fills.insert(layer, fill);
            depth = depth.max(layer + 1);
            occupied += 1;
            mass = mass + Dd::from_f64(native_mass) * Dd::from_f64(mass_unit);
        }
        if depth == 0
            || depth >= size[axis]
            || depth.checked_mul(columns) != Some(occupied)
            || fills
                .iter()
                .any(|(&layer, &fill)| layer + 1 < depth && fill != 1.0)
        {
            return Err(FreeSurfaceSlabHeatError::Invalid(
                "slab heat needs complete columns and a gas layer",
            ));
        }
        let h0 = h0.ok_or(FreeSurfaceSlabHeatError::Invalid("slab heat has no liquid"))?;
        if self
            .volume_slab
            .as_ref()
            .is_some_and(|v| v.specific_enthalpy_j_kg.to_bits() != h0.to_bits())
        {
            return Err(FreeSurfaceSlabHeatError::Invalid(
                "slab heat changed its retained enthalpy",
            ));
        }
        let initial = self
            .caloric
            .state_at_specific_enthalpy(h0)
            .map_err(Error::from)?;
        if initial.phase() != SolidLiquidPhase::Liquid {
            return Err(FreeSurfaceSlabHeatError::Invalid(
                "slab heat requires pure liquid",
            ));
        }
        let terminal = fills[&(depth - 1)];
        let domain = planar_domain(size, axis, upper, depth, terminal, dx, gate, cx)?;
        let height = if upper {
            domain.face_coord(axis, size[axis]) - domain.face_coord(axis, size[axis] - depth)
        } else {
            domain.face_coord(axis, depth)
        };
        let unit_rate = height.recip();
        if !height.is_finite() || height <= 0.0 || !unit_rate.is_finite() || unit_rate == 0.0 {
            return Err(FreeSurfaceSlabHeatError::Invalid(
                "slab heat has unrepresentable height",
            ));
        }
        let unit = affine_predictor(&domain, axis, upper, unit_rate, gate, cx)?;
        let mut density = vec![0.0; count];
        let mut target = vec![0.0; count];
        let mut openings = Vec::new();
        let free_face = Face3::ALL[2 * axis + usize::from(!upper)];
        for c in 0..count {
            checkpoint(gate, cx)?;
            if !domain.is_fluid(c) {
                continue;
            }
            let coordinates = domain.coords(c);
            let location = FreeSurface3::transport_location(coordinates).ok_or(
                FreeSurfaceSlabHeatError::Invalid("slab heat lost a cell address"),
            )?;
            let cell_mass = masses[&location] * mass_unit;
            density[c] = cell_mass / domain.volume(c);
            target[c] = unit.net_outflow(&domain, c);
            if domain
                .neighbor(c, free_face as usize)
                .is_some_and(|n| !domain.is_fluid(n))
            {
                openings.push(InternalPressureFace {
                    cell: c,
                    face: free_face,
                    pressure_pa: 0.0,
                });
            }
        }
        let measured = unit
            .project_density_volume_source_with_internal_pressure_faces(
                &domain,
                [FlowFace::Wall; 6],
                &openings,
                &density,
                &target,
                controls.volume.volume.projection,
                gate,
            )
            .map_err(Error::from)?;
        let inertia = 2.0 * measured.report.kinetic_energy_before_j;
        if !inertia.is_finite() || inertia <= 0.0 {
            return Err(FreeSurfaceSlabHeatError::Invalid(
                "slab heat native face inertia is unrepresentable",
            ));
        }
        let incoming_rate = self
            .volume_slab
            .as_ref()
            .map_or(0.0, |v| v.affine_rate_s_inv);
        let incoming_speed = slab::product_ratio([incoming_rate, height, 1.0, 1.0], [1.0, 1.0]);
        let reference_density = self.flow.scaling.material().density_kg_m3();
        let rho0 = initial.bulk_density_kg_m3();
        let solution = balance::solve(
            mass,
            supplied_heat_j,
            controls.net_specific_enthalpy_bracket_j_kg,
            controls.maximum_heat_defect_j,
            controls.maximum_source_evaluations,
            |increment| -> Result<balance::Response, Error> {
                checkpoint(gate, cx)?;
                let next = self.caloric.advance_specific_energy(initial, increment)?;
                if next.phase() != SolidLiquidPhase::Liquid {
                    return Err(Error::Invalid(
                        "slab heat bracket leaves the original liquid phase",
                    ));
                }
                let rho1 = next.bulk_density_kg_m3();
                // Same stable specific-volume difference as the native demand.
                let relative = (rho0 - rho1) / rho0.max(rho1);
                let speed = slab::product_ratio(
                    [reference_density, height, relative, 1.0],
                    [rho0.min(rho1), dt],
                );
                let slope = self
                    .caloric
                    .specific_volume_derivative_at_specific_enthalpy(
                        next.specific_enthalpy_j_kg(),
                    )?;
                let speed_slope =
                    slab::product_ratio([reference_density, height, slope, 1.0], [dt, 1.0]);
                if !speed.is_finite()
                    || !speed_slope.is_finite()
                    || !incoming_speed.is_finite()
                    || (rho0 != rho1 && speed == 0.0)
                    || (slope != 0.0 && speed_slope == 0.0)
                {
                    return Err(Error::Invalid(
                        "slab heat density/velocity response is unrepresentable",
                    ));
                }
                Ok(balance::pressure_response(
                    inertia,
                    speed,
                    incoming_speed,
                    speed_slope,
                ))
            },
        )
        .map_err(|failure| match failure {
            balance::Failure::Invalid(what) => FreeSurfaceSlabHeatError::Invalid(what),
            balance::Failure::Evaluation(error) => FreeSurfaceSlabHeatError::Native(error),
            balance::Failure::Balance {
                evaluations,
                residual_j,
            } => FreeSurfaceSlabHeatError::Balance {
                evaluations,
                residual_j,
            },
        })?;
        let endpoint = self.law.prepare_caloric_endpoint(
            self.caloric,
            initial,
            solution.increment_j_kg,
            controls.volume.volume.maximum_source_density_defect_kg_m3,
        )?;
        checkpoint(gate, cx)?;
        Ok(FreeSurfaceSlabHeatPlan {
            original: self,
            endpoint,
            anchored_face,
            controls,
            supplied_heat_j,
            solution,
            surface_speed_inertia_kg: inertia,
        })
    }
}
