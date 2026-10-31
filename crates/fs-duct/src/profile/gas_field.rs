//! Continuous, prescribed gas fields on the existing Webster/Magnus owner.
//!
//! The material field retains primitive T/composition and source receipts.
//! A cold binding evaluates that field at the two actual quadrature positions
//! of every Magnus step. Frequency sweeps then reuse those admitted states;
//! no density, sound-speed or transport coefficient is interpolated. Local
//! two-ports compose through the existing stable scattering owner.

use super::{MAX_RICCATI_STEPS, RiccatiSteps, bad, local_generator, magnus_step, positive};
use crate::compose::{
    ComposeOptions, CompositionMode, InputReflection, PreparedChain, PreparedOpenHole, TwoPort,
    compose_chain, compose_chain_reflection, open_hole_flows, prepare_open_hole,
};
use crate::{Duct, DuctError, DuctResponse, LossModel, OpenHoleFlow, Segment, Termination};
use fs_material::MaterialError;
use fs_material::gas::GasState;
use fs_material::gas_field::{AxialGasField, AxialGasSampling};
use fs_math::{c64::C64, det};
use fs_phs::WallPin;
use crate::WallModel;

/// Original source or propagation refusal, without a fallback material law.
#[derive(Debug)]
pub enum GasFieldDuctError {
    /// Primitive state, composition, source domain or phase refusal.
    Material(MaterialError),
    /// Geometry, integration budget or propagation refusal.
    Duct(DuctError),
}
impl core::fmt::Display for GasFieldDuctError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Material(error) => write!(f, "gas field: {error}"),
            Self::Duct(error) => write!(f, "gas-field duct: {error}"),
        }
    }
}
impl std::error::Error for GasFieldDuctError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Material(error) => Some(error),
            Self::Duct(error) => Some(error),
        }
    }
}
impl From<MaterialError> for GasFieldDuctError {
    fn from(error: MaterialError) -> Self {
        Self::Material(error)
    }
}
impl From<DuctError> for GasFieldDuctError {
    fn from(error: DuctError) -> Self {
        Self::Duct(error)
    }
}

#[derive(Clone, Debug)]
struct QuadratureStep {
    /// Negative because the physical ABCD map propagates outlet to inlet.
    h: f64,
    radii: [f64; 2],
    /// Absolute positions of the original material/geometry quadrature.
    positions_m: [f64; 2],
    /// Physical wall area per axial area at the same two geometry samples.
    slants: [f64; 2],
    gases: [GasState; 2],
}

#[derive(Clone, Debug)]
struct ToneHoleStep {
    /// Keep the zero-length junction at its exact authored axial station.
    before_step: usize,
    segment: usize,
    axial_position_m: f64,
    gas: GasState,
    extra_series_m: f64,
    bore_one_way_time_s: f64,
}

/// Read-only geometry/material sample shared with the molecular constitutive
/// binding. The original field owner alone chooses knots and quadrature.
pub(crate) enum MolecularFieldStep<'a> {
    Magnus { h: f64, radii: [f64; 2], positions_m: [f64; 2], slants: [f64; 2] },
    Hole { geometry: &'a Segment, axial_position_m: f64 },
}

/// Physical open chimney and its travel time through the original material
/// field. These are observation coordinates, not an extra propagation model.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasFieldOpening {
    /// Axial coordinate from the inlet [m].
    pub axial_position_m: f64,
    /// Radiating radius including the declared open area fraction [m].
    pub effective_radius_m: f64,
    /// Integral of 1/c to this junction on the original propagation quadrature.
    pub bore_one_way_time_s: f64,
    /// Bore travel plus physical/inner-matching chimney length at its local c.
    pub one_way_time_s: f64,
}

/// Acoustic binding of existing geometry to the one shared material field.
///
/// This is a frozen, no-mean-flow plane-wave operator. An explicit wall pin
/// adds the shared locally reacting material law at each quadrature point.
/// Temperature
/// and species transport, condensation and transverse acoustic modes remain
/// with their respective owners. A uniform field delegates the original duct
/// verbatim, including its existing profile/cone axes and operation order.
#[derive(Clone, Debug)]
pub struct GasFieldDuct {
    duct: Duct,
    field: AxialGasField,
    quadrature: Vec<QuadratureStep>,
    holes: Vec<ToneHoleStep>,
    inlet: GasState,
    outlet: GasState,
    one_way_time_s: f64,
}

impl GasFieldDuct {
    /// Bind the complete field to a cylinder/cone/smooth-profile chain with
    /// optional local tone-hole junctions.
    ///
    /// For a nonuniform field, steps_per_span controls the shared Magnus
    /// integration on EACH smooth interval. Authored geometry boundaries,
    /// primitive field stations, and spline knots are exact step boundaries;
    /// derivative jumps are never hidden inside a fourth-order step. The total
    /// work is bounded by MAX_RICCATI_STEPS, and steps_used exposes it.
    /// This explicit field axis replaces the uniform profile's step axis when
    /// the field varies; a uniform field needs no field discretization.
    ///
    /// All local material states are admitted before this binding is returned.
    /// Each tone hole uses the exact gas at its axial station for its junction
    /// and homogeneous chimney, as in the regional-duct model. Closed holes may
    /// lie anywhere in the field. Multiple open/vented holes require a constant
    /// primitive state throughout their connecting span, preserving the existing
    /// homogeneous mutual-hole law. A gradient across that span refuses; no
    /// averaged medium or omitted mutual correction is substituted.
    pub fn new(
        duct: &Duct,
        field: &AxialGasField,
        steps_per_span: RiccatiSteps,
    ) -> Result<Self, GasFieldDuctError> {
        Self::new_controlled(duct, field, steps_per_span, &mut || false)
    }

    /// As new, with cancellation before each pair of source evaluations.
    /// A refused or cancelled construction publishes no partial binding.
    pub fn new_controlled(
        duct: &Duct,
        field: &AxialGasField,
        steps_per_span: RiccatiSteps,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, GasFieldDuctError> {
        steps_per_span.validate()?;
        if duct.segments.is_empty() {
            return Err(DuctError::EmptyDuct.into());
        }
        if duct.segments.len() > MAX_RICCATI_STEPS {
            return Err(bad("gas-field geometry exceeds the integration work cap").into());
        }
        for segment in &duct.segments {
            segment.validate()?;
        }
        let length: f64 = duct.segments.iter().map(Segment::axial_length).sum();
        if !positive(length) || length != field.length_m() {
            return Err(bad("gas field must cover the exact complete duct length").into());
        }
        if cancelled() {
            return Err(DuctError::Cancelled { completed_steps: 0 }.into());
        }
        let inlet = *field.state_at(0.0)?.state();
        let outlet = *field.state_at(length)?.state();
        let mut result = Self {
            duct: duct.clone(),
            field: field.clone(),
            quadrature: Vec::new(),
            holes: Vec::new(),
            inlet,
            outlet,
            one_way_time_s: length / inlet.sound_speed,
        };
        if field.uniform_state().is_some() {
            return Ok(result);
        }
        let mut open_sites = Vec::new();
        let mut position = 0.0;
        for segment in &duct.segments {
            if let Segment::ToneHole {
                hole_radius,
                bore_radius,
                state,
                ..
            } = *segment
            {
                let sigma = crate::hole_sigma(state);
                if sigma > 0.0 {
                    open_sites.push((position, hole_radius * sigma.sqrt(), bore_radius));
                }
            }
            position += segment.axial_length();
        }
        if let (Some(first), Some(last)) = (open_sites.first(), open_sites.last()) {
            // Equal endpoint coefficients alone are insufficient: temperature
            // or composition could change between them. Inspect the actual
            // reconstruction on every intersecting span. A cellwise interval
            // holds its RIGHT station; its left neighbor need not be identical.
            let first_index = field.stations_m().partition_point(|x| *x < first.0);
            let first_state = &field.primitive_states()[first_index];
            for (positions, states) in field
                .stations_m()
                .windows(2)
                .zip(field.primitive_states().windows(2))
            {
                let varies = match field.sampling() {
                    AxialGasSampling::LinearPrimitives => states[0] != states[1],
                    AxialGasSampling::CellwiseConstant => &states[1] != first_state,
                };
                if first.0 < last.0
                    && positions[1] > first.0
                    && positions[0] < last.0
                    && varies
                {
                    return Err(bad("open tone holes across a varying gas field need an unavailable mutual-hole law").into());
                }
            }
        }
        let root3 = det::sqrt(3.0);
        let (c1, c2) = (0.5 - root3 / 6.0, 0.5 + root3 / 6.0);
        let mut offset = 0.0;
        let mut travel = 0.0;
        let mut open_index = 0;
        for (segment_index, segment) in duct.segments.iter().enumerate() {
            if let Segment::ToneHole { state, .. } = *segment {
                if result.steps_used() == MAX_RICCATI_STEPS {
                    return Err(bad("gas-field integration exceeds MAX_RICCATI_STEPS").into());
                }
                if cancelled() {
                    return Err(DuctError::Cancelled {
                        completed_steps: result.steps_used(),
                    }
                    .into());
                }
                let mut extra_series_m = 0.0;
                if crate::hole_sigma(state) > 0.0 {
                    let (xi, bi, ai) = open_sites[open_index];
                    for &(xj, bj, _) in &open_sites {
                        extra_series_m +=
                            fs_phs::side_hole_mutual_length(bi, bj, ai, (xj - xi).abs());
                    }
                    open_index += 1;
                }
                result.holes.push(ToneHoleStep {
                    before_step: result.quadrature.len(),
                    segment: segment_index,
                    axial_position_m: offset,
                    gas: *field.state_at(offset)?.state(),
                    extra_series_m,
                    bore_one_way_time_s: travel,
                });
                continue;
            }
            let length = segment.axial_length();
            let end = offset + length;
            let mut boundaries = vec![0.0, length];
            boundaries.extend(
                field
                    .stations_m()
                    .iter()
                    .filter_map(|&x| (x > offset && x < end).then_some(x - offset)),
            );
            if let Segment::Profile { profile, .. } = segment {
                if let Some(stations) = profile.stations() {
                    let origin = stations[0].0;
                    boundaries.extend(
                        stations[1..stations.len() - 1]
                            .iter()
                            .map(|&(x, _)| x - origin),
                    );
                }
            }
            boundaries.sort_by(f64::total_cmp);
            boundaries.dedup();
            let count = (boundaries.len() - 1)
                .checked_mul(steps_per_span.0)
                .and_then(|count| result.quadrature.len().checked_add(count))
                .filter(|&count| count <= MAX_RICCATI_STEPS - result.holes.len())
                .ok_or_else(|| bad("gas-field integration exceeds MAX_RICCATI_STEPS"))?;
            result.quadrature.reserve(count - result.quadrature.len());
            let radius = |x: f64| -> Result<f64, DuctError> {
                match segment {
                    Segment::Cylinder { radius, .. } => Ok(*radius),
                    Segment::Cone {
                        inlet_radius,
                        outlet_radius,
                        ..
                    } => Ok(*inlet_radius + (*outlet_radius - *inlet_radius) * (x / length)),
                    Segment::Profile { profile, .. } => profile.radius(x),
                    Segment::ToneHole { .. } => unreachable!("axial geometry admitted"),
                }
            };
            let slant = |x: f64| -> Result<f64, DuctError> {
                let slope = match segment {
                    Segment::Cylinder { .. } => 0.0,
                    Segment::Cone { inlet_radius, outlet_radius, .. } => {
                        (outlet_radius - inlet_radius) / length
                    }
                    Segment::Profile { profile, .. } => profile.slope(x)?,
                    Segment::ToneHole { .. } => unreachable!("axial geometry admitted"),
                };
                let value = det::sqrt(1.0 + slope * slope);
                if !positive(value) {
                    return Err(bad("gas-field wall surface slope is not representable"));
                }
                Ok(value)
            };
            for span in boundaries.windows(2) {
                for index in 0..steps_per_span.0 {
                    if cancelled() {
                        return Err(DuctError::Cancelled {
                            completed_steps: result.steps_used(),
                        }
                        .into());
                    }
                    let left =
                        span[0] + (span[1] - span[0]) * (index as f64 / steps_per_span.0 as f64);
                    let right = if index + 1 == steps_per_span.0 {
                        span[1]
                    } else {
                        span[0]
                            + (span[1] - span[0]) * ((index + 1) as f64 / steps_per_span.0 as f64)
                    };
                    let h = left - right;
                    if !positive(-h) {
                        return Err(bad(
                            "gas-field integration interval is numerically unresolved",
                        )
                        .into());
                    }
                    let positions = [right + c1 * h, right + c2 * h];
                    let gases = [
                        *field.state_at(offset + positions[0])?.state(),
                        *field.state_at(offset + positions[1])?.state(),
                    ];
                    let radii = [radius(positions[0])?, radius(positions[1])?];
                    let slants = [slant(positions[0])?, slant(positions[1])?];
                    travel += -0.5 * h * (1.0 / gases[0].sound_speed + 1.0 / gases[1].sound_speed);
                    result.quadrature.push(QuadratureStep {
                        h, radii, slants, gases,
                        positions_m: [offset + positions[0], offset + positions[1]],
                    });
                }
            }
            offset = end;
        }
        if !positive(travel) {
            return Err(bad("gas-field travel time must be positive and finite").into());
        }
        result.one_way_time_s = travel;
        Ok(result)
    }

    /// Original source-backed primitive field, retained without resampling.
    pub const fn field(&self) -> &AxialGasField {
        &self.field
    }
    /// Original geometry, including authored profile knots.
    pub const fn duct(&self) -> &Duct {
        &self.duct
    }
    /// Exact source state at the inlet, also used by the characteristic port.
    pub const fn inlet_gas(&self) -> &GasState {
        &self.inlet
    }
    /// Exact source state at the outlet, the default exterior radiation gas.
    pub const fn outlet_gas(&self) -> &GasState {
        &self.outlet
    }
    /// Integral of 1/c using the same field quadrature as propagation.
    pub const fn one_way_time_s(&self) -> f64 {
        self.one_way_time_s
    }

    /// Open and vented mouths in authored order, with no closed-hole sources.
    /// Travel times use the same admitted primitive field and quadrature as the
    /// feedback operator; they never use an averaged or exterior sound speed.
    pub fn open_hole_sites(&self) -> Vec<GasFieldOpening> {
        let mut sites = Vec::new();
        let mut axial_position_m = 0.0;
        let mut steps = self.holes.iter();
        for segment in &self.duct.segments {
            if let Segment::ToneHole { hole_radius, chimney_height, bore_radius, state } = *segment {
                let (gas, time) = if self.field.uniform_state().is_some() {
                    (&self.inlet, axial_position_m / self.inlet.sound_speed)
                } else {
                    let step = steps.next().expect("each admitted junction retains its local gas");
                    (&step.gas, step.bore_one_way_time_s)
                };
                let sigma = crate::hole_sigma(state);
                if sigma > 0.0 {
                    let inner = fs_phs::side_hole_inner_length(hole_radius, bore_radius);
                    sites.push(GasFieldOpening {
                        axial_position_m,
                        effective_radius_m: hole_radius * sigma.sqrt(),
                        bore_one_way_time_s: time,
                        one_way_time_s: time + (chimney_height + inner) / gas.sound_speed,
                    });
                }
            }
            axial_position_m += segment.axial_length();
        }
        sites
    }
    /// Actual field integration blocks, including zero-length tone-hole
    /// junctions; zero for the original uniform path.
    pub fn steps_used(&self) -> usize {
        self.quadrature.len() + self.holes.len()
    }

    /// Traverse original quadrature and exact authored junctions without
    /// copying or inventing a second field-integration mesh.
    pub(crate) fn molecular_steps(&self) -> impl Iterator<Item = MolecularFieldStep<'_>> {
        let mut axial = 0;
        let mut hole = 0;
        std::iter::from_fn(move || {
            if let Some(junction) = self.holes.get(hole).filter(|h| h.before_step == axial) {
                hole += 1;
                return Some(MolecularFieldStep::Hole {
                    geometry: &self.duct.segments[junction.segment],
                    axial_position_m: junction.axial_position_m,
                });
            }
            let step = self.quadrature.get(axial)?;
            axial += 1;
            Some(MolecularFieldStep::Magnus {
                h: step.h, radii: step.radii, positions_m: step.positions_m,
                slants: step.slants,
            })
        })
    }
    /// Exact inlet rho*c/area, with no independently supplied medium knob.
    pub fn inlet_impedance(&self) -> f64 {
        let radius = self.duct.segments[0].inlet_radius();
        self.inlet.density * self.inlet.sound_speed / (core::f64::consts::PI * radius * radius)
    }

    /// Resolve local two-ports in inlet-to-outlet order. Every finite step
    /// shares the existing Magnus and exact-ZK owners; composition never forms
    /// a global growing ABCD matrix. Callers may refine the field step axis if
    /// an individual exponential is unresolved.
    pub fn chain_controlled(
        &self,
        omega: f64,
        loss: LossModel,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedChain, DuctError> {
        self.chain_wall_controlled(omega, loss, None, cancelled)
    }

    /// The same stable local two-ports with an explicit per-area wall law.
    /// Geometry, gas and wall admittance are evaluated together; gas-field
    /// changes do not discard the material wall or change its constitutive law.
    /// None retains the rigid-wall operation order.
    pub fn chain_wall_controlled(
        &self,
        omega: f64,
        loss: LossModel,
        wall: Option<&WallPin>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedChain, DuctError> {
        self.chain_with_wall_controlled(omega, loss, WallModel::from_pin(wall), cancelled)
    }

    /// Complete wall-frequency law at every retained material/geometry sample.
    pub fn chain_with_wall_controlled(
        &self,
        omega: f64,
        loss: LossModel,
        wall: WallModel<'_>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedChain, DuctError> {
        self.chain_with_wall_exterior_controlled(omega, loss, None, wall, cancelled)
    }

    /// The same field and chimney propagation with an independently declared
    /// exterior radiation gas at every exposed hole. `None` preserves the
    /// original local-gas radiation convention; no gas is averaged or replaced.
    pub fn chain_with_wall_exterior_controlled(
        &self,
        omega: f64,
        loss: LossModel,
        exterior: Option<&GasState>,
        wall: WallModel<'_>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedChain, DuctError> {
        self.chain_with_holes_controlled(omega, loss, exterior, wall, false, cancelled)
            .map(|(chain, _)| chain)
    }

    fn chain_with_holes_controlled(
        &self,
        omega: f64,
        loss: LossModel,
        exterior: Option<&GasState>,
        wall: WallModel<'_>,
        observe_holes: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(PreparedChain, Vec<PreparedOpenHole>), DuctError> {
        if let Some(exterior) = exterior {
            crate::exterior::admit(&self.inlet, exterior)?;
        }
        if self.field.uniform_state().is_some() {
            return crate::compose::duct_chain_in_exterior(&self.duct, &self.inlet,
                exterior.unwrap_or(&self.inlet), omega, loss, wall,
                ComposeOptions::default(), true, cancelled);
        }
        if !positive(omega) {
            return Err(bad("gas-field frequency must be positive and finite"));
        }
        let reference = self.inlet_impedance();
        if !positive(reference) {
            return Err(bad(
                "gas-field reference impedance must be positive and finite",
            ));
        }
        let root3 = det::sqrt(3.0);
        let mut blocks = Vec::with_capacity(self.steps_used());
        let mut observed_holes = Vec::new();
        let mut min_shear = f64::INFINITY;
        let open_sites = crate::collect_open_hole_sites(&self.duct);
        let mut axial_index = 0;
        let mut hole_index = 0;
        let mut open_hole_index = 0;
        for index in 0..self.steps_used() {
            if cancelled() {
                return Err(DuctError::Cancelled {
                    completed_steps: index,
                });
            }
            if let Some(hole) = self
                .holes
                .get(hole_index)
                .filter(|hole| hole.before_step == axial_index)
            {
                let Segment::ToneHole {
                    hole_radius,
                    chimney_height,
                    bore_radius,
                    state,
                } = self.duct.segments[hole.segment]
                else {
                    unreachable!("admitted tone hole")
                };
                if matches!(loss, LossModel::WideTube | LossModel::Bessel) {
                    min_shear = min_shear.min(
                        hole_radius
                            * (hole.gas.density * omega / hole.gas.dynamic_viscosity).sqrt(),
                    );
                }
                if crate::hole_sigma(state) > 0.0 {
                    if observe_holes {
                        observed_holes.push(prepare_open_hole(
                            &hole.gas, exterior.unwrap_or(&hole.gas), omega, loss, wall,
                            &self.duct.segments[hole.segment], open_sites[open_hole_index],
                            blocks.len(), hole.extra_series_m, true,
                        )?);
                    }
                    open_hole_index += 1;
                }
                let matrix = crate::tone_hole_t_junction_in_exterior_with_wall(
                    &hole.gas,
                    exterior.unwrap_or(&hole.gas),
                    hole_radius,
                    chimney_height,
                    state,
                    omega,
                    loss,
                    bore_radius,
                    hole.extra_series_m,
                    wall,
                )?;
                blocks.push(TwoPort::from_abcd(matrix, reference)?);
                hole_index += 1;
                continue;
            }
            let step = &self.quadrature[axial_index];
            let (a, shear_a) = local_generator(
                &step.gases[0], step.radii[0], omega, loss, reference, wall, step.slants[0],
            )?;
            let (b, shear_b) = local_generator(
                &step.gases[1], step.radii[1], omega, loss, reference, wall, step.slants[1],
            )?;
            let mut matrix = magnus_step(a, b, step.h, root3)?;
            matrix[1] = matrix[1].scale(reference);
            matrix[2] = matrix[2].scale(1.0 / reference);
            blocks.push(TwoPort::from_abcd(matrix, reference)?);
            min_shear = min_shear.min(shear_a).min(shear_b);
            axial_index += 1;
        }
        Ok((PreparedChain {
            blocks,
            reference_impedance: reference,
            min_shear_number: min_shear,
        }, observed_holes))
    }

    fn exterior<'a>(&'a self, gas: Option<&'a GasState>) -> Result<&'a GasState, DuctError> {
        let gas = gas.unwrap_or(&self.outlet);
        gas.validate_coherence()
            .map_err(|_| bad("gas-field exterior must be a coherent material state"))?;
        if gas.pressure != self.inlet.pressure {
            return Err(bad("gas field and exterior need a common static pressure"));
        }
        Ok(gas)
    }

    /// Input impedance and physical outlet transfer with the existing analytic
    /// termination owner. Exterior replaces the homogeneous exterior only;
    /// omitted means the exact x=L material, never a guessed air state.
    pub fn response(
        &self,
        omega: f64,
        loss: LossModel,
        termination: Termination,
        exterior: Option<&GasState>,
    ) -> Result<DuctResponse, DuctError> {
        self.response_wall(omega, loss, termination, exterior, None)
    }

    /// Input and outlet response including the existing compliant-wall law.
    /// The wall state is independent of the gas mixture; changing either data
    /// source recomputes their coupled local acoustic generator.
    pub fn response_wall(
        &self,
        omega: f64,
        loss: LossModel,
        termination: Termination,
        exterior: Option<&GasState>,
        wall: Option<&WallPin>,
    ) -> Result<DuctResponse, DuctError> {
        self.response_with_wall(omega, loss, termination, exterior, WallModel::from_pin(wall))
    }

    /// Input and outlet response with the complete source-owned wall law.
    pub fn response_with_wall(
        &self, omega: f64, loss: LossModel, termination: Termination,
        exterior: Option<&GasState>, wall: WallModel<'_>,
    ) -> Result<DuctResponse, DuctError> {
        self.response_with_holes(omega, loss, termination, exterior, wall, false)
            .map(|(response, _)| response)
    }

    /// The same feedback response plus physical open/vented chimney flows.
    /// Each flow is reconstructed from its actual place in the Magnus cascade,
    /// retaining local material, wall, mutual correction and exterior loading.
    /// Closed holes remain acoustic loads and never become radiating sources.
    pub fn response_with_wall_and_holes(
        &self, omega: f64, loss: LossModel, termination: Termination,
        exterior: Option<&GasState>, wall: WallModel<'_>,
    ) -> Result<(DuctResponse, Vec<OpenHoleFlow>), DuctError> {
        self.response_with_holes(omega, loss, termination, exterior, wall, true)
    }

    fn response_with_holes(
        &self, omega: f64, loss: LossModel, termination: Termination,
        exterior: Option<&GasState>, wall: WallModel<'_>, observe_holes: bool,
    ) -> Result<(DuctResponse, Vec<OpenHoleFlow>), DuctError> {
        let radiation = self.exterior(exterior)?;
        let radius = self
            .duct
            .segments
            .last()
            .expect("admitted duct")
            .outlet_radius();
        if self.field.uniform_state().is_some() {
            if observe_holes {
                return crate::exterior::input_impedance_and_holes(&self.duct, &self.inlet,
                    radiation, omega, loss, termination, wall);
            }
            return crate::exterior::input_impedance(&self.duct, &self.inlet, radiation,
                omega, loss, termination, wall).map(|response| (response, Vec::new()));
        }
        let (chain, holes) = self.chain_with_holes_controlled(
            omega, loss, exterior, wall, observe_holes, &mut || false,
        )?;
        let (load, mouth_ka) = crate::exterior::termination_impedance(
            termination, &self.outlet, radiation, radius, omega)?;
        let solved = compose_chain(&chain.blocks, load, ComposeOptions::default())?;
        let flows = open_hole_flows(&chain, &holes, load, ComposeOptions::default())?;
        Ok((DuctResponse {
            omega,
            impedance: solved.z_in,
            min_shear_number: chain.min_shear_number,
            mouth_ka,
            u_mouth_over_u_in: solved.u_out_over_u_in,
            p_mouth_over_u_in: if load.is_none() {
                C64::ZERO
            } else {
                solved.p_out_over_u_in
            },
            condition_estimate: solved.condition_estimate,
            composition: CompositionMode::ImpedanceRecursion,
        }, flows))
    }

    /// Direct input reflection for characteristic-line construction. The
    /// existing stable owner does not divide through an input impedance pole.
    pub fn reflection(
        &self,
        omega: f64,
        loss: LossModel,
        termination: Termination,
        exterior: Option<&GasState>,
    ) -> Result<InputReflection, DuctError> {
        self.reflection_wall(omega, loss, termination, exterior, None)
    }

    /// Wall-aware direct reflection, avoiding an intermediate input-impedance
    /// pole just as the rigid-wall path does.
    pub fn reflection_wall(
        &self,
        omega: f64,
        loss: LossModel,
        termination: Termination,
        exterior: Option<&GasState>,
        wall: Option<&WallPin>,
    ) -> Result<InputReflection, DuctError> {
        self.reflection_with_wall(omega, loss, termination, exterior, WallModel::from_pin(wall))
    }

    /// Direct reflection with the same complete material wall as the response.
    pub fn reflection_with_wall(
        &self, omega: f64, loss: LossModel, termination: Termination,
        exterior: Option<&GasState>, wall: WallModel<'_>,
    ) -> Result<InputReflection, DuctError> {
        let radiation = self.exterior(exterior)?;
        let chain = self.chain_with_wall_exterior_controlled(omega, loss, exterior, wall, &mut || false)?;
        let radius = self
            .duct
            .segments
            .last()
            .expect("admitted duct")
            .outlet_radius();
        let (load, _) = crate::exterior::termination_impedance(
            termination, &self.outlet, radiation, radius, omega)?;
        compose_chain_reflection(&chain.blocks, load, ComposeOptions::default())
    }
}

#[cfg(test)]
pub(crate) mod tests;
