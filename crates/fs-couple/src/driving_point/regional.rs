//! Prescribed gas fields connected at TMM pressure/volume-flow ports.
//!
//! Each constant-medium run is evaluated by `fs-duct`, including its cones,
//! wall law and tone-hole junctions. Its downstream load is the next run's input
//! impedance. Pressure and volume flow are continuous at interfaces: no averaged
//! gas, artificial reflection coefficient or additional propagation law is used.
//! A continuous binding uses fs-duct's shared Magnus field operator through
//! these same source, outlet and characteristic-line ports.
//! This is frozen, no-mean-flow acoustics, not thermal/composition transport.

/// Physical outlet observation of the existing regional transfer.
pub mod outlet;
/// Independent exterior gas at the terminal pressure/volume-flow port.
pub mod exterior;
/// Primitive gas-field sampling through the existing regional propagation law.
pub mod field;
/// Source-bound membrane walls on the same regional/continuous propagation.
pub mod material_wall;
mod baked;

use super::{
    DrivingPointError, LineRealization, PoleState, bulk_delay_pole_state,
    characteristic_line_response_with_poles, duct_reflectance, impedance_reflectance,
    sealed_volume_compliance_gain,
};
use fs_duct::{
    Duct, DuctError, DuctResponse, HoleState, LossModel, MouthLoad, OpenHoleFlow, Segment,
    Termination, WallModel,
};
use fs_material::gas::GasState;
use fs_phs::WallPin;
use fs_vfit::discretize::DelayedFilter;

/// Cold-path cap, including axial segments and zero-length tone-hole records.
/// It admits a complete measured bore table with its split stations and chimneys.
pub const MAX_REGIONAL_SEGMENTS: usize = 256;

#[derive(Debug, Clone)]
struct Region {
    duct: Duct,
    gas: GasState,
}

/// Static regional duct composed from the existing uniform TMM owners.
///
/// Adjacent identical gases coalesce, preserving the exact uniform calculation.
/// Every open/vented tone hole must lie in the SAME constant-medium run so that
/// `fs-duct` still evaluates every mutual-hole correction. Cross-region mutual
/// hole laws are not available and refuse rather than silently disappearing.
#[derive(Debug, Clone)]
pub struct RegionalDuct {
    regions: Vec<Region>,
    gases: Vec<GasState>,
    one_way_time_s: f64,
    exterior_gas: Option<GasState>,
    inlet_gas_override: Option<GasState>,
    continuous: Option<fs_duct::profile::gas_field::GasFieldDuct>,
    closing_pressure_pa: f64,
}

fn invalid(what: &'static str) -> DuctError { DuctError::BadParameter { what } }

impl RegionalDuct {
    /// Bind exactly one resolved state to every authored geometry record.
    ///
    /// States must be coherent and at exactly the same static pressure. Species
    /// and temperature may vary. A tone-hole record uses its assigned gas for
    /// both its compact bore junction and chimney. Every exposed mouth uses
    /// one homogeneous exterior, defaulting to the actual outlet gas; an
    /// explicit binding replaces only that exterior. Geometry/radiation/loss validity remains checked by `fs-duct`
    /// at each requested frequency. Construction does not advance sound state.
    ///
    /// # Errors
    /// Empty/oversized or incomplete binding, incoherent gases, unequal pressure,
    /// nonpositive axial lengths, nonaxial ends, or unsupported mutual-hole scope.
    pub fn new(duct: &Duct, gases: &[GasState]) -> Result<Self, DuctError> {
        if duct.segments.is_empty() { return Err(DuctError::EmptyDuct); }
        if duct.segments.len() > MAX_REGIONAL_SEGMENTS || gases.len() != duct.segments.len() {
            return Err(invalid("regional duct requires 1..=256 segments and one gas per segment"));
        }
        if matches!(duct.segments.first(), Some(Segment::ToneHole { .. }))
            || matches!(duct.segments.last(), Some(Segment::ToneHole { .. })) {
            return Err(invalid("regional duct ends must be explicit axial segments"));
        }
        let mut regions: Vec<Region> = Vec::new();
        let mut open_hole_region = None;
        for (segment, gas) in duct.segments.iter().zip(gases) {
            gas.validate_coherence().map_err(|_| invalid("regional gas is not a coherent material state"))?;
            if gas.pressure != gases[0].pressure {
                return Err(invalid("regional duct needs a common static pressure; mean flow is not modeled"));
            }
            if regions.last().is_none_or(|r| r.gas != *gas) {
                regions.push(Region { duct: Duct { segments: Vec::new() }, gas: *gas });
            }
            match *segment {
                Segment::Profile { .. } => {}
                Segment::Cylinder { length, .. } | Segment::Cone { length, .. } => {
                    if !length.is_finite() || length <= 0.0 {
                        return Err(invalid("regional axial lengths must be positive and finite"));
                    }
                }
                Segment::ToneHole { state, .. } => {
                    let open = match state {
                        HoleState::Open => true,
                        HoleState::Closed => false,
                        HoleState::Vent(s) => {
                            if !s.is_finite() || !(0.0..=1.0).contains(&s) {
                                return Err(invalid("regional vent fraction must be finite and in [0,1]"));
                            }
                            s > 0.0
                        }
                    };
                    if open {
                        let index = regions.len() - 1;
                        if open_hole_region.is_some_and(|previous| previous != index) {
                            return Err(invalid("open tone holes in different gas regions need an unavailable mutual-hole law"));
                        }
                        open_hole_region = Some(index);
                    }
                }
            }
            regions.last_mut().expect("run just admitted").duct.segments.push(segment.clone());
        }
        let one_way_time_s = regions.iter().map(|region| {
            let length: f64 = region.duct.segments.iter().map(Segment::axial_length).sum();
            length / region.gas.sound_speed
        }).sum::<f64>();
        if !one_way_time_s.is_finite() || one_way_time_s <= 0.0 {
            return Err(invalid("regional travel time must be positive and finite"));
        }
        Ok(Self { regions, gases: gases.to_vec(), one_way_time_s, exterior_gas: None,
            inlet_gas_override: None, continuous: None,
            closing_pressure_pa: fs_duct::DEFAULT_REED_CLOSING_PRESSURE_PA })
    }

    /// Bind the static reed closing pressure $P_M = K_a y_m$ [Pa] used for the
    /// Raman/Dalmont threshold acoustic pressure amplitude
    /// $p_{\text{ac,peak}} = (\sqrt{2}/3)\, P_M$ in the open-hole vortex
    /// describing function (`Dalmont, Gilbert & Ollivier 2003` Eq. 12).
    #[must_use]
    pub const fn with_closing_pressure_pa(mut self, closing_pressure_pa: f64) -> Self {
        self.closing_pressure_pa = closing_pressure_pa;
        self
    }

    /// Static reed closing pressure $P_M$ [Pa] bound to this regional duct.
    #[must_use]
    pub const fn closing_pressure_pa(&self) -> f64 {
        self.closing_pressure_pa
    }

    pub(super) fn compose_options(&self) -> fs_duct::compose::ComposeOptions {
        fs_duct::compose::ComposeOptions {
            closing_pressure_pa: self.closing_pressure_pa,
            ..fs_duct::compose::ComposeOptions::default()
        }
    }

    /// Resolved media in authored segment order, including hole records.
    /// For a continuous field these are exact authored-segment INLET samples;
    /// propagation uses the field operator's quadrature, not these summaries.
    #[must_use]
    pub fn segment_gases(&self) -> &[GasState] { &self.gases }

    /// Actual inlet medium for Bernoulli density and the characteristic port.
    #[must_use]
    pub fn inlet_gas(&self) -> &GasState { self.inlet_gas_override.as_ref().unwrap_or(&self.gases[0]) }

    /// Actual internal medium at the terminal axial segment. Radiation uses
    /// [`Self::exterior_gas`], which defaults to this state.
    #[must_use]
    pub fn outlet_gas(&self) -> &GasState {
        self.continuous.as_ref().map_or_else(
            || self.gases.last().expect("nonempty binding"),
            |field| field.outlet_gas(),
        )
    }

    /// Sum of constant-medium travel times, or the continuous operator's
    /// quadrature integral of local 1/c along its authored geometry.
    /// This excludes side-branch chimney length, as does the uniform line owner.
    #[must_use]
    pub const fn one_way_time_s(&self) -> f64 { self.one_way_time_s }

    /// Evaluate one frequency with the original TMM loss, wall and load laws.
    ///
    /// # Errors
    /// The existing duct owner refuses geometry, frequency, loss or radiation;
    /// or an intermediate port is nonfinite (for example an exact lossless pole).
    pub fn response(&self, omega: f64, loss: LossModel, termination: Termination,
        wall: Option<&WallPin>) -> Result<DuctResponse, DuctError> {
        self.response_with_wall(omega, loss, termination, WallModel::from_pin(wall))
    }

    /// Complete local wall response with the same gas, geometry and terminal load.
    pub fn response_with_wall(&self, omega: f64, loss: LossModel, termination: Termination,
        wall: WallModel<'_>) -> Result<DuctResponse, DuctError> {
        if !omega.is_finite() || omega <= 0.0 {
            return Err(invalid("frequency must be positive and finite"));
        }
        let exterior = self.exterior_gas();
        if let Some(field) = &self.continuous {
            return field.response_with_wall(omega, loss, termination, Some(exterior), wall);
        }
        let options = self.compose_options();
        let mut runs = self.regions.iter().rev();
        let terminal = runs.next().expect("nonempty binding");
        let mut response = self.terminal_response(terminal, omega, loss, termination, wall)?;
        for region in runs {
            let upstream = fs_duct::exterior::input_impedance_load_options(
                &region.duct,
                &region.gas,
                exterior,
                omega,
                loss,
                response.impedance,
                wall,
                options,
            )?;
            // A downstream one-port receives the upstream run's output flow.
            // Keep the REAL mouth's ka and load; this interface is not another
            // radiating termination. Units are pressure and VOLUME flow on both
            // sides, so no density, area or power-wave renormalization is needed.
            response.u_mouth_over_u_in = response.u_mouth_over_u_in * upstream.u_mouth_over_u_in;
            response.p_mouth_over_u_in = response.p_mouth_over_u_in * upstream.u_mouth_over_u_in;
            response.impedance = upstream.impedance;
            response.min_shear_number = response.min_shear_number.min(upstream.min_shear_number);
        }
        Ok(response)
    }

    /// [`Self::response_with_wall`] plus outward volume-velocity transfers at every
    /// open or vented tone hole's flanged mouth, solved from the same regional chain.
    pub fn response_with_wall_and_holes(
        &self,
        omega: f64,
        loss: LossModel,
        termination: Termination,
        wall: WallModel<'_>,
    ) -> Result<(DuctResponse, Vec<OpenHoleFlow>), DuctError> {
        if !omega.is_finite() || omega <= 0.0 {
            return Err(invalid("frequency must be positive and finite"));
        }
        let exterior = self.exterior_gas();
        if let Some(field) = &self.continuous {
            return field.response_with_wall_and_holes(
                omega, loss, termination, Some(exterior), wall,
            );
        }
        let options = self.compose_options();
        let mut offsets = Vec::with_capacity(self.regions.len());
        let mut axial_x = 0.0;
        for region in &self.regions {
            offsets.push(axial_x);
            axial_x += region
                .duct
                .segments
                .iter()
                .map(Segment::axial_length)
                .sum::<f64>();
        }
        let mut runs = self.regions.iter().enumerate().rev();
        let (terminal_idx, terminal) = runs.next().expect("nonempty binding");
        let (mut response, mut holes) =
            self.terminal_response_with_holes(terminal, omega, loss, termination, wall)?;
        for hole in &mut holes {
            hole.axial_position_m += offsets[terminal_idx];
        }
        for (region_idx, region) in runs {
            let (upstream, mut upstream_holes) =
                fs_duct::exterior::input_impedance_load_and_holes_options(
                    &region.duct,
                    &region.gas,
                    exterior,
                    omega,
                    loss,
                    response.impedance,
                    wall,
                    options,
                )?;
            for hole in &mut holes {
                hole.u_hole_over_u_in = hole.u_hole_over_u_in * upstream.u_mouth_over_u_in;
                hole.u_chimney_in_over_u_in =
                    hole.u_chimney_in_over_u_in * upstream.u_mouth_over_u_in;
            }
            if !upstream_holes.is_empty() {
                for hole in &mut upstream_holes {
                    hole.axial_position_m += offsets[region_idx];
                }
                upstream_holes.extend(holes);
                holes = upstream_holes;
            }
            response.u_mouth_over_u_in = response.u_mouth_over_u_in * upstream.u_mouth_over_u_in;
            response.p_mouth_over_u_in = response.p_mouth_over_u_in * upstream.u_mouth_over_u_in;
            response.impedance = upstream.impedance;
            response.min_shear_number = response.min_shear_number.min(upstream.min_shear_number);
        }
        Ok((response, holes))
    }

    /// Reflection at the original inlet port, with the complete local wall.
    /// Continuous fields and uniform regions use their stable reflection owner.
    pub fn reflection_with_wall(&self, omega: f64, loss: LossModel,
        termination: Termination, wall: WallModel<'_>) -> Result<fs_math::c64::C64, DuctError> {
        let zc = self.inlet_impedance();
        let exterior = self.exterior_gas();
        if let Some(field) = &self.continuous {
            return field.reflection_with_wall(omega, loss, termination,
                Some(exterior), wall)?.referenced_to(zc);
        }
        if self.regions.len() == 1 {
            let region = &self.regions[0];
            return fs_duct::exterior::input_reflection_options(
                &region.duct,
                &region.gas,
                exterior,
                omega,
                loss,
                termination,
                wall,
                self.compose_options(),
            )?
            .referenced_to(zc);
        }
        self.response_with_wall(omega, loss, termination, wall)
            .map(|response| impedance_reflectance(&response, zc))
    }

    /// Realize the regional driving-point reflectance using the SAME FIR owner.
    ///
    /// The physical return time uses local sound speeds. Port impedance uses
    /// the inlet gas and inlet area. `n_fft_override` has exactly the uniform
    /// line's meaning, budgets and finite-history approximation.
    ///
    /// # Errors
    /// As [`Self::response`], plus original history/Fourier/filter admission.
    pub fn characteristic_line<'a>(&self, termination: impl Into<MouthLoad<'a>>, sample_rate_hz: u32,
        n: usize, wall: Option<&WallPin>, n_fft_override: Option<usize>)
        -> Result<DelayedFilter, DrivingPointError> {
        self.characteristic_line_realized(termination, sample_rate_hz, n, wall, n_fft_override)
            .map(|(line, _)| line)
    }

    /// Realize the regional driving-point reflectance and return the [`LineRealization`] receipt.
    ///
    /// # Errors
    /// As [`Self::characteristic_line`].
    pub fn characteristic_line_realized<'a>(&self, termination: impl Into<MouthLoad<'a>>, sample_rate_hz: u32,
        n: usize, wall: Option<&WallPin>, n_fft_override: Option<usize>)
        -> Result<(DelayedFilter, LineRealization), DrivingPointError> {
        self.characteristic_line_with_wall_realized(termination, sample_rate_hz, n,
            WallModel::from_pin(wall), n_fft_override)
    }

    /// Original FIR owner with one complete source-owned wall binding.
    pub fn characteristic_line_with_wall<'a>(&self, termination: impl Into<MouthLoad<'a>>, sample_rate_hz: u32,
        n: usize, wall: WallModel<'_>, n_fft_override: Option<usize>)
        -> Result<DelayedFilter, DrivingPointError> {
        self.characteristic_line_with_wall_realized(termination, sample_rate_hz, n,
            wall, n_fft_override).map(|(line, _)| line)
    }

    /// Complete wall-frequency response and the original line-realization receipt.
    pub fn characteristic_line_with_wall_realized<'a>(&self, termination: impl Into<MouthLoad<'a>>, sample_rate_hz: u32,
        n: usize, wall: WallModel<'_>, n_fft_override: Option<usize>)
        -> Result<(DelayedFilter, LineRealization), DrivingPointError> {
        let termination = match termination.into() {
            MouthLoad::Analytic(t) => t,
            MouthLoad::Tabulated(table) => return self.tabulated_line(table, sample_rate_hz, n, wall, n_fft_override),
            MouthLoad::Matched => return Err(super::baked::matched_refusal()),
        };
        // Delegate constant pins verbatim, including travel-time arithmetic.
        // A closed/ideal-open terminal does not remove radiation at side holes.
        if let WallModel::Rigid | WallModel::Pin(_) = wall
            && self.continuous.is_none() && self.regions.len() == 1 && (self.exterior_gas() == self.outlet_gas()
            || (matches!(termination, Termination::Closed | Termination::IdealOpen)
                && !self.regions[0].duct.segments.iter().any(|segment| {
                    matches!(segment, Segment::ToneHole { state, .. }
                        if matches!(state, HoleState::Open)
                            || matches!(state, HoleState::Vent(fraction) if *fraction > 0.0))
                }))) {
            let pin = match wall { WallModel::Pin(pin) => Some(pin), _ => None };
            let r = &self.regions[0];
            let length = super::duct_axial_length(&r.duct);
            let dt = 1.0 / f64::from(sample_rate_hz);
            let geo_delay = 2.0 * length / r.gas.sound_speed / dt;
            let mut pole_states: Vec<PoleState> = Vec::new();
            if let Some(gain) = sealed_volume_compliance_gain(&r.duct, &r.gas, termination, pin) {
                pole_states.push(PoleState::Compliance { gain });
            }
            pole_states.push(bulk_delay_pole_state(geo_delay));
            let zc = self.inlet_impedance();
            return characteristic_line_response_with_poles(
                geo_delay,
                termination,
                sample_rate_hz,
                n,
                zc,
                n_fft_override,
                pole_states,
                |omega| {
                    let loss = LossModel::Bessel;
                    if self.closing_pressure_pa == fs_duct::DEFAULT_REED_CLOSING_PRESSURE_PA {
                        duct_reflectance(&r.duct, &r.gas, omega, loss, termination, pin, zc)
                    } else {
                        let options = self.compose_options();
                        match input_impedance_with_wall_options(&r.duct, &r.gas, omega, loss, termination, wall, options) {
                            Ok(response) => Ok(impedance_reflectance(&response, zc)),
                            Err(DuctError::Singular) => {
                                fs_duct::compose::input_reflection_with_wall_options(&r.duct, &r.gas, omega, loss, termination, wall, options)?.referenced_to(zc)
                            }
                            Err(e) => Err(e),
                        }
                    }
                },
            );
        }
        let dt = 1.0 / f64::from(sample_rate_hz);
        let geo_delay = 2.0 * self.one_way_time_s / dt;
        let pole_states = vec![bulk_delay_pole_state(geo_delay)];
        let zc = self.inlet_impedance();
        // The continuous operator supplies R directly through its stable local
        // two-ports. The retained piecewise path still composes through Z.
        characteristic_line_response_with_poles(
            geo_delay,
            termination,
            sample_rate_hz,
            n,
            zc,
            n_fft_override,
            pole_states,
            |omega| {
                self.reflection_with_wall(omega, LossModel::Bessel, termination, wall)
            },
        )
    }

    /// Plane-wave impedance rho*c/S at the original inlet geometry.
    #[must_use]
    pub fn inlet_impedance(&self) -> f64 {
        let radius = self.regions[0].duct.segments[0].inlet_radius();
        let area = core::f64::consts::PI * radius * radius;
        self.inlet_gas().density * self.inlet_gas().sound_speed / area
    }
}

#[cfg(test)]
mod tests;
