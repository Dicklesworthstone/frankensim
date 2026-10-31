//! Spatially resolved friction across a rigid or material-fiber ribbon.
//!
//! The width is a physical interval along the string, discretized by midpoint
//! quadrature. Each cell has its OWN slip and traction in the existing joint
//! solve; the string velocity is never averaged before applying friction.
//! Every contact uses the original material-derived transverse/torsional map.
//! `with_uniform_bow_ribbon` supplies the rigid, uniformly loaded limit;
//! `with_tensioned_fiber_ribbon` retains reciprocal fiber compliance. Neither
//! supplies thermal rosin or a validated violin. Refine cells, modes and time.
//! Finite-width and point contacts need not produce the same motion: see
//! Pitteroff & Woodhouse, Acta Acustica 84 (1998), 744–757.

use super::{MaterialBowContact, MultiBowStringBank, MultiBowStringBankSource};
use super::super::{ReciprocalStringBank, StringBowLocation};
use crate::bowed_string::runtime::reciprocal::{ReciprocalStringError, invalid, poll};
use crate::render::schedule::force::coupled::contact::friction::{
    FrictionDrive, multiple::MultiFrictionConfig,
};
use fs_exec::CancelGate;
use crate::fiber_ribbon::{FiberDirection, TensionedFiberRibbon};
use crate::render::schedule::force::coupled::ModalAttachment;
use crate::render::schedule::force::coupled::contact::friction::multiple::surfaces::FrictionSurfaceAttachment;

fn fiber_memory_error(error: crate::fiber_ribbon::FiberRibbonError) -> ReciprocalStringError {
    match error {
        crate::fiber_ribbon::FiberRibbonError::Invalid(what) => invalid(what),
        crate::fiber_ribbon::FiberRibbonError::Cancelled => ReciprocalStringError::Network(
            crate::render::schedule::force::coupled::ModalCouplingError::Cancelled),
        error => ReciprocalStringError::Material(
            crate::acoustic_realize::AcousticRealizeError::Nonlinear(error.to_string())),
    }
}

/// Translating fiber supports, their inertial loads and moving material contact.
pub mod motion;

/// Material fibers present as mechanical components of the SAME string bank.
#[derive(Clone, Debug)]
pub struct MaterialFiberRibbon {
    model: TensionedFiberRibbon,
    component_start: usize,
    string_count: usize,
    station_fraction: f64,
    normal_shapes: Vec<f64>,
}
impl MaterialFiberRibbon {
    /// Original geometry, material sources and two independent resolutions.
    pub const fn model(&self) -> &TensionedFiberRibbon { &self.model }
    /// First actual fiber component after the original string/body/fluid states.
    pub const fn component_start(&self) -> usize { self.component_start }
    /// Contact material-point station along each fiber's supported length.
    pub const fn station_fraction(&self) -> f64 { self.station_fraction }
    /// Normal maps for `MaterialStringStop::body`, in its string-relative body
    /// indexing. Use these with the SAME ribbon cells and contact orientation;
    /// `with_bow_contact_loads` then consumes the solved normal reactions.
    /// A normal law, clearance and direction are explicit independent inputs.
    pub fn normal_stop_receivers(&self) -> Vec<ModalAttachment> {
        (0..self.model.resolution().cells).map(|cell| ModalAttachment {
            component: self.component_start + cell - self.string_count,
            shapes: self.normal_shapes.clone(),
        }).collect()
    }
}

/// Frozen uniform line load and its actual spatial contact quadrature.
/// Contact zero is the original centre; other cells follow from left to right.
#[derive(Clone, Debug, PartialEq)]
pub struct UniformBowRibbon {
    centre: StringBowLocation,
    width_m: f64,
    stations: Vec<f64>,
}

impl UniformBowRibbon {
    /// Resolve physical width into centre-first midpoint contact stations.
    /// This same projection is used by source-file contact bindings and runtime
    /// preparation, so a ribbon cell is never reconstructed from rounded text.
    pub fn try_new(
        centre: StringBowLocation,
        length_m: f64,
        width_m: f64,
        cells: usize,
    ) -> Result<Self, ReciprocalStringError> {
        // Odd quadrature retains the already bound physical centre as port zero.
        if !width_m.is_finite() || width_m <= 0.0 || !length_m.is_finite()
            || length_m <= 0.0 || !(1..=63).contains(&cells) || cells % 2 == 0
        {
            return Err(invalid("ribbon needs positive finite width and an odd cell count in 1..=63"));
        }
        let half = 0.5 * width_m / length_m;
        let left = centre.station_fraction - half;
        let right = centre.station_fraction + half;
        if !left.is_finite() || !right.is_finite() || left <= 0.0 || right >= 1.0
            || !(left < centre.station_fraction && centre.station_fraction < right)
        {
            return Err(invalid("the complete finite ribbon must lie strictly inside its actual string"));
        }
        let step = width_m / length_m / cells as f64;
        let mid = cells / 2;
        let mut stations = Vec::with_capacity(cells);
        stations.push(centre.station_fraction);
        let mut previous = left;
        for cell in 0..cells {
            let station = centre.station_fraction + (cell as f64 - mid as f64) * step;
            if !(previous < station && station < right) {
                return Err(invalid("ribbon cells must have distinct representable interior stations"));
            }
            previous = station;
            if cell != mid { stations.push(station); }
        }
        Ok(Self { centre, width_m, stations })
    }

    /// Actual material string and fixed ribbon centre.
    #[must_use]
    pub const fn centre(&self) -> StringBowLocation { self.centre }
    /// Physical contact width along that string [m], not a tonal parameter.
    #[must_use]
    pub const fn width_m(&self) -> f64 { self.width_m }
    /// Actual contact stations, in the joint solver's port order.
    #[must_use]
    pub fn stations(&self) -> &[f64] { &self.stations }
    /// Equal midpoint-quadrature share of the total normal force.
    #[must_use]
    pub fn load_fraction(&self) -> f64 { 1.0 / self.stations.len() as f64 }

    /// Convert ONE total ribbon load to the original solver's per-cell drives.
    /// The prescribed speed is common; each cell still solves its own slip.
    /// Output storage is caller-owned and must cover every spatial cell.
    pub fn project_drive(
        &self, drive: FrictionDrive, output: &mut [FrictionDrive],
    ) -> Result<(), ReciprocalStringError> {
        if output.len() != self.stations.len() || !drive.speed_m_s.is_finite()
            || !drive.normal_force_n.is_finite() || drive.normal_force_n < 0.0
        {
            return Err(invalid("ribbon drive requires finite speed, nonnegative total load and one output per cell"));
        }
        output.fill(FrictionDrive { speed_m_s: drive.speed_m_s,
            normal_force_n: drive.normal_force_n * self.load_fraction() });
        Ok(())
    }
}

impl MultiBowStringBankSource {
    /// Finite-width source when this bank was constructed from a uniform ribbon.
    #[must_use]
    pub const fn ribbon(&self) -> Option<&UniformBowRibbon> { self.ribbon.as_ref() }
    /// Original source of the reciprocally moving fiber surfaces, if present.
    pub const fn fiber_ribbon(&self) -> Option<&MaterialFiberRibbon> { self.fiber_ribbon.as_ref() }
}

impl ReciprocalStringBank {
    /// Replace rigid ribbon surfaces with source-derived fiber motion.
    ///
    /// Every width cell receives its own axial and normal continuum, retained
    /// in the original exact mechanical owner. Opposite tangential reactions
    /// stretch the fibers and feed back into each string contact in the SAME
    /// joint solve; their material losses enter the existing energy account.
    /// Selected fiber Young spectra append real first-order material histories
    /// to that same owner in both axial and bending directions. Existing
    /// string/body/cavity histories remain an unchanged prefix.
    /// No post-filter or independently stepped bow is used. The original body,
    /// cavity, radiation, material memory and nonlinear strain storage remain
    /// attached. Appended fiber modes leave the original string and retained
    /// longitudinal coordinates in place and join their constitutive solve.
    /// An explicitly bound fiber strain law couples its normal deflection and
    /// axial extension through that same owner. It retains all existing string
    /// nonlinear storage and requires the same numerical solve configuration.
    /// Selected Maxwell spectra then follow that complete strain, including
    /// mean geometric stretch, while bending keeps its independent histories.
    ///
    /// Bind cold. Contact station is fixed along the
    /// fibers: this is a local frozen-station bowing model in a uniformly
    /// translating support frame, not a complete stroke/moving-frog model.
    /// `FrictionDrive::speed_m_s` is that frame's prescribed constant speed.
    /// Nonuniform support acceleration and contact migration are not inferred.
    /// Explicit `with_fiber_support_motion` adds these from endpoint velocities
    /// and original material mass for the prescribed-load dry-friction rung.
    /// A performance's `with_fiber_support_program` supplies continuous velocity
    /// knots to that same owner; laboratory work includes its momentum correction.
    /// Normal load is still prescribed unless the returned normal receivers
    /// are explicitly composed with `with_string_stops` and
    /// `with_bow_contact_loads`. Dry friction is not a thermal rosin law.
    pub fn with_tensioned_fiber_ribbon(
        mut self, ribbon: TensionedFiberRibbon, station_fraction: f64,
        joint: MultiFrictionConfig, gate: &CancelGate,
    ) -> Result<MultiBowStringBank, ReciprocalStringError> {
        poll(gate)?;
        let cells = ribbon.resolution().cells;
        let extra_modes = ribbon.cell_model().modes().len().checked_mul(cells)
            .and_then(|n| n.checked_add(self.system.mode_count()));
        if extra_modes.is_none_or(|n| n > self.source.options.network.max_modes)
            || cells > joint.max_ports
            || ribbon.cell_model().sample_period_s().to_bits() != self.system.sample_period_s().to_bits()
        {
            return Err(invalid("material ribbon must fit the original physical clock and network/contact budgets"));
        }
        let axial = ribbon.displacement_map(station_fraction, FiberDirection::Axial)
            .map_err(|_| invalid("invalid axial fiber contact map"))?;
        let normal_shapes = ribbon.displacement_map(station_fraction, FiberDirection::Normal)
            .map_err(|_| invalid("invalid normal fiber contact map"))?;
        let component_start = self.system.components().len();
        let first_fiber_mode = self.system.mode_count();
        let string_count = self.cards.len();
        self.system = self.system.with_internal_components(ribbon.cell_models(), gate)
            .map_err(ReciprocalStringError::Friction)?;
        if let Some(memory) = ribbon.material_relaxation() {
            let branches = memory.branches(first_fiber_mode, self.system.mode_count(), gate)
                .map_err(fiber_memory_error)?;
            self.system = self.system.with_additional_material_relaxation(branches, None, gate)
                .map_err(ReciprocalStringError::Friction)?;
        }
        if let Some((bindings, config)) = ribbon.strain_bindings(component_start) {
            self.system = self.system.with_mapped_nonlinear_potentials(bindings, config, gate)
                .map_err(ReciprocalStringError::Friction)?;
        }
        if let Some(memory) = ribbon.material_relaxation() {
            let bindings = memory.strain_bindings(component_start);
            if !bindings.is_empty() {
                // The source budget already covers all ribbon bending and
                // strain states. Preserve earlier nonlinear histories as well;
                // the joint owner independently enforces its 256-state ceiling.
                let count = self.system.nonlinear_material_memory().iter().map(|m| m.len()).sum::<usize>()
                    .checked_add(bindings.iter().map(|b| b.law.memory_count()).sum::<usize>())
                    .ok_or_else(|| invalid("complete fiber material history count overflows"))?;
                self.system = self.system.with_mapped_material_histories(bindings, count, gate)
                    .map_err(ReciprocalStringError::Friction)?;
            }
            memory.validate_band(self.system.linear_frequency_upper_bound_rad_s(),
                self.system.sample_period_s()).map_err(fiber_memory_error)?;
            if self.source.longitudinal_waves.is_some() {
                super::super::nonlinear::longitudinal::relaxation::finite_strain::validate_sources(
                    &self.source, &self.system, gate).map_err(|error| match error {
                        crate::render::plate::impact::ImpactError::Cancelled =>
                            ReciprocalStringError::Network(crate::render::schedule::force::coupled::ModalCouplingError::Cancelled),
                        error => ReciprocalStringError::Material(
                            crate::acoustic_realize::AcousticRealizeError::Nonlinear(error.to_string())),
                    })?;
            } else {
                self.source.validate_body_memory_frequencies(&self.system, gate)?;
                self.source.validate_material_bridge_frequencies(&self.system, gate)?;
            }
        }
        let mut bank = self.with_uniform_bow_ribbon(ribbon.geometry().width_m, cells, joint, gate)?;
        let attachments = (0..cells).map(|port| FrictionSurfaceAttachment { port,
            surface: ModalAttachment { component: component_start + port, shapes: axial.clone() },
        }).collect();
        bank.system = bank.system.with_surface_attachments(attachments, gate)
            .map_err(ReciprocalStringError::Friction)?;
        bank.source.fiber_ribbon = Some(MaterialFiberRibbon {
            model: ribbon, component_start, string_count, station_fraction, normal_shapes,
        });
        Ok(bank)
    }

    /// Resolve a uniform finite-width ribbon on the existing bow's string.
    ///
    /// The original material/body/cavity state, clock and observer survive.
    /// Every new point port participates in the SAME simultaneous implicit
    /// friction solve. The width does not filter sound or prescribe oscillation.
    /// The supplied odd cell count is a numerical resolution, not a hair count.
    ///
    /// An already bound interface refuses: its nominal contact area belongs to
    /// the whole original contact. Bind explicit per-cell interface areas after
    /// this construction with `with_bow_interface`, retaining the original law.
    /// Hair compliance and a mechanically resolved normal-load distribution are
    /// separate physics; this method prescribes uniform loading and rigid speed.
    pub fn with_uniform_bow_ribbon(
        self, width_m: f64, cells: usize, joint: MultiFrictionConfig, gate: &CancelGate,
    ) -> Result<MultiBowStringBank, ReciprocalStringError> {
        poll(gate)?;
        if self.bow_interface().is_some() {
            return Err(invalid("bind ribbon interfaces per cell after resolving the contact area"));
        }
        if cells > joint.max_ports {
            return Err(invalid("ribbon cells exceed the declared joint contact budget"));
        }
        let centre = self.source.options.bow;
        let ribbon = UniformBowRibbon::try_new(centre, self.cards[centre.string].length_m, width_m, cells)?;
        let law = self.system.port().law.clone();
        let friction = self.source.options.friction;
        let additional = ribbon.stations[1..].iter().map(|&station_fraction| MaterialBowContact {
            location: StringBowLocation { string: centre.string, station_fraction },
            law: law.clone(), friction,
        }).collect();
        let mut bank = self.with_additional_bows(additional, joint, gate)?;
        bank.source.ribbon = Some(ribbon);
        Ok(bank)
    }
}
