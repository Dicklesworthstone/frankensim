//! Bind one physical baffled exterior to original DG endpoint and side ports.
//!
//! This reuses the existing geometry/gas BEM sampling, rational bake and common
//! load admission. A circular outlet at the origin of an infinite rigid baffle
//! is explicit source physics, not a free-field bell/body approximation.

use super::{DgBoundary, DgBranch, DgDuct, DgPortLoad, invalid};
use crate::bernoulli_aperture::radiation::environment::multiport::input::MimoRadiationInputError;
use crate::bernoulli_aperture::radiation::environment::multiport::source::ComputedBaffledInput;
use crate::pcm_wav::baffled::openings::CircularOpening;
use crate::pcm_wav::baffled::{BaffledPressure, CircularOutletReceiver};
use crate::pcm_wav::baffled::gas::{GasPropagation, GasPropagationReport};
use crate::acoustic_realize::AcousticRealizeError;
use fs_exec::CancelGate;
use fs_duct::dg::PortId;
use fs_duct::profile::gas_field::GasFieldDuctError;
use fs_material::gas_species::input::GasRecord;
use fs_math::det;
use fs_vfit::mimo::load::{MimoLoadLimits, PreparedMimoLoad};

impl DgPortLoad {
    /// Compute a causal radiation load from the actual outlet area and explicit
    /// exterior material, preserving the complete source record. This is cold
    /// preparation; a material substitution recomputes source and runtime once.
    ///
    /// The bore and exterior have one equilibrium static pressure. Different
    /// compositions/temperatures are allowed, but a pressure jump would require
    /// a mean-flow/interface model and is refused by this quiescent adapter.
    ///
    /// The original source bake must resolve the full mechanical Nyquist band.
    /// Its numerical fit and GridOnly checks retain their actual report classes.
    /// Selecting a microphone later does not recompute or alter this reaction.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_computed_baffled(
        duct: DgDuct,
        mechanics_rate_hz: u32,
        exponential_tolerance: f64,
        maximum_dofs: usize,
        exterior: &GasRecord,
        source: ComputedBaffledInput,
        gate: &CancelGate,
    ) -> Result<Self, MimoRadiationInputError> {
        let DgBoundary::Port(port) = duct.outlet() else {
            return Err(MimoRadiationInputError::Physics(invalid(
                "computed DG outlet radiation requires its original port",
            )));
        };
        let radius_m = det::sqrt(duct.outlet_area_m2() / core::f64::consts::PI);
        Self::new_with_computed_baffled_ports(
            duct, mechanics_rate_hz, exponential_tolerance, maximum_dofs,
            exterior, source, &[(port, CircularOpening { center_m: [0.0; 2], radius_m })], gate,
        )
    }

    /// Couple every non-driven DG port to one complete baffled radiation matrix.
    ///
    /// Each tuple binds an original port identity to its actual circular opening
    /// in a shared infinite rigid baffle. Tuple order is retained by the matrix,
    /// accepted terminal frames and microphone. The original BEM/VF owner keeps
    /// all mutual entries; they are not replaced by independent self loads.
    ///
    /// An axial outlet must retain its original bore area. Interior ports are
    /// explicit zero-volume flush openings: their areas are authored separately
    /// from bore area, with no inferred chimney or junction correction. Each
    /// consumes its actual face gas for equilibrium admission. All non-inlet
    /// ports must be supplied exactly once before the expensive source bake.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_computed_baffled_ports(
        duct: DgDuct,
        mechanics_rate_hz: u32,
        exponential_tolerance: f64,
        maximum_dofs: usize,
        exterior: &GasRecord,
        source: ComputedBaffledInput,
        openings: &[(PortId, CircularOpening)],
        gate: &CancelGate,
    ) -> Result<Self, MimoRadiationInputError> {
        Self::new_with_computed_baffled_branches(
            duct, mechanics_rate_hz, exponential_tolerance, maximum_dofs,
            exterior, source, openings, Vec::new(), gate,
        )
    }

    /// Add finite original DG branches between declared bore side ports and
    /// the same complete exterior matrix. Opening geometry and accepted flow
    /// belong to each branch's physical OUTLET. The internal junction is never
    /// substituted as a radiation source. Branch and bore gas may differ at
    /// common equilibrium pressure; their original material states are retained.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_computed_baffled_branches(
        duct: DgDuct,
        mechanics_rate_hz: u32,
        exponential_tolerance: f64,
        maximum_dofs: usize,
        exterior: &GasRecord,
        source: ComputedBaffledInput,
        openings: &[(PortId, CircularOpening)],
        branches: Vec<DgBranch>,
        gate: &CancelGate,
    ) -> Result<Self, MimoRadiationInputError> {
        let dt = 1.0 / f64::from(mechanics_rate_hz);
        Self::new_with_computed_baffled_branches_clock(
            duct, dt, mechanics_rate_hz, exponential_tolerance, maximum_dofs,
            exterior, source, openings, branches, gate,
        )
    }

    /// Retain an arbitrary original physical step for the same baffled load.
    ///
    /// The integer-Hz source baker covers Nyquist using `ceil(1 / dt)` only as
    /// its band ceiling. Both the DG exponential and terminal realization use
    /// the exact supplied `dt`, including their actual sampled-passivity check.
    /// No geometry, gas state or physical clock is rounded to the bake rate.
    /// Exact reciprocal integer clocks retain the existing bake and arithmetic.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_computed_baffled_branches_at_dt(
        duct: DgDuct,
        dt: f64,
        exponential_tolerance: f64,
        maximum_dofs: usize,
        exterior: &GasRecord,
        source: ComputedBaffledInput,
        openings: &[(PortId, CircularOpening)],
        branches: Vec<DgBranch>,
        gate: &CancelGate,
    ) -> Result<Self, MimoRadiationInputError> {
        if gate.is_requested() { return Err(MimoRadiationInputError::Cancelled); }
        let rate = 1.0 / dt;
        let nearest = rate.round();
        let bake_rate = if (1.0..=f64::from(u32::MAX)).contains(&nearest)
            && dt.to_bits() == (1.0 / nearest).to_bits()
        { nearest } else { rate.ceil() };
        if !dt.is_finite() || dt <= 0.0
            || !bake_rate.is_finite() || !(1.0..=f64::from(u32::MAX)).contains(&bake_rate)
        {
            return Err(MimoRadiationInputError::Physics(invalid(
                "computed DG radiation requires a finite positive physical step within its source bake band budget",
            )));
        }
        Self::new_with_computed_baffled_branches_clock(
            duct, dt, bake_rate as u32, exponential_tolerance, maximum_dofs,
            exterior, source, openings, branches, gate,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn new_with_computed_baffled_branches_clock(
        duct: DgDuct,
        dt: f64,
        bake_rate_hz: u32,
        exponential_tolerance: f64,
        maximum_dofs: usize,
        exterior: &GasRecord,
        source: ComputedBaffledInput,
        openings: &[(PortId, CircularOpening)],
        branches: Vec<DgBranch>,
        gate: &CancelGate,
    ) -> Result<Self, MimoRadiationInputError> {
        if gate.is_requested() { return Err(MimoRadiationInputError::Cancelled); }
        if bake_rate_hz == 0 || !dt.is_finite() || dt <= 0.0
            || !matches!(duct.inlet(), DgBoundary::Port(_))
            || duct.total_dofs() > maximum_dofs
            || duct.viscothermal_band_rad_s().is_some_and(|band|
                band[1] * dt > core::f64::consts::PI)
        {
            return Err(MimoRadiationInputError::Physics(invalid(
                "computed DG radiation requires owned ports, a resolved clock/state budget and equal equilibrium pressures",
            )));
        }
        admit_branches(&duct, &branches, dt, maximum_dofs)
            .map_err(MimoRadiationInputError::Physics)?;
        admit_openings(&duct, &branches, exterior, openings)
            .map_err(MimoRadiationInputError::Physics)?;
        let ports: Vec<_> = openings.iter().map(|(port, _)| *port).collect();
        let apertures: Vec<_> = openings.iter().map(|(_, aperture)| *aperture).collect();
        let (model, record) = source.build_apertures(
            exterior, &apertures,
            bake_rate_hz, gate,
        )?;
        let limits = source.limits();
        // The source covers this band, but its integer-Hz ceiling is not the
        // physical step. Admit the actual runtime through its own Nyquist.
        let terminal = PreparedMimoLoad::new(&model.impedance, dt, MimoLoadLimits {
            band_rad_s: [
                core::f64::consts::TAU * model.band_hz[0],
                core::f64::consts::PI / dt,
            ],
            sampled_passivity_points: limits.sampled_passivity_points,
            max_states: limits.max_states,
            max_memory_bytes: limits.max_memory_bytes,
        }).map_err(|error| MimoRadiationInputError::Physics(invalid(error.0)))?;
        if gate.is_requested() { return Err(MimoRadiationInputError::Cancelled); }
        let mut load = Self::new_with_terminal_branches(
            duct, dt, exponential_tolerance, maximum_dofs, ports, terminal, branches,
        ).map_err(MimoRadiationInputError::Physics)?;
        load.radiation_record = Some(record);
        if gate.is_requested() { return Err(MimoRadiationInputError::Cancelled); }
        Ok(load)
    }
}

#[cfg(test)]
mod tests;

fn admit_openings(
    duct: &DgDuct,
    branches: &[DgBranch],
    exterior: &GasRecord,
    openings: &[(PortId, CircularOpening)],
) -> Result<(), AcousticRealizeError> {
    use crate::bernoulli_aperture::radiation::environment::multiport::bake::MAX_BAFFLED_MIMO_PORTS;
    let sealed = branches
        .iter()
        .filter(|branch| branch.duct.outlet() == DgBoundary::Closed)
        .count();
    let expected = duct.side_ports().len()
        + usize::from(matches!(duct.outlet(), DgBoundary::Port(_))) - sealed;
    if openings.is_empty() || openings.len() > MAX_BAFFLED_MIMO_PORTS
        || openings.len() != expected
    {
        return Err(invalid("computed DG radiation must bind every unsealed non-inlet port within the original aperture budget"));
    }
    for (index, (port, aperture)) in openings.iter().enumerate() {
        if duct.inlet() == DgBoundary::Port(*port)
            || branches.iter().any(|branch| branch.port == *port && branch.duct.outlet() == DgBoundary::Closed)
            || openings[..index].iter().any(|(other, _)| other == port)
            || !aperture.radius_m.is_finite() || aperture.radius_m <= 0.0
            || !aperture.center_m.iter().all(|value| value.is_finite())
        {
            return Err(invalid("computed DG openings need distinct owned terminal IDs and finite positive geometry"));
        }
        let endpoint = branches.iter().find(|branch| branch.port == *port)
            .map(|branch| &branch.duct)
            .or_else(|| (duct.outlet() == DgBoundary::Port(*port)).then_some(duct));
        let gas = if let Some(endpoint) = endpoint {
            let area = core::f64::consts::PI * aperture.radius_m * aperture.radius_m;
            let actual = endpoint.outlet_area_m2();
            if !area.is_finite() || (area - actual).abs() > 64.0 * f64::EPSILON * actual {
                return Err(invalid("computed DG opening must retain its actual bore or chimney outlet area"));
            }
            *endpoint.outlet_gas()
        } else {
            duct.side_port_source(*port).map_err(|error| match error {
                GasFieldDuctError::Duct(source) => AcousticRealizeError::Duct(source),
                GasFieldDuctError::Material(source) => AcousticRealizeError::Ambient(source.to_string()),
            })?.gas
        };
        if gas.pressure.to_bits() != exterior.state().pressure.to_bits() {
            return Err(invalid("computed DG openings and their exterior require equal equilibrium pressures"));
        }
    }
    Ok(())
}

// Cheap source admission precedes the original expensive radiation bake.
// Runtime construction performs the same ownership/energy admission on the
// exact prepared branches; this check does not build another propagator.
fn admit_branches(
    duct: &DgDuct,
    branches: &[DgBranch],
    dt: f64,
    maximum_dofs: usize,
) -> Result<(), AcousticRealizeError> {
    super::branch::admit_sources(duct, branches, dt, maximum_dofs)
}

impl DgPortLoad {
    /// Prepare an observer of this computed exterior's accepted terminal flows.
    /// Source geometry, full gas state, clock and radiation/wall bands remain
    /// owned by this load; a microphone never substitutes or re-bakes them.
    /// The original Rayleigh owner admits quadrature, delay/history and its
    /// same-clock differentiating image. This adds no acoustic reaction.
    pub fn baffled_receiver(
        &self,
        mechanics_rate_hz: u32,
        receiver: CircularOutletReceiver,
    ) -> Result<BaffledPressure, AcousticRealizeError> {
        self.admit_baffled_receiver(mechanics_rate_hz, receiver)?;
        let source = self.radiation_record().expect("admitted exterior source");
        BaffledPressure::circular_openings_multirate_differentiating(
            source.apertures(), mechanics_rate_hz, mechanics_rate_hz,
            receiver, source.gas_record().state().into(),
        ).map_err(AcousticRealizeError::Nonlinear)
    }

    /// Observe the same accepted terminal flows with the complete bound exterior
    /// material on every source-to-microphone path. The selected gas law belongs
    /// to the shared Rayleigh owner, never to a reed or jet. It does not change
    /// the physical terminal reaction or add bulk absorption inside the bore.
    pub fn baffled_receiver_with_gas(
        &self, mechanics_rate_hz: u32, receiver: CircularOutletReceiver,
        propagation: &GasPropagation, gate: &CancelGate,
    ) -> Result<(BaffledPressure, GasPropagationReport), AcousticRealizeError> {
        self.admit_baffled_receiver(mechanics_rate_hz, receiver)?;
        let source = self.radiation_record().expect("admitted exterior source");
        BaffledPressure::circular_openings_multirate_gas_differentiating(
            source.apertures(), mechanics_rate_hz, mechanics_rate_hz,
            receiver, source.gas_record().state(), propagation,
            &mut || gate.is_requested(),
        ).map_err(AcousticRealizeError::Nonlinear)
    }

    /// Observe conservatively projected ORIGINAL terminal flows on an
    /// independent regular clock. The caller must perform that explicit time
    /// projection before `step_trace`; this constructor never retimes the DG
    /// reaction or interprets raw physical-step samples at another rate.
    /// Receiver bandwidth must fit both physical and output Nyquist, plus the
    /// retained radiation and material-wall bands. Projection/aliasing error
    /// remains a caller-owned refinement axis, not an observer accuracy claim.
    pub fn baffled_receiver_projected(
        &self, output_rate_hz: u32, receiver: CircularOutletReceiver,
    ) -> Result<BaffledPressure, AcousticRealizeError> {
        self.admit_baffled_receiver_clock(output_rate_hz, receiver, true)?;
        let source = self.radiation_record().expect("admitted exterior source");
        BaffledPressure::circular_openings_multirate_differentiating(
            source.apertures(), output_rate_hz, output_rate_hz,
            receiver, source.gas_record().state().into(),
        ).map_err(AcousticRealizeError::Nonlinear)
    }

    /// The same explicit terminal-flow projection with original gas absorption
    /// and dispersion on every retained source-to-receiver path.
    pub fn baffled_receiver_projected_with_gas(
        &self, output_rate_hz: u32, receiver: CircularOutletReceiver,
        propagation: &GasPropagation, gate: &CancelGate,
    ) -> Result<(BaffledPressure, GasPropagationReport), AcousticRealizeError> {
        self.admit_baffled_receiver_clock(output_rate_hz, receiver, true)?;
        let source = self.radiation_record().expect("admitted exterior source");
        BaffledPressure::circular_openings_multirate_gas_differentiating(
            source.apertures(), output_rate_hz, output_rate_hz,
            receiver, source.gas_record().state(), propagation,
            &mut || gate.is_requested(),
        ).map_err(AcousticRealizeError::Nonlinear)
    }

    fn admit_baffled_receiver(
        &self, mechanics_rate_hz: u32, receiver: CircularOutletReceiver,
    ) -> Result<(), AcousticRealizeError> {
        self.admit_baffled_receiver_clock(mechanics_rate_hz, receiver, false)
    }

    fn admit_baffled_receiver_clock(
        &self, mechanics_rate_hz: u32, receiver: CircularOutletReceiver,
        projected: bool,
    ) -> Result<(), AcousticRealizeError> {
        if mechanics_rate_hz == 0 || self.state().step_count() != 0
            || (!projected && self.propagator().dt().to_bits()
                != (1.0 / f64::from(mechanics_rate_hz)).to_bits())
        {
            return Err(invalid("DG receiver requires a cold source and its original physical clock"));
        }
        let source = self.radiation_record().ok_or_else(|| invalid(
            "baffled DG observation requires the original computed radiation source",
        ))?;
        // The private source record and ordered terminal basis were bound
        // together at construction; observation cannot substitute geometry.
        if source.driven_mouth().is_some_and(|mouth|
            receiver.radial_rings != mouth.terminal_radial_rings
                || receiver.angular_points != mouth.terminal_angular_points)
        {
            return Err(invalid("joint DG receiver must retain the original terminal surface quadrature"));
        }
        if source.apertures().is_empty()
            || source.apertures().len() != self.terminal_ports().len()
            || (!projected && source.band_hz()[1].to_bits()
                != (0.5 * f64::from(mechanics_rate_hz)).to_bits())
        {
            return Err(invalid(
                "baffled DG observation must use its original terminal geometry and physical clock",
            ));
        }
        // Exact integer source clocks retain their original endpoint even
        // when the reciprocal dt rounds by an ULP (for example 8001 Hz).
        // A fractional physical clock still uses its actual Nyquist; this
        // equality does not admit a source baked below that band.
        let physical_nyquist = if self.propagator().dt().to_bits()
            == (1.0 / (2.0 * source.band_hz()[1])).to_bits()
        { source.band_hz()[1] } else { 0.5 / self.propagator().dt() };
        if projected && (source.band_hz()[1] < physical_nyquist
            || receiver.maximum_frequency_hz > physical_nyquist
            || receiver.maximum_frequency_hz > 0.5 * f64::from(mechanics_rate_hz))
        {
            return Err(invalid("projected DG receiver must fit original physical, source-bake and output Nyquist bands"));
        }
        // Compare in the wall owner's native angular-frequency arithmetic.
        // Converting its endpoint back to Hz can reject the identical endpoint
        // solely because multiplication and division round differently.
        if receiver.maximum_frequency_hz > source.band_hz()[1]
            || self.minimum_wall_band_upper_rad_s().is_some_and(|upper|
                receiver.maximum_frequency_hz * core::f64::consts::TAU > upper)
        {
            return Err(invalid("baffled DG receiver exceeds its radiation or wall-loss band"));
        }
        Ok(())
    }
}
