//! Wall radiation reaction on original pressure/volume-flow network ports.
//!
//! Each unflanged compact outlet keeps the native duct self impedance. The
//! original 3D monopole Green function supplies every free-field mutual term,
//! and the original spherical impedance-wall owner supplies every reflected
//! term, including the finite self return. The COMPLETE load enters the native
//! network source solve; pressure-driven outlet flows therefore react to the
//! wall and to each other. A reflected correction alone need not be passive.
//!
//! This is the compact unflanged-mouth approximation, not a finite-body BEM
//! solve. Caller-declared ka and radius/distance limits bound its application;
//! they are not error certificates. Mouth orientation, finite-source multipoles,
//! pipe-body scattering, mean flow and nonlinear source changes are unresolved.
//! No fitted delay, absorption filter or replacement source impedance is added.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::refraction::full_wave::FullWaveMonopoleLimits;
use fs_aeroac::curle3d::refraction::full_wave::wall::{
    ImpedanceWallDiagnostics, ImpedanceWallFieldKind, PlanarImpedanceWall,
    PreparedImpedanceWallCompactSourceObserver, PreparedImpedanceWallMonopole,
};
use fs_aeroac::curle3d::{
    CoherentCompactField, CompactSourcePhasor, PreparedCompactSourceObserver,
};
use fs_duct::fluid::FluidDuctError;
use fs_duct::network::{
    CoupledTermination, DuctNetwork, MAX_NETWORK_PORTS, NetworkCoherentResponse,
    NetworkSourceResponse, PortSource, PortTermination, network_drive_sources,
};
use fs_duct::{DuctError, LossModel, Termination, compact_radiation_impedance_from_material};
use fs_material::fluid_acoustics::FluidAcoustics;
use fs_material::gas::GasState;
use fs_math::{c64::C64, dd::Dd};

/// One original network outlet, with its unflanged circular-mouth radius.
#[derive(Clone, Debug, PartialEq)]
pub struct CompactWallOutlet {
    pub label: String,
    pub position_m: [f64; 3],
    pub radius_m: f64,
}

/// Explicit compact-model admission and aggregate reflected-spectral work.
#[derive(Clone, Copy, Debug)]
pub struct CompactWallRadiationLimits {
    pub maximum_outlets: usize,
    /// Sum over independent reciprocal pairs, including diagonal self returns.
    pub maximum_total_evaluations: usize,
    /// Positive and at most 0.5; choose a smaller bound for the required accuracy.
    pub maximum_ka: f64,
    /// Positive and at most 0.25. Bounds radius / wall clearance and the sum
    /// of two outlet radii / centre separation, keeping mouths disjoint.
    pub maximum_radius_over_distance: f64,
    pub monopole: FullWaveMonopoleLimits,
}

#[derive(Debug)]
pub enum CompactWallRadiationError {
    Invalid(&'static str),
    Cancelled,
    Acoustic(AeroacError),
    Duct(DuctError),
    Fluid(FluidDuctError),
}
impl core::fmt::Display for CompactWallRadiationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "compact wall radiation refused: {self:?}")
    }
}
impl std::error::Error for CompactWallRadiationError {}

/// Complete reciprocal exterior load and its unchanged source/field owners.
#[derive(Clone, Debug)]
pub struct PreparedCompactWallRadiation {
    wall: PlanarImpedanceWall,
    frequency_hz: f64,
    outlets: Vec<CompactWallOutlet>,
    free_space: CoupledTermination,
    reflected: Vec<C64>,
    termination: CoupledTermination,
    diagnostics: Vec<ImpedanceWallDiagnostics>,
    evaluations: usize,
}

impl PreparedCompactWallRadiation {
    /// Prepare one frequency on the wall owner's original exp(-i omega t)
    /// peak-phasor clock. `source_id` identifies the caller's outlet assembly;
    /// original radius, position, order, fluid and wall data remain retained.
    ///
    /// Self impedances use `Termination::UnflangedOpen`; off-diagonal free-field
    /// terms and ALL wall returns use the same physical Q [m³/s] normalization.
    /// Reciprocity evaluates each pair once. Every Q kernel shares the aggregate
    /// work cap; failure or cancellation publishes no partial load. Full-matrix
    /// passivity is checked after combination without clipping or repair.
    pub fn new(
        wall: PlanarImpedanceWall,
        frequency_hz: f64,
        outlets: &[CompactWallOutlet],
        source_id: impl Into<String>,
        limits: CompactWallRadiationLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, CompactWallRadiationError> {
        poll(cancelled)?;
        let source_id = source_id.into();
        let n = outlets.len();
        let omega = core::f64::consts::TAU * frequency_hz;
        let k = omega / wall.fluid().sound_speed_m_s();
        if n == 0
            || n > limits.maximum_outlets
            || n > MAX_NETWORK_PORTS
            || !positive(omega)
            || !positive(k)
            || source_id.trim().is_empty()
            || limits.maximum_total_evaluations == 0
            || !positive(limits.maximum_ka)
            || limits.maximum_ka > 0.5
            || !positive(limits.maximum_radius_over_distance)
            || limits.maximum_radius_over_distance > 0.25
        {
            return Err(invalid(
                "need bounded compact outlets, positive frequency and explicit compact/work limits",
            ));
        }
        let mut free_space = CoupledTermination {
            omega,
            labels: outlets.iter().map(|p| p.label.clone()).collect(),
            z: vec![C64::ZERO; n * n],
            source_id: format!("{source_id}:unflanged-free-space"),
        };
        let mut reflected = free_space.z.clone();
        let mut diagnostics = Vec::with_capacity(n * (n + 1) / 2);
        let mut evaluations = 0;
        let exact = wall.impedance_pa_s_m().is_none() || wall.impedance_pa_s_m() == Some(C64::ZERO);
        for (i, outlet) in outlets.iter().enumerate() {
            let clearance = -dot_distance(outlet.position_m, wall.point_m(), wall.outward_normal());
            if outlet.label.trim().is_empty()
                || outlets[..i].iter().any(|p| p.label == outlet.label)
                || outlet.position_m.iter().any(|x| !x.is_finite())
                || !positive(outlet.radius_m)
                || k * outlet.radius_m > limits.maximum_ka
                || !positive(clearance)
                || outlet.radius_m / clearance > limits.maximum_radius_over_distance
            {
                return Err(invalid(
                    "outlet label, radius, ka or original wall clearance is outside the compact model",
                ));
            }
            free_space.z[i * n + i] = compact_radiation_impedance_from_material(
                Termination::UnflangedOpen,
                &wall.fluid(),
                outlet.radius_m,
                omega,
            )
            .map_err(CompactWallRadiationError::Fluid)?;
            for (j, other) in outlets.iter().enumerate().take(i) {
                let d = distance(outlet.position_m, other.position_m);
                if !positive(d)
                    || (outlet.radius_m / d + other.radius_m / d)
                        > limits.maximum_radius_over_distance
                {
                    return Err(invalid(
                        "outlet separation is outside the compact mutual-radiation model",
                    ));
                }
                let value = PreparedCompactSourceObserver::new(
                    wall.fluid(),
                    frequency_hz,
                    outlet.position_m,
                    &[other.position_m],
                    1,
                )
                .and_then(|p| p.evaluate_harmonic_line(&[CompactSourcePhasor::monopole(C64::ONE)]))
                .map_err(CompactWallRadiationError::Acoustic)?
                .pressure_pa;
                free_space.z[i * n + j] = value;
                free_space.z[j * n + i] = value;
            }
        }
        for i in 0..n {
            for j in i..n {
                poll(cancelled)?;
                let remaining = limits.maximum_total_evaluations - evaluations;
                if remaining == 0 && !exact {
                    return Err(invalid(
                        "compact wall load exhausted its aggregate spectral budget",
                    ));
                }
                let response = PreparedImpedanceWallMonopole::new(
                    wall,
                    frequency_hz,
                    outlets[j].position_m,
                    outlets[i].position_m,
                    FullWaveMonopoleLimits {
                        maximum_evaluations: limits
                            .monopole
                            .maximum_evaluations
                            .min(remaining.max(1)),
                        ..limits.monopole
                    },
                    cancelled,
                )
                .map_err(CompactWallRadiationError::Acoustic)?;
                let z = response
                    .evaluate_reflected_harmonic_line(C64::ONE)
                    .map_err(CompactWallRadiationError::Acoustic)?
                    .pressure_pa;
                reflected[i * n + j] = z;
                reflected[j * n + i] = z;
                evaluations += response.diagnostics().evaluations;
                diagnostics.push(*response.diagnostics());
            }
        }
        free_space
            .validate()
            .map_err(CompactWallRadiationError::Duct)?;
        let mut termination = free_space.clone();
        termination.source_id = format!("{source_id}:unflanged-impedance-wall");
        for (z, &returning) in termination.z.iter_mut().zip(&reflected) {
            *z = C64::new(
                (Dd::from_f64(z.re) + Dd::from_f64(returning.re)).to_f64(),
                (Dd::from_f64(z.im) + Dd::from_f64(returning.im)).to_f64(),
            );
        }
        termination
            .validate()
            .map_err(CompactWallRadiationError::Duct)?;
        poll(cancelled)?;
        Ok(Self {
            wall,
            frequency_hz,
            outlets: outlets.to_vec(),
            free_space,
            reflected,
            termination,
            diagnostics,
            evaluations,
        })
    }

    pub const fn wall(&self) -> PlanarImpedanceWall {
        self.wall
    }
    pub const fn frequency_hz(&self) -> f64 {
        self.frequency_hz
    }
    pub fn outlets(&self) -> &[CompactWallOutlet] {
        &self.outlets
    }
    pub const fn termination(&self) -> &CoupledTermination {
        &self.termination
    }
    /// Original self impedances and complete free-field mutual terms.
    pub const fn free_space_termination(&self) -> &CoupledTermination {
        &self.free_space
    }
    /// Signed reciprocal correction in row-major order; NOT a standalone load.
    pub fn reflected_impedance(&self) -> &[C64] {
        &self.reflected
    }
    pub const fn evaluations(&self) -> usize {
        self.evaluations
    }
    /// Upper-triangular pair order (0,0),(0,1),...,(1,1),...; unit p/Q evidence.
    pub fn reflected_diagnostics(&self) -> &[ImpedanceWallDiagnostics] {
        &self.diagnostics
    }

    /// Close actual pressure/flow generators against this COMPLETE load through
    /// the original network solver. Outlet labels must be passive ports in that
    /// network; all remaining ports retain explicit source or load conditions.
    /// The caller binds each physical outlet radius/location to its geometry.
    pub fn drive_network(
        &self,
        network: &DuctNetwork,
        gas: &GasState,
        loss: LossModel,
        sources: &[PortSource],
        other_loads: &[PortTermination],
    ) -> Result<NetworkSourceResponse, CompactWallRadiationError> {
        let fluid = gas
            .acoustic_state()
            .map_err(|e| CompactWallRadiationError::Fluid(FluidDuctError::Material(e)))?;
        if fluid != self.wall.fluid() {
            return Err(invalid(
                "wall reaction and network exterior require the same original fluid",
            ));
        }
        network_drive_sources(
            network,
            gas,
            self.termination.omega,
            loss,
            sources,
            other_loads,
            Some(&self.termination),
        )
        .map_err(CompactWallRadiationError::Duct)
    }

    /// Observe the ACTUAL solved outward outlet flows through a matching native
    /// total-field cache. Network currents point inward, hence Q=-U. The source
    /// clock, original wall and ordered locations must match before evaluation.
    /// No power, phase, wave propagation or source amplitudes are reconstructed.
    pub fn observe(
        &self,
        response: &NetworkCoherentResponse,
        observer: &PreparedImpedanceWallCompactSourceObserver,
    ) -> Result<CoherentCompactField, CompactWallRadiationError> {
        let wall = observer.wall();
        if response.omega.to_bits() != self.termination.omega.to_bits()
            || observer.frequency_hz().to_bits() != self.frequency_hz.to_bits()
            || observer.field_kind() != ImpedanceWallFieldKind::Total
            || wall.fluid() != self.wall.fluid()
            || wall.point_m() != self.wall.point_m()
            || wall.outward_normal() != self.wall.outward_normal()
            || wall.impedance_pa_s_m() != self.wall.impedance_pa_s_m()
            || observer.source_count() != self.outlets.len()
            || response.ports.len() > MAX_NETWORK_PORTS
            || response.ports.len() != response.port_inflows.len()
        {
            return Err(invalid(
                "wall observation needs the same solved frequency, fluid, wall and outlets",
            ));
        }
        let mut sources = [CompactSourcePhasor::monopole(C64::ZERO); MAX_NETWORK_PORTS];
        for (i, outlet) in self.outlets.iter().enumerate() {
            if observer.source_position_m(i) != Some(outlet.position_m) {
                return Err(invalid(
                    "wall observer source order differs from the loaded outlets",
                ));
            }
            let mut matching = response
                .ports
                .iter()
                .enumerate()
                .filter(|(_, label)| **label == outlet.label);
            let index = matching
                .next()
                .ok_or_else(|| invalid("loaded outlet missing from network response"))?
                .0;
            if matching.next().is_some() {
                return Err(invalid("network response repeats an outlet label"));
            }
            sources[i].volume_flow_m3_s = -response.port_inflows[index];
        }
        observer
            .evaluate(&sources[..self.outlets.len()])
            .map_err(CompactWallRadiationError::Acoustic)
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), CompactWallRadiationError> {
    if cancelled() {
        Err(CompactWallRadiationError::Cancelled)
    } else {
        Ok(())
    }
}
fn invalid(what: &'static str) -> CompactWallRadiationError {
    CompactWallRadiationError::Invalid(what)
}
fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn distance(a: [f64; 3], b: [f64; 3]) -> f64 {
    (a[0] - b[0]).hypot(a[1] - b[1]).hypot(a[2] - b[2])
}
fn dot_distance(a: [f64; 3], b: [f64; 3], normal: [f64; 3]) -> f64 {
    (0..3)
        .fold(Dd::from_f64(0.0), |sum, i| {
            sum + (Dd::from_f64(a[i]) - Dd::from_f64(b[i])) * Dd::from_f64(normal[i])
        })
        .to_f64()
}

#[cfg(test)]
mod tests;
