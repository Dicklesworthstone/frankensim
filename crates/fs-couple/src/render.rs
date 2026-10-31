//! Block render API (music bead `frankensim-music-v8-root-3ez8g.2.1`).
//!
//! Every performance image eventually runs in a callback-shaped world:
//! AUDIO-RATE BLOCKS advance composed voices; CONTROL-RATE deltas apply
//! BETWEEN blocks (doctrine D17: state lifts at control boundaries so
//! charts move without clicks); no allocation inside a block on the
//! admitted no-alloc voices; cancellation is polled at block boundaries
//! (request -> complete the current block -> stop; blocks are
//! complete on success, so a cancelled render is always a whole number of
//! blocks and resumes bitwise-identically). A physics refusal may leave
//! partial state and permanently poisons the context; it cannot be resumed.
//!
//! This module HOSTS the existing steppers — the exact-FIR characteristic
//! line, the scalar reed islands, the exact-ZOH modal runtime — it never
//! adds a fourth integration scheme (consolidation ruling: exactly three
//! steppers exist by design). Determinism contract: same construction
//! inputs + same control schedule => bitwise-identical pascals regardless
//! of block-size choice, because a block boundary performs no arithmetic
//! (the per-sample operation sequence is invariant; only the loop
//! partition moves).
//!
//! Allocation honesty (measured, not aspirational): the massless-reed
//! voice with an empty plate bank is allocation-free per block and is
//! gated by the counting-allocator test. The massive-reed lay path calls
//! `fs_dcontact::Obstacle::dissipative_modal_forces`, which returns a
//! `Vec` per sample, and the pHS plate path allocates inside
//! `fs_phs::step` — both are DISCLOSED allocating voices (fusion
//! candidates for bead 3ez8g.15), not silently admitted to the no-alloc
//! set.

/// Sample-accurate, bounded control schedules over the existing voice steppers.
pub mod schedule;
/// Geometry-derived compact plate voices with physical force controls.
pub mod plate;
/// Realized resolution and explicit solver-tolerance study seams.
pub mod ladder;
mod scattering;
pub use plate::{CompactPlateVoice, PlateVoiceConfig};
use schedule::force::coupled::{ModalCouplingError, render::{CoupledModalVoice, contact::ContactModalVoice}};
use schedule::force::coupled::render::contact::multiple::MultiContactModalVoice;
use schedule::force::coupled::render::friction::FrictionModalVoice;
use schedule::force::coupled::contact::friction::{FrictionDrive, ModalFrictionError};

use crate::acoustic_realize::AcousticRealizeError;
use crate::bowed_string::runtime::schedule::{BowedScheduleError, ScheduledBowedRenderer};
use crate::modal_acoustic_time::{
    ModalAcousticTimeError, ModalAcousticTimeModel, ModalAcousticWorkspace,
};
use crate::reed_bore::{
    FastSolveStats, ReedSolverMode, blowing_envelope, reed_structural, solve_reed_wave,
};
use crate::thin_plate::PlateBank;
use fs_dcontact::Obstacle;
use fs_duct::{Duct, MouthLoad, Termination};
use fs_exec::CancelGate;
use fs_material::gas::GasState;
use fs_scenario::BeatingReed;

#[cfg(test)]
mod gas_field_tests;

#[cfg(test)]
mod inlet_geometry_tests {
    use super::*;
    use fs_duct::Segment;
    use fs_material::gas::GasSpec;

    /// G0: the source impedance uses the local inlet area, not the downstream
    /// area of a tapered segment. Cylinders preserve their previous value.
    #[test]
    fn reed_source_uses_the_upstream_cross_section() {
        let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap();
        for segment in [
            Segment::Cylinder {
                radius: 0.002,
                length: 0.5,
            },
            Segment::Cone {
                inlet_radius: 0.002,
                outlet_radius: 0.0015,
                length: 0.5,
            },
            Segment::Cone {
                inlet_radius: 0.002,
                outlet_radius: 0.0022,
                length: 0.5,
            },
        ] {
            assert_eq!(segment.inlet_radius(), 0.002);
            let mut voice = ReedBoreVoice::new(
                &Duct {
                    segments: vec![segment],
                },
                &gas,
                BeatingReed {
                    rest_opening_m: 4e-4,
                    width_m: 0.013,
                    closing_pressure_pa: 6000.0,
                    blowing_pressure_pa: 2800.0,
                    attack_s: 0.008,
                    mass_kg: 0.0,
                    stiffness_n_m: 0.0,
                    damping_ratio: 0.35,
                },
                Termination::UnflangedOpen,
                PlateBank::default(),
                1.0,
                48_000,
                4096,
                None,
            )
            .unwrap();
            let area = core::f64::consts::PI * 0.002 * 0.002;
            assert_eq!(voice.area_bore, area);
            assert_eq!(voice.zc, gas.density * gas.sound_speed / area);
            let mut pressure = [0.0; 256];
            voice.step_block(&mut pressure).unwrap();
            assert!(pressure.iter().all(|p| p.is_finite()));
            assert!(pressure.iter().any(|p| p.abs() > 0.0));
        }
    }

    /// G0: a one-edge network bore drives the reed like the same duct: the
    /// active port's area sets the source, and the network line sits on the
    /// duct line's ladder rung. The 2 mm mouth keeps the compact radiation
    /// load at ka <= 1 up to the 24 kHz Nyquist bin.
    #[test]
    fn reed_voice_on_a_one_edge_network_matches_the_duct_voice() {
        use crate::driving_point::LineLadder;
        use crate::driving_point::network::NetworkDrivingPoint;
        use fs_duct::network::{
            DuctNetwork, JunctionCorrections, NetworkEdge, NetworkLoad, NetworkNode, NodeBoundary,
        };
        use fs_ports::{CoordinateBinding, PortKind, PortOrientation, PortTimestamp, StableId};
        let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap();
        let duct = Duct {
            segments: vec![Segment::Cylinder {
                radius: 0.002,
                length: 0.40,
            }],
        };
        let id = |s: &str| StableId::new(s).unwrap();
        let schema = PortKind::FluidPressureFlux
            .scalar_seed_schema(
                id("reed-network-test"),
                CoordinateBinding::new(
                    id("pressure-volume-flow"),
                    id("mouthpiece-plane"),
                    PortOrientation::OutwardFromOwner,
                ),
                PortTimestamp::new(id("declared-frequency-sweep"), 0),
            )
            .unwrap();
        let network = DuctNetwork {
            nodes: vec![
                NetworkNode {
                    id: "in".into(),
                    boundary: NodeBoundary::Port { label: "in".into() },
                },
                NetworkNode {
                    id: "bell".into(),
                    boundary: NodeBoundary::Termination(NetworkLoad::Radiation {
                        radius: 0.002,
                        termination: Termination::UnflangedOpen,
                    }),
                },
            ],
            edges: vec![NetworkEdge {
                id: "bore".into(),
                from: 0,
                to: 1,
                duct: duct.clone(),
                corrections: JunctionCorrections::Disabled,
            }],
        };
        let point = NetworkDrivingPoint::new(network, "in".into(), Vec::new(), schema).unwrap();
        let reed = BeatingReed {
            rest_opening_m: 4e-4,
            width_m: 0.013,
            closing_pressure_pa: 6000.0,
            blowing_pressure_pa: 2800.0,
            attack_s: 0.008,
            mass_kg: 0.0,
            stiffness_n_m: 0.0,
            damping_ratio: 0.35,
        };
        let ladder = LineLadder::default();
        let mut by_duct = ReedBoreVoice::new_with_ladder(
            &duct,
            &gas,
            reed,
            Termination::UnflangedOpen,
            PlateBank::default(),
            1.0,
            48_000,
            65_536,
            None,
            Some(&ladder),
        )
        .unwrap();
        let mut by_network = ReedBoreVoice::new_with_network(
            &point,
            &gas,
            fs_duct::LossModel::Bessel,
            reed,
            PlateBank::default(),
            1.0,
            48_000,
            &ladder,
            None,
        )
        .unwrap();
        assert_eq!(by_network.zc.to_bits(), by_duct.zc.to_bits());
        assert_eq!(
            by_network.line_realization.n_fft,
            by_duct.line_realization.n_fft
        );
        let (mut a, mut b) = (vec![0.0; 4800], vec![0.0; 4800]);
        by_duct.step_block(&mut a).unwrap();
        by_network.step_block(&mut b).unwrap();
        let peak = a.iter().fold(0.0_f64, |m, p| m.max(p.abs()));
        let error = a
            .iter()
            .zip(&b)
            .fold(0.0_f64, |m, (x, y)| m.max((x - y).abs()))
            / peak;
        eprintln!(
            "{{\"suite\":\"fs-couple-render\",\"case\":\"reed-network-one-edge\",\"n_fft\":{},\"max_relative_difference\":{error:.3e}}}",
            by_duct.line_realization.n_fft
        );
        assert!(
            peak > 0.0 && error <= 1.0e-6,
            "network voice differs from the duct voice by {error:.3e}"
        );
    }
}

/// A reed-lay contact that replaces the default slit lay of a massive reed:
/// the power law `stiffness * p^alpha` [N] in the penetration `p = onset_m - y`
/// of the opening `y`, with Hunt-Crossley loss `loss_s_per_m` [s/m]. A real
/// mouthpiece lay is curved, so the reed meets it at the positive opening
/// `onset_m` and the contact force is already carrying load when the slit shuts.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct ReedLay {
    /// Power-law stiffness [N/m^alpha].
    pub(crate) stiffness: f64,
    /// Power-law exponent, at least 1.
    pub(crate) alpha: f64,
    /// Opening [m] below which the reed touches the lay.
    pub(crate) onset_m: f64,
    /// Hunt-Crossley internal loss [s/m].
    pub(crate) loss_s_per_m: f64,
}

/// The two parts of one [`ReedBoreVoice`] observer sample. The voice writes
/// `internal_pa + plate + jet` as one number; an observer law must act on the
/// radiated part only, because the bore pressure travels no distance.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct ReedObserverParts {
    /// Bore pressure at the reed: an internal observation, never propagated.
    pub internal_pa: f64,
    /// Plate radiation plus the jet monopole, both evaluated at `listener_m`.
    pub radiated_pa: f64,
}

/// One reed-on-a-characteristic-line voice: the `realize_reed_bore`
/// physics, restructured so a caller can advance it block by block. The
/// one-shot realizer is now a thin wrapper over this type, so the two
/// paths cannot drift.
pub struct ReedBoreVoice {
    line: fs_vfit::discretize::DelayedFilter,
    line_realization: crate::driving_point::LineRealization,
    molecular_driving_point: Option<Box<crate::driving_point::molecular::MolecularDrivingPoint>>,
    register_shunt: Option<crate::driving_point::regional::outlet::RegisterHoleShuntLine>,
    plates: PlateBank,
    reed: BeatingReed,
    lay: Option<Obstacle>,
    inlet_gas: GasState,
    viscous_slit: Option<crate::bernoulli_aperture::viscous::ViscousSlit>,
    rho: f64,
    sound_speed: f64,
    zc: f64,
    area_bore: f64,
    listener_m: f64,
    dt: f64,
    reed_y: f64,
    reed_v: f64,
    closing_force_n: f64,
    p_plus_prev: f64,
    f_jet_prev: f64,
    last_parts: ReedObserverParts,
    sample_index: u64,
    /// Aperture-junction solver mode. `Strict` is the certification
    /// default; `FastNewton` is an explicit opt-in declared fast mode
    /// (bead frankensim-2s4i5).
    solver_mode: ReedSolverMode,
    /// Fast-mode fallback counters; idle while in Strict mode.
    fast_stats: FastSolveStats,
    junction_step_tolerance: f64,
    junction_residual_max: Option<f64>,
    ports: crate::port_adoption::ModulePorts,
}

impl ReedBoreVoice {
    /// Build the voice exactly as `realize_reed_bore` builds its state
    /// (same admission checks, same line realization, same lay matching,
    /// same `p_plus = 5.0` priming push), then hand stepping to the
    /// caller.
    ///
    /// # Errors
    /// Domain, TMM, or solver refusals — identical to the one-shot path.
    #[allow(clippy::too_many_arguments)] // mirrors the realizer signature
    pub fn new<'a>(
        physics: &Duct,
        gas: &GasState,
        reed: BeatingReed,
        termination: impl Into<fs_duct::MouthLoad<'a>>,
        plates: PlateBank,
        listener_m: f64,
        sample_rate_hz: u32,
        line_samples: usize,
        wall: Option<&fs_phs::WallPin>,
    ) -> Result<Self, AcousticRealizeError> {
        Self::new_with_ladder(
            physics,
            gas,
            reed,
            termination,
            plates,
            listener_m,
            sample_rate_hz,
            line_samples,
            wall,
            None,
        )
    }

    /// Build a reed-bore voice with an explicit characteristic-line FFT ladder.
    ///
    /// # Errors
    /// Domain, TMM, ladder-convergence, or solver refusals.
    #[allow(clippy::too_many_arguments)] // mirrors the realizer signature
    pub fn new_with_ladder<'a>(
        physics: &Duct,
        gas: &GasState,
        reed: BeatingReed,
        termination: impl Into<fs_duct::MouthLoad<'a>>,
        plates: PlateBank,
        listener_m: f64,
        sample_rate_hz: u32,
        line_samples: usize,
        wall: Option<&fs_phs::WallPin>,
        ladder: Option<&crate::driving_point::LineLadder>,
    ) -> Result<Self, AcousticRealizeError> {
        let termination = termination.into();
        Self::with_line(physics, gas, reed, plates, listener_m, sample_rate_hz, |zc| {
            crate::driving_point::characteristic_line_realized_with_ladder(
                physics,
                gas,
                termination,
                sample_rate_hz,
                line_samples,
                zc,
                wall,
                None,
                ladder,
            )
            .map_err(crate::acoustic_realize::map_drive)
        })
    }

    /// Build the same reed mechanics with one fixed gas for each duct record.
    ///
    /// The inlet medium supplies the Bernoulli law and existing compact-jet
    /// proxy; each downstream region supplies its own propagation and losses.
    /// This is not a calibrated exterior observer or an evolving composition.
    /// All mechanical initialization and per-sample stepping are shared with
    /// [`Self::new`]. No uniform trial duct is constructed to admit a regional one.
    ///
    /// # Errors
    /// Original reed/line refusals, incomplete or unequal-pressure gas binding,
    /// or unsupported cross-region mutual-hole interactions.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_segment_gases<'a>(
        physics: &Duct, gases: &[GasState], reed: BeatingReed,
        termination: impl Into<fs_duct::MouthLoad<'a>>, plates: PlateBank, listener_m: f64,
        sample_rate_hz: u32, line_samples: usize, wall: Option<&fs_phs::WallPin>,
    ) -> Result<Self, AcousticRealizeError> {
        let termination = termination.into();
        let regional = crate::driving_point::regional::RegionalDuct::new(physics, gases)
            .map_err(|e| crate::acoustic_realize::map_drive(crate::driving_point::DrivingPointError::Duct(e)))?;
        Self::new_with_regional(&regional, reed, termination, plates, listener_m,
            sample_rate_hz, line_samples, wall)
    }

    /// Bind the original reed to an already-admitted regional duct, including
    /// its independent exterior radiation load. Source area, source gas and
    /// propagation all come from that ONE geometry/material binding.
    ///
    /// A physical receiver consumes [`Self::outgoing_pressure_pa`] after each
    /// accepted sample. The voice's ordinary output remains its legacy proxy.
    /// This uses the regional owner's existing default convergence ladder.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_regional<'a>(
        regional: &crate::driving_point::regional::RegionalDuct,
        reed: BeatingReed,
        termination: impl Into<fs_duct::MouthLoad<'a>>,
        plates: PlateBank,
        listener_m: f64,
        sample_rate_hz: u32,
        line_samples: usize,
        wall: Option<&fs_phs::WallPin>,
    ) -> Result<Self, AcousticRealizeError> {
        Self::new_with_regional_wall(regional, reed, termination, plates,
            listener_m, sample_rate_hz, line_samples, fs_duct::WallModel::from_pin(wall))
    }

    /// Bind a shared local wall law to the original regional reed feedback.
    ///
    /// A source-resolved material wall supplies the same geometry-dependent
    /// impedance used by the physical outlet observer. Elastic, viscous and
    /// hereditary material response remain with the wall owner; neither the
    /// reed nor this adapter copies coefficients into a second wall model.
    /// Every sampled radius and frequency must satisfy the original source
    /// domain. This frozen membrane image does not add exterior wall radiation.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_regional_wall<'a>(
        regional: &crate::driving_point::regional::RegionalDuct,
        reed: BeatingReed,
        termination: impl Into<fs_duct::MouthLoad<'a>>,
        plates: PlateBank,
        listener_m: f64,
        sample_rate_hz: u32,
        line_samples: usize,
        wall: fs_duct::WallModel<'_>,
    ) -> Result<Self, AcousticRealizeError> {
        let termination = termination.into();
        Self::with_inlet_line(|| Ok(regional.inlet_radius_m()), regional.inlet_gas(),
            reed, plates, listener_m, sample_rate_hz, |_| {
                regional.characteristic_line_with_wall_realized(termination, sample_rate_hz,
                    line_samples, wall, None).map_err(crate::acoustic_realize::map_drive)
            })
    }

    /// Render the same reed mechanics against a prescribed continuous
    /// temperature/composition field. Local NASA-9 and material transport
    /// evaluations feed the existing regional characteristic line; the exact
    /// inlet state supplies source density and characteristic impedance.
    /// `cells_per_segment` is the explicit midpoint field-refinement axis.
    ///
    /// This keeps the existing reed and compact-jet observer model. The field
    /// is fixed during rendering and does not simulate breath heat transfer,
    /// species advection, a player or a calibrated complete instrument.
    ///
    /// # Errors
    /// Existing field/source, geometry/subdivision, reed or line refusals.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_gas_field<'a>(
        physics: &Duct,
        field: &fs_material::gas_field::AxialGasField,
        cells_per_segment: usize,
        reed: BeatingReed,
        termination: impl Into<MouthLoad<'a>>,
        plates: PlateBank,
        listener_m: f64,
        sample_rate_hz: u32,
        line_samples: usize,
        wall: Option<&fs_phs::WallPin>,
    ) -> Result<Self, AcousticRealizeError> {
        use crate::driving_point::regional::{RegionalDuct, field::GasFieldDuctError};
        let termination = termination.into();
        let regional = RegionalDuct::from_gas_field(physics, field, cells_per_segment)
            .map_err(|error| match error {
                GasFieldDuctError::Material(error) => AcousticRealizeError::Ambient(error.to_string()),
                GasFieldDuctError::Duct(error) => AcousticRealizeError::Duct(error),
            })?;
        Self::with_line(physics, regional.inlet_gas(), reed, plates, listener_m, sample_rate_hz, |_| {
            regional.characteristic_line_realized(termination, sample_rate_hz, line_samples, wall, None)
                .map_err(crate::acoustic_realize::map_drive)
        })
    }

    /// Render against a continuous primitive gas field and smooth radius
    /// geometry using the same Magnus and stable two-port propagation owners.
    /// The exact inlet medium sets reed flow and the characteristic source
    /// port; local gas properties set propagation/loss and the exact outlet
    /// medium sets terminal loading. The supplied integration axis is refinable.
    ///
    /// This is the existing reed and compact-source observer with a frozen
    /// axial field, not a species/thermal transport or complete-player model.
    /// Tone-hole junctions use their local gas. Multiple active holes require
    /// a homogeneous connecting field for the existing mutual-hole law.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_continuous_gas_field<'a>(
        physics: &Duct,
        field: &fs_material::gas_field::AxialGasField,
        steps_per_span: fs_duct::RiccatiSteps,
        reed: BeatingReed,
        termination: impl Into<MouthLoad<'a>>,
        plates: PlateBank,
        listener_m: f64,
        sample_rate_hz: u32,
        line_samples: usize,
        wall: Option<&fs_phs::WallPin>,
    ) -> Result<Self, AcousticRealizeError> {
        use crate::driving_point::regional::{RegionalDuct, field::GasFieldDuctError};
        let termination = termination.into();
        let regional = RegionalDuct::from_continuous_gas_field(physics, field, steps_per_span)
            .map_err(|error| match error {
                GasFieldDuctError::Material(error) => AcousticRealizeError::Ambient(error.to_string()),
                GasFieldDuctError::Duct(error) => AcousticRealizeError::Duct(error),
            })?;
        Self::with_line(physics, regional.inlet_gas(), reed, plates, listener_m,
            sample_rate_hz, |_| {
                regional.characteristic_line_realized(
                    termination, sample_rate_hz, line_samples, wall, None,
                ).map_err(crate::acoustic_realize::map_drive)
            })
    }

    /// Build the same reed mechanics on a duct network's driving point. The
    /// network line ([`crate::driving_point::network::NetworkDrivingPoint::realize_line`])
    /// replaces the single-duct line, and the source impedance uses the active
    /// port's cross-section. Coupled external ports, for example a
    /// multi-aperture radiation matrix, are supplied per sampled frequency by
    /// `coupled_at`. This is the existing reed and compact-jet observer.
    ///
    /// # Errors
    /// Reed or listener refusals, an active port without one incident edge,
    /// and network-line refusals.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_network(
        driving_point: &crate::driving_point::network::NetworkDrivingPoint,
        gas: &GasState,
        loss: fs_duct::LossModel,
        reed: BeatingReed,
        plates: PlateBank,
        listener_m: f64,
        sample_rate_hz: u32,
        ladder: &crate::driving_point::LineLadder,
        coupled_at: Option<crate::driving_point::network::CoupledSource<'_>>,
    ) -> Result<Self, AcousticRealizeError> {
        Self::with_inlet_line(
            || {
                driving_point
                    .active_port_radius()
                    .map_err(AcousticRealizeError::Duct)
            },
            gas,
            reed,
            plates,
            listener_m,
            sample_rate_hz,
            |zc| {
                driving_point
                    .realize_line(gas, loss, sample_rate_hz, zc, ladder, coupled_at)
                    .map_err(crate::acoustic_realize::map_drive)
            },
        )
    }

    /// Build the same pressure/structure junction with the original gas's
    /// molecular reservoirs participating in bore feedback. The prepared
    /// constitutive source and exterior load remain available on the voice;
    /// no observer attenuation is substituted for the interior material law.
    ///
    /// This retains the existing reed mechanics, contact and compact-jet
    /// observation, whose gas is the inlet state. A separately bound exterior
    /// belongs to terminal loading and the physical baffled outlet receiver.
    /// Source refinement and the FIR ladder remain independent
    /// requirements; it is not a new nonlinear-gas or bell-radiation model.
    pub fn new_with_molecular_driving_point(
        driving_point: crate::driving_point::molecular::MolecularDrivingPoint,
        termination: fs_duct::Termination,
        reed: BeatingReed,
        plates: PlateBank,
        listener_m: f64,
        sample_rate_hz: u32,
        ladder: &crate::driving_point::LineLadder,
    ) -> Result<Self, AcousticRealizeError> {
        let mut voice = Self::with_inlet_line(
            || Ok(driving_point.inlet_radius_m()),
            driving_point.inlet_gas(), reed, plates, listener_m, sample_rate_hz,
            |zc| driving_point.realize_line(termination, sample_rate_hz, zc, ladder)
                .map_err(crate::acoustic_realize::map_drive),
        )?;
        voice.molecular_driving_point = Some(Box::new(driving_point));
        Ok(voice)
    }

    /// Source gas, molecular reservoirs, geometry and exterior retained by the
    /// material-backed characteristic path. Legacy paths retain their own laws.
    #[must_use]
    pub fn molecular_driving_point(
        &self,
    ) -> Option<&crate::driving_point::molecular::MolecularDrivingPoint> {
        self.molecular_driving_point.as_deref()
    }

    fn with_line(
        physics: &Duct,
        gas: &GasState,
        reed: BeatingReed,
        plates: PlateBank,
        listener_m: f64,
        sample_rate_hz: u32,
        build_line: impl FnOnce(
            f64,
        ) -> Result<
            (fs_vfit::discretize::DelayedFilter, crate::driving_point::LineRealization),
            AcousticRealizeError,
        >,
    ) -> Result<Self, AcousticRealizeError> {
        let inlet_radius = || -> Result<f64, AcousticRealizeError> {
            Ok(physics
                .segments
                .first()
                .ok_or(AcousticRealizeError::InvalidDescription {
                    what: "duct has no segments",
                })?
                .inlet_radius())
        };
        Self::with_inlet_line(
            inlet_radius,
            gas,
            reed,
            plates,
            listener_m,
            sample_rate_hz,
            build_line,
        )
    }

    // One mechanical initial state, lay-contact admission and priming operation
    // for every physical line provider. No provider advances the voice.
    fn with_inlet_line(
        inlet_radius: impl FnOnce() -> Result<f64, AcousticRealizeError>,
        gas: &GasState,
        reed: BeatingReed,
        plates: PlateBank,
        listener_m: f64,
        sample_rate_hz: u32,
        build_line: impl FnOnce(
            f64,
        ) -> Result<
            (fs_vfit::discretize::DelayedFilter, crate::driving_point::LineRealization),
            AcousticRealizeError,
        >,
    ) -> Result<Self, AcousticRealizeError> {
        if !crate::reed_bore::reed_parameters_valid(reed) {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "reed parameters must be physical and finite",
            });
        }
        if !listener_m.is_finite() || listener_m <= 0.0 {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "listener distance must be positive and finite",
            });
        }
        let inlet_r = inlet_radius()?;
        let area_bore = core::f64::consts::PI * inlet_r * inlet_r;
        let zc = gas.density * gas.sound_speed / area_bore;
        let (mut line, line_realization) = build_line(zc)?;
        let lay = if reed.mass_kg > 0.0 {
            let k_lay = 1.0e7 * reed.width_m;
            let (_k, r_damp) = reed_structural(reed);
            let chi = r_damp / (k_lay * reed.rest_opening_m * reed.rest_opening_m).max(1.0e-18);
            Some(
                crate::unilateral_contact::slit_lay(k_lay, 2.0)
                    .and_then(|o| o.with_internal_loss(chi))
                    .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?,
            )
        } else {
            None
        };
        let p_plus_prev = 5.0;
        let _ = line
            .push(p_plus_prev)
            .map_err(|_| AcousticRealizeError::Reed {
                what: "characteristic line left the finite set",
            })?;
        Ok(Self {
            line,
            line_realization,
            molecular_driving_point: None,
            register_shunt: None,
            plates,
            reed,
            lay,
            inlet_gas: *gas,
            viscous_slit: None,
            rho: gas.density,
            sound_speed: gas.sound_speed,
            zc,
            area_bore,
            listener_m,
            dt: 1.0 / f64::from(sample_rate_hz),
            reed_y: reed.rest_opening_m,
            reed_v: 0.0,
            closing_force_n: 0.0,
            p_plus_prev,
            f_jet_prev: 0.0,
            last_parts: ReedObserverParts::default(),
            sample_index: 0,
            solver_mode: ReedSolverMode::Strict,
            fast_stats: FastSolveStats::default(),
            junction_step_tolerance: 1.0e-9,
            junction_residual_max: None,
            ports: crate::port_adoption::ModulePorts::new("reed-bore-line", vec![
                ("inlet", fs_ports::PortKind::FluidPressureFlux),
            ]),
        })
    }

    /// Total samples rendered since construction.
    #[must_use]
    pub const fn samples_rendered(&self) -> u64 {
        self.sample_index
    }

    /// Last accepted outgoing characteristic pressure [Pa] at the original
    /// inlet plane. Immediately after construction this is the existing priming
    /// sample; a downstream observer must consume it once before sample zero.
    #[must_use]
    pub const fn outgoing_pressure_pa(&self) -> f64 { self.p_plus_prev }

    /// Characteristic-line realization receipt (`n_fft`, rungs tried, tail energy,
    /// peak shift, and explicit pole states).
    #[must_use]
    pub const fn line_realization(&self) -> &crate::driving_point::LineRealization {
        &self.line_realization
    }

    /// Mutable access to the plate bank (the one-shot realizer swaps the
    /// caller's bank in and out through this seam so its error semantics
    /// stay byte-identical to the pre-render-API code).
    pub const fn plate_bank_mut(&mut self) -> &mut PlateBank {
        &mut self.plates
    }

    pub(crate) fn set_register_shunt(
        &mut self,
        shunt: crate::driving_point::regional::outlet::RegisterHoleShuntLine,
    ) {
        self.register_shunt = Some(shunt);
    }

    /// Select the shared viscous channel law from physical streamwise length
    /// [m] and this voice's exact inlet gas. Density and viscosity therefore
    /// remain bound to the SAME regional or axial material state as the source
    /// impedance. This law enters every trial contact/flow junction.
    ///
    /// Omitting this selection preserves the inviscid source. The selected
    /// model is fully developed, quasisteady parallel-plate flow plus separated
    /// jet loss; it does not resolve unsteady channel inertia or turbulence.
    /// Admission is allowed only before the first rendered sample.
    pub fn set_viscous_slit_length(&mut self, length_m: f64) -> Result<(), AcousticRealizeError> {
        if self.sample_index != 0 {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "viscous reed channel selection requires an unadvanced voice",
            });
        }
        let slit = crate::bernoulli_aperture::viscous::ViscousSlit::from_gas(&self.inlet_gas, length_m)?;
        self.viscous_slit = Some(slit);
        Ok(())
    }

    /// Selected material-gas channel and physical length, if explicitly supplied.
    #[must_use]
    pub const fn viscous_slit(&self) -> Option<&crate::bernoulli_aperture::viscous::ViscousSlit> {
        self.viscous_slit.as_ref()
    }

    /// Replace the default slit lay by a measured one ([`ReedLay`]). The voice
    /// must not have advanced, the reed must be massive, and the lay must start
    /// below the rest opening so the reed's initial state is out of contact.
    ///
    /// # Errors
    /// A massless reed, an advanced voice, an onset at or above the rest
    /// opening, or contact admission (negative stiffness or loss, exponent < 1).
    pub(crate) fn set_lay(&mut self, lay: ReedLay) -> Result<(), AcousticRealizeError> {
        if self.reed.mass_kg <= 0.0 || self.sample_index != 0 {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "a measured lay needs a massive reed and a voice that has not advanced",
            });
        }
        if lay.onset_m.is_nan() || lay.onset_m >= self.reed.rest_opening_m {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "the lay onset must lie below the reed's rest opening",
            });
        }
        let obstacle =
            crate::unilateral_contact::slit_lay_with_onset(lay.stiffness, lay.alpha, lay.onset_m)
                .and_then(|o| o.with_internal_loss(lay.loss_s_per_m))
                .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
        self.lay = Some(obstacle);
        if self.closing_force_n > 0.0 {
            self.reed_y = self.static_equilibrium_opening();
            self.reed_v = 0.0;
        }
        Ok(())
    }

    fn static_equilibrium_opening(&self) -> f64 {
        let face = crate::reed_bore::reed_pressure_face(self.reed);
        let (stiffness, _damping) = reed_structural(self.reed);
        let downward_n = face * self.reed.blowing_pressure_pa + self.closing_force_n;
        if self.reed.mass_kg <= 0.0
            || !downward_n.is_finite()
            || downward_n <= 0.0
            || !stiffness.is_finite()
            || stiffness <= 0.0
        {
            return self.reed.rest_opening_m;
        }
        let y_free = self.reed.rest_opening_m - downward_n / stiffness;
        let Some(lay) = self.lay.as_ref() else {
            return y_free;
        };
        let contact_at = |y: f64| -> f64 {
            crate::unilateral_contact::distributed::ApertureContactStep::new(lay, y)
                .and_then(|step| step.response(y, 0.0))
                .map_or(0.0, |resp| resp.force)
        };
        if contact_at(y_free) <= 0.0 {
            return y_free;
        }
        let mut lo = y_free;
        let mut hi = self.reed.rest_opening_m;
        for _ in 0..64 {
            let mid = f64::midpoint(lo, hi);
            let residual = stiffness * (mid - self.reed.rest_opening_m) + downward_n - contact_at(mid);
            if residual < 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        f64::midpoint(lo, hi)
    }

    /// Replace the reed's held blowing pressure. Applied between blocks
    /// only (the render context enforces the boundary); the change is a
    /// pure input-parameter move — no stored state depends on it, so the
    /// D17 lift map is empty and the transition is click-bounded by the
    /// reed's own dynamics.
    pub fn set_blowing_pressure(&mut self, pressure_pa: f64) {
        self.reed.blowing_pressure_pa = pressure_pa;
        if self.sample_index == 0 && self.closing_force_n > 0.0 {
            self.reed_y = self.static_equilibrium_opening();
            self.reed_v = 0.0;
        }
    }

    /// Held generalized closing force [N] on the reed's opening coordinate.
    #[must_use]
    pub const fn closing_force_n(&self) -> f64 {
        self.closing_force_n
    }

    /// Replace the held generalized closing force [N] on a massive reed's
    /// opening coordinate (such as tongue-reed contact prior to and during a
    /// tongued release). At sample zero, a positive closing force initializes
    /// the reed at its static contact equilibrium so pre-existing tongue
    /// contact does not inject a synthetic slap transient.
    pub fn set_closing_force(&mut self, force_n: f64) {
        self.closing_force_n = force_n;
        if self.sample_index == 0 && self.closing_force_n > 0.0 {
            self.reed_y = self.static_equilibrium_opening();
            self.reed_v = 0.0;
        }
    }

    /// Opt in to the declared fast-mode aperture solver
    /// ([`ReedSolverMode::FastNewton`], bead frankensim-2s4i5). The
    /// strict bisection stays the certification default and the
    /// constructor's value; this is an explicit, per-voice upgrade.
    pub const fn set_solver_mode(&mut self, mode: ReedSolverMode) {
        self.solver_mode = mode;
    }

    /// Fast-mode fallback counters. Zero `newton + fallback` means no
    /// solver-work sample was seen (strict mode or closed-branch only).
    #[must_use]
    pub const fn fast_solver_stats(&self) -> FastSolveStats {
        self.fast_stats
    }

    /// The internal and radiated parts of the last sample `step_block` wrote.
    /// Their sum in the order `(internal + plate) + jet` is that sample's
    /// bits; `radiated_pa` is `plate + jet`, so a caller that recombines them
    /// rounds once more.
    #[must_use]
    pub const fn last_observer_parts(&self) -> ReedObserverParts {
        self.last_parts
    }

    /// Advance exactly `out.len()` samples, writing observer pascals.
    ///
    /// Source pressure, flow and reflected wave use the same sample of the
    /// original line, including instantaneous reflection. The one-shot path
    /// shares this loop; a block boundary performs no arithmetic.
    ///
    /// # Errors
    /// Solver or finite-set refusals; on refusal the voice must be
    /// discarded (mid-sample state is not rewound).
    pub fn step_block(&mut self, out: &mut [f64]) -> Result<(), AcousticRealizeError> {
        for slot in out.iter_mut() {
            let t = self.sample_index as f64 * self.dt;
            let p_m = blowing_envelope(self.reed, t);
            let mut boundary = self.line.boundary_response().map_err(|_|
                AcousticRealizeError::Reed { what: "characteristic boundary left the finite set" })?;
            let shunt_incoming = self
                .register_shunt
                .as_ref()
                .map_or(
                    0.0,
                    crate::driving_point::regional::outlet::RegisterHoleShuntLine::s12_history,
                );
            boundary.history += shunt_incoming;
            let u_body = self.plates.volume_velocity();
            let p_plus = if self.reed.mass_kg > 0.0 {
                let (pp, y, v) = self.step_massive_scattering(boundary, p_m, u_body)?;
                self.reed_y = y;
                self.reed_v = v;
                pp
            } else if let Some(slit) = self.viscous_slit {
                let outgoing = crate::reed_bore::solve_reed_wave_viscous(
                    self.reed, self.rho, self.zc, boundary.direct, boundary.history,
                    p_m, self.p_plus_prev, u_body, slit,
                )?;
                if self.solver_mode == ReedSolverMode::FastNewton {
                    // The old Newton Jacobian describes inviscid flow only.
                    // Count this explicit handoff to the shared strict solver.
                    self.fast_stats.fallback_samples += 1;
                }
                outgoing
            } else {
                match self.solver_mode {
                    ReedSolverMode::Strict => solve_reed_wave(
                        self.reed,
                        self.rho,
                        self.zc,
                        boundary.direct,
                        boundary.history,
                        p_m,
                        self.p_plus_prev,
                        u_body,
                    )?,
                    ReedSolverMode::FastNewton => crate::reed_bore::solve_reed_wave_fast_tolerance(
                        self.reed,
                        self.rho,
                        self.zc,
                        boundary.direct,
                        boundary.history,
                        p_m,
                        self.p_plus_prev,
                        u_body,
                        &mut self.fast_stats,
                        self.junction_step_tolerance,
                    )?,
                }
            };
            self.p_plus_prev = p_plus;
            if self.reed.mass_kg == 0.0 {
                let residual = match self.viscous_slit {
                    Some(slit) => crate::reed_bore::reed_viscous_flow_mismatch(
                        self.reed, self.rho, self.zc, boundary.direct, boundary.history,
                        p_m, p_plus, u_body, slit,
                    )?,
                    None => crate::reed_bore::reed_flow_mismatch(
                        self.reed, self.rho, self.zc, boundary.direct, boundary.history,
                        p_m, p_plus, u_body,
                    ),
                }.abs();
                self.junction_residual_max = Some(self.junction_residual_max.unwrap_or(0.0).max(residual));
            }
            let p_minus_now = self
                .line
                .push(p_plus)
                .map_err(|_| AcousticRealizeError::Reed {
                    what: "characteristic line left the finite set",
                })?
                + shunt_incoming;
            if let Some(shunt) = &mut self.register_shunt {
                shunt.step(p_plus).map_err(|_| AcousticRealizeError::Reed {
                    what: "register-hole orifice shunt left the finite set",
                })?;
            }
            let p_bore = p_plus + p_minus_now;
            let plate_pa = self.plates.drive_and_radiate(
                p_bore * self.area_bore,
                self.dt,
                self.rho,
                self.listener_m,
            )?;
            let mut p_obs = p_bore;
            p_obs += plate_pa;
            let opening = if self.reed.mass_kg > 0.0 {
                self.reed_y.max(0.0)
            } else {
                crate::reed_bore::aperture_of(self.reed)
                    .opening_m(p_m - p_bore)
                    .max(0.0)
            };
            let f_jet = (p_m - p_bore) * self.reed.width_m * opening;
            let dfdt = (f_jet - self.f_jet_prev) / self.dt;
            self.f_jet_prev = f_jet;
            let jet_pa = self.rho * dfdt
                / (4.0 * core::f64::consts::PI * self.listener_m * self.sound_speed);
            p_obs += jet_pa;
            *slot = p_obs;
            self.last_parts = ReedObserverParts { internal_pa: p_bore, radiated_pa: plate_pa + jet_pa };
            self.sample_index += 1;
            // Current accepted incoming wave, matching the affine boundary used
            // by the source up to convolution summation roundoff. This is port
            // work, not an inferred line storage or filter-dissipation account.
            let volume = self.dt*(p_plus-p_minus_now)/self.zc;
            self.ports.accept(crate::port_adoption::StepAccount {
                tick:self.sample_index,duration_s:self.dt,storage_change_j:None,
                loss_j:None,included_work_j:0.0,entropy_change_j_k:0.0,tolerance_j:None,
                unavailable:Some("ReedBoreVoice/DelayedFilter has no physical line storage/loss owner or energy tolerance [J]"),
            },[crate::port_adoption::PortTransfer::work(p_bore*volume,true)]);
        }
        Ok(())
    }
}

crate::port_adoption::impl_ported_module!(ReedBoreVoice, ports);

/// A modal string voice hosting the exact-ZOH runtime: per-sample free
/// decay (or a held generalized force) with the model's own observer
/// pressure. Candidate states and energy/loss diagnostics reuse a workspace
/// allocated at construction; the callback does not construct per-sample
/// vectors. All exact-ZOH budget and passivity gates remain active.
pub struct ModalStringVoice {
    model: ModalAcousticTimeModel,
    workspace: ModalAcousticWorkspace,
    held_force: Vec<f64>,
    sample_index: u64,
}

impl ModalStringVoice {
    /// Wrap an admitted model with a held (zero for free decay)
    /// generalized-force vector.
    ///
    /// # Errors
    /// Force-length mismatch against the model's mode count.
    pub fn new(
        model: ModalAcousticTimeModel,
        held_force: Vec<f64>,
    ) -> Result<Self, ModalAcousticTimeError> {
        if held_force.len() != model.modes().len() {
            return Err(ModalAcousticTimeError::ForceCountMismatch {
                expected: model.modes().len(),
                found: held_force.len(),
            });
        }
        let workspace = ModalAcousticWorkspace::new(&model);
        Ok(Self {
            model,
            workspace,
            held_force,
            sample_index: 0,
        })
    }

    /// Total samples rendered since construction.
    #[must_use]
    pub const fn samples_rendered(&self) -> u64 {
        self.sample_index
    }

    /// Advance exactly `out.len()` samples of observer pressure.
    ///
    /// # Errors
    /// Model refusals (budget ceilings, non-finite states).
    pub fn step_block(&mut self, out: &mut [f64]) -> Result<(), ModalAcousticTimeError> {
        for slot in out.iter_mut() {
            let frame = self
                .model
                .step_into(&self.held_force, &mut self.workspace)?;
            *slot = frame.observer_pressure_pa;
            self.sample_index += 1;
        }
        Ok(())
    }
}

/// Typed control deltas applied BETWEEN blocks (control rate). v1 carries
/// the minimal set the current voices support; gesture schedules
/// (bead 3ez8g.2.3) lower into these, and image-swap lifts land with
/// their track beads.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ControlDelta {
    /// Held physical tangential speed and compressive load on a friction port.
    SetFrictionDrive {
        /// Friction-network voice slot.
        voice: usize,
        /// Signed surface speed [m/s].
        speed_m_s: f64,
        /// Nonnegative prescribed normal load [N]; zero releases the port.
        normal_force_n: f64,
    },
    /// Total signed force [N] on a plate reduction's declared footprint.
    SetPlateForce {
        /// Plate voice slot.
        voice: usize,
        /// Zero releases the load while preserving vibration.
        force_n: f64,
    },
    /// Replace one retained mode's held generalized force without changing state.
    /// Zero releases this input; existing vibration continues to decay.
    SetModalForce {
        /// Voice index in the context.
        voice: usize,
        /// Zero-based index in the voice's admitted modal basis.
        mode: usize,
        /// Signed mass-normalized generalized force [N/sqrt(kg)].
        force_n_per_sqrt_kg: f64,
    },
    /// Replace a reed voice's held blowing pressure [Pa].
    SetBlowingPressure {
        /// Voice index in the context.
        voice: usize,
        /// New held pressure [Pa].
        pressure_pa: f64,
    },
    /// Replace a massive reed voice's held generalized closing force [N].
    SetReedClosingForce {
        /// Voice index in the context.
        voice: usize,
        /// Nonnegative generalized closing force [N].
        force_n: f64,
    },
}

/// One applied-control record: what changed at which block boundary, and
/// what state lift it required (empty = pure input-parameter move). This
/// is the D17 lift map made inspectable.
#[derive(Clone, Debug, PartialEq)]
pub struct ControlRecord {
    /// Block index the delta applied BEFORE.
    pub block_index: u64,
    /// The delta, verbatim.
    pub delta: ControlDelta,
    /// Human-readable lift description (empty when no state moved).
    pub lift: &'static str,
}

/// Typed refusal from the render context.
#[derive(Debug)]
pub enum RenderError {
    /// The implicit friction/network owner refused its physical sample.
    Friction(ModalFrictionError),
    /// A previous voice refusal left partial state; construct a new context.
    Poisoned,
    /// A voice refused mid-block; the context is poisoned.
    Voice(AcousticRealizeError),
    /// The modal voice refused mid-block; the context is poisoned.
    Modal(ModalAcousticTimeError),
    /// A retained bow performance refused; its incomplete callback is discarded.
    Bowed(BowedScheduleError),
    /// A two-way coupled modal network refused its complete sample transaction.
    Coupled(ModalCouplingError),
    /// A control named a voice that does not exist or cannot accept it.
    Control {
        /// What was wrong.
        what: &'static str,
    },
    /// A control named a voice slot beyond the admitted voices.
    UnknownVoice {
        /// The zero-based slot the control named.
        index: usize,
    },
    /// A caller-supplied size exceeds a documented pre-sized budget (block
    /// length versus construction-time scratch, or a block-count product).
    /// This is an admission failure of the request shape, not a control
    /// defect.
    Sizing {
        /// What was wrong.
        what: &'static str,
    },
    /// `block` was called with an empty output slice.
    EmptyBlock,
}

impl core::fmt::Display for RenderError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Friction(e) => write!(f, "friction voice refusal: {e}"),
            Self::Poisoned => write!(
                f,
                "render context has partial state after a voice refusal; construct a new context"
            ),
            Self::Voice(e) => write!(f, "voice refusal: {e:?}"),
            Self::Modal(e) => write!(f, "modal voice refusal: {e:?}"),
            Self::Bowed(e) => write!(f, "bowed voice refusal: {e}"),
            Self::Coupled(e) => write!(f, "coupled modal voice refusal: {e}"),
            Self::Control { what } => write!(f, "control refusal: {what}"),
            Self::UnknownVoice { index } => {
                write!(f, "control names unknown voice slot {index}")
            }
            Self::Sizing { what } => write!(f, "sizing refusal: {what}"),
            Self::EmptyBlock => write!(f, "block output slice is empty"),
        }
    }
}

impl std::error::Error for RenderError {}

/// One voice slot: the admitted performance images this context can host.
#[allow(clippy::large_enum_variant)] // a handful of voice slots; boxing would ripple the API
pub enum RenderVoice {
    /// Two-way mechanical network with an implicit driven friction port.
    FrictionModal(Box<FrictionModalVoice>),
    /// A physical bow performance and its compact-body pressure observer.
    /// Bow controls remain owned by the existing sample-clock scheduler; outer
    /// callbacks and other voices' event splits never reset the string or body.
    BowedString(Box<ScheduledBowedRenderer>),
    /// Reed on a characteristic line (wind-reed filling).
    ReedBore(ReedBoreVoice),
    /// Exact-ZOH modal string (string filling).
    ModalString(ModalStringVoice),
    /// Geometry/material-derived compact plate, using the existing body stepper.
    CompactPlate(CompactPlateVoice),
    /// Modal components exchanging two-way physical connection forces.
    CoupledModal(Box<CoupledModalVoice>),
    /// Modal network with an implicit two-body unilateral contact.
    ContactModal(Box<ContactModalVoice>),
    /// Modal network with a jointly solved set of compliant normal contacts.
    MultiContactModal(Box<MultiContactModalVoice>),
}

/// The block render context: a set of voices summed into one observer
/// channel, advanced a block at a time, with controls applied only at
/// block boundaries and cancellation polled only between blocks.
pub struct RenderContext {
    voices: Vec<RenderVoice>,
    /// Scratch buffer reused across blocks so summation allocates nothing
    /// after construction.
    scratch: Vec<f64>,
    blocks_rendered: u64,
    samples_rendered: u64,
    controls_applied: Vec<ControlRecord>,
    poisoned: bool,
}

impl RenderContext {
    /// Build a context over admitted voices, pre-sizing the scratch buffer
    /// to `max_block` samples (larger blocks refuse rather than allocate).
    #[must_use]
    pub fn new(voices: Vec<RenderVoice>, max_block: usize) -> Self {
        Self {
            voices,
            scratch: vec![0.0; max_block],
            blocks_rendered: 0,
            samples_rendered: 0,
            controls_applied: Vec::new(),
            poisoned: false,
        }
    }

    /// Blocks rendered so far.
    #[must_use]
    pub const fn blocks_rendered(&self) -> u64 {
        self.blocks_rendered
    }

    /// Samples completed since this context was constructed. Unlike the block
    /// counter, this clock is independent of the callback partition.
    #[must_use]
    pub const fn samples_rendered(&self) -> u64 {
        self.samples_rendered
    }

    /// Construction-time scratch capacity in samples.
    #[must_use]
    pub fn max_block_len(&self) -> usize {
        self.scratch.len()
    }

    /// The applied-control log (the inspectable D17 record).
    #[must_use]
    pub fn control_log(&self) -> &[ControlRecord] {
        &self.controls_applied
    }

    /// Inspect a hosted bow's state, pending inputs and finite horizon without
    /// allowing its clock to advance behind the shared renderer.
    #[must_use]
    pub fn bowed_performance(&self, index: usize) -> Option<&ScheduledBowedRenderer> {
        match self.voices.get(index) {
            Some(RenderVoice::BowedString(bow)) => Some(bow.as_ref()),
            _ => None,
        }
    }

    /// Inspect the accepted friction mechanics without bypassing the shared clock.
    #[must_use]
    pub fn friction_voice(&self, index: usize) -> Option<&FrictionModalVoice> {
        match self.voices.get(index) {
            Some(RenderVoice::FrictionModal(voice)) => Some(voice),
            _ => None,
        }
    }

    /// Inspect a hosted reed-bore voice without bypassing the shared clock.
    #[must_use]
    pub fn reed_voice(&self, index: usize) -> Option<&ReedBoreVoice> {
        match self.voices.get(index) {
            Some(RenderVoice::ReedBore(voice)) => Some(voice),
            _ => None,
        }
    }

    /// Check an entire control batch without changing inputs, state or logs.
    ///
    /// # Errors
    /// A poisoned context, unknown voice index, or a delta the voice kind
    /// cannot accept.
    pub fn validate_controls(&self, deltas: &[ControlDelta]) -> Result<(), RenderError> {
        if self.poisoned {
            return Err(RenderError::Poisoned);
        }
        for voice in &self.voices {
            if let RenderVoice::BowedString(bow) = voice {
                if !bow.state().has_radiation() {
                    return Err(RenderError::Control {
                        what: "a mixed bowed voice requires a physical pressure observer",
                    });
                }
            }
        }
        // Validate everything first: control application is transactional.
        for delta in deltas {
            match delta {
                ControlDelta::SetFrictionDrive { voice, speed_m_s, normal_force_n } => {
                    match self.voices.get(*voice) {
                        Some(RenderVoice::FrictionModal(network)) => network.validate_drive(FrictionDrive {
                            speed_m_s: *speed_m_s, normal_force_n: *normal_force_n,
                        })?,
                        Some(_) => return Err(RenderError::Control { what: "friction controls require a friction-network voice" }),
                        None => return Err(RenderError::UnknownVoice { index: *voice }),
                    }
                }
                ControlDelta::SetPlateForce { voice, force_n } => {
                    match self.voices.get(*voice) {
                        Some(RenderVoice::CompactPlate(plate)) => plate.validate_force(*force_n)?,
                        Some(_) => return Err(RenderError::Control {
                            what: "a plate footprint force requires a compact plate voice",
                        }),
                        None => return Err(RenderError::UnknownVoice { index: *voice }),
                    }
                }
                ControlDelta::SetModalForce {
                    voice,
                    mode,
                    force_n_per_sqrt_kg,
                } => {
                    if !force_n_per_sqrt_kg.is_finite() {
                        return Err(RenderError::Control {
                            what: "modal force must be finite N/sqrt(kg)",
                        });
                    }
                    match self.voices.get(*voice) {
                        Some(RenderVoice::FrictionModal(network)) => network.validate_force(*mode, *force_n_per_sqrt_kg)?,
                        Some(RenderVoice::MultiContactModal(network)) => network.validate_force(*mode, *force_n_per_sqrt_kg)?,
                        Some(RenderVoice::ContactModal(network)) => network.validate_force(*mode, *force_n_per_sqrt_kg)?,
                        Some(RenderVoice::CoupledModal(network)) => network.validate_force(*mode, *force_n_per_sqrt_kg)?,
                        Some(RenderVoice::ModalString(string)) => {
                            if *mode >= string.held_force.len() {
                                return Err(RenderError::Control {
                                    what: "modal force index is outside the admitted basis",
                                });
                            }
                        }
                        Some(RenderVoice::CompactPlate(_)) => {
                            return Err(RenderError::Control {
                                what: "use a physical footprint force for a compact plate voice",
                            });
                        }
                        Some(RenderVoice::BowedString(_) | RenderVoice::ReedBore(_)) => {
                            return Err(RenderError::Control {
                                what: "this voice has no retained modal-force input",
                            });
                        }
                        None => return Err(RenderError::UnknownVoice { index: *voice }),
                    }
                }
                ControlDelta::SetBlowingPressure { voice, pressure_pa } => {
                    if !pressure_pa.is_finite() || *pressure_pa < 0.0 {
                        return Err(RenderError::Control {
                            what: "blowing pressure must be finite and non-negative",
                        });
                    }
                    match self.voices.get(*voice) {
                        Some(RenderVoice::ReedBore(_)) => {}
                        Some(RenderVoice::BowedString(_)) => return Err(RenderError::Control {
                            what: "a bowed string has physical bow controls, not blowing pressure",
                        }),
                        Some(RenderVoice::FrictionModal(_) | RenderVoice::CoupledModal(_) | RenderVoice::ContactModal(_) | RenderVoice::MultiContactModal(_)) => return Err(RenderError::Control {
                            what: "a coupled modal network has no blowing pressure",
                        }),
                        Some(RenderVoice::CompactPlate(_)) => {
                            return Err(RenderError::Control {
                                what: "a compact plate voice has no blowing pressure",
                            });
                        }
                        Some(RenderVoice::ModalString(_)) => {
                            return Err(RenderError::Control {
                                what: "a modal string voice has no blowing pressure",
                            });
                        }
                        None => {
                            return Err(RenderError::UnknownVoice { index: *voice });
                        }
                    }
                }
                ControlDelta::SetReedClosingForce { voice, force_n } => {
                    if !force_n.is_finite() || *force_n < 0.0 {
                        return Err(RenderError::Control {
                            what: "reed closing force must be finite and non-negative",
                        });
                    }
                    match self.voices.get(*voice) {
                        Some(RenderVoice::ReedBore(reed)) if reed.reed.mass_kg > 0.0 => {}
                        Some(RenderVoice::ReedBore(_)) => {
                            return Err(RenderError::Control {
                                what: "reed closing force requires a massive reed voice",
                            });
                        }
                        Some(_) => {
                            return Err(RenderError::Control {
                                what: "reed closing force requires a reed-bore voice",
                            });
                        }
                        None => {
                            return Err(RenderError::UnknownVoice { index: *voice });
                        }
                    }
                }
            }
        }
        Ok(())
    }

    /// Apply control deltas at the CURRENT block boundary (before the next
    /// `block` call). Refusals leave every voice untouched.
    ///
    /// # Errors
    /// A poisoned context, unknown voice index, or incompatible/invalid input.
    pub fn apply_controls(&mut self, deltas: &[ControlDelta]) -> Result<(), RenderError> {
        self.validate_controls(deltas)?;
        for delta in deltas {
            match delta {
                ControlDelta::SetFrictionDrive { voice, speed_m_s, normal_force_n } => {
                    if let Some(RenderVoice::FrictionModal(network)) = self.voices.get_mut(*voice) {
                        network.set_drive_admitted(FrictionDrive {
                            speed_m_s: *speed_m_s, normal_force_n: *normal_force_n,
                        });
                    }
                }
                ControlDelta::SetPlateForce { voice, force_n } => {
                    if let Some(RenderVoice::CompactPlate(plate)) = self.voices.get_mut(*voice) {
                        plate.set_force_admitted(*force_n);
                    }
                }
                ControlDelta::SetModalForce {
                    voice,
                    mode,
                    force_n_per_sqrt_kg,
                } => {
                    match self.voices.get_mut(*voice) {
                        Some(RenderVoice::FrictionModal(network)) => network.set_force_admitted(*mode, *force_n_per_sqrt_kg),
                        Some(RenderVoice::ModalString(string)) => string.held_force[*mode] = *force_n_per_sqrt_kg,
                        Some(RenderVoice::MultiContactModal(network)) => network.set_force_admitted(*mode, *force_n_per_sqrt_kg),
                        Some(RenderVoice::ContactModal(network)) => network.set_force_admitted(*mode, *force_n_per_sqrt_kg),
                        Some(RenderVoice::CoupledModal(network)) => network.set_force_admitted(*mode, *force_n_per_sqrt_kg),
                        _ => unreachable!("complete control batch was admitted"),
                    }
                }
                ControlDelta::SetBlowingPressure { voice, pressure_pa } => {
                    if let Some(RenderVoice::ReedBore(reed)) = self.voices.get_mut(*voice) {
                        reed.set_blowing_pressure(*pressure_pa);
                    }
                }
                ControlDelta::SetReedClosingForce { voice, force_n } => {
                    if let Some(RenderVoice::ReedBore(reed)) = self.voices.get_mut(*voice) {
                        reed.set_closing_force(*force_n);
                    }
                }
            }
            self.controls_applied.push(ControlRecord {
                block_index: self.blocks_rendered,
                delta: *delta,
                lift: "", // pure input-parameter move; no state lifted
            });
        }
        Ok(())
    }

    /// Render exactly one block: zero the output, advance every voice
    /// `out.len()` samples, sum. No allocation occurs for blocks within
    /// the pre-sized maximum and no-alloc voices; larger blocks refuse.
    ///
    /// # Errors
    /// Voice refusals poison the context (mid-sample state is not
    /// rewound); discard that block's partial output. Later calls refuse
    /// without touching output or state. An oversized or empty block
    /// refuses before any state moves and does not poison a healthy context.
    pub fn block(&mut self, out: &mut [f64]) -> Result<(), RenderError> {
        self.validate_block_len(out.len())?;
        out.fill(0.0);
        for voice in &mut self.voices {
            let scratch = &mut self.scratch[..out.len()];
            let result = match voice {
                RenderVoice::FrictionModal(network) => network.step_block(scratch),
                RenderVoice::BowedString(bow) => {
                    bow.pressure_block(scratch).map_err(RenderError::Bowed)
                }
                RenderVoice::MultiContactModal(network) => network.step_block(scratch),
                RenderVoice::ContactModal(network) => network.step_block(scratch),
                RenderVoice::CoupledModal(network) => network.step_block(scratch),
                RenderVoice::CompactPlate(plate) => plate.step_block(scratch),
                RenderVoice::ReedBore(reed) => reed.step_block(scratch).map_err(RenderError::Voice),
                RenderVoice::ModalString(string) => {
                    string.step_block(scratch).map_err(RenderError::Modal)
                }
            };
            if let Err(error) = result {
                self.poisoned = true;
                return Err(error);
            }
            for (accumulator, sample) in out.iter_mut().zip(scratch.iter()) {
                *accumulator += *sample;
            }
        }
        self.blocks_rendered += 1;
        // The conversion and both additions were admitted before stepping.
        self.samples_rendered += out.len() as u64;
        Ok(())
    }

    fn validate_block_len(&self, len: usize) -> Result<(), RenderError> {
        if self.poisoned {
            return Err(RenderError::Poisoned);
        }
        if len == 0 {
            return Err(RenderError::EmptyBlock);
        }
        if len > self.scratch.len() {
            return Err(RenderError::Sizing {
                what: "block exceeds the pre-sized maximum; grow max_block at construction",
            });
        }
        let samples = u64::try_from(len).map_err(|_| RenderError::Sizing {
            what: "block length cannot be represented by the sample clock",
        })?;
        if self.samples_rendered.checked_add(samples).is_none()
            || self.blocks_rendered.checked_add(1).is_none()
        {
            return Err(RenderError::Sizing {
                what: "render sample or block clock would overflow",
            });
        }
        // Admit every bow's complete callback before another voice advances or
        // output is cleared. A horizon/capacity refusal is not a physics failure.
        for voice in &self.voices {
            if let RenderVoice::BowedString(bow) = voice {
                if !bow.state().has_radiation() {
                    return Err(RenderError::Control {
                        what: "a mixed bowed voice requires a physical pressure observer",
                    });
                }
                if len > bow.state().max_block_len() || samples > bow.remaining_samples() {
                    return Err(RenderError::Sizing {
                        what: "callback exceeds a bowed voice's capacity or performance horizon",
                    });
                }
            }
        }
        Ok(())
    }
}

/// Outcome of a gated multi-block render.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GatedRenderOutcome {
    /// Every requested block rendered.
    Completed {
        /// Blocks rendered.
        blocks: u64,
    },
    /// Cancellation was observed at a block boundary; the render stopped
    /// cleanly after a whole number of blocks (request -> drain the
    /// current block -> stop).
    Cancelled {
        /// Whole blocks rendered before stopping.
        blocks: u64,
    },
}

/// Render `blocks` blocks of `block_len` samples into `out` under a
/// cancellation gate, polling ONLY at block boundaries: an in-flight
/// block always completes (drain semantics), so the rendered prefix is a
/// whole number of blocks and a resumed context continues
/// bitwise-identically from the boundary.
///
/// # Errors
/// Voice refusals, oversized blocks, or an output slice smaller than the
/// requested render.
pub fn render_under_gate(
    context: &mut RenderContext,
    gate: &CancelGate,
    out: &mut [f64],
    block_len: usize,
    blocks: usize,
) -> Result<GatedRenderOutcome, RenderError> {
    let required_samples = block_len.checked_mul(blocks).ok_or(RenderError::Sizing {
        what: "block_len * blocks overflows usize",
    })?;
    if block_len == 0 || out.len() < required_samples {
        return Err(RenderError::Sizing {
            what: "output slice must hold block_len * blocks samples",
        });
    }
    for index in 0..blocks {
        if gate.is_requested() {
            return Ok(GatedRenderOutcome::Cancelled {
                blocks: index as u64,
            });
        }
        let start = index * block_len;
        context.block(&mut out[start..start + block_len])?;
    }
    Ok(GatedRenderOutcome::Completed {
        blocks: blocks as u64,
    })
}
