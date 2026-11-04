//! A physical duct microphone on the ordinary material-assembly entry point.
//!
//! The moving-reed junction, regional TMM terminal transfer or native spatial
//! ODE terminal flow, retarded Rayleigh receiver and causal atmosphere filter
//! keep their original owners. Interior bore pressure and the compact jet proxy
//! are not microphone signals. Observation never rebuilds or retunes the reed,
//! its geometry, or either medium.

use super::{
    AcousticRealizeError, BoreRealization, BoreRealizationReceipt, RealizeOptions,
    RealizedAssembly, gas_state, physics_duct, realize_reed_on_duct_with_receiver, sample_count,
    selected_observer_law, wall_pin,
};
use crate::air_path::{GasAbsorptionModel, ObserverFilter, ObserverFilterReceipt};
use crate::driving_point::regional::RegionalDuct;
use crate::driving_point::regional::material_wall::MaterialWallLayout;
use crate::pcm_wav::baffled::{BaffledPressure, CircularOutletReceiver, RayleighMedium};
use crate::render::ReedBoreVoice;
use crate::thin_plate::PlateBank;
use fs_duct::{Termination, WallModel};
use fs_material::gas::GasState;
use fs_scenario::{AcousticAssembly, WaveguideEnd};

/// Explicit regional materials and exterior receiver for one reed-driven bore.
///
/// This bounded image admits an unbranched cylinder/cone chain with its existing
/// local wall law and flanged termination. Strings, separate body/plate/cavity
/// members, other excitation and tone holes refuse; none is silently dropped.
/// The explicit lumped ladder uses the original compression-volume/volume-flow
/// frustum owner, including for a single widening cone, with regional interior
/// gases, an independent pressure-matched exterior and its original numerical
/// wall. The exact same prepared chain supplies passive-pitch admission, reed
/// feedback and terminal flow. The legacy spherical psi-line is not used by
/// this physical image.
/// It observes the actual admitted loaded terminal flow;
/// spatial propagation remains in the LC state, without a separate TMM transfer.
/// The original compact/Foster mouth-load approximation remains; this receiver
/// neither replaces it nor certifies its continuum radiation match.
/// Interior cell storage/losses use each local gas and the reed uses the inlet; only the
/// terminal radiation load, receiver and atmosphere use the exterior. The
/// passive-pitch reference uses the same assignment through the original TMM.
/// Nonuniform interiors keep at least two cells per authored region within the
/// exact requested total, refining the greatest local cell travel time first.
/// Their spatial screening estimate does not certify material-interface mode
/// error: the measured regional TMM pitch gate remains decisive. Source-bound
/// walls still require the characteristic line. Side-hole and general body microphones
/// retain their existing dedicated owners and are not inferred from this input.
///
/// Frozen gas regions share one static pressure. This is not advective mixing,
/// evolving breath, wall radiation or experimental qualification of an instrument.
/// Any supplied uniform, regional or exterior gas requires an explicit observer
/// law. Stokes-Kirchhoff/lossless selections retain the labelled molecular-loss
/// omission; ISO-air or sourced molecular relaxation retains its own admission.
/// Absorption/dispersion uses the receiver's range from the outlet center; the
/// original finite-aperture Rayleigh paths retain that center-range approximation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BaffledDuctObserver<'a> {
    /// One state per authored axial segment, in geometry order. `None` uses
    /// `RealizeOptions::gas`, or the scenario air, throughout the interior.
    /// Supplying both this field and `RealizeOptions::gas` refuses ambiguity.
    pub segment_gases: Option<&'a [GasState]>,
    /// Independent homogeneous exterior at the same static pressure. If absent,
    /// the actual last interior region supplies radiation loading and observation.
    /// Its density/speed change the outlet load as well as the microphone.
    pub exterior_gas: Option<GasState>,
    /// Original uniform or complete segment material layout, including each
    /// state, thickness, axial restraint and local membrane spectrum. Pressure
    /// feedback and terminal flow share that same assignment. Every sampled
    /// radius and frequency retains its source's admission. The scenario's authored
    /// `duct.wall` must be absent when this binding is supplied.
    /// This does not add radiation from the wall or equilibrate its temperature
    /// with the independently specified interior or exterior gas.
    pub material_wall: Option<&'a MaterialWallLayout>,
    /// Receiver position and quadrature in the outlet's outward baffle frame.
    /// Its actual range replaces the legacy scalar listener distance.
    pub receiver: CircularOutletReceiver,
}

/// Physical media and delay accounting retained by the realized pressure.
#[derive(Clone, Debug, PartialEq)]
pub struct BaffledDuctObservation {
    /// Actual interior gas states in the original axial-segment order.
    pub segment_gases: Vec<GasState>,
    /// Same homogeneous exterior used by the load and the pressure receiver.
    pub exterior_gas: GasState,
    /// Exact source-backed layout used by both feedback and observation, when
    /// explicitly selected. Every original material identity, state and loss
    /// spectrum remains inspectable in authored segment order.
    pub material_wall: Option<MaterialWallLayout>,
    /// Original receiver geometry and quadrature, not an output gain.
    pub receiver: CircularOutletReceiver,
    /// Original causal atmospheric law. Its declared excess/realization delay
    /// and the receiver's delay together count the center-reference range/c
    /// once. This is not individual dispersive propagation on each surface ray.
    pub filter: ObserverFilterReceipt,
    /// Whole samples of axial propagation factored out of the characteristic
    /// terminal transfer. Zero for the lumped ladder: its propagation lives in
    /// the actual spatial state, with the reported spatial approximation error.
    pub axial_delay_samples: usize,
    /// Retarded Rayleigh kernel's minimum and maximum sample delays.
    pub exterior_delay_samples: (usize, usize),
}

fn invalid(what: &'static str) -> AcousticRealizeError {
    AcousticRealizeError::InvalidDescription { what }
}

pub(super) fn realize(
    assembly: &AcousticAssembly,
    has_chart: bool,
    options: &RealizeOptions<'_>,
    binding: BaffledDuctObserver<'_>,
) -> Result<RealizedAssembly, AcousticRealizeError> {
    if has_chart
        || assembly.string.is_some()
        || assembly.plate.is_some()
        || assembly.soundboard.is_some()
        || !assembly.body_modes.is_empty()
        || assembly.cavity.is_some()
        || assembly.pluck.is_some()
        || assembly.bow.is_some()
        || assembly.blow.is_some()
        || !assembly.obstacles.is_empty()
        || assembly.contact_texture.is_some()
    {
        return Err(invalid(
            "physical duct observation requires only the reed and bore; other coupled members or excitations need their own physical receiver",
        ));
    }
    let is_lumped = matches!(options.bore, BoreRealization::LumpedLadder { .. });
    if !is_lumped && options.bore != BoreRealization::default() {
        return Err(invalid(
            "physical duct observation requires the default characteristic line or an explicit lumped ladder",
        ));
    }
    let duct = assembly
        .duct
        .as_ref()
        .ok_or_else(|| invalid("physical duct observation requires a duct"))?;
    let reed = assembly
        .reed
        .ok_or_else(|| invalid("physical duct observation requires a reed excitation"))?;
    if !duct.tone_holes.is_empty() {
        return Err(invalid(
            "physical assembly observation currently requires an unbranched bore; side mouths need individual receiver paths",
        ));
    }
    if duct.termination != WaveguideEnd::FlangedOpen {
        return Err(invalid(
            "baffled duct observation requires the same flanged load in pressure feedback",
        ));
    }
    if binding.segment_gases.is_some() && options.gas.is_some() {
        return Err(invalid(
            "supply either regional duct gases or the uniform gas option, not both",
        ));
    }
    if binding.material_wall.is_some() && duct.wall.is_some() {
        return Err(invalid(
            "supply either a source-bound material wall or an authored duct wall, not both",
        ));
    }
    let physical = physics_duct(duct)?;
    let gases = match binding.segment_gases {
        Some(gases) => gases.to_vec(),
        None => {
            let gas = match options.gas {
                Some(gas) => gas,
                None => gas_state(assembly.ambient)?,
            };
            vec![gas; physical.segments.len()]
        }
    };
    let regional = RegionalDuct::new(&physical, &gases).map_err(AcousticRealizeError::Duct)?;
    let regional = match binding.exterior_gas {
        Some(gas) => regional
            .with_exterior_gas(gas)
            .map_err(AcousticRealizeError::Duct)?,
        None => regional,
    };
    let regional = match binding.material_wall {
        Some(MaterialWallLayout::Segments(walls)) => regional
            .with_material_wall_regions(walls)
            .map_err(AcousticRealizeError::Duct)?,
        _ => regional,
    };
    let exterior = *regional.exterior_gas();
    let interior = *regional.inlet_gas();
    if is_lumped && binding.material_wall.is_some() {
        return Err(invalid(
            "physical lumped outlet observation requires the original numerical wall; source-bound walls need their own spatial owner",
        ));
    }
    let [x, y, z] = binding.receiver.position_m;
    let range_m = x.hypot(y).hypot(z);
    if !range_m.is_finite() || range_m <= 0.0 {
        return Err(invalid(
            "physical receiver range must be finite and positive",
        ));
    }
    let explicit_medium =
        options.gas.is_some() || binding.segment_gases.is_some() || binding.exterior_gas.is_some();
    if explicit_medium && options.observer_absorption.is_none() {
        return Err(invalid(
            "physical duct material overrides require an explicit observer absorption law; select or explicitly omit molecular relaxation",
        ));
    }
    let law = selected_observer_law(assembly.ambient, &exterior, options.observer_absorption);
    if options.observer_relaxation.is_some() && law != GasAbsorptionModel::MolecularRelaxation {
        return Err(invalid(
            "an observer relaxation state needs the MolecularRelaxation observer law",
        ));
    }
    let mut observer = ObserverFilter::for_law(
        law,
        &exterior,
        options.observer_relaxation,
        range_m,
        assembly.sample_rate_hz,
    )?;
    let filter = *observer.receipt();
    let radius_m = physical
        .segments
        .last()
        .expect("admitted nonempty duct")
        .outlet_radius();
    let near_m = (x.hypot(y) - radius_m).max(0.0).hypot(z);
    let rate = f64::from(assembly.sample_rate_hz);
    let available_samples = near_m / exterior.sound_speed * rate;
    let needed_samples = near_m * filter.low_frequency_delay_s / range_m * rate + 1.0;
    if filter.delay_speed_m_s.is_none() || available_samples < needed_samples {
        return Err(AcousticRealizeError::ObserverLatency {
            range_m,
            available_samples,
            needed_samples,
            budget_energy: filter.budget_energy,
        });
    }
    let delay_speed_m_s = filter
        .delay_speed_m_s
        .expect("admitted positive receiver delay");
    let count = sample_count(assembly.sample_rate_hz, assembly.duration_s)?;
    if is_lumped {
        let mut microphone = BaffledPressure::circular_outlet(
            radius_m,
            assembly.sample_rate_hz,
            binding.receiver,
            RayleighMedium {
                density: exterior.density,
                sound_speed: delay_speed_m_s,
            },
        )
        .map_err(|error| {
            AcousticRealizeError::DuctObserver(
                crate::driving_point::regional::outlet::OutletError::Receiver(error),
            )
        })?;
        let exterior_delay_samples = microphone.delay_samples;
        let (mut pressure_pa, bore_receipt) = realize_reed_on_duct_with_receiver(
            duct,
            reed,
            &mut PlateBank::default(),
            &interior,
            range_m,
            assembly.sample_rate_hz,
            count,
            options.bore,
            Some((&mut microphone, regional)),
        )?;
        for pressure in &mut pressure_pa {
            *pressure = observer.process(*pressure);
            if !pressure.is_finite() {
                return Err(invalid("physical lumped duct pressure left the finite set"));
            }
        }
        return Ok(RealizedAssembly {
            sample_rate_hz: assembly.sample_rate_hz,
            pressure_pa,
            gas: exterior,
            bore_receipt: Some(bore_receipt),
            observer_absorption: filter.absorption,
            duct_observer: Some(BaffledDuctObservation {
                segment_gases: gases,
                exterior_gas: exterior,
                material_wall: None,
                receiver: binding.receiver,
                filter,
                axial_delay_samples: 0,
                exterior_delay_samples,
            }),
        });
    }
    let authored_wall = wall_pin(duct.wall);
    let wall = match binding.material_wall {
        Some(MaterialWallLayout::Uniform(material)) => WallModel::Law(material),
        Some(MaterialWallLayout::Segments(_)) => WallModel::Rigid,
        None => WallModel::from_pin(authored_wall.as_ref()),
    };
    let mut voice = ReedBoreVoice::new_with_regional_wall(
        &regional,
        reed,
        Termination::FlangedOpen,
        PlateBank::default(),
        range_m,
        assembly.sample_rate_hz,
        count,
        wall,
    )?;
    let line = voice.line_realization().clone();
    let mut microphone = regional
        .baffled_outlet_delayed_at(
            assembly.sample_rate_hz,
            line.n_fft,
            binding.receiver,
            wall,
            delay_speed_m_s,
        )
        .map_err(AcousticRealizeError::DuctObserver)?;
    let receiver = microphone.info();
    // Sample -1 is already in the original feedback history. The physical
    // terminal consumes exactly that same sample once before source sample zero.
    microphone
        .step(voice.outgoing_pressure_pa())
        .map_err(AcousticRealizeError::DuctObserver)?;
    let mut pressure_pa = Vec::with_capacity(count);
    for _ in 0..count {
        voice.step_block(&mut [0.0])?;
        let emitted = microphone
            .step(voice.outgoing_pressure_pa())
            .map_err(AcousticRealizeError::DuctObserver)?;
        let pressure = observer.process(emitted);
        if !pressure.is_finite() {
            return Err(invalid("physical duct pressure left the finite set"));
        }
        pressure_pa.push(pressure);
    }
    Ok(RealizedAssembly {
        sample_rate_hz: assembly.sample_rate_hz,
        pressure_pa,
        gas: exterior,
        bore_receipt: Some(BoreRealizationReceipt::CharacteristicLine(line)),
        observer_absorption: filter.absorption,
        duct_observer: Some(BaffledDuctObservation {
            segment_gases: gases,
            exterior_gas: exterior,
            material_wall: binding.material_wall.cloned(),
            receiver: binding.receiver,
            filter,
            axial_delay_samples: receiver.axial_delay_samples,
            exterior_delay_samples: receiver.exterior_delay_samples,
        }),
    })
}
