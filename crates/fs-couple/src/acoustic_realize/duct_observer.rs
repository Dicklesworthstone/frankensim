//! A physical duct microphone on the ordinary material-assembly entry point.
//!
//! The regional TMM, moving-reed junction, terminal volume-flow transfer,
//! retarded Rayleigh receiver and causal atmosphere filter keep their original
//! owners. Interior bore pressure and the compact jet proxy are not microphone
//! signals. This adapter observes the accepted outgoing wave without rebuilding
//! or retuning the reed, its geometry, or either medium.

use super::{
    AcousticRealizeError, BoreRealization, BoreRealizationReceipt, RealizeOptions,
    RealizedAssembly, gas_state, physics_duct, sample_count, selected_observer_law, wall_pin,
};
use crate::air_path::{GasAbsorptionModel, ObserverFilter, ObserverFilterReceipt};
use crate::driving_point::regional::RegionalDuct;
use crate::pcm_wav::baffled::CircularOutletReceiver;
use crate::render::ReedBoreVoice;
use crate::thin_plate::PlateBank;
use fs_duct::{Termination, WallModel};
use fs_material::gas::GasState;
use fs_scenario::{AcousticAssembly, WaveguideEnd};

/// Explicit regional materials and exterior receiver for one reed-driven bore.
///
/// This bounded image admits an unbranched cylinder/cone chain with its existing
/// local wall law and flanged termination. Strings, separate body/plate/cavity
/// members, other excitation, tone holes and nondefault bore discretizations
/// refuse; none is silently dropped. Side-hole and general body microphones
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
    /// Original receiver geometry and quadrature, not an output gain.
    pub receiver: CircularOutletReceiver,
    /// Original causal atmospheric law. Its declared excess/realization delay
    /// and the receiver's delay together count the center-reference range/c
    /// once. This is not individual dispersive propagation on each surface ray.
    pub filter: ObserverFilterReceipt,
    /// Whole samples of axial propagation retained by the terminal transfer.
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
    if options.bore != BoreRealization::default() {
        return Err(invalid(
            "physical regional duct observation requires the original default characteristic-line ladder",
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
    let exterior = *regional.exterior_gas();
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
    let wall = wall_pin(duct.wall);
    let mut voice = ReedBoreVoice::new_with_regional(
        &regional,
        reed,
        Termination::FlangedOpen,
        PlateBank::default(),
        range_m,
        assembly.sample_rate_hz,
        count,
        wall.as_ref(),
    )?;
    let line = voice.line_realization().clone();
    let mut microphone = regional
        .baffled_outlet_delayed_at(
            assembly.sample_rate_hz,
            line.n_fft,
            binding.receiver,
            WallModel::from_pin(wall.as_ref()),
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
            receiver: binding.receiver,
            filter,
            axial_delay_samples: receiver.axial_delay_samples,
            exterior_delay_samples: receiver.exterior_delay_samples,
        }),
    })
}
