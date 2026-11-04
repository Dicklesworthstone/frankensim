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
/// This bounded image admits a cylinder/cone chain with its existing local wall
/// law and flanged termination. The physical lumped route also admits native
/// tone holes with uniform gas, the same exterior and rigid walls. Strings,
/// separate body/plate/cavity members and other excitation retain their refusals.
/// The explicit lumped ladder uses the original compression-volume/volume-flow
/// frustum owner, including for a single widening cone, with regional interior
/// gases, an independent pressure-matched exterior and either the original
/// numerical wall or source-bound local material walls. The exact same prepared chain supplies passive-pitch admission, reed
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
/// error: the measured regional TMM pitch gate remains decisive. Source walls
/// retain their full equilibrium, Kelvin and Maxwell law at each native cell's
/// actual radius and slanted strip area; the same source supplies the TMM gate.
/// This is frozen local membrane motion, without axial shell coupling, static
/// pressure prestress or wall radiation. Source applicability and the original
/// exact-propagator frequency admission remain binding. Open/vented side mouths
/// use their actual native terminal flow (after all long-chimney compression),
/// effective open radius and individual receiver/atmosphere history. Closed pads
/// retain their load but emit no mouth flow. Their geometry uses the existing
/// parallel-baffle convention: every opening faces +z and a hole at axial station
/// x has receiver coordinate receiver.x + bore_length - x. This does not infer
/// an instrument's hole azimuths, arbitrary normals or unmodeled mutual radiation.
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
    /// Active side mouths in original bore order, excluding closed pads. Each
    /// uses its own range, delay and atmospheric filter before pressure summation.
    pub openings: Vec<BaffledDuctOpeningObservation>,
}

/// One original open/vented mouth in the explicitly selected parallel baffle.
#[derive(Clone, Debug, PartialEq)]
pub struct BaffledDuctOpeningObservation {
    /// Authored physical axial station [m], before numerical cell alignment.
    pub axial_position_m: f64,
    /// Native aperture radius including sqrt(open_fraction) [m].
    pub effective_radius_m: f64,
    /// Bell receiver translated by bore_length - axial_position along x;
    /// every mouth uses the same +z normal. No separate chimney delay is added.
    pub receiver: CircularOutletReceiver,
    /// This mouth's own center-range atmosphere and delay allocation.
    pub filter: ObserverFilterReceipt,
    /// Physical Rayleigh delay support on the source sample clock.
    pub exterior_delay_samples: (usize, usize),
}

struct SpatialMouthReceiver {
    microphone: BaffledPressure,
    atmosphere: ObserverFilter,
}

impl SpatialMouthReceiver {
    fn step(&mut self, flow: f64) -> Result<f64, AcousticRealizeError> {
        let pressure = self.microphone.step(&[flow]).map_err(receiver_error)?;
        Ok(self.atmosphere.process(pressure))
    }
}

// Original per-mouth receivers only; no source transfer or extra axial delay.
pub(super) struct SpatialDuctReceivers {
    bell: SpatialMouthReceiver,
    openings: Vec<SpatialMouthReceiver>,
}

impl SpatialDuctReceivers {
    pub(super) fn bell_only(microphone: BaffledPressure, atmosphere: ObserverFilter) -> Self {
        Self {
            bell: SpatialMouthReceiver {
                microphone,
                atmosphere,
            },
            openings: Vec::new(),
        }
    }

    pub(super) fn opening_count(&self) -> usize {
        self.openings.len()
    }

    pub(super) fn observe(
        &mut self,
        bell: &fs_phs::acoustic_outlet::AcousticOutletFlow,
        openings: &[fs_phs::acoustic_outlet::AcousticTapFlow],
        state: &[f64],
    ) -> Result<f64, AcousticRealizeError> {
        if openings.len() != self.openings.len() {
            return Err(invalid(
                "physical opening rows and receivers must remain aligned",
            ));
        }
        let flow = bell
            .volume_flow_m3_s(state)
            .map_err(super::render_phs_error)?;
        let mut pressure = self.bell.step(flow)?;
        for (opening, receiver) in openings.iter().zip(&mut self.openings) {
            let flow = opening
                .flow
                .volume_flow_m3_s(state)
                .map_err(super::render_phs_error)?;
            pressure += receiver.step(flow)?;
        }
        if !pressure.is_finite() {
            return Err(invalid("physical lumped duct pressure left the finite set"));
        }
        Ok(pressure)
    }
}

fn receiver_error(error: String) -> AcousticRealizeError {
    AcousticRealizeError::DuctObserver(
        crate::driving_point::regional::outlet::OutletError::Receiver(error),
    )
}

fn receiver_filter(
    law: GasAbsorptionModel,
    exterior: &GasState,
    options: &RealizeOptions<'_>,
    receiver: CircularOutletReceiver,
    radius_m: f64,
    sample_rate_hz: u32,
) -> Result<ObserverFilter, AcousticRealizeError> {
    let [x, y, z] = receiver.position_m;
    let range_m = x.hypot(y).hypot(z);
    if !range_m.is_finite() || range_m <= 0.0 {
        return Err(invalid(
            "physical receiver range must be finite and positive",
        ));
    }
    let observer = ObserverFilter::for_law(
        law,
        exterior,
        options.observer_relaxation,
        range_m,
        sample_rate_hz,
    )?;
    let filter = observer.receipt();
    let near_m = (x.hypot(y) - radius_m).max(0.0).hypot(z);
    let rate = f64::from(sample_rate_hz);
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
    Ok(observer)
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
    let has_holes = !duct.tone_holes.is_empty();
    if has_holes && !is_lumped {
        return Err(invalid(
            "physical assembly side-mouth observation requires the explicit lumped ladder",
        ));
    }
    if has_holes && (duct.wall.is_some() || binding.material_wall.is_some()) {
        return Err(invalid(
            "physical side mouths require rigid walls; chimney wall partition is not represented",
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
            vec![gas; duct.segments.len()]
        }
    };
    let regional_gases = if has_holes {
        if gases.len() != duct.segments.len()
            || gases.first().is_none()
            || gases.iter().any(|gas| Some(gas) != gases.first())
        {
            return Err(invalid(
                "physical side mouths require one uniform interior gas; regional tap ownership is not represented",
            ));
        }
        // The same uniform state also owns native zero-length hole records.
        vec![gases[0]; physical.segments.len()]
    } else {
        gases.clone()
    };
    let regional =
        RegionalDuct::new(&physical, &regional_gases).map_err(AcousticRealizeError::Duct)?;
    let regional = match binding.exterior_gas {
        Some(gas) => regional
            .with_exterior_gas(gas)
            .map_err(AcousticRealizeError::Duct)?,
        None => regional,
    };
    let regional = match binding.material_wall {
        Some(MaterialWallLayout::Segments(walls)) if !is_lumped => regional
            .with_material_wall_regions(walls)
            .map_err(AcousticRealizeError::Duct)?,
        _ => regional,
    };
    let exterior = *regional.exterior_gas();
    let interior = *regional.inlet_gas();
    if has_holes && exterior != interior {
        return Err(invalid(
            "physical side mouths require the same exterior gas as their native chimney load",
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
    let radius_m = physical
        .segments
        .last()
        .expect("admitted nonempty duct")
        .outlet_radius();
    let mut observer = receiver_filter(
        law,
        &exterior,
        options,
        binding.receiver,
        radius_m,
        assembly.sample_rate_hz,
    )?;
    let filter = *observer.receipt();
    let delay_speed_m_s = filter
        .delay_speed_m_s
        .expect("admitted positive receiver delay");
    let count = sample_count(assembly.sample_rate_hz, assembly.duration_s)?;
    if is_lumped {
        let microphone = BaffledPressure::circular_outlet(
            radius_m,
            assembly.sample_rate_hz,
            binding.receiver,
            RayleighMedium {
                density: exterior.density,
                sound_speed: delay_speed_m_s,
            },
        )
        .map_err(receiver_error)?;
        let exterior_delay_samples = microphone.delay_samples;
        let mut receivers = SpatialDuctReceivers::bell_only(microphone, observer);
        let mut openings = Vec::new();
        let length_m = regional.total_axial_length_m();
        for site in regional.open_hole_sites() {
            let mut receiver = binding.receiver;
            receiver.position_m[0] += length_m - site.axial_position_m;
            let atmosphere = receiver_filter(
                law,
                &exterior,
                options,
                receiver,
                site.effective_radius_m,
                assembly.sample_rate_hz,
            )?;
            let filter = *atmosphere.receipt();
            let microphone = BaffledPressure::circular_outlet(
                site.effective_radius_m,
                assembly.sample_rate_hz,
                receiver,
                RayleighMedium {
                    density: exterior.density,
                    sound_speed: filter.delay_speed_m_s.expect("admitted side-mouth delay"),
                },
            )
            .map_err(receiver_error)?;
            openings.push(BaffledDuctOpeningObservation {
                axial_position_m: site.axial_position_m,
                effective_radius_m: site.effective_radius_m,
                receiver,
                filter,
                exterior_delay_samples: microphone.delay_samples,
            });
            receivers.openings.push(SpatialMouthReceiver {
                microphone,
                atmosphere,
            });
        }
        let (pressure_pa, bore_receipt) = realize_reed_on_duct_with_receiver(
            duct,
            reed,
            &mut PlateBank::default(),
            &interior,
            range_m,
            assembly.sample_rate_hz,
            count,
            options.bore,
            Some((&mut receivers, regional, binding.material_wall)),
        )?;
        return Ok(RealizedAssembly {
            sample_rate_hz: assembly.sample_rate_hz,
            pressure_pa,
            gas: exterior,
            bore_receipt: Some(bore_receipt),
            observer_absorption: filter.absorption,
            duct_observer: Some(BaffledDuctObservation {
                segment_gases: gases,
                exterior_gas: exterior,
                material_wall: binding.material_wall.cloned(),
                receiver: binding.receiver,
                filter,
                axial_delay_samples: 0,
                exterior_delay_samples,
                openings,
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
            openings: Vec::new(),
        }),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn opening_observation_filters_each_native_flow_before_pressure_summation() {
        use fs_phs::{AcousticSection, AcousticTap, ChainOutlet, MouthFlange};
        let gas = GasState::try_new(
            &fs_material::gas::GasSpec::dry_air_ussa1976(),
            293.15,
            101_325.0,
        )
        .unwrap();
        let rate = 16_000;
        let sections = [AcousticSection {
            length: 0.3,
            radius: 0.004,
            outlet_radius: 0.004,
            cells: 6,
        }];
        let taps = [
            AcousticTap {
                station: 0.3,
                neck_length: 0.003,
                neck_radius: 0.001,
                open_fraction: 0.25,
            },
            AcousticTap {
                station: 0.5,
                neck_length: 0.01,
                neck_radius: 0.001,
                open_fraction: 0.0,
            },
            AcousticTap {
                station: 0.7,
                neck_length: 0.05,
                neck_radius: 0.00125,
                open_fraction: 1.0,
            },
        ];
        let (source, bell, holes) = fs_phs::acoustic_outlet::acoustic_chain_with_opening_flows(
            &sections,
            gas.density,
            gas.sound_speed,
            ChainOutlet::Radiating(MouthFlange::Flanged),
            1,
            &taps,
            None,
            None,
        )
        .unwrap();
        let receiver = CircularOutletReceiver {
            position_m: [0.25, 0.0, 0.75],
            radial_rings: 1,
            angular_points: 8,
            maximum_frequency_hz: 400.0,
        };
        let settings = RealizeOptions::default();
        let law = GasAbsorptionModel::StokesKirchhoff;
        let mut geometry = vec![(0.004, receiver)];
        for hole in &holes {
            let mut position = receiver;
            position.position_m[0] += sections[0].length * (1.0 - hole.station);
            geometry.push((hole.effective_radius_m, position));
        }
        let actual_mouth = |(radius, receiver): (f64, CircularOutletReceiver)| {
            let atmosphere = receiver_filter(law, &gas, &settings, receiver, radius, rate).unwrap();
            let microphone = BaffledPressure::circular_outlet(
                radius,
                rate,
                receiver,
                RayleighMedium {
                    density: gas.density,
                    sound_speed: atmosphere.receipt().delay_speed_m_s.unwrap(),
                },
            )
            .unwrap();
            SpatialMouthReceiver {
                microphone,
                atmosphere,
            }
        };
        let mut actual = SpatialDuctReceivers {
            bell: actual_mouth(geometry[0]),
            openings: geometry[1..].iter().copied().map(actual_mouth).collect(),
        };
        // Independent original receiver/filter instances use each physical
        // center range. Applying the bell's filter after summing is different.
        let mut expected: Vec<_> = geometry
            .iter()
            .map(|&(radius, receiver)| {
                let [x, y, z] = receiver.position_m;
                let filter =
                    ObserverFilter::for_law(law, &gas, None, x.hypot(y).hypot(z), rate).unwrap();
                let microphone = BaffledPressure::circular_outlet(
                    radius,
                    rate,
                    receiver,
                    RayleighMedium {
                        density: gas.density,
                        sound_speed: filter.receipt().delay_speed_m_s.unwrap(),
                    },
                )
                .unwrap();
                (microphone, filter)
            })
            .collect();
        assert_ne!(
            expected[0].1.receipt().range_m,
            expected[1].1.receipt().range_m
        );
        let mut state = vec![0.0; source.state_dim()];
        let mut amplitude = 0.0_f64;
        for sample in 0..128 {
            state = fs_phs::step(
                &source,
                &state,
                &[if sample < 20 { 1e-4 } else { 0.0 }],
                1.0 / f64::from(rate),
            )
            .unwrap()
            .x;
            let effort = source.effort(&state);
            // Actual physical effort coordinates: bell; compact vent; long
            // terminal after its two compression cells. No neck substitution.
            let flows = [effort[11], -effort[12], effort[17]];
            let mut independent = 0.0;
            for ((microphone, atmosphere), flow) in expected.iter_mut().zip(flows) {
                independent += atmosphere.process(microphone.step(&[flow]).unwrap());
            }
            let observed = actual.observe(&bell, &holes, &state).unwrap();
            amplitude = amplitude.max(independent.abs());
            assert!((observed - independent).abs() < 1e-12 * independent.abs().max(1e-12));
        }
        assert!(amplitude > 0.0);
    }
}
