//! Material-bound, complete-field pressure acoustics with exact held-input steps.
//!
//! cargo run -p fs-couple --example heterogeneous_acoustic_exact -- cold-hot-air 4 160 0.00025 100
//! cargo run -p fs-couple --example heterogeneous_acoustic_exact -- air-water 8 320 0.000125 100
//! Volume-only drive with a 1e-6 m^3/s peak injection into the left fluid:
//! cargo run -p fs-couple --example heterogeneous_acoustic_exact -- cold-hot-air 4 160 0.00025 100 0 1e-6
//! Arguments: CASE NX STEPS DT_S FREQUENCY_HZ [INLET_VELOCITY_M_S VOLUME_FLOW_M3_S].
//! The optional amplitudes must appear together and be finite and signed;
//! zero disables that input. Defaults are 1e-4 m/s and zero m^3/s.
//! NX is even so the material interface at x=0.5 m is an actual tetrahedral face.
//! The 1 x 0.1 x 0.1 m tube
//! retains every P1 pressure/potential coordinate, including the potential gauge.
//! Preparation uses a bounded dense exponential; this is a small-domain example,
//! not a modal model or a sparse large-domain time solver.
//!
//! Dry air at 100 kPa, 300 K occupies the left half; the right half is either
//! dry air at 600 K or IF97 region-1 water at 100 kPa, 300 K. These original
//! material states supply both density and sound speed. The side walls are rigid.
//! Each end has the local real admittance 1/(rho*c) of its own fluid. This load
//! matches normally incident plane waves, not arbitrary three-dimensional waves.
//! A sine velocity source acts at the inlet in parallel with that load.
//! The optional distributed source has a fixed shape of 200 m^-3 throughout
//! the original left half (0.005 m^3), zero elsewhere, and unit volume integral.
//! Its volume-flow amplitude [m^3/s] prescribes injection into continuity;
//! no mass-transfer or heat-release law is inferred. Both amplitudes multiply
//! the SAME sine. Its midpoint value is HELD over each interval; propagation does not
//! remove this source approximation or the FE spatial error. Initial fields
//! are at rest, so the output includes the physical startup transient.
//!
//! Stdout is CSV of pressure, particle velocity and instantaneous intensity at
//! the SAME interface point approached from each fluid. The original one-sided
//! gradients are reported; their normal velocities are not averaged or forced
//! to agree. Stderr reports physical stored energy and cumulative total source
//! work, with separate signed boundary and volume contributions from their
//! integrated conjugate pressures. It also reports independently integrated
//! boundary absorption and the signed work defect. Refine NX and DT_S separately
//! at fixed duration/frequency; samples are not an antialiased audio recording.
//! Media/geometry are frozen, bulk loss is absent, and no heat/mass/phase transfer
//! or measured/converged acoustic accuracy is claimed by this demonstration.

use std::error::Error;

use fs_couple::cavity_modal::interior::heterogeneous::regions::FluidRegionsInterior;
use fs_feec::acoustics::helmholtz::heterogeneous::HeterogeneousError;
use fs_feec::acoustics::helmholtz::heterogeneous::time::{
    ResistiveBoundary,
    exact::{ExactTimeLimits, ExactVelocityPort},
};
use fs_feec::acoustics::helmholtz::{HelmholtzFe, Medium};
use fs_material::fluid_acoustics::FluidAcoustics;
use fs_material::gas::{GasSpec, GasState};
use fs_math::det;
use fs_phs::linear_exact::MAX_LINEAR_EXPONENTIAL_BLOCK_DIM;
use fs_qty::{Pressure, Temperature};
use fs_thermochem::water::if97::{If97Region, If97State};

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args().collect();
    if args.len() != 6 && args.len() != 8 {
        return Err(
            "usage: heterogeneous_acoustic_exact air-water|cold-hot-air NX STEPS DT_S FREQUENCY_HZ [INLET_VELOCITY_M_S VOLUME_FLOW_M3_S]"
                .into(),
        );
    }
    let nx: usize = args[2].parse()?;
    let steps: usize = args[3].parse()?;
    let dt: f64 = args[4].parse()?;
    let frequency: f64 = args[5].parse()?;
    let (inlet_velocity_m_s, volume_flow_m3_s): (f64, f64) = if args.len() == 8 {
        (args[6].parse()?, args[7].parse()?)
    } else {
        (1e-4, 0.0)
    };
    if !inlet_velocity_m_s.is_finite() || !volume_flow_m3_s.is_finite() {
        return Err(
            "inlet velocity [m/s] and volume flow [m^3/s] must be finite signed amplitudes".into(),
        );
    }
    if !(2..=14).contains(&nx)
        || nx % 2 != 0
        || !(1..=200_000).contains(&steps)
        || !dt.is_finite()
        || dt <= 0.0
        || !(dt * steps as f64).is_finite()
        || !frequency.is_finite()
        || frequency <= 0.0
        || !(frequency * dt < 0.5)
    {
        return Err("require even NX in 2..=14, STEPS in 1..=200000, finite positive duration/frequency and FREQUENCY_HZ*DT_S < 0.5".into());
    }
    let limits = ExactTimeLimits::default();
    let pressure_dofs = 4 * (nx + 1);
    let input_count = usize::from(inlet_velocity_m_s != 0.0) + usize::from(volume_flow_m3_s != 0.0);
    // The original pHS cap counts potential, momentum, each input and affine
    // coordinate, then doubles that block for the flow integrals. Never trim it.
    if pressure_dofs > limits.max_dofs
        || 2 * (2 * pressure_dofs + input_count + 1) > MAX_LINEAR_EXPONENTIAL_BLOCK_DIM
    {
        return Err("requested full field exceeds the original exact preparation capacity".into());
    }
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 100_000.0)?;
    let right: Box<dyn FluidAcoustics> = match args[1].as_str() {
        "cold-hot-air" => Box::new(GasState::try_new(
            &GasSpec::dry_air_ussa1976(),
            600.0,
            100_000.0,
        )?),
        "air-water" => Box::new(If97State::at_pressure_temperature(
            Pressure::new(100_000.0),
            Temperature::new(300.0),
            If97Region::Region1,
        )?),
        _ => return Err("CASE must be air-water or cold-hot-air".into()),
    };
    let (complex, mut positions) = fs_feec::fixtures::masked_cube_grid(nx, 1, 1, &|_, _, _| true);
    for point in &mut positions {
        point[0] /= nx as f64;
        point[1] *= 0.1;
        point[2] *= 0.1;
    }
    let left = gas.acoustic_state()?;
    let geometry = HelmholtzFe::new(
        &complex,
        &positions,
        1,
        Medium {
            density_kg_m3: left.density_kg_m3(),
            sound_speed_m_s: left.sound_speed_m_s(),
        },
    )?;
    let regions: Vec<_> = geometry
        .tetrahedra()
        .iter()
        .map(|tet| {
            let x = tet.iter().map(|&v| positions[v as usize][0]).sum::<f64>() / 4.0;
            usize::from(x > 0.5)
        })
        .collect();
    let acoustic =
        FluidRegionsInterior::new(&geometry, &[&gas, right.as_ref()], &regions, 3, &mut || {
            false
        })?;
    let end_faces = |x: f64| {
        geometry
            .boundary_faces()
            .into_iter()
            .filter(|face| face.iter().all(|&v| positions[v as usize][0] == x))
            .collect::<Vec<_>>()
    };
    let ports = if inlet_velocity_m_s != 0.0 {
        vec![ExactVelocityPort {
            profile: end_faces(0.0)
                .into_iter()
                .map(|face| (face, [1.0; 3]))
                .collect(),
        }]
    } else {
        Vec::new()
    };
    let volume_source = if volume_flow_m3_s != 0.0 {
        Some(acoustic.operator().volume_velocity_source(
            3,
            &mut |element, _| {
                Ok::<_, HeterogeneousError>(if regions[element] == 0 { 200.0 } else { 0.0 })
            },
            &mut || false,
        )?)
    } else {
        None
    };
    let mut boundaries = Vec::with_capacity(2);
    for region in 0..2 {
        let medium = acoustic
            .region_state(region)
            .ok_or("unused example fluid")?;
        boundaries.push(ResistiveBoundary {
            faces: end_faces(region as f64),
            admittance_m_pa_s: 1.0 / medium.characteristic_impedance_pa_s_m(),
        });
        eprintln!(
            "# region={region},rho_kg_m3={},c_m_s={}",
            medium.density_kg_m3(),
            medium.sound_speed_m_s()
        );
    }
    let prepared = if let Some(source) = &volume_source {
        acoustic.operator().prepare_exact_time_with_volume_sources(
            dt,
            &boundaries,
            &ports,
            &[source],
            limits,
            &mut || false,
        )?
    } else {
        acoustic
            .operator()
            .prepare_exact_time(dt, &boundaries, &ports, limits, &mut || false)?
    };
    // Both probes use one fixed interface triangle, with each element's own
    // material density and original potential gradient. No interface averaging.
    let face = geometry
        .tetrahedra()
        .iter()
        .find_map(|tet| {
            let face: Vec<_> = tet
                .iter()
                .copied()
                .filter(|&v| positions[v as usize][0] == 0.5)
                .collect();
            (face.len() == 3).then_some(face)
        })
        .ok_or("mesh has no conforming interface face")?;
    let mut probes = Vec::with_capacity(2);
    for region in 0..2 {
        let (element, tet) = geometry
            .tetrahedra()
            .iter()
            .enumerate()
            .find(|(element, tet)| {
                regions[*element] == region && face.iter().all(|vertex| tet.contains(vertex))
            })
            .ok_or("interface lacks an original neighboring fluid element")?;
        let barycentric = tet.map(|vertex| {
            if face.contains(&vertex) {
                1.0 / 3.0
            } else {
                0.0
            }
        });
        let probe = acoustic.probe(element, barycentric)?;
        eprintln!(
            "# probe_region={region},position_m={:?}",
            probe.position_m()
        );
        probes.push(probe);
    }
    eprintln!(
        "# full_P1_pressure_dofs={},tetrahedra={},dt_s={dt},source_hz={frequency}",
        geometry.ndof(),
        geometry.tetrahedra().len()
    );
    eprintln!("# inlet_velocity_m_s={inlet_velocity_m_s},volume_flow_m3_s={volume_flow_m3_s}");
    eprintln!(
        "time_s,kinetic_j,compression_j,source_work_j,absorbed_j,cumulative_defect_j,step_defect_j,boundary_source_work_j,volume_source_work_j"
    );
    println!("time_s,region,pressure_pa,vx_m_s,vy_m_s,vz_m_s,ix_w_m2,iy_w_m2,iz_w_m2");
    let mut state = prepared.zero_state();
    let initial_energy = prepared.energy(&state)?.total_j();
    let mut supplied = 0.0;
    let mut boundary_supplied = 0.0;
    let mut volume_supplied = 0.0;
    let mut absorbed = 0.0;
    for index in 0..steps {
        let phase = core::f64::consts::TAU * (frequency * dt) * (index as f64 + 0.5);
        let sine = det::sin(phase);
        let inlet_amplitudes = [inlet_velocity_m_s * sine];
        let inlet_amplitudes = if ports.is_empty() {
            &[][..]
        } else {
            &inlet_amplitudes[..]
        };
        let volume_amplitude = volume_flow_m3_s * sine;
        let step = if volume_source.is_some() {
            prepared.step_with_volume_sources(
                &state,
                inlet_amplitudes,
                &[volume_amplitude],
                &mut || false,
            )?
        } else {
            prepared.step(&state, inlet_amplitudes, &mut || false)?
        };
        // Independently use the returned conjugate-pressure integral. Its
        // multiplication by each actual held amplitude is signed source work.
        let boundary_work: f64 = inlet_amplitudes
            .iter()
            .zip(&step.integrated_port_pressure_pa_m2_s)
            .map(|(amplitude, integral)| amplitude * integral)
            .sum();
        let volume_work = if volume_source.is_some() {
            volume_amplitude * step.integrated_volume_source_pressure_pa_s[0]
        } else {
            0.0
        };
        boundary_supplied += boundary_work;
        volume_supplied += volume_work;
        supplied += boundary_work + volume_work;
        absorbed += step.work.boundary_dissipation_j;
        let energy = prepared.energy(&step.state)?;
        let time = (index + 1) as f64 * dt;
        let defect = energy.total_j() - initial_energy + absorbed - supplied;
        eprintln!(
            "{time},{},{},{supplied},{absorbed},{defect},{},{boundary_supplied},{volume_supplied}",
            energy.kinetic_j,
            energy.compression_j,
            step.work.balance_residual_j()
        );
        for (region, probe) in probes.iter().enumerate() {
            let sample = probe.transient(&step.state)?;
            let [vx, vy, vz] = sample.particle_velocity_m_s;
            let [ix, iy, iz] = sample.instantaneous_intensity_w_m2;
            println!(
                "{time},{region},{},{vx},{vy},{vz},{ix},{iy},{iz}",
                sample.pressure_pa
            );
        }
        state = step.state;
    }
    Ok(())
}
