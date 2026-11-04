//! Material-bound, complete-field pressure acoustics with exact held-input steps.
//!
//! cargo run -p fs-couple --example heterogeneous_acoustic_exact -- cold-hot-air 4 160 0.00025 100
//! cargo run -p fs-couple --example heterogeneous_acoustic_exact -- air-water 8 320 0.000125 100
//! Arguments: CASE NX STEPS DT_S FREQUENCY_HZ. NX is even so the material
//! interface at x=0.5 m is an actual tetrahedral face. The 1 x 0.1 x 0.1 m tube
//! retains every P1 pressure/potential coordinate, including the potential gauge.
//! Preparation uses a bounded dense exponential; this is a small-domain example,
//! not a modal model or a sparse large-domain time solver.
//!
//! Dry air at 100 kPa, 300 K occupies the left half; the right half is either
//! dry air at 600 K or IF97 region-1 water at 100 kPa, 300 K. These original
//! material states supply both density and sound speed. The side walls are rigid.
//! Each end has the local real admittance 1/(rho*c) of its own fluid. This load
//! matches normally incident plane waves, not arbitrary three-dimensional waves.
//! A 1e-4 m/s sine velocity source acts at the inlet in parallel with that load.
//! Its midpoint value is HELD over each interval; exact propagation does not
//! remove this source approximation or the FE spatial error. Initial fields
//! are at rest, so the output includes the physical startup transient.
//!
//! Stdout is CSV of pressure, particle velocity and instantaneous intensity at
//! the SAME interface point approached from each fluid. The original one-sided
//! gradients are reported; their normal velocities are not averaged or forced
//! to agree. Stderr reports physical stored energy and cumulative source work
//! from the integrated conjugate pressure, with independently integrated
//! boundary absorption and the signed work defect. Refine NX and DT_S separately
//! at fixed duration/frequency; samples are not an antialiased audio recording.
//! Media/geometry are frozen, bulk loss is absent, and no heat/mass/phase transfer
//! or measured/converged acoustic accuracy is claimed by this demonstration.

use std::error::Error;

use fs_couple::cavity_modal::interior::heterogeneous::regions::FluidRegionsInterior;
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
    if args.len() != 6 {
        return Err(
            "usage: heterogeneous_acoustic_exact air-water|cold-hot-air NX STEPS DT_S FREQUENCY_HZ"
                .into(),
        );
    }
    let nx: usize = args[2].parse()?;
    let steps: usize = args[3].parse()?;
    let dt: f64 = args[4].parse()?;
    let frequency: f64 = args[5].parse()?;
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
    // The original pHS cap counts potential, momentum, the source and affine
    // coordinate, then doubles that block for the flow integrals. Never trim it.
    if pressure_dofs > limits.max_dofs
        || 2 * (2 * pressure_dofs + 2) > MAX_LINEAR_EXPONENTIAL_BLOCK_DIM
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
    let source = ExactVelocityPort {
        profile: end_faces(0.0)
            .into_iter()
            .map(|face| (face, [1.0; 3]))
            .collect(),
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
    let prepared =
        acoustic
            .operator()
            .prepare_exact_time(dt, &boundaries, &[source], limits, &mut || false)?;
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
    eprintln!(
        "time_s,kinetic_j,compression_j,source_work_j,absorbed_j,cumulative_defect_j,step_defect_j"
    );
    println!("time_s,region,pressure_pa,vx_m_s,vy_m_s,vz_m_s,ix_w_m2,iy_w_m2,iz_w_m2");
    let mut state = prepared.zero_state();
    let initial_energy = prepared.energy(&state)?.total_j();
    let mut supplied = 0.0;
    let mut absorbed = 0.0;
    for index in 0..steps {
        let phase = core::f64::consts::TAU * (frequency * dt) * (index as f64 + 0.5);
        let amplitude = 1e-4 * det::sin(phase);
        let step = prepared.step(&state, &[amplitude], &mut || false)?;
        // Independently use the returned conjugate-pressure integral. Its
        // multiplication by the actual held velocity is source work, not heat.
        supplied += amplitude * step.integrated_port_pressure_pa_m2_s[0];
        absorbed += step.work.boundary_dissipation_j;
        let energy = prepared.energy(&step.state)?;
        let time = (index + 1) as f64 * dt;
        let defect = energy.total_j() - initial_energy + absorbed - supplied;
        eprintln!(
            "{time},{},{},{supplied},{absorbed},{defect},{}",
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
