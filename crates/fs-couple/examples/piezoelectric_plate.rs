//! Nominal PZT-5A plate: the actual 3D source/tensor/electrode solve.
//!
//! Run: cargo run -p fs-couple --example piezoelectric_plate
//! Frequency response: append `-- --harmonic` for open/resistive/RL shunts.
//! Acoustic RLC sensing: append `-- --acoustic` for a periodic guide-driven run.
//! Data and approximation boundaries: examples/PIEZOELECTRIC_PLATE.md.

use fs_blake3::ContentHash;
use fs_couple::piezoelectric::spatial::{
    ElectrodeCondition, PiezoelectricElectrode, TetPiezoelectricControls, TetPiezoelectricModel,
    harmonic::{
        HarmonicElectrode, HarmonicElectrodeCondition, HarmonicElectrodeLoad,
        HarmonicPiezoelectricControls,
        acoustic::{HarmonicAcousticPort, PlaneWavePortControls},
    },
    transient::TetPiezoelectricTransient,
    transient::acoustic::TransientAcousticPort,
    transient::load::reactive::{
        ReactiveElectrodeLoad, ReactiveElectrodeState, ReactivePiezoelectricTransient,
    },
    transient::load::{ElectrodeDrive, ResistiveElectrodeLoad},
};
use fs_evidence::ValidityDomain;
use fs_exec::{CancelGate, Cx};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::conductor::{ELECTRICAL_RESISTIVITY_DIMS, resolve_uniform_conductor};
use fs_material::field_storage::{
    ELECTRIC_PERMITTIVITY_DIMS, ELECTRIC_PERMITTIVITY_PROPERTY,
    MAGNETIC_PERMEABILITY_DIMS, MAGNETIC_PERMEABILITY_PROPERTY, ResolvedCapacitor, ResolvedInductor,
};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_material::piezoelectric::tensor::{TensorPiezoelectric, waves::PiezoelectricWaveCondition};
use fs_material::state_point::MaterialPropertySelection;
use fs_math::c64::C64;
use fs_solid::linear3::TetAssemblyBudget;
use std::sync::Arc;

const SOURCE: &str =
    "https://info.piezo.com/hubfs/Data-Sheets/piezo-material-properties-data-sheet-20201112.pdf";
const FRAME: ContentHash = ContentHash([59; 32]);
const EPSILON_0: f64 = 8.854_187_818_8e-12;

fn nominal_material() -> Result<TensorPiezoelectric, Box<dyn std::error::Error>> {
    // Explicit ideal transversely isotropic poled-ceramic model. The source
    // reports nominal small-signal coefficients, not a full measured tensor
    // with covariance. All symmetry-filled entries below are model assumptions.
    let mut compliance = [[0.0; 6]; 6];
    compliance[0][0] = 15.1e-12;
    compliance[1][1] = 15.1e-12;
    compliance[2][2] = 18.6e-12;
    compliance[0][1] = -4.8e-12;
    compliance[1][0] = -4.8e-12;
    for i in 0..2 {
        compliance[i][2] = -7.6e-12;
        compliance[2][i] = -7.6e-12;
    }
    compliance[3][3] = (15.1 + 4.8) * 1e-12; // s66_eng/2 = s11-s12
    compliance[4][4] = 40.0e-12 / 2.0;
    compliance[5][5] = 40.0e-12 / 2.0;
    let mut d = [[0.0; 6]; 3];
    d[0][5] = 460e-12 * core::f64::consts::FRAC_1_SQRT_2;
    d[1][4] = 460e-12 * core::f64::consts::FRAC_1_SQRT_2;
    d[2][0] = -190e-12;
    d[2][1] = -190e-12;
    d[2][2] = 390e-12;
    let epsilon = [1600.0 * EPSILON_0, 1600.0 * EPSILON_0, 1800.0 * EPSILON_0];
    let mut claims = ClaimSet::new();
    for role in TensorPiezoelectric::strain_charge_requirements()? {
        let name = role.name();
        let value = if name == "density" {
            7950.0
        } else {
            let suffix = name
                .rsplit('_')
                .next()
                .expect("matrix source schema")
                .as_bytes();
            let (i, j) = ((suffix[0] - b'1') as usize, (suffix[1] - b'1') as usize);
            if name.starts_with("piezoelectric_compliance") {
                compliance[i][j]
            } else if name.starts_with("piezoelectric_strain_charge") {
                d[i][j]
            } else if i == j {
                epsilon[i]
            } else {
                0.0
            }
        };
        let dims = role.dims();
        claims.insert_claim(PropertyClaim {key:PropertyKey::new(name,dims),
            value:PropertyValue::Scalar{value,dims},
            // Nominal frozen research model. This does not invent a measured
            // temperature/frequency/amplitude validity band for the datasheet.
            validity:ValidityDomain::unconstrained(),uncertainty:UncertaintyModel::Unstated,
            interpolation:InterpolationPolicy::ConstantWithinValidity,observations:vec![],
            provenance:Provenance {source:format!("{SOURCE}; nominal coefficients; authored transverse-isotropy and quasistatic approximation"),
                license:"Factual nominal values; source document rights retained by publisher".into(),artifact:None}})?;
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "nominal PZT-5A".into(),
            phase: "solid".into(),
            process: "authored ideal poled-ceramic small-signal model from public typical values"
                .into(),
            revision: 0,
        },
        claims,
        vec![],
    )?;
    Ok(TensorPiezoelectric::resolve_strain_charge(
        &card,
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
        FRAME,
    )?)
}

fn declared_load() -> Result<ResistiveElectrodeLoad, Box<dyn std::error::Error>> {
    // An explicitly authored 1 kohm circuit load, not a measured ceramic loss.
    // rho=1e-3 ohm m, length=10 mm and area=0.01 mm^2 define R=rho L/A.
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new("electrical_resistivity", ELECTRICAL_RESISTIVITY_DIMS),
        value: PropertyValue::Scalar {
            value: 1e-3,
            dims: ELECTRICAL_RESISTIVITY_DIMS,
        },
        validity: ValidityDomain::unconstrained(),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![],
        provenance: Provenance {
            source: "authored fixed Ohmic load for the piezoelectric plate example".into(),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    })?;
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "authored external resistor".into(),
            phase: "solid".into(),
            process: "declared constant resistivity; no temperature model".into(),
            revision: 0,
        },
        claims,
        vec![],
    )?;
    Ok(ResistiveElectrodeLoad {
        electrodes: [1, 0],
        conductor: std::sync::Arc::new(resolve_uniform_conductor(
            &card,
            &QueryPoint::new(),
            MaterialPropertySelection::SingleClaimOnly,
            0.01,
            1e-8,
        )?),
    })
}

fn declared_winding() -> Result<ResolvedInductor, Box<dyn std::error::Error>> {
    // Authored ideal uniform-field winding: 1000 turns, 100 mm magnetic path,
    // 100 mm^2 cross-section and nominal vacuum permeability. This is an
    // approximation, not a measured coil, core-loss or parasitic model.
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(MAGNETIC_PERMEABILITY_PROPERTY, MAGNETIC_PERMEABILITY_DIMS),
        value: PropertyValue::Scalar {
            value: 4.0 * core::f64::consts::PI * 1e-7,
            dims: MAGNETIC_PERMEABILITY_DIMS,
        },
        validity: ValidityDomain::unconstrained(),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![],
        provenance: Provenance {
            source:
                "authored ideal uniform-field vacuum winding for the piezoelectric plate example"
                    .into(),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    })?;
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "authored ideal vacuum winding".into(),
            phase: "uniform magnetic field".into(),
            process: "fixed permeability and geometry; no parasitic or loss model".into(),
            revision: 0,
        },
        claims,
        vec![],
    )?;
    Ok(ResolvedInductor::resolve(
        &card,
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
        1e-4,
        0.1,
        1000,
    )?)
}

fn declared_capacitor() -> Result<ResolvedCapacitor, Box<dyn std::error::Error>> {
    // Separate authored dielectric: epsilon=1e-8 F/m, A=100 mm^2, gap=1 mm.
    // Its 1 nF is additional external storage, not the ceramic capacitance.
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(ELECTRIC_PERMITTIVITY_PROPERTY, ELECTRIC_PERMITTIVITY_DIMS),
        value: PropertyValue::Scalar { value: 1e-8, dims: ELECTRIC_PERMITTIVITY_DIMS },
        validity: ValidityDomain::unconstrained(), uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity, observations: vec![],
        provenance: Provenance {
            source: "authored ideal fixed dielectric capacitor for the piezoelectric plate example".into(),
            license: "CC0-1.0".into(), artifact: None,
        },
    })?;
    let card = MaterialCard::assemble(MaterialStateId {
        chemistry: "authored ideal dielectric".into(), phase: "solid".into(),
        process: "uniform fixed field; no leakage or dielectric loss".into(), revision: 0,
    }, claims, vec![])?;
    Ok(ResolvedCapacitor::resolve(&card, &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly, 1e-4, 1e-3)?)
}

/// A balanced mechanical drive powers the actual electrical load. Each row is
/// an independent steady periodic solve of the same finite plate, not a sweep
/// of transient samples or a fitted terminal resonance.
#[allow(clippy::too_many_arguments)]
fn harmonic_response(
    model: &TetPiezoelectricModel<'_>,
    supports: &[(usize, f64)],
    bottom: &[usize],
    top: &[usize],
    nodes: usize,
    controls: TetPiezoelectricControls,
    cx: &Cx<'_>,
) -> Result<(), Box<dyn std::error::Error>> {
    let resistor = declared_load()?;
    let winding = std::sync::Arc::new(declared_winding()?);
    let load_cases = [
        ("open", vec![]),
        (
            "resistor",
            vec![HarmonicElectrodeLoad::Resistive(resistor.clone())],
        ),
        (
            "series_rl",
            vec![HarmonicElectrodeLoad::Reactive(
                ReactiveElectrodeLoad::SeriesRl {
                    electrodes: resistor.electrodes,
                    inductor: winding.clone(),
                    conductor: Some(resistor.conductor.clone()),
                },
            )],
        ),
    ];
    let prescribed: Vec<_> = supports
        .iter()
        .map(|&(i, x)| (i, C64::from_re(x)))
        .collect();
    let electrodes = [
        HarmonicElectrode {
            nodes: bottom,
            condition: HarmonicElectrodeCondition::Voltage(C64::ZERO),
        },
        HarmonicElectrode {
            nodes: top,
            condition: HarmonicElectrodeCondition::Current(C64::ZERO),
        },
    ];
    let mut forces = vec![[C64::ZERO; 3]; nodes];
    // A 10 mN peak resultant, distributed equally over each face's nodes.
    // This is a declared nodal-force pattern, not a uniform-traction claim.
    for (face, sign) in [(top, 1.0), (bottom, -1.0)] {
        for &node in face {
            forces[node][2] = C64::from_re(sign * 0.01 / face.len() as f64);
        }
    }
    let charges = vec![C64::ZERO; nodes];
    let controls = HarmonicPiezoelectricControls {
        field: controls,
        power_tolerance_w: 1e-10,
        maximum_real_coordinates: 8 * nodes,
    };
    println!("Nominal PZT-5A: 10 x 10 x 2 mm finite plate; peak phasors use exp(-i omega t).");
    println!("Source: {SOURCE}");
    println!(
        "Authored frozen lossless ceramic approximation; 10 mN balanced nodal force, no measured-device accuracy claim."
    );
    println!(
        "External R={:.6e} ohm; ideal L={:.6e} H; only the external resistor produces heat.",
        resistor.conductor.resistance_ohm(),
        winding.inductance_h()
    );
    println!(
        "frequency_hz,load,peak_voltage_v,peak_thickness_motion_m,peak_load_current_a,mean_resistor_heat_w,mean_mechanical_power_w,power_residual_w"
    );
    for frequency_hz in [1e3, 2e4, 5e4, 1e5, 2e5] {
        for (name, loads) in &load_cases {
            let response = model.solve_harmonic(
                core::f64::consts::TAU * frequency_hz,
                &prescribed,
                &electrodes,
                loads,
                &forces,
                &charges,
                controls,
                cx,
            )?;
            let mean_motion = |face: &[usize]| {
                face.iter()
                    .fold(C64::ZERO, |sum, &i| sum + response.displacement_m[i][2])
                    .scale(1.0 / face.len() as f64)
            };
            let thickness = mean_motion(top) - mean_motion(bottom);
            let voltage = response.electrodes[1].voltage_v - response.electrodes[0].voltage_v;
            let load_current = response
                .loads
                .first()
                .map_or(0.0, |load| load.current_a.abs());
            println!(
                "{frequency_hz:.6e},{name},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.3e}",
                voltage.abs(),
                thickness.abs(),
                load_current,
                response.average_load_heat_w,
                response.average_mechanical_power_w,
                response.power_balance_residual_w
            );
        }
    }
    Ok(())
}

/// A single in-band acoustic tone drives the complete retained R/L/C system.
/// Initialize its periodic history rather than launching an undeclared broadband
/// mechanical startup. Midpoint's frequency map makes this initialization agree
/// with the actual time-discrete equations; it is not a fitted forcing gain.
#[allow(clippy::too_many_arguments)]
fn acoustic_response(
    model: &TetPiezoelectricModel<'_>, supports: &[(usize, f64)],
    bottom: &[usize], top: &[usize], tets: &[[usize; 4]],
    zeros: &[[f64; 3]], zero_charge: &[f64],
    controls: TetPiezoelectricControls, cx: &Cx<'_>,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut faces = std::collections::BTreeSet::new();
    for tet in tets {
        for omit in 0..4 {
            let face: Vec<_> = (0..4).filter(|&i| i != omit).map(|i| tet[i]).collect();
            if face.iter().all(|i| top.contains(i)) {
                let mut face = [face[0], face[1], face[2]];
                face.sort_unstable();
                faces.insert(face);
            }
        }
    }
    let port = model.plane_wave_acoustic_port(&faces.into_iter().collect::<Vec<_>>(), 1e-4,
        PlaneWavePortControls {
            maximum_surface_dofs: 3 * top.len(), coplanarity_tolerance: 1e-10,
            maximum_wavenumber_rad_m: 0.8 * core::f64::consts::PI / 0.01,
            maximum_relative_amplitude: 1e-3,
        }, cx)?;
    // Declared constant fluid coefficients; not an empirically validated water model.
    let medium = FluidAcousticState::try_new(1000.0, 1480.0)?;
    let resistor = declared_load()?;
    let winding = Arc::new(declared_winding()?);
    let capacitor = Arc::new(declared_capacitor()?);
    let loads = vec![
        ReactiveElectrodeLoad::SeriesRl { electrodes: [1, 0], inductor: winding.clone(),
            conductor: Some(resistor.conductor.clone()) },
        ReactiveElectrodeLoad::Capacitor { electrodes: [1, 0], capacitor: capacitor.clone() },
    ];
    let frequency_hz = 5000.0;
    let dt = 1.0 / (64.0 * frequency_hz);
    let omega = core::f64::consts::TAU * frequency_hz;
    let half_phase = 0.5 * omega * dt;
    let discrete_omega = 2.0 * half_phase.tan() / dt;
    let peak_pressure = 1000.0;
    let controls = TetPiezoelectricControls { energy_tolerance_j: 1e-18,
        charge_tolerance_c: 1e-18, ..controls };
    let prescribed: Vec<_> = supports.iter().map(|&(i, v)| (i, C64::from_re(v))).collect();
    let helectrodes = [
        HarmonicElectrode { nodes: bottom, condition: HarmonicElectrodeCondition::Voltage(C64::ZERO) },
        HarmonicElectrode { nodes: top, condition: HarmonicElectrodeCondition::Current(C64::ZERO) },
    ];
    let hloads: Vec<_> = loads.iter().cloned().map(HarmonicElectrodeLoad::Reactive).collect();
    let steady = model.solve_harmonic_with_acoustic_ports(discrete_omega,
        &[HarmonicAcousticPort { port: &port, medium: &medium,
            incident_pressure_pa: C64::from_re(peak_pressure / half_phase.cos()) }],
        None, &prescribed, &helectrodes, &hloads,
        &vec![[C64::ZERO; 3]; zeros.len()], &vec![C64::ZERO; zeros.len()],
        HarmonicPiezoelectricControls { field: controls, power_tolerance_w: 1e-12,
            maximum_real_coordinates: 8 * zeros.len() }, cx)?;
    let displacement: Vec<_> = steady.displacement_m.iter().map(|u| u.map(|v| v.re)).collect();
    let velocity: Vec<_> = steady.velocity_m_s.iter().map(|u| u.map(|v| v.re)).collect();
    let body = TetPiezoelectricTransient::new(model, supports, &[
        PiezoelectricElectrode { nodes: bottom, condition: ElectrodeCondition::Voltage(0.0) },
        PiezoelectricElectrode { nodes: top, condition: ElectrodeCondition::Charge(steady.electrodes[1].charge_c.re) },
    ], &displacement, &velocity, zero_charge, controls, cx)?;
    let flux = winding.inductance_h() * steady.loads[0].current_a.re;
    let charge = capacitor.capacitance_f() * steady.loads[1].voltage_v.re;
    let initial_energy = body.state().total_energy_j()
        + winding.stored_energy_j(flux)? + capacitor.stored_energy_j(charge)?;
    let mut coupled = ReactivePiezoelectricTransient::new(body, loads,
        vec![ReactiveElectrodeState::FluxLinkage(flux), ReactiveElectrodeState::Charge(charge)],
        vec![], controls, cx)?;
    println!("Nominal PZT-5A acoustic sensor: 5 kHz, 1000 Pa peak, periodic initial history, 64 steps/cycle.");
    println!("Source: {SOURCE}");
    println!("Authored lossless fluid rho=1000 kg/m^3, c=1480 m/s; 10 mm square uniform guide, first transverse cutoff=74 kHz.");
    println!("External R={:.6e} ohm, L={:.6e} H, C={:.6e} F; no external electrical or nodal-force supply.",
        resistor.conductor.resistance_ohm(), winding.inductance_h(), capacitor.capacitance_f());
    println!("time_s,voltage_v,mean_thickness_motion_m,external_reactive_energy_j,conductor_heat_j,incident_energy_j,outgoing_energy_j,net_outward_work_j,balance_j");
    let (mut heat, mut incident, mut outgoing, mut outward) = (0.0, 0.0, 0.0, 0.0);
    for i in 0..128 {
        let pressure = peak_pressure * (omega * (f64::from(i) + 0.5) * dt).cos();
        let candidate = coupled.step_with_acoustic_ports_candidate(dt,
            &[TransientAcousticPort { port: &port, medium: &medium,
                midpoint_incident_pressure_pa: pressure,
                maximum_angular_frequency_rad_s: 2.0 * omega }], None,
            &[ElectrodeDrive::Voltage(0.0), ElectrodeDrive::Current(0.0)],
            zeros, zero_charge, controls, cx)?;
        let step = candidate.step();
        let wave = &step.acoustic_ports[0];
        heat += step.joule_heat_j;
        incident += wave.incident_energy_j;
        outgoing += wave.outgoing_energy_j;
        outward += wave.net_outward_work_j;
        if (i + 1) % 16 == 0 {
            let state = candidate.model().body().state();
            let mean = |face: &[usize]| face.iter().map(|&n| state.displacement_m()[n][2]).sum::<f64>() / face.len() as f64;
            let balance = state.total_energy_j() + step.reactive_energy_j[1] + heat + outward - initial_energy;
            println!("{:.6e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.3e}",
                state.time_s(), state.electrodes()[1].voltage_v, mean(top) - mean(bottom),
                step.reactive_energy_j[1], heat, incident, outgoing, outward, balance);
        }
        coupled = candidate.into_parts().0;
    }
    Ok(())
}

/// Drive the actual finite plate, then disconnect its driven electrode without
/// changing its charge, displacement or velocity. The open electrode voltage
/// subsequently follows the computed motion, rather than an imposed waveform.
fn transient_response(
    model: &TetPiezoelectricModel<'_>,
    supports: &[(usize, f64)],
    bottom: &[usize],
    top: &[usize],
    zeros: &[[f64; 3]],
    zero_charge: &[f64],
    controls: TetPiezoelectricControls,
    cx: &Cx<'_>,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut body = TetPiezoelectricTransient::new(
        model,
        supports,
        &[
            PiezoelectricElectrode {
                nodes: bottom,
                condition: ElectrodeCondition::Voltage(0.0),
            },
            PiezoelectricElectrode {
                nodes: top,
                condition: ElectrodeCondition::Voltage(0.0),
            },
        ],
        zeros,
        zeros,
        zero_charge,
        controls,
        cx,
    )?;
    let thickness = |u: &[[f64; 3]]| {
        top.iter().map(|&i| u[i][2]).sum::<f64>() / top.len() as f64
            - bottom.iter().map(|&i| u[i][2]).sum::<f64>() / bottom.len() as f64
    };
    let dt = 2e-8;
    let mut supplied = 0.0;
    println!("Transient: 400 ns linear 0-to-100 V ramp; native consistent mass, dt=20 ns.");
    for i in 1..=20 {
        let candidate = body.step_candidate(
            dt,
            &[
                ElectrodeCondition::Voltage(0.0),
                ElectrodeCondition::Voltage(-5.0 * f64::from(i)),
            ],
            zeros,
            zero_charge,
            controls,
            cx,
        )?;
        supplied += candidate.step().electrical_work_j + candidate.step().mechanical_work_j;
        body = candidate.into_parts().0;
    }
    println!(
        "End of ramp: thickness change={:.6} nm, kinetic energy={:.9e} J, accumulated work defect={:.3e} J",
        thickness(body.state().displacement_m()) * 1e9,
        body.state().kinetic_energy_j(),
        body.state().total_energy_j() - supplied,
    );

    // Rebind with the actual accepted total charge. This is a matched open
    // circuit, not a projection to zero charge or a sudden short circuit.
    let charge = body.state().electrodes()[1].charge_c;
    let open = [
        ElectrodeCondition::Voltage(0.0),
        ElectrodeCondition::Charge(charge),
    ];
    let mut isolated = TetPiezoelectricTransient::new(
        model,
        supports,
        &[
            PiezoelectricElectrode {
                nodes: bottom,
                condition: open[0],
            },
            PiezoelectricElectrode {
                nodes: top,
                condition: open[1],
            },
        ],
        body.state().displacement_m(),
        body.state().velocity_m_s(),
        zero_charge,
        controls,
        cx,
    )?;
    let initial_energy = isolated.state().total_energy_j();
    if (initial_energy - body.state().total_energy_j()).abs() > controls.energy_tolerance_j {
        return Err("charge-preserving electrode disconnection changed stored energy".into());
    }
    let mut loaded = isolated.clone();
    let load = declared_load()?;
    let drives = [ElectrodeDrive::Voltage(0.0), ElectrodeDrive::Current(0.0)];
    let mut resistor_heat = 0.0;
    println!(
        "Isolated sensor: time since disconnection [us], thickness change [nm], voltage [V], energy drift [J]"
    );
    for i in 1..=100 {
        isolated = isolated
            .step_candidate(dt, &open, zeros, zero_charge, controls, cx)?
            .into_parts()
            .0;
        let (candidate, report) = loaded.step_with_loads_candidate(
            dt,
            &drives,
            std::slice::from_ref(&load),
            zeros,
            zero_charge,
            controls,
            cx,
        )?;
        resistor_heat += report.joule_heat_j;
        loaded = candidate.into_parts().0;
        if i % 20 == 0 {
            let state = isolated.state();
            println!(
                "{:.3}, {:.6}, {:.9}, {:.3e}",
                state.time_s() * 1e6,
                thickness(state.displacement_m()) * 1e9,
                state.electrodes()[1].voltage_v,
                state.total_energy_j() - initial_energy
            );
        }
    }
    println!(
        "Same initial state with authored 1 kohm load: thickness change={:.6} nm, voltage={:.9} V, resistor heat={:.9e} J, energy defect={:.3e} J",
        thickness(loaded.state().displacement_m()) * 1e9,
        loaded.state().electrodes()[1].voltage_v,
        resistor_heat,
        loaded.state().total_energy_j() + resistor_heat - initial_energy,
    );
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mode = match std::env::args().skip(1).collect::<Vec<_>>().as_slice() {
        [] => "default",
        [mode] if mode == "--harmonic" => "harmonic",
        [mode] if mode == "--acoustic" => "acoustic",
        _ => return Err("usage: piezoelectric_plate [--harmonic|--acoustic]".into()),
    };
    let material = nominal_material()?;
    let sources = [material];
    let (mesh, mut nodes) = fs_feec::fixtures::kuhn_cube(2);
    for p in &mut nodes {
        p[0] *= 0.01;
        p[1] *= 0.01;
        p[2] *= 0.002;
    }
    let tets = mesh
        .tets
        .iter()
        .map(|t| t.map(|v| v as usize))
        .collect::<Vec<_>>();
    let bottom = nodes
        .iter()
        .enumerate()
        .filter_map(|(i, p)| (p[2] == 0.0).then_some(i))
        .collect::<Vec<_>>();
    let top = nodes
        .iter()
        .enumerate()
        .filter_map(|(i, p)| (p[2] == 0.002).then_some(i))
        .collect::<Vec<_>>();
    let origin = nodes.iter().position(|p| *p == [0.0, 0.0, 0.0]).unwrap();
    let along_x = nodes.iter().position(|p| *p == [0.01, 0.0, 0.0]).unwrap();
    let along_y = nodes.iter().position(|p| *p == [0.0, 0.01, 0.0]).unwrap();
    // Six rigid-motion supports, allowing all free d31/d33 strains.
    let supports = [
        (3 * origin, 0.0),
        (3 * origin + 1, 0.0),
        (3 * origin + 2, 0.0),
        (3 * along_x + 1, 0.0),
        (3 * along_x + 2, 0.0),
        (3 * along_y + 2, 0.0),
    ];
    let zero_forces = vec![[0.0; 3]; nodes.len()];
    let zero_charge = vec![0.0; nodes.len()];
    let controls = TetPiezoelectricControls {
        relative_tolerance: 1e-11,
        maximum_iterations: 4000,
        force_tolerance_n: 1e-6,
        charge_tolerance_c: 1e-15,
        energy_tolerance_j: 1e-12,
        maximum_strain_norm: 1e-3,
        maximum_electric_field_v_m: 1e5,
    };
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena|->Result<(),Box<dyn std::error::Error>> {
        let cx=Cx::new(&gate,arena,fs_exec::StreamKey {seed:59,kernel_id:59,tile:0,iteration:0},
            fs_exec::Budget::INFINITE,fs_exec::ExecMode::Deterministic);
        let model=TetPiezoelectricModel::assemble(&nodes,&tets,&sources,FRAME,TetAssemblyBudget::standard(),&cx)?;
        if mode == "harmonic" {
            return harmonic_response(&model,&supports,&bottom,&top,nodes.len(),controls,&cx);
        }
        if mode == "acoustic" {
            return acoustic_response(&model,&supports,&bottom,&top,&tets,&zero_forces,&zero_charge,controls,&cx);
        }
        let driven=model.solve(&supports,&[
            PiezoelectricElectrode {nodes:&bottom,condition:ElectrodeCondition::Voltage(0.0)},
            PiezoelectricElectrode {nodes:&top,condition:ElectrodeCondition::Voltage(-100.0)}],
            &zero_forces,&zero_charge,controls,&cx)?;
        let extension=top.iter().map(|&i|driven.displacement_m[i][2]).sum::<f64>()/top.len()as f64;
        let charge=driven.electrodes[1].charge_c;
        let isolated=model.solve(&supports,&[
            PiezoelectricElectrode {nodes:&bottom,condition:ElectrodeCondition::Voltage(0.0)},
            PiezoelectricElectrode {nodes:&top,condition:ElectrodeCondition::Charge(charge)}],
            &zero_forces,&zero_charge,controls,&cx)?;
        println!("Nominal PZT-5A: 10 x 10 x 2 mm, {} tetrahedra",tets.len());
        println!("Source: {SOURCE}");
        println!("Authored frozen, lossless, transversely isotropic approximation; no measured-device accuracy claim.");
        println!("100 V drive: thickness extension = {:.6} nm (nominal d33 oracle: 39 nm)",extension*1e9);
        println!("Top electrode charge = {:.6} nC; free capacitance = {:.6} nF",charge*1e9,-charge/100.0*1e9);
        println!("Isolated same-charge electrode: solved voltage = {:.9} V",isolated.electrodes[1].voltage_v);
        println!("Stored energy = {:.9e} J; source-work residual = {:.3e} J",driven.stored_energy_j,driven.energy_residual_j);
        println!("MINRES iterations = {}; force residual = {:.3e} N; charge residual = {:.3e} C",
            driven.iterations,driven.maximum_force_residual_n,driven.maximum_charge_residual_c);
        for condition in [PiezoelectricWaveCondition::FixedField,PiezoelectricWaveCondition::ChargeFree] {
            let waves=sources[0].directional_waves([0.0,0.0,1.0],condition,1e-3)?;
            println!("Frozen local {condition:?} phase speeds [m/s]: {:?}",waves.mechanical.modes.map(|m|m.phase_speed_m_s));
            println!("Frozen local {condition:?} energy velocities [m/s]: {:?}",waves.mechanical.modes.map(|m|m.energy_velocity_m_s));
        }
        transient_response(&model,&supports,&bottom,&top,&zero_forces,&zero_charge,controls,&cx)?;
        Ok(())
    })
}
