//! Source stiffness/inertia/viscosity -> native frame -> prepared coupled motion.
//! `cargo run -p fs-couple --example frame_vibration > frame.csv`
//! Optional arguments: density_kg_m3 young_modulus_pa axial_eta_pa_s bending_eta_pa_s.
//! All four are manufactured constant source values, not measured material claims.
//! Geometry stays a 1 m cantilever with eight members and a 0.02 m solid radius.
use fs_couple::frame_vibration::{FrameVibrationBasis, FrameVibrationBudget};
use fs_couple::modal_acoustic_time::ModalAcousticTimeBudget;
use fs_couple::render::schedule::force::coupled::{ModalCouplingConfig, PreparedModalBudget};
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use fs_qty::{Density, DynViscosity, Pressure};
use fs_solid::frame2d::{Frame2d, FrameMassFormulation, FrameViscousSection, Geometry};
use std::io::{BufWriter, Write};

type Error = Box<dyn std::error::Error>;

fn section(values: [f64; 4]) -> Result<FrameViscousSection, Error> {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for ((name, dims), value) in [
        ("density", Density::DIMS),
        ("young_modulus", Pressure::DIMS),
        ("axial_eta", DynViscosity::DIMS),
        ("bending_eta", DynViscosity::DIMS),
    ]
    .into_iter()
    .zip(values)
    {
        let mut validity = ValidityDomain::unconstrained().with("T", 300.0, 300.0);
        if name.ends_with("_eta") {
            validity = validity.with("omega", 0.0, 1e6);
        }
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Scalar { value, dims },
            validity,
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: "explicit manufactured frame example input".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })?;
        requirements.push(ScalarPropertyRequirement::try_new(
            name,
            dims,
            ScalarAdmissibility::Finite,
        )?);
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured solid".into(),
            phase: "solid".into(),
            process: "constant linear elastic and viscous example".into(),
            revision: 0,
        },
        claims,
        vec![],
    )?;
    let material = resolve_material_state_point(
        &card,
        &QueryPoint::new().with("T", 300.0)?.with("omega", 0.0)?,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )?;
    let elastic = ConcentricElasticSection::try_new(
        0.0,
        vec![ConcentricMaterialLayer {
            outer_radius_m: 0.02,
            material,
        }],
    )?;
    Ok(FrameViscousSection::try_new(
        elastic,
        "axial_eta",
        "bending_eta",
    )?)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    let mut values = [1000.0, 2e9, 2e5, 4e6];
    if !args.is_empty() {
        if args.len() != 4 {
            return Err("usage: frame_vibration [density_kg_m3 young_modulus_pa axial_eta_pa_s bending_eta_pa_s]".into());
        }
        for (value, argument) in values.iter_mut().zip(args) {
            *value = argument.parse()?;
        }
    }
    let section = section(values)?;
    let mut frame =
        Frame2d::with_mass_formulation(Geometry::Linear, FrameMassFormulation::ConsistentReference);
    let elements = 8;
    for i in 0..=elements {
        frame.add_node(i as f64 / elements as f64, 0.0);
    }
    frame.fix(0, [true; 3]);
    for i in 0..elements {
        frame.add_material_viscous_element(i, i + 1, 3, section.clone())?;
    }
    // Resolve the COMPLETE native spectrum before choosing the sampled clock;
    // six retained modes cannot conceal an unresolved high-frequency member mode.
    let maximum_omega = frame
        .modal_analysis()?
        .modes
        .iter()
        .filter_map(|mode| mode.angular_frequency())
        .fold(0.0_f64, f64::max);
    let rate = (2.0 * maximum_omega / std::f64::consts::PI)
        .ceil()
        .max(4000.0);
    if !rate.is_finite() || rate > 192_000.0 {
        return Err("complete frame spectrum exceeds this example's 192 kHz work budget".into());
    }
    let rate = rate as u32;
    let basis = FrameVibrationBasis::from_frame(
        &frame,
        6,
        FrameVibrationBudget {
            max_frame_dofs: 27,
            max_dense_entries: 729,
        },
        || false,
    )?;
    let component = basis.coupled_component(
        rate,
        &[0.0; 6],
        ModalAcousticTimeBudget::audible_reference(),
        6 * 6 * 3 * elements,
        || false,
    )?;
    let config = ModalCouplingConfig {
        max_modes: 6,
        max_connections: 3 * elements,
        max_setup_terms: 10_000,
        nyquist_guard_fraction: 0.9,
        maximum_total_energy_j: 1e3,
        maximum_abs_pressure_pa: 1.0,
        maximum_abs_connection_force_n: 1e5,
        solve_relative_tolerance: 1e-10,
        energy_absolute_tolerance_j: 1e-9,
        energy_relative_tolerance: 1e-8,
    };
    let gate = CancelGate::new();
    let (basis, mut system) = component.into_system(
        config,
        Some(PreparedModalBudget {
            max_modes: 6,
            max_dense_bytes: 4 * 1024 * 1024,
            max_multiply_terms: 200_000_000,
        }),
        &gate,
    )?;
    eprintln!(
        "manufactured rho/E/eta_axial/eta_bending={values:?}; {rate} Hz; {} of {} modes; vertical effective/full mass={}/{} kg",
        basis.modes().len(),
        basis.available_mode_count(),
        basis.retained_effective_mass_kg()[1],
        basis.base_influence_mass_kg()[1]
    );
    let mut load = vec![0.0; frame.ndof()];
    let mut generalized = [0.0; 6];
    let (mut displacement, mut velocity) = (vec![0.0; frame.ndof()], vec![0.0; frame.ndof()]);
    let (mut work, mut material_loss) = (0.0, 0.0);
    let stdout = std::io::stdout();
    let mut out = BufWriter::new(stdout.lock());
    writeln!(
        out,
        "time_s,tip_y_m,tip_vy_m_s,energy_j,cumulative_nodal_work_j,cumulative_material_loss_j,step_energy_residual_j"
    )?;
    for sample in 0..2 * rate {
        // One newton held at the tip for the first 20 ms; then free decay.
        load[3 * elements + 1] = if f64::from(sample) / f64::from(rate) < 0.02 {
            1.0
        } else {
            0.0
        };
        basis.project_force_into(&load, [0.0; 2], &mut generalized, || false)?;
        let report = system.step_under_gate(&generalized, &gate)?;
        work += report.external_work_j;
        material_loss += report.connection_dissipation_j;
        let (energy, residual) = (report.modal_energy_j, report.energy_residual_j);
        basis.reconstruct_into(
            system.components()[0].states(),
            &mut displacement,
            &mut velocity,
            || false,
        )?;
        writeln!(
            out,
            "{:.9},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e}",
            f64::from(sample + 1) / f64::from(rate),
            displacement[3 * elements + 1],
            velocity[3 * elements + 1],
            energy,
            work,
            material_loss,
            residual
        )?;
    }
    out.flush()?;
    Ok(())
}
