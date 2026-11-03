//! Nominal PZT-5A plate: the actual 3D source/tensor/electrode solve.
//!
//! Run: cargo run -p fs-couple --example piezoelectric_plate
//! Data and approximation boundaries: examples/PIEZOELECTRIC_PLATE.md.

use fs_blake3::ContentHash;
use fs_couple::piezoelectric::spatial::{
    ElectrodeCondition, PiezoelectricElectrode, TetPiezoelectricControls, TetPiezoelectricModel,
};
use fs_evidence::ValidityDomain;
use fs_exec::{CancelGate, Cx};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::piezoelectric::tensor::{TensorPiezoelectric, waves::PiezoelectricWaveCondition};
use fs_material::state_point::MaterialPropertySelection;
use fs_solid::linear3::TetAssemblyBudget;

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

fn main() -> Result<(), Box<dyn std::error::Error>> {
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
        Ok(())
    })
}
