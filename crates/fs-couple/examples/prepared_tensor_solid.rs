//! Complete TET10 motion with consistent inertia and original tensor memory.
//!
//! Run `cargo run -p fs-couple --example prepared_tensor_solid --release > motion.csv`.
//! Standard output is CSV; source/capacity information and elapsed times go to
//! standard error. All 30 nodal coordinates are retained by an identity basis;
//! no modal reduction, mass lumping, pole fitting or empirical damping is used.
//!
//! This is an authored anisotropic material fixture, not measured material data
//! or a mesh/time convergence qualification. The reference tetrahedron is free
//! to translate and rotate within small-strain kinematics. Its initial affine
//! extension and unrelaxed Maxwell memory are specified initial conditions.
//! A balanced axial force pair acts for 24 intervals, followed by free release.
//! The elapsed times include this example's observations and CSV formatting;
//! they are not a performance benchmark or a real-time guarantee.

use std::io::{BufWriter, Write};
use std::time::Instant;

use fs_couple::tensor_material::{
    InitialMemory, TensorMaterialBody, TensorMaterialBudget, prepared::PreparedTensorControls,
    solid::prepared::PreparedTensorSolidBody,
};
use fs_exec::{CancelGate, Cx};
use fs_material::visco::tensor::{Mandel6, PronyBranch, TensorProny};
use fs_solid::linear3::{
    TetAssemblyBudget,
    tet10::{TetQuadrature, promote_p1_to_p2},
    visco::{
        SolidMeshRef, SolidProjectionBudget, TetViscoField, TetViscoMaterial, ViscoElements,
        ViscoStateBudget,
    },
};

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

fn diagonal(values: [f64; 6]) -> Mandel6 {
    core::array::from_fn(|i| core::array::from_fn(|j| if i == j { values[i] } else { 0.0 }))
}

fn material() -> TetViscoMaterial {
    let equilibrium = diagonal([120e6, 90e6, 120e6, 60e6, 50e6, 40e6]);
    // Distinct directional ratios make Kelvin and Maxwell response part of
    // the coupled tensor problem, rather than one authored modal loss factor.
    let maxwell = diagonal([48e6, 27e6, 60e6, 21e6, 22.5e6, 14e6]);
    let kelvin = diagonal([12000.0, 14000.0, 8000.0, 6500.0, 3500.0, 6000.0]);
    TetViscoMaterial {
        law: TensorProny::new_mandel(
            equilibrium,
            vec![PronyBranch::new(maxwell, 0.004).unwrap()],
            0.01,
        )
        .unwrap()
        .with_kelvin_viscosity_mandel(kelvin)
        .unwrap(),
        density_kg_m3: 1000.0,
        material_state_identity: [97; 32],
    }
}

fn controls() -> PreparedTensorControls {
    PreparedTensorControls {
        maximum_coordinates: 30,
        maximum_material_points: 64,
        maximum_state_coordinates: 90,
        maximum_heat_groups: 1,
        maximum_preparation_bytes: 512 * 1024 * 1024,
        maximum_preparation_multiply_terms: 1_000_000_000_000,
        maximum_step_multiply_terms: 100_000_000,
        maximum_relative_operator_error: 1e-9,
        maximum_work_error_j: 1e-10,
        maximum_steps: 480,
    }
}

fn run(cx: &Cx<'_>, gate: &CancelGate) -> Result<()> {
    let cold_start = Instant::now();
    let corners = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 0.2, 0.0],
        [0.0, 0.0, 0.2],
    ];
    let (nodes, cells) = promote_p1_to_p2(&corners, &[[0, 1, 2, 3]], None, cx)?;
    let material = material();
    let source = ViscoElements::new(
        SolidMeshRef::Tet10 {
            nodes_m: &nodes,
            elements: &cells,
            quadrature: TetQuadrature {
                stiffness_points_per_dir: 3,
                mass_points_per_dir: 4,
            },
        },
        TetViscoField::Uniform(&material),
        &[],
        TetAssemblyBudget::standard(),
        ViscoStateBudget {
            max_internal_variables: 1000,
        },
        cx,
    )?;
    let coordinates = 3 * nodes.len();
    let identity: Vec<Vec<f64>> = (0..coordinates)
        .map(|column| {
            (0..coordinates)
                .map(|row| if row == column { 1.0 } else { 0.0 })
                .collect()
        })
        .collect();
    let projection = source.project_basis(
        &identity,
        SolidProjectionBudget {
            max_coordinates: 30,
            max_basis_scalars: 900,
            max_material_points: 64,
            max_scalar_terms: 20_000_000,
        },
        cx,
    )?;
    let point_count = projection.material_points().len();
    let original_memory_rows: usize = projection
        .material_points()
        .iter()
        .flat_map(|point| point.branches())
        .map(|branch| branch.rows().len())
        .sum();
    let original_states = 2 * coordinates + original_memory_rows;
    let original_block = 2 * (original_states + coordinates + 1);
    match TensorMaterialBody::from_mass_matrix(
        coordinates,
        projection.mass_matrix_kg().to_vec(),
        projection.material_points().to_vec(),
        TensorMaterialBudget {
            max_state_coordinates: 256,
            max_material_points: 64,
        },
        gate,
    ) {
        Err(error) => eprintln!("original_full_memory_refusal: {error}"),
        Ok(_) => return Err("the original full-memory capacity refusal was not reproduced".into()),
    }
    let dt = 1.0 / 48000.0;
    let body = PreparedTensorSolidBody::new(projection, dt, controls(), gate)?;
    if coordinates != 30 || point_count != 27 || body.material_body().state_dimension() != 90 {
        return Err(
            "the complete TET10 fixture did not retain its expected source dimensions".into(),
        );
    }
    let initial: Vec<f64> = nodes
        .iter()
        .flat_map(|position| [1e-5 * position[0], 0.0, 0.0])
        .collect();
    let mut history = body.history(&initial, &[0.0; 30], InitialMemory::Instantaneous, gate)?;
    let mut stiffness_times_initial = vec![0.0; coordinates];
    source
        .equilibrium()
        .stiffness
        .spmv(&initial, &mut stiffness_times_initial);
    let initial_equilibrium = 0.5
        * initial
            .iter()
            .zip(stiffness_times_initial)
            .map(|(u, ku)| u * ku)
            .sum::<f64>();
    let initial_memory: f64 = body
        .projection()
        .material_points()
        .iter()
        .flat_map(|point| point.branches())
        .flat_map(|branch| branch.rows())
        .map(|row| {
            let spring = row.iter().zip(&initial).map(|(a, u)| a * u).sum::<f64>();
            0.5 * spring * spring
        })
        .sum();
    eprintln!(
        "synthetic TET10: coordinates={coordinates}, material_points={point_count}, original_states={original_states}, original_exponential_block={original_block}, prepared_states={}, relaxation_time_s=0.004, reference_mass_kg={:.9e}",
        body.material_body().state_dimension(),
        source.equilibrium().total_mass_kg,
    );
    eprintln!(
        "source_and_flow_preparation_ms={:.3}; force_pair=+/-0.1N for 24 intervals; dt_s={dt:.9e}; tip_extension=x_displacement(corner1)-x_displacement(corner0)",
        cold_start.elapsed().as_secs_f64() * 1000.0,
    );

    let stdout = std::io::stdout();
    let mut csv = BufWriter::new(stdout.lock());
    writeln!(
        csv,
        "step,time_s,tip_extension_m,equilibrium_energy_j,kinetic_energy_j,maxwell_memory_energy_j,cumulative_external_work_j,cumulative_kelvin_maxwell_heat_j,interval_external_work_j,interval_heat_j,interval_balance_error_j,cumulative_balance_error_j"
    )?;
    writeln!(
        csv,
        "0,0.000000000000e0,{:.12e},{initial_equilibrium:.12e},0.000000000000e0,{initial_memory:.12e},0.000000000000e0,0.000000000000e0,0.000000000000e0,0.000000000000e0,0.000000000000e0,0.000000000000e0",
        initial[3] - initial[0]
    )?;
    let warm_start = Instant::now();
    let mut maximum_balance = 0.0_f64;
    let mut maximum_nodal_work_error = 0.0_f64;
    for interval in 0..480 {
        let mut force = [[0.0; 3]; 10];
        if interval < 24 {
            force[0][0] = -0.1;
            force[1][0] = 0.1;
        }
        let step = history.advance(&force, gate)?;
        let displacement = history.nodal_displacement_m(gate)?;
        let energy = step.material.energy;
        maximum_balance = maximum_balance.max(energy.cumulative_balance_error_j.abs());
        maximum_nodal_work_error = maximum_nodal_work_error.max(step.nodal_work_error_j.abs());
        writeln!(
            csv,
            "{},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e}",
            step.material.step,
            step.material.time_s,
            displacement[1][0] - displacement[0][0],
            energy.equilibrium_after_j,
            energy.kinetic_after_j,
            energy.memory_after_j,
            energy.cumulative_work_j,
            energy.cumulative_dissipated_j,
            energy.supplied_work_j,
            energy.dissipated_j,
            energy.balance_error_j,
            energy.cumulative_balance_error_j,
        )?;
    }
    csv.flush()?;
    eprintln!(
        "completed_steps=480, observed_run_and_csv_ms={:.3}, maximum_cumulative_balance_error_j={maximum_balance:.9e}, maximum_nodal_work_error_j={maximum_nodal_work_error:.9e}",
        warm_start.elapsed().as_secs_f64() * 1000.0,
    );
    Ok(())
}

fn main() -> Result<()> {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        run(
            &Cx::new(
                &gate,
                arena,
                fs_exec::StreamKey {
                    seed: 97,
                    kernel_id: 44,
                    tile: 0,
                    iteration: 0,
                },
                fs_exec::Budget::INFINITE,
                fs_exec::ExecMode::Deterministic,
            ),
            &gate,
        )
    })
}
