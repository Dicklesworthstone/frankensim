//! Causal material memory must preserve resolved nonlinear force and work.
//! These are numerical checks of the existing Maxwell owner, not a new loss law.

use fs_nlmodal::{SosModalStorage, StressChannel, assemble_storage_with_damping};
use fs_phs::linear_exact::{LinearSplit, LinearSplitBudget, split_quadratic, step_render};
use fs_phs::{PortHamiltonian, RelaxationBranch, Storage};

struct EnergyDatum {
    base: SosModalStorage,
    datum: f64,
}

impl Storage for EnergyDatum {
    fn hamiltonian(&self, state: &[f64]) -> f64 {
        self.datum + self.base.hamiltonian(state)
    }

    fn gradient(&self, state: &[f64], output: &mut [f64]) {
        self.base.gradient(state, output);
    }

    fn quadratic_part(&self) -> Option<fs_phs::linear_exact::QuadraticPart> {
        self.base.quadratic_part()
    }

    fn midpoint_energy_remainder(&self, a: &[f64], b: &[f64], mid: &[f64]) -> Option<f64> {
        self.base.midpoint_energy_remainder(a, b, mid)
    }
}

fn base(datum: f64) -> PortHamiltonian {
    let storage = EnergyDatum {
        base: SosModalStorage {
            omegas: vec![5.0, 7.0],
            channels: vec![StressChannel {
                coefficient: 1.4,
                coupling: vec![1.0, 0.3, 0.3, 0.8],
            }],
        },
        datum,
    };
    assemble_storage_with_damping(
        2,
        &[5.0, 7.0],
        &[0.3, -0.1, -0.1, 0.2],
        1,
        vec![0.0, 0.8, 0.0, -0.25],
        Box::new(storage),
    )
    .unwrap()
}

fn arms() -> Vec<RelaxationBranch> {
    vec![
        RelaxationBranch {
            projection: vec![2.0, 0.0, -1.0, 0.0],
            stiffness: 1.0,
            relaxation_time_s: 0.12,
        },
        RelaxationBranch {
            projection: vec![-0.6, 0.0, 1.3, 0.0],
            stiffness: 1.0,
            relaxation_time_s: 0.04,
        },
    ]
}

fn split(system: &PortHamiltonian) -> LinearSplit {
    split_quadratic(
        system,
        LinearSplitBudget {
            max_modes: 6,
            max_omega_dt: 3.0,
        },
    )
    .unwrap()
}

#[test]
fn maxwell_memory_preserves_nonlinear_response_under_a_large_energy_datum() {
    let reference = base(0.0).with_relaxation_branches(arms()).unwrap();
    let shifted = base(1.0e24).with_relaxation_branches(arms()).unwrap();
    let reference_split = split(&reference);
    let shifted_split = split(&shifted);
    let initial = [0.3, 0.2, -0.15, 0.08, -0.2, 0.11];
    let mut a = initial.to_vec();
    let mut b = initial.to_vec();
    let mut loss = 0.0;
    for _ in 0..40 {
        let ra = step_render(&reference, &reference_split, &a, &[0.4], 0.001).unwrap();
        let rb = step_render(&shifted, &shifted_split, &b, &[0.4], 0.001).unwrap();
        // An energy datum changes neither motion nor constitutive memory. The
        // complete endpoint H subtraction cannot resolve ANY work here.
        assert_eq!(shifted.hamiltonian(&b), shifted.hamiltonian(&rb.x));
        assert_eq!(ra.x, rb.x);
        assert_eq!(ra.delta_h, rb.delta_h);
        assert_eq!(ra.dissipated, rb.dissipated);
        assert_eq!(ra.supplied, rb.supplied);
        assert!(ra.balance_residual().abs() < 2.0e-12);
        assert!(ra.dissipated >= 0.0);
        loss += ra.dissipated;
        a = ra.x;
        b = rb.x;
    }
    assert!(loss > 0.05, "the material histories must actually relax");
    assert!((a[1] - initial[1]).abs() > 0.1);
    assert!((a[4] - initial[4]).abs() > 0.1);
}

// Independent polarization of the original modal and material energies. No
// subtraction of total H and no reference to the discrete-gradient formula.
fn energy_change(a: &[f64], b: &[f64]) -> f64 {
    let mut change = 0.0;
    for (mode, omega2) in [25.0, 49.0].iter().enumerate() {
        let q = 2 * mode;
        let p = q + 1;
        change +=
            0.5 * (a[p] + b[p]) * (b[p] - a[p]) + 0.5 * omega2 * (a[q] + b[q]) * (b[q] - a[q]);
    }
    let form = |x: &[f64]| x[0] * x[0] + 0.6 * x[0] * x[2] + 0.8 * x[2] * x[2];
    let form_change = (a[0] + b[0]) * (b[0] - a[0])
        + 0.3 * ((a[0] + b[0]) * (b[2] - a[2]) + (a[2] + b[2]) * (b[0] - a[0]))
        + 0.8 * (a[2] + b[2]) * (b[2] - a[2]);
    change += 0.35 * (form(a) + form(b)) * form_change;
    for (arm, row) in [[2.0, -1.0], [-0.6, 1.3]].iter().enumerate() {
        let before = row[0] * a[0] + row[1] * a[2] - a[4 + arm];
        let after = row[0] * b[0] + row[1] * b[2] - b[4 + arm];
        let extension_change =
            row[0] * (b[0] - a[0]) + row[1] * (b[2] - a[2]) - (b[4 + arm] - a[4 + arm]);
        change += 0.5 * (before + after) * extension_change;
    }
    change
}

#[test]
fn nested_maxwell_storage_keeps_every_arm_and_its_independent_work() {
    let combined = base(0.0).with_relaxation_branches(arms()).unwrap();
    let mut branches = arms();
    let mut second = branches.pop().unwrap();
    // The second attachment sees the first memory coordinate, whose zero
    // projection is explicit. Both paths describe the same original mechanics.
    second.projection.push(0.0);
    let nested = base(1.0e24)
        .with_relaxation_branches(branches)
        .unwrap()
        .with_relaxation_branches(vec![second])
        .unwrap();
    let combined_split = split(&combined);
    let nested_split = split(&nested);
    let mut a = vec![0.3, 0.2, -0.15, 0.08, -0.2, 0.11];
    let mut b = a.clone();
    let dt = 0.001;
    for sample in 0..40 {
        let force = if sample < 20 { 0.4 } else { -0.2 };
        let ra = step_render(&combined, &combined_split, &a, &[force], dt).unwrap();
        let rb = step_render(&nested, &nested_split, &b, &[force], dt).unwrap();
        for (one, two) in ra.x.iter().zip(&rb.x) {
            assert!((one - two).abs() < 2.0e-12);
        }
        let displacement = [rb.x[0] - b[0], rb.x[2] - b[2]];
        let velocity = [displacement[0] / dt, displacement[1] / dt];
        let work = force * (0.8 * displacement[0] - 0.25 * displacement[1]);
        let viscous_power = 0.3 * velocity[0] * velocity[0] - 0.2 * velocity[0] * velocity[1]
            + 0.2 * velocity[1] * velocity[1];
        // z' = -effort_z/tau; therefore the accepted Maxwell heat is
        // sum tau*(z1-z0)^2/dt, separately from Kelvin velocity damping.
        let heat = dt * viscous_power
            + 0.12 * (rb.x[4] - b[4]).powi(2) / dt
            + 0.04 * (rb.x[5] - b[5]).powi(2) / dt;
        let stored_change = energy_change(&b, &rb.x);
        assert!((rb.delta_h - stored_change).abs() < 2.0e-12);
        assert!((rb.dissipated - heat).abs() < 2.0e-12);
        assert!((rb.supplied - work).abs() < 2.0e-12);
        assert!((stored_change + heat - work).abs() < 2.0e-12);
        a = ra.x;
        b = rb.x;
    }
}
