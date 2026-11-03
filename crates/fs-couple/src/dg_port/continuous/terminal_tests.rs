//! Complete original DG/material states, mixed terminal constraints and work.
use super::*;
use fs_duct::dg::{DgMesh, DgSidePort};
use fs_duct::{Duct, Segment};
use fs_exec::CancelGate;
use fs_material::gas::{GasSpec, GasState};
use fs_math::c64::C64;
use fs_phs::linear_exact::prepared::PreparedQuadraticFlow;
use fs_phs::linear_exact::{DenseExponentialBudget, LinearSplitBudget, split_quadratic};
use fs_phs::{PortHamiltonian, QuadraticStorage, WallPin};
use fs_vfit::mimo::{MimoModel, load::MimoLoadLimits};
use fs_vfit::model::PoleTerm;

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn medium() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap()
}
fn network(wall: bool) -> DgContinuousNetwork {
    let mut duct = DgDuct::new(
        &Duct {
            segments: vec![Segment::Cylinder {
                radius: 0.01,
                length: 0.12,
            }],
        },
        &medium(),
        DgMesh {
            elements: 2,
            order: 2,
        },
        DgBoundary::Port(PortId(0)),
        DgBoundary::Port(PortId(1)),
    )
    .unwrap()
    .with_side_ports(&[DgSidePort {
        id: PortId(2),
        downstream_element: 1,
    }])
    .unwrap();
    if wall {
        duct = duct
            .with_compliant_wall(WallPin {
                surface_density: 0.2,
                stiffness_per_area: 1e8,
                resistance: 350.0,
            })
            .unwrap();
    }
    let view = duct.continuous_operator();
    let trace = |id| {
        view.ports()
            .iter()
            .find(|(p, _)| *p == PortId(id))
            .unwrap()
            .1
            .clone()
    };
    let n = view.dimension();
    let a = view.generator_per_s().to_vec();
    let mut loss = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            loss[i * n + j] = -0.5 * (a[i * n + j] + a[j * n + i]);
        }
    }
    DgContinuousNetwork {
        dimension: n,
        generator: a,
        inlet_column: trace(0),
        inlet_pressure_column: trace(0),
        terminal_pressure_columns: vec![trace(1), trace(2)],
        terminal_input_columns: vec![trace(1), trace(2)],
        owners: vec![network::DgContinuousOwner {
            offset: 0,
            dimension: n,
            loss,
        }],
        junctions: vec![],
    }
}

fn source(kind: usize) -> MimoModel {
    let impedance =
        medium().density * medium().sound_speed / (core::f64::consts::PI * 0.01_f64.powi(2));
    let d = match kind {
        0 => vec![0.0; 4],
        // Rank one, with signed mutual coupling. Memory uses different axes.
        1 => vec![0.25, -0.5, -0.5, 1.0],
        _ => vec![0.7, -0.15, -0.15, 0.4],
    };
    MimoModel {
        ports: 2,
        poles: vec![PoleTerm::Real {
            pole: -400.0,
            residue: 0.0,
        }],
        residues: vec![
            [0.8, 0.18, 0.18, 0.6]
                .iter()
                .map(|r| C64::from_re(r * impedance * 400.0))
                .collect(),
        ],
        d: d.iter().map(|d| d * impedance).collect(),
        e: vec![0.0; 4],
        reciprocal: true,
    }
}
fn load(model: &MimoModel) -> PreparedMimoLoad {
    PreparedMimoLoad::new(
        model,
        1e-5,
        MimoLoadLimits {
            band_rad_s: [1.0, 2e6],
            sampled_passivity_points: 33,
            max_states: 16,
            max_memory_bytes: 1 << 20,
        },
    )
    .unwrap()
}
fn sym_product(a: &[f64], b: &[f64]) -> Vec<f64> {
    let n = a.len();
    let mut out = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            out[i * n + j] = 0.5 * (a[i] * b[j] + b[i] * a[j]);
        }
    }
    out
}
fn extend(row: &[f64], last: f64) -> Vec<f64> {
    let mut out = row.to_vec();
    out.push(last);
    out
}

// The production exact-flow owner integrates the ORIGINAL joint generator and
// independently authored physical work forms. The held inlet is explicit and
// stores no energy. The fixture's relaxation law has known H=z^T K z/2,
// D=400*z^T K z+U^T D U; this does not infer storage for arbitrary rational fits.
fn flow(
    original: &DgContinuousNetwork,
    terminal: &DgContinuousTerminal,
    model: &MimoModel,
    dt: f64,
) -> PreparedQuadraticFlow {
    let nd = original.dimension;
    let n = terminal.inlet_column.len();
    let total = n + 1;
    let mut j = vec![0.0; nd * nd];
    let mut r = j.clone();
    let mut q = j.clone();
    for i in 0..nd {
        q[i * nd + i] = 1.0;
        for k in 0..nd {
            j[i * nd + k] = 0.5 * (original.generator[i * nd + k] - original.generator[k * nd + i]);
            r[i * nd + k] =
                -0.5 * (original.generator[i * nd + k] + original.generator[k * nd + i]);
        }
    }
    let system = PortHamiltonian::new(
        nd,
        1,
        j,
        r,
        original.inlet_column.clone(),
        Box::new(QuadraticStorage::new(q, nd).unwrap()),
    )
    .unwrap();
    let split = split_quadratic(
        &system,
        LinearSplitBudget {
            max_modes: 128,
            max_omega_dt: 1.0,
        },
    )
    .unwrap();
    let mut generator = vec![0.0; total * total];
    for i in 0..n {
        generator[i * total..i * total + n]
            .copy_from_slice(&terminal.generator[i * n..(i + 1) * n]);
        generator[i * total + n] = terminal.inlet_column[i];
    }
    let mut inlet_pressure = vec![0.0; total];
    inlet_pressure[..nd].copy_from_slice(&original.inlet_pressure_column);
    let mut held = vec![0.0; total];
    held[n] = 1.0;
    let mut forms = vec![sym_product(&inlet_pressure, &held)];
    let mut linear = Vec::new();
    let mut fluxes = Vec::new();
    for p in 0..2 {
        let flux = extend(&terminal.flow_columns[p], terminal.inlet_flow_factors[p]);
        let pressure = extend(&terminal.pressure_columns[p], 0.0);
        let load_pressure = extend(
            &terminal.load_pressure_columns[p],
            terminal.load_pressure_inlet_factors[p],
        );
        forms.push(sym_product(&pressure, &flux));
        forms.push(sym_product(&load_pressure, &flux));
        linear.push(pressure);
        linear.push(load_pressure);
        fluxes.push(flux);
    }
    let mut material_heat = vec![0.0; total * total];
    for i in 0..nd {
        for k in 0..nd {
            material_heat[i * total + k] = original.owners[0].loss[i * nd + k];
        }
    }
    forms.push(material_heat);
    let mut load_heat = vec![0.0; total * total];
    for i in 0..2 {
        for k in 0..2 {
            load_heat[(nd + i) * total + nd + k] += 400.0 * model.residues[0][i * 2 + k].re;
            for a in 0..total {
                for b in 0..total {
                    load_heat[a * total + b] += fluxes[i][a] * model.d[i * 2 + k] * fluxes[k][b];
                }
            }
        }
    }
    forms.push(load_heat);
    let mut scales = vec![
        model.residues[0]
            .iter()
            .map(|x| x.re.abs())
            .fold(0.0, f64::max)
            .sqrt();
        2
    ];
    scales.push(1.0);
    PreparedQuadraticFlow::new_extended_controlled(
        &system,
        &split,
        &generator,
        &scales,
        dt,
        &forms,
        &linear,
        DenseExponentialBudget {
            max_dimension: 128,
            max_dense_bytes: 64 << 20,
            max_multiply_terms: 4_000_000_000,
        },
        || false,
    )
    .unwrap()
}
fn stored_load(model: &MimoModel, z: &[f64]) -> f64 {
    let mut sum = 0.0;
    for i in 0..2 {
        for j in 0..2 {
            sum += 0.5 * z[i] * model.residues[0][i * 2 + j].re * z[j];
        }
    }
    sum
}

#[test]
fn mixed_rank_original_dg_and_wall_states_keep_pressure_and_separate_integrated_work() {
    let gate = CancelGate::new_clock_free();
    for kind in [0, 1, 2] {
        let original = network(true);
        let nd = original.dimension;
        assert!(nd > 12, "wall coordinates must be live");
        let model = source(kind);
        let boundary = load(&model);
        let prepared = DgContinuousTerminal::new(
            &boundary,
            &[PortId(1), PortId(2)],
            &original,
            8 << 20,
            &gate,
        )
        .unwrap();
        assert_eq!(prepared.generator.len(), (nd + 2) * (nd + 2));
        assert_eq!(
            prepared.constraints.len(),
            match kind {
                0 => 2,
                1 => 1,
                _ => 0,
            }
        );
        let dt = 2e-6;
        let stepper = flow(&original, &prepared, &model, dt);
        let n = nd + 2;
        let mut state = vec![0.0; n + 1];
        state[n] = 2e-6;
        let mut supplied = 0.0;
        let mut duct_work = 0.0;
        let mut load_work = 0.0;
        let mut material_heat = 0.0;
        let mut terminal_heat = 0.0;
        // Run past the bore transit so both stored load motion and its heat
        // are observable in the original full-state trajectory.
        for _ in 0..200 {
            let step = stepper.evaluate(&state, &[0.0]).unwrap();
            supplied += step.integrals[0];
            duct_work += step.integrals[1] + step.integrals[3];
            load_work += step.integrals[2] + step.integrals[4];
            material_heat += step.integrals[5];
            terminal_heat += step.integrals[6];
            for p in 0..2 {
                let pressure = dot(&prepared.pressure_columns[p], &step.state[..n]);
                let effort = dot(&prepared.load_pressure_columns[p], &step.state[..n])
                    + prepared.load_pressure_inlet_factors[p] * step.state[n];
                let scale = prepared.pressure_columns[p]
                    .iter()
                    .zip(&step.state[..n])
                    .map(|(a, b)| (a * b).abs())
                    .sum::<f64>()
                    + prepared.load_pressure_columns[p]
                        .iter()
                        .zip(&step.state[..n])
                        .map(|(a, b)| (a * b).abs())
                        .sum::<f64>()
                    + 1e-20;
                assert!(
                    (pressure - effort).abs() < 2e-10 * scale,
                    "kind={kind},p={p},pressure={pressure},effort={effort}"
                );
            }
            for row in &prepared.constraints {
                let scale = row
                    .iter()
                    .zip(&step.state[..n])
                    .map(|(a, b)| (a * b).abs())
                    .sum::<f64>()
                    + 1e-20;
                assert!(dot(row, &step.state[..n]).abs() < 2e-10 * scale);
            }
            state = step.state;
        }
        let hd = 0.5 * dot(&state[..nd], &state[..nd]);
        let hl = stored_load(&model, &state[nd..n]);
        let scale = supplied.abs()
            + duct_work.abs()
            + load_work.abs()
            + hd
            + hl
            + material_heat
            + terminal_heat;
        assert!(material_heat > 0.0 && terminal_heat > 0.0 && hl > 0.0);
        assert!(
            (hd + material_heat + duct_work - supplied).abs() < 2e-10 * scale,
            "kind={kind}, H_duct={hd:e}, heat={material_heat:e}, work={duct_work:e}, input={supplied:e}, scale={scale:e}"
        );
        assert!(
            (hl + terminal_heat - load_work).abs() < 2e-10 * scale,
            "kind={kind}, H_load={hl:e}, heat={terminal_heat:e}, work={load_work:e}, scale={scale:e}"
        );
        assert!((duct_work - load_work).abs() < 2e-10 * scale, "kind={kind}");
    }
}

#[test]
fn mixed_rank_full_state_evolution_has_the_same_fixed_input_semigroup_across_clocks() {
    let gate = CancelGate::new_clock_free();
    let original = network(false);
    let model = source(1);
    let boundary = load(&model);
    let prepared = DgContinuousTerminal::new(
        &boundary,
        &[PortId(1), PortId(2)],
        &original,
        8 << 20,
        &gate,
    )
    .unwrap();
    let mut reference: Option<Vec<f64>> = None;
    for divisions in [1, 2, 4] {
        let stepper = flow(&original, &prepared, &model, 8e-6 / divisions as f64);
        let mut state = vec![0.0; prepared.inlet_column.len() + 1];
        *state.last_mut().unwrap() = 1e-6;
        for _ in 0..divisions {
            state = stepper.evaluate(&state, &[0.0]).unwrap().state;
        }
        if let Some(expected) = &reference {
            let scale = expected.iter().map(|x| x * x).sum::<f64>().sqrt();
            let error = state
                .iter()
                .zip(expected)
                .map(|(a, b)| (a - b).powi(2))
                .sum::<f64>()
                .sqrt();
            assert!(error < 2e-11 * scale, "error={error},scale={scale}");
        } else {
            reference = Some(state);
        }
    }
}

#[test]
fn physical_terminal_permutation_preserves_original_joint_operator_and_flows() {
    let gate = CancelGate::new_clock_free();
    let mut original = network(false);
    let mut model = source(1);
    let first = DgContinuousTerminal::new(
        &load(&model),
        &[PortId(1), PortId(2)],
        &original,
        8 << 20,
        &gate,
    )
    .unwrap();
    original.terminal_pressure_columns.swap(0, 1);
    original.terminal_input_columns.swap(0, 1);
    model.d = vec![model.d[3], model.d[2], model.d[1], model.d[0]];
    let r = &model.residues[0];
    model.residues = vec![vec![r[3], r[2], r[1], r[0]]];
    let second = DgContinuousTerminal::new(
        &load(&model),
        &[PortId(2), PortId(1)],
        &original,
        8 << 20,
        &gate,
    )
    .unwrap();
    let nd = original.dimension;
    let n = nd + 2;
    let at = |i| if i < nd { i } else { nd + 1 - (i - nd) };
    for i in 0..n {
        for j in 0..n {
            let a = first.generator[i * n + j];
            let b = second.generator[at(i) * n + at(j)];
            assert!((a - b).abs() < 3e-12 * (a.abs() + b.abs() + 1.0));
        }
        for p in 0..2 {
            let a = first.flow_columns[p][i];
            let b = second.flow_columns[1 - p][at(i)];
            assert!((a - b).abs() < 3e-12 * (a.abs() + b.abs() + 1e-20));
        }
    }
}

#[test]
fn mixed_rank_preparation_admits_complete_resources_and_never_mutates_original_history() {
    let gate = CancelGate::new_clock_free();
    let original = network(false);
    let boundary = load(&source(1));
    let before = boundary.pole_states().to_vec();
    let n = original.dimension + 2;
    assert!(
        DgContinuousTerminal::new(
            &boundary,
            &[PortId(1), PortId(2)],
            &original,
            256 * n * n - 1,
            &gate
        )
        .is_err()
    );
    assert!(DgContinuousTerminal::new(&boundary, &[PortId(1)], &original, 8 << 20, &gate).is_err());
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(
        DgContinuousTerminal::new(
            &boundary,
            &[PortId(1), PortId(2)],
            &original,
            8 << 20,
            &cancelled
        )
        .is_err()
    );
    assert_eq!(boundary.pole_states(), before);
    assert_eq!(
        DgContinuousTerminal::new(
            &boundary,
            &[PortId(1), PortId(2)],
            &original,
            8 << 20,
            &gate
        )
        .unwrap()
        .ports,
        [PortId(1), PortId(2)]
    );
}
