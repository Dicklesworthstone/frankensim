//! Feature-gated conformance battery for proposal-only generators.
#![cfg(feature = "proposal-gen")]

use std::collections::BTreeSet;

use fs_gen::{GraphGenerator, ModelCard, MutationKernel, Proposal, ShapePrior, acquire};

fn card() -> ModelCard {
    ModelCard {
        model: "test-model".to_string(),
        corpus_hash: "test-corpus".to_string(),
        determinism: "deterministic".to_string(),
    }
}

fn corpus() -> Vec<Vec<f64>> {
    vec![
        vec![0.0, 0.0, 1.0],
        vec![1.0, 0.5, 1.5],
        vec![2.0, 1.0, 2.0],
        vec![3.0, 1.5, 2.5],
    ]
}

#[test]
fn proposal_payload_exits_only_through_validator() {
    let proposal = Proposal::new(vec![1.0, 2.0, 3.0], card());
    let rejected = proposal
        .clone()
        .promote(|_| Err("missing certificate"))
        .unwrap_err();
    assert_eq!(rejected.reason, "missing certificate");

    let promoted = proposal
        .promote(|x| (x.len() == 3).then_some(()).ok_or("wrong dimension"))
        .unwrap();
    assert_eq!(promoted, vec![1.0, 2.0, 3.0]);
}

#[test]
fn shape_prior_is_deterministic_and_dimension_checked() {
    let prior = ShapePrior::fit(&corpus());
    let a = prior.propose(11).promote(|x| {
        (x.len() == 3 && x.iter().all(|v| v.is_finite()))
            .then_some(())
            .ok_or("bad proposal")
    });
    let b = prior.propose(11).promote(|x| {
        (x.len() == 3 && x.iter().all(|v| v.is_finite()))
            .then_some(())
            .ok_or("bad proposal")
    });
    assert_eq!(a.unwrap(), b.unwrap());
    assert!(prior.density(&[1.0, 0.5, 1.5]).is_finite());
    assert!(prior.density(&[1.0, 0.5, 1.5]) > 0.0);
}

#[test]
#[should_panic(expected = "shape prior needs non-empty design vectors")]
fn shape_prior_rejects_zero_dimensional_corpus() {
    let _ = ShapePrior::fit(&[Vec::new()]);
}

#[test]
#[should_panic(expected = "mutation kernel corpus rows must share one finite")]
fn mutation_kernel_rejects_ragged_corpus() {
    let _ = MutationKernel::fit(&[vec![0.0, 1.0], vec![2.0]]);
}

#[test]
fn mutation_kernel_is_deterministic() {
    let kernel = MutationKernel::fit(&corpus());
    let a = kernel.propose(17).promote(|x| {
        (x.len() == 3 && x.iter().all(|v| v.is_finite()))
            .then_some(())
            .ok_or("bad mutation")
    });
    let b = kernel.propose(17).promote(|x| {
        (x.len() == 3 && x.iter().all(|v| v.is_finite()))
            .then_some(())
            .ok_or("bad mutation")
    });
    assert_eq!(a.unwrap(), b.unwrap());
}

#[test]
fn mutation_kernel_handles_degenerate_covariance() {
    let kernel = MutationKernel::fit(&[vec![1.0, 2.0, 3.0]]);
    let direction = kernel.propose(19).promote(|x| {
        (x.len() == 3 && x.iter().all(|v| v.is_finite()))
            .then_some(())
            .ok_or("bad degenerate mutation")
    });
    assert_eq!(direction.unwrap().len(), 3);
}

#[test]
fn graph_generator_proposes_simple_edges() {
    let generator = GraphGenerator::fit(4, &[vec![(0, 1), (1, 2)], vec![(2, 3)]]);
    let edges = generator.propose(23, 3).promote(|edges| {
        let unique: BTreeSet<(usize, usize)> = edges.iter().copied().collect();
        (edges.len() == 3
            && unique.len() == edges.len()
            && edges.iter().all(|&(a, b)| a < b && a < 4 && b < 4))
        .then_some(())
        .ok_or("not a simple graph")
    });
    assert_eq!(edges.unwrap().len(), 3);
}

#[test]
fn graph_generator_fills_requested_possible_edge_count() {
    let generator = GraphGenerator::fit(4, &[]);
    let edges = generator.propose(0, 6).promote(|edges| {
        let unique: BTreeSet<(usize, usize)> = edges.iter().copied().collect();
        (edges.len() == 6 && unique.len() == 6)
            .then_some(())
            .ok_or("missing possible edge")
    });
    assert_eq!(edges.unwrap().len(), 6);
}

#[test]
#[should_panic(expected = "graph generator needs at least one node")]
fn graph_generator_rejects_zero_nodes() {
    let _ = GraphGenerator::fit(0, &[]);
}

#[test]
#[should_panic(expected = "requested edge count exceeds")]
fn graph_generator_rejects_impossible_edge_count() {
    let generator = GraphGenerator::fit(2, &[]);
    let _ = generator.propose(29, 2);
}

#[test]
fn acquisition_is_deterministic_and_bounded_by_candidate_count() {
    let prior = ShapePrior::fit(&corpus());
    let archive = vec![vec![0.0, 0.0, 1.0], vec![3.0, 1.5, 2.5]];
    let a = acquire(&prior, &archive, 5, 8, 31);
    let b = acquire(&prior, &archive, 5, 8, 31);
    assert_eq!(a.len(), 5);
    assert_eq!(a.len(), b.len());
    for (p, q) in a.iter().zip(&b) {
        assert_eq!(p.inspect_for_logging_only(), q.inspect_for_logging_only());
        assert_eq!(p.card, q.card);
    }
}

#[test]
#[should_panic(expected = "archive rows must match the prior dimension")]
fn acquisition_rejects_malformed_archive_rows() {
    let prior = ShapePrior::fit(&corpus());
    let _ = acquire(&prior, &[vec![0.0, 1.0]], 1, 1, 37);
}
