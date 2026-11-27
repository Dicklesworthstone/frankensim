//! Exactly zero factor rows leave the dense system and decay on their own.
//! Goldens are the f64 bits of the full dense evolution, taken before the
//! change (bead ks7e6).
use super::*;

/// Body momenta, memory amplitudes, [heat, memory change] and planned terms.
type Step = (Vec<f64>, Vec<f64>, Vec<f64>, u64);
fn step(
    states: &[RigidBodyState],
    links: &[RelaxationNetworkLink],
    memory: &[Vec<[f64; 6]>],
    dt: f64,
) -> Step {
    let (properties, mut limits) = (line_properties(), memory_limits());
    limits.maximum_branches = 16;
    let poses: Vec<_> = states.iter().map(|s| s.pose()).collect();
    let body = &properties[..states.len()];
    let mut source =
        PreparedRelaxationNetwork::new(&poses, body, links, dt, limits, || false).unwrap();
    let efforts = vec![[0.0; 6]; links.len()];
    let out = source
        .advance(states, memory, &efforts, controls(), || false)
        .unwrap();
    let momenta = out.states_after.iter().flat_map(|s| {
        let (p, h) = (s.linear_momentum_world(), s.angular_momentum_body());
        [p.x, p.y, p.z, h.x, h.y, h.z]
    });
    let memory = out.memory_after.iter().flatten().flatten();
    let energy = vec![out.dissipated_energy_j, out.memory_energy_change_j];
    let terms = source.planned_preparation_terms();
    (momenta.collect(), memory.copied().collect(), energy, terms)
}
/// Each quantity within 1e-12 of its largest golden value; all bits on failure.
fn golden(quantities: &[(&[f64], &[u64])]) {
    let mut close = true;
    for &(actual, bits) in quantities {
        let golden: Vec<f64> = bits.iter().map(|&b| f64::from_bits(b)).collect();
        let scale = golden.iter().fold(0.0, |s: f64, g| s.max(g.abs()));
        close &= actual.len() == bits.len();
        for (a, g) in actual.iter().zip(&golden) {
            close &= (a - g).abs() <= 1e-12 * scale;
        }
    }
    let found: Vec<Vec<u64>> = quantities
        .iter()
        .map(|(actual, _)| actual.iter().map(|x| x.to_bits()).collect())
        .collect();
    assert!(close, "{found:#x?}");
}

#[test]
fn six_stated_axes_keep_the_dense_path_bit_for_bit() {
    let i = UnitQuaternion::IDENTITY;
    let states = [
        state(X.scale(-1.0), i, v([0.4, -0.2, 0.3]), v([0.1, 0.2, -0.3])),
        state(Y.scale(0.2), i, v([-0.1, 0.5, 0.0]), v([0.0, -0.4, 0.2])),
    ];
    // Lower triangular with a nonzero diagonal: every row is stated.
    let rows = core::array::from_fn(|i| {
        core::array::from_fn(|j| if j > i { 0.0 } else { 1.0 + j as f64 })
    });
    let links = [RelaxationNetworkLink {
        viscous: ViscousNetworkLink {
            bodies: [0, 1],
            arms_body_m: [v([0.1, 0.2, 0.0]), v([-0.1, 0.0, 0.3])],
            damping: diagonal([0.5, 0.4, 0.3, 0.2, 0.1, 0.05]),
        },
        branches: vec![SpatialRelaxationBranch::from_gram_rows(rows, 0.4).unwrap()],
    }];
    let memory = [vec![[0.1, -0.2, 0.05, 0.0, 0.3, -0.1]]];
    let (momenta, memory, energy, _) = step(&states, &links, &memory, 0.2);
    let all = [momenta, memory, energy].concat();
    let found: Vec<u64> = all.iter().map(|x| x.to_bits()).collect();
    #[rustfmt::skip]
    let bits: [u64; 20] = [
        0x3fd76a0280c50ea3, 0xbfc85c42ab87e250, 0x3fd489c1f9502cda, 0x3fb950bd19776c8d,
        0x3fc77b2d4631bcf4, 0xbfdab65e1acf9f1c, 0xbfb0db3d36476dbe, 0x3fdf615488f7245b,
        0xbf9568ec61cf9a70, 0x3f75ae851d2f4c5b, 0xbfd9e0f23602a4ee, 0x3fd55e8e30696359,
        0x3f85c698dc00d4a5, 0xbfa1efdd9bdde510, 0x3fa8584a5b18e59a, 0xbf9d445639473ba4,
        0xbf8c277784de3309, 0xbfbbaa35580cf239, 0x3fad3db19c917ecf, 0xbfb16f7b534f7cf6,
    ];
    assert_eq!(found, bits, "{found:#x?}");
}

#[test]
fn a_thirteen_branch_x_only_mount_drops_its_dead_axes_from_the_plan() {
    let states = &line_states([-0.5, 0.5, 0.0])[..2];
    let branch = |i: usize| {
        let mut rows = [[0.0; 6]; 6];
        rows[0][0] = (10.0 * (2 + i) as f64).sqrt();
        let tau = 10f64.powi(i as i32 - 6);
        SpatialRelaxationBranch::from_gram_rows(rows, tau).unwrap()
    };
    let links = [RelaxationNetworkLink {
        viscous: axial_link([0, 1], 0.3),
        branches: (0..13).map(branch).collect(),
    }];
    let x_only = |i: usize| [0.01 * (i + 1) as f64, 0.0, 0.0, 0.0, 0.0, 0.0];
    let mut memory = vec![(0..13).map(x_only).collect::<Vec<_>>()];
    // Stored energy on two unstated axes: it exerts nothing and only decays.
    memory[0][6][3] = 0.2;
    memory[0][12][5] = -0.1;
    let dt = 0.125;
    let (momenta, after, energy, terms) = step(states, &links, &memory, dt);
    let x = |a: u64| [a, 0, 0, 0, 0, 0];
    #[rustfmt::skip]
    let mut memory_bits = [
        0xbec9f6d2915e8a6d, 0xbf03e12c91ee23ba, 0xbf3cc3bb5ce35a7d, 0xbf74949f17f2aea0,
        0xbfb031d7e16a422e, 0xbfd0540d050d1b43, 0xbfcbf8a572745e3f, 0xbfcb03c25a69fb38,
        0xbfcb851cc0abd64b, 0xbfcc1eca60f9986e, 0xbfcca838bc6dd2d1, 0xbfcd1f356e77d299,
        0xbfcd85620efb7efc,
    ]
    .map(x);
    memory_bits[6][3] = 0x3fc697881d8ace82;
    memory_bits[12][5] = 0xbfb9999963e9b49d;
    let x_momenta = [momenta[0], momenta[6]];
    golden(&[
        (&x_momenta, &[0x3fe625882c68a7da, 0xbfe625882c68a7da]),
        (&after, memory_bits.as_flattened()),
        (&energy, &[0x3fb9b374049c8fff, 0x3fc47eb7793b9aea]),
    ]);
    near(after[6 * 6 + 3], 0.2 * (-dt).exp(), 1e-15);
    near(after[12 * 6 + 5], -0.1 * (-dt / 1e6).exp(), 1e-15);
    // Twelve body and thirteen stated coordinates (all 90 before: 16327580316).
    assert_eq!(terms, 413_362_646);
}

#[test]
fn only_an_exactly_zero_row_is_dead() {
    let states = &line_states([-0.5, 0.5, 0.0])[..2];
    let terms = |root: f64| {
        let mut rows = [[0.0; 6]; 6];
        rows[0][0] = 1.0;
        rows[1][1] = root;
        let branch = SpatialRelaxationBranch::from_gram_rows(rows, 0.5).unwrap();
        let links = [RelaxationNetworkLink {
            viscous: axial_link([0, 1], 0.3),
            branches: vec![branch],
        }];
        step(states, &links, &[vec![[0.0; 6]]], 0.1).3
    };
    assert_eq!(terms(1e-300), terms(1.0));
    assert_eq!(terms(-0.0), terms(0.0));
    assert!(terms(0.0) < terms(1.0));
}
