//! G0: an owner that declares its stored entropy closes it. The entropy its
//! ports carry in plus the entropy its receipt says it generated equals the
//! measured change, to roundoff relative to the inventories it reports; an
//! owner that declares nothing is not audited.

#![allow(clippy::float_cmp)]

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cancelled, Cx, ExecMode, StreamKey};
use fs_ports::StableId;
use fs_ports::component::conformance::{
    Case, Failure, Observed, UnitExponents, check_component, check_unit_rescaling, unit_factor,
};
use fs_ports::component::{
    Class, Component, ContractViolation, EntropyStepAudit, PortLayout, PortRole, PortSpec,
    Proposal, Receipt, audit_entropy_step, audit_step,
};
use fs_qty::Dims;

fn heat(name: &str, role: PortRole) -> PortSpec {
    PortSpec::heat(StableId::new(name).expect("id"), role).expect("heat port")
}

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let key = StreamKey {
            seed: 1,
            kernel_id: 2,
            tile: 3,
            iteration: 4,
        };
        let cx = Cx::new(gate, arena, key, Budget::INFINITE, ExecMode::Deterministic);
        f(&cx)
    })
}

/// One heat port taking 2000 W at 400 K for 0.5 s: 2.5 J/K in, exactly.
fn one_port(role: PortRole) -> ([PortSpec; 1], [f64; 1], [f64; 1]) {
    let (input, output) = if role == PortRole::FlowIn {
        ([2000.0], [400.0])
    } else {
        ([400.0], [2000.0])
    };
    ([heat("heat", role)], input, output)
}

fn audited(
    ports: &[PortSpec],
    inputs: &[f64],
    outputs: &[f64],
    generated: f64,
    old: Option<f64>,
    new: Option<f64>,
) -> Result<Option<EntropyStepAudit>, ContractViolation> {
    let layout = PortLayout::new(ports).expect("layout");
    audit_entropy_step(ports, &layout, 0.5, inputs, outputs, generated, old, new)
}

#[test]
fn port_entropy_and_the_booked_generation_close_the_measured_change() {
    for role in [PortRole::FlowIn, PortRole::EffortIn] {
        let (ports, input, output) = one_port(role);
        let audit = audited(&ports, &input, &output, 1.5, Some(0.0), Some(4.0))
            .expect("entropy closes")
            .expect("measured inventory");
        assert_eq!(audit.port_entropy_in_j_per_k, 2.5);
        assert_eq!(audit.entropy_generated_j_per_k, 1.5);
        assert_eq!(audit.stored_change_j_per_k, 4.0);
        assert_eq!(audit.balance_defect_j_per_k, 0.0);
        // Heat leaving takes its entropy out; a negative inventory is legal,
        // since the reference is the owner's.
        let leaving = [-2000.0];
        let (input, output) = if role == PortRole::FlowIn {
            (leaving, output)
        } else {
            (input, leaving)
        };
        let audit = audited(&ports, &input, &output, 0.5, Some(-2.0), Some(-4.0))
            .expect("entropy closes")
            .expect("measured inventory");
        assert_eq!(audit.port_entropy_in_j_per_k, -2.5);
        assert_eq!(audit.balance_defect_j_per_k, 0.0);
    }
    // A stateless conductor passing 1000 W from 400 K to 200 K for 0.5 s
    // takes 1.25 J/K in and puts 2.5 J/K out: it must book the 1.25 J/K
    // difference as generated, and a claim of a reversible step is refused.
    let ports = [
        heat("hot", PortRole::EffortIn),
        heat("cold", PortRole::EffortIn),
    ];
    let (input, output) = ([400.0, 200.0], [1000.0, -1000.0]);
    let audit = audited(&ports, &input, &output, 1.25, Some(0.0), Some(0.0))
        .expect("the conductor closes")
        .expect("zero inventory is available");
    assert_eq!(audit.port_entropy_in_j_per_k, -1.25);
    assert_eq!(audit.balance_defect_j_per_k, 0.0);
    // The tolerance counts the generation with the absolute transfers:
    // 1.25 + 2.5 + 1.25 J/K.
    assert_eq!(audit.tolerance_j_per_k, 128.0 * f64::EPSILON * 5.0);
    assert!(matches!(
        audited(&ports, &input, &output, 0.0, Some(0.0), Some(0.0)),
        Err(ContractViolation::EntropyBalance {
            port_entropy_in_j_per_k,
            entropy_generated_j_per_k,
            ..
        }) if port_entropy_in_j_per_k == -1.25 && entropy_generated_j_per_k == 0.0
    ));
}

#[test]
fn the_tolerance_is_strict_and_relative_to_the_reported_inventories() {
    let (ports, input, output) = one_port(PortRole::FlowIn);
    let closed = audited(&ports, &input, &output, 1.5, Some(0.0), Some(4.0))
        .expect("closes")
        .expect("inventory");
    // 128 EPSILON of the largest term: the 4 J/K inventory, which equals the
    // 2.5 J/K carried in plus the 1.5 J/K generated.
    let tolerance = closed.tolerance_j_per_k;
    assert_eq!(tolerance, 128.0 * f64::EPSILON * 4.0);
    let within = audited(
        &ports,
        &input,
        &output,
        1.5,
        Some(0.0),
        Some(4.0 + tolerance / 4.0),
    )
    .expect("a quarter of the tolerance is admitted")
    .expect("inventory");
    assert!(
        within.balance_defect_j_per_k != 0.0,
        "accepted defects remain visible"
    );
    for off in [2.0 * tolerance, -2.0 * tolerance, 1.0e-3] {
        assert!(matches!(
            audited(&ports, &input, &output, 1.5, Some(0.0), Some(4.0 + off)),
            Err(ContractViolation::EntropyBalance { .. })
        ));
    }
    // The same 1e-3 J/K defect against inventories whose reference sits
    // 1e12 J/K away is inside 128 EPSILON of them: the gate loosens with the
    // reference the owner states, as the method's documentation says.
    let offset = 1.0e12;
    let candidate = offset + 4.0 + 1.0e-3;
    let loose = audited(&ports, &input, &output, 1.5, Some(offset), Some(candidate))
        .expect("admitted at an offset reference")
        .expect("inventory");
    assert_eq!(loose.tolerance_j_per_k, 128.0 * f64::EPSILON * candidate);
    assert!(loose.balance_defect_j_per_k.abs() > 9.0e-4);
}

#[test]
fn unavailable_or_invalid_accounts_never_pass_as_a_balance() {
    let (ports, input, output) = one_port(PortRole::FlowIn);
    let layout = PortLayout::new(&ports).expect("layout");
    assert_eq!(
        audited(&ports, &input, &output, 0.5, None, None),
        Ok(None),
        "unavailable is not a passing balance"
    );
    let audit = |dt_s: f64, ports: &[PortSpec], inputs: &[f64], outputs: &[f64], generated| {
        audit_entropy_step(ports, &layout, dt_s, inputs, outputs, generated, None, None)
    };
    assert!(matches!(
        audit(0.0, &ports, &input, &output, 0.5),
        Err(ContractViolation::Step { .. })
    ));
    assert!(matches!(
        audit(0.5, &[], &input, &output, 0.5),
        Err(ContractViolation::VectorLength {
            vector: "ports",
            ..
        })
    ));
    assert!(matches!(
        audit(0.5, &ports, &[], &output, 0.5),
        Err(ContractViolation::VectorLength {
            vector: "inputs",
            ..
        })
    ));
    assert!(matches!(
        audit(0.5, &ports, &input, &[f64::NAN], 0.5),
        Err(ContractViolation::NonFinite {
            quantity: "entropy port values"
        })
    ));
    for bad in [f64::NAN, f64::INFINITY] {
        assert!(matches!(
            audit(0.5, &ports, &input, &output, bad),
            Err(ContractViolation::NonFinite {
                quantity: "entropy_generated_j_per_k"
            })
        ));
    }
    assert!(matches!(
        audit(0.5, &ports, &input, &output, -1.0e-300),
        Err(ContractViolation::NegativeProduction {
            quantity: "entropy_generated_j_per_k",
            ..
        })
    ));
    // Heat at 0 K has no entropy; heat over f64::MAX seconds has no finite one.
    for (dt_s, kelvin) in [(0.5, 0.0), (f64::MAX, 400.0)] {
        assert!(matches!(
            audit(dt_s, &ports, &input, &[kelvin], 0.5),
            Err(ContractViolation::NonFinite {
                quantity: "port entropy transfer"
            })
        ));
    }
    for (old, new) in [(Some(0.0), None), (None, Some(3.0))] {
        assert!(matches!(
            audited(&ports, &input, &output, 0.5, old, new),
            Err(ContractViolation::EntropyAvailability { old_available, new_available })
                if old_available == old.is_some() && new_available == new.is_some()
        ));
    }
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        for (old, new, expected) in [(bad, 3.0, "accepted"), (0.0, bad, "candidate")] {
            assert!(matches!(
                audited(&ports, &input, &output, 0.5, Some(old), Some(new)),
                Err(ContractViolation::StoredEntropy { endpoint, .. }) if endpoint == expected
            ));
        }
    }
    // Finite terms that balance, 1e308 J/K carried in and 1e308 J/K
    // generated, around an inventory change that is not finite.
    assert!(matches!(
        audited(
            &ports,
            &[1.0e308],
            &[0.5],
            1.0e308,
            Some(-1.0e308),
            Some(1.0e308)
        ),
        Err(ContractViolation::NonFinite {
            quantity: "entropy account"
        })
    ));
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum Fault {
    None,
    /// Its stored entropy reads a millionth too large.
    Skewed,
    /// Its receipt books no generation.
    Reversible,
    /// Its stored entropy is unavailable above 350 K.
    Lapsing,
}

/// A body of heat capacity `capacity_j_k` at one temperature [K], its state,
/// warmed through one heat port by a backward-Euler step: the heat rate `q`
/// it takes in raises it to `T1 = T0 + dt q / C`, the temperature the port
/// returns. Its stored entropy, `C ln(T / reference_k)`, is zero at
/// `reference_k`. Heat entering at the end temperature carries `dt q / T1`,
/// which `ln x >= 1 - 1/x` keeps at or below the change `C ln(T1 / T0)`; the
/// body books the rest, `C ln(T1 / T0) - dt q / T1`, as generated, a number
/// that does not depend on the reference.
struct Body {
    ports: [PortSpec; 1],
    capacity_j_k: f64,
    reference_k: f64,
    fault: Fault,
}

impl Body {
    fn new(fault: Fault) -> Self {
        Self {
            ports: [heat("heat", PortRole::FlowIn)],
            capacity_j_k: 10.0,
            reference_k: 300.0,
            fault,
        }
    }

    fn entropy(&self, kelvin: f64) -> f64 {
        self.capacity_j_k * (kelvin / self.reference_k).ln()
    }
}

impl Component for Body {
    type State = f64;
    type Error = Cancelled;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }
    fn class(&self) -> Class {
        Class::ExactLinear
    }
    fn stored_energy_j(&self, kelvin: &f64) -> Result<f64, Cancelled> {
        Ok(self.capacity_j_k * kelvin)
    }
    fn stored_entropy_j_per_k(&self, kelvin: &f64) -> Result<Option<f64>, Cancelled> {
        Ok(match self.fault {
            Fault::Skewed => Some(self.entropy(*kelvin) * (1.0 + 1.0e-6)),
            Fault::Lapsing if *kelvin > 350.0 => None,
            _ => Some(self.entropy(*kelvin)),
        })
    }

    fn propose(
        &self,
        old: &f64,
        dt_s: f64,
        inputs: &[f64],
        new: &mut f64,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, Cancelled> {
        cx.checkpoint()?;
        let heat_j = dt_s * inputs[0];
        *new = old + heat_j / self.capacity_j_k;
        outputs[0] = *new;
        let generated = if self.fault == Fault::Reversible {
            0.0
        } else {
            self.capacity_j_k * (*new / *old).ln() - heat_j / *new
        };
        Ok(Proposal::Step(Receipt {
            entropy_generated_j_per_k: generated,
            ..Receipt::lossless(heat_j, 0.0)
        }))
    }
}

/// The body with the default `stored_entropy_j_per_k`.
struct Undeclared(Body);

impl Component for Undeclared {
    type State = f64;
    type Error = Cancelled;

    fn ports(&self) -> &[PortSpec] {
        self.0.ports()
    }
    fn class(&self) -> Class {
        self.0.class()
    }
    fn stored_energy_j(&self, kelvin: &f64) -> Result<f64, Cancelled> {
        self.0.stored_energy_j(kelvin)
    }

    fn propose(
        &self,
        old: &f64,
        dt_s: f64,
        inputs: &[f64],
        new: &mut f64,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, Cancelled> {
        self.0.propose(old, dt_s, inputs, new, outputs, cx)
    }
}

/// 2000 W into the body at 300 K for 0.5 s: it ends at 400 K.
fn case() -> Case<'static, f64> {
    Case {
        state: &300.0,
        stale: &7.0,
        dt_s: 0.5,
        inputs: &[2000.0],
    }
}

#[test]
fn the_conformance_kit_audits_a_declared_entropy_and_refuses_each_fault() {
    let gate = CancelGate::new_clock_free();
    let seen = with_cx(&gate, |cx| {
        check_component(&Body::new(Fault::None), &case(), cx)
    })
    .expect("the body closes its entropy");
    let Observed::Stepped { receipt, .. } = &seen else {
        panic!("unexpected shortening: {seen:?}")
    };
    // 10 ln(4/3) J/K stored against 2.5 J/K carried in at 400 K.
    assert!((receipt.entropy_generated_j_per_k - 0.376_820_724_517_809).abs() < 1.0e-12);
    for (fault, refused) in [
        (Fault::Skewed, "balance"),
        (Fault::Reversible, "balance"),
        (Fault::Lapsing, "availability"),
    ] {
        let body = Body::new(fault);
        let failure = with_cx(&gate, |cx| check_component(&body, &case(), cx))
            .expect_err("each fault is refused");
        assert!(
            matches!(
                (refused, &failure),
                (
                    "balance",
                    Failure::Contract(ContractViolation::EntropyBalance { .. })
                ) | (
                    "availability",
                    Failure::Contract(ContractViolation::EntropyAvailability {
                        old_available: true,
                        new_available: false,
                    }),
                )
            ),
            "{fault:?}: {failure}"
        );
    }
    // The first law cannot see a missing generation: the reversible body's
    // energy audit passes.
    let body = Body::new(Fault::Reversible);
    let (mut new, mut output) = (0.0, [0.0]);
    let Proposal::Step(receipt) = with_cx(&gate, |cx| {
        body.propose(&300.0, 0.5, &[2000.0], &mut new, &mut output, cx)
    })
    .expect("proposal") else {
        panic!("unexpected shortening")
    };
    let layout = PortLayout::new(body.ports()).expect("layout");
    audit_step(
        body.ports(),
        &layout,
        0.5,
        &[2000.0],
        &output,
        &receipt,
        body.stored_energy_j(&300.0).expect("energy"),
        body.stored_energy_j(&new).expect("energy"),
    )
    .expect("the energy audit admits the reversible claim");
    // Undeclared, the same reversible body is not audited: the default is
    // unavailable, not zero storage.
    let undeclared = Undeclared(Body::new(Fault::Reversible));
    assert!(matches!(
        undeclared.stored_entropy_j_per_k(&300.0),
        Ok(None)
    ));
    with_cx(&gate, |cx| check_component(&undeclared, &case(), cx))
        .expect("an undeclared entropy is not audited");
}

/// Each base dimension alone, then all six together, as base-2 exponents in
/// the order `[m, kg, s, K, A, mol]`.
const UNIT_CHANGES: [UnitExponents; 7] = [
    [3, 0, 0, 0, 0, 0],
    [0, -2, 0, 0, 0, 0],
    [0, 0, 5, 0, 0, 0],
    [0, 0, 0, 2, 0, 0],
    [0, 0, 0, 0, -3, 0],
    [0, 0, 0, 0, 0, 1],
    [3, -2, 5, 2, -3, 1],
];

/// The body in the units `exponents` name: capacity [J/K] and reference [K].
fn body_in_units(exponents: &UnitExponents) -> Body {
    let mut body = Body::new(Fault::None);
    body.capacity_j_k *= unit_factor(Dims([2, 1, -2, -1, 0, 0]), exponents);
    body.reference_k *= unit_factor(Dims([0, 0, 0, 1, 0, 0]), exponents);
    body
}

/// The body's state in those units. The kit's `scale_state` takes the state
/// by reference.
#[allow(clippy::trivially_copy_pass_by_ref)]
fn kelvin_in_units(kelvin: &f64, exponents: &UnitExponents) -> f64 {
    kelvin * unit_factor(Dims([0, 0, 0, 1, 0, 0]), exponents)
}

#[test]
fn a_declared_entropy_rescales_with_its_reference_to_the_bit() {
    let gate = CancelGate::new_clock_free();
    for exponents in &UNIT_CHANGES {
        with_cx(&gate, |cx| {
            check_unit_rescaling(&body_in_units, &case(), &kelvin_in_units, exponents, cx)
        })
        .unwrap_or_else(|failure| panic!("units {exponents:?}: {failure}"));
    }
    // A reference left in the old units changes no step, only the entropy
    // the start reports.
    let kept = |exponents: &UnitExponents| {
        let mut body = body_in_units(exponents);
        body.reference_k = Body::new(Fault::None).reference_k;
        body
    };
    let failure = with_cx(&gate, |cx| {
        check_unit_rescaling(&kept, &case(), &kelvin_in_units, &UNIT_CHANGES[6], cx)
    })
    .expect_err("refused");
    let Failure::Law { law, detail } = failure else {
        panic!("expected a law failure, got {failure}");
    };
    assert_eq!(
        law,
        "the rescaled start's energy, mass, entropy and step limits are the base ones in the new units"
    );
    assert!(detail.contains("stored_entropy_j_per_k"), "{detail}");
}
