# Accelerator pilot doctrine

Status: ratified schema v1 for
`frankensim-extreal-program-f85xj.15.1`.

FrankenSim does not ship an accelerator backend. No GPU runtime has been
admitted as a production dependency, no device evidence record exists, and no
speedup or energy benefit is claimed. This doctrine exists so profiling can
refuse speculative accelerator work before silicon-specific code enters the
product path.

The policy is deliberately conditional. The end-to-end profiling campaign
owned by `frankensim-extreal-program-f85xj.15.2` must retain every workflow
phase, identify the top three kernels, and evaluate the named falsifier. Only
then may `frankensim-extreal-program-f85xj.15.3` consider one feature-gated
`[M]` kernel, and only after separate dependency admission and moonshot
displacement decisions. A refused-with-evidence pilot is a successful policy
outcome.

## Canonical policy, evidence class, and candidates

The section between the code-derived markers is rendered from
`fs_govern::ACCELERATOR_DOCTRINE`, `BACKEND_EVIDENCE_FIELDS`, and
`ACCELERATOR_CANDIDATES`, then checked byte-for-byte by the crate's G0 tests.
Edit the Rust registry and this section together.

<!-- BEGIN CODE-DERIVED ACCELERATOR DOCTRINE -->
| Policy field | Canonical value |
| --- | --- |
| Ambition | `[M]` |
| Profiling gate | `frankensim-extreal-program-f85xj.15.2` |
| Conditional pilot | `frankensim-extreal-program-f85xj.15.3` |
| Dependency ruling | `frankensim-extreal-program-f85xj.11.1` |
| Moonshot displacement | `frankensim-extreal-program-f85xj.16.3` |
| Thresholds | top-three wall: 50.0%; top-three energy where measured: 50.0%; selected kernel wall: 15.0% |
| Cancellation | kernel-batch boundary followed by request, drain, and finalize evidence |
| Permanent CPU reference | required |
| Named falsifier | refuse when the top three kernels jointly account for less than 50% of end-to-end wall time, jointly account for less than 50% of measured energy where energy is available, or contain no transfer-and-synchronization-suitable candidate with at least 15% of end-to-end wall time |
| Refusal | close the conditional pilot as refused-with-evidence; retain the CPU path and do not admit an accelerator dependency |
| No claim | governance schema only; no accelerator backend, dependency admission, device execution, speedup, energy saving, numerical equivalence, cancellation completion, or production authority is established |

| ID | Required field | Record mapping | Status | Requirement and boundary |
| --- | --- | --- | --- | --- |
| `AE-01` | device identity | `AcceleratorEnvironmentReceipt` in `fs-roofline`<br>no implementation locator | `explicitly-new` | Requirement: vendor, architecture, model, stable device identifier, memory topology, and capability fingerprint<br>No claim: a device name does not establish which binary ran or whether results are equivalent |
| `AE-02` | driver and runtime version | `AcceleratorEnvironmentReceipt` in `fs-roofline`<br>no implementation locator | `explicitly-new` | Requirement: exact driver, userspace runtime, backend API, and enabled feature versions<br>No claim: version strings do not authenticate the loaded implementation |
| `AE-03` | compiler and backend identity | `AcceleratorEnvironmentReceipt` in `fs-roofline`<br>no implementation locator | `explicitly-new` | Requirement: compiler, code-generation backend, flags, target features, build identity, and dependency decision receipt<br>No claim: compiler metadata does not prove semantic preservation or production admission |
| `AE-04` | kernel source identity | `fs_blake3::ContentHash` in `fs-blake3`<br>crates/fs-blake3/src/lib.rs | `existing` | Requirement: content identity of canonical kernel source and generated or embedded device binary<br>No claim: a content hash identifies bytes but does not authenticate their author or execution |
| `AE-05` | kernel and workload identity | `fs_roofline::KernelSpec` in `fs-roofline`<br>crates/fs-roofline/src/lib.rs | `existing` | Requirement: versioned kernel, dimensions, units, arithmetic-intensity model, dataset, phase, and denominator<br>No claim: a static kernel specification is not a timed production run |
| `AE-06` | CPU machine identity and measured axes | `fs_roofline::MachineAxes` in `fs-roofline`<br>crates/fs-roofline/src/axes.rs | `existing` | Requirement: CPU topology fingerprint plus measured bandwidth and compute axes before interpreting a comparison<br>No claim: CPU axes do not describe accelerator throughput, transfer cost, or energy |
| `AE-07` | reduction and determinism policy | `AcceleratorReductionPolicyReceipt` in `fs-roofline`<br>no implementation locator | `explicitly-new` | Requirement: fixed-order reduction topology or a versioned tolerance policy, including accumulation type and tie breaking<br>No claim: declaring a policy does not prove that a device kernel followed it |
| `AE-08` | CPU-versus-device numerical comparison | `AcceleratorEquivalenceReceipt` in `fs-evidence`<br>no implementation locator | `explicitly-new` | Requirement: per-QoI CPU and device values with units, numerical uncertainty, acceptance envelope, and corpus identity<br>No claim: no exact equivalence receipt exists today, and one future accepted metric would not establish equivalence outside its corpus, QoI, or envelope |
| `AE-09` | permanent CPU reference | `fs_roofline::RecordedProductionRun and FreshProductionEvidence` in `fs-roofline`<br>crates/fs-roofline/src/production.rs | `existing` | Requirement: retain the CPU implementation and exact comparison run; revalidate freshness whenever it is cited positively<br>No claim: a recorded operation is not fresh positive evidence until live authority revalidation succeeds |
| `AE-10` | cancellation, drain, and finalize outcome | `fs_exec::DrainFinalizeReport` in `fs-exec`<br>crates/fs-exec/src/cx.rs | `existing` | Requirement: executor-minted proof that every admitted kernel batch observed request, drained, and finalized<br>No claim: a caller-authored cancelled flag or dropped device handle is not drain evidence |
| `AE-11` | phase wall-time and energy attribution | `PipelineAttributionReceipt` in `fs-roofline`<br>no implementation locator | `explicitly-new` | Requirement: end-to-end phase totals, kernel shares, transfer and synchronization costs, profiling overhead, energy provenance, and named gaps<br>No claim: a kernel microbenchmark cannot stand in for workflow-level speed or energy benefit |
| `AE-12` | go or no-go decision | `AcceleratorPilotDecisionReceipt` in `fs-govern`<br>no implementation locator | `explicitly-new` | Requirement: bind the exact profile, thresholds, top-three ranking, suitability findings, dependency ruling, moonshot displacement, and terminal decision<br>No claim: a governance decision does not prove scientific correctness or future production fitness |

| ID | Candidate family | CPU source | Profiling hypothesis | Known falsifier pressure |
| --- | --- | --- | --- | --- |
| `AK-01` | D3Q19 LBM collide and stream | `fs-lbm`<br>crates/fs-lbm/src/d3q19/sparse.rs | regular active-tile arithmetic may expose enough parallel work per launch | halo traffic, sparse activity, boundary handling, and publication barriers may dominate |
| `AK-02` | sparse matrix-vector multiplication | `fs-sparse`<br>crates/fs-sparse/src/lib.rs | large repeated sparse operators may benefit from device memory bandwidth | irregular gathers and problem sizes may remain bandwidth-bound and transfer-dominated |
| `AK-03` | FFT and batched transforms | `fs-fft`<br>crates/fs-fft/src/lib.rs | regular batched transforms can offer substantial parallel work | pencil transposes, full-array passes, and host-device movement may erase kernel gains |
| `AK-04` | batched constitutive evaluation | `fs-material`<br>crates/fs-material/src/lib.rs | many independent material points may form efficient batches | no dedicated batch API exists today; branching laws, internal-state traffic, and small batches may underfill the device |
| `AK-05` | spectral path tracing | `fs-render`<br>crates/fs-render/src/tracer.rs | independent pixel samples provide abundant parallelism | divergent paths, deterministic accumulation, scene transfer, and feature-gated maturity may dominate |
<!-- END CODE-DERIVED ACCELERATOR DOCTRINE -->

## Evidence ownership

`existing` means the named Rust type and source locator exist today. It does not
mean that type already records a device run. `explicitly-new` means the future
pilot would need a new type at the named owner and that no implementation
locator is allowed in this schema. The mapping prevents adjacent CPU evidence
from being presented as if a GPU-specific receipt already existed.

The CPU reference is permanent. A pilot may add a device path but may not
replace the reference implementation, relax its comparison corpus, or cite a
bare recorded operation as fresh positive evidence. Cancellation likewise
requires executor-minted request, drain, and finalize evidence at kernel-batch
boundaries; a dropped future or device handle is not a cancellation receipt.

## Workflow-level decision boundary

The profiling denominator is the complete representative workflow, including
phases accelerators do not help. Kernel-only timing cannot satisfy the gate.
Energy is evaluated where the platform exposes credible data and its provenance
must be retained; an unavailable energy source is a named gap, not a fabricated
zero. Transfer, synchronization, compilation, setup, reporting, and ledger I/O
remain in the accounting.

Passing the numeric screen does not open a pilot by itself. The candidate must
also be one of the profiled top three, be suitable after transfer and
synchronization costs, receive an explicit production-dependency ruling, and
displace an existing `[M]` slot under the fixed-size moonshot policy. Even then,
the pilot is one feature-gated kernel and reports workflow-level Amdahl impact,
not a generalized FrankenSim GPU claim.
