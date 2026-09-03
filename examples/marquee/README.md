# Marquee 2D Topology Optimization Example (`bracket-2d.fsim`)

This directory contains the canonical tracked example for FrankenSim's **Marquee 2D Topology Optimization** pipeline (Journey B, bead `frankensim-rc-root-q61wp.20`).

The study defines a 2D elastic plate subject to fixed clamping on the left boundary and traction on the right boundary. The optimizer adjusts hole radii using projected gradient descent with Armijo backtracking line search while strictly preserving the material volume budget.

---

## Section-by-Section Anatomy

Every `.fsim` study document adheres strictly to the **Five Explicits** and the frozen schema specification (`STUDY_FSIM_VERSION = 1`):

### 1. Root & Version (`:version 1`)
Defines the schema version envelope. Admitted documents must match `STUDY_FSIM_VERSION` (currently `1`).

### 2. `metadata`
Carries user-facing context, intended decision, decision gate (`scoping-estimate`, `design-selection`, `compliance-signoff`), and consequence class (`advisory`, `reliability`, `safety-critical`).

### 3. `versions`
Explicit engine and schema versions preventing silent semantic drift across toolchain upgrades.

### 4. `seeds`
Counter-based pseudo-random generator seed (`:rng 1337`). Ensures bit-level reproducibility of all stochastic or initial operations on the same ISA.

### 5. `budgets`
Declares execution budgets:
- `:wall-time 60 s`: Maximum execution wall clock time.
- `:memory 1073741824 B`: Memory consumption ceiling.
- `:max-iterations 8`: Maximum optimization iterations permitted.

### 6. `capabilities`
Names the required execution capabilities:
- `"optimization.marquee-topopt"`: 2D topology optimization driver.
- `"geometry.sdf"`: Signed distance field representations.
- `"physics.cutfem"`: Unfitted CutFEM elasticity solver.

### 7. `units`
Declares the units doctrine (`:system "SI"`). Omission of declared units triggers an immediate admission refusal (`project-undeclared-units`).

### 8. `domain`
Specifies the geometric domain:
- `:type sdf-plate-with-holes`: Plate domain with circular voids.
- `:bounds ((0.0 0.0) (1.0 1.0))`: Bounding box of the unit square $[0, 1]^2$.
- `:initial-holes`: Centers and initial radii of voids in the plate.

### 9. `physics`
Discretization parameters for CutFEM elasticity:
- `:type elasticity-2d`: 2D plane stress/strain elasticity.
- `:mesh-level 4`: Quadtree refinement level for background grid cells.

### 10. `scenario`
Supports and external loads:
- `:fixed-boundary left`: Clamped Dirichlet boundary condition along $x = 0$.
- `:load-region right`: Traction load along $x = 1$. Non-boundary load declarations trigger refusal (`study-load-non-boundary`).

### 11. `objective`
Optimization goal:
- `:type compliance`: Structural compliance functional $J(u) = \int_{\Gamma_N} f \cdot u \, ds$.
- `:sense minimize`: Minimization objective.
- `:unit "J"`: Energy units $[M \cdot L^2 \cdot T^{-2}]$. Supplying incompatible units (e.g. `"W"` or `"m"`) triggers refusal (`study-objective-dimension-mismatch`).

### 12. `constraints`
Design constraints:
- `:volume-fraction 0.853`: Target solid volume fraction. Values outside $(0, 1)$ trigger refusal (`study-volume-fraction-out-of-bounds`).

### 13. `optimizer`
Hyperparameters for the optimization loop:
- `:type projected-gradient`: Gradient projection onto the volume equality constraint.
- `:step-size 1.0`: Initial gradient descent step size.
- `:r-min 0.08`, `:r-max 0.20`: Lower and upper bounds on hole radii.
- `:steps 8`: Number of optimization steps.

---

## Running the Study

Run the study end-to-end using the `frankensim study` CLI verb:

```bash
cargo run -p fs-cli --bin frankensim -- study examples/marquee/bracket-2d.fsim ledger.db
```

### Running with a Budget Limit

To run with an iteration budget (e.g., stopping after 2 iterations):

```bash
cargo run -p fs-cli --bin frankensim -- study examples/marquee/bracket-2d.fsim ledger.db --budget 2
```

When an iteration budget is exhausted mid-loop:
1. The process exits with exit code `6` (`exit::BUDGET`).
2. The last certified iterate is durable and retained.
3. The run receipt and report mark status as `"budget-exhausted"` with complete lineage.
