# CONTRACT: fs-couple

> **Layout (restored 2026-11-09).** This contract is an index plus part files. Once this file grew
> past 1 MiB (1,048,576 bytes), some commits replaced the whole contract with their one new
> section; no edit below that size ever did. The cause has not been identified: the cutting commits
> are ordinary multi-file commits, so a tool that cannot read past 1 MiB and then writes the file
> whole is the likely kind of cause.
> This file was cut to a single section 18 times between 2026-10-30 and 2026-11-09; the text written before 2026-11-05 19:45 (2.3 MB) was absent from then until this restore.
>
> - All earlier contract text is in `contracts/part-NN.md`, unchanged. A higher part number is newer;
>   inside a part, sections are newest first, as they were in the single file.
> - Add a new section under "New sections" below, newest first.
> - There is no size limit on this file or on the parts (owner ruling, 2026-11-09). If a tool shows
>   you this file empty or cut short, do not write it back whole: add your section with a patch, and
>   check that `git diff --stat` shows no large deletion here before you commit.
> - The ten required sections are in the parts named under their headings below.

## New sections (newest first)

### An assembled body's modal solve kept between renders (`MaterialShellAssembly::reduce_with_modal_damping`'s `store`)

`reduce_with_modal_damping` and `ReciprocalStringBank::from_material_shell_assembly_prepared` take an optional store, the `SkinLoadStore` that skin loads use;
`reduce` and `from_material_shell_assembly` pass none. With a store, the tied pencil's solve is keyed by blake3 over its own inputs: both matrices (their dimensions,
each row's columns and the bits of its values), the window, and the `SliceOptions` by their `Debug` text, which names every field and writes each finite `f64` so
that it reads back to the same bits. A key the store holds is read back; any other is solved and kept. Each mode of the window (its place, λ, residual, interval and
every coordinate's bits) is kept under a key made from the solve's key and its place, and then the report's own numbers under the solve's key, so that this head
exists only once every mode it counts was kept. A keep encodes, and a load decodes, one mode's entry at a time beside the report. Every entry carries the solve's key
and ends in a hash of every preceding byte. An entry that is damaged, short, missing or another solve's, or whose length is not the pencil's order, refuses (`kept
shell modes are damaged, incomplete or another solve's`), and nothing is solved behind it. The rest of the reduction reads the same modes, so a kept solve gives
every bit of the image a solve gives. Only the bank's constructor reduces its body, so no rebuild solves again. `MaterialShellAssemblyImage::kept_modes` gives the
key and whether the solve was read. The key cannot see a change to the eigensolver, so a store holds the solves of one build; the example names each file with the
executable's hash.

The example's `--skin-load-store DIR` now applies after the card is parsed. It keeps `.shell-modes` entries for a v4 body, and skin loads only for a performance
with `body-radiation skin`; with neither it refuses. A render, a body response and a loaded response print `shell_modes_store_entry`, `shell_modes_loaded` and
`skin_load_store_keeps` (`shell-modes` alone without a skin load).

Evidence. `tests/material_string_bank/shell/assembly.rs` reduces a one-part body without a store, into an empty store and from it: `kept_modes` is `None`, then
`(key, false)`, then `(key, true)`; the modes, force ports, damping and intervals are the same bits all three times, the read adds no entry, another window has
another key, and a flipped head byte, a mode cut to 8 bytes, a removed mode and another solve's head under the key each refuse. Two control cards, each rendered
by one binary beside its own store-free render at the same time, give the same wav bytes and the same log but for the store's line: viola D♭4 (261 modes of
57,876 coordinates, 120,881,039 bytes) keeps its solve in 436.39 s against 432.48 s and reads it in 295.05 s against 400.91 s; double bass A1 at 3 workers (1,104
modes of 55,098 coordinates, 486,776,978 bytes) keeps it in 1230.13 s against 1219.92 s and reads it in 483.72 s against 1119.75 s, at 2.43 GiB of peak memory
against 3.36 GiB.

Not tested: a store that fails to write, and a render stopped while it keeps a solve, whose short mode file is never replaced, so a later read of it refuses.

### The prepared step reads its rows from copies by blocks of eight (`upper_rows::UpperBlocks`, `upper_rows::RowBlocks`)

The loss sweep by rows and the product on one thread advance eight consecutive rows together, each with its own sum. In a row-major square the eight entries of one
column lie a row apart, so each column took eight scalar loads. Each kernel now reads a copy in which they lie together. `UpperBlocks` holds a Gramian's upper triangle by
blocks of eight rows: the block's own square of eight rows and columns, row by row, with +0.0 below its diagonal and past the side; then, for a whole block, each column
right of it as eight entries. `RowBlocks` holds the transition's computed rows the same way over every column, and the rows after the last whole block as they are. A row
starts from the same value and adds the same terms in the same order as before, so every sum keeps its bits: the copies change where an entry is read, not what is added.

`sum_by_rows` makes the Gramians' copies once it has held them to their transposes. The first product on one thread makes the transition's (a `OnceLock`); a product that
`apply` splits over threads reads the row-major rows as before, and an operator whose products all run there never makes it. The row-major arrays stay: `apply` reads the
transition, the chain's sweep reads the Gramians, and the energy identity's check reads both when the operator is prepared; once an operator sums by rows, its sweep only
asks whether each Gramian exists. A Gramian's copy is about half its bytes: 64 entries a block, and eight a column right of each whole block, which is `s² / 2 + 4 s`
entries at a side `s` of whole blocks; the transition's is `product_rows` rows of the side. The admission (`admit_with_memory`) is unchanged: it counts `losses + 128`
squares of the side, and the copies add at most about two and a half (half of each of three Gramians and a transition of every row). A scratch build printed each copy's
bytes: the viola Db4 control card (338 modes, side 1014) holds two Gramian copies of 4,146,424 bytes and the transition's of 5,483,712; the violin D4 card with longitudinal
waves and the violin A4 control card (586 modes, side 1758) 12,420,448 twice and 16,483,008. The guitar F5 control card (379 modes, side 1137, 8 workers), the cello D3
control card (795 modes, side 2385, 7 workers) and the double-bass E1 control card (641 modes, side 1923, 7 workers) run their products on threads and hold no transition
copy: 5,209,080, 22,832,088 and 14,855,560 bytes twice. That is about 1.67 squares of the side on the first three cards and about 1.0 on the last three.

No path writes the arrays a copy holds once the copy exists. Their fields are private to `exact/prepared.rs`, where the one method of `PreparedNetwork` that takes
`&mut self` is `sum_by_rows`, which writes `upper_blocks`, `sums_by_rows` and the energy memo. Every other change builds a new operator, which makes its own copies and
starts with no kept sums, and replaces the old one whole. A changed link law, changed relaxation times or more material arms, a radiation load and a constant-loss binding
go through `preparation::prepare` (`coefficients.rs:342`, `relaxation.rs:223`, `radiation.rs:118`, `coupled.rs:591`), which calls `sum_by_rows` on the new operator when
the system sums by rows (`preparation.rs:124-129`). A prepared budget is bound through `prepare` too (`coupled.rs:540-545`), and `with_time_propagation` binds it after a
`StructureFirst` build (`coupled.rs:567`); the bank then asks for sums by rows (`bowed_string/runtime/reciprocal/bank.rs:396` and `:417`). The friction owner's network for
added components is a new system (`coupled/contact/friction/nonlinear.rs:505-506`) that takes the order of the one it replaces (`:508-509`).

Evidence. `tests/upper_row_forms.rs` gives the kernels the copies and holds them to the plain loops over the row-major square, bit for bit, at sides 1, 7, 8, 9, 16, 17 and
100 (whole blocks, a partial one, and 67 rows of side 100); a square whose lower triangle is moved has a copy with the same form. `tests/coupled_coefficient_update.rs`: a
system that sweeps and multiplies under one law and then takes another has, for the next three samples, the frames and the state bits of a system constructed with the
second law. Rendered beside 2819b9ab89's binary at the same time, each with the same frames and the same log: the viola Db4 control card takes 450.54 s against 533.46 s,
its wav byte-equal to the reference; the violin D4 card with longitudinal waves steps in 4.956 ms a frame against 6.296 ms over 600 s (the 67,712 frames both wrote), and
the violin A4 control card in 5.722 ms against 7.364 ms (52,864 frames).

Not tested: a material Gramian's copy against a reference in a system test; an operator whose declared workers change between its products, so that its copy is made
after products on threads.

### A loss sweep continues from the sums kept for its state (`upper_rows::upper_row_form_from`)

The sweeps of one sample share their state and differ in the held force, as its products do: a card with longitudinal waves computes about three a sample, one for each
nonlinear trial. A row of a loss form sums its columns right of the diagonal from left to right, and `z = [D q, v, memory | u]` has the state's `2 n + arms` entries
first. For each Gramian of the sweep by rows (component, connection, material) the operator keeps, for the last state it swept, every row's sum and magnitude over the
state's columns, with the bits of the state's entries (`upper_rows::StateForms`). A sweep whose state has the kept bits in every entry continues each row over the held
force's `n` columns; any other sweep computes the state's sums first and keeps them. With no arms a continued form adds 2.5 n² terms where the whole form adds 4.5 n².

The kernel is one pass over a range of columns (`rows_extended`): it adds to each row its terms for the columns of the range right of its diagonal, in the form's order
(the columns inside the row's block of eight, then the columns a whole block shares). Two passes over adjoining ranges add what one pass over their union adds, in its
order. `upper_row_form` is the pass over all the columns, and the continued form is the pass over the state's columns, kept, and then the pass over the force's. So a row
adds the same terms in the same order either way, the rows are added in row order as before, and each form keeps its bits.

The sums belong to one operator and live in it, as a product's kept sums do: every preparation builds a new operator, which starts with none, and a changed law that
keeps the scale `D` keeps the key. The slots are behind one mutex, held for the three forms of a sweep; a poisoned lock sums every row whole. The chain's sweep, which a
system that does not sum by rows computes, keeps nothing and is as it was.

Evidence. `tests/upper_row_forms.rs` holds the continued form to the plain loop of the form's order, bit for bit, at sides 1, 7, 8, 9, 16, 17 and 100, with a state of
0, 7, 8, 9, two thirds of the side, all but one and all of the entries: under three forces on one state, then on states that differ in the first entry or in the last,
by its last bit and by a half (one bit of one term can vanish in a row's sum, and a key that ignored the entry would then go unseen). The whole form is held to the same
loop by the test it had. `tests/coupled_coefficient_update.rs`: a refused step sweeps at a state, the law changes to one with the same own stiffnesses and another
coupling, and the next three steps have the frames of a system that never swept under the first law.

Not tested: a material Gramian's kept sums in a system test (the system test has no memory arm, and the kernel's test sums no particular Gramian).

### A bank built without its operator and prepared once (`ModalPreparation::StructureFirst`, `ReciprocalStringBank::with_time_propagation`)

A prepared bank computes a whole operator at its construction and again at every rebuild (the series exponential, a dense product for each doubling, the same again for
the finer check, and two products a level for each loss Gramian, all of side `3 n + arms`), and the friction owner takes its held-force line from that operator at once. A
card with residual flexibility builds two operators and a card that adds axial modes three, and only the last is ever stepped.

`ModalPreparation::StructureFirst(budget)` builds the complete bank as `HarmonicOnly` does, without a time propagator, and admits it as `Prepared` does. Every
construction and rebuild runs what a prepared construction runs before it computes its operator, in its order: the held-load pressure preparation, its reserve, and the
operator's mode, byte and work admission (`PreparedModalBudget::admit_structure`, which `PreparedNetwork::new` runs before it allocates). The harmonic admission follows.
Without an air bank its figures are the smaller ones (about 12 n² entries against at least 1152 n², 9 n³ terms against at least 27 n³ times 272); its one stricter limit
is 1024 modes. Air poles count toward that limit and enter its workspace squared, 8 (n + poles)² entries, where a prepared air bank's grow linearly: beside more than
about 11 n poles its bytes can refuse first. `material_strings` binds no air bank before the call. A system keeps its choice, so its rebuilds are admitted the same way.

`with_time_propagation` is the step from there to playback: a cold harmonic-only owner gains its propagator from the budget it holds.
- `CoupledModalSystem::with_time_propagation` hands the budget to the cold budget binding (`with_prepared_linear_propagation`), which prepares one operator and its
  pressure preparation from the system's present structure by `preparation::prepare`, the function that prepares a system constructed prepared.
- `FrictionModalSystem::with_time_propagation` then computes the held-force line, which a harmonic-only owner leaves at zero.
- `ReciprocalStringBank::with_time_propagation` then opts in to sums by rows and records the bank as prepared.

Each refuses a started owner and one that is not harmonic-only, so a second call is refused.

`material_strings` builds a playback's bank this way when the card's `prepared-linear` budget allows at most `MAX_HARMONIC_UNKNOWNS` (1024) modes, the capacity of a
harmonic-only construction. A larger budget is constructed prepared, as before. The call comes after the builders that only rebuild structure (residual flexibility, bridge
interfaces, body connections, cylinder gas loading, geometric stretching, longitudinal waves, material memory) and before computed radiation and the multi-bow conversion,
which treat a prepared bank and a harmonic-only one differently (`bank/radiation/rayleigh.rs`, `bank/radiation/bem.rs`, `contact/friction/multiple.rs`). The builders
before it read components and their states, memory, the port, pressure weights, the interface and the frequency bound, which is structure. None reads a response, a
compliance or the held-force line.

What changes for a refused card. A budget that cannot hold the bank's operator is refused where it was and with the words it had: by the constructor or by the rebuild
that first passes it, for its mode count, its bytes, its work or the pressure preparation's share. What moves is a refusal of the operator's own computation (a memory
arm that is not a passive Maxwell arm, a transition or a loss integral that fails its check, the energy identity, a held-force line that is not positive, a storage that
cannot be summed by rows): the call reports it, with the same words, where the constructor or a rebuild did. A builder that runs in between and refuses is heard first.

Evidence. `tests/material_string_nonlinear/bank/longitudinal.rs`: for five builders (none, cylinder gas loading, residual flexibility, geometric stretching, longitudinal
waves) a bank constructed structure-first, taken through the builder and then given the call sums by rows, records the choice, computes the same product rows and has,
over 48 bowed samples, the frames and the states of the bank constructed prepared and taken through the same builder; a second call is refused.
`tests/coupled_coefficient_update.rs`: two oscillators whose pressure has a quadrature part, joined by links with rest extensions. Built structure-first with a budget of
one byte, of too few bytes for the operator, of one term of work or of one mode, the system is refused with the words of the prepared construction. Built with room and
given the call, it has the 64 frames of the system constructed prepared, the pressure about its held equilibrium among them; a system constructed prepared, or already
given the call, is refused. The viola Db4 control card renders the same wav and the same log, every line. That card with a work budget of 1,000,000 prints the refusal
it printed (`Network(Budget { .. })`).

Not tested by a test of its own: the bridge-interface, body-connection, memory and cavity builders on a structure-first bank. The example's own tests reach the call with
some of them: its test files for torsion, bending and axial memory and for body relaxation set a `prepared-linear` budget of 64 modes. A budget above 1024 modes is
constructed prepared and runs no line of this section.

### A prepared product continues from the sums kept for its state (`upper_rows::product_from`)

The products of one sample share their state and differ in the held force: `z = [D q, v, memory | u]`, the state's `2 n + arms` entries first. A row's sum runs over the
columns in order, so its terms over the state's columns are the same for each of those products. The operator keeps them: for the last state it multiplied by, each
computed row's sum over the state's columns, with the bits of the state's entries (`upper_rows::StateSums`). A product whose state has the kept bits in every entry starts
each row from its kept sum and adds the held force's `n` columns. Any other product computes the state's sums first and keeps them. Either way a row adds the terms of one
left-to-right sum in that sum's order, so each entry keeps its bits, and a row that leaves the finite range is refused as before.

The sums belong to one operator and live in it. Every preparation (a changed link law, changed relaxation times, a radiation load, more material arms) builds a new
operator, which starts with none. The key alone would not tell two operators apart: it is the scaled state, and a changed law that keeps each mode's own stiffness keeps
the scale `D`.

The slot is behind a mutex, as the energy memo is, and a product holds it from the comparison to the last row. Two threads that multiplied by one operator would wait for
one another, and each would get the plain product's bits whichever state's sums it found. The kept sums are used on one thread: a product that `apply` splits over
threads (a side from 1024 with more than one declared worker) is computed by `apply` as before and keeps nothing.

Evidence. `tests/upper_row_forms.rs` holds `product_from` to one plain left-to-right sum a row, bit for bit, at sides 1, 7, 8, 9, 16, 17 and 100: on all the rows of a
square and on its first two thirds; with no state, a third of the columns, all but the last and all of them; under three forces on one state; then on a state that
differs in the last bit of its first entry and on one that differs in the last bit of its last. `tests/coupled_coefficient_update.rs`: a refused step multiplies at a
state, the law changes to one with the same own stiffnesses and another coupling, and the next three steps have the bits of a system that never multiplied under the
first law. The viola Db4 control card, and the violin D4 card with longitudinal waves over the first 900 s of its render, write the same frames as before.

Not tested: a product on more than one thread beside kept sums (no test target has a side of 1024).

### The prepared operator's product stops before the held force's rows (`CoupledModalSystem::prepared_product_rows`)

A prepared operator multiplies its transition by `z = [D q, v, memory | u]` (`exact/prepared.rs`). The transition's side is `3 n + arms`: the state's `2 n + arms` rows,
then one row for each entry of the held force `u`. The generator's rows of the held force are zero, so those rows of the transition are the identity's and their product is
`u` again. No caller reads it: `propagate` reads the first `2 n + arms` entries of the product, `response` each mode's displacement entry, and `velocity_response` each
mode's displacement and velocity entries. When an operator is prepared, `upper_rows::identity_rows` holds the last `n` rows to the identity's bit for bit: 1.0 on the
row's own column and +0.0 on every other. If they are, the product computes the state's rows only (two thirds of the terms when there are no arms). If one bit differs, it
computes every row as before. Each entry of the product is its own row's sum in column order, on any number of threads, so the entries that are read keep their bits.
A row that is left out could refuse only a `z` that is not finite. `state_vector` refuses such an `x` or `u` before any product, and every row that is kept would refuse
such a `z` too (a zero entry times an infinity is not a number).

Evidence. `tests/upper_row_forms.rs` holds `identity_rows` at sides 1, 7, 8, 9, 16, 17 and 100, with the held rows starting at the first row, at the middle and at the
last: true on the identity's rows under mixed rows; false with the one moved by a unit in its last place, with a zero made the smallest positive number, with the negative
zero, and at a wrong side; a row before the first held row is not read. `prepared_product_rows()` reports the count on the system and on the two friction owners; only
tests call it. `tests/coupled_coefficient_update.rs`: a prepared system of two modes computes 4 of its 6 rows before and after a changed link law, and 5 of 7 with one
material arm; without the prepared operator it is `None`. `tests/measured_string_bank.rs`: a prepared bank computes two rows a mode and one a material arm.
The viola Db4 control card renders the same wav and the same log as before.

Not tested: that the product does stop (no value shows it; the render times do); the fall back to every row through a system (no public call builds a transition whose
held rows are not the identity's; a build with one such entry planted computed 6 of 6 rows and passed every other test of its target).

### The prepared operator's storage by rows, and the dashpot norms it holds (`CoupledModalSystem::with_sums_by_rows`)

`with_upper_row_losses` and `sums_losses_by_upper_rows` are now `with_sums_by_rows` and `sums_by_rows`: one opt-in sets two orders. The first is the loss sweep's, described
in the section below. The second is the order of the operator's storage (`PreparedNetwork::hamiltonian`). By default that is one chain over each mode's kinetic term and
its `n` terms `0.5 q_i K_ij q_j`. By rows, `s_i` is the sum of `K_ij q_j` over every column in column order, eight rows advance together, and the rows are added in row
order as `0.5 v_i v_i + 0.5 q_i s_i` (`upper_rows::row_storage`); the material arms follow as before. The whole square is read. `mechanical_stiffness` adds a link as
`k * c_i * c_j`, which is `(k c_i) c_j` at (i, j) and `(k c_j) c_i` at (j, i), so the stiffness need not equal its transpose bit for bit and no half of it is relied on.

The storage differs from the chain's in its last bits. It is only compared, in two places that refuse. `advance` holds `(after - before) + dissipated - supplied` to 1e-11
of the largest of its terms. `exact::stage` hands the start state's storage to the bound on each link's force, which is held to `maximum_abs_connection_force_n`. No state,
path, sub-step count or scale reads it. `material_strings` prints `energy_order=rows` after `loss_order=upper-rows`.

The dashpot norm is not an opt-in and changes no bit. That bound takes the norm of the column a link's dashpot acts along: the damping column of a constant-loss link, the
link's own column otherwise. A step computed it for every link at every sample. The prepared operator computes it once for every link when it is prepared
(`prepared::column_norm`, the step's expression) and the bound reads it there; an owner without the prepared operator computes it in the step as before.
`connection_dashpot_column` and `connection_dashpot_norm` return the column and the held norm (`None` without the prepared operator). Only the test calls them.

Evidence. `tests/upper_row_forms.rs` holds `row_storage` to a plain loop of its order bit for bit and to the chain within 1e-13, at sides 1, 7, 8, 9, 16, 17 and 100, on
stiffnesses built as `mechanical_stiffness` builds them; from side 7 each differs from its transpose. `tests/coupled_coefficient_update.rs`: the two-orders test of the
section below runs both orders; on a fixture with one constant-loss link the held norm of every link has the bits of `det::sqrt` of its column's squares before and after
a changed law, and the constant-loss link's norm moves with the law. The viola Db4 and viola E3 control cards render the same wav, and the same log lines but for the added
`energy_order=rows`. A scratch build that ran both orders on those renders found the largest residual at 0.35 % and 0.37 % of its tolerance by rows (1.4 % with the chain's
storage) and the largest force bound at 0.90 % and 1.40 % of its limit in either order.

Not tested: which order the storage ran in through a system. No public value reads the storage; the oracle holds the kernel, and the opt-in is the flag whose loss half
the two-orders test observes. In the scratch build 437,162 of viola Db4's 438,751 evaluations of the storage, and 437,087 of viola E3's 438,751, had other bits in the two
orders. Not covered: the two transition products of a step.

### The prepared operator's loss sweep by rows over the upper triangle (`CoupledModalSystem::with_sums_by_rows`)

A step of a prepared operator evaluates three loss quadratics `z^T W z`: component, summed connection and material. By default each is one chain of additions over the whole
square in row-major order. `CoupledModalSystem::with_sums_by_rows` changes the order for that system. Every loss Gramian is exactly symmetric, so row `i` contributes
`z_i (W_ii z_i + 2 s_i)`, where `s_i` is the sum of `W_ij z_j` over the columns `j > i` in column order; eight rows advance together and the rows are added in row order
(`exact/prepared/upper_rows.rs`). That is half the terms and no single chain. The order is a function of the side alone, takes no thread count and runs on one thread. The
three losses differ from the chain's in their last bits. In a step of the linear operator, and under the friction owners, whose balances only refuse, the states, the
observer pressure and the stored energies do not read them and keep their bits. The trials of a system with nonlinear potentials are the exception: the paragraph on them
is below. The per-link port ledger and the air poles' forms keep the chain.

The opt-in is refused after the first sample, without the prepared operator, and when a loss Gramian differs from its transpose in any bit: the order never reads the lower
triangle. The system holds the order. `preparation::prepare` applies it to every operator prepared for the system afterwards (a changed link law, changed relaxation times, a
radiation load, more material arms) and holds each new operator's Gramians to their transposes again. `sums_by_rows()` reports it.

The bowed-string bank opts in whenever its operator is prepared: in `assemble_sources`, the one function of the bank's code that constructs a network, and in
`ReciprocalStringBank::with_prepared_linear_propagation`. Two builders add components through the friction owner instead, and it constructs the new network itself
(`FrictionModalSystem::with_internal_components_and_connections`): `with_longitudinal_waves_selected`, which the three other `with_longitudinal_waves` builders and the
`longitudinal-waves` row of `material_strings` reach, and `with_tensioned_fiber_ribbon`, which its `bow_fibers` row reaches. That network takes the order of the one it
replaces: by rows when that one summed by rows, the chain when it did not. `material_strings` prints `loss_order=upper-rows` for a bank that sums by rows. No other owner
opts in, and their losses keep the chain's bits.

A system with nonlinear potentials bound to it (a bank's longitudinal waves bind them through `with_mapped_nonlinear_potentials`) solves each sample by trials, and there
one more comparison reads the losses. `NonlinearFrictionState::iterate` accepts a trial when its force residual is within tolerance and its energy balance is within
`energy_absolute_tolerance_j + energy_relative_tolerance * scale`; the balance holds the staged frame's losses. A trial over that tolerance is iterated again, and a step
whose trials run out with the force residual within its tolerance is refused as an energy imbalance. So the order of the losses could change a state only at a sample whose
balance stands within the two orders' difference of that tolerance.

Evidence. `tests/upper_row_forms.rs` compiles the kernel's file a second time (the file names nothing of its crate). At sides 1, 7, 8, 9, 16, 17 and 100, on mixed-sign
Gramians, it holds the kernel to a plain loop of the same order bit for bit and to the chain within 1e-13, and holds `bitwise_symmetric` to one unit in the last place of one
lower entry, to 0.0 against -0.0 and to a wrong side. `tests/coupled_coefficient_update.rs`: two prepared systems that differ only in the order hold the same states,
pressure and stored energy bit for bit over 400 samples, through a changed link law at sample 200; their connection losses agree within 1e-13 and have other bits before and
after the change; the refusals. `tests/measured_string_bank.rs` and `tests/material_string_bank/shell/coating.rs` assert the order on a bank prepared after its construction
and on one constructed prepared. `tests/material_string_nonlinear/bank/longitudinal.rs` and `tests/material_string_relaxation/fiber_ribbon.rs` assert it on a prepared bank
before and after its axial modes, and a fiber ribbon's cells, are added; `tests/material_string_relaxation.rs` asserts that a bank without the prepared operator has not
chosen it after its axial modes are added. The viola Db4 and viola E3 control cards render the same wav as before, and the same log lines but for the added `loss_order=upper-rows`.
The violin D4 card with longitudinal waves, stopped 900 s into its render, has written the same frames as it did with the chain, and its log has both order lines.

Not tested: the refusal of an asymmetric Gramian through a system (no public call builds one); the component and material channels of a system (the fixture has neither
loss); sides from 1024, where the default order may run each Gramian on its own thread. Not covered: the two transition products of a step, which are most of what is
left of a string render's step. The prepared Hamiltonian is the section above.

### One evaluation of a coupled network's total energy for each accepted state (`CoupledModalSystem::total_energy_j`)

A step asks `CoupledModalSystem::total_energy_j` for the energy of its accepted state more than once: `exact::stage` asks, the candidate's gates ask, and a friction or contact
owner asks before it stages. The time-domain operator (`exact::ExactNetwork`) now keeps the last evaluation with the words it read from the state: each mode's displacement and
velocity, then the material memory, then the radiation state. A request whose words have the same bits is served the kept value. That value is the one the evaluation returns
for those words, so no reported number changes. One state is kept. A refusal is not kept. A system with no time-domain operator evaluates every request.
`CoupledModalSystem::total_energy_evaluations()` counts the evaluations since the operator was prepared, served requests excluded. It exists to observe the memo; nothing
physical and no receipt reads it.

The constants the evaluation reads are not in the key: the modes' frequencies, the link columns and laws, the relaxation branches and the radiation model. They are the inputs
of `preparation::prepare`, and the memo is a field of the operator prepared from them. Every place that replaces one of them installs, in the same function, the operator
prepared from the new ones (grep of the non-test files under `render/schedule/force/coupled` at 372d94a611): `new_with_columns_inner`, `with_prepared_linear_propagation`,
`with_constant_loss_connections`, `with_radiation_load`, `with_additional_material_relaxation`, `install_connection_coefficients`, and `exchange` of
`coefficients/rates.rs`, which swaps the models, the branches and the operator together. The link columns are written by the constructor only, and fs-modal's
`ModalAcousticTimeModel` has no method that changes a mode. A new place that changed one of these without installing its operator would be served the energy of the constants
before.

Evidence. `tests/modal_multi_tangential.rs`: over 60 steps of the two-contact owner, which asks once before the network's stage asks twice, the network evaluates its energy
once for each accepted state. `tests/modal_material_composition.rs`: an arm with a zero projection relaxes while the one mode stays at rest, so its memory is the only state
word that moves; the energy asked after each of 8 steps is the one before times `exp(-2 dt / tau)`. Existing tests guard the rest: `tests/coupled_coefficient_update.rs` asserts
the energy change across a changed link law on an unchanged state, and `tests/coupled_impulse.rs` moves velocities alone. A memo that never serves makes the first test
count 180 evaluations for its 60 steps, three a step. The viola Db4 control card (432000 samples at 144 kHz) renders the same wav and the same log lines as before.

Not tested: the radiation state moving alone (a cold attach starts it at zero, and it then moves only with a mode); the system with no time-domain operator. Not covered: the
other energy evaluations of a step (`PreparedNetwork::hamiltonian` has its own memo) and the loss sweep, about half of that card's time for a step.

### The linear impact description on the composition kernel (`render::plate::impact::graph`, feature `kernel`)

`ImpactGraph::new` takes the inputs `LinearImpactSystem::new` takes (bodies, contacts, volumes, config) and a sink temperature, declares them as a graph on fs-system
and seals it. A body that no volume touches is an fs-nlmodal `ModalStructure`. Each set of bodies that volumes join is one fs-system `ExactPortHamiltonianComponent`, with
storage `diag(omega^2) + sum k c c^T` over its interleaved `(q, p)` (`k = K A^2 / V`, `c` the signed modal areas over `A`, the module's volume spring) and each mode's rate on
its momentum. The one contact point is an fs-dcontact `ContactInterface` on the load line between the two parts its row touches; each part takes it through its port
`contact`, weighted by minus the row at end `a` and plus the row at end `b`, and every part publishes its ends. Each part's heat goes to a sink at the given temperature.
`step` reads its frame from the kernel's ledger: stored energy is the build's plus the ledger's stored change, dissipation the sinks' heat over the step, the residual the
step's change of the uncorrected first-law residual. `state` holds the interleaved coordinates in the description's order, so a caller rebuilds (a catch, a launch) from it
as it does from the module's.

Refused by name: a nonlinear body; more than one contact point (distributed or several contacts); a contact row that does not join exactly two parts; a contact with
internal loss (`ContactInterface` has no heat port); a sink temperature that is not positive; a step the kernel shortens. Not declared: dampers, external generalized
forces, the module's energy ceiling. Feature `kernel` (`dep:fs-system`, `fs-system/phs`, `fs-dcontact/component`) is off by default; a build without it compiles none of this.

Evidence (`tests/impact_graph.rs`, required-features `kernel`): the second graph of `tests/drum_strike_graph.rs` (a stick strikes a two-mode head joined by its air to a
second two-mode head), here declared from the module's own inputs. First touch at sample 11, 4 strikes; every coordinate within -270.84 dB of its full scale of the module's
run (the hand-built graph's figure); stored energy within 8.3e-17 J of the module's (start 1.28e-2 J), summed residual -4.3e-18 J, heat equal to the module's to ten digits;
after a catch, both rebuilt from their own head states, within -270.01 dB; no allocation in a warm step; the refusals above. Used by engine-drum's `"mechanics":"kernel"`.
Not covered: radiation, several contacts, dampers, speed (engine-drum's bench).

### The nonlinear bore with its own gas (`bore_gas` beside `mouth_load`)

A card may now give the nonlinear bore (`bore nonlinear`) a gas of its own inside and a different room outside (MB 24783). `bore_gas`, the gas the bore, cup and lip hold,
stands beside `mouth_load modal` on `bore nonlinear` (it used to be refused there, as it needs the room named by the flanged load): the record must be `thermochemical` with its
packs, because the nonlinear bore evolves the `bore_gas` row's own NASA-9 laws (the gas source, the wall gas, the duct fill state); the `ambient` row is then the room, which the
kernel of `mouth_load` names as its exterior (the two-gas flanged law of fs-duct, with the BORE gas's characteristic impedance as the kernel's reference, which is the bore end's
own), and which the receiver and the observer keep; the `ambient` row may be any gas. Refused by name: `bore_gas` on `bore nonlinear` without `mouth_load` (a baffled-piston end would
radiate into the bore gas), a `bore_gas` that is not `thermochemical` on `bore nonlinear`, a static pressure other than the ambient's, and `modal_baffle` beside `bore nonlinear`
(the multimodal bank is the linear route's). Without the row nothing changes: equal gases give the same kernel bit for bit (fs-duct tests, digests read on the unchanged builder), and
a card without `bore_gas` renders the same bytes (the digest of a 2400-sample nonlinear render, piston end and `mouth_load` end, was read on the unchanged sources and on the changed ones
on one host and compiler; the values are in the commit message).

What it does not do, and what it costs: the receiver is still the baffled-piston integral of the plane mouth flow (not the finite-distance integral over every mouth mode of the linear
route); the load is the plane-mode equivalent of the flanged load, which sits sharp of the multimodal bore above about 400 Hz on the tenor trombone (fs-duct contract); above the evaluated
top the realized load is the bore gas's characteristic impedance where the radiation resistance tends to the exterior's: 1.3 % for the 27 C saturated breath state against dry 20 C air
(`rho c` 408.03 against 413.28 Pa s/m, a reflection of 0.64 % in amplitude), and the realized load's worst error against the two-gas load between its bins is 3.98e-4 `Z_inf` up to
3.6 kHz of a 4 kHz claim (fs-duct test, gases of that `rho c` ratio); and the walls carry their own 128-pole memory, not the linear bore's.

### Reed opening tap

`ReedBoreVoice::record_reed_opening` and `ReedPerformance::record_reed_opening` give a massive reed an opt-in tap: at the end of every solver sample one text line `y v` (the blade's opening in metres and velocity in m/s, as the closure stop and the lay see them) goes to a caller-supplied writer. It reads state after the step and changes no audio sample (test `reed_trace_tests`: the pressure bits with and without the tap are equal, and with the stop on `y` never goes below zero, with `closure_stop off` it can). The `music_render wind` command routes the existing `--aperture-flows FILE` to it on a `frankensim-reed-performance-v1` file (`--labium` stays flue-only); on any other file type the flag keeps its earlier meaning and refusals.

### The reed-noise row against the print's power (measured)

A held, passive reed (damping ratio 30, it cannot oscillate) on an 8 m cylinder (radius 6.935 mm; the characteristic-line ladder refuses this card at 9, 10, 15 and 25 m) with the blow held at 1800 Pa was run with `constriction_noise` on and off, and the outgoing wave p+ read in the window of samples 800 to 2200 (the first reflection returns at 2239), where the bore carries only the wave the slit launches. The print gives, per direction, W = K(St)^2 A dp^2 / (4 rho c) (Oldham and Ukpoho as printed in Guedel, AIVC CR02 eq. 3), an outgoing wave of rms K dp / 2 with dp the static drop across the slit (1368 Pa here). Measured per third octave from 630 to 8000 Hz (12 bands) the row's wave is below that by 6.0, 13.5, 8.2, 9.1, 7.4, 10.9, 7.9, 9.0, 11.4, 9.3, 8.9 and 9.7 dB, band mean 9.29 dB (6.3 dB if the printed power is the total over both sides of the duct, which is not settled). The likely reason is that the row is a series pressure source on a slit whose small-signal impedance is a few times the bore's, so the launched wave is only part of the source (not measured separately). The test that measured it (`noise_power_tests`, band mean within 1.5 dB) fails by that amount and is not part of this change; it returns when the print's definition (one- or two-sided power, a dipole force or a monopole volume source) is settled. No row, value or engine path changes.

### Reed residual flow as a card row

Almeida, Vergez and Causse 2007 (JASA 121(1) 536-546), Sec. III.1 and Fig. 9, measure a natural cane oboe reed on an
artificial mouth: the flow "first stabilises at a certain minimum value" instead of vanishing at the beating pressure
("very hard to completely close the reed"; "narrow channels remaining between the two blades ... behaving like rigid
capillary ducts"), and Fig. 10 shows a 1/2 power dependence of that residual flow. The engine's reed flow went to exactly zero
at zero opening. The reed card now accepts `residual_flow AREA_M2` right after `closure_stop` (or `contraction`,
`staple_recovery`, `slit_length`, the reed row), before `tongue_contact`, for a massive reed on the inviscid law (refused with
`slit_length` and with `staple_recovery`). The jet area becomes `C_c w max(h, 0) + AREA` (`C_c` the `contraction` value, else 1): the
voice hands the junction the opening `max(h, 0) + AREA / (C_c w)` through the existing flow-opening argument of
`step_with_junction`, so no kernel and no function signature changes; the blade's mechanics, lay contact and closure stop keep
the blade opening `h` alone. The row has no default and carries a measured quantity: one reed's Fig. 9 read at 450 dpi gives
0.094 to 0.100 mm^2 (mean 0.097e-6 m^2, 22 to 27 ml/s at 32 to 47 kPa), about 3 % of that reed's fitted rest area, and it is
one reed's. A card without the row renders as before, bit for bit (the flow-opening argument is unchanged when the voice has no
residual). `ReedBoreVoice::set_residual_flow_area` and `residual_flow_area_m2` are the voice-level switch. Tests (G0 to G2):
`render::schedule::reed::residual_flow_tests`.

`AREA` is the contracted JET area, with discharge coefficient 1: a shut reed passes `AREA sqrt(2 dp / rho)` and the junction sees the floor opening `AREA / (C_c w)`, 27.3 um for 9.33e-8 m^2 on a 5.6 mm reed with `C_c` 0.611 (`AREA / w` would be 16.7 um). A shut reed with no row takes the junction's exact closed-valve branch (`characteristic_junction_with_recovery`, `a == 0`); any positive `AREA`, however small, takes the general branch, which reproduces it to rounding (about 1e-16 relative). So a vanishing area agrees with no row to rounding and not bit for bit: `g1` pins it at 1e-20 (the pressure within 1e-9 of the peak), and on the Bb4 oboe card `residual_flow 1e-300` moves the first closure by 1.5e-21 m and leaves the clip level equal to five digits (3185.7 LSB, rms of 1 to 2 s). Measured on that card (rest jet area 1.222 mm^2, blow 7.4 kPa, rho 1.2 assumed): 9.33e-8 m^2 is 7.64 % of the rest jet area, 11 to 12 % of the mean flow, and cancels 30 % of the reed's quasi-static negative slope (reed -2.317e-9, leak +7.0e-10 m^3/s/Pa). It raises the oscillation threshold. Clip level (rms of 1 to 2 s) in dB against AREA 0: heavy-reed Bb4 at 7.4 kPa 0, -0.25, -0.76, -1.57, -2.45 for AREA 0, 1e-8, 3e-8, 6e-8, 9.33e-8; F4 heavy reed at 6.865 kPa 0, -0.20, -0.56 and silent at 9.33e-8; light reed (7.69e-6 kg, damping ratio 0.32) Bb4 at 7.4 kPa 0, -0.21, -0.66, -1.03 (4.5e-8) and silent at 6e-8, where the same card sounds at blows of 8.5, 9.5 and 10.5 kPa (-0.36, -0.14, -0.97 dB): the leak raises the blowing pressure the reed needs to start, and a card that holds a leak the print names says so.

### Reed closure stop as a card row

Commit a857c6f17 made the reed step's rigid closure stop at zero opening unconditional (a reed that the free implicit
midpoint solve would carry below zero opening ends the step at exactly zero with end velocity 0; no restitution is printed
for a reed on its lay). The reed card now accepts `closure_stop on|off` right after `contraction` (or `staple_recovery`,
`slit_length`, the reed row), before `tongue_contact`, for a massive reed. `on` is the default and a card without the row
renders as before, bit for bit (`step_with_contact_response` still passes the stop). `off` selects
`OpeningPropagator::MidpointFree`, the same implicit midpoint step without the stop, so a reed that goes below zero keeps that
opening and the lay contact alone closes it; a step whose free solution stays at or above zero is the same step either way.
The row carries no value: it is an A/B selection, not a physical claim. `ReedBoreVoice::set_closure_stop` and `closure_stop`
are the voice-level switch (massive reed, voice not advanced). The loaded finite-gas path keeps passing `false` as before.
Tests (G0, G1): `render::schedule::reed::closure_stop_tests`.

### A prescribed-flow step of the DG port load

`DgPortLoad::step(volume_flow_m3_s)` advances the load alone by one sample
under a prescribed bore volume flow: the crate-private trial of
`preview_step`, accepted at once. It adds no numerics and no state; a refusal
publishes nothing, as for the trial. It is the public counterpart of
`step_with_mouth_flow` for a load with no driven mouth, and it exists so that
a test outside the crate can march the load as its own reference:
`tests/dg_finite_gas/molecular.rs` called the crate-private trial and the
`dg_finite_gas` target did not compile from 2026-10-29 until 2026-11-10. No
claim: unlike `step_with_mouth_flow` it takes no cancellation gate, and it is
not a coupling path (nonlinear participants keep using the trial).

### Thick isotropic wall: Lame above thickness/radius 0.05

`LocalMaterialWall::base_patch` (graph, DG and acoustic-cell walls all pass through it) keeps the
membrane expressions at thickness/radius <= 0.05, bit for bit. Above that an ISOTROPIC wall takes
Lame's thick cylinder (Timoshenko and Goodier, thick-walled cylinder under internal pressure; free
outer face; plane stress for a stress-free axial restraint, plane strain for a strain-fixed one):
the stiffness per inner area is `p / u(a)` and the surface density is the kinetic energy of the same
radial field `u = A r + B / r` per inner area (a Rayleigh estimate, below the through-thickness
resonances). The Poisson ratio is the wall's own resolved material's; the Kelvin-Voigt resistance and
Maxwell arms scale with `stiffness / E` as before. An ORTHOTROPIC wall over the limit is still
refused. At the gate for nu = 1/3 the two laws differ by 4.06 % in stiffness and 0.68 % in mass
(stress-free), 4.82 % and 0.12 % (strain-fixed); the tests carry those amounts in their names.
This replaces the refusal `local membrane wall requires positive geometry and thickness/radius <= 0.05`
for isotropic walls (test `g3_the_same_input_over_a_twentieth_of_the_radius_now_takes_the_thick_wall_law`).

### Adaptive laminar enthalpy and external hydraulic continuation

`LaminarEnthalpy::advance_windkessel_adaptive` reuses the original coupled
BDF1 coarse/two-fine controller. Coarse and first-half solves borrow the same
accepted flow/compliance history; the second half borrows the actual first-half
pair. Rejected comparisons advance no external state. The existing ordinary,
Boussinesq and driven entries retain their original semantics and leave the
new optional hydraulic observations empty. No second time controller or
extrapolated physical field is introduced.

An explicit hydraulic norm joins the original caloric, velocity and work
criteria. It checks each port's compliance and added pressure, actual integrated
native Q, separate C/L storage, and signed plus resolved-magnitude resistor,
reservoir and port work. Signed and magnitude Q*dt*depth throughput have their
own cubic-metre tolerance, so equal endpoint flow cannot conceal different
transported volume. Pressure scaling removes the held reservoir offset and
uses the start compliance/Ohmic state; it invents no initial acceleration.
Storage changes and errors use difference of squares before dimensional
multiplication. All spans and local corrections come from the original
first-order controller, including exact target clipping and one-ULP cap tails.

Only accepted fine receipts contribute compensated per-port work and signed/
absolute transported volume totals. Initial Q is independently obtained from
the original accepted BDM edge transport; initial/final C*x²/2 and L*Q²/2
changes use the paired physical states and declared depth. Each port and their
aggregate must satisfy port work minus reservoir work, both resistor losses,
both storage changes and the separate native BDF defect within the explicit
window allowance. This external-circuit identity does not add heat, replace
the independent native caloric balance, or certify whole-system passivity.

All coarse/rejected work consumes the shared native-step and thermal-endpoint
ceilings. Source, material, momentum, stream and thermal refusals remain
terminal. Cancellation or a failed final circuit account rejects the complete
window, including any internally accepted prefix. The final step's retained
compliance coordinates continue with the returned original velocity and clock.
Circuit parameters/references and physical boundary data remain held; this
entry adds no scheduled circuit events or Boussinesq/load combination.

Native regression definitions cover symmetric compliance discharge with zero
fluid/thermal change, hydraulic-only refinement and an independent accepted BE
recurrence; signed equilibrium all-port plug transport and exact work/storage;
nonzero resistor losses excluded from heat; decimal target completion; and
late source, work-budget and cancellation rollback after an accepted prefix.
All seven changed Rust files passed the pinned parser, formatter and
`rustfmt --check`. Bounded independent source review found no paired-history
or accounting blocker. Full native FEM/thermal compilation and these native
integration test executions remain unverified in the current partial restored
environment; no runtime pass is inferred from those checks.

## Parts

| Part | Bytes | `##` sections | First (newest) section | Last (oldest) section |
|---|---|---|---|---|
| [`contracts/part-08.md`](contracts/part-08.md) | 449466 | 145 | Driven full-curve Hertz friction with original thermal sources | CST/DKT hereditary material memory on the original modal basis |
| [`contracts/part-07.md`](contracts/part-07.md) | 442691 | 141 | Thermal porous follower pressure through enthalpy feedback | Electrical contact heat in thermally live actuator supports |
| [`contracts/part-06.md`](contracts/part-06.md) | 446940 | 130 | Complete material and acoustic memory in exponential impact | Shared material cavitation, films and bubble acoustic history |
| [`contracts/part-05.md`](contracts/part-05.md) | 420894 | 27 | Finite solvent heat during gas absorption and exsolution | Conformance tests |
| [`contracts/part-04.md`](contracts/part-04.md) | 449410 | 65 | No-claim boundaries | Live gas changes through the shared acoustic port |
| [`contracts/part-03.md`](contracts/part-03.md) | 444523 | 95 | Finite thermochemical feedback through the complete acoustic load | Nonlinear spatial current with complete material heat exchange |
| [`contracts/part-02.md`](contracts/part-02.md) | 425077 | 113 | Source-resolved Hertz islands on shared finite thermal cells — v0 | Complete sliding-film layer heat on the impact clock |
| [`contracts/part-01.md`](contracts/part-01.md) | 285661 | 57 | Generalized-Newtonian film resistance in impact mechanics | Native joint friction with finite-cell thermal feedback |

## Purpose and layer

Full text: [`contracts/part-05.md`](contracts/part-05.md), under the same heading.

## Public types and semantics

Full text: [`contracts/part-05.md`](contracts/part-05.md), under the same heading.

## Invariants

Full text: [`contracts/part-05.md`](contracts/part-05.md), under the same heading.

## Error model

Full text: [`contracts/part-05.md`](contracts/part-05.md), under the same heading.

## Determinism class

Full text: [`contracts/part-05.md`](contracts/part-05.md), under the same heading.

## Cancellation behavior

Full text: [`contracts/part-05.md`](contracts/part-05.md), under the same heading.

## Unsafe boundary

Full text: [`contracts/part-05.md`](contracts/part-05.md), under the same heading.

## Feature flags

Full text: [`contracts/part-05.md`](contracts/part-05.md), under the same heading.

## Conformance tests

Full text: [`contracts/part-05.md`](contracts/part-05.md), under the same heading.

## No-claim boundaries

Full text: [`contracts/part-04.md`](contracts/part-04.md), under the same heading.
