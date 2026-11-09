# Chronomètre Souverain · Calibre FS-1903 Tourbillon (`apps/pocket-watch`)

An interactive real-time 3D + WebAssembly multiphysics showcase for **FrankenSim**, simulating a Patek-style Grand Complication mechanical pocket watch constructed from sourced material cards in `data/matdb/seed-v1`.

## Quick Start

```bash
cd apps/pocket-watch
npm run wasm      # Rebuilds the Rust -> WASM engine (engine-watch)
npm test          # Runs Node WASM integration tests across all multiphysics regimes
npm run build     # Typechecks (tsc --noEmit) and builds the production Vite bundle
npm run dev       # Launches the interactive 3D + 48 kHz WebAudio laboratory
```

## Horological & Multiphysics Architecture

- **Horological Core ([`horology.rs`](file:///Users/jemanuel/projects/frankensim/apps/pocket-watch/engine-watch/src/horology.rs))**:
  - **Mainspring Barrel**: $78$ teeth ($R = 7.80\text{ mm}$), housing a $520\text{ mm}$ cold-drawn high-carbon steel ribbon (`music-wire-spring-moduli-fuchs-1968`) with a $42\text{ h}$ power reserve ($8.5$ arbor turns) and a $78$-tooth winding ratchet + recoil click.
  - **Going Train**: Center wheel ($80$ teeth / pinion $10$, $1\text{ rev/h}$), third wheel ($75$ teeth / pinion $10$), fourth wheel ($80$ teeth / pinion $10$, $1\text{ rev/min}$), and $12:1$ motion works + power-reserve differential + $29.53\text{ d}$ moonphase complication.
  - **1-Minute Titanium Tourbillon & Swiss Lever Escapement**: Grade 5 Titanium carriage (`titanium-ti6al4v-nist-thermal-reference`) rotating at $1\text{ rpm}$ around the fixed fourth wheel, driving an $8$-leaf escape pinion and $15$-tooth club-tooth steel escape wheel (`aisi-52100-cvm-hot-hardness`), a Swiss lever pallet fork with monocrystalline synthetic ruby entry/exit pallet stones (`sapphire-al2o3-crystran-reference`, $\pm 5.25^\circ$ banking lift), and a coaxial Glucydur balance wheel (`beryllium-copper-c17200-ngk-reference`) with $16$ gold poising screws and a $14$-turn Breguet overcoil hairspring ($18{,}000\text{ vph} = 2.500000\text{ Hz}$).

- **12-Node Coupled Multiphysics ([`multiphysics.rs`](file:///Users/jemanuel/projects/frankensim/apps/pocket-watch/engine-watch/src/multiphysics.rs) & [`materials.rs`](file:///Users/jemanuel/projects/frankensim/apps/pocket-watch/engine-watch/src/materials.rs))**:
  - **Virtual Oxy-Butane Blowtorch ($2250\text{ K}$)**: Heats targeted watch sectors via a 3D Gaussian plume, driving thermal expansion $\epsilon_{\text{th}}(T)$ of the balance rim ($I(T) = I_0(1+\epsilon_{\text{th}})^2$), thermoelastic Young's modulus drift $E(T)/E_0$ of the hairspring, Arrhenius $\text{Fe}_3\text{O}_4$ oxide thin-film growth, and Planck blackbody incandescence.
  - **18k Gold Case Melting & Liquid Metal Fouling**: Sustained blowtorch heating pushes the 18k gold case (`gold-pure-nbs-c447`, $T_{\text{liq}} = 1337.33\text{ K}$, $L_f = 63.7\text{ kJ/kg}$) across its enthalpy phase-change plateau (`MushyPhaseLaw`). As Carman-Kozeny mobility $\mathcal{M}(\phi_{\text{liq}})$ activates, the case wall slumps and molten gold flows inward over the Tourbillon, escapement, and gear train, bogging down and halting the movement.
  - **Cryogenic Liquid Nitrogen ($\text{LN}_2$, $77.36\text{ K}$)**: Extracts heat through Leidenfrost film boiling and nucleate boiling (`nitrogen-transport-lemmon2004`), condensing atmospheric frost. Pouring $\text{LN}_2$ onto a melting watch immediately quenches $\phi_{\text{liq}} \to 0$ and freezes the molten gold in its deformed shape (`LOCKED — ENCASED IN SOLIDIFIED GOLD`). Pouring $\text{LN}_2$ onto an unmelted watch spikes oil viscosity by $>10^4\times$ and pinches brass-to-steel pivot clearances.
  - **Watch Oil Tribology & Capillary Bridging**: Models Vogel-Fulcher-Tammann viscosity $\mu(T)$ and surface tension $\gamma(T)$ (`pennzane-shf-x-2000-bearing-oil`), transitioning from boundary friction to EHL at $1.6\text{ }\mu\text{m}$, while flooding with oil creates Couette drag and capillary meniscus bridges across adjacent Breguet hairspring coils.
  - **Sintered $\text{NdFeB N42}$ Super-Magnet ($B_r = 1.32\text{ T}$)**: Evaluates 3D dipolar magnetic flux density $\mathbf{B}(\mathbf{r})$ (`jinshan-n42-pristine-temperature`), inducing hairspring inter-coil attraction, beat-error bias, eddy-current braking, ferromagnetic remanence $M_{\text{rem}}$ (cleared via the 60 Hz AC degaussing coil), and pallet-fork magnetic pinning above $22\text{ mT}$.

- **48 kHz 3D Virtual Microphone & Modal Acoustics ([`acoustics.rs`](file:///Users/jemanuel/projects/frankensim/apps/pocket-watch/engine-watch/src/acoustics.rs))**:
  - Synthesizes the 3-part horological beat signature (*unlock impact*, *ruby pallet impulse slide*, and *escape tooth drop / banking lock*), winding ratchet clicks, blowtorch combustion roar, and $\text{LN}_2$ Leidenfrost sizzle through a 12-mode structural resonator bank attenuated to the movable 3D studio microphone probe.
