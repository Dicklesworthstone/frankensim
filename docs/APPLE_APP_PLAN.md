# FrankenSim for Apple Platforms

Status: implementation plan for the first native iPhone, iPad, and Mac release.

## Product promise

FrankenSim is a native simulation studio, not a website wrapper and not a gallery
of canned animations. Every runnable experiment is produced locally by the same
bounded Rust kernels exposed by `crates/fs-wasm`. SwiftUI owns navigation,
controls, accessibility, and platform adaptation; Swift Canvas owns interactive
rendering. The product must preserve FrankenSim's central distinction between a
number and a justified claim.

The app does **not** claim that the whole FrankenSim workspace is a finished
general-purpose simulator. It exposes the current laboratory and campaign
surfaces with their real maturity and no-claim boundaries.

## Information architecture

1. **Studio** — one focused experiment with immediate parameter controls, an
   animated native visualization, runtime statistics, and an evidence card.
2. **Laboratory** — the website's 30 kernels grouped into Foundations, Frontier,
   and Deep Kernel. Search and filters remain visible on wide layouts.
3. **Campaigns** — ten evidence-bearing end-to-end campaigns plus the ornithoid,
   frame, and vessel flagships. These are presented as composed workflows, not
   as stronger maturity claims than their source contracts allow.
4. **Atlas** — concise native explanations of the seven layers, Three Colors,
   Five Explicits, Gauntlet, glossary, and architecture.

## Platform behavior

- **iPhone:** single-column focus. The run controls, receipt, and evidence
  boundary follow the live canvas in one predictable reading order. Primary
  controls never require a toolbar menu.
- **iPad portrait:** canvas above a two-column experiment/inspector region.
  Landscape uses a persistent three-column navigation split.
- **Mac Catalyst:** freely resizable window with a native sidebar, central live
  canvas, and persistent inspector. Minimum size prevents unusable layouts but
  there is no maximum size.
- Animations pause when the scene is inactive and honor Reduce Motion. Every
  color-coded claim also has a textual label.

## Native engine boundary

`ios/rust` builds a static library over `fs-wasm`. The C ABI is deliberately
small: run a bounded catalog entry, query the immutable result length, and copy
individual values before the next run. Calls are panic-contained and results
are thread-local, so no Rust allocation or ownership crosses the ABI.

The result packet begins with six finite `f64` values:

`[schema, experiment_id, shape, width, height, frames, payload…]`

Unknown schemas, malformed dimensions, non-finite metadata, and over-budget
payloads are refused in Swift before rendering.

## Visual language

The family resemblance is the friendly Franken character, deep laboratory
black, restrained emerald/cyan/violet energy, rounded Apple typography, and
small monospaced labels only where code or measurements justify them. Color is
used to explain simulation state and evidence—not as decorative rainbow noise.

The App Store icon uses the suite's large upper-right letter badge (`S`) and
turns the character into a simulation scientist holding a wireframe geometry
and optimized lattice.

## Acceptance gates

- The generated Xcode project builds for arm64 iPhone, arm64 simulator, and Mac
  Catalyst.
- Every catalog entry either runs a real Rust kernel or is visibly marked as
  explanatory content; there are no fake progress animations or random data.
- App launch and experiment switching remain responsive because computation is
  off the main actor.
- iPhone portrait, iPad portrait/landscape, and small/large Mac windows are
  visually inspected.
- Privacy manifest declares no collection and the app performs no network
  requests.

## Fresh-eyes review before implementation

The first pass was narrowed after comparing the engine README, plan, `fs-wasm`
contract, website catalog, and the existing Franken Apple apps:

- A single flashy flagship would hide the real breadth, so the complete 30-lab
  catalog and all campaign identities are first-class from day one.
- Recreating kernels in Swift would create a second implementation and violate
  the website's strongest promise. The Rust bridge is therefore foundational.
- A giant scrolling dashboard performs poorly on iPhone and Mac. The app uses a
  focused studio plus adaptive sidebars/inspectors instead.
- Evidence colors need text, provenance, and no-claim context. A glowing badge
  alone would misrepresent scientific maturity.
