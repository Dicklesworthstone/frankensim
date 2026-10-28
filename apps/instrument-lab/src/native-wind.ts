import "./native-wind.css";
// Canonical native fixtures, bundled verbatim by Vite; no copied physics defaults.
import reedExampleSource from "../../../crates/fs-couple/examples/reed-two-medium-tone-hole.performance?raw";
import flueExampleSource from "../../../crates/fs-couple/examples/flue-nolle.performance?raw";
import flueDgSource from "../../../crates/fs-couple/examples/flue-dg-baffled.performance?raw";
import materialReedSource from "../../../crates/fs-couple/examples/material_modal_valve/paired-walls.performance?raw";
import materialReedDgSource from "../../../crates/fs-couple/examples/material_modal_valve/paired-dg-chimney.performance?raw";
import stringExampleSource from "../../../crates/fs-couple/examples/material_strings/pluck-and-bow.material-strings?raw";
import observedStringSource from "../../../crates/fs-couple/examples/material_strings/rayleigh-observed.material-strings?raw";
import layeredStringSource from "../../../crates/fs-couple/examples/material_strings/relaxing-layered-strings.material-strings?raw";
import shearStringSource from "../../../crates/fs-couple/examples/material_strings/shear-layer-bow.material-strings?raw";
import { nativeWindStatus, renderNativePerformance, type NativeWindStatus, type NativeWindAttachment, type NativeWindGasSources, type NativeStringProfiles } from "./native-wind-client.mjs";
import { el } from "./util";
import { createNativeGasControls } from "./native-wind-gases";
import { createNativeAtmosphereControls } from "./native-atmosphere";
import { createNativeStringProfiles } from "./native-string-profiles";
import { createNativeComparisonPanel } from "./native-comparison";
import { MAX_REQUEST_BYTES, editableAttachmentText, editedAttachmentBytes, serializeNativeExperiment, validateInput } from "./native-experiment.mjs";

const EXAMPLES = [
  {
    name: "reed-two-medium-tone-hole.performance",
    label: "Reed — two-medium tone hole",
    renderer: "wind",
    source: reedExampleSource,
    decimate: false,
    intro: "Starter: two-medium tone hole — a short boundary diagnostic with a reed and three tube sections in different gas regions. Change the gas mixtures, temperature, geometry or gestures in the description.",
    help: "ambient sets the inlet and chimney gas; a segment’s gas entry overrides its local medium; exterior sets the surrounding gas. The interior and exterior must share a static pressure because this boundary model omits mean flow. This short numerical fixture is not a measured clarinet or oboe.",
  },
  {
    name: "flue-nolle.performance",
    label: "Flute — Nolle/Price jet",
    renderer: "wind",
    source: flueExampleSource,
    decimate: false,
    intro: "Starter: flute jet — a short blowing-pressure phrase through a cylindrical tube, using the native Nolle/Price jet transfer. Change the gas, tube, mouth geometry or pressure gesture in the description.",
    help: "gas owns the internal medium; exterior sets the surrounding gas. jet_transfer selects the short-channel Nolle/Price source. The native description owns the mechanical and output clocks, and the jet and acoustic mouth are solved together. This synthetic example is not a measured flute.",
  },
  {
    name: "flue-dg-baffled.performance",
    label: "Flute — resolved duct and material atmosphere",
    renderer: "wind",
    source: flueDgSource,
    decimate: true,
    intro: "Starter: a flute jet drives a spatial duct with a radiating outlet and mouth. Change the bore and exterior gas compositions independently; the native pressure observer propagates each source through the declared exterior atmosphere.",
    help: "dg_exterior supplies density, sound speed, viscosity, conductivity and heat capacities to the shared Stokes–Kirchhoff propagation law. observer_gas selects that law and its numerical budgets; lossless is an explicit comparison. Native decimation observes the unchanged 96 kHz mechanics at 48 kHz. This short integration fixture uses an interim jet model and is not a measured flute.",
  },
  {
    name: "material_modal_valve/paired-walls.performance",
    label: "Material reeds — paired blades and compliant duct walls",
    renderer: "wind",
    source: materialReedSource,
    decimate: false,
    intro: "Starter: two material-derived reed blades coupled to an upstream tract, downstream bore and material walls. Attach blade.matpack and wall.matpack from the canonical material-modal example. Change each blade region or wall source independently.",
    help: "The original meshes, thickness regions, material stiffness and loss determine blade motion, aperture flow and external microphone pressure. This source already observes its 768 kHz mechanics at 48 kHz; leave additional decimation off. Its manufactured materials are not measured cane or a qualified clarinet/oboe specimen.",
  },
  {
    name: "material_modal_valve/paired-dg-chimney.performance",
    label: "Material reeds — resolved duct, chimney and regional gases",
    renderer: "wind",
    source: materialReedDgSource,
    decimate: true,
    intro: "Starter: the same material-derived paired blades drive a spatially resolved duct with a side chimney. Attach blade.matpack. The bore, chimney and exterior have separately declared gas sources and feed the original radiating boundary.",
    help: "This canonical source runs at 96 kHz, with native decimation selected for 48 kHz playback. observer_gas propagates the outlet and chimney sound through dg_exterior using the shared Stokes–Kirchhoff law. Change the material pack or gas data in the source; the same native solver derives the pressure output.",
  },
  {
    name: "material_strings/pluck-and-bow.material-strings",
    label: "Strings — pluck and bow with a shared body",
    renderer: "strings",
    source: stringExampleSource,
    decimate: false,
    intro: "Starter: one material-derived string, body and cavity, with a physical preload and bow gesture. Attach wire.fsmatpk and body.fsmatpk using the source-pack instructions linked below, or open your own complete material-string performance.",
    help: "The string performer accepts the same full material cards as the reed and wall solvers, along with legacy claim packs. Replace an attached source at its declared path to change the solved material without browser coefficients. Optional profiles add body geometry, fingers/frets and contact heat. This is a numerical integration fixture, not a measured violin, cello or bass.",
  },
  {
    name: "material_strings/rayleigh-observed.material-strings",
    label: "Strings — material plate and exterior microphone",
    renderer: "strings",
    source: observedStringSource,
    decimate: false,
    intro: "Starter: a plucked material string drives a regional orthotropic plate, with reciprocal acoustic loading and a microphone at its declared position. Attach wire.fsmatpk, directional-body.fsmatpk and orthotropic.body, then select orthotropic.body as the body profile below.",
    help: "The native observer integrates signed plate motion over source-to-microphone paths through the original gas. Change material sources or atmosphere data and render again. The body uses an infinite baffle and synthetic specimen data; this flat-plate fixture does not represent a complete measured violin, cello or bass.",
  },
  {
    name: "material_strings/relaxing-layered-strings.material-strings",
    label: "Strings — bonded layers with original relaxation spectra",
    renderer: "strings",
    source: layeredStringSource,
    decimate: false,
    intro: "Starter: layered strings with independent source relaxation spectra. Attach each original material pack named in the performance. Change layer materials and geometry in the data; the shared constitutive solver retains each layer's history.",
    help: "This perfectly bonded concentric-section model retains original per-layer moduli and relaxation times. It does not add sliding helical windings. Keep the stated end-load or displacement comparison and pressure scale when comparing material substitutions.",
  },
  {
    name: "material_strings/shear-layer-bow.material-strings",
    label: "Strings — material shear layer and distributed bow contact",
    renderer: "strings",
    source: shearStringSource,
    decimate: false,
    intro: "Starter: a material shear layer drives distributed bow contact on the original string/body/cavity system. Attach wire.fsmatpk, body.fsmatpk and manufactured-shear-layer.fsmatpk from the canonical material-string example.",
    help: "Layer shear modulus, yield stress and plastic viscosity control retained elastic/plastic history and the actual pressure output. The native instructions describe geometric stretching, longitudinal waves and optional contact heat. Heat feedback requires the same materials to supply caloric and transport data; this isothermal fixture does not invent them.",
  },
];
const exampleUrl = (name: string): string =>
  `https://github.com/Dicklesworthstone/frankensim/blob/main/crates/fs-couple/examples/${name}`;
const LOCAL_FILE_LIMIT = 1_048_576; // Bound file reading even on a static host without a renderer.
const message = (error: unknown): string => error instanceof Error ? error.message : String(error);

/** A source editor and native render client. Synthesis stays in the Rust renderer. */
export function createNativeWindPanel(): HTMLElement {
  let selectedExample = EXAMPLES[0]!;
  const renderer = el("select", {id: "native-performance-renderer", class: "btn"},
    el("option", {value: "wind"}, "Wind — reeds, jets and material valves"),
    el("option", {value: "strings"}, "Strings — materials, bow/contact and shared body"));
  const selectedRenderer = (): "wind" | "strings" => renderer.value === "strings" ? "strings" : "wind";
  const rendererAvailable = (next: NativeWindStatus | null): boolean => Boolean(next?.available
    && (selectedRenderer() === "strings" ? next.renderers?.strings === true : next.renderers?.wind !== false));
  const starter = el("select", { id: "native-wind-starter", class: "btn" },
    ...EXAMPLES.map((example) => el("option", { value: example.name }, example.label)));
  const introduction = el("p", { class: "nw-intro" }, selectedExample.intro);
  const sourceHelp = el("p", { id: "native-wind-source-help" }, selectedExample.help);
  const canonical = el("a", { href: exampleUrl(selectedExample.name), target: "_blank", rel: "noopener" }, "Canonical example");
  const source = el("textarea", { id: "native-wind-source", class: "nw-source", rows: "15",
    spellcheck: "false", autocapitalize: "off", autocomplete: "off", wrap: "off",
    "aria-describedby": "native-wind-source-help" });
  source.value = selectedExample.source;
  // A textarea normalizes CRLF. Keep the loaded UTF-8 text until the user edits it.
  let loadedSource = selectedExample.source;
  let loadedEditorValue = source.value;
  const currentSource = (): string => source.value === loadedEditorValue ? loadedSource : source.value;
  const sourceInfo = el("span", { class: "small muted" });
  const connection = el("span", { class: "nw-connection", role: "status", "aria-live": "polite" }, "Checking local renderer…");
  const result = el("p", { class: "nw-result small", role: "status", "aria-live": "polite", "aria-atomic": "true" });
  const render = el("button", { type: "button", class: "btn primary" }, "Render WAV");
  const cancel = el("button", { type: "button", class: "btn" }, "Cancel render");
  const check = el("button", { type: "button", class: "btn" }, "Check renderer");
  const open = el("button", { type: "button", class: "btn" }, "Open performance…");
  const openExperiment = el("button", {type:"button", class:"btn"}, "Open experiment…");
  const saveExperiment = el("button", {type:"button", class:"btn"}, "Save experiment");
  const experimentFile = el("input", {type:"file", accept:".json,application/json",
    "aria-label":"Open a complete native experiment with its material sources"});
  experimentFile.hidden = true;
  const reset = el("button", { type: "button", class: "btn" }, "Reset example");
  const attach = el("button", { type: "button", class: "btn" }, "Attach source files…");
  const attachmentFile = el("input", { type: "file", multiple: "",
    "aria-label": "Attach material packs or other referenced native source files" });
  attachmentFile.hidden = true;
  const decimate = el("input", { type: "checkbox", id: "native-wind-decimate" });
  const decimation = el("label", { for: "native-wind-decimate", class: "small nw-decimation" }, decimate,
    " Use native decimation for a performance above 48 kHz (mechanics keep their original clock)");
  const attachmentList = el("div", { class: "nw-attachments" });
  const attachmentInfo = el("p", { class: "small muted" });
  const attachmentPanel = el("section", { "aria-label": "Attached source files" },
    attachmentInfo, attachmentList);
  type Attached = NativeWindAttachment & { sha256: string };
  let attachments: Attached[] = [];
  let replacing: Attached | null = null;
  const replacementFile = el("input", {type:"file", "aria-label":"Replace an attached source without changing its path"});
  replacementFile.hidden = true;
  let editingAttachment: Attached | null = null;
  let attachmentEditorInitialValue = "";
  const attachmentEditorTitle = el("label", {for:"native-attachment-source"});
  const attachmentSource = el("textarea", {id:"native-attachment-source", class:"nw-source", rows:"12",
    spellcheck:"false", autocapitalize:"off", autocomplete:"off", wrap:"off"});
  const applyAttachment = el("button", {type:"button", class:"btn primary"}, "Apply source edit");
  const cancelAttachment = el("button", {type:"button", class:"btn"}, "Cancel edit");
  const attachmentEditor = el("section", {"aria-label":"Edit attached physical source"},
    attachmentEditorTitle,
    el("p", {class:"small muted"}, "Edit original profile, geometry or control records. Apply keeps this file's path and all profile associations; the native renderer checks the physical data. Binary material, interface and model packs use Replace."),
    attachmentSource, el("div", {class:"btn-row"}, applyAttachment, cancelAttachment));
  attachmentEditor.hidden = true;
  const attachmentIdentity = (): string => JSON.stringify(attachments.map(file => [file.path, file.sha256]));
  const gas = createNativeGasControls(() => edited());
  const profiles = createNativeStringProfiles(() => edited());
  const atmosphere = createNativeAtmosphereControls(currentSource, next => {
    source.value = next; loadedSource = next; loadedEditorValue = source.value;
    edited();
  }, () => controls());
  const file = el("input", { type: "file", accept: ".performance,.material-strings,.txt,text/plain",
    "aria-label": "Open a native physical performance" });
  file.hidden = true;
  const audio = el("audio", { controls: "", preload: "metadata", "aria-label": "Your native physical render" });
  const download = el("a", { class: "btn", download: "wind-render.wav" }, "Download WAV");
  const renderedFrom = el("p", { class: "small muted" });
  const renderedFiles = el("p", { class: "small muted nw-file-identities" });
  const playback = el("div", { class: "nw-playback" }, audio, download, renderedFrom, renderedFiles);
  playback.hidden = true;
  cancel.hidden = true;
  const comparison = createNativeComparisonPanel(() => audio.pause());
  audio.addEventListener("play", () => comparison.stop());

  let status: NativeWindStatus | null = null;
  let checking: AbortController | null = null;
  let rendering: AbortController | null = null;
  let loadingFile = false;
  let previewUrl: string | null = null;
  let renderedSource: string | null = null;
  let renderedAttachmentIdentity: string | null = null;
  let renderedDecimate = false;
  let renderedGasIdentity: string | null = null;
  let renderedProfileIdentity: string | null = null;
  let renderedRenderer: string | null = null;
  let sourceName = selectedExample.name;

  const controls = (): void => {
    const active = rendering !== null;
    const otherEdit = active || loadingFile || editingAttachment !== null;
    const locked = otherEdit || atmosphere.hasDraft();
    render.disabled = locked || checking !== null || !rendererAvailable(status) || Boolean(status?.busy) || !source.value.trim();
    cancel.hidden = !active;
    check.disabled = active || checking !== null;
    open.disabled = locked;
    openExperiment.disabled = locked;
    saveExperiment.disabled = locked || !source.value.trim();
    reset.disabled = locked;
    starter.disabled = locked;
    renderer.disabled = locked;
    attach.disabled = locked;
    decimate.disabled = locked;
    gas.setDisabled(locked);
    atmosphere.setDisabled(otherEdit);
    profiles.element.hidden = selectedRenderer() !== "strings";
    profiles.setDisabled(locked || selectedRenderer() !== "strings");
    attachmentList.querySelectorAll<HTMLInputElement | HTMLButtonElement>("input,button")
      .forEach(control => { control.disabled = locked; });
    source.readOnly = locked;
    source.setAttribute("aria-busy", String(active));
    attachmentSource.readOnly = active || loadingFile;
    applyAttachment.disabled = active || loadingFile;
    cancelAttachment.disabled = active || loadingFile;
    render.textContent = active ? "Rendering…" : "Render WAV";
  };
  const edited = (): void => {
    atmosphere.refresh(currentSource());
    sourceInfo.textContent = `${sourceName} · ${new TextEncoder().encode(currentSource()).byteLength.toLocaleString()} bytes`;
    if (renderedSource !== null) {
      renderedFrom.textContent = currentSource() === renderedSource && attachmentIdentity() === renderedAttachmentIdentity && decimate.checked === renderedDecimate && gas.identity() === renderedGasIdentity && selectedRenderer() === renderedRenderer && profiles.identity() === renderedProfileIdentity
        ? "This WAV was rendered from the performance and source files currently shown."
        : "This WAV is from your previous performance or source files. Render again to hear your edits.";
    }
    controls();
  };
  const clearPreview = (): void => {
    audio.pause();
    audio.removeAttribute("src");
    audio.load();
    download.removeAttribute("href");
    playback.hidden = true;
    renderedSource = null;
    renderedAttachmentIdentity = null;
    renderedGasIdentity = null;
    renderedProfileIdentity = null;
    renderedRenderer = null;
    renderedFiles.textContent = "";
    if (previewUrl !== null) URL.revokeObjectURL(previewUrl);
    previewUrl = null;
  };
  const showResult = (text: string, refused = false): void => {
    result.textContent = text;
    result.classList.toggle("nw-refused", refused);
  };

  const checkConnection = async (): Promise<void> => {
    checking?.abort();
    const controller = new AbortController();
    checking = controller;
    status = null;
    connection.textContent = "Checking local renderer…";
    controls();
    const timeout = window.setTimeout(() => controller.abort(), 5000);
    try {
      const next = await nativeWindStatus(controller.signal);
      if (checking !== controller) return;
      status = next;
      connection.textContent = !rendererAvailable(next) ? "Selected native renderer unavailable"
        : next.busy ? "Local renderer is busy" : "Local renderer ready";
      if (!rendererAvailable(next)) {
        showResult("Configure the selected native executable in the local Instrument Lab server to render and listen. You can edit the source here.");
      } else if (next.busy) {
        showResult("Another render is running. Check the renderer again when it finishes.");
      } else {
        showResult(`Ready to render. The local render limit is ${Math.ceil(next.timeoutMs / 1000)} seconds.`);
      }
    } catch {
      if (checking !== controller) return;
      connection.textContent = "Native renderer unavailable";
      showResult("This host has no native renderer connection. Open the local Instrument Lab server to render and listen.");
    } finally {
      window.clearTimeout(timeout);
      if (checking === controller) {
        checking = null;
        controls();
      }
    }
  };

  render.addEventListener("click", () => {
    if (!status || !rendererAvailable(status) || status.busy || rendering || checking || loadingFile) return;
    const connected = status;
    const submitted = currentSource();
    const submittedFiles = attachments.map(file => ({ ...file }));
    const submittedAttachmentIdentity = attachmentIdentity();
    const submittedDecimate = decimate.checked;
    const submittedRenderer = selectedRenderer();
    let submittedGas: NativeWindGasSources | undefined;
    let submittedProfiles: NativeStringProfiles | undefined;
    try { submittedGas = gas.read(); submittedProfiles = submittedRenderer === "strings" ? profiles.read() : undefined; }
    catch (error) { showResult(message(error), true); return; }
    const submittedGasIdentity = gas.identity();
    const submittedProfileIdentity = profiles.identity();
    const controller = new AbortController();
    rendering = controller;
    audio.pause();
    comparison.stop();
    showResult("Rendering your performance from its material, gas and geometry data…");
    controls();
    let timedOut = false;
    const timeout = window.setTimeout(() => {
      timedOut = true;
      controller.abort();
    }, connected.timeoutMs + 5000);
    void (async () => {
      try {
        const native = await renderNativePerformance(submitted, submittedFiles, connected, {renderer:submittedRenderer,
          decimate:submittedDecimate, gasSources:submittedGas, stringProfiles:submittedProfiles}, controller.signal);
        if (rendering !== controller || controller.signal.aborted) return;
        let details = submittedFiles.length
          ? "Rendered source files: " + submittedFiles.map(file => `${file.path} · ${file.bytes.byteLength.toLocaleString()} bytes · SHA-256 ${file.sha256}`).join("\n")
          : "This render used no attached source files.";
        if (submittedGas) details += "\nGas associations: " + submittedGas.components
          .map(component => `${component.speciesPath} / ${component.nasa9Path} / ${component.transportPath}`).join("; ")
          + ` · phase ${submittedGas.phase} · water ${submittedGas.waterSpecies ?? "none"}`;
        if (submittedProfiles) details += "\nString profiles: " + Object.entries(submittedProfiles)
          .map(([role, path]) => `${role}: ${path}`).join("; ");
        await comparison.setCurrent(native, sourceName,
          `Performer: ${submittedRenderer} · native decimation: ${submittedDecimate}\n${details}\n\nOriginal performance:\n${submitted}`, controller.signal);
        if (rendering !== controller || controller.signal.aborted) return;
        const nextUrl = URL.createObjectURL(native.wave);
        clearPreview();
        previewUrl = nextUrl;
        renderedSource = submitted;
        renderedAttachmentIdentity = submittedAttachmentIdentity;
        renderedDecimate = submittedDecimate;
        renderedGasIdentity = submittedGasIdentity;
        renderedProfileIdentity = submittedProfileIdentity;
        renderedRenderer = submittedRenderer;
        renderedFiles.textContent = details;
        audio.src = nextUrl;
        download.href = nextUrl;
        download.download = `${sourceName.replace(/\.[^.]*$/, "") || "physical"}-render.wav`;
        playback.hidden = false;
        edited();
        showResult("Render complete. Listen above or keep this result as reference A to compare a material or atmosphere change below.");
      } catch (error) {
        if (rendering !== controller) return;
        showResult(controller.signal.aborted
          ? timedOut ? "The render timed out. Check the local renderer before trying again." : "Render cancelled."
          : message(error), !controller.signal.aborted || timedOut);
      } finally {
        window.clearTimeout(timeout);
        if (rendering === controller) {
          rendering = null;
          controls();
        }
      }
    })();
  });
  cancel.addEventListener("click", () => {
    rendering?.abort();
    showResult("Cancelling render…");
  });
  check.addEventListener("click", () => void checkConnection());
  source.addEventListener("input", edited);
  open.addEventListener("click", () => file.click());
  const loadExample = (): void => {
    source.value = selectedExample.source;
    loadedSource = selectedExample.source;
    loadedEditorValue = source.value;
    sourceName = selectedExample.name.split("/").pop()!;
    renderer.value = selectedExample.renderer;
    introduction.textContent = selectedExample.intro;
    sourceHelp.textContent = selectedExample.help;
    canonical.href = exampleUrl(selectedExample.name);
    canonical.hidden = false;
    attachments = [];
    decimate.checked = selectedExample.decimate;
    gas.reset();
    profiles.reset();
    drawAttachments();
    edited();
    showResult("Example loaded. Edit the source or render it.");
  };
  reset.addEventListener("click", () => { loadExample(); void checkConnection(); });
  starter.addEventListener("change", () => {
    if (rendering || loadingFile) return;
    const selected = EXAMPLES.find((example) => example.name === starter.value);
    if (!selected) return;
    selectedExample = selected;
    loadExample();
    void checkConnection();
  });
  renderer.addEventListener("change", () => { edited(); void checkConnection(); });
  const drawAttachments = (): void => {
    attachmentInfo.textContent = attachments.length
      ? `${attachments.length} source files · ${attachments.reduce((sum, file) => sum + file.bytes.byteLength, 0).toLocaleString()} bytes. Paths must match the performance exactly.`
      : "Attach the original material packs or artifacts referenced by your performance. Keep their source paths; edit a path below when the description uses a subdirectory.";
    attachmentList.replaceChildren();
    gas.setAttachments(attachments.map(file => file.path));
    profiles.setAttachments(attachments.map(file => file.path));
    for (const attached of attachments) {
      const path = el("input", { type: "text", value: attached.path, spellcheck: "false",
        "aria-label": `Source path for ${attached.path}` });
      path.addEventListener("input", () => {
        attached.path = path.value; gas.setAttachments(attachments.map(file => file.path));
        profiles.setAttachments(attachments.map(file => file.path)); edited();
      });
      const editText = el("button", {type:"button", class:"btn"}, "Edit text…");
      editText.addEventListener("click", () => {
        if (rendering || loadingFile || editingAttachment) return;
        try {
          const text = editableAttachmentText(attached.bytes);
          attachmentSource.value = text;
          attachmentEditorInitialValue = attachmentSource.value;
          editingAttachment = attached;
          attachmentEditorTitle.textContent = `Editing ${attached.path}`;
          attachmentEditor.hidden = false;
          controls(); attachmentSource.focus();
          showResult("Apply or cancel the attached-source edit before rendering or saving.");
        } catch (error) { showResult(message(error), true); }
      });
      const replace = el("button", {type:"button", class:"btn"}, "Replace…");
      replace.addEventListener("click", () => {
        if (rendering || loadingFile) return;
        replacing = attached; replacementFile.click();
      });
      const remove = el("button", { type: "button", class: "btn" }, "Remove");
      remove.addEventListener("click", () => {
        attachments = attachments.filter(file => file !== attached);
        drawAttachments(); edited();
      });
      attachmentList.append(el("div", { class: "nw-attachment" }, path,
        el("span", { class: "small muted", title: `SHA-256 ${attached.sha256}` },
          `${attached.bytes.byteLength.toLocaleString()} bytes · ${attached.sha256.slice(0, 12)}…`), editText, replace, remove));
    }
    controls();
  };
  const replaceAttachedSource = async (target: Attached, bytes: Uint8Array): Promise<void> => {
    const maximumFile = Math.min(8 * 1024 * 1024, status?.maximumAttachmentFileBytes ?? 8 * 1024 * 1024);
    const maximumTotal = Math.min(32 * 1024 * 1024, status?.maximumAttachmentBytes ?? 32 * 1024 * 1024);
    const otherBytes = attachments.reduce((sum, file) => sum + (file === target ? 0 : file.bytes.byteLength), 0);
    if (!bytes.length || bytes.length > maximumFile || otherBytes + bytes.length > maximumTotal)
      throw new Error("Replacement exceeds the source-file limit. The original source remains attached.");
    const digest = await crypto.subtle.digest("SHA-256", new Uint8Array(bytes));
    const sha256 = Array.from(new Uint8Array(digest), byte => byte.toString(16).padStart(2, "0")).join("");
    const index = attachments.indexOf(target);
    if (index === -1) throw new Error("The source is no longer attached.");
    // Publish bytes and identity together; paths, associations and prior audio stay intact.
    attachments[index] = {path:target.path, bytes, sha256};
    drawAttachments(); edited();
  };
  const closeAttachmentEditor = (): void => {
    editingAttachment = null;
    attachmentEditorInitialValue = "";
    attachmentSource.value = "";
    attachmentEditor.hidden = true;
    controls();
  };
  cancelAttachment.addEventListener("click", () => {
    if (rendering || loadingFile) return;
    closeAttachmentEditor();
    showResult("Edit cancelled. The original attached source is unchanged.");
  });
  applyAttachment.addEventListener("click", () => {
    const target = editingAttachment;
    if (!target || rendering || loadingFile) return;
    let bytes: Uint8Array;
    try {
      bytes = attachmentSource.value === attachmentEditorInitialValue ? target.bytes
        : editedAttachmentBytes(target.bytes, attachmentSource.value);
    }
    catch (error) { showResult(message(error), true); return; }
    if (bytes === target.bytes) {
      closeAttachmentEditor();
      showResult("No source change. The original bytes and line endings are retained.");
      return;
    }
    loadingFile = true; controls();
    void replaceAttachedSource(target, bytes).then(() => {
      closeAttachmentEditor();
      showResult(`Updated ${target.path}. Render again to hear the new physical configuration; reference A is unchanged.`);
    }).catch(error => {
      showResult(`Edit was not applied: ${message(error)}`, true);
    }).finally(() => { loadingFile = false; controls(); });
  });
  decimate.addEventListener("change", edited);
  attach.addEventListener("click", () => attachmentFile.click());
  openExperiment.addEventListener("click", () => experimentFile.click());
  saveExperiment.addEventListener("click", () => {
    if (rendering || loadingFile) return;
    try {
      const packet = serializeNativeExperiment(currentSource(), attachments, {
        renderer:selectedRenderer(), decimate:decimate.checked, gasSources:gas.read(),
        stringProfiles:selectedRenderer() === "strings" ? profiles.read() : undefined,
      });
      // Use the service's same envelope checks before saving an unusable packet.
      // Native physics admission remains at render time.
      validateInput(new TextEncoder().encode(packet), true);
      const url = URL.createObjectURL(new Blob([packet], {type:"application/json"}));
      const link = el("a", {href:url, download:sourceName.replace(/(?:\.native)?\.json$|\.(?:performance|material-strings)$/, "") + ".native.json"});
      link.click();
      window.setTimeout(() => URL.revokeObjectURL(url), 1000);
      showResult("Experiment saved with the original performance, source files, profiles and gas associations. Open it to continue this material comparison.");
    } catch (error) { showResult(message(error), true); }
  });
  experimentFile.addEventListener("change", () => {
    const selected = experimentFile.files?.[0];
    experimentFile.value = "";
    if (!selected || rendering || loadingFile) return;
    if (selected.size > MAX_REQUEST_BYTES) {
      showResult(`Choose an experiment no larger than ${MAX_REQUEST_BYTES.toLocaleString()} bytes.`, true);
      return;
    }
    loadingFile = true; controls();
    void selected.arrayBuffer().then(async buffer => {
      const experiment = validateInput(new Uint8Array(buffer), true);
      const text = new TextDecoder("utf-8", {fatal:true, ignoreBOM:true}).decode(experiment.source);
      const nextAttachments = await Promise.all(experiment.attachments.map(async attached => {
        const digest = await crypto.subtle.digest("SHA-256", new Uint8Array(attached.bytes));
        const sha256 = Array.from(new Uint8Array(digest), byte => byte.toString(16).padStart(2, "0")).join("");
        return {...attached, sha256};
      }));
      // Commit the complete imported experiment only after every byte and
      // association was admitted. A bad file leaves the current model intact.
      attachments = nextAttachments;
      source.value = text; loadedSource = text; loadedEditorValue = source.value;
      sourceName = selected.name;
      renderer.value = experiment.renderer ?? "wind";
      decimate.checked = experiment.decimate;
      drawAttachments();
      gas.load(experiment.gasSources);
      profiles.load(experiment.stringProfiles);
      starter.value = ""; canonical.hidden = true;
      introduction.textContent = `Loaded ${selected.name}: the complete performance and ${attachments.length} original source files.`;
      sourceHelp.textContent = "This experiment keeps its original performer, mechanical clock, pressure scale, material paths, specimen profiles and gas-source associations. Replace an attached source or edit the performance to run another physical configuration.";
      edited();
      await checkConnection();
      showResult("Experiment loaded. Render to hear it; earlier comparison audio is unchanged.");
    }).catch(error => {
      showResult(`Experiment was not loaded: ${message(error)} Your current sources are unchanged.`, true);
    }).finally(() => { loadingFile = false; controls(); });
  });
  attachmentFile.addEventListener("change", () => {
    const selected = Array.from(attachmentFile.files ?? []);
    attachmentFile.value = "";
    if (!selected.length || rendering || loadingFile) return;
    const maxFiles = Math.min(32, status?.maximumAttachments ?? 32);
    const maxFileBytes = Math.min(8 * 1024 * 1024, status?.maximumAttachmentFileBytes ?? 8 * 1024 * 1024);
    const maxBytes = Math.min(32 * 1024 * 1024, status?.maximumAttachmentBytes ?? 32 * 1024 * 1024);
    const total = attachments.reduce((sum, file) => sum + file.bytes.byteLength, 0)
      + selected.reduce((sum, file) => sum + file.size, 0);
    if (attachments.length + selected.length > maxFiles || total > maxBytes
        || selected.some(file => !file.size || file.size > maxFileBytes)) {
      showResult(`Choose at most ${maxFiles} nonempty files, ${maxFileBytes.toLocaleString()} bytes per file and ${maxBytes.toLocaleString()} bytes combined.`, true);
      return;
    }
    loadingFile = true; controls();
    void Promise.all(selected.map(async file => {
      const buffer = await file.arrayBuffer();
      const digest = await crypto.subtle.digest("SHA-256", buffer);
      const sha256 = Array.from(new Uint8Array(digest), byte => byte.toString(16).padStart(2, "0")).join("");
      return {path:file.name, bytes:new Uint8Array(buffer), sha256};
    })).then(added => {
      attachments.push(...added);
      drawAttachments(); edited();
      showResult("Source files attached. Check that each path matches the performance, then render.");
    }).catch(() => {
      showResult("Source files could not be read. Your previous attachments are unchanged.", true);
    }).finally(() => { loadingFile = false; controls(); });
  });
  replacementFile.addEventListener("change", () => {
    const selected = replacementFile.files?.[0], target = replacing;
    replacementFile.value = ""; replacing = null;
    if (!selected || !target || !attachments.includes(target) || rendering || loadingFile) return;
    const maximumFile = Math.min(8 * 1024 * 1024, status?.maximumAttachmentFileBytes ?? 8 * 1024 * 1024);
    const maximumTotal = Math.min(32 * 1024 * 1024, status?.maximumAttachmentBytes ?? 32 * 1024 * 1024);
    const otherBytes = attachments.reduce((sum, file) => sum + (file === target ? 0 : file.bytes.byteLength), 0);
    if (!selected.size || selected.size > maximumFile || otherBytes + selected.size > maximumTotal) {
      showResult("Replacement exceeds the source-file limit. The original material remains attached.", true);
      return;
    }
    loadingFile = true; controls();
    void selected.arrayBuffer().then(buffer => replaceAttachedSource(target, new Uint8Array(buffer)))
      .then(() => showResult(`Replaced ${target.path}. Render again to hear the new source; reference A is unchanged.`))
      .catch(error => showResult(`The replacement was not applied: ${message(error)}`, true))
      .finally(() => { loadingFile = false; controls(); });
  });
  file.addEventListener("change", () => {
    const selected = file.files?.[0];
    file.value = "";
    if (!selected || rendering || loadingFile) return;
    const limit = Math.min(LOCAL_FILE_LIMIT, status?.maximumSourceBytes ?? LOCAL_FILE_LIMIT);
    if (selected.size > limit) {
      showResult(`Choose a performance file no larger than ${limit.toLocaleString()} bytes.`, true);
      return;
    }
    loadingFile = true;
    controls();
    void selected.arrayBuffer().then((bytes) => {
      const text = new TextDecoder("utf-8", { fatal: true, ignoreBOM: true }).decode(bytes);
      source.value = text;
      loadedSource = text;
      loadedEditorValue = source.value;
      sourceName = selected.name;
      starter.value = ""; canonical.hidden = true;
      introduction.textContent = `Loaded performance: ${selected.name}. Existing source attachments and profile choices are retained.`;
      showResult("Performance loaded. Edit the source or render it.");
      edited();
    }).catch(() => {
      showResult("This performance could not be read as UTF-8 text. Your previous source is still in the editor.", true);
    }).finally(() => {
      loadingFile = false;
      controls();
    });
  });
  window.addEventListener("pagehide", () => {
    checking?.abort();
    rendering?.abort();
    clearPreview();
    comparison.dispose();
  });
  window.addEventListener("pageshow", (event) => {
    if (event.persisted) void checkConnection();
  });

  const panel = el("section", { class: "panel nw-panel", "aria-labelledby": "native-wind-title" },
    el("div", { class: "panel-head" },
      el("div", {}, el("h3", { id: "native-wind-title" }, "Render a material-driven performance"),
        el("p", { class: "small muted" }, "Edit the source, run the native physics, and listen to your result.")),
      connection),
    el("div", {class: "nw-source-head"}, el("label", {for: "native-performance-renderer"}, "Physical performer"), renderer),
    el("div", { class: "nw-source-head" }, el("label", { for: "native-wind-starter" }, "Starter performance"), starter),
    introduction,
    el("details", { class: "nw-help small" },
      el("summary", {}, "About the source and local renderer"),
      sourceHelp,
      el("p", {}, "Wind and material-string descriptions use their original native loaders. Attach the exact material packs and artifacts named by the source, including subdirectory paths. Changing a pack changes its physical properties; missing data and unsupported states remain native refusals. Keep the same declared pressure scale when comparing renders."),
      el("p", {}, "Rendering uses the local music_render or material_strings executable. The hosted review page can show and edit the source; audio is available when you open it through the local renderer server. String contact-thermal profiles retain their original mechanical clock when pressure is decimated to 48 kHz."),
      el("p", {}, "Open experiment loads a complete native source packet, including its original material files, gas associations and specimen profiles. Save experiment keeps the current configuration in that same format; neither action renders sound. The material-string launcher can prepare these packets directly from a physical assembly."),
      el("p", {}, canonical,
        " · ", el("a", {href: "https://github.com/Dicklesworthstone/frankensim/tree/main/crates/fs-couple/examples/material_strings", target: "_blank", rel: "noopener"}, "String source-pack instructions"),
        " · ", el("a", { href: "https://github.com/Dicklesworthstone/frankensim/tree/main/apps/instrument-lab#native-wind-and-string-renderer", target: "_blank", rel: "noopener" }, "Local renderer setup"))),
    atmosphere.element,
    el("div", { class: "nw-source-head" }, el("label", { for: "native-wind-source" }, "Performance source"),
      sourceInfo),
    source, attachmentPanel, attachmentEditor, profiles.element, gas.element, decimation,
    el("div", { class: "btn-row nw-actions" }, render, cancel, openExperiment, saveExperiment, open, attach, reset, check,
      file, experimentFile, attachmentFile, replacementFile),
    result, playback, comparison.element);
  drawAttachments();
  edited();
  void checkConnection();
  return panel;
}
