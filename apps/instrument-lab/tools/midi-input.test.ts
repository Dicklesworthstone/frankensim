import assert from "node:assert/strict";
import { test } from "node:test";
import { decodeMidiMessage, LiveMidiInput, type MidiChannelMessage, type MidiInputState } from "../src/live/midiInput.ts";

function deferred<T>() {
  let resolve!: (value: T) => void, reject!: (reason: Error) => void;
  const promise = new Promise<T>((yes, no) => { resolve = yes; reject = no; });
  return { promise, resolve, reject };
}

class Input extends EventTarget {
  state = "connected";
  connection = "closed";
  openCalls = 0;
  closeCalls = 0;
  opening: Promise<void> | null = null;
  onmidimessage: ((event: Event) => void) | null = null;
  constructor(readonly id: string, readonly name = id) { super(); }
  async open() {
    ++this.openCalls;
    if (this.opening) await this.opening;
    this.connection = "open";
    return this;
  }
  async close() { ++this.closeCalls; this.connection = "closed"; return this; }
  send(data: number[], timeStamp = 42.25) {
    const event = new Event("midimessage");
    Object.defineProperties(event, { data: { value: new Uint8Array(data) }, timeStamp: { value: timeStamp } });
    this.onmidimessage?.(event);
    this.dispatchEvent(event);
  }
}

class Access extends EventTarget {
  inputs = new Map<string, Input>();
  onstatechange: (() => void) | null = null;
  constructor(...inputs: Input[]) { super(); inputs.forEach((input) => this.inputs.set(input.id, input)); }
  changed() { this.onstatechange?.(); this.dispatchEvent(new Event("statechange")); }
}

function observe(access: Access, requestAccess = async () => access as unknown as MIDIAccess) {
  const messages: { message: MidiChannelMessage; timeStamp: number; inputId: string }[] = [];
  const states: MidiInputState[] = [];
  let resets = 0;
  const input = new LiveMidiInput({ requestAccess,
    onMessage: (message, timeStamp, inputId) => messages.push({ message, timeStamp, inputId }),
    onReset: () => { ++resets; }, onState: (state) => states.push(state) });
  return { input, messages, states, get resets() { return resets; }, get state() { return states.at(-1)!; } };
}

// Let the permission/open promises run without timers or a real browser device.
async function tick() { for (let i = 0; i < 8; ++i) await Promise.resolve(); }

test("decoder preserves channels, release velocity and controller values", () => {
  assert.deepEqual(decodeMidiMessage([0x95, 60, 127]), { kind: "noteon", channel: 5, key: 60, velocity: 127 });
  assert.deepEqual(decodeMidiMessage(new Uint8Array([0x82, 21, 46])), { kind: "noteoff", channel: 2, key: 21, velocity: 46 });
  assert.deepEqual(decodeMidiMessage([0x9f, 127, 0]), { kind: "noteoff", channel: 15, key: 127, velocity: 0 });
  assert.deepEqual(decodeMidiMessage([0xb9, 4, 91]), { kind: "cc", channel: 9, controller: 4, value: 91 });
});

test("decoder refuses malformed, truncated, running-status and unsupported messages", () => {
  for (const data of [null, [], [0x90, 60], [0x90, 60, 1, 2], [60, 90, 1], [0x90, 128, 1],
    [0x90, 60, 128], [0x190, 60, 1], [-1, 60, 1], [0x90, -1, 1], [0x90, 60, -1],
    [0x90, 60.5, 1], [0x90, 60, NaN], [NaN, 60, 1], [0x90, Infinity, 1],
    [0xf8], [0xff, 0, 0], [0xf0, 1, 0xf7], [0xe0, 0, 64], [0xa0, 60, 64], [0xc0, 1], [0xd0, 1]]) {
    assert.equal(decodeMidiMessage(data), null, String(data));
  }
});

test("connect chooses one input, shares pending requests and delivers its receipt timestamp", async () => {
  const a = new Input("a", "Keyboard"), b = new Input("b", "Pads"), access = new Access(a, b);
  let requests = 0;
  const h = observe(access, async () => { ++requests; return access as unknown as MIDIAccess; });
  assert.equal(h.state.status, "disconnected", "constructor publishes a usable initial state");
  const first = h.input.connect();
  assert.equal(h.input.connect(), first);
  await first;
  await h.input.connect();
  assert.equal(requests, 1);
  assert.equal(a.openCalls, 1);
  assert.equal(b.openCalls, 0);
  assert.equal(h.state.status, "connected");
  assert.equal(h.state.selectedId, "a");
  assert.deepEqual(h.state.inputs, [{ id: "a", name: "Keyboard" }, { id: "b", name: "Pads" }]);
  a.send([0x93, 60, 90], 123.5);
  b.send([0x90, 61, 90]);
  a.send([0x90, 60]);
  assert.deepEqual(h.messages, [{ message: { kind: "noteon", channel: 3, key: 60, velocity: 90 }, timeStamp: 123.5, inputId: "a" }]);
  h.input.dispose();
});

test("selection resets held input once and never listens to the old device", async () => {
  const a = new Input("a"), b = new Input("b"), h = observe(new Access(a, b));
  await h.input.connect();
  a.send([0x90, 60, 90]);
  h.input.select("b");
  await h.input.connect();
  assert.equal(h.resets, 1);
  assert.equal(h.state.selectedId, "b");
  a.send([0x80, 60, 0]);
  b.send([0x91, 61, 80]);
  assert.deepEqual(h.messages.map((event) => event.inputId), ["a", "b"]);
  h.input.select("b");
  await h.input.connect();
  assert.equal(h.resets, 1);
  h.input.dispose();
});

test("unplug waits for the selected identity and reconnects its replacement", async () => {
  const a = new Input("a"), b = new Input("b"), access = new Access(a, b), h = observe(access);
  await h.input.connect();
  a.state = "disconnected";
  access.changed();
  assert.equal(h.resets, 1);
  assert.equal(h.state.status, "waiting");
  assert.equal(h.state.selectedId, "a");
  assert.deepEqual(h.state.inputs, [{ id: "b", name: "b" }]);
  b.send([0x90, 60, 80]);
  assert.equal(h.messages.length, 0);
  const replacement = new Input("a", "Keyboard reconnected");
  access.inputs.set("a", replacement);
  access.changed();
  await h.input.connect();
  a.send([0x90, 60, 80]);
  replacement.send([0x90, 62, 80]);
  assert.equal(h.messages.length, 1);
  assert.equal(h.state.status, "connected");
  assert.equal(h.resets, 1);
  h.input.dispose();
});

test("a connected port object replaced under the same ID releases the old source", async () => {
  const a = new Input("a"), access = new Access(a), h = observe(access);
  await h.input.connect();
  const replacement = new Input("a");
  access.inputs.set("a", replacement);
  access.changed();
  await h.input.connect();
  assert.equal(h.resets, 1);
  a.send([0x90, 60, 80]);
  replacement.send([0x90, 61, 80]);
  assert.equal(h.messages.length, 1);
  assert.equal(replacement.openCalls, 1);
  h.input.dispose();
});

test("empty inputs wait for hotplug; an explicit missing selection is not replaced", async () => {
  const access = new Access(), h = observe(access);
  await h.input.connect();
  assert.equal(h.state.status, "waiting");
  assert.equal(h.state.selectedId, null);
  const a = new Input("a");
  access.inputs.set("a", a);
  access.changed();
  await h.input.connect();
  assert.equal(h.state.selectedId, "a");
  h.input.select("missing");
  assert.equal(h.state.status, "waiting");
  h.input.disconnect();
  await h.input.connect();
  assert.equal(h.state.selectedId, "missing");
  assert.equal(h.state.status, "waiting");
  h.input.dispose();
});

test("permission denial can retry without leaving an active subscription", async () => {
  const a = new Input("a"), access = new Access(a);
  let requests = 0;
  const h = observe(access, async () => {
    if (++requests === 1) throw new Error("Permission denied");
    return access as unknown as MIDIAccess;
  });
  await h.input.connect();
  assert.equal(h.state.status, "error");
  assert.match(h.state.message, /Permission denied/);
  assert.equal(a.openCalls, 0);
  await h.input.connect();
  assert.equal(h.state.status, "connected");
  assert.equal(requests, 2);
  h.input.dispose();
});

for (const action of ["disconnect", "dispose"] as const) {
  test(`${action} settles a pending permission request and ignores its late result`, async () => {
    const a = new Input("a"), access = new Access(a), permission = deferred<MIDIAccess>();
    const h = observe(access, () => permission.promise);
    const connecting = h.input.connect();
    await tick();
    h.input[action]();
    await connecting;
    const stateCount = h.states.length;
    permission.resolve(access as unknown as MIDIAccess);
    await tick();
    access.changed();
    a.send([0x90, 60, 80]);
    assert.equal(a.openCalls, 0);
    assert.equal(h.states.length, stateCount);
    assert.equal(h.messages.length, 0);
    assert.equal(h.resets, 1);
    h.input.dispose();
  });
}

test("late permission from a cancelled connection cannot replace a newer connection", async () => {
  const old = new Input("old"), current = new Input("current"), permission = deferred<MIDIAccess>();
  let requests = 0;
  const h = observe(new Access(), () => ++requests === 1 ? permission.promise : Promise.resolve(new Access(current) as unknown as MIDIAccess));
  const first = h.input.connect();
  await tick();
  h.input.disconnect();
  await first;
  await h.input.connect();
  permission.resolve(new Access(old) as unknown as MIDIAccess);
  await tick();
  assert.equal(h.state.selectedId, "current");
  assert.equal(old.openCalls, 0);
  current.send([0x90, 60, 80]);
  assert.equal(h.messages[0]?.inputId, "current");
  h.input.dispose();
});

test("changing input cancels the old open wait and ignores its eventual completion", async () => {
  const a = new Input("a"), b = new Input("b"), opening = deferred<void>();
  a.opening = opening.promise;
  const h = observe(new Access(a, b));
  const connecting = h.input.connect();
  await tick();
  assert.equal(a.openCalls, 1);
  h.input.select("b");
  await connecting;
  await h.input.connect();
  assert.equal(h.state.selectedId, "b");
  const stateCount = h.states.length;
  opening.resolve();
  await tick();
  a.send([0x90, 60, 80]);
  b.send([0x90, 61, 80]);
  assert.equal(h.states.length, stateCount);
  assert.equal(h.messages[0]?.inputId, "b");
  assert.equal(h.resets, 1);
  h.input.dispose();
});

for (const action of ["disconnect", "dispose"] as const) {
  test(`${action} cancels an in-flight open without closing another consumer's port`, async () => {
    const a = new Input("a"), opening = deferred<void>();
    a.opening = opening.promise;
    const h = observe(new Access(a));
    const connecting = h.input.connect();
    await tick();
    assert.equal(a.openCalls, 1);
    h.input[action]();
    await connecting;
    const stateCount = h.states.length;
    opening.resolve();
    await tick();
    a.send([0x90, 60, 80]);
    assert.equal(h.states.length, stateCount);
    assert.equal(h.messages.length, 0);
    assert.equal(h.resets, 1);
    assert.equal(a.closeCalls, 0);
    h.input.dispose();
  });
}

test("open failure is reported once and an explicit retry reuses the permission", async () => {
  const a = new Input("a"), access = new Access(a);
  let requests = 0;
  const h = observe(access, async () => { ++requests; return access as unknown as MIDIAccess; });
  a.opening = Promise.reject(new Error("Device busy"));
  await h.input.connect();
  assert.equal(h.state.status, "error");
  assert.match(h.state.message, /Device busy/);
  access.changed();
  await tick();
  assert.equal(a.openCalls, 1, "hotplug notifications must not retry a rejected device in a loop");
  a.opening = null;
  await h.input.connect();
  assert.equal(h.state.status, "connected");
  assert.equal(a.openCalls, 2);
  assert.equal(requests, 1);
  h.input.dispose();
});

test("two subscriptions and pre-existing handlers survive each other's disconnect", async () => {
  const a = new Input("a"), access = new Access(a), one = observe(access), two = observe(access);
  let externalMessages = 0, externalChanges = 0;
  const messageHandler = () => { ++externalMessages; }, changeHandler = () => { ++externalChanges; };
  a.onmidimessage = messageHandler;
  access.onstatechange = changeHandler;
  a.addEventListener("midimessage", messageHandler);
  await Promise.all([one.input.connect(), two.input.connect()]);
  a.send([0x90, 60, 80]);
  one.input.disconnect();
  one.input.disconnect();
  a.send([0x80, 60, 0]);
  access.changed();
  assert.equal(one.messages.length, 1);
  assert.equal(two.messages.length, 2);
  assert.equal(one.resets, 1);
  assert.equal(two.resets, 0);
  assert.equal(externalMessages, 4);
  assert.equal(externalChanges, 1);
  assert.equal(a.onmidimessage, messageHandler);
  assert.equal(access.onstatechange, changeHandler);
  assert.equal(a.closeCalls, 0);
  two.input.dispose();
  one.input.dispose();
});

test("default permission request explicitly disables system-exclusive access", async () => {
  const previous = Object.getOwnPropertyDescriptor(globalThis, "navigator"), access = new Access();
  const requests: unknown[] = [];
  Object.defineProperty(globalThis, "navigator", { configurable: true, value: {
    requestMIDIAccess: async (options: unknown) => { requests.push(options); return access; },
  } });
  const input = new LiveMidiInput({ onMessage() {}, onReset() {}, onState() {} });
  try {
    assert.equal(requests.length, 0, "construction never asks for permission");
    await input.connect();
    assert.deepEqual(requests, [{ sysex: false }]);
  } finally {
    input.dispose();
    if (previous) Object.defineProperty(globalThis, "navigator", previous);
    else Reflect.deleteProperty(globalThis, "navigator");
  }
});
