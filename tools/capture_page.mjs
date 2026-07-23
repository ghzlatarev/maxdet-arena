#!/usr/bin/env node

import fs from "node:fs";

const [endpoint, url, output, rawWidth, rawHeight, ...flags] =
  process.argv.slice(2);
if (!endpoint || !url || !output || !rawWidth || !rawHeight) {
  console.error(
    "usage: capture_page.mjs CDP_ENDPOINT URL OUTPUT WIDTH HEIGHT [--full-page] [--scroll=Y]",
  );
  process.exit(2);
}

const width = Number(rawWidth);
const height = Number(rawHeight);
const fullPage = flags.includes("--full-page");
const scrollFlag = flags.find((flag) => flag.startsWith("--scroll="));
const scrollY = scrollFlag ? Number(scrollFlag.slice("--scroll=".length)) : 0;
if (
  !Number.isInteger(width) ||
  !Number.isInteger(height) ||
  width <= 0 ||
  height <= 0 ||
  !Number.isFinite(scrollY) ||
  scrollY < 0
) {
  throw new Error("width, height, and scroll offset must be non-negative numbers");
}
const targets = await fetch(`${endpoint}/json`).then((response) => response.json());
const target = targets.find((item) => item.type === "page");
if (!target) throw new Error("no Chrome page target found");

const socket = new WebSocket(target.webSocketDebuggerUrl);
await new Promise((resolve, reject) => {
  socket.addEventListener("open", resolve, { once: true });
  socket.addEventListener("error", reject, { once: true });
});

let sequence = 0;
const pending = new Map();
socket.addEventListener("message", (event) => {
  const message = JSON.parse(event.data);
  if (!message.id) return;
  const handler = pending.get(message.id);
  if (!handler) return;
  pending.delete(message.id);
  if (message.error) handler.reject(new Error(message.error.message));
  else handler.resolve(message.result);
});

function call(method, params = {}) {
  const id = ++sequence;
  socket.send(JSON.stringify({ id, method, params }));
  return new Promise((resolve, reject) => {
    pending.set(id, { resolve, reject });
  });
}

await call("Page.enable");
await call("Emulation.setDeviceMetricsOverride", {
  width,
  height,
  deviceScaleFactor: 1,
  mobile: width < 700,
});
await call("Page.navigate", { url });
await new Promise((resolve) => setTimeout(resolve, 800));
if (scrollY) {
  await call("Runtime.evaluate", {
    expression: `document.documentElement.style.scrollBehavior = "auto"; scrollTo(0, ${JSON.stringify(scrollY)})`,
  });
  await new Promise((resolve) => setTimeout(resolve, 200));
}

let clip;
if (fullPage) {
  const metrics = await call("Page.getLayoutMetrics");
  clip = {
    x: 0,
    y: 0,
    width,
    height: Math.ceil(metrics.cssContentSize.height),
    scale: 1,
  };
}

const screenshot = await call("Page.captureScreenshot", {
  format: "png",
  fromSurface: true,
  captureBeyondViewport: Boolean(clip),
  ...(clip ? { clip } : {}),
});
fs.writeFileSync(output, Buffer.from(screenshot.data, "base64"));

const dimensions = await call("Runtime.evaluate", {
  expression:
    "JSON.stringify({innerWidth,innerHeight,scrollX,scrollY,scrollWidth:document.documentElement.scrollWidth,scrollHeight:document.documentElement.scrollHeight})",
  returnByValue: true,
});
console.log(dimensions.result.value);
socket.close();
