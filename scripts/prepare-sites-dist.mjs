import { access, copyFile, mkdir, readFile } from "node:fs/promises";
import path from "node:path";

const root = process.cwd();
const serverEntry = path.join(root, "dist", "server", "index.js");
const hostingSource = path.join(root, ".openai", "hosting.json");
const hostingDirectory = path.join(root, "dist", ".openai");
const hostingTarget = path.join(hostingDirectory, "hosting.json");

await access(serverEntry);

const hosting = JSON.parse(await readFile(hostingSource, "utf8"));
if (
  typeof hosting.project_id !== "string" ||
  hosting.project_id.length === 0
) {
  throw new Error(".openai/hosting.json has no valid Sites project_id");
}

await mkdir(hostingDirectory, { recursive: true });
await copyFile(hostingSource, hostingTarget);

console.log("Sites bundle prepared at dist/server/index.js");
