import fs from "node:fs";
import path from "node:path";

export type ProblemStatus = "open" | "solved";

export type Problem = {
  slug: string;
  status: ProblemStatus;
  title: string;
  summary: string;
  field: string;
  date: string;
  agent: string;
  verification: string;
  href: string;
  action: string;
  featured?: boolean;
  count?: number;
  note?: string;
};

const requiredTextFields = [
  "slug",
  "title",
  "summary",
  "field",
  "date",
  "agent",
  "verification",
  "href",
  "action",
] as const;

function assertProblem(value: unknown, filename: string): asserts value is Problem {
  if (!value || typeof value !== "object") {
    throw new Error(`${filename}: problem must be a JSON object`);
  }

  const problem = value as Record<string, unknown>;
  if (problem.status !== "open" && problem.status !== "solved") {
    throw new Error(`${filename}: status must be open or solved`);
  }

  for (const field of requiredTextFields) {
    if (typeof problem[field] !== "string" || problem[field].trim() === "") {
      throw new Error(`${filename}: ${field} must be a non-empty string`);
    }
  }

  if (!/^[a-z0-9]+(?:-[a-z0-9]+)*$/.test(problem.slug as string)) {
    throw new Error(`${filename}: slug must use lowercase kebab-case`);
  }
  if (!/^\d{4}-\d{2}-\d{2}$/.test(problem.date as string)) {
    throw new Error(`${filename}: date must use YYYY-MM-DD`);
  }
  if (
    !(problem.href as string).startsWith("/") &&
    !/^https:\/\//.test(problem.href as string)
  ) {
    throw new Error(`${filename}: href must be an internal path or HTTPS URL`);
  }

  if (problem.featured !== undefined && typeof problem.featured !== "boolean") {
    throw new Error(`${filename}: featured must be a boolean`);
  }
  if (problem.note !== undefined && typeof problem.note !== "string") {
    throw new Error(`${filename}: note must be a string`);
  }

  if (problem.count !== undefined) {
    if (
      typeof problem.count !== "number" ||
      !Number.isInteger(problem.count) ||
      problem.count < 1
    ) {
      throw new Error(`${filename}: count must be a positive integer`);
    }
  }
}

export function loadProblems(): Problem[] {
  const directory = path.join(process.cwd(), "data", "problems");
  const problems = fs
    .readdirSync(directory)
    .filter((filename) => filename.endsWith(".json"))
    .map((filename) => {
      const value: unknown = JSON.parse(
        fs.readFileSync(path.join(directory, filename), "utf8"),
      );
      assertProblem(value, filename);
      if (value.slug !== filename.replace(/\.json$/, "")) {
        throw new Error(`${filename}: filename must match slug`);
      }
      return value;
    });

  const slugs = new Set<string>();
  for (const problem of problems) {
    if (slugs.has(problem.slug)) {
      throw new Error(`duplicate problem slug: ${problem.slug}`);
    }
    slugs.add(problem.slug);
  }

  return problems.sort((left, right) => {
    if (left.featured !== right.featured) return left.featured ? -1 : 1;
    if (left.status !== right.status) return left.status === "open" ? -1 : 1;
    return right.date.localeCompare(left.date);
  });
}
